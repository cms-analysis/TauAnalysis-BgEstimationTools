#include "TauAnalysis/BgEstimationTools/plugins/TemplateBgEstFit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include "TauAnalysis/BgEstimationTools/interface/BgEstMean.h"
#include "TauAnalysis/BgEstimationTools/interface/BgEstMedian.h"
#include "TauAnalysis/BgEstimationTools/interface/BgEstCovMatrix.h"
#include "TauAnalysis/BgEstimationTools/interface/templateBgEstFitAuxFunctions.h"
#include "TauAnalysis/BgEstimationTools/interface/histogramAuxFunctions.h"

#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include <TEllipse.h>
#include <TColor.h>
#include <TMarker.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TArrayD.h>
#include <TLegend.h>
#include <TF1.h>
#include <TStyle.h>

#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooProdPdf.h>
#include <RooArgSet.h>
#include <RooGlobalFunc.h>
#include <RooPlot.h>
#include <RooFit.h>
#include <RooLinkedList.h>
#include <RooCmdArg.h>
#include <RooTFnPdfBinding.h>
#include <RooParametricStepFunction.h>

#include <iostream>
#include <fstream>

const int defaultCanvasSizeX = 800;
const int defaultCanvasSizeY = 600;

const double maxNorm = 1.e+6;
const double normStartValue = 1.e+3;

const std::string fluctMode_coherent = "coherent";
const std::string fluctMode_incoherent = "incoherent";

const double epsilon = 1.e-3;
const double epsilon_integral = 5.e-2;

const int fitStatus_converged = 0;

typedef std::pair<double, double> double_pair;

TemplateBgEstFit::data1dType::data1dType(const std::string& processName, const std::string& varName, 
					 const std::string& meName, RooRealVar* x, bool cutUnfittedRegion)
  : processName_(processName),
    varName_(varName),
    meName_(meName),
    xRef_(x),
    cutUnfittedRegion_(cutUnfittedRegion),
    histogram_(0),
    fluctHistogram_(0),
    error_(0)
{}

TemplateBgEstFit::data1dType::~data1dType()
{
  delete histogram_;
  delete fluctHistogram_;
}

void TemplateBgEstFit::data1dType::initialize()
{
  //std::cout << "<data1dType::initialize>:" << std::endl;

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  me_ = dqmStore.get(meName_);
  if ( !me_ ) {
    edm::LogError ("data1dType") << " Failed to access dqmMonitorElement = " << meName_ << " !!";
    error_ = 1;
    return;
  }

//--- make sure that x-axis range excluded from fit
//    does really have no effect on fit results
//    (differences between the event selection criteria
//     applied in the final analysis and the background enriched samples
//     may cause sizeable deviations between the templates obtained from the background enriched samples 
//     and the distributions observed in the final analysis;
//     even in case regions with deviations are not within the fitted region,
//     they may affect the fit results via differences in the fractions of events 
//     that are outside of the region included in the fit,
//     because the normalization of PDFs depends on these fractions)
  histogram_ = (TH1*)me_->getTH1()->Clone();
  if ( cutUnfittedRegion_ ) {
    const RooAbsBinning& xRange = xRef_->getBinning();
    double xMin = xRange.lowBound();
    double xMax = xRange.highBound();
    int numBins = histogram_->GetNbinsX();
    for ( int iBin = 0; iBin < (numBins + 2); ++iBin ) {
      double binCenter = histogram_->GetBinCenter(iBin);
      
      if ( !(binCenter > xMin && binCenter < xMax) ) {
	histogram_->SetBinContent(iBin, 0.);
	histogram_->SetBinError(iBin, 0.);
      }
    }
  }

  fluctHistogram_ = (TH1*)histogram_->Clone();
}

void TemplateBgEstFit::data1dType::fluctuate(bool, bool)
{  
  sampleHistogram_stat(histogram_, fluctHistogram_);
  makeHistogramPositive(fluctHistogram_);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::dataNdType::dataNdType(bool cutUnfittedRegion)
  : numDimensions_(0),
    cutUnfittedRegion_(cutUnfittedRegion),
    normCorrFactor_(1.),
    error_(0)
{}

TemplateBgEstFit::dataNdType::~dataNdType()
{
  //std::cout << "<dataNdType::~dataNdType>:" << std::endl;

  for ( std::map<std::string, data1dType*>::iterator it = data1dEntries_.begin();
	it != data1dEntries_.end(); ++it ) {
    delete it->second;
  }
}

void TemplateBgEstFit::dataNdType::addVar(const std::string& varName, RooRealVar* x, const std::string& meName)
{
  varNames_.push_back(varName);
  
  data1dType* data1dEntry = new data1dType("data", varName, meName, x, cutUnfittedRegion_);
  if ( data1dEntry->error_ ) error_ = 1;
  data1dEntries_[varName] = data1dEntry;
  
  ++numDimensions_;
}

void TemplateBgEstFit::dataNdType::initialize()
{
  std::cout << "<dataNdType::initialize>:" << std::endl;

  for ( std::map<std::string, data1dType*>::iterator data1dEntry = data1dEntries_.begin();
	data1dEntry != data1dEntries_.end(); ++data1dEntry ) { 
    data1dEntry->second->initialize();
    if ( data1dEntry->second->error_ ) error_ = 1;

    //std::cout << "processName = " << data1dEntry->second->processName_ << std::endl;

    const RooAbsBinning& xRange = data1dEntry->second->xRef_->getBinning();
    double xMin = xRange.lowBound();
    double xMax = xRange.highBound();
    double integral_fitted = getIntegral(data1dEntry->second->histogram_, xMin, xMax);
    //std::cout << " integral_fitted = " << integral_fitted << std::endl;
    
    double integral = getIntegral(data1dEntry->second->me_->getTH1());
    //std::cout << " integral = " << integral << std::endl;
    
    if ( integral_fitted > 0. ) normCorrFactor_ *= (integral/integral_fitted);
  }

  std::cout << "--> using normCorrFactor = " << normCorrFactor_ << std::endl;
}

void TemplateBgEstFit::dataNdType::fluctuate(bool, bool)
{
  for ( std::map<std::string, data1dType*>::iterator data1dEntry = data1dEntries_.begin();
	data1dEntry != data1dEntries_.end(); ++data1dEntry ) {
    data1dEntry->second->fluctuate(true, false);
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::model1dType::model1dType(const std::string& processName, const std::string& varName,
					   const std::string& meName, RooRealVar* x, bool cutUnfittedRegion,
					   bool applySmoothing, const edm::ParameterSet& cfgSmoothing)
  : data1dType(processName, varName, meName, x, cutUnfittedRegion),
    pdfName_(std::string(processName).append("_").append(varName).append("_").append("pdf")),
    pdf_(0),
    applySmoothing_(applySmoothing), cfgSmoothing_(cfgSmoothing),
    auxTF1Wrapper_(0),
    pdfBinning_(0),
    pdfCoeffCollection_(0),
    pdfCoeffArgs_(0)
{}

TemplateBgEstFit::model1dType::~model1dType()
{
  delete pdf_;
  delete auxTF1Wrapper_;
  delete pdfBinning_;
  delete pdfCoeffCollection_;
  delete pdfCoeffArgs_;
}

void TemplateBgEstFit::model1dType::initialize()
{
  //std::cout << "<model1dType::initialize>:" << std::endl;

  data1dType::initialize();

  if ( error_ ) return;

  buildPdf();

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  for ( std::vector<sysErrFluctType>::iterator sysErrFluctuation = sysErrFluctuations_.begin();
	sysErrFluctuation != sysErrFluctuations_.end(); ++sysErrFluctuation ) {
    sysErrFluctuation->me_ = dqmStore.get(sysErrFluctuation->meName_);
    if ( !sysErrFluctuation->me_ ) {
      edm::LogError ("model1dType") << " Failed to access dqmMonitorElement = " << sysErrFluctuation->meName_ << " !!";
      error_ = 1;
      continue;
    }

//--- check that histograms representing systematic uncertainties have the same binning
//    as that representing expectation
    if ( !isCompatibleBinning(histogram_, sysErrFluctuation->me_->getTH1()) ) {
      edm::LogError ("model1dType") << " Incompatible binning of histograms " << meName_ 
				    << " and " << sysErrFluctuation->meName_ << " !!";
      error_ = 1;
      continue;
    }
  }
}

void TemplateBgEstFit::model1dType::fluctuate(bool fluctStat, bool fluctSys)
{
  if ( fluctStat ) {
    sampleHistogram_stat(histogram_, fluctHistogram_);
  } else {
    int numBins = histogram_->GetNbinsX();
    for ( int iBin = 0; iBin < (numBins + 2); ++iBin ) {
      double binContent = histogram_->GetBinContent(iBin);
      double binError = histogram_->GetBinError(iBin);

      fluctHistogram_->SetBinContent(iBin, binContent);
      fluctHistogram_->SetBinError(iBin, binError);
    }
  }

  if ( fluctSys ) {
    for ( std::vector<sysErrFluctType>::const_iterator sysErrFluctuation = sysErrFluctuations_.begin();
	  sysErrFluctuation != sysErrFluctuations_.end(); ++sysErrFluctuation ) {
      sampleHistogram_sys(fluctHistogram_, sysErrFluctuation->me_->getTH1(),  
			  sysErrFluctuation->pullRMS_, sysErrFluctuation->pullMin_, sysErrFluctuation->pullMax_, 
			  sysErrFluctuation->fluctMode_);
    }
  }

  makeHistogramPositive(fluctHistogram_);
  
  buildPdf();
}

void TemplateBgEstFit::model1dType::buildPdf()
{
  std::cout << "<model1dType::buildPdf>:" << std::endl;

  if ( applySmoothing_ ) {
    bool isFirstFit = (!auxTF1Wrapper_);
    
    if ( isFirstFit ) {
      std::string pluginTypeTF1Wrapper = cfgSmoothing_.getParameter<std::string>("pluginType");
      auxTF1Wrapper_ = TF1WrapperPluginFactory::get()->create(pluginTypeTF1Wrapper, cfgSmoothing_);
    } else {
      auxTF1Wrapper_->reinitializeTF1Parameter();
    }

    std::string fitOption = ( isFirstFit ) ? "RB0" : "RB0Q";
    
    fluctHistogram_->Fit(auxTF1Wrapper_->getTF1(), fitOption.data());

    delete pdf_;
    pdf_ = new RooTFnPdfBinding(pdfName_.data(), pdfName_.data(), auxTF1Wrapper_->getTF1(), RooArgList(*xRef_));
  } else {
    bool isFirstFit = (!pdfCoeffCollection_);

    if ( isFirstFit ) {
      pdfBinning_ = new TArrayD(getBinning(fluctHistogram_));
      pdfCoeffCollection_ = new TObjArray();

      const RooAbsBinning& xRange = xRef_->getBinning();
      double xMin = xRange.lowBound();
      double xMax = xRange.highBound();
      std::cout << "pdfName = " << pdfName_ << ": xMin = " << xMin << ", xMax = " << xMax << std::endl;

      unsigned numBins = pdfBinning_->GetSize() - 1;
      for ( unsigned iBin = 0; iBin < numBins; ++iBin ) {
	std::ostringstream pdfCoeffName;
	pdfCoeffName << processName_ << "_" << varName_ << "_coeff" << iBin;	

	RooAbsReal* pdfCoeff = 0;
	if ( pdfBinning_->At(iBin + 1) > xMin && pdfBinning_->At(iBin) < xMax ) {
	  pdfCoeff = new RooRealVar(pdfCoeffName.str().data(), pdfCoeffName.str().data(), 0., 1.);
	} else {
	  pdfCoeff = new RooConstVar(pdfCoeffName.str().data(), pdfCoeffName.str().data(), fluctHistogram_->GetBinContent(iBin));
	}

	pdfCoeffCollection_->Add(pdfCoeff);
      }
      
      pdfCoeffCollection_->SetOwner();

      std::string pdfCoeffArgName = std::string(processName_).append("_").append(varName_).append("_pdfCoeffArgs");
      pdfCoeffArgs_ = new RooArgList(*pdfCoeffCollection_, pdfCoeffArgName.data());
    }

    unsigned numBins = pdfCoeffCollection_->GetEntries();
    for ( unsigned iBin = 0; iBin < numBins; ++iBin ) {
      RooAbsRealLValue* pdfCoeff = dynamic_cast<RooAbsRealLValue*>(pdfCoeffCollection_->At(iBin));
      if ( pdfCoeff && !pdfCoeff->isConstant() ) {
	double pdfCoeffValue = 0.25*(3.*fluctHistogram_->GetBinContent(iBin) + 1./numBins);
	pdfCoeff->setVal(pdfCoeffValue);
      }
    }
    
    delete pdf_;
    pdf_ = new RooParametricStepFunction(pdfName_.data(), pdfName_.data(), *xRef_, *pdfCoeffArgs_, *pdfBinning_, numBins);
  } 
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::modelNdType::modelNdType(const std::string& processName, const edm::ParameterSet& cfgProcess,
					   bool cutUnfittedRegion, 
					   bool applyNormConstraint, double valueNormConstraint, double sigmaNormConstraint)
  : processName_(processName),
    applyNormConstraint_(applyNormConstraint),
    pdfNormConstraint_(0), 
    meanNormConstraint_(0),
    sigmaNormConstraint_(0),
    norm_(0),
    cutUnfittedRegion_(cutUnfittedRegion),
    normCorrFactor_(1.),
    numDimensions_(0),
    error_(0)
{
  std::string normName = std::string(processName_).append("_").append("norm");
  double norm_initial = ( cfgProcess.exists("norm") ) ? 
    cfgProcess.getParameter<edm::ParameterSet>("norm").getParameter<double>("initial") : normStartValue;
  norm_ = new RooRealVar(normName.data(), normName.data(), norm_initial, 0., maxNorm);

  if ( applyNormConstraint_ ) {
    std::cout << "<modelNdType>: constraining norm = " << valueNormConstraint << " +/- " << sigmaNormConstraint << ","
	      << " process = " << processName_ << std::endl;
    
    std::string meanNormConstraintName = std::string(processName_).append("_").append("meanNormConstraint");
    meanNormConstraint_ = new RooConstVar(meanNormConstraintName.data(), meanNormConstraintName.data(), valueNormConstraint);
    std::string sigmaNormConstraintName = std::string(processName_).append("_").append("sigmaNormConstraint");
    sigmaNormConstraint_ = new RooConstVar(sigmaNormConstraintName.data(), sigmaNormConstraintName.data(), sigmaNormConstraint);
    
    std::string pdfNormConstraintName = std::string(processName_).append("_").append("pdfNormConstraint");
    pdfNormConstraint_ = new RooGaussian(pdfNormConstraintName.data(), pdfNormConstraintName.data(), 
					 *norm_, *meanNormConstraint_, *sigmaNormConstraint_);
  }

  edm::ParameterSet cfgDrawOptions_process = cfgProcess.getParameter<edm::ParameterSet>("drawOptions");
  lineColor_ = cfgDrawOptions_process.getParameter<int>("lineColor");
  lineStyle_ = cfgDrawOptions_process.getParameter<int>("lineStyle");
  lineWidth_ = cfgDrawOptions_process.getParameter<int>("lineWidth");
}

TemplateBgEstFit::modelNdType::~modelNdType()
{
  for ( std::map<std::string, model1dType*>::iterator it = model1dEntries_.begin();
	it != model1dEntries_.end(); ++it ) {
    delete it->second;
  }

  delete pdfNormConstraint_;
  delete meanNormConstraint_;
  delete sigmaNormConstraint_;

  delete norm_;
}

void TemplateBgEstFit::modelNdType::addVar(const std::string& varName, RooRealVar* x, const std::string& meName,
					   bool applySmoothing, const edm::ParameterSet& cfgSmoothing)
{
  varNames_.push_back(varName);
  
  model1dType* model1dEntry = new model1dType(processName_, varName, meName, x, cutUnfittedRegion_,
					      applySmoothing, cfgSmoothing);
  if ( model1dEntry->error_ ) error_ = 1;
  model1dEntries_[varName] = model1dEntry;
  
  ++numDimensions_;
}

void TemplateBgEstFit::modelNdType::initialize()
{
  std::cout << "<modelNdType::initialize>:" << std::endl;

  for ( std::map<std::string, model1dType*>::iterator model1dEntry = model1dEntries_.begin();
	model1dEntry != model1dEntries_.end(); ++model1dEntry ) {
    model1dEntry->second->initialize();
    if ( model1dEntry->second->error_ ) error_ = 1;

    std::cout << "processName = " << model1dEntry->second->processName_ << std::endl;
      
    const RooAbsBinning& xRange = model1dEntry->second->xRef_->getBinning();
    double xMin = xRange.lowBound();
    double xMax = xRange.highBound();
    double integral_fitted = getIntegral(model1dEntry->second->histogram_, xMin, xMax);
    std::cout << " integral_fitted = " << integral_fitted << std::endl;
    
    double integral = getIntegral(model1dEntry->second->me_->getTH1());
    std::cout << " integral = " << integral << std::endl;
    
    if ( integral_fitted > 0. ) normCorrFactor_ *= (integral/integral_fitted);
  }

  std::cout << "--> using normCorrFactor = " << normCorrFactor_ << std::endl;
}

void TemplateBgEstFit::modelNdType::fluctuate(bool fluctStat, bool fluctSys)
{
  for ( std::map<std::string, model1dType*>::iterator model1dEntry = model1dEntries_.begin();
	model1dEntry != model1dEntries_.end(); ++model1dEntry ) {
    model1dEntry->second->fluctuate(fluctStat, fluctSys);
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::TemplateBgEstFit(const edm::ParameterSet& cfg)
  : fitData_(0),
    fitModel_(0),
    fitCategories_(0),
    fitResult_(0),
    error_(0)
{
  std::cout << "<TemplateBgEstFit::TemplateBgEstFit>:" << std::endl;

  edm::ParameterSet cfgFit = cfg.getParameter<edm::ParameterSet>("fit");

  std::map<std::string, bool> applyNormConstraints;
  std::map<std::string, double> meanNormConstraints;
  std::map<std::string, double> sigmaNormConstraints;
  if ( cfgFit.exists("constraints") ) {
    edm::ParameterSet cfgProcesses = cfgFit.getParameter<edm::ParameterSet>("constraints");
    vstring processNames = cfgProcesses.getParameterNamesForType<edm::ParameterSet>();
    for ( vstring::const_iterator processName = processNames.begin(); 
	  processName != processNames.end(); ++processName ) {
      edm::ParameterSet cfgProcess = cfgProcesses.getParameter<edm::ParameterSet>(*processName);
      
      edm::ParameterSet cfgNorm = cfgProcess.getParameter<edm::ParameterSet>("norm");

      applyNormConstraints[*processName] = true;
      meanNormConstraints[*processName] = cfgNorm.getParameter<double>("value");
      sigmaNormConstraints[*processName] = cfgNorm.getParameter<double>("uncertainty");
    }
  }

  cutUnfittedRegion_ = ( cfgFit.exists("cutUnfittedRegion") ) ? cfgFit.getParameter<bool>("cutUnfittedRegion") : false;

  printLevel_ = ( cfgFit.exists("printLevel") ) ? cfgFit.getParameter<int>("printLevel") : 1;
  printWarnings_  = ( cfgFit.exists("printWarnings") ) ? cfgFit.getParameter<bool>("printWarnings") : true;

//--- read list of variables to be used in fit
//    (for each variable: name, title and range to be fitted)
  edm::ParameterSet cfgVariables = cfgFit.getParameter<edm::ParameterSet>("variables");
  vstring varNames = cfgVariables.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator varName = varNames.begin(); 
	varName != varNames.end(); ++varName ) {
    edm::ParameterSet cfgVariable = cfgVariables.getParameter<edm::ParameterSet>(*varName);

    std::cout << "--> including variable = " << (*varName) << " in fit." << std::endl;
    varNames_.push_back(*varName);

    std::string name = cfgVariable.getParameter<std::string>("name");
    std::string title = cfgVariable.getParameter<std::string>("title");
    double xMin = cfgVariable.getParameter<double>("xMin");
    double xMax = cfgVariable.getParameter<double>("xMax");

    x_[*varName] = new RooRealVar(name.data(), title.data(), xMin, xMax);
  }

//--- read configuration parameters specifying signal and background templates
//    fitted to the distributions observed in (pseudo)data
//    (for each process: name; for each combination of variable and process: name of DQM MonitorElement holding template histogram)
  edm::ParameterSet cfgProcesses = cfg.getParameter<edm::ParameterSet>("processes");
  vstring processNames = cfgProcesses.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator processName = processNames.begin(); 
	processName != processNames.end(); ++processName ) {
    edm::ParameterSet cfgProcess = cfgProcesses.getParameter<edm::ParameterSet>(*processName);

    processNames_.push_back(*processName);

    modelNdType* modelEntry = new modelNdType(*processName, cfgProcess, cutUnfittedRegion_,
					      applyNormConstraints[*processName], 
					      meanNormConstraints[*processName], sigmaNormConstraints[*processName]);
    
    edm::ParameterSet cfgTemplates = cfgProcess.getParameter<edm::ParameterSet>("templates");
    for ( vstring::const_iterator varName = varNames_.begin();
	  varName != varNames_.end(); ++varName ) {
      if ( !cfgTemplates.exists(*varName) ) {
	edm::LogError ("TemplateBgEstFit") << " No Template of variable = " << (*varName) 
					   << " defined for process = " << (*processName) << " !!";
	error_ = 1;
	continue;
      }

      edm::ParameterSet cfgTemplate = cfgTemplates.getParameter<edm::ParameterSet>(*varName);

      std::string meName = cfgTemplate.getParameter<std::string>("meName");

      bool applySmoothing = cfgTemplate.exists("smoothing");
      edm::ParameterSet cfgSmoothing = ( applySmoothing ) ? 
	cfgTemplate.getParameter<edm::ParameterSet>("smoothing") : edm::ParameterSet();

      modelEntry->addVar(*varName, x_[*varName], meName, applySmoothing, cfgSmoothing);
    }

    modelEntries_[*processName] = modelEntry;
  }
  
//--- read configuration parameters specifying distributions observed in (pseudo)data 
//    that are to be fitted
//    (for each variable: name of DQM MonitorElements holding template histogram)
//
//    WARNING: dataNdType::addElement and modelNdType::addElement need to be called
//              with variable names sorted in the **exact** same order !!
//
  edm::ParameterSet cfgData = cfg.getParameter<edm::ParameterSet>("data");

  dataEntry_ = new dataNdType(cutUnfittedRegion_);

  edm::ParameterSet cfgDistributions = cfgData.getParameter<edm::ParameterSet>("distributions");
  for ( vstring::const_iterator varName = varNames_.begin();
        varName != varNames_.end(); ++varName ) {
    if ( !cfgDistributions.exists(*varName) ) {
      edm::LogError ("TemplateBgEstFit") << " No Template of variable = " << (*varName) << " defined for data !!";
      error_ = 1;
      continue;
    }

    edm::ParameterSet cfgDistribution = cfgDistributions.getParameter<edm::ParameterSet>(*varName);

    std::string meName = cfgDistribution.getParameter<std::string>("meName");

    dataEntry_->addVar(*varName, x_[*varName], meName);
  }

//--- read configuration parameters specifying options for making control plots
  edm::ParameterSet cfgControlPlots = cfg.getParameter<edm::ParameterSet>("output").getParameter<edm::ParameterSet>("controlPlots");
  controlPlotsFileName_ = cfgControlPlots.getParameter<std::string>("fileName");

//--- read configuration parameters specifying how statistical and systematic uncertainties 
//    on normalization factors determined by fit get estimated
  edm::ParameterSet cfgStatErr = cfg.getParameter<edm::ParameterSet>("estStatUncertainties");
  statErrNumSamplings_ = cfgStatErr.getParameter<int>("numSamplings");
  statErrChi2redMax_ = cfgStatErr.getParameter<double>("chi2redMax");  
  statErrPrintLevel_ = cfgStatErr.getParameter<edm::ParameterSet>("verbosity").getParameter<int>("printLevel");
  statErrPrintWarnings_ = cfgStatErr.getParameter<edm::ParameterSet>("verbosity").getParameter<bool>("printWarnings");

  edm::ParameterSet cfgSysErr = cfg.getParameter<edm::ParameterSet>("estSysUncertainties");
  edm::ParameterSet cfgSysFluctuations = cfgSysErr.getParameter<edm::ParameterSet>("fluctuations");
  vstring sysErrFluctNames = cfgSysFluctuations.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator sysErrFluctName = sysErrFluctNames.begin(); 
	sysErrFluctName != sysErrFluctNames.end(); ++sysErrFluctName ) {
    edm::ParameterSet cfgSysFluct = cfgSysFluctuations.getParameter<edm::ParameterSet>(*sysErrFluctName);

    double pullRMS = cfgSysFluct.getParameter<double>("pullRMS");
    double pullMin = cfgSysFluct.getParameter<double>("pullMin");
    double pullMax = cfgSysFluct.getParameter<double>("pullMax");

    std::string fluctMode_string = cfgSysFluct.getParameter<std::string>("mode");
    int fluctMode_int = -1;
    if ( fluctMode_string == fluctMode_coherent ) {
      fluctMode_int = kCoherent;
    } else if ( fluctMode_string == fluctMode_incoherent ) {
      fluctMode_int = kIncoherent;
    } else {
      edm::LogError ("TemplateBgEstFit::TemplateBgEstFit") << " Invalid 'mode' parameter = " << fluctMode_string << " !!";
      error_ = 1;
    }

    edm::ParameterSet cfgProcesses = cfgSysFluct.getParameter<edm::ParameterSet>("meNames");

    for ( vstring::iterator processName = processNames_.begin();
	  processName != processNames_.end(); ++processName ) {
      if ( !cfgProcesses.exists(*processName) ) {
	edm::LogError ("TemplateBgEstFit") << " No Estimate of systematic Uncertainty = " << (*sysErrFluctName) 
					   << " defined for process = " << (*processName) << " !!";
	error_ = 1;
	continue;
      }

      edm::ParameterSet cfgProcess = cfgProcesses.getParameter<edm::ParameterSet>(*processName);

      for ( vstring::const_iterator varName = varNames_.begin();
	    varName != varNames_.end(); ++varName ) {
	if ( !cfgProcess.exists(*varName) ) {
	  edm::LogError ("TemplateBgEstFit") << " No Estimate of systematic Uncertainty = " << (*sysErrFluctName) 
					     << " on variable = " << (*varName) << " defined for process = " << (*processName) << " !!";
	  error_ = 1;
	}

	std::string meName = cfgProcess.getParameter<std::string>(*varName);

	model1dType::sysErrFluctType sysErrFluct;
	sysErrFluct.fluctName_ = (*sysErrFluctName);
	sysErrFluct.meName_ = meName;
	sysErrFluct.pullRMS_ = pullRMS;
	sysErrFluct.pullMin_ = pullMin;
	sysErrFluct.pullMax_ = pullMax;
	sysErrFluct.fluctMode_ = fluctMode_int;

	modelEntries_[*processName]->model1dEntries_[*varName]->sysErrFluctuations_.push_back(sysErrFluct);
      }
    }
  }
  sysErrNumSamplings_ = cfgSysErr.getParameter<int>("numSamplings");
  sysErrChi2redMax_ = cfgSysErr.getParameter<double>("chi2redMax");  
  sysErrPrintLevel_ = cfgSysErr.getParameter<edm::ParameterSet>("verbosity").getParameter<int>("printLevel");
  sysErrPrintWarnings_ = cfgSysErr.getParameter<edm::ParameterSet>("verbosity").getParameter<bool>("printWarnings");
}

TemplateBgEstFit::~TemplateBgEstFit()
{
  std::cout << "<TemplateBgEstFit::~TemplateBgEstFit>:" << std::endl;

  delete dataEntry_;
  delete fitData_;

  for ( modelEntryMap::iterator it = modelEntries_.begin();
	it != modelEntries_.end(); ++it ) {
    delete it->second;
  }

  for ( pdfMap::iterator it = pdfModelSums_.begin();
	it != pdfModelSums_.end(); ++it ) {
    delete it->second;
  }
  
  for ( normMap::iterator it = normTemplateShapes_.begin();
	it != normTemplateShapes_.end(); ++it ) {
    delete it->second;
  }
  
  for ( pdfMap::iterator it = pdfTemplateShapeSums_.begin();
	it != pdfTemplateShapeSums_.end(); ++it ) {
    delete it->second;
  }
  
  delete fitModel_;
  
  delete fitCategories_;
  
  for ( realVarMap::iterator it = x_.begin();
	it != x_.end(); ++it ) {
    delete it->second;
  }
  
  delete fitResult_;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

std::string getCategoryName_data(const std::string& varName)
{
  return std::string(varName).append("_data");
}

std::string getCategoryName_template(const std::string& processName, const std::string& varName)
{
  return std::string(processName).append("_").append(varName).append("_template");
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void TemplateBgEstFit::buildFitData()
{
  std::string fitDataName = "fitData";

  if ( fitCategories_->numTypes() == 1 ) {
    delete fitData_;
    RooRealVar* x = x_[varNames_.front()];
    TH1* hist = dataEntry_->data1dEntries_[varNames_.front()]->fluctHistogram_;
    fitData_ = new RooDataHist(fitDataName.data(), fitDataName.data(), *x, hist);
  } else {
    TObjArray xCollection;
    
    std::map<std::string, TH1*> histMap;
    
    for ( vstring::const_iterator varName = varNames_.begin();
	  varName != varNames_.end(); ++varName ) {
      
      xCollection.Add(x_[*varName]);
      
      std::string categoryName_data = getCategoryName_data(*varName);
      histMap[categoryName_data] = dataEntry_->data1dEntries_[*varName]->fluctHistogram_;
      
      for ( vstring::const_iterator processName = processNames_.begin();
	    processName != processNames_.end(); ++processName ) {
	if ( !modelEntries_[*processName]->model1dEntries_[*varName]->applySmoothing_ ) {
	  std::string categoryName_template = getCategoryName_template(*processName, *varName);
	  histMap[categoryName_template] = modelEntries_[*processName]->model1dEntries_[*varName]->fluctHistogram_;
	}
      }
    }
    
    for ( std::map<std::string, TH1*>::const_iterator histMapEntry = histMap.begin();
	  histMapEntry != histMap.end(); ++histMapEntry ) {
      std::cout << "histMap[" << histMapEntry->first << "] = " << histMapEntry->second->GetName() << std::endl;
    }
    
    delete fitData_;
    fitData_ = new RooDataHist(fitDataName.data(), fitDataName.data(), RooArgList(xCollection), *fitCategories_, histMap); 
  }
}

void TemplateBgEstFit::buildFitModel()
{
  for ( vstring::const_iterator varName = varNames_.begin();
	varName != varNames_.end(); ++varName ) {
    std::string pdfModelSumName = std::string(*varName).append("_pdfModelSum");
    
    TObjArray pdfModelSum_pdfCollection;
    TObjArray pdfModelSum_normCollection;
    for ( vstring::const_iterator processName = processNames_.begin();
	  processName != processNames_.end(); ++processName ) {
      pdfModelSum_pdfCollection.Add(modelEntries_[*processName]->model1dEntries_[*varName]->pdf_);
      pdfModelSum_normCollection.Add(modelEntries_[*processName]->norm_);
    }
    
    //std::cout << "--> creating RooAddPdf with name = " << pdfModelSumName << std::endl;
    //std::cout << " #pdfs = " << RooArgList(pdfModelSum_pdfCollection).getSize() << std::endl;
    //std::cout << " #coefficients = " << RooArgList(pdfModelSum_normCollection).getSize() << std::endl;
    delete pdfModelSums_[*varName];
    pdfModelSums_[*varName] = new RooAddPdf(pdfModelSumName.data(), pdfModelSumName.data(), 
					    RooArgList(pdfModelSum_pdfCollection), RooArgList(pdfModelSum_normCollection));
  }
  
  std::string fitModelName = "fitModel";

  if ( fitCategories_->numTypes() == 1 ) {
    fitModel_ = pdfModelSums_[varNames_.front()];
  } else {
    delete fitModel_;
    fitModel_ = new RooSimultaneous(fitModelName.data(), fitModelName.data(), *fitCategories_);
    
    for ( vstring::const_iterator varName = varNames_.begin();
	  varName != varNames_.end(); ++varName ) {
      
      std::string categoryName_data = getCategoryName_data(*varName);
      ((RooSimultaneous*)fitModel_)->addPdf(*pdfModelSums_[*varName], categoryName_data.data());
      
      for ( vstring::const_iterator processName = processNames_.begin();
	    processName != processNames_.end(); ++processName ) {
	if ( !modelEntries_[*processName]->model1dEntries_[*varName]->applySmoothing_ ) {
	  std::string key = std::string(*processName).append("_").append(*varName);
	  
	  std::string normTemplateShapeName = std::string(*processName).append("_").append(*varName).append("_normTemplateShape");
	  delete normTemplateShapes_[key];
	  RooRealVar* normTemplateShape = new RooRealVar(normTemplateShapeName.data(), normTemplateShapeName.data(), 0., maxNorm);
	  normTemplateShape->setVal(getIntegral(modelEntries_[*processName]->model1dEntries_[*varName]->fluctHistogram_));
	  normTemplateShapes_[key] = normTemplateShape;
	  
	  RooAbsPdf* pdfTemplateShape = modelEntries_[*processName]->model1dEntries_[*varName]->pdf_;
	  
	  std::string pdfTemplateShapeSumName = std::string(*processName).append("_").append(*varName).append("_pdfTemplateShapeSum");
	  
	  //std::cout << "--> creating RooAddPdf with name = " << pdfTemplateShapeSumName << std::endl;
	  //std::cout << " #pdfs = " << RooArgList(*pdfTemplateShape).getSize() << std::endl;
	  //std::cout << " #coefficients = " << RooArgList(*normTemplateShape).getSize() << std::endl;
	  delete pdfTemplateShapeSums_[key];
	  RooAbsPdf* pdfTemplateShapeSum = new RooAddPdf(pdfTemplateShapeSumName.data(), pdfTemplateShapeSumName.data(), 
							 RooArgList(*pdfTemplateShape), RooArgList(*normTemplateShape));
	  pdfTemplateShapeSums_[key] = pdfTemplateShapeSum;
	  
	  std::string categoryName_template = getCategoryName_template(*processName, *varName);
	  ((RooSimultaneous*)fitModel_)->addPdf(*pdfTemplateShapeSum, categoryName_template.data());
	}
      }
    }
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void TemplateBgEstFit::fit(bool saveFitResult, int printLevel, int printWarnings) 
{
//--- fit model to data

//--- build list of fit options
  RooLinkedList fitOptions;

//--- check if "external" constraints exist on normalization factors to be determined by fit
//    (specified by Gaussian probability density functions with mean and sigma obtained
//     e.g. by level of agreement between Monte Carlo simulation and number of events observed in background enriched samples)
  TObjArray normConstraints_pdfCollection;
  for ( modelEntryMap::iterator modelEntry = modelEntries_.begin();
	modelEntry != modelEntries_.end(); ++modelEntry ) {
    if ( modelEntry->second->applyNormConstraint_ ) normConstraints_pdfCollection.Add(modelEntry->second->pdfNormConstraint_);
  }

  if ( normConstraints_pdfCollection.GetEntries() > 0 ) {
    std::string normConstraints_pdfArgName = std::string("normConstraints").append("_pdfArgs");
    RooArgSet normConstraints_pdfArgs(normConstraints_pdfCollection, normConstraints_pdfArgName.data());
    
    fitOptions.Add(new RooCmdArg(RooFit::ExternalConstraints(normConstraints_pdfArgs)));
  }

//--- check if results of fit are to be saved for later analysis
//    (only necessary when fitting for "central values", **not** when estimating statistical or systematic uncertainties)
  if ( saveFitResult ) fitOptions.Add(new RooCmdArg(RooFit::Save(true)));

//--- stop Minuit from printing lots of information 
//    about progress on fit and warnings
  fitOptions.Add(new RooCmdArg(RooFit::PrintLevel(printLevel)));
  fitOptions.Add(new RooCmdArg(RooFit::PrintEvalErrors(printWarnings)));
  fitOptions.Add(new RooCmdArg(RooFit::Warnings(printWarnings_)));

  RooFitResult* fitResult = fitModel_->fitTo(*fitData_, fitOptions);
  if ( saveFitResult ) fitResult_ = fitResult;

//--- delete fit option objects
  int numCmdArgs = fitOptions.GetSize();
  for ( int iCmdArg = 0; iCmdArg < numCmdArgs; ++iCmdArg ) {
    delete fitOptions.At(iCmdArg);
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void TemplateBgEstFit::endJob()
{
  std::cout << "<TemplateBgEstFit::endJob>:" << std::endl;

//--- check that DQMStore service is available
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("TemplateBgEstFit::endJob") << " Failed to access dqmStore" 
					       << " --> histogram will NOT be fitted !!";
    return;
  }

//--- check that configuration parameters contain no errors,
//    retrieve MonitorElements from DQMStore
//    and check that all DQM MonitorElements have successfully been retrieved
  for ( modelEntryMap::iterator modelEntry = modelEntries_.begin();
	modelEntry != modelEntries_.end(); ++modelEntry ) {
    modelEntry->second->initialize();
    if ( modelEntry->second->error_ ) error_ = 1;
  }

  dataEntry_->initialize();
  if ( dataEntry_->error_ ) error_ = 1;
  
  if ( error_ ) {
    edm::LogError ("TemplateBgEstFit::endJob") << " Error in Configuration ParameterSet" 
					       << " --> histogram will NOT be fitted !!";
    return;
  }

//--- configure RooFit structure;
//    print-out structure once configuration finished
  fitCategories_ = new RooCategory("fitCategories", "fitCategories");
  for ( vstring::const_iterator varName = varNames_.begin();
	varName != varNames_.end(); ++varName ) {
    std::string categoryName_data = getCategoryName_data(*varName);
    fitCategories_->defineType(categoryName_data.data());
    
    for ( vstring::const_iterator processName = processNames_.begin();
	  processName != processNames_.end(); ++processName ) {
      if ( !modelEntries_[*processName]->model1dEntries_[*varName]->applySmoothing_ ) {
	std::string categoryName_template = getCategoryName_template(*processName, *varName);
	fitCategories_->defineType(categoryName_template.data());
      }
    }
  }

  std::cout << "number of Categories = " << fitCategories_->numTypes() << std::endl;
  unsigned numCategories = fitCategories_->numTypes();
  for ( unsigned iCategory = 0; iCategory < numCategories; ++iCategory ) {
    std::cout << "Category[" << iCategory << "] = " << fitCategories_->lookupType(iCategory)->GetName() << std::endl;
  }

  buildFitData();

  buildFitModel();

  std::cout << ">>> RootFit model used for Template method Fit <<<" << std::endl;
  fitModel_->printCompactTree();

  std::cout << ">>> RootFit Parameters <<<" << std::endl;
  for ( vstring::const_iterator varName = varNames_.begin();
	varName != varNames_.end(); ++varName ) {
    std::cout << "for Variable = " << (*varName) << ":" << std::endl;
    RooAbsPdf* pdfModelSum = pdfModelSums_[*varName];
    RooArgSet* pdfModelSumParameters = pdfModelSum->getParameters(*pdfModelSum->getComponents());
    pdfModelSumParameters->Print("v");
    delete pdfModelSumParameters;
  }

  std::cout << ">>> RootFit Observables <<<" << std::endl;
  for ( vstring::const_iterator varName = varNames_.begin();
	varName != varNames_.end(); ++varName ) {
    std::cout << "for Variable = " << (*varName) << ":" << std::endl;
    RooAbsPdf* pdfModelSum = pdfModelSums_[*varName];
    RooArgSet* pdfModelSumObservables = pdfModelSum->getObservables(*pdfModelSum->getComponents());
    pdfModelSumObservables->Print("v");
    delete pdfModelSumObservables;
  }

//--- fit template shapes of different signal and background processes
//    to distribution observed in (pseudo)data
  fit(true, printLevel_, printWarnings_);
  assert(fitResult_ != NULL);

//--- print-out fit results
  std::cout << ">>> Fit Results <<<" << std::endl;
  std::cout << " fitStatus = " << fitResult_->status() << std::endl;
  std::cout << " Chi2red = " << compChi2red() << std::endl;
  print(std::cout);

//--- save normalization factors determined and covariance matrix estimated by fit
  unpackFitResult(fitResult_, fitResultMean_, fitResultCov_);

//--- produce plot of different signal and background processes
//    using scale factors determined by fit
//    compared to distribution of (pseudo)data
  if ( controlPlotsFileName_ != "" ) makeControlPlots();

//--- estimate statistical uncertainties
  std::cout << ">>> Statistical Uncertainties <<<" << std::endl;
  estimateUncertainties(true, false, statErrNumSamplings_, statErrChi2redMax_, 
			"estStatUncertainties", statErrPrintLevel_, statErrPrintWarnings_);
 
//--- estimate systematic uncertainties
  std::cout << ">>> Systematic Uncertainties <<<" << std::endl;
  estimateUncertainties(false, true, sysErrNumSamplings_, sysErrChi2redMax_,
			"estSysUncertainties", sysErrPrintLevel_, sysErrPrintWarnings_);

//--- estimate total (statistical + systematic) uncertainties
  std::cout << ">>> Total (statistical + systematic) Uncertainties <<<" << std::endl;
  double chi2redMax = TMath::Max(statErrChi2redMax_, sysErrChi2redMax_);
  int totErrPrintLevel = TMath::Min(statErrPrintLevel_, sysErrPrintLevel_);
  bool totErrPrintWarnings = (statErrPrintWarnings_ && sysErrPrintWarnings_);
  estimateUncertainties(true, true, sysErrNumSamplings_, chi2redMax,
			"estTotUncertainties", totErrPrintLevel, totErrPrintWarnings);
  
  std::cout << "done." << std::endl;
}

void TemplateBgEstFit::print(std::ostream& stream)
{
  stream << "Fit Parameter:" << std::endl;
  for ( modelEntryMap::iterator modelEntry = modelEntries_.begin();
	modelEntry != modelEntries_.end(); ++modelEntry ) {
    const std::string& processName = modelEntry->first;
    RooRealVar* processNorm = modelEntry->second->norm_;

//--- correct normalization factor determined by fit
//    for events events passing final analysis selection criteria
//    that are outside fitted region
    double normCorrFactor = dataEntry_->normCorrFactor_;

    stream << " " << processName << ": normalization = " << normCorrFactor*processNorm->getVal()
	   << " (within fitted region = " << processNorm->getVal() << ")";

    if ( processNorm->hasAsymError() ) {
      stream << " + " << normCorrFactor*processNorm->getAsymErrorHi() << " (" << processNorm->getAsymErrorHi() << ")"
	     << " - " << fabs(normCorrFactor*processNorm->getAsymErrorLo()) << " (" << fabs(processNorm->getAsymErrorLo()) << ")";
    } else if ( processNorm->hasError() ) {
      stream << " +/- " << normCorrFactor*processNorm->getError() << " (" << processNorm->getError() << ")";
    }

    stream << std::endl;
  }
}

void TemplateBgEstFit::makeControlPlots()
{
//--- produce control plot of distribution observed in (pseudo)data
//    versus sum of signal and background templates using normalization determined by fit

//--- stop ROOT from opening X-window for canvas output
//    (in order to be able to run in batch mode) 
  gROOT->SetBatch(true);

//--- produce control plots of smoothing functions
//    fitted to shape templates in order to reduce effect of statistical fluctuations
//    of shape templates on fit results
//    (and in particular effect of different statistical precision with which shape templates
//     are determined for different background processes in background enriched samples)
  makeControlPlotsSmoothing();

//--- produce control plots of one and two sigma error contours 
//    showing correlation of estimated normalization factors
  const RooArgList& fitParameter = fitResult_->floatParsFinal();
  int numFitParameter = fitParameter.getSize();
  vstring labels(numFitParameter);
  for ( int iX = 0; iX < numFitParameter; ++iX ) {
    const RooAbsArg* paramX_arg = fitParameter.at(iX);
    labels[iX] = paramX_arg->GetName();
  }

  makeControlPlotsCovariance(fitResultMean_, fitResultMean_, fitResultCov_, labels, controlPlotsFileName_, "");

//--- produce control plots of distributions observed in (pseudo)data
//    compared to sum of signal and background templates
//    scaled by normalization factors determined by the fit
  makeControlPlotsObsDistribution();
}

void TemplateBgEstFit::makeControlPlotsSmoothing()
{
  TCanvas canvas("TemplateBgEstFit", "TemplateBgEstFit", defaultCanvasSizeX, defaultCanvasSizeY);
  canvas.SetFillColor(10);
  canvas.SetFrameFillColor(10);

  int defStyle_optStat = gStyle->GetOptStat();
  int defStyle_optFit = gStyle->GetOptStat();
  float defStyle_labelSizeX = gStyle->GetLabelSize("x");
  float defStyle_labelSizeY = gStyle->GetLabelSize("y");

  //gStyle->SetOptStat(1111);
  //gStyle->SetOptFit(111);
  gStyle->SetLabelSize(0.03, "x");
  gStyle->SetLabelSize(0.03, "y");
  
  for ( modelEntryMap::iterator modelEntry = modelEntries_.begin();
	modelEntry != modelEntries_.end(); ++modelEntry ) {
    const std::string& processName = modelEntry->first;
    
    for ( vstring::const_iterator varName = varNames_.begin(); 
	  varName != varNames_.end(); ++varName ) {

      model1dType* model1dEntry = modelEntry->second->model1dEntries_[*varName];

      if ( !model1dEntry->applySmoothing_ ) continue;

      std::string histogramName = std::string(model1dEntry->fluctHistogram_->GetName()).append("_cloned");
      TH1* histogram_cloned = (TH1*)model1dEntry->fluctHistogram_->Clone(histogramName.data());
      histogram_cloned->SetStats(false);
      histogram_cloned->GetXaxis()->SetTitle(varName->data());
      histogram_cloned->SetLineStyle(1);
      histogram_cloned->SetLineColor(1);
      histogram_cloned->SetLineWidth(2);

      std::string tf1Name = std::string(model1dEntry->auxTF1Wrapper_->getTF1()->GetName()).append("_cloned");
      TF1* tf1_cloned = (TF1*)model1dEntry->auxTF1Wrapper_->getTF1()->Clone(tf1Name.data());
      tf1_cloned->SetLineStyle(1);
      tf1_cloned->SetLineColor(2);
      tf1_cloned->SetLineWidth(2);
      
      double yMax = histogram_cloned->GetMaximum();
      histogram_cloned->SetMaximum(1.3*yMax);

      TLegend legend(0.63, 0.38, 0.89, 0.54);
      legend.SetBorderSize(0);
      legend.SetFillColor(0);

      legend.AddEntry(histogram_cloned, "Template Histogram", "p");
      legend.AddEntry(tf1_cloned, "Smoothing Function", "l");

      histogram_cloned->Draw("e1p");
      tf1_cloned->Draw("lsame");
      
      //legend.Draw();

      canvas.Update();

      int errorFlag = 0;
      std::string fileNameParam = std::string(processName).append("_").append(*varName).append("_smoothing");
      std::string fileName = replace_string(controlPlotsFileName_, plotKeyword, fileNameParam, 1, 1, errorFlag);
      if ( !errorFlag ) {
	canvas.Print(fileName.data());
      } else {
	edm::LogError("TemplateBgEstFit::makeControlPlotsSmoothing") << " Failed to decode controlPlotsFileName = " 
								     << controlPlotsFileName_ << " --> skipping !!";
	return;
      }

      delete histogram_cloned;
      delete tf1_cloned;
    }
  }

  gStyle->SetOptStat(defStyle_optStat);
  gStyle->SetOptFit(defStyle_optFit);
  gStyle->SetLabelSize(defStyle_labelSizeX, "x");
  gStyle->SetLabelSize(defStyle_labelSizeY, "y");
}

void drawErrorEllipse(double x0, double y0, double errX0, double errY0, double Sxx, double Sxy, double Syy, 
		      const char* labelX, const char* labelY, const char* fileName)
{
//--- draw best fit estimate (x0,y0) 
//    plus one and two sigma error contours centered at (errX0,errY0) 
//    and with (correlated) uncertainties estimated by elements Sxx, Sxy, Syy of covariance matrix 
//    passed as function arguments
//    (note that since the covariance matrix is symmetric, 
//     there is no need to pass element Syx of the covariance matrix)

  TCanvas canvas("drawErrorEllipse", "drawErrorEllipse", 600, 600);
  canvas.SetFillColor(10);
  canvas.SetFrameFillColor(10);

//--- compute angle between first principal axis of error ellipse
//    and x-axis
  double alpha = 0.5*TMath::ATan2(-2*Sxy, Sxx - Syy);

  //std::cout << "alpha = " << alpha*180./TMath::Pi() << std::endl;

  double sinAlpha = TMath::Sin(alpha);
  double cosAlpha = TMath::Cos(alpha);

//--- compute covariance axis in coordinate system
//    defined by principal axes of error ellipse
  double Suu = Sxx*sinAlpha*sinAlpha + 2*Sxy*sinAlpha*cosAlpha + Syy*cosAlpha*cosAlpha;
  double Svv = Sxx*cosAlpha*cosAlpha + 2*Sxy*sinAlpha*cosAlpha + Syy*sinAlpha*sinAlpha;
  
//--- resolve ambiguity which axis represents the first principal axis
//    and which represents the second principal axis
//
//    NOTE: in case Sxy > 0. (correlation of variables X and Y), 
//          the principal axis needs to point in direction of either the first or the third quadrant;
//          in case Sxy < 0. (anti-correlation of variables X and Y), 
//          the principal axis needs to point in direction of either the second or the fourth quadrant.
  double sigmaX_transformed = 0.;
  double sigmaY_transformed = 0.;
  if ( (Sxy >= 0. && TMath::Abs(alpha) <= 0.5*TMath::Pi()) || 
       (Sxy <  0. && TMath::Abs(alpha) >  0.5*TMath::Pi()) ) {
    sigmaX_transformed = TMath::Sqrt(TMath::Max(Suu, Svv));
    sigmaY_transformed = TMath::Sqrt(TMath::Min(Suu, Svv));
  } else {
    sigmaX_transformed = TMath::Sqrt(TMath::Min(Suu, Svv));
    sigmaY_transformed = TMath::Sqrt(TMath::Max(Suu, Svv));
  }

  TEllipse oneSigmaErrorEllipse(errX0, errY0, sigmaX_transformed*1., sigmaY_transformed*1., 0., 360., alpha*180./TMath::Pi()); 
  oneSigmaErrorEllipse.SetFillColor(5);
  oneSigmaErrorEllipse.SetLineColor(44);
  oneSigmaErrorEllipse.SetLineWidth(1);
  TEllipse twoSigmaErrorEllipse(errX0, errY0, sigmaX_transformed*2., sigmaY_transformed*2., 0., 360., alpha*180./TMath::Pi()); 
  TSeqCollection* colors = gROOT->GetListOfColors();
  if ( colors && colors->At(42) ) {
    TColor* orange = (TColor*)colors->At(42);
    orange->SetRGB(1.00,0.80,0.00);
  } else {
    edm::LogWarning ("drawErrorEllipse") << " Failed to access list of Colors from gROOT object"
					 << " --> skipping definition of Color 'orange' !!";
  }
  twoSigmaErrorEllipse.SetFillColor(42);
  twoSigmaErrorEllipse.SetLineColor(44);
  twoSigmaErrorEllipse.SetLineWidth(1);

  TMarker centralValueMarker(x0, y0, 5);
  centralValueMarker.SetMarkerSize(2);

//--- create dummy histogram  
//    defining region to be plotted
  double minX = TMath::Min(errX0 - 2.2*(TMath::Abs(cosAlpha*sigmaX_transformed) + TMath::Abs(sinAlpha*sigmaY_transformed)),
			   x0 - 0.2*(TMath::Abs(cosAlpha*sigmaX_transformed) + TMath::Abs(sinAlpha*sigmaY_transformed)));
  double maxX = TMath::Max(errX0 + 2.8*(TMath::Abs(cosAlpha*sigmaX_transformed) + TMath::Abs(sinAlpha*sigmaY_transformed)),
			   x0 + 0.2*(TMath::Abs(cosAlpha*sigmaX_transformed) + TMath::Abs(sinAlpha*sigmaY_transformed)));
  double minY = TMath::Min(errY0 - 2.2*(TMath::Abs(sinAlpha*sigmaX_transformed) + TMath::Abs(cosAlpha*sigmaY_transformed)),
			   y0 - 0.2*(TMath::Abs(sinAlpha*sigmaX_transformed) + TMath::Abs(cosAlpha*sigmaY_transformed)));
  double maxY = TMath::Max(errY0 + 2.8*(TMath::Abs(sinAlpha*sigmaX_transformed) + TMath::Abs(cosAlpha*sigmaY_transformed)),
			   y0 + 0.2*(TMath::Abs(sinAlpha*sigmaX_transformed) + TMath::Abs(cosAlpha*sigmaY_transformed)));
			   
  if ( TMath::Abs(maxX - minX) < epsilon || 
       TMath::Abs(maxY - minY) < epsilon ) {
    if ( TMath::Abs(maxX - minX) < epsilon ) edm::LogWarning ("drawErrorEllipse") << " Invalid x-range: minX = maxX = " << minX;
    if ( TMath::Abs(maxY - minY) < epsilon ) edm::LogWarning ("drawErrorEllipse") << " Invalid y-range: minY = maxY = " << minY;
    edm::LogWarning ("drawErrorEllipse") << " --> skipping drawing of Error ellipse for labelX = " << labelX << ","
					 << " labelY = " << labelY << " !!";
    return;
  }

//--- create dummy histogram  
  TH2F dummyHistogram("dummyHistogram", "dummyHistogram", 5, minX, maxX, 5, minY, maxY);
  dummyHistogram.SetTitle("");
  dummyHistogram.SetStats(false);
  dummyHistogram.SetXTitle(labelX);
  dummyHistogram.SetYTitle(labelY);
  dummyHistogram.SetTitleOffset(1.35, "Y");

  dummyHistogram.Draw("AXIS");
  
  twoSigmaErrorEllipse.Draw();
  oneSigmaErrorEllipse.Draw();

  centralValueMarker.Draw();

  TLegend legend(0.70, 0.70, 0.89, 0.89);
  legend.SetBorderSize(0);
  legend.SetFillColor(0);
  legend.AddEntry(&centralValueMarker, "best Fit value", "p");
  legend.AddEntry(&oneSigmaErrorEllipse, "1#sigma Contour", "f");
  legend.AddEntry(&twoSigmaErrorEllipse, "2#sigma Contour", "f");

  legend.Draw();

  canvas.Print(fileName);
}

void TemplateBgEstFit::makeControlPlotsCovariance(TVectorD estimate, TVectorD errMean, TMatrixD errCov, const vstring& labels, 
						  const std::string& controlPlotsFileName, const char* type)
{
  int numFitParameter = estimate.GetNoElements();
  assert(errMean.GetNoElements() == numFitParameter);
  assert(errCov.GetNrows() == numFitParameter && errCov.GetNcols() == numFitParameter);
  for ( int iX = 0; iX < numFitParameter; ++iX ) {
    double x0 = estimate(iX);
    double errX0 = errMean(iX);
    double Sxx = errCov(iX, iX);
    const char* labelX = labels[iX].data();

    for ( int iY = 0; iY < iX; ++iY ) {
      double y0 = estimate(iY);
      double errY0 = errMean(iY);
      double Syy = errCov(iY, iY);
      const char* labelY = labels[iY].data();

      double Sxy = errCov(iX, iY);
      std::string fileNameParam = std::string("corr_").append(labelX).append("_vs_").append(labelY);
      if ( type != "" ) fileNameParam.append("_").append(type);
      
      int errorFlag = 0;
      std::string fileName = replace_string(controlPlotsFileName, plotKeyword, fileNameParam, 1, 1, errorFlag);
      if ( !errorFlag ) {
	drawErrorEllipse(x0, y0, errX0, errY0, Sxx, Sxy, Syy, labelX, labelY, fileName.data());
      } else {
	edm::LogError("drawErrorEllipses") << " Failed to decode controlPlotsFileName = " << controlPlotsFileName 
					   << " --> skipping !!";
	return;
      }
    }
  }
}

void TemplateBgEstFit::makeControlPlotsObsDistribution()
{
  TCanvas canvas("TemplateBgEstFit", "TemplateBgEstFit", defaultCanvasSizeX, defaultCanvasSizeY);
  canvas.SetFillColor(10);
  canvas.SetFrameFillColor(10);

  for ( vstring::const_iterator varName = varNames_.begin(); 
	varName != varNames_.end(); ++varName ) {
    TH1* fittedHistogram_sum = 0;

    std::vector<TH1*> fittedHistograms;

    TLegend legend(0.67, 0.63, 0.89, 0.89);
    legend.SetBorderSize(0);
    legend.SetFillColor(0);

    for ( modelEntryMap::iterator modelEntry = modelEntries_.begin();
	  modelEntry != modelEntries_.end(); ++modelEntry ) {
      const std::string& processName = modelEntry->first;
      double processNorm = modelEntry->second->norm_->getVal();

//--- correct normalization factor determined by fit
//    for events events passing final analysis selection criteria
//    that are outside fitted region
      double normCorrFactor = dataEntry_->normCorrFactor_;

      TH1* processShapeHist = modelEntry->second->model1dEntries_[*varName]->me_->getTH1();

      std::string fittedHistogramName_process = std::string(processShapeHist->GetName()).append("_cloned");
      TH1* fittedHistogram_process = (TH1*)processShapeHist->Clone(fittedHistogramName_process.data());

      if ( getIntegral(fittedHistogram_process) > 0. ) {
	fittedHistogram_process->Scale(normCorrFactor*processNorm/getIntegral(fittedHistogram_process));
      }

      fittedHistogram_process->SetLineColor(modelEntry->second->lineColor_);
      fittedHistogram_process->SetLineStyle(modelEntry->second->lineStyle_);
      fittedHistogram_process->SetLineWidth(modelEntry->second->lineWidth_);
      
      fittedHistograms.push_back(fittedHistogram_process);

      legend.AddEntry(fittedHistogram_process, processName.data(), "l");

      if ( !fittedHistogram_sum ) {
	std::string fittedHistogramName_sum = std::string(*varName).append("_histogram_fittedSum");
	fittedHistogram_sum = (TH1*)fittedHistogram_process->Clone(fittedHistogramName_sum.data());
	fittedHistogram_sum->SetStats(false);
	fittedHistogram_sum->GetXaxis()->SetTitle(varName->data());
	fittedHistogram_sum->SetLineColor(1); // black
	fittedHistogram_sum->SetLineStyle(1); // solid
	fittedHistogram_sum->SetLineWidth(fittedHistogram_process->GetLineWidth());
      } else {
	fittedHistogram_sum->Add(fittedHistogram_process);
      }
    }

    TH1* dataHistogram = dataEntry_->data1dEntries_[*varName]->me_->getTH1();
    dataHistogram->SetMarkerStyle(8);
    legend.AddEntry(dataHistogram, "final Evt. Sel.", "p");

    double yMax = TMath::Max(fittedHistogram_sum->GetMaximum(),
			     dataHistogram->GetMaximum());
    fittedHistogram_sum->SetMaximum(1.3*yMax);
    legend.AddEntry(fittedHistogram_sum, "fitted #Sigma", "l");

    fittedHistogram_sum->Draw("hist");

    for ( std::vector<TH1*>::const_iterator fittedHistogram_process = fittedHistograms.begin();
	  fittedHistogram_process != fittedHistograms.end(); ++fittedHistogram_process ) {
      (*fittedHistogram_process)->Draw("histsame");
    }

    dataHistogram->Draw("e1psame");

    legend.Draw();
  
    canvas.Update();

    int errorFlag = 0;
    std::string fileName = replace_string(controlPlotsFileName_, plotKeyword, *varName, 1, 1, errorFlag);
    if ( !errorFlag ) {
      canvas.Print(fileName.data());
    } else {
      edm::LogError("TemplateBgEstFit::makeControlPlotsObsDistribution") << " Failed to decode controlPlotsFileName = " 
									 << controlPlotsFileName_ << " --> skipping !!";
      return;
    }

    delete fittedHistogram_sum;
    for ( std::vector<TH1*>::iterator it = fittedHistograms.begin();
	  it != fittedHistograms.end(); ++it ) {
      delete (*it);
    }
  }
}

double TemplateBgEstFit::compChi2red() 
{
  //std::cout << "<TemplateBgEstFit::compChi2red>:" << std::endl;
  
  double chi2 = 0.;
  int numDoF = 0;

  for ( vstring::const_iterator varName = varNames_.begin(); 
	varName != varNames_.end(); ++varName ) {
    const TH1* histogramData = dataEntry_->data1dEntries_[*varName]->me_->getTH1();
    const RooAbsBinning& xRange = x_[*varName]->getBinning();
    double xMin = xRange.lowBound();
    double xMax = xRange.highBound();
    int numBins = histogramData->GetNbinsX();
    for ( int iBin = 0; iBin < (numBins + 2); ++iBin ) {
      double dataBinCenter = histogramData->GetBinCenter(iBin);

//--- restrict computation of chi^2 to region included in fit
      if ( !(dataBinCenter > xMin && dataBinCenter < xMax) ) continue;

      double dataBinContent = histogramData->GetBinContent(iBin);
      double dataBinError = histogramData->GetBinError(iBin);
      
      double fitBinContent = 0.;
      double fitBinError2 = 0.;
      for ( modelEntryMap::iterator modelEntry = modelEntries_.begin();
	    modelEntry != modelEntries_.end(); ++modelEntry ) {
	const TH1* histogramProcess = modelEntry->second->model1dEntries_[*varName]->me_->getTH1();

	double processBinContent = histogramProcess->GetBinContent(iBin);
	double processBinError = histogramProcess->GetBinError(iBin);

	double processNorm = modelEntry->second->norm_->getVal();
	double processNormCorrFactor = modelEntry->second->normCorrFactor_;
	//std::cout << "processName = " << modelEntry->first << ": processNormCorrFactor = " << processNormCorrFactor << std::endl;

	double processBinContent_scaled = processNormCorrFactor*processNorm*processBinContent;
	double processBinError_scaled = processNormCorrFactor*processNorm*processBinError;

	fitBinContent += processBinContent_scaled;
	fitBinError2 += processBinError_scaled*processBinError_scaled;
      }
      
      //std::cout << "iBin = " << iBin << ": dataBinContent = " << dataBinContent << ", fitBinContent = " << fitBinContent << std::endl;

      double diffBinContent2 = (dataBinContent - fitBinContent)*(dataBinContent - fitBinContent);
      double diffBinError2 = fitBinError2 + dataBinError*dataBinError;
      
      if ( diffBinError2 > 0. ) {
	chi2 += (diffBinContent2/diffBinError2);
	++numDoF;
      }
    }
  }

//--- correct number of degrees of freedom
//    for number of fitted parameters
  numDoF -= modelEntries_.size();

  //std::cout << "chi2 = " << chi2 << std::endl;
  //std::cout << "numDoF = " << numDoF << std::endl;
  
  if ( numDoF > 0 ) {
    return (chi2/numDoF);
  } else {
    edm::LogWarning ("compChi2red") << " numDoF = " << numDoF << " must not be negative"
				    << " returning Chi2red = 1.e+3 !!";
    return 1.e+3;
  }
}

void TemplateBgEstFit::estimateUncertainties(bool fluctStat, bool fluctSys, int numSamplings, 
					     double chi2redMax, const char* type, int printLevel, bool printWarnings)
{
  int numProcesses = modelEntries_.size();
  TVectorD fitValues(numProcesses);
  BgEstMean mean(numProcesses);
  BgEstMedian median(numProcesses);
  BgEstCovMatrix cov(numProcesses);

  unsigned numTotFits = 0;
  unsigned numGoodFits = 0;

  for ( int iRndFluct = 0; iRndFluct < numSamplings; ++iRndFluct ) {

    std::cout << "<TemplateBgEstFit::estimateUncertainties>: iRndFluct = " << iRndFluct << std::endl;

//--- fluctuate distributions observed in (pseudo)data
    dataEntry_->fluctuate(true, false);

//--- fluctuate template histograms fitted to the (pseudo)data  
//
//    CV: treat statistical uncertainties on template histograms
//        as systematic uncertainties of background estimation
//
    for ( modelEntryMap::iterator modelEntry = modelEntries_.begin();
	  modelEntry != modelEntries_.end(); ++modelEntry ) {
      //modelEntry->second->fluctuate(fluctStat, fluctSys);
      modelEntry->second->fluctuate(fluctSys, fluctSys);
    }

    buildFitData();

    buildFitModel();

    fit(false, printLevel, printWarnings);
    ++numTotFits;
    
    for ( int iProcess = 0; iProcess < numProcesses; ++iProcess ) {
      const std::string& processName = processNames_[iProcess];
      fitValues(iProcess) = modelEntries_[processName]->norm_->getVal();
      //std::cout << " fitValue(iProcess = " << iProcess << ", processName = " << processName << ")"
      //	    << " = " << fitValues(iProcess) << std::endl;
    }
    
    double chi2red = compChi2red();
    //std::cout << "Chi2red = " << chi2red << std::endl;
    if ( !(chi2red < chi2redMax) ) continue;
    
    mean.update(fitValues);
    median.update(fitValues);
    cov.update(fitValues);
    
    ++numGoodFits;
  }

  double badFitFraction = (numTotFits - numGoodFits)/((double)numTotFits);
  std::cout << "fraction of Samplings discarded due to bad Fit quality = " << badFitFraction << std::endl;

  std::cout << "Mean:" << std::endl;
  mean.print(std::cout, &processNames_);
  std::cout << "Median:" << std::endl;
  median.print(std::cout, &processNames_);
  std::cout << "Covariance Matrix:" << std::endl;
  cov.print(std::cout, &processNames_);
  
  if ( controlPlotsFileName_ != "" ) {
    //makeControlPlotsCovariance(fitResultMean_, median(), cov(), processNames_, controlPlotsFileName_, type);
    makeControlPlotsCovariance(fitResultMean_, fitResultMean_, cov(), processNames_, controlPlotsFileName_, type);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TemplateBgEstFit);
