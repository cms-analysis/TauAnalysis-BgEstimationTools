#include "TauAnalysis/BgEstimationTools/plugins/TemplateBgEstFit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"

#include "TauAnalysis/BgEstimationTools/interface/histogramAuxFunctions.h"
#include "TauAnalysis/BgEstimationTools/interface/BgEstMean.h"
#include "TauAnalysis/BgEstimationTools/interface/BgEstCovMatrix.h"

#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include <TEllipse.h>
#include <TColor.h>
#include <TMarker.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TVectorD.h>
#include <TMatrixD.h>

#include <RooArgSet.h>
#include <RooGlobalFunc.h>
#include <RooPlot.h>
#include <RooFit.h>

#include <iostream>

const int defaultCanvasSizeX = 800;
const int defaultCanvasSizeY = 600;

const float maxNorm = 1.e+6;

const std::string fluctDirection_up = "up";
const std::string fluctDirection_down = "down";
const std::string fluctDirection_bidirectional = "bidirectional";

enum { kUp, kDown, kBidirectional };

const std::string fluctMode_coherent = "coherent";
const std::string fluctMode_incoherent = "incoherent";

enum { kCoherent, kIncoherent };

TRandom3 gRndNum;

const double epsilon = 1.e-3;

const int fitStatus_converged = 0;

void drawErrorEllipse(double, double, double, double, double, const char*, const char*, const char*);
void makeCovariancePlots(TVectorD, TMatrixD, const std::vector<std::string>&, const std::string&, const char*);

double getSampledPull(int fluctDirection)
{
  assert(fluctDirection == kUp || fluctDirection == kDown || fluctDirection == kBidirectional);
  
  double fluctPull = 0.;
  bool fluctPull_isValid = false;

  while ( !fluctPull_isValid ) {
    double x = gRndNum.Gaus(0., 1.);
    if ( (fluctDirection == kUp   && x >= 0.) ||
	 (fluctDirection == kDown && x <= 0.) ||
	  fluctDirection == kBidirectional    ) {
      fluctPull = x;
      fluctPull_isValid = true;
    }
  }

  return fluctPull;
}

void sampleStatTH1(TH1* origHistogram, TH1* fluctHistogram)
{
//--- fluctuate bin-contents by Gaussian distribution
//    with zero mean and standard deviation given by bin-errors

  int numBins = origHistogram->GetNbinsX();
  for ( int iBin = 0; iBin <= numBins; ++iBin ) {
    double origBinContent = origHistogram->GetBinContent(iBin);
    double origBinError = origHistogram->GetBinError(iBin);

    double fluctPull = getSampledPull(kBidirectional);
    double fluctBinContent = origBinContent + fluctPull*origBinError;
    
    fluctHistogram->SetBinContent(iBin, fluctBinContent);
    fluctHistogram->SetBinError(iBin, origBinError);
  }
}

void sampleSysTH1(TH1* fluctHistogram, TH1* sysHistogram, int fluctDirection, int fluctMode)
{
  assert(fluctMode == kCoherent || fluctMode == kIncoherent);

//--- fluctuate bin-contents by Gaussian distribution
//    with zero mean and standard deviation given by bin-errors
  double sampledPull = getSampledPull(fluctDirection);

  int numBins = fluctHistogram->GetNbinsX();
  for ( int iBin = 0; iBin <= numBins; ++iBin ) {
    double fluctBinContent = fluctHistogram->GetBinContent(iBin);
    double fluctBinError = fluctHistogram->GetBinError(iBin);
    
    double sysBinContent = sysHistogram->GetBinContent(iBin);
    double sysBinError = sysHistogram->GetBinError(iBin);

    double modBinContent = fluctBinContent + sampledPull*sysBinContent;
    double modBinError = TMath::Sqrt(fluctBinError*fluctBinError + sysBinError*sysBinError);

    fluctHistogram->SetBinContent(iBin, modBinContent);
    fluctHistogram->SetBinError(iBin, modBinError);

    if ( fluctMode == kIncoherent ) sampledPull = getSampledPull(fluctDirection);
  }
}

void makePositiveTH1(TH1* fluctHistogram)
{
  int numBins = fluctHistogram->GetNbinsX();
  for ( int iBin = 0; iBin <= numBins; ++iBin ) {
    if ( fluctHistogram->GetBinContent(iBin) < 0. ) fluctHistogram->SetBinContent(iBin, 0.);
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::dataEntryType::dataEntryType(const std::string& processName, const std::string& meName)
  : histogram_(0), 
    cfgError_(0)
{
  //std::cout << "<dataEntryType::dataEntryType>:" << std::endl;
  //std::cout << " processName = " << processName << std::endl;
  //std::cout << " meName = " << meName << std::endl;

  processName_ = processName;

  meName_ = meName;
}

TemplateBgEstFit::dataEntryType::~dataEntryType()
{
  delete histogram_;
  delete fluctHistogram_;
}

void TemplateBgEstFit::dataEntryType::initialize(DQMStore& dqmStore, RooRealVar* x, int& error)
{
  //std::cout << "<dataEntryType::initialize>:" << std::endl;

  me_ = dqmStore.get(meName_);
  if ( !me_ ) {
    edm::LogError ("dataEntryType::initialize") << " Failed to access dqmMonitorElement = " << meName_ << " !!";
    error = 1;
    return;
  }

  x_ = x;

  histogramName_ = std::string(processName_).append("_histogram");
  histogram_ = new RooDataHist(histogramName_.data(), histogramName_.data(), *x_, me_->getTH1());

  fluctHistogram_ = (TH1*)me_->getTH1()->Clone();
}

void TemplateBgEstFit::dataEntryType::fluctuate(bool, bool)
{  
  sampleStatTH1(me_->getTH1(), fluctHistogram_);

  makePositiveTH1(fluctHistogram_);

  delete histogram_;
  histogram_ = new RooDataHist(histogramName_.data(), histogramName_.data(), *x_, fluctHistogram_);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::processEntryType::processEntryType(const std::string& processName, const std::string& meName,
						     int lineColor, int lineStyle, int lineWidth)
  : dataEntryType(processName, meName), 
    pdf_(0), norm_(0)
{
  lineColor_ = lineColor;
  lineStyle_ = lineStyle;
  lineWidth_ = lineWidth;
}

TemplateBgEstFit::processEntryType::~processEntryType()
{
  delete pdf_;
  delete norm_;
}

void TemplateBgEstFit::processEntryType::initialize(DQMStore& dqmStore, RooRealVar* x, int& error)
{
  //std::cout << "<processEntryType::initialize>:" << std::endl;

  dataEntryType::initialize(dqmStore, x, error);

  if ( error ) return;
 
  pdfName_ = std::string(processName_).append("_").append("pdf");
  pdf_ = new RooHistPdf(pdfName_.data(), pdfName_.data(), *x_, *histogram_);

  std::string normName = std::string(processName_).append("_").append("norm");
  norm_ = new RooRealVar(normName.data(), normName.data(), 0., maxNorm);

  for ( std::vector<sysFluctEntryType>::iterator sysErrFluctuation = sysErrFluctuations_.begin();
	sysErrFluctuation != sysErrFluctuations_.end(); ++sysErrFluctuation ) {
    //std::cout << " sysErrFluctuation->meName = " << sysErrFluctuation->meName_ << std::endl;
    sysErrFluctuation->me_ = dqmStore.get(sysErrFluctuation->meName_);
    if ( !sysErrFluctuation->me_ ) {
      edm::LogError ("processEntryType::initialize") << " Failed to access dqmMonitorElement = " << sysErrFluctuation->meName_ << " !!";
      error = 1;
      return;
    }

//--- check that histograms representing systematic uncertainties have the same binning
//    as that representing expectation
    if ( !isCompatibleBinning(me_->getTH1(), sysErrFluctuation->me_->getTH1()) ) {
      edm::LogError ("processEntryType::initialize") << " Incompatible binning of histograms " << meName_ 
						     << " and " << sysErrFluctuation->meName_ << " !!";
      error = 1;
      return;
    }
  }
}

void TemplateBgEstFit::processEntryType::fluctuate(bool fluctStat, bool fluctSys)
{
  if ( fluctStat ) {
    sampleStatTH1(me_->getTH1(), fluctHistogram_);
  } else {
    TH1* origHistogram = me_->getTH1();
    int numBins = origHistogram->GetNbinsX();
    for ( int iBin = 0; iBin <= numBins; ++iBin ) {
      double origBinContent = origHistogram->GetBinContent(iBin);
      double origBinError = origHistogram->GetBinError(iBin);

      fluctHistogram_->SetBinContent(iBin, origBinContent);
      fluctHistogram_->SetBinError(iBin, origBinError);
    }
  }

  if ( fluctSys ) {
    for ( std::vector<sysFluctEntryType>::const_iterator sysErrFluctuation = sysErrFluctuations_.begin();
	  sysErrFluctuation != sysErrFluctuations_.end(); ++sysErrFluctuation ) {
      sampleSysTH1(fluctHistogram_, sysErrFluctuation->me_->getTH1(),  sysErrFluctuation->direction_, sysErrFluctuation->mode_);
    }
  }

  makePositiveTH1(fluctHistogram_);

  delete histogram_;
  histogram_ = new RooDataHist(histogramName_.data(), histogramName_.data(), *x_, fluctHistogram_);
  delete pdf_;
  pdf_ = new RooHistPdf(pdfName_.data(), pdfName_.data(), *x_, *histogram_);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

TemplateBgEstFit::TemplateBgEstFit(const edm::ParameterSet& cfg)
  : cfgError_(0)
{
  std::cout << "<TemplateBgEstFit::TemplateBgEstFit>:" << std::endl;

  edm::ParameterSet cfgProcesses = cfg.getParameter<edm::ParameterSet>("processes");
  std::vector<std::string> processNames = cfgProcesses.getParameterNamesForType<edm::ParameterSet>();
  for ( std::vector<std::string>::const_iterator processName = processNames.begin(); 
	processName != processNames.end(); ++processName ) {
    edm::ParameterSet cfgProcess = cfgProcesses.getParameter<edm::ParameterSet>(*processName);

    std::string meName_process = cfgProcess.getParameter<std::string>("meName");

    edm::ParameterSet cfgDrawOptions_process = cfgProcess.getParameter<edm::ParameterSet>("drawOptions");
    int lineColor_process = cfgDrawOptions_process.getParameter<int>("lineColor");
    int lineStyle_process = cfgDrawOptions_process.getParameter<int>("lineStyle");
    int lineWidth_process = cfgDrawOptions_process.getParameter<int>("lineWidth");

    processNames_.push_back(*processName);

    processEntryType* processEntry = new processEntryType(*processName, meName_process, 
							  lineColor_process, lineStyle_process, lineWidth_process);
    processEntries_.push_back(processEntry);
  }

  edm::ParameterSet cfgData = cfg.getParameter<edm::ParameterSet>("data");
  std::string meName_data = cfgData.getParameter<std::string>("meName");
  dataEntry_ = new dataEntryType("data", meName_data);
  
  edm::ParameterSet cfgFit = cfg.getParameter<edm::ParameterSet>("fit");
  variableName_ = ( cfgFit.exists("variableName") ) ? cfgFit.getParameter<std::string>("variableName") : "";
  variableTitle_ = ( cfgFit.exists("variableTitle") ) ? cfgFit.getParameter<std::string>("variableTitle") : "";
  xMin_ = cfgFit.getParameter<double>("xMin");
  xMax_ = cfgFit.getParameter<double>("xMax");
  printLevel_ = ( cfgFit.exists("printLevel") ) ? cfgFit.getParameter<int>("printLevel") : 1;
  printWarnings_  = ( cfgFit.exists("printWarnings") ) ? cfgFit.getParameter<bool>("printWarnings") : true;

  edm::ParameterSet cfgControlPlots = cfg.getParameter<edm::ParameterSet>("output").getParameter<edm::ParameterSet>("controlPlots");
  controlPlotsFileName_ = cfgControlPlots.getParameter<std::string>("fileName");

  edm::ParameterSet cfgStatErr = cfg.getParameter<edm::ParameterSet>("estStatUncertainties");
  statErrNumSamplings_ = cfgStatErr.getParameter<edm::ParameterSet>("numSamplings").getParameter<int>("stat");
  statErrChi2redMax_ = cfgStatErr.getParameter<double>("chi2redMax");  
  statErrPrintLevel_ = cfgStatErr.getParameter<edm::ParameterSet>("verbosity").getParameter<int>("printLevel");
  statErrPrintWarnings_ = cfgStatErr.getParameter<edm::ParameterSet>("verbosity").getParameter<bool>("printWarnings");

  edm::ParameterSet cfgSysErr = cfg.getParameter<edm::ParameterSet>("estSysUncertainties");
  edm::ParameterSet cfgSysFluctuations = cfgSysErr.getParameter<edm::ParameterSet>("fluctuations");
  std::vector<std::string> sysFluctNames = cfgSysFluctuations.getParameterNamesForType<edm::ParameterSet>();
  for ( std::vector<std::string>::const_iterator sysFluctName = sysFluctNames.begin(); 
	sysFluctName != sysFluctNames.end(); ++sysFluctName ) {
    edm::ParameterSet cfgSysFluct = cfgSysFluctuations.getParameter<edm::ParameterSet>(*sysFluctName);

    std::string fluctDirection_string = cfgSysFluct.getParameter<std::string>("direction");
    int fluctDirection_int = -1;
    if ( fluctDirection_string == fluctDirection_up ) {
      fluctDirection_int = kUp;
    } else if ( fluctDirection_string == fluctDirection_down ) {
      fluctDirection_int = kDown;
    } else if ( fluctDirection_string == fluctDirection_bidirectional ) {
      fluctDirection_int = kBidirectional;
    } else {
      edm::LogError ("TemplateBgEstFit::TemplateBgEstFit") << " Invalid 'direction' parameter = " << fluctDirection_string << " !!";
      cfgError_ = 1;
    }

    std::string fluctMode_string = cfgSysFluct.getParameter<std::string>("mode");
    int fluctMode_int = -1;
    if ( fluctMode_string == fluctMode_coherent ) {
      fluctMode_int = kCoherent;
    } else if ( fluctMode_string == fluctMode_incoherent ) {
      fluctMode_int = kIncoherent;
    } else {
      edm::LogError ("TemplateBgEstFit::TemplateBgEstFit") << " Invalid 'mode' parameter = " << fluctMode_string << " !!";
      cfgError_ = 1;
    }

    edm::ParameterSet cfgProcesses = cfgSysFluct.getParameter<edm::ParameterSet>("meNames");

    for ( std::vector<processEntryType*>::iterator processEntry = processEntries_.begin();
	  processEntry != processEntries_.end(); ++processEntry ) {
      const std::string& processName = (*processEntry)->processName_;
      if ( !cfgProcesses.exists(processName) ) {
	edm::LogError ("TemplateBgEstFit::TemplateBgEstFit") << " No MonitorElement defined for estimating systematic Uncertainty = " 
							     << (*sysFluctName) << " of process = " << processName << " !!";
	cfgError_ = 1;
	continue;
      }

      std::string meName = cfgProcesses.getParameter<std::string>(processName);

      sysFluctEntryType sysFluctEntry;
      sysFluctEntry.fluctName_ = (*sysFluctName);
      sysFluctEntry.meName_ = meName;
      sysFluctEntry.direction_ = fluctDirection_int;
      sysFluctEntry.mode_ = fluctMode_int;

      (*processEntry)->sysErrFluctuations_.push_back(sysFluctEntry);
    }
  }
  sysErrNumStatSamplings_ = cfgSysErr.getParameter<edm::ParameterSet>("numSamplings").getParameter<int>("stat");
  sysErrNumSysSamplings_ = cfgSysErr.getParameter<edm::ParameterSet>("numSamplings").getParameter<int>("sys");
  sysErrChi2redMax_ = cfgSysErr.getParameter<double>("chi2redMax");  
  sysErrPrintLevel_ = cfgSysErr.getParameter<edm::ParameterSet>("verbosity").getParameter<int>("printLevel");
  sysErrPrintWarnings_ = cfgSysErr.getParameter<edm::ParameterSet>("verbosity").getParameter<bool>("printWarnings");

  x_ = 0;
  model_ = 0;
  fitResult_ = 0;
}

TemplateBgEstFit::~TemplateBgEstFit()
{
  //std::cout << "<TemplateBgEstFit::~TemplateBgEstFit>:" << std::endl;

  for ( std::vector<processEntryType*>::iterator it = processEntries_.begin();
	it != processEntries_.end(); ++it ) {
    delete (*it);
  }

  delete dataEntry_;

  delete x_;
  delete model_;
  delete fitResult_;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void TemplateBgEstFit::buildModel()
{
  std::string modelName = std::string("pdfModel");
  std::string modelTitle = modelName;

  TObjArray model_pdfCollection;
  TObjArray model_normCollection;
  for ( std::vector<processEntryType*>::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    model_pdfCollection.Add((*processEntry)->pdf_);
    model_normCollection.Add((*processEntry)->norm_);
  }

  std::string model_pdfArgName = std::string("pdfModel").append("_pdfArgs");
  RooArgList model_pdfArgs(model_pdfCollection, model_pdfArgName.data());
  std::string model_normArgName = std::string("pdfModel").append("_normArgs");
  RooArgList model_normArgs(model_normCollection, model_normArgName.data());

  //std::cout << "--> creating RooAddPdf with name = " << modelName << std::endl;
  model_ = new RooAddPdf(modelName.data(), modelTitle.data(), model_pdfArgs, model_normArgs);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void TemplateBgEstFit::endJob()
{
  std::cout << "<TemplateBgEstFit::endJob>:" << std::endl;

//--- check that configuration parameters contain no errors
  for ( std::vector<processEntryType*>::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    if ( (*processEntry)->cfgError_ ) cfgError_ = 1;
  }

  if ( dataEntry_->cfgError_ ) cfgError_ = 1;
  
  if ( cfgError_ ) {
    edm::LogError ("TemplateBgEstFit::endJob") << " Error in Configuration ParameterSet" 
					       << " --> histogram will NOT be fitted !!";
    return;
  }

//--- check that DQMStore service is available
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("TemplateBgEstFit::endJob") << " Failed to access dqmStore" 
					       << " --> histogram will NOT be fitted !!";
    return;
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());

  //dqmStore.showDirStructure();
  
//--- configure RooFit structure
  x_ = new RooRealVar("x", "x", xMin_, xMax_);

  int error = 0;
  for ( std::vector<processEntryType*>::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    (*processEntry)->initialize(dqmStore, x_, error);
  }

  dataEntry_->initialize(dqmStore, x_, error);
  
  if ( error ) {
    edm::LogError ("TemplateBgEstFit::endJob") << " Failed to access dqmMonitorElement(s) --> histograms will NOT be fitted !!";
    return;
  }
  
//--- fit template shapes of different signal and background processes
//    to distribution observed in (pseudo)data
  buildModel();

  std::cout << ">>> RootFit model used for generalized Matrix method Fit <<<" << std::endl;
  model_->printCompactTree();

  std::cout << ">>> RootFit Variables <<<" << std::endl;
  model_->getVariables()->Print("v");

  std::cout << ">>> RootFit Parameters <<<" << std::endl;
  model_->getParameters(dataEntry_->histogram_)->Print("v");

  std::cout << ">>> RootFit Observables <<<" << std::endl;
  model_->getObservables(dataEntry_->histogram_)->Print("v");

//--- NOTE: RooFit::Warnings not implemented in RooFit version 
//          included in ROOT 5.18/00a linked against CMSSW_2_2_13
  fitResult_ = model_->fitTo(*dataEntry_->histogram_, RooFit::Extended(), RooFit::Save(true), 
			     RooFit::PrintLevel(printLevel_) /* , RooFit::Warnings(printWarnings_) */ );
  assert(fitResult_ != NULL);

//--- print-out fit results
  std::cout << ">>> Fit Results <<<" << std::endl;
  std::cout << " fitStatus = " << fitResult_->status() << std::endl;
  std::cout << " chi2red = " << compChi2red() << std::endl;
  print(std::cout);

//--- produce plot of different signal and background processes
//    using scale factors determined by fit
//    compared to distribution of (pseudo)data
  if ( controlPlotsFileName_ != "" ) makeControlPlots();

//--- estimate statistical uncertainties
  std::cout << ">>> Statistical Uncertainties <<<" << std::endl;
  estimateUncertainties(true, statErrNumSamplings_, false, 1, statErrChi2redMax_, 
			"estStatUncertainties", statErrPrintLevel_, statErrPrintWarnings_);

//--- estimate systematic uncertainties
  std::cout << ">>> Systematic Uncertainties <<<" << std::endl;
  estimateUncertainties(false, sysErrNumStatSamplings_, true, sysErrNumSysSamplings_, sysErrChi2redMax_,
			"estSysUncertainties", sysErrPrintLevel_, sysErrPrintWarnings_);

//--- estimate total (statistical + systematic) uncertainties
  std::cout << ">>> Total (statistical + systematic) Uncertainties <<<" << std::endl;
  double chi2redMax = TMath::Max(statErrChi2redMax_, sysErrChi2redMax_);
  int totErrPrintLevel = TMath::Min(statErrPrintLevel_, sysErrPrintLevel_);
  bool totErrPrintWarnings = (statErrPrintWarnings_ && sysErrPrintWarnings_);
  estimateUncertainties(true, sysErrNumStatSamplings_, true, sysErrNumSysSamplings_, chi2redMax,
			"estTotUncertainties", totErrPrintLevel, totErrPrintWarnings);

  std::cout << "done." << std::endl;
}

void TemplateBgEstFit::print(std::ostream& stream)
{
  stream << "Fit Parameter:" << std::endl;
  for ( std::vector<processEntryType*>::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    stream << " " << (*processEntry)->processName_ << ": normalization = " << (*processEntry)->norm_->getVal();
    if ( (*processEntry)->norm_->hasAsymError() ) {
      stream << " + " << (*processEntry)->norm_->getAsymErrorHi() << " - " << fabs((*processEntry)->norm_->getAsymErrorLo());
    } else if ( (*processEntry)->norm_->hasError() ) {
      stream << " +/- " << (*processEntry)->norm_->getError();
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

  TCanvas canvas("TemplateBgEstFit", "TemplateBgEstFit", defaultCanvasSizeX, defaultCanvasSizeY);
  canvas.SetFillColor(10);

  RooPlot* plotFrame = x_->frame();
  std::string plotTitle = std::string("TemplateBgEstFit - Results for fitting ").append(variableName_);
  plotFrame->SetTitle(plotTitle.data());
  dataEntry_->histogram_->plotOn(plotFrame, RooFit::MarkerColor(kBlack), RooFit::MarkerStyle(2));
  for ( std::vector<processEntryType*>::const_iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    std::string componentName = (*processEntry)->pdf_->GetName();
    model_->plotOn(plotFrame, RooFit::Components(componentName.data()), 
		   RooFit::LineColor((*processEntry)->lineColor_), 
		   RooFit::LineStyle((*processEntry)->lineStyle_), 
		   RooFit::LineWidth((*processEntry)->lineWidth_));
  }
  model_->plotOn(plotFrame, RooFit::LineColor(kBlack), RooFit::LineStyle(kSolid), RooFit::LineWidth(2));
  plotFrame->SetXTitle(variableTitle_.data());
  //plotFrame->SetTitleOffset(1.2, "X");
  plotFrame->SetYTitle("");
  //plotFrame->SetTitleOffset(1.2, "Y");
  plotFrame->Draw();
  
  canvas.Update();

  int errorFlag = 0;
  std::string fileName = replace_string(controlPlotsFileName_, plotKeyword, std::string("normalization"), 1, 1, errorFlag);
  if ( !errorFlag ) {
    canvas.Print(fileName.data());
  } else {
    edm::LogError("TemplateBgEstFit::makeControlPlots") << " Failed to decode controlPlotsFileName = " << controlPlotsFileName_ 
							<< " --> skipping !!";
    return;
  }

//--- produce control plots of one and two sigma error contours 
//    showing correlation of estimated normalization parameters 
  const RooArgList& fitParameter = fitResult_->floatParsFinal();
  int numFitParameter = fitParameter.getSize();
  TVectorD mean(numFitParameter);
  TMatrixD cov(numFitParameter, numFitParameter);
  std::vector<std::string> labels(numFitParameter);
  for ( int iX = 0; iX < numFitParameter; ++iX ) {
    const RooAbsArg* paramX_arg = fitParameter.at(iX);
    const RooRealVar* paramX = dynamic_cast<const RooRealVar*>(paramX_arg);
    assert(paramX != NULL);
    mean(iX) = paramX->getVal();
    double sigmaX = paramX->getError();    
    for ( int iY = 0; iY < numFitParameter; ++iY ) {
      if ( iY == iX ) {
	cov(iX, iX) = sigmaX*sigmaX;
      } else {
	const RooAbsArg* paramY_arg = fitParameter.at(iY);
	const RooRealVar* paramY = dynamic_cast<const RooRealVar*>(paramY_arg);
	assert(paramY != NULL);
	double sigmaY = paramY->getError();
	double corrXY = fitResult_->correlation(*paramX_arg, *paramY_arg);
	cov(iX, iY) = sigmaX*sigmaY*corrXY;
      }
    }
    labels[iX] = paramX_arg->GetName();
  }

  makeCovariancePlots(mean, cov, labels, controlPlotsFileName_, "");
}

double TemplateBgEstFit::compChi2red() const
{
  //std::cout << "<TemplateBgEstFit::compChi2red>:" << std::endl;
/*
//--- check that template histograms for all processes
//    are properly normalized to unit area
  for ( std::vector<processEntryType*>::const_iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    TH1* histogramProcess = (*processEntry)->me_->getTH1();

    std::cout << "processName = " << (*processEntry)->processName_ << ": integral = " << histogramProcess->Integral() << std::endl;
  }
 */
  double chi2 = 0.;
  int numDoF = 0;
  TH1* histogramData = dataEntry_->me_->getTH1();
  int numBins = histogramData->GetNbinsX();
  for ( int iBin = 1; iBin < numBins; ++iBin ) {
    double dataBinCenter = histogramData->GetBinCenter(iBin);

//--- restrict computation of chi^2 to range included in fit
    if ( !(dataBinCenter > xMin_ && dataBinCenter < xMax_) ) continue;

    //std::cout << "iBin = " << iBin << ":" << std::endl;

    double dataBinContent = histogramData->GetBinContent(iBin);
    //std::cout << " dataBinContent = " << dataBinContent << std::endl;
    double dataBinError = histogramData->GetBinError(iBin);
    //std::cout << " dataBinError = " << dataBinError << std::endl;

    double fitBinContent = 0.;
    double fitBinError2 = 0.;
    for ( std::vector<processEntryType*>::const_iterator processEntry = processEntries_.begin();
	    processEntry != processEntries_.end(); ++processEntry ) {
      TH1* histogramProcess = (*processEntry)->me_->getTH1();

      double processBinContent = histogramProcess->GetBinContent(iBin);
      double processBinError = histogramProcess->GetBinError(iBin);

      //std::cout << " processName = " << (*processEntry)->processName_ << std::endl;

      double processNorm = (*processEntry)->norm_->getVal();
      //std::cout << "  processNorm = " << processNorm << std::endl;

      double processBinContent_scaled = processNorm*processBinContent;
      //std::cout << "  processBinContent_scaled = " << processBinContent_scaled << std::endl;
      double processBinError_scaled = processNorm*processBinError;
      //std::cout << "  processBinError_scaled = " << processBinError_scaled << std::endl;

      fitBinContent += processBinContent_scaled;
      fitBinError2 += processBinError_scaled*processBinError_scaled;
    }

    double diffBinContent2 = (dataBinContent - fitBinContent)*(dataBinContent - fitBinContent);
    //std::cout << " diffBinContent2 = " << diffBinContent2 << std::endl;
    double diffBinError2 = fitBinError2 + dataBinError*dataBinError;
    //std::cout << " diffBinError2 = " << diffBinError2 << std::endl;

    if ( diffBinError2 > 0. ) {
      chi2 += (diffBinContent2/diffBinError2);
      ++numDoF;
    }
  }

//--- correct number of degrees of freedom
//    for number of fitted parameters
  numDoF -= processEntries_.size();

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

void TemplateBgEstFit::estimateUncertainties(bool fluctStat, int numStatSamplings, bool fluctSys, int numSysSamplings, 
					     double chi2redMax, const char* type, int printLevel, bool printWarnings)
{
  int numProcesses = processEntries_.size();
  TVectorD fitValues(numProcesses);
  BgEstMean mean(numProcesses);
  BgEstCovMatrix cov(numProcesses);

  unsigned numTotFits = 0;
  unsigned numGoodFits = 0;

  for ( int iRndStat = 0; iRndStat < numStatSamplings; ++iRndStat ) {
    for ( int iRndSys = 0; iRndSys < numSysSamplings; ++iRndSys ) {

      //std::cout << "iRndStat = " << iRndStat << ", iRndSys = " << iRndSys << ":" << std::endl;

//--- fluctuate distribution observed in (pseudo)data
      dataEntry_->fluctuate(true, false);

//--- fluctuate template histograms fitted to the (pseudo)data  
      for ( std::vector<processEntryType*>::iterator processEntry = processEntries_.begin();
	    processEntry != processEntries_.end(); ++processEntry ) {
	(*processEntry)->fluctuate(fluctStat, fluctSys);
      }

      delete model_;
      buildModel();

      delete fitResult_;
      fitResult_ = model_->fitTo(*dataEntry_->histogram_, RooFit::Extended(), RooFit::Save(true), 
				 RooFit::PrintLevel(printLevel) /* , RooFit::Warnings(printWarnings) */ );

      ++numTotFits;

      for ( int iProcess = 0; iProcess < numProcesses; ++iProcess ) {
	fitValues(iProcess) = processEntries_[iProcess]->norm_->getVal();
	//std::cout << " fitValue(iProcess = " << iProcess << ", processName = " << processEntries_[iProcess]->processName_ << ")"
	//	    << " = " << fitValues(iProcess) << std::endl;
      }

      int fitStatus = fitResult_->status();
      //std::cout << " fitStatus = " << fitStatus << std::endl;

      double chi2red = compChi2red();
      //std::cout << " chi2red = " << chi2red << std::endl;

      if ( !(fitStatus == fitStatus_converged && chi2red < chi2redMax) ) continue;

      mean.update(fitValues);
      cov.update(fitValues);

      ++numGoodFits;
    }
  }

  double badFitFraction = (numTotFits - numGoodFits)/((double)numTotFits);
  std::cout << "fraction of Samplings discarded due to bad Fit quality = " << badFitFraction << std::endl;

  std::cout << "Mean:" << std::endl;
  mean.print(std::cout, &processNames_);
  std::cout << "Covariance Matrix:" << std::endl;
  cov.print(std::cout, &processNames_);
  
  if ( controlPlotsFileName_ != "" ) makeCovariancePlots(mean(), cov(), processNames_, controlPlotsFileName_, type);
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void drawErrorEllipse(double x0, double y0, double Sxx, double Sxy, double Syy, 
		      const char* labelX, const char* labelY, const char* fileName)
{
  //std::cout << "<drawErrorEllipse>:" << std::endl;
  //std::cout << " x0 = " << x0 << std::endl;
  //std::cout << " y0 = " << y0 << std::endl;
  //if ( Sxy >  0. ) std::cout << "variables are correlated" << std::endl;
  //if ( Sxy == 0. ) std::cout << "variables are uncorrelated" << std::endl;
  //if ( Sxy <  0. ) std::cout << "variables are anti-correlated" << std::endl;
  //std::cout << " -2*Sxy = " << -2*Sxy << std::endl;
  //std::cout << " Sxx - Syy = " << Sxx - Syy << std::endl;
  //std::cout << " labelX = " << labelX << std::endl;
  //std::cout << " labelY = " << labelY << std::endl;
  //std::cout << " fileName = " << fileName << std::endl;

//--- draw one and two sigma error contours 
//    centered at fit results (x0,y0) and with (correlated) uncertainties 
//    estimated by elements Sxx, Sxy, Syy of covariance matrix passed as function arguments
//    (note that since the covariance matrix is symmetric, 
//     there is no need to pass element Syx of the covariance matrix)

  TCanvas canvas("drawErrorEllipse", "drawErrorEllipse", 600, 600);
  canvas.SetFillColor(10);

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

  TEllipse oneSigmaErrorEllipse(x0, y0, sigmaX_transformed*1., sigmaY_transformed*1., 0., 360., alpha*180./TMath::Pi()); 
  oneSigmaErrorEllipse.SetFillColor(5);
  oneSigmaErrorEllipse.SetLineColor(44);
  oneSigmaErrorEllipse.SetLineWidth(1);
  TEllipse twoSigmaErrorEllipse(x0, y0, sigmaX_transformed*2., sigmaY_transformed*2., 0., 360., alpha*180./TMath::Pi()); 
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
  double minX = x0 - 2.2*TMath::Abs(cosAlpha*sigmaX_transformed + sinAlpha*sigmaY_transformed);
  double maxX = x0 + 2.8*TMath::Abs(cosAlpha*sigmaX_transformed + sinAlpha*sigmaY_transformed);
  double minY = y0 - 2.2*TMath::Abs(sinAlpha*sigmaX_transformed + cosAlpha*sigmaY_transformed);
  double maxY = y0 + 2.8*TMath::Abs(sinAlpha*sigmaX_transformed + cosAlpha*sigmaY_transformed);

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

void makeCovariancePlots(TVectorD mean, TMatrixD cov, const std::vector<std::string>& labels, 
			 const std::string& controlPlotsFileName, const char* type)
{
  int numFitParameter = mean.GetNoElements();
  assert(cov.GetNrows() == numFitParameter && cov.GetNcols() == numFitParameter);
  for ( int iX = 0; iX < numFitParameter; ++iX ) {
    double x0 = mean(iX);
    double Sxx = cov(iX, iX);
    const char* labelX = labels[iX].data();

    for ( int iY = 0; iY < iX; ++iY ) {
      double y0 = mean(iY);
      double Syy = cov(iY, iY);
      const char* labelY = labels[iY].data();

      double Sxy = cov(iX, iY);
      std::string fileNameParam = std::string("corr_").append(labelX).append("_vs_").append(labelY);
      if ( type != "" ) fileNameParam.append("_").append(type);
      
      int errorFlag = 0;
      std::string fileName = replace_string(controlPlotsFileName, plotKeyword, fileNameParam, 1, 1, errorFlag);
      if ( !errorFlag ) {
	drawErrorEllipse(x0, y0, Sxx, Sxy, Syy, labelX, labelY, fileName.data());
      } else {
	edm::LogError("drawErrorEllipses") << " Failed to decode controlPlotsFileName = " << controlPlotsFileName 
					   << " --> skipping !!";
	return;
      }
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TemplateBgEstFit);
