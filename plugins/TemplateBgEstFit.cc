#include "TauAnalysis/BgEstimationTools/plugins/TemplateBgEstFit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"

#include <TCanvas.h>
#include <TROOT.h>

#include <RooArgSet.h>
#include <RooGlobalFunc.h>
#include <RooPlot.h>
#include <RooFit.h>

#include <iostream>

const int defaultCanvasSizeX = 800;
const int defaultCanvasSizeY = 600;

const float maxNorm = 1.e+6;

TemplateBgEstFit::dataEntryType::dataEntryType(const std::string& processName, const std::string& meName)
  : histogram_(0),
    cfgError_(0)
{
  processName_ = processName;

  if ( meName.find("/") != std::string::npos ) {
    size_t posSeparator = meName.find_last_of("/");

    if ( posSeparator < (meName.length() - 1) ) {
      dqmDirectory_store_ = std::string(meName, 0, posSeparator);
      meName_ = std::string(meName, posSeparator + 1);
    } else {
      edm::LogError ("processEntryType") << " Invalid Configuration Parameter 'meName' = " << meName << " !!";
      cfgError_ = 1;
    }
  } else {
    dqmDirectory_store_ = "";
    meName_ = meName;
  }
}

TemplateBgEstFit::dataEntryType::~dataEntryType()
{
  delete histogram_;
}

void TemplateBgEstFit::dataEntryType::initialize(DQMStore& dqmStore, RooRealVar* x, int& error)
{
  std::cout << "<dataEntryType::initialize>:" << std::endl;

  std::string meName_full = dqmDirectoryName(dqmDirectory_store_).append(meName_);
  std::cout << " meName_full = " << meName_full << std::endl;
  
  MonitorElement* me = dqmStore.get(meName_full);
  if ( !me ) {
    edm::LogError ("dataEntryType::initialize") << " Failed to access dqmMonitorElement = " << meName_full << " !!";
    error = 1;
    return;
  }

  std::string histogramName = std::string(processName_).append("_histogram");
  histogram_ = new RooDataHist(histogramName.data(), histogramName.data(), *x, me->getTH1());
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
  std::cout << "<processEntryType::initialize>:" << std::endl;

  dataEntryType::initialize(dqmStore, x, error);

  if ( error ) return;
 
  std::string pdfName = std::string(processName_).append("_").append("pdf");
  pdf_ = new RooHistPdf(pdfName.data(), pdfName.data(), *x, *histogram_);

  std::string normName = std::string(processName_).append("_").append("norm");
  norm_ = new RooRealVar(normName.data(), normName.data(), 0., maxNorm);
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

  edm::ParameterSet cfgControlPlots = cfg.getParameter<edm::ParameterSet>("output").getParameter<edm::ParameterSet>("controlPlots");
  controlPlotsFileName_ = cfgControlPlots.getParameter<std::string>("fileName");
}

TemplateBgEstFit::~TemplateBgEstFit()
{
  std::cout << "<TemplateBgEstFit::~TemplateBgEstFit>:" << std::endl;

  for ( std::vector<processEntryType*>::iterator it = processEntries_.begin();
	it != processEntries_.end(); ++it ) {
    delete (*it);
  }

  delete dataEntry_;

  delete x_;

  delete model_;
}

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
    edm::LogError ("endJob") << " Error in Configuration ParameterSet" 
			     << " --> histogram will NOT be fitted !!";
    return;
  }

//--- check that DQMStore service is available
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("endJob") << " Failed to access dqmStore" 
			     << " --> histogram will NOT be fitted !!";
    return;
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());

  dqmStore.showDirStructure();
  
//--- configure RooFit structure
  x_ = new RooRealVar("x", "x", xMin_, xMax_);

  int error = 0;
  for ( std::vector<processEntryType*>::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    (*processEntry)->initialize(dqmStore, x_, error);
  }

  dataEntry_->initialize(dqmStore, x_, error);
  
  if ( error ) {
    edm::LogError ("beginJob") << " Failed to access dqmMonitorElement(s) --> histograms will NOT be fitted !!";
    return;
  }
  
//--- fit template shapes of different signal and background processes
//    to distribution observed in (pseudo)data
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
  std::cout << "--> creating RooAddPdf with name = " << modelName << std::endl;
  model_ = new RooAddPdf(modelName.data(), modelTitle.data(), model_pdfArgs, model_normArgs);

  std::cout << ">>> RootFit model used for generalized Matrix method Fit <<<" << std::endl;
  model_->printCompactTree();

  std::cout << ">>> RootFit Variables <<<" << std::endl;
  model_->getVariables()->Print("v");

  std::cout << ">>> RootFit Parameters <<<" << std::endl;
  model_->getParameters(dataEntry_->histogram_)->Print("v");

  std::cout << ">>> RootFit Observables <<<" << std::endl;
  model_->getObservables(dataEntry_->histogram_)->Print("v");

  model_->fitTo(*dataEntry_->histogram_, RooFit::Extended());

//--- print-out fit results
  print(std::cout);

//--- produce plot of different signal and background processes
//    using scale factors determined by fit
//    compared to distribution of (pseudo)data
  if ( controlPlotsFileName_ != "" ) makeControlPlots();
}

void TemplateBgEstFit::print(std::ostream& stream)
{
  stream << "<TemplateBgEstFit::print>:" << std::endl;

  stream << "Fit Parameter:" << std::endl;
  for ( std::vector<processEntryType*>::iterator processEntry = processEntries_.begin();
	processEntry != processEntries_.end(); ++processEntry ) {
    stream << (*processEntry)->processName_ << ": normalization = " << (*processEntry)->norm_->getVal() << std::endl;
  }
}

void TemplateBgEstFit::makeControlPlots()
{
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
  
  canvas.Print(controlPlotsFileName_.data());
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(TemplateBgEstFit);
