#ifndef TauAnalysis_BgEstimationTools_TemplateBgEstFit_h  
#define TauAnalysis_BgEstimationTools_TemplateBgEstFit_h

/** \class TemplateBgEstFit
 *
 * Estimate contribution of signal and background processes
 * to final event sample by fitting shape "templates" for different processes
 * to distribution observed in data
 * (class implements "template" method for data-driven background estimation)
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: TemplateBgEstFit.h,v 1.2 2009/07/15 09:03:47 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooRealVar.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>

#include <string>

class TemplateBgEstFit : public edm::EDAnalyzer
{
  struct dataEntryType
  {
    dataEntryType(const std::string&, const std::string&);
    virtual ~dataEntryType();
    virtual void initialize(DQMStore&, RooRealVar*, int&);
    std::string processName_;
    std::string dqmDirectory_store_;
    std::string meName_;
    RooDataHist* histogram_;
    int cfgError_;
  };

  struct processEntryType : public dataEntryType
  {
    processEntryType(const std::string&, const std::string&, int, int, int);
    virtual ~processEntryType();
    virtual void initialize(DQMStore&, RooRealVar*, int&);
    RooHistPdf* pdf_;
    RooRealVar* norm_;
    int lineColor_;
    int lineStyle_;
    int lineWidth_;
  };

 public:
  
  explicit TemplateBgEstFit(const edm::ParameterSet&);
  ~TemplateBgEstFit();
  
 private:

  void beginJob(const edm::EventSetup&) {}
  void analyze(const edm::Event&, const edm::EventSetup&) {}
  void endJob();

//--- private auxiliary functions
  void print(std::ostream& stream);
  void makeControlPlots();

//--- configuration parameters
  std::string variableName_;
  std::string variableTitle_;
  double xMin_;
  double xMax_;

  std::string controlPlotsFileName_;

//--- internal data-members for handling histograms
//    and performing fit
  std::vector<processEntryType*> processEntries_;
  
  dataEntryType* dataEntry_;

  RooRealVar* x_;

  RooAddPdf* model_;
  RooFitResult* fitResult_;

  int cfgError_;
};

#endif  


