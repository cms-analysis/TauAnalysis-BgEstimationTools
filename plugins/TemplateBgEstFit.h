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
 * \version $Revision: 1.3 $
 *
 * $Id: TemplateBgEstFit.h,v 1.3 2009/08/04 14:36:59 veelken Exp $
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
    virtual void fluctuate(bool, bool);
    std::string processName_;
    std::string meName_;
    MonitorElement* me_;
    RooRealVar* x_;
    std::string histogramName_;
    RooDataHist* histogram_;
    TH1* fluctHistogram_;
    int cfgError_;
  };

  struct sysFluctEntryType
  {
    std::string fluctName_;
    std::string meName_;
    MonitorElement* me_;
    int direction_;
    int mode_;
  };

  struct processEntryType : public dataEntryType
  {
    processEntryType(const std::string&, const std::string&, int, int, int);
    virtual ~processEntryType();
    void initialize(DQMStore&, RooRealVar*, int&);
    void fluctuate(bool, bool);
    std::string pdfName_;
    RooHistPdf* pdf_;
    RooRealVar* norm_;
    int lineColor_;
    int lineStyle_;
    int lineWidth_;
    std::vector<sysFluctEntryType> sysErrFluctuations_;
  };

 public:
  
  explicit TemplateBgEstFit(const edm::ParameterSet&);
  ~TemplateBgEstFit();
  
 private:

  void beginJob(const edm::EventSetup&) {}
  void analyze(const edm::Event&, const edm::EventSetup&) {}
  void endJob();

//--- private auxiliary functions
  void buildModel();
  void print(std::ostream& stream);
  void makeControlPlots();
  double compChi2red() const;
  void estimateUncertainties(bool, int, bool, int, double, const char*, int, bool);

//--- configuration parameters
  std::string variableName_;
  std::string variableTitle_;
  double xMin_;
  double xMax_;
  int printLevel_;
  bool printWarnings_;

  std::string controlPlotsFileName_;

  int statErrNumSamplings_;
  double statErrChi2redMax_;
  int statErrPrintLevel_;
  bool statErrPrintWarnings_;

  int sysErrNumStatSamplings_;
  int sysErrNumSysSamplings_;
  double sysErrChi2redMax_;
  int sysErrPrintLevel_;
  bool sysErrPrintWarnings_;
  
//--- internal data-members for handling histograms
//    and performing fit
  std::vector<std::string> processNames_;
  std::vector<processEntryType*> processEntries_;
  
  dataEntryType* dataEntry_;

  RooRealVar* x_;

  RooAddPdf* model_;
  RooFitResult* fitResult_;

  int cfgError_;
};

#endif  


