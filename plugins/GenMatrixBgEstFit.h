#ifndef TauAnalysis_BgEstimationTools_GenMatrixBgEstFit_h  
#define TauAnalysis_BgEstimationTools_GenMatrixBgEstFit_h

/** \class GenMatrixBgEstFit
 *
 * Estimate contribution of signal and background processes
 * to final event sample by fitting number of events 
 * observed in different control regions in data
 * to signal and background models,
 * based on the assumption that observables are uncorrelated
 * (class implements "generalized matrix" method for data-driven background estimation)
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: GenMatrixBgEstFit.h,v 1.1 2009/06/11 07:23:28 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/GenMatrixBgEst1dPdf.h"
#include "TauAnalysis/BgEstimationTools/interface/BinGrid.h"

#include <TTree.h>
#include <TChain.h>

#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooAbsPdf.h>
#include <RooDataSet.h>

#include <string>
#include <vector>
#include <iostream>

class GenMatrixBgEstFit : public edm::EDAnalyzer
{
  struct objVarEntryType
  {
    objVarEntryType(const std::string&, float, float, float);
    ~objVarEntryType();
    std::string name_;
    RooRealVar* x_;
    float xBoundary_;
    float xMin_;
    float xMax_;  
  };

  struct pdf1dEntryType
  {
    pdf1dEntryType(const std::string&, objVarEntryType*, float, bool);
    ~pdf1dEntryType();
    std::string name_;
    GenMatrixBgEst1dPdf* pdf1d_;
    bool fixProb_;
    RooAbsReal* prob1_;
    RooAbsReal* prob2_;
  };

  struct pdfSingleProcessEntryType
  {
    pdfSingleProcessEntryType(const std::string&, const std::vector<objVarEntryType*>&, 
			      const std::vector<float>&, const std::vector<bool>&, float, bool);
    ~pdfSingleProcessEntryType();
    std::string name_;
    RooAbsPdf* pdfSingleProcess_;
    std::vector<pdf1dEntryType*> pdf1dEntries_;
    bool fixNorm_;
    RooAbsReal* norm_;
    int lineColor_;
    int lineStyle_;
    int lineWidth_;
  };

  struct pdfProcessSumEntryType
  {
    pdfProcessSumEntryType(const std::string&, std::vector<pdfSingleProcessEntryType*>&);
    ~pdfProcessSumEntryType();
    std::string name_;
    RooAbsPdf* pdfProcessSum_;
    std::vector<pdfSingleProcessEntryType*> pdfSingleProcessEntries_;
  };

  struct processTreeEntryType
  {
    processTreeEntryType(const std::string&, const std::vector<std::string>&, 
			 const std::string&, const std::string&, 
			 std::vector<objVarEntryType*>&, const std::string&);
    ~processTreeEntryType();
    std::string name_;
    TChain* allEventsTree_;
    TTree* selEventsTree_;
    TTree* goodEventsTree_;
    std::vector<Float_t> objVarValues_; 
    Float_t eventWeight_;
  };

 public:
  
  explicit GenMatrixBgEstFit(const edm::ParameterSet&);
  ~GenMatrixBgEstFit();
  
 private:

  void beginJob(const edm::EventSetup&) {}
  void analyze(const edm::Event&, const edm::EventSetup&) {}
  void endJob();

//--- private auxiliary functions
  double compNorm(processTreeEntryType*);
  double compProb(processTreeEntryType*, unsigned);
  void print(std::ostream& stream);
  void makeControlPlot(const RooRealVar*, const std::string&, const std::string&, const std::string&);
  void makeScaleFactors();

//--- configuration parameters
  typedef std::vector<std::string> vstring;
  vstring dataFileNames_;
  std::string treeName_;
  std::string treeSelection_;

  vstring branchNames_;
  std::string branchName_eventWeight_;
  RooRealVar* eventWeight_;

  std::string controlPlotsFileName_;

  typedef std::vector<double> vdouble;
  std::string scaleFactorFileName_;
  vdouble scaleFactorSignalRegion_;

//--- internal data-members for handling branches
  std::vector<objVarEntryType*> objVarEntries_;
  unsigned numVar_;

  std::map<std::string, processTreeEntryType*> treeData_;

  RooDataSet* dataSet_; 
  pdfProcessSumEntryType* pdfProcessSumEntry_;

  BinGrid* binGrid_;

  int cfgError_;
};

#endif  


