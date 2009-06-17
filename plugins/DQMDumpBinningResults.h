#ifndef TauAnalysis_BgEstimationTools_DQMDumpBinningResults_h
#define TauAnalysis_BgEstimationTools_DQMDumpBinningResults_h

/** \class DQMDumpBinningResults
 *  
 *  Class to print-out binning results information contained in objects inheriting from BinningBase class
 *
 *  $Date: 2009/04/17 18:59:54 $
 *  $Revision: 1.1 $
 *  \author Christian Veelken, UC Davis
 */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/BinningBase.h"

#include <string>
#include <vector>
#include <map>

class DQMDumpBinningResults : public edm::EDAnalyzer
{
 public:
  explicit DQMDumpBinningResults(const edm::ParameterSet&);
  virtual ~DQMDumpBinningResults();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();  

private:
  BinningBase* loadBinning(const std::string&);

  typedef std::vector<std::string> vstring;
  vstring processes_;

  std::map<std::string, std::string> dqmDirectories_;

  std::map<std::string, BinningBase*> binningResults_;

  int cfgError_;
};

#endif


