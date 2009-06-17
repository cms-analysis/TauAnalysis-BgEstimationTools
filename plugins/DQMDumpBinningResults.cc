#include "TauAnalysis/BgEstimationTools/plugins/DQMDumpBinningResults.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"
#include "TauAnalysis/DQMTools/interface/generalAuxFunctions.h"
#include "TauAnalysis/BgEstimationTools/interface/binningAuxFunctions.h"
#include "TauAnalysis/BgEstimationTools/interface/DataBinning.h"

DQMDumpBinningResults::DQMDumpBinningResults(const edm::ParameterSet& cfg)
{
  //std::cout << "<DQMDumpBinningResults::DQMDumpBinningResults>:" << std::endl;

  cfgError_ = 0;

  edm::ParameterSet dqmDirectoryEntries = cfg.getParameter<edm::ParameterSet>("dqmDirectories");
  vstring dqmDirectoryEntryNames = dqmDirectoryEntries.getParameterNamesForType<std::string>();
  for ( vstring::const_iterator dqmDirectoryEntryName = dqmDirectoryEntryNames.begin(); 
	dqmDirectoryEntryName != dqmDirectoryEntryNames.end(); ++dqmDirectoryEntryName ) {
    std::string dqmDirectoryEntry = dqmDirectoryEntries.getParameter<std::string>(*dqmDirectoryEntryName);

    processes_.push_back(*dqmDirectoryEntryName);

    dqmDirectories_[*dqmDirectoryEntryName] = dqmDirectoryEntry;
  }

  if ( processes_.size() == 0 ) {
    edm::LogError("DQMDumpBinningResults") << " Configuration Parameter dqmDirectories contains no Entries --> skipping !!";
    cfgError_ = 1;
  }
}

DQMDumpBinningResults::~DQMDumpBinningResults() 
{
  for ( std::map<std::string, BinningBase*>::iterator it = binningResults_.begin();
	it != binningResults_.end(); ++it ) {
    delete it->second;
  }
}

void DQMDumpBinningResults::analyze(const edm::Event&, const edm::EventSetup&)
{
//--- nothing to be done yet
}

BinningBase* DQMDumpBinningResults::loadBinning(const std::string& dqmDirectory)
{
//--- load all monitor elements present in dqmDirectory given as function argument

  std::vector<std::string> buffer;

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory);

  std::vector<std::string> meNames = dqmStore.getMEs();
  for ( std::vector<std::string>::const_iterator meName = meNames.begin();
	meName != meNames.end(); ++meName ) {
    std::string meName_full = dqmDirectoryName(dqmDirectory).append(*meName);
    //std::cout << " meName_full = " <<  meName_full << std::endl;

    dqmStore.setCurrentFolder(dqmDirectory);
    MonitorElement* me = dqmStore.get(meName_full);
    //std::cout << " me = " << me << std::endl;
    if ( !me ) {
      edm::LogError ("loadBinning") << " Failed to access meName = " << (*meName) << " in DQMStore" 
				    << " --> skipping !!";
      continue;
    }

//--- skip "invalid" MonitorElements
    if ( me->kind() == MonitorElement::DQM_KIND_INVALID ) {
      edm::LogWarning ("loadBinning") << " MonitorElement meName = " << (*meName) << " marked as invalid" 
				      << " --> skipping !!";
      continue;
    }

    if ( me->kind() == MonitorElement::DQM_KIND_STRING ) {
      std::string meValue = me->getStringValue();
      
      std::string entry = encodeBinningStringRep(*meName, "string", meValue);
      buffer.push_back(entry);
    } else if ( me->kind() == MonitorElement::DQM_KIND_REAL ) {
      std::ostringstream meValue;
      meValue << std::setprecision(3) << std::fixed << me->getFloatValue();
     
      std::string entry = encodeBinningStringRep(*meName, "float", meValue.str());
      buffer.push_back(entry);
    } else if ( me->kind() == MonitorElement::DQM_KIND_INT ) {
      std::ostringstream meValue;
      meValue << me->getIntValue();

      std::string entry = encodeBinningStringRep(*meName, "int", meValue.str());
      buffer.push_back(entry);
    } else {
      edm::LogError ("loadBinning") << " MonitorElement meName = " << (*meName) << " of unsupported type" 
				    << " --> skipping !!";
      continue;
    } 
  }

//--- only DataBinning objects supported so far...
  DataBinning* binning = new DataBinning();
  buffer << (*binning);

  return binning;
}

void DQMDumpBinningResults::endJob()
{
  //std::cout << "<DQMDumpBinningResults::endJob>:" << std::endl;

//--- check that configuration parameters contain no errors
  if ( cfgError_ ) {
    edm::LogError ("endjob") << " Error in Configuration ParameterSet --> FilterStatisticsTables will NOT be printed-out !!";
    return;
  }

//--- check that DQMStore service is available
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("endJob") << " Failed to access dqmStore --> FilterStatisticsTables will NOT be printed-out !!";
    return;
  }

//--- load objects inhertiting from BinningBase class from DQM directories
  for ( vstring::const_iterator process = processes_.begin();
	process != processes_.end(); ++process ) {
    const std::string& dqmDirectory = dqmDirectories_[*process];
    std::cout << "retrieving Binning results for process = " << (*process) 
	      << " from dqmDirectory = " << dqmDirectory << "..." << std::endl;

    BinningBase* binningResult = loadBinning(dqmDirectory); 

    if ( binningResult ) {
      binningResults_[*process] = binningResult;
    } else {
      edm::LogError ("DQMDumpBinningResults") << " Failed to load Binning result from dqmDirectory = " << dqmDirectory
					      << " --> Binning results will NOT be printed-out !!";
      return;
    }
  }

//--- print objects inhertiting from BinningBase class 
  for ( std::vector<std::string>::const_iterator process = processes_.begin();
	process != processes_.end(); ++process ) {
    BinningBase* binningResult = binningResults_[*process];
    std::cout << "Binning result for process = " << (*process) << ":" << std::endl;
    binningResult->print(std::cout);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DQMDumpBinningResults);
