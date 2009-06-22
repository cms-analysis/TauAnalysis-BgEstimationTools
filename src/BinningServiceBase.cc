#include "TauAnalysis/BgEstimationTools/interface/BinningServiceBase.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "TauAnalysis/DQMTools/interface/dqmAuxFunctions.h"
#include "TauAnalysis/BgEstimationTools/interface/binningAuxFunctions.h"

#include <iostream>
#include <iomanip>
#include <algorithm>

const std::string meNameSeparator = ". :";

BinningServiceBase::BinningServiceBase(const edm::ParameterSet&)
{}

BinningServiceBase::~BinningServiceBase()
{
//--- nothing to be done yet...
}

bool operator<(const class BinningServiceBase::meEntry& entry1, const class BinningServiceBase::meEntry& entry2)
{
  return entry1.id_ < entry2.id_;
} 

BinningBase* BinningServiceBase::loadBinningResults(const std::string& dqmDirectory) const
{
//--- check if DQMStore is available;
//    print an error message and return if not
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("loadBinningResults") << " Failed to access dqmStore !!";
    return 0;
  }

  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  dqmStore.setCurrentFolder(dqmDirectory);

  std::vector<meEntry> meEntries;

  std::vector<std::string> meNames = dqmStore.getMEs();
  for ( std::vector<std::string>::const_iterator meName = meNames.begin();
	meName != meNames.end(); ++meName ) {
    std::string meName_full = dqmDirectoryName(dqmDirectory).append(*meName);
    std::cout << " meName_full = " <<  meName_full << std::endl;

    dqmStore.setCurrentFolder(dqmDirectory);
    MonitorElement* me = dqmStore.get(meName_full);
    std::cout << " me = " << me << std::endl;
    if ( !me ) {
      edm::LogError ("loadBinningResults") << " Failed to access meName = " << (*meName) << " in DQMStore" 
					   << " --> skipping !!";
      continue;
    }

//--- skip "invalid" MonitorElements
    if ( me->kind() == MonitorElement::DQM_KIND_INVALID ) {
      edm::LogWarning ("loadBinningResults") << " MonitorElement meName = " << (*meName) << " marked as invalid" 
					     << " --> skipping !!";
      continue;
    }

//--- decode id of MonitorElement,
//    so that MonitorElements loaded from DQMStore
//    can be sorted in the order in which they have been saved
    size_t posSeparator = meName->find(meNameSeparator);
    int meId = atoi(std::string(*meName, 0, posSeparator).data());
    std::string meName_decoded = std::string(*meName, posSeparator + meNameSeparator.length());
    std::cout << "meName_decoded = " << meName_decoded << std::endl;

    std::string meType, meValue;

    if ( me->kind() == MonitorElement::DQM_KIND_STRING ) {
      meType = "string";
      meValue = me->getStringValue();
      std::cout << "meValue(string) = " << meValue << std::endl;
    } else if ( me->kind() == MonitorElement::DQM_KIND_REAL ) {
      meType = "float";
      std::ostringstream meValue_ostringstream;
      meValue_ostringstream << std::setprecision(3) << std::fixed << me->getFloatValue();
      meValue = meValue_ostringstream.str();
      std::cout << "meValue(float) = " << meValue << std::endl;
    } else if ( me->kind() == MonitorElement::DQM_KIND_INT ) {
      meType = "int";
      std::ostringstream meValue_ostringstream;
      meValue_ostringstream << me->getIntValue();
      meValue = meValue_ostringstream.str();
      std::cout << "meValue(int) = " << meValue << std::endl;
    } else {
      edm::LogError ("loadBinningResults") << " MonitorElement meName = " << (*meName) << " of unsupported type" 
					   << " --> skipping !!";
      continue;
    } 

    meEntries.push_back(meEntry(meId, meName_decoded, meType, meValue));
  }

//--- sort meEntry objects representing MonitorElements
//    in the order by which MonitorElements have been saved
  std::sort(meEntries.begin(), meEntries.end()); 

//--- fill meEntry objects into buffer in sorted order
  std::vector<std::string> buffer;
  for ( std::vector<meEntry>::const_iterator meEntry_i = meEntries.begin();
	meEntry_i != meEntries.end(); ++meEntry_i ) {
    std::string entry = encodeBinningStringRep(meEntry_i->name_, meEntry_i->type_, meEntry_i->value_);
    std::cout << "entry = " << entry << std::endl;
    buffer.push_back(entry);
  }

//--- create binning object and initialize with content of buffer
  BinningBase* binning = createBinning();
  buffer >> (*binning);

  return binning;
}

void BinningServiceBase::saveBinningResults(const std::string& dqmDirectory, const BinningBase* binning) const
{
//--- check if DQMStore is available;
//    print an error message and return if not
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("saveBinningResults") << " Failed to access dqmStore --> binning results will NOT be saved !!";
    return;
  }
  
  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  
  dqmStore.setCurrentFolder(dqmDirectory);

  std::vector<std::string> buffer;
  buffer << (*binning);

  int id = 1;
  for ( std::vector<std::string>::const_iterator entry = buffer.begin();
	entry != buffer.end(); ++entry ) {
    std::string meName, meType, meValue;
    int error = 0;
    decodeBinningStringRep(*entry, meName, meType, meValue, error);

    if ( error ) {
      edm::LogError ("saveBinningResults") << " Error in parsing string = " << (*entry) << " --> skipping !!";
      continue;
    }

//--- encode id of MonitorElement,
//    so that MonitorElements loaded from DQMStore
//    can be sorted in the order in which they have been saved
    std::ostringstream meName_encoded_ostringstream;
    meName_encoded_ostringstream << std::setw(3) << id << meNameSeparator << meName;
    std::string meName_encoded = meName_encoded_ostringstream .str();
    //std::cout << "meName_encoded = " << meName_encoded << std::endl;

    if ( meType == "string" ) {
      dqmStore.bookString(meName_encoded, meValue);
    } else if ( meType == "float" ) {
      MonitorElement* me = dqmStore.bookFloat(meName_encoded);
      me->Fill(atof(meValue.data()));
    } else if ( meType == "int" ) {
      MonitorElement* me = dqmStore.bookInt(meName_encoded);
      me->Fill(atoi(meValue.data()));
    } else {
      edm::LogError ("saveBinningResults") << " Undefined meType = " << meType << " --> skipping !!";
      continue;
    }

    ++id;
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(BinningServicePluginFactory, "BinningServicePluginFactory");
