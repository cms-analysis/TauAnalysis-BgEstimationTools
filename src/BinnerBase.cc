#include "TauAnalysis/BgEstimationTools/interface/BinnerBase.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"
#include "TauAnalysis/BgEstimationTools/interface/binningAuxFunctions.h"

BinnerBase::BinnerBase(const edm::ParameterSet& cfg)
{
  //std::cout << "<BinnerBase::BinnerBase>:" << std::endl; 

  edm::ParameterSet cfgBinning = cfg.getParameter<edm::ParameterSet>("binning");
  objValExtractor_ = ObjValVectorExtractorPluginFactory::get()->create("MultiObjValExtractor", cfgBinning);

  binning_ = 0;

  dqmDirectory_store_ = cfg.getParameter<std::string>("dqmDirectory_store");
  //std::cout << " dqmDirectory_store = " << dqmDirectory_store_ << std::endl;
}

BinnerBase::~BinnerBase()
{
  delete objValExtractor_;
  delete binning_;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void BinnerBase::analyze(const edm::Event& evt, const edm::EventSetup& es) 
{ 
  bin(evt, es); 
}

void BinnerBase::endJob() 
{ 
  if ( binning_ ) {
    binning_->print(std::cout); 
    saveBinning(); 
  }
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

void BinnerBase::bin(const edm::Event& evt, const edm::EventSetup& es) 
{
  if ( !binning_ ) {
    edm::LogError ("BinnerBase::bin") << " No binning object defined --> skipping !!";
    return;
  }

  std::vector<double> x = (*objValExtractor_)(evt);

  binning_->bin(x);
}

void BinnerBase::saveBinning() const
{
  if ( !edm::Service<DQMStore>().isAvailable() ) {
    edm::LogError ("saveBinning") << " Failed to access dqmStore --> binning results will NOT be saved !!";
    return;
  }
  
  DQMStore& dqmStore = (*edm::Service<DQMStore>());
  
  dqmStore.setCurrentFolder(dqmDirectory_store_);

  std::vector<std::string> buffer;
  buffer << (*binning_);

  for ( std::vector<std::string>::const_iterator entry = buffer.begin();
	entry != buffer.end(); ++entry ) {
    std::string meName, meType, meValue;
    int error = 0;
    decodeBinningStringRep(*entry, meName, meType, meValue, error);

    if ( error ) {
      edm::LogError ("saveBinning") << " Error in parsing string = " << (*entry) << " --> skipping !!";
      continue;
    }

    if ( meType == "string" ) {
      dqmStore.bookString(meName, meValue);
    } else if ( meType == "float" ) {
      MonitorElement* me = dqmStore.bookFloat(meName);
      me->Fill(atof(meValue.data()));
    } else if ( meType == "int" ) {
      MonitorElement* me = dqmStore.bookInt(meName);
      me->Fill(atoi(meValue.data()));
    } else {
      edm::LogError ("saveBinning") << " Undefined meType = " << meType << " --> skipping !!";
      continue;
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(BinnerPluginFactory, "BinnerPluginFactory");


