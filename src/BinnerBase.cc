#include "TauAnalysis/BgEstimationTools/interface/BinnerBase.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

BinnerBase::BinnerBase(const edm::ParameterSet& cfg)
{
  //std::cout << "<BinnerBase::BinnerBase>:" << std::endl; 

  objValExtractor_ = ObjValVectorExtractorPluginFactory::get()->create("MultiObjValExtractor", cfg);

  binning_ = 0;
}

BinnerBase::~BinnerBase()
{
  delete objValExtractor_;
  delete binning_;
}

void BinnerBase::bin(const edm::Event& evt, const edm::EventSetup& es) 
{
  if ( !binning_ ) {
    edm::LogError ("BinnerBase::bin") << " No binning object defined --> skipping !!";
    return;
  }

  std::vector<double> x = (*objValExtractor_)(evt);
  binning_->bin(x);
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(BinnerPluginFactory, "BinnerPluginFactory");

