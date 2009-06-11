#include "TauAnalysis/BgEstimationTools/interface/BinningBase.h"
 
BinningBase::BinningBase(const edm::ParameterSet& cfg)
{
  name_ = cfg.getParameter<std::string>("name");

  binGrid_ = new BinGrid(cfg);
}

BinningBase::~BinningBase()
{
  delete binGrid_;
}
