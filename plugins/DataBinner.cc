#include "TauAnalysis/BgEstimationTools/plugins/DataBinner.h"

#include "TauAnalysis/BgEstimationTools/interface/DataBinning.h"

#include <iostream>

DataBinner::DataBinner(const edm::ParameterSet& cfg)
  : BinnerBase(cfg)
{
  //std::cout << "<DataBinner::DataBinner>:" << std::endl; 

  binning_ = new DataBinning(cfg);
}

DataBinner::~DataBinner()
{
//--- nothing to be done yet...
}
