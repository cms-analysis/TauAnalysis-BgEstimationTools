#include "TauAnalysis/BgEstimationTools/interface/DataBinning.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <iomanip>
 
DataBinning::DataBinning(const edm::ParameterSet& cfg) 
  : BinningBase(cfg)
{
  numBins_ = binGrid_->numBins();

  binContents_.resize(numBins_);
  binSumw2_.resize(numBins_);
  for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
    binContents_[iBin] = 0.;
    binSumw2_[iBin] = 0.;
  }
}

DataBinning::~DataBinning()
{
//--- nothing to be done yet...
}

void DataBinning::bin(const std::vector<double>& x, double weight)
{
  int iBin = binGrid_->binNumber(x);
  if ( iBin >= 0 && iBin < (int)numBins_ ) {
    binContents_[iBin] += weight;
    binSumw2_[iBin] += weight*weight;
  }
}

void DataBinning::print(std::ostream& stream) const
{
  stream << "<DataBinning::print>:" << std::endl;
  stream << " name = " << name_ << std::endl;
  
  const std::vector<std::string>& objVarNames = binGrid_->objVarNames();
  
  for ( unsigned iBin = 0; iBin < numBins_; ++iBin ) {
    stream << " bin " << std::setw(2) << iBin << "(center: ";

    vdouble binCenter = binGrid_->binCenter(iBin);
    if ( binCenter.size() != objVarNames.size() ) {
      edm::LogError ("DataBinning::print") << "Invalid dimension of bin-center vector !!";
      return;
    }

    unsigned numObjVarNames = objVarNames.size();
    for ( unsigned iObjVar = 0; iObjVar < numObjVarNames; ++iObjVar ) {
      stream << objVarNames[iObjVar] << " = " << std::setprecision(3) << std::fixed << binCenter[iObjVar];
      if ( iObjVar < (numObjVarNames - 1) ) stream << ", ";
    }

    stream << "): " << std::setprecision(3) << std::fixed << binContents_[iBin] << " +/- " << binSumw2_[iBin] << std::endl;
  }
}

