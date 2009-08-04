#include "TauAnalysis/BgEstimationTools/interface/BgEstMean.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

BgEstMean::BgEstMean(unsigned d)
  : mean_(d), auxDelta_(d)
{
  numVar_ = d;

  iValue_ = 0;
}

void BgEstMean::update(const TVectorD& value)
{
  if ( value.GetNoElements() != (int)numVar_ ) {
    edm::LogError("BgEstMean::update") << "Given value has invalid dimension = " << value.GetNoElements() << "," 
				       << " expected = " << numVar_ << " --> mean value will NOT be updated !!";
    return;
  }

  if ( iValue_ == 0 ) {
    for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
      mean_(iVar) = value(iVar);
    }
  } else {
    for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
      auxDelta_(iVar) = value(iVar) - mean_(iVar);
      mean_(iVar) += auxDelta_(iVar)/(iValue_ + 1.);
    }
  }

  ++iValue_;
}



