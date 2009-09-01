#ifndef TauAnalysis_BgEstimationTools_TF1landauXgausWrapper_h
#define TauAnalysis_BgEstimationTools_TF1landauXgausWrapper_h

/** \class TF1landauXgausWrapper
 *
 * Numerical implementation of Landau density convoluted with Gaussian
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: TF1landauXgausWrapper.h,v 1.1 2009/06/11 07:23:28 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/TF1WrapperBase.h"

class TF1landauXgausWrapper : public TF1WrapperBase
{
 public:
  // constructor 
  explicit TF1landauXgausWrapper(const edm::ParameterSet& cfg)
    : TF1WrapperBase(cfg) {}
  
  // destructor
  ~TF1landauXgausWrapper() {}

 protected:
  // auxiliary function to create TF1 object
  void makeTF1(const edm::ParameterSet&);
};

#endif  

