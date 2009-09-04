#ifndef TauAnalysis_BgEstimationTools_TF1Wrapper_h
#define TauAnalysis_BgEstimationTools_TF1Wrapper_h

/** \class TF1Wrapper
 *
 * Auxiliary class for creating TF1 objects via plugin factory
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: TF1Wrapper.h,v 1.1 2009/09/01 14:17:42 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/TF1WrapperBase.h"

class TF1Wrapper : public TF1WrapperBase
{
 public:
  // constructor 
  explicit TF1Wrapper(const edm::ParameterSet& cfg);
  
  // destructor
  virtual ~TF1Wrapper() {}
};

#endif  

