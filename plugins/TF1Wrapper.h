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
 * $Id: TF1Wrapper.h,v 1.1 2009/06/11 07:23:28 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/TF1WrapperBase.h"

class TF1Wrapper : public TF1WrapperBase
{
 public:
  // constructor 
  explicit TF1Wrapper(const edm::ParameterSet& cfg)
    : TF1WrapperBase(cfg) {}
  
  // destructor
  virtual ~TF1Wrapper() {}
  
 protected:
  // auxiliary function to create TF1 object
  void makeTF1(const edm::ParameterSet&);
};

#endif  

