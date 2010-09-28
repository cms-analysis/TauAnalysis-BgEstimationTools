#ifndef TauAnalysis_BgEstimationTools_StringObjVectorValExtractor_h  
#define TauAnalysis_BgEstimationTools_StringObjVectorValExtractor_h

/** \class StringObjVectorValExtractor
 *
 * Auxiliary class for extracting scalar values from PAT objects
 * (used for Ntuple filling)
 *
 * NOTE: the values are extracted from the PAT object
 *       specified by the "index" configuration parameter (**first** PAT object in case "index" is not specified)
 *       contained in the collection specified by the "src" configuration parameter;
 *       in case the collection of PAT objects is empty, 
 *       a substitute value of -1. is returned by operator()
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: StringObjVectorValExtractor.h,v 1.3 2010/05/21 07:53:52 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"

template<typename T>
class StringObjVectorValExtractor : public ObjValVectorExtractorBase
{
 public:
  
  explicit StringObjVectorValExtractor(const edm::ParameterSet&);
  ~StringObjVectorValExtractor();
    
  std::vector<double> operator()(const edm::Event&) const;

 private:

//--- configuration parameters
  edm::InputTag src_;

  StringObjectFunction<T> stringObjFunction_;  
};

#endif  


