#ifndef TauAnalysis_BgEstimationTools_BinningService_h  
#define TauAnalysis_BgEstimationTools_BinningService_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "TauAnalysis/BgEstimationTools/interface/BinningServiceBase.h"

#include <string>
#include <vector>

template<typename T>
class BinningService : public BinningServiceBase
{
 public: 
  explicit BinningService(const edm::ParameterSet&);
  ~BinningService();

 private:
  virtual BinningBase* createBinning() const;
};

#endif  


