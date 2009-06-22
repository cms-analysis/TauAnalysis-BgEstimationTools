#ifndef TauAnalysis_BgEstimationTools_BinnerBase_h  
#define TauAnalysis_BgEstimationTools_BinnerBase_h

/** \class BinnerBase
 *
 * Pure virtual base-class for extracting observables to be binned from the event
 * and filling them into object derrived from BinningBase,
 * using the bin-grid associated with the binning object.
 * 
 * \author Christian Veelken, UC Davis
 *         (inspired by code written for H1 by Paul Laycock, University of Liverpool)
 *
 * \version $Revision: 1.2 $
 *
 * $Id: BinnerBase.h,v 1.2 2009/06/17 12:39:10 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/Core/interface/AnalyzerPluginBase.h"
#include "TauAnalysis/BgEstimationTools/interface/ObjValVectorExtractorBase.h"
#include "TauAnalysis/BgEstimationTools/interface/BinningBase.h"
#include "TauAnalysis/BgEstimationTools/interface/BinningServiceBase.h"

#include <vector>
#include <string>

class BinnerBase : public AnalyzerPluginBase
{
 public: 
  explicit BinnerBase(const edm::ParameterSet&);
  virtual ~BinnerBase();

  const BinningBase* getBinning() const { return binning_; }

  virtual void beginJob() {}
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

 protected:
  virtual void bin(const edm::Event&, const edm::EventSetup&);

  ObjValVectorExtractorBase* objValExtractor_;

  BinningBase* binning_;

  std::string dqmDirectory_store_;

  BinningServiceBase* binningService_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<BinnerBase* (const edm::ParameterSet&)> BinnerPluginFactory;

#endif  


