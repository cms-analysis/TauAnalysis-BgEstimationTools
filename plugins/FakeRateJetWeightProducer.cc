#include "TauAnalysis/BgEstimationTools/plugins/FakeRateJetWeightProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/TauReco/interface/PFTau.h"

#include "DataFormats/Common/interface/ValueMap.h" 
#include "DataFormats/PatCandidates/interface/LookupTableRecord.h"

typedef edm::ValueMap<pat::LookupTableRecord> LookupTableMap;

FakeRateJetWeightProducer::FakeRateJetWeightProducer(const edm::ParameterSet& cfg)
  : FakeRateWeightProducerBase(cfg)
{
  produces<LookupTableMap>();
}

FakeRateJetWeightProducer::~FakeRateJetWeightProducer()
{
//--- nothing to be done yet...
}

void FakeRateJetWeightProducer::produce(edm::Event& evt, const edm::EventSetup&) 
{ 
  //std::cout << "<FakeRateJetWeightProducer::produce>:" << std::endl;

  if ( cfgError_ ) return;

  edm::Handle<edm::View<reco::BaseTau> > tauJets;
  evt.getByLabel(tauJetSource_, tauJets);

  std::vector<pat::LookupTableRecord> fakeRateJetWeights;

  unsigned numTauJets = tauJets->size();
  for ( unsigned iTauJet = 0; iTauJet < numTauJets; ++iTauJet ) {
    edm::RefToBase<reco::BaseTau> tauJetRef = tauJets->refAt(iTauJet);

    double tauJetIdEff = 1.;
    double qcdJetFakeRate = 1.;

    bool tauJetDiscr_passed = true;

    getTauJetProperties(evt, tauJetRef, iTauJet, tauJetIdEff, qcdJetFakeRate, tauJetDiscr_passed);

    double fakeRateJetWeight = 0.;
    if ( tauJetIdEff > qcdJetFakeRate ) {
      fakeRateJetWeight = ( tauJetDiscr_passed ) ? 
	-qcdJetFakeRate*(1 - tauJetIdEff)/(tauJetIdEff - qcdJetFakeRate) : qcdJetFakeRate*tauJetIdEff/(tauJetIdEff - qcdJetFakeRate);
    }

    //std::cout << " --> jet weight = " << fakeRateJetWeight << std::endl;

    fakeRateJetWeights.push_back(pat::LookupTableRecord(fakeRateJetWeight));
  }

  std::auto_ptr<LookupTableMap> fakeRateJetWeightMap(new LookupTableMap());
  LookupTableMap::Filler valueMapFiller(*fakeRateJetWeightMap);
  valueMapFiller.insert(tauJets, fakeRateJetWeights.begin(), fakeRateJetWeights.end());
  valueMapFiller.fill();

  evt.put(fakeRateJetWeightMap);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(FakeRateJetWeightProducer);
