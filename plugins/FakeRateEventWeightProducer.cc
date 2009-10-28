#include "TauAnalysis/BgEstimationTools/plugins/FakeRateEventWeightProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/TauReco/interface/PFTau.h"

FakeRateEventWeightProducer::FakeRateEventWeightProducer(const edm::ParameterSet& cfg)
  : FakeRateWeightProducerBase(cfg)
{
  produces<double>();
}

FakeRateEventWeightProducer::~FakeRateEventWeightProducer()
{
//--- nothing to be done yet...
}

void FakeRateEventWeightProducer::produce(edm::Event& evt, const edm::EventSetup&) 
{ 
  //std::cout << "<FakeRateEventWeightProducer::produce>:" << std::endl;

  if ( cfgError_ ) return;

  edm::Handle<edm::View<reco::BaseTau> > allTauJets;
  evt.getByLabel(allTauJetSource_, allTauJets);

  edm::Handle<edm::View<reco::Candidate> > preselTauJets;
  evt.getByLabel(preselTauJetSource_, preselTauJets);

  double tauJetIdEffSum = 0.;
  double qcdJetFakeRateSum = 0.;
  unsigned numTauJetDiscrPassed = 0;

  double qcdJetFakeRateProduct_complement = 1.;

  unsigned numTauJets = allTauJets->size();
  for ( unsigned iTauJet = 0; iTauJet < numTauJets; ++iTauJet ) {
    edm::RefToBase<reco::BaseTau> tauJetRef = allTauJets->refAt(iTauJet);

    double tauJetIdEff = 1.;
    double qcdJetFakeRate = 1.;

    bool tauJetDiscr_passed = true;

    getTauJetProperties(evt, tauJetRef, iTauJet, preselTauJets, tauJetIdEff, qcdJetFakeRate, tauJetDiscr_passed);

    //std::cout << "Pt = " << tauJetRef->pt() << ", eta = " << tauJetRef->eta() << ", tau id. discr = " << tauJetDiscr_passed << ":" 
    //	        << " tau id. eff = " << tauJetIdEff << ", fake-rate = " << qcdJetFakeRate << std::endl;

    tauJetIdEffSum += tauJetIdEff;
    qcdJetFakeRateSum += qcdJetFakeRate;

    if ( tauJetDiscr_passed ) ++numTauJetDiscrPassed;
    
    qcdJetFakeRateProduct_complement *= (1 - qcdJetFakeRate);
  }

  double fakeRateEventWeight = 0.;
  
  if ( method_ == "simple" ) {
    fakeRateEventWeight = (1 - qcdJetFakeRateProduct_complement);
  } else if ( method_ == "CDF" ) {
    double effTerm = 0.5*(tauJetIdEffSum + qcdJetFakeRateSum);
    double frTerm = qcdJetFakeRateSum;
    
    if ( effTerm > frTerm ) {
      fakeRateEventWeight = ( numTauJetDiscrPassed > 0 ) ?
	-frTerm*(1 - effTerm)/(effTerm - frTerm) : frTerm*effTerm/(effTerm - frTerm);
    }
  }
  
  //std::cout << " --> event weight = " << fakeRateEventWeight << std::endl;
  
  std::auto_ptr<double> fakeRateEventWeightPtr(new double(fakeRateEventWeight));
  
  evt.put(fakeRateEventWeightPtr);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(FakeRateEventWeightProducer);
