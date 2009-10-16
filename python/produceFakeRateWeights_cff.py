import FWCore.ParameterSet.Config as cms

bgEstFakeRateJetWeights = cms.EDProducer("FakeRateJetWeightProducer",
    tauJetSource = cms.InputTag("shrinkingConePFTauProducer"),                         
    tauJetDiscriminators = cms.VPSet(
        cms.PSet(
            tauJetIdEffSource = cms.InputTag("shrinkingConeEfficienciesProducerFromFile", "effByECALIsolationZtautausim"),
            qcdJetFakeRateSource = cms.InputTag("shrinkingConeEfficienciesProducerFromFile", "frByECALIsolationMuEnrichedQCDsim"),
            tauJetDiscrSource = cms.InputTag("shrinkingConePFTauDiscriminationByIsolation")
        )
    )
)

bgEstFakeRateEventWeights = cms.EDProducer("FakeRateEventWeightProducer",
    tauJetSource = cms.InputTag('shrinkingConePFTauProducer'),                                              
    tauJetDiscriminators = bgEstFakeRateJetWeights.tauJetDiscriminators                       
)

produceFakeRateWeights = cms.Sequence( bgEstFakeRateJetWeights * bgEstFakeRateEventWeights )
