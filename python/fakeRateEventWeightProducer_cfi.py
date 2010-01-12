import FWCore.ParameterSet.Config as cms

from TauAnalysis.BgEstimationTools.fakeRateJetWeightProducer_cfi import *

bgEstFakeRateEventWeights = cms.EDProducer("FakeRateEventWeightProducer",
    method = bgEstFakeRateJetWeights.method,
                                           
    allTauJetSource = bgEstFakeRateJetWeights.allTauJetSource,
    preselTauJetSource = bgEstFakeRateJetWeights.preselTauJetSource,
                                           
    dRmatch = bgEstFakeRateJetWeights.dRmatch,

    frTypes = bgEstFakeRateJetWeights.frTypes
)
