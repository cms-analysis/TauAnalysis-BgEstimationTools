import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.BgEstimationTools.fakeRateJetWeightProducer_cfi import *
from TauAnalysis.BgEstimationTools.fakeRateEventWeightProducer_cfi import *

#
# NOTE: preselection of reco::PFTaus need to correspond **exactly** to selection of pat::Taus
#       defined for Z --> muon + tau-jet channel in TauAnalysis/RecoTools/python/patLeptonSelection_cff.py !!
#
preselPFTausForMuTauFakeRateWeights = cms.EDFilter("PFTauSelector",
    src = cms.InputTag('shrinkingConePFTauProducer'),
    discriminators = cms.VPSet(
        cms.PSet(
            discriminator = cms.InputTag("shrinkingConePFTauDiscriminationAgainstMuon"),
            selectionCut = cms.double(0.5)
        )
    ),
    cut = cms.string("abs(eta) < 2.1 & pt > 20."), 
    filter = cms.bool(False)
)

bgEstFakeRateJetWeightsForMuTau = copy.deepcopy(bgEstFakeRateJetWeights)
bgEstFakeRateJetWeightsForMuTau.allTauJetSource = preselPFTausForMuTauFakeRateWeights.src
bgEstFakeRateJetWeightsForMuTau.preselTauJetSource = cms.InputTag('preselPFTausForMuTauFakeRateWeights')

bgEstFakeRateEventWeightsForMuTau = copy.deepcopy(bgEstFakeRateEventWeights)
bgEstFakeRateEventWeightsForMuTau.allTauJetSource = bgEstFakeRateJetWeightsForMuTau.allTauJetSource
bgEstFakeRateEventWeightsForMuTau.preselTauJetSource = bgEstFakeRateJetWeightsForMuTau.preselTauJetSource

produceFakeRateWeightsForMuTau = cms.Sequence(
    preselPFTausForMuTauFakeRateWeights
   * bgEstFakeRateJetWeightsForMuTau * bgEstFakeRateEventWeightsForMuTau
)
