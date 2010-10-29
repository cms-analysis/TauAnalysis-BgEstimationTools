import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

#--------------------------------------------------------------------------------
# select loosely selected Z --> tau+ tau- --> muon + tau-jet event sample
# for measuring tau id. efficiency
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------  
# produce collection of pat::Muons
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patMuonSelection_cfi import *

muonsForTauIdEffZtoMuTauCombinedFitCombRelIsoLoose = cms.EDFilter("PATMuonSelector",
    cut = cms.string('(userIsolation("pat::TrackIso") + userIsolation("pat::EcalIso")) < (0.24*pt)'),
    filter = cms.bool(False)
)

muonsForTauIdEffZtoMuTauCombinedFitPionVetoLooseIso = copy.deepcopy(selectedPatMuonsPionVeto)

muonsForTauIdEffZtoMuTauCombinedFitTrkLooseIso = copy.deepcopy(selectedPatMuonsTrk)

muonsForTauIdEffZtoMuTauCombinedFitTrkIPlooseIso = copy.deepcopy(selectedPatMuonsTrkIP)

muonSelConfiguratorTauIdEffZtoMuTauCombinedFitLooseIso = objSelConfigurator(
    [ muonsForTauIdEffZtoMuTauCombinedFitCombRelIsoLoose,
      muonsForTauIdEffZtoMuTauCombinedFitPionVetoLooseIso,
      muonsForTauIdEffZtoMuTauCombinedFitTrkLooseIso,
      muonsForTauIdEffZtoMuTauCombinedFitTrkIPlooseIso ],
    src = "selectedPatMuonsPt10Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuonsForTauIdEffZtoMuTauCombinedFitLooseIso = muonSelConfiguratorTauIdEffZtoMuTauCombinedFitLooseIso.configure(pyNameSpace = locals())

muonsForTauIdEffZtoMuTauCombinedFitCombRelIsoTight = cms.EDFilter("PATMuonSelector",
    cut = cms.string('(userIsolation("pat::TrackIso") + userIsolation("pat::EcalIso")) < (0.06*pt)'),
    filter = cms.bool(False)
)

muonsForTauIdEffZtoMuTauCombinedFitPionVetoTightIso = copy.deepcopy(selectedPatMuonsPionVeto)

muonsForTauIdEffZtoMuTauCombinedFitTrkTightIso = copy.deepcopy(selectedPatMuonsTrk)

muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIso = copy.deepcopy(selectedPatMuonsTrkIP)

muonSelConfiguratorTauIdEffZtoMuTauCombinedFitTightIso = objSelConfigurator(
    [ muonsForTauIdEffZtoMuTauCombinedFitCombRelIsoTight,
      muonsForTauIdEffZtoMuTauCombinedFitPionVetoTightIso,
      muonsForTauIdEffZtoMuTauCombinedFitTrkTightIso,
      muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIso ],
    src = "selectedPatMuonsPt10Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuonsForTauIdEffZtoMuTauCombinedFitTightIso = muonSelConfiguratorTauIdEffZtoMuTauCombinedFitTightIso.configure(pyNameSpace = locals())

selectMuonsForTauIdEffZtoMuTauCombinedFit = cms.Sequence(
    selectMuonsForTauIdEffZtoMuTauCombinedFitLooseIso
   * selectMuonsForTauIdEffZtoMuTauCombinedFitTightIso
)    

#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patPFTauSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *

tausForTauIdEffZtoMuTauCombinedFitMuonVeto = copy.deepcopy(selectedPatTausMuonVeto)

tauSelConfiguratorTauIdEffZtoMuTauCombinedFit = objSelConfigurator(
    [ tausForTauIdEffZtoMuTauCombinedFitMuonVeto ],
    src = "selectedPatTausForMuTauLeadTrkPtCumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectTausForTauIdEffZtoMuTauCombinedFit = tauSelConfiguratorTauIdEffZtoMuTauCombinedFit.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of Z --> muon+ muon- hypotheses
#--------------------------------------------------------------------------------

# require muon candidates considered for Z --> mu+ mu- hypothesis
# to be reconstructed in muon system
# (with or without a track reconstructed in Pixel/SiStrip tracking detectors linked to it)
muonsLooseForZmumuHypothesesTauIdEffZtoMuTauCombinedFit = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string('isGlobalMuon | isStandAloneMuon | pt > 5.'),
    filter = cms.bool(False)
)

muonsTightForZmumuHypothesesTauIdEffZtoMuTauCombinedFit = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string('isGlobalMuon | isStandAloneMuon'),
    filter = cms.bool(False)
)

allDiMuPairZmumuHypothesesTauIdEffZtoMuTauCombinedFit = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsTightForZmumuHypothesesTauIdEffZtoMuTauCombinedFit'),
    srcLeg2 = cms.InputTag('muonsLooseForZmumuHypothesesTauIdEffZtoMuTauCombinedFit'),
    dRmin12 = cms.double(0.5),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                        
    verbosity = cms.untracked.int32(0)
)

selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauCombinedFit = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag("allDiMuPairZmumuHypothesesTauIdEffZtoMuTauCombinedFit"),                                   
    cut = cms.string('p4Vis.mass > 70. & p4Vis.mass < 110.'),
    filter = cms.bool(False)
)

produceDiMuPairsTauIdEffZtoMuTauCombinedFit = cms.Sequence(
    muonsLooseForZmumuHypothesesTauIdEffZtoMuTauCombinedFit * muonsTightForZmumuHypothesesTauIdEffZtoMuTauCombinedFit
   * allDiMuPairZmumuHypothesesTauIdEffZtoMuTauCombinedFit
   * selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauCombinedFit
)

#--------------------------------------------------------------------------------  
# produce collection of muon + tau-jet combinations
#--------------------------------------------------------------------------------

from TauAnalysis.CandidateTools.muTauPairProduction_cff import *
from TauAnalysis.CandidateTools.resolutions_cfi import *

muTauPairsTauIdEffZtoMuTauCombinedFit = allMuTauPairs.clone(
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIsoCumulative'),
    srcLeg2 = cms.InputTag('tausForTauIdEffZtoMuTauCombinedFitMuonVetoCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('patMETs'),
    recoMode = cms.string(""),
    doSVreco = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauCombinedFitBackToBack = cms.EDFilter("PATMuTauPairSelector",
    cut = cms.string('dPhi12 > 2.793'),
    filter = cms.bool(False)
)

muTauPairConfiguratorTauIdEffZtoMuTauCombinedFit = objSelConfigurator(
    [ muTauPairsTauIdEffZtoMuTauCombinedFitBackToBack ],
    src = "muTauPairsTauIdEffZtoMuTauCombinedFit",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuTauPairsTauIdEffZtoMuTauCombinedFit = \
  muTauPairConfiguratorTauIdEffZtoMuTauCombinedFit.configure(pyNameSpace = locals())

muTauPairsTauIdEffZtoMuTauCombinedFitWplusJets = copy.deepcopy(muTauPairsTauIdEffZtoMuTauCombinedFit)

muTauPairsTauIdEffZtoMuTauCombinedFitMt1MEtWplusJets = cms.EDFilter("PATMuTauPairSelector",
    cut = cms.string('mt1MET > 50.'),
    filter = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffWplusJets = cms.EDFilter("PATMuTauPairSelector",
    cut = cms.string('(pZeta - 1.5*pZetaVis) < -20.'),
    filter = cms.bool(False)
)

muTauPairConfiguratorTauIdEffZtoMuTauCombinedFitWplusJets = objSelConfigurator(
    [ muTauPairsTauIdEffZtoMuTauCombinedFitMt1MEtWplusJets,
      muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffWplusJets ],
    src = "muTauPairsTauIdEffZtoMuTauCombinedFitWplusJets",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuTauPairsTauIdEffZtoMuTauCombinedFitWplusJets = \
  muTauPairConfiguratorTauIdEffZtoMuTauCombinedFitWplusJets.configure(pyNameSpace = locals())

muTauPairsTauIdEffZtoMuTauCombinedFitQCD = muTauPairsTauIdEffZtoMuTauCombinedFit.clone(
    srcLeg1 = cms.InputTag('muonsForTauIdEffZtoMuTauCombinedFitTrkIPlooseIsoCumulative')
)

muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackQCD = copy.deepcopy(muTauPairsTauIdEffZtoMuTauCombinedFitBackToBack)

muTauPairsTauIdEffZtoMuTauCombinedFitMt1MEtQCD = muTauPairsTauIdEffZtoMuTauCombinedFitMt1MEtWplusJets.clone(
    cut = cms.string('mt1MET < 30.')
)

muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffQCD = muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffWplusJets.clone(
    cut = cms.string('(pZeta - 1.5*pZetaVis) > -20.')
)

muTauPairConfiguratorTauIdEffZtoMuTauCombinedFitQCD = objSelConfigurator(
    [ muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackQCD,
      muTauPairsTauIdEffZtoMuTauCombinedFitMt1MEtQCD,
      muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffQCD ],
    src = "muTauPairsTauIdEffZtoMuTauCombinedFitQCD",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuTauPairsTauIdEffZtoMuTauCombinedFitQCD = \
  muTauPairConfiguratorTauIdEffZtoMuTauCombinedFitQCD.configure(pyNameSpace = locals())

produceMuTauPairsTauIdEffZtoMuTauCombinedFit = cms.Sequence(
    muTauPairsTauIdEffZtoMuTauCombinedFit * selectMuTauPairsTauIdEffZtoMuTauCombinedFit
   * muTauPairsTauIdEffZtoMuTauCombinedFitWplusJets * selectMuTauPairsTauIdEffZtoMuTauCombinedFitWplusJets
   * muTauPairsTauIdEffZtoMuTauCombinedFitQCD * selectMuTauPairsTauIdEffZtoMuTauCombinedFitQCD 
)    

#--------------------------------------------------------------------------------  
# produce collection of central jets not overlapping with muon or tau-jet
#--------------------------------------------------------------------------------

jetsTauIdEffZtoMuTauCombinedFitEta25 = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) < 2.5'),
    filter = cms.bool(False)
)

jetsTauIdEffZtoMuTauCombinedFitEt15 = cms.EDFilter("PATJetSelector",
    cut = cms.string('et > 15.'), 
    filter = cms.bool(False)
)

jetsTauIdEffZtoMuTauCombinedFitAntiOverlapWithLeptonsVeto = cms.EDFilter("PATJetAntiOverlapSelector",
    srcNotToBeFiltered = cms.VInputTag(
        "muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIsoCumulative",
        "tausForTauIdEffZtoMuTauCombinedFitMuonVetoCumulative"
    ),                                                           
    dRmin = cms.double(0.7),
    filter = cms.bool(False)                                           
)

jetSelConfiguratorTauIdEffZtoMuTauCombinedFit = objSelConfigurator(
    [ jetsTauIdEffZtoMuTauCombinedFitEta25,
      jetsTauIdEffZtoMuTauCombinedFitEt15,
      jetsTauIdEffZtoMuTauCombinedFitAntiOverlapWithLeptonsVeto ],
    src = "cleanPatJets",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectJetsForTauIdEffZtoMuTauCombinedFit = jetSelConfiguratorTauIdEffZtoMuTauCombinedFit.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *

muonHistManagerTauIdEffZtoMuTauCombinedFit = copy.deepcopy(muonHistManager)
muonHistManagerTauIdEffZtoMuTauCombinedFit.pluginName = 'muonHistManagerTauIdEffZtoMuTauCombinedFit'
muonHistManagerTauIdEffZtoMuTauCombinedFit.muonSource = 'muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIsoCumulative'

tauHistManagerTauIdEffZtoMuTauCombinedFit = copy.deepcopy(tauHistManager)
tauHistManagerTauIdEffZtoMuTauCombinedFit.pluginName = 'tauHistManagerTauIdEffZtoMuTauCombinedFit'
tauHistManagerTauIdEffZtoMuTauCombinedFit.tauSource = 'tausForTauIdEffZtoMuTauCombinedFitMuonVetoCumulative'

diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit.pluginName = \
  'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit'
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit.diTauCandidateSource = \
  'muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackCumulative'
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit.visMassHypothesisSource = cms.InputTag('')
from TauAnalysis.Core.diTauCandidateEventActivityHistManager_cfi import *
diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit = copy.deepcopy(diTauCandidateEventActivityHistManager)
diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit.pluginName = \
  'diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit'
diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit.pluginType = 'PATMuTauPairEventActivityHistManager'
diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit.diTauCandidateSource = \
  'muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackCumulative'

diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets.pluginName = \
  'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets'
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets.diTauCandidateSource = \
  'muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffWplusJetsCumulative'

diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD.pluginName = \
  'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD'
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD.diTauCandidateSource = \
  'muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffQCDcumulative'

jetHistManagerTauIdEffZtoMuTauCombinedFit = copy.deepcopy(jetHistManager)
jetHistManagerTauIdEffZtoMuTauCombinedFit.pluginName = 'jetHistManagerTauIdEffZtoMuTauCombinedFit'
jetHistManagerTauIdEffZtoMuTauCombinedFit.jetSource = 'jetsTauIdEffZtoMuTauCombinedFitAntiOverlapWithLeptonsVetoCumulative'

from TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauHistManager_cfi import *
tauIdEffZtoMuTauHistManagerCombinedFit = copy.deepcopy(tauIdEffZtoMuTauHistManager)
tauIdEffZtoMuTauHistManagerCombinedFit.pluginName = 'tauIdEffZtoMuTauHistManagerCombinedFit'
tauIdEffZtoMuTauHistManagerCombinedFit.muonSource = 'muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIsoCumulative'
tauIdEffZtoMuTauHistManagerCombinedFit.tauSource = 'tausForTauIdEffZtoMuTauCombinedFitMuonVetoCumulative'
tauIdEffZtoMuTauHistManagerCombinedFit.diTauSource = 'muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackCumulative'
tauIdEffZtoMuTauHistManagerCombinedFit.diTauChargeSignExtractor.src = tauIdEffZtoMuTauHistManagerCombinedFit.diTauSource

tauIdEffZtoMuTauHistManagerCombinedFitWplusJets = copy.deepcopy(tauIdEffZtoMuTauHistManagerCombinedFit)
tauIdEffZtoMuTauHistManagerCombinedFitWplusJets.pluginName = 'tauIdEffZtoMuTauHistManagerCombinedFitWplusJets'
tauIdEffZtoMuTauHistManagerCombinedFitWplusJets.diTauSource = 'muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffWplusJetsCumulative'
tauIdEffZtoMuTauHistManagerCombinedFitWplusJets.diTauChargeSignExtractor.src = tauIdEffZtoMuTauHistManagerCombinedFitWplusJets.diTauSource

tauIdEffZtoMuTauHistManagerCombinedFitQCD = copy.deepcopy(tauIdEffZtoMuTauHistManagerCombinedFit)
tauIdEffZtoMuTauHistManagerCombinedFitQCD.pluginName ='tauIdEffZtoMuTauHistManagerCombinedFitQCD'
tauIdEffZtoMuTauHistManagerCombinedFitQCD.diTauSource = 'muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffQCDcumulative'
tauIdEffZtoMuTauHistManagerCombinedFitQCD.diTauChargeSignExtractor.src = tauIdEffZtoMuTauHistManagerCombinedFitQCD.diTauSource

from TauAnalysis.BgEstimationTools.bgEstBinGridZtoMuTau_cfi import *

dataBinnerTauIdEffZtoMuTauCombinedFit = copy.deepcopy(dataBinner)
dataBinnerTauIdEffZtoMuTauCombinedFit.pluginName = 'dataBinnerTauIdEffZtoMuTauCombinedFit'

binningTauIdEffZtoMuTauCombinedFit_ewkTauId = copy.deepcopy(binning_ewkTauId)
binningTauIdEffZtoMuTauCombinedFit_ewkTauId.extractor.src = 'tausForTauIdEffZtoMuTauCombinedFitMuonVetoCumulative'

binningTauIdEffZtoMuTauCombinedFit_relMuonIso = copy.deepcopy(binning_relMuonIso)
binningTauIdEffZtoMuTauCombinedFit_relMuonIso.extractor.src = 'muonsForTauIdEffZtoMuTauCombinedFitTrkIPlooseIsoCumulative'
binningTauIdEffZtoMuTauCombinedFit_relMuonIso.binning = cms.PSet(
    boundaries = cms.vdouble(0.06, 0.12),
    min = cms.double(-0.01),
    max = cms.double(0.24)
)

binningTauIdEffZtoMuTauCombinedFit_diTauAbsCharge = copy.deepcopy(binning_diTauAbsCharge)
binningTauIdEffZtoMuTauCombinedFit_diTauAbsCharge.extractor.src = \
  'muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackCumulative'

binningTauIdEffZtoMuTauCombinedFitWplusJets_diTauAbsCharge = copy.deepcopy(binning_diTauAbsCharge)
binningTauIdEffZtoMuTauCombinedFitWplusJets_diTauAbsCharge.extractor.src = \
  'muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffWplusJetsCumulative'

binningTauIdEffZtoMuTauCombinedFitQCD_diTauAbsCharge = copy.deepcopy(binning_diTauAbsCharge)
binningTauIdEffZtoMuTauCombinedFitQCD_diTauAbsCharge.extractor.src = \
  'muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffQCDcumulative'

tauIdEffBinningZtoMuTau_comb2d = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_comb2d"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauCombinedFit_ewkTauId,
        binningTauIdEffZtoMuTauCombinedFit_diTauAbsCharge
    )
)

tauIdEffBinningZtoMuTau_comb2dWplusJets = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_comb2dWplusJets"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauCombinedFit_ewkTauId,
        binningTauIdEffZtoMuTauCombinedFitWplusJets_diTauAbsCharge
    )
)

tauIdEffBinningZtoMuTau_comb3dQCD = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_comb3dQCD"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauCombinedFit_ewkTauId,
        binningTauIdEffZtoMuTauCombinedFit_relMuonIso,
        binningTauIdEffZtoMuTauCombinedFitQCD_diTauAbsCharge
    )
)

analyzeEventsTauIdEffZtoMuTauCombinedFit = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('TauIdEffAnalyzerZtoMuTauCombinedFit'), 
                            
    filters = cms.VPSet(
        evtSelGenPhaseSpace,
        evtSelTrigger,
        evtSelPrimaryEventVertex,
        evtSelPrimaryEventVertexQuality,
        evtSelPrimaryEventVertexPosition,
        evtSelGlobalMuon,
        evtSelMuonEta,
        evtSelMuonPt,
        cms.PSet(
            pluginName = cms.string('muonCombRelIsoLooseCutTauIdEffZtoMuTauCombinedFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauCombinedFitCombRelIsoLooseCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muonAntiPionCutTauIdEffZtoMuTauCombinedFitLooseIso'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauCombinedFitPionVetoLooseIsoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muonTrkIPcutTauIdEffZtoMuTauCombinedFitLooseIso'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauCombinedFitTrkLooseIsoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muonCombRelIsoTightCutTauIdEffZtoMuTauCombinedFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauCombinedFitCombRelIsoTightCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muonAntiPionCutTauIdEffZtoMuTauCombinedFitTightIso'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauCombinedFitPionVetoTightIsoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muonTrkIPcutTauIdEffZtoMuTauCombinedFitTightIso'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauCombinedFitTrkTightIsoCumulative'),
            minNumber = cms.uint32(1)
        ),
        evtSelTauAntiOverlapWithMuonsVeto,
        evtSelTauEta,
        evtSelTauPt,
        evtSelTauLeadTrk,
        evtSelTauLeadTrkPt,
        cms.PSet(
            pluginName = cms.string('tauMuonVetoTauIdEffZtoMuTauCombinedFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('tausForTauIdEffZtoMuTauCombinedFitMuonVetoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauCombinedFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFit'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairBackToBackTauIdEffZtoMuTauCombinedFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauCombinedFitWplusJets'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitWplusJets'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairMt1MEtTauIdEffZtoMuTauCombinedFitWplusJets'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitMt1MEtWplusJetsCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauCombinedFitWplusJets'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffWplusJetsCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauCombinedFitQCD'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitQCD'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairBackToBackTauIdEffZtoMuTauCombinedFitQCD'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackQCDcumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairMt1MEtTauIdEffZtoMuTauCombinedFitQCD'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitMt1MEtQCDcumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauCombinedFitQCD'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffQCDcumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('diMuPairZmumuHypothesisVeto'),
            pluginType = cms.string('PATCandViewMaxEventSelector'),
            src = cms.InputTag('selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauCombinedFit'),
            maxNumber = cms.uint32(0)
        ),
        ##cms.PSet(
        ##    pluginName = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTauCombinedFit'),
        ##    pluginType = cms.string('PATCandViewMaxEventSelector'),
        ##    src = cms.InputTag('tausForTauIdEffZtoMuTauCombinedFitMuonVeto'),
        ##    maxNumber = cms.uint32(0)
        ##),
        cms.PSet(
            pluginName = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'),
            pluginType = cms.string('PATCandViewMaxEventSelector'),
            src = cms.InputTag('muonsTightForZmumuHypothesesTauIdEffZtoMuTauCombinedFit'),
            maxNumber = cms.uint32(1)
        )
    ),
  
    analyzers = cms.VPSet(
        muonHistManagerTauIdEffZtoMuTauCombinedFit,
        tauHistManagerTauIdEffZtoMuTauCombinedFit,
        diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit,
        diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets,
        diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD,
        diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit,
        jetHistManagerTauIdEffZtoMuTauCombinedFit,
        caloMEtHistManager,
        pfMEtHistManager,
        tauIdEffZtoMuTauHistManagerCombinedFit,
        tauIdEffZtoMuTauHistManagerCombinedFitWplusJets,
        tauIdEffZtoMuTauHistManagerCombinedFitQCD,
        dataBinnerTauIdEffZtoMuTauCombinedFit,
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinner2d'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_comb2d,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResults2d')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManager2d'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_comb2d,
            histManagers = cms.VPSet(
                muonHistManagerTauIdEffZtoMuTauCombinedFit,
                tauHistManagerTauIdEffZtoMuTauCombinedFit,
                diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit,
                tauIdEffZtoMuTauHistManagerCombinedFit
            ),
            dqmDirectory_store = cms.string('tauIdEffHistograms2d')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinnerComb2dWplusJets'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_comb2dWplusJets,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResultsComb2dWplusJets')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManagerComb2dWplusJets'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_comb2dWplusJets,
            histManagers = cms.VPSet(
                muonHistManagerTauIdEffZtoMuTauCombinedFit,
                tauHistManagerTauIdEffZtoMuTauCombinedFit,
                diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets,
                tauIdEffZtoMuTauHistManagerCombinedFitWplusJets
            ),
            dqmDirectory_store = cms.string('tauIdEffHistogramsComb2dWplusJets')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinnerComb3dQCD'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_comb3dQCD,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResultsComb3dQCD')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManagerComb3dQCD'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_comb3dQCD,
            histManagers = cms.VPSet(
                muonHistManagerTauIdEffZtoMuTauCombinedFit,
                tauHistManagerTauIdEffZtoMuTauCombinedFit,
                diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD,
                tauIdEffZtoMuTauHistManagerCombinedFitQCD
            ),
            dqmDirectory_store = cms.string('tauIdEffHistogramsComb3dQCD')
        )
    ),

    eventDumps = cms.VPSet(),

    analysisSequence = cms.VPSet(
    
        # generator level phase-space selection
        # (NOTE: (1) to be used in case of Monte Carlo samples
        #            overlapping in simulated phase-space only !!
        #        (2) genPhaseSpaceCut needs to be **always** the first entry in the list of cuts
        #           - otherwise the script submitToBatch.csh for submission of cmsRun jobs
        #            to the CERN batch system will not work !!)
        cms.PSet(
            filter = cms.string('genPhaseSpaceCut'),
            title = cms.string('gen. Phase-Space')
        ),
        cms.PSet(
            filter = cms.string('evtSelTrigger'),
            title = cms.string('Muon Trigger')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertex'),
            title = cms.string('Vertex')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertexQuality'),
            title = cms.string('p(chi2Vertex) > 0.01')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertexPosition'),
            title = cms.string('-25 < zVertex < +25 cm')
        ),
        cms.PSet(
            filter = cms.string('evtSelGlobalMuon'),
            title = cms.string('global Muon')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonEta'),
            title = cms.string('-2.1 < eta(Muon) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonPt'),
            title = cms.string('Pt(Muon) > 15 GeV')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauAntiOverlapWithMuonsVeto'),
            title = cms.string('Tau not overlapping w. Muon'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauEta'),
            title = cms.string('-2.1 < eta(Tau) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauPt'),
            title = cms.string('Pt(Tau) > 20 GeV')
        ),        
        cms.PSet(
            filter = cms.string('muonCombRelIsoTightCutTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('Muon tight Track + ECAL relative iso.')
        ),
        cms.PSet(
            filter = cms.string('muonAntiPionCutTauIdEffZtoMuTauCombinedFitTightIso'),
            title = cms.string('Muon pi-Veto')
        ),
        cms.PSet(
            filter = cms.string('muonTrkIPcutTauIdEffZtoMuTauCombinedFitTightIso'),
            title = cms.string('Muon Track IP')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrk'),
            title = cms.string('Tau lead. Track find.')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrkPt'),
            title = cms.string('Tau lead. Track Pt')
        ),
        cms.PSet(
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairBackToBackTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('dPhi(Muon,Tau) > 160 deg.')
        ),
        cms.PSet(
            filter = cms.string('diMuPairZmumuHypothesisVeto'),
            title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ),        
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('num. global Muons < 2')
        ),
        ##cms.PSet(
        ##    filter = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTauCombinedFit'),
        ##    title = cms.string('num. Tau-Jet Candidates < 2')
        ##),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauCombinedFit',
                'tauHistManagerTauIdEffZtoMuTauCombinedFit',
                'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit',
                'diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit',
                'jetHistManagerTauIdEffZtoMuTauCombinedFit',
                'caloMEtHistManager',
                'pfMEtHistManager',
                'tauIdEffZtoMuTauHistManagerCombinedFit',
                'dataBinnerTauIdEffZtoMuTauCombinedFit',
                'tauIdEffDataBinner2d',
                'tauIdEffBinGridHistManager2d'
            ),
            replace = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauCombinedFit.muonSource = ' \
               + 'muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIsoCumulative',
                'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit.diTauCandidateSource = ' \
               + 'muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackCumulative',
                'diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit.diTauCandidateSource = ' \
               + 'muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackCumulative'
            )
        )
    )
)

analyzeEventsTauIdEffZtoMuTauCombinedFitWplusJets = analyzeEventsTauIdEffZtoMuTauCombinedFit.clone(
  
    name = cms.string('TauIdEffAnalyzerZtoMuTauCombinedFitWplusJets'), 

    analysisSequence = cms.VPSet(
    
        # generator level phase-space selection
        # (NOTE: (1) to be used in case of Monte Carlo samples
        #            overlapping in simulated phase-space only !!
        #        (2) genPhaseSpaceCut needs to be **always** the first entry in the list of cuts
        #           - otherwise the script submitToBatch.csh for submission of cmsRun jobs
        #            to the CERN batch system will not work !!)
        cms.PSet(
            filter = cms.string('genPhaseSpaceCut'),
            title = cms.string('gen. Phase-Space')
        ),
        cms.PSet(
            filter = cms.string('evtSelTrigger'),
            title = cms.string('mu15 || isoMu11 Trigger')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertex'),
            title = cms.string('Vertex')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertexQuality'),
            title = cms.string('p(chi2Vertex) > 0.01')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertexPosition'),
            title = cms.string('-25 < zVertex < +25 cm')
        ),
        cms.PSet(
            filter = cms.string('evtSelGlobalMuon'),
            title = cms.string('global Muon')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonEta'),
            title = cms.string('-2.1 < eta(Muon) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonPt'),
            title = cms.string('Pt(Muon) > 15 GeV')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauAntiOverlapWithMuonsVeto'),
            title = cms.string('Tau not overlapping w. Muon'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauEta'),
            title = cms.string('-2.1 < eta(Tau) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauPt'),
            title = cms.string('Pt(Tau) > 20 GeV')
        ),        
        cms.PSet(
            filter = cms.string('muonCombRelIsoTightCutTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('Muon tight Track + ECAL relative iso.')
        ),
        cms.PSet(
            filter = cms.string('muonAntiPionCutTauIdEffZtoMuTauCombinedFitTightIso'),
            title = cms.string('Muon pi-Veto')
        ),
        cms.PSet(
            filter = cms.string('muonTrkIPcutTauIdEffZtoMuTauCombinedFitTightIso'),
            title = cms.string('Muon Track IP')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrk'),
            title = cms.string('Tau lead. Track find.')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrkPt'),
            title = cms.string('Tau lead. Track Pt')
        ),
        cms.PSet(
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTauCombinedFitWplusJets'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairMt1MEtTauIdEffZtoMuTauCombinedFitWplusJets'),
            title = cms.string('M_{T}(Muon-MET) > 50 GeV')
        ),
        cms.PSet(
            filter = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauCombinedFitWplusJets'),
            title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} < -20 GeV')
        ),
        cms.PSet(
            filter = cms.string('diMuPairZmumuHypothesisVeto'),
            title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ),        
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('num. global Muons < 2')
        ),
        ##cms.PSet(
        ##    filter = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTauCombinedFit'),
        ##    title = cms.string('num. Tau-Jet Candidates < 2')
        ##),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauCombinedFit',
                'tauHistManagerTauIdEffZtoMuTauCombinedFit',
                'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets',
                'tauIdEffZtoMuTauHistManagerCombinedFitWplusJets',
                'dataBinnerTauIdEffZtoMuTauCombinedFit',
                'tauIdEffDataBinnerComb2dWplusJets',
                'tauIdEffBinGridHistManagerComb2dWplusJets'
            ),
            replace = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauCombinedFit.muonSource = ' \
               + 'muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIsoCumulative',
                'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets.diTauCandidateSource = ' \
               + 'muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffWplusJetsCumulative'
            )
        )
    )
)

analyzeEventsTauIdEffZtoMuTauCombinedFitQCD = analyzeEventsTauIdEffZtoMuTauCombinedFit.clone(
  
    name = cms.string('TauIdEffAnalyzerZtoMuTauCombinedFitQCD'), 

   analysisSequence = cms.VPSet(
    
        # generator level phase-space selection
        # (NOTE: (1) to be used in case of Monte Carlo samples
        #            overlapping in simulated phase-space only !!
        #        (2) genPhaseSpaceCut needs to be **always** the first entry in the list of cuts
        #           - otherwise the script submitToBatch.csh for submission of cmsRun jobs
        #            to the CERN batch system will not work !!)
        cms.PSet(
            filter = cms.string('genPhaseSpaceCut'),
            title = cms.string('gen. Phase-Space')
        ),
        cms.PSet(
            filter = cms.string('evtSelTrigger'),
            title = cms.string('mu15 || isoMu11 Trigger')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertex'),
            title = cms.string('Vertex')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertexQuality'),
            title = cms.string('p(chi2Vertex) > 0.01')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertexPosition'),
            title = cms.string('-25 < zVertex < +25 cm')
        ),
        cms.PSet(
            filter = cms.string('evtSelGlobalMuon'),
            title = cms.string('global Muon')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonEta'),
            title = cms.string('-2.1 < eta(Muon) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonPt'),
            title = cms.string('Pt(Muon) > 15 GeV')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauAntiOverlapWithMuonsVeto'),
            title = cms.string('Tau not overlapping w. Muon'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauEta'),
            title = cms.string('-2.1 < eta(Tau) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauPt'),
            title = cms.string('Pt(Tau) > 20 GeV')
        ),        
        cms.PSet(
            filter = cms.string('muonCombRelIsoLooseCutTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('Muon loose Track + ECAL relative iso.')
        ),
        cms.PSet(
            filter = cms.string('muonAntiPionCutTauIdEffZtoMuTauCombinedFitLooseIso'),
            title = cms.string('Muon pi-Veto')
        ),
        cms.PSet(
            filter = cms.string('muonTrkIPcutTauIdEffZtoMuTauCombinedFitLooseIso'),
            title = cms.string('Muon Track IP')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrk'),
            title = cms.string('Tau lead. Track find.')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrkPt'),
            title = cms.string('Tau lead. Track Pt')
        ),
        cms.PSet(
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTauCombinedFitQCD'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairBackToBackTauIdEffZtoMuTauCombinedFitQCD'),
            title = cms.string('dPhi(Muon,Tau) > 160 deg.')
        ),
        cms.PSet(
            filter = cms.string('muTauPairMt1MEtTauIdEffZtoMuTauCombinedFitQCD'),
            title = cms.string('M_{T}(Muon-MET) < 30 GeV')
        ),
        cms.PSet(
            filter = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauCombinedFitQCD'),
            title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV')
        ),
        cms.PSet(
            filter = cms.string('diMuPairZmumuHypothesisVeto'),
            title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ),        
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'),
            title = cms.string('num. global Muons < 2')
        ),
        ##cms.PSet(
        ##    filter = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTauCombinedFit'),
        ##    title = cms.string('num. Tau-Jet Candidates < 2')
        ##),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauCombinedFit',
                'tauHistManagerTauIdEffZtoMuTauCombinedFit',
                'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD',
                'tauIdEffZtoMuTauHistManagerCombinedFitQCD',
                'dataBinnerTauIdEffZtoMuTauCombinedFit',
                'tauIdEffDataBinnerComb3dQCD',
                'tauIdEffBinGridHistManagerComb3dQCD'
            ),
            replace = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauCombinedFit.muonSource = ' \
               + 'muonsForTauIdEffZtoMuTauCombinedFitTrkIPlooseIsoCumulative',
                'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD.diTauCandidateSource = ' \
               + 'muTauPairsTauIdEffZtoMuTauCombinedFitPzetaDiffQCDcumulative'
            )
        )
    )
)

#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstTauIdEffZtoMuTauCombinedFitAnalysisSequence = cms.Sequence(
    selectMuonsForTauIdEffZtoMuTauCombinedFit 
   + selectTausForTauIdEffZtoMuTauCombinedFit
   + produceDiMuPairsTauIdEffZtoMuTauCombinedFit
   + produceMuTauPairsTauIdEffZtoMuTauCombinedFit
   + selectJetsForTauIdEffZtoMuTauCombinedFit
   + analyzeEventsTauIdEffZtoMuTauCombinedFit
   + analyzeEventsTauIdEffZtoMuTauCombinedFitWplusJets + analyzeEventsTauIdEffZtoMuTauCombinedFitQCD
)
