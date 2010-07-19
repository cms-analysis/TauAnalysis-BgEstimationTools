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
    src = "selectedPatMuonsPt15Cumulative",
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
    src = "selectedPatMuonsPt15Cumulative",
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

from TauAnalysis.CandidateTools.resolutions_cfi import *

muTauPairsTauIdEffZtoMuTauCombinedFitLooseMuonIso = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsForTauIdEffZtoMuTauCombinedFitTrkIPlooseIsoCumulative'),
    srcLeg2 = cms.InputTag('tausForTauIdEffZtoMuTauCombinedFitMuonVetoCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('patMETs'),
    recoMode = cms.string(""),
    collinearApproxMassCompatibility = cms.PSet(
        mZ = cms.PSet(
            resonanceMass = cms.double(91.2),
            resonanceWidth = cms.double(2.5),
            metResolutionPx = pfMEtResolutionPx,
            metResolutionPy = pfMEtResolutionPy
        )
    ),                                                       
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                        
    verbosity = cms.untracked.int32(0)
)

muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackLooseMuonIso = cms.EDFilter("PATMuTauPairSelector",
    cut = cms.string('dPhi12 > 2.793'),
    filter = cms.bool(False)
)

muTauPairConfiguratorTauIdEffZtoMuTauCombinedFitLooseMuonIso = objSelConfigurator(
    [ muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackLooseMuonIso ],
    src = "muTauPairsTauIdEffZtoMuTauCombinedFitLooseMuonIso",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuTauPairsTauIdEffZtoMuTauCombinedFitLooseMuonIso = muTauPairConfiguratorTauIdEffZtoMuTauCombinedFitLooseMuonIso.configure(pyNameSpace = locals())

muTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso = copy.deepcopy(muTauPairsTauIdEffZtoMuTauCombinedFitLooseMuonIso)
muTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso.srcLeg1 = cms.InputTag('muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIsoCumulative')

muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackTightMuonIso = copy.deepcopy(muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackLooseMuonIso)

muTauPairsTauIdEffZtoMuTauCombinedFitZeroChargeTightMuonIso = cms.EDFilter("PATMuTauPairSelector",
    cut = cms.string('(leg1.charge + leg2.leadTrack.charge) = 0'),
    filter = cms.bool(False)                                                               
)

muTauPairConfiguratorTauIdEffZtoMuTauCombinedFitTightMuonIso = objSelConfigurator(
    [ muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackTightMuonIso,
      muTauPairsTauIdEffZtoMuTauCombinedFitZeroChargeTightMuonIso ],
    src = "muTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso = muTauPairConfiguratorTauIdEffZtoMuTauCombinedFitTightMuonIso.configure(pyNameSpace = locals())

produceMuTauPairsTauIdEffZtoMuTauCombinedFit = cms.Sequence(
    muTauPairsTauIdEffZtoMuTauCombinedFitLooseMuonIso * selectMuTauPairsTauIdEffZtoMuTauCombinedFitLooseMuonIso
   * muTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso * selectMuTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso
)    

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
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit.pluginName = 'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit'
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit.diTauCandidateSource = 'muTauPairsTauIdEffZtoMuTauCombinedFitZeroChargeTightMuonIsoCumulative'
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit.visMassHypothesisSource = cms.InputTag('')
from TauAnalysis.Core.diTauCandidateEventActivityHistManager_cfi import *
diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit = copy.deepcopy(diTauCandidateEventActivityHistManager)
diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit.pluginName = 'diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit'
diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit.pluginType = 'PATMuTauPairEventActivityHistManager'
diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit.diTauCandidateSource = 'muTauPairsTauIdEffZtoMuTauCombinedFitZeroChargeTightMuonIsoCumulative'

diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets.pluginName = 'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets'
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets.diTauCandidateSource = 'muTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso'

diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD.pluginName = 'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD'
diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD.diTauCandidateSource = 'muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackLooseMuonIsoCumulative'

from TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauHistManager_cfi import *
tauIdEffZtoMuTauHistManagerCombinedFit = copy.deepcopy(tauIdEffZtoMuTauHistManager)
tauIdEffZtoMuTauHistManagerCombinedFit.pluginName = 'tauIdEffZtoMuTauHistManagerCombinedFit'
tauIdEffZtoMuTauHistManagerCombinedFit.muonSource = 'muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIsoCumulative'
tauIdEffZtoMuTauHistManagerCombinedFit.tauSource = 'tausForTauIdEffZtoMuTauCombinedFitMuonVetoCumulative'
tauIdEffZtoMuTauHistManagerCombinedFit.diTauSource = 'muTauPairsTauIdEffZtoMuTauCombinedFitZeroChargeTightMuonIsoCumulative'
tauIdEffZtoMuTauHistManagerCombinedFit.diTauChargeSignExtractor.src = tauIdEffZtoMuTauHistManagerCombinedFit.diTauSource

tauIdEffZtoMuTauHistManagerCombinedFitWplusJets = copy.deepcopy(tauIdEffZtoMuTauHistManagerCombinedFit)
tauIdEffZtoMuTauHistManagerCombinedFitWplusJets.pluginName = 'tauIdEffZtoMuTauHistManagerCombinedFitWplusJets'
tauIdEffZtoMuTauHistManagerCombinedFitWplusJets.diTauSource = 'muTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso'
tauIdEffZtoMuTauHistManagerCombinedFitWplusJets.diTauChargeSignExtractor.src = tauIdEffZtoMuTauHistManagerCombinedFitWplusJets.diTauSource

tauIdEffZtoMuTauHistManagerCombinedFitQCD = copy.deepcopy(tauIdEffZtoMuTauHistManagerCombinedFit)
tauIdEffZtoMuTauHistManagerCombinedFitQCD.pluginName ='tauIdEffZtoMuTauHistManagerCombinedFitQCD'
tauIdEffZtoMuTauHistManagerCombinedFitQCD.diTauSource = 'muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackLooseMuonIsoCumulative'
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

binningTauIdEffZtoMuTauCombinedFit_diTauDPhi12 = copy.deepcopy(binning_diTauDPhi12)
binningTauIdEffZtoMuTauCombinedFit_diTauDPhi12.extractor.src = 'muTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso'

binningTauIdEffZtoMuTauCombinedFit_diTauAbsChargeLooseMuonIso = copy.deepcopy(binning_diTauAbsCharge)
binningTauIdEffZtoMuTauCombinedFit_diTauAbsChargeLooseMuonIso.extractor.src = 'muTauPairsTauIdEffZtoMuTauCombinedFitLooseMuonIso'

binningTauIdEffZtoMuTauCombinedFit_diTauAbsChargeTightMuonIso = copy.deepcopy(binning_diTauAbsCharge)
binningTauIdEffZtoMuTauCombinedFit_diTauAbsChargeTightMuonIso.extractor.src = 'muTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso'

tauIdEffBinningZtoMuTau_comb1d = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_genMatrix1d"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauCombinedFit_ewkTauId
    )
)

tauIdEffBinningZtoMuTau_comb3dWplusJets = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_comb3dWplusJets"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauCombinedFit_ewkTauId,
        binningTauIdEffZtoMuTauCombinedFit_diTauDPhi12,
        binningTauIdEffZtoMuTauCombinedFit_diTauAbsChargeTightMuonIso
    )
)

tauIdEffBinningZtoMuTau_comb3dQCD = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_comb3dQCD"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauCombinedFit_ewkTauId,
        binningTauIdEffZtoMuTauCombinedFit_relMuonIso,
        binningTauIdEffZtoMuTauCombinedFit_diTauAbsChargeLooseMuonIso
    )
)

analyzeEventsTauIdEffZtoMuTauCombinedFit = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('TauIdEffAnalyzerZtoMuTauCombinedFit'), 
                            
    filters = cms.VPSet(
        genPhaseSpaceCut,
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
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauCombinedFitLooseMuonIso'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitLooseMuonIso'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairBackToBackTauIdEffZtoMuTauCombinedFitLooseMuonIso'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackLooseMuonIsoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauCombinedFitTightMuonIso'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairBackToBackTauIdEffZtoMuTauCombinedFitTightMuonIso'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackTightMuonIsoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairZeroChargeTauIdEffZtoMuTauCombinedFitTightMuonIso'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauCombinedFitZeroChargeTightMuonIsoCumulative'),
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
        caloMEtHistManager,
        pfMEtHistManager,
        tauIdEffZtoMuTauHistManagerCombinedFit,
        tauIdEffZtoMuTauHistManagerCombinedFitWplusJets,
        tauIdEffZtoMuTauHistManagerCombinedFitQCD,
        dataBinnerTauIdEffZtoMuTauCombinedFit,
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinner2regions'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_comb1d,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResults2regions')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManager2regions'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_comb1d,
            histManagers = cms.VPSet(
                muonHistManagerTauIdEffZtoMuTauCombinedFit,
                tauHistManagerTauIdEffZtoMuTauCombinedFit,
                diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit,
                tauIdEffZtoMuTauHistManagerCombinedFit
            ),
            dqmDirectory_store = cms.string('tauIdEffHistograms2regions')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinnerComb3dWplusJets'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_comb3dWplusJets,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResultsComb3dWplusJets')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManagerComb3dWplusJets'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_comb3dWplusJets,
            histManagers = cms.VPSet(
                diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets,
                tauIdEffZtoMuTauHistManagerCombinedFitWplusJets
            ),
            dqmDirectory_store = cms.string('tauIdEffHistogramsComb3dWplusJets')
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
            filter = cms.string('muTauPairTauIdEffZtoMuTauCombinedFitTightMuonIso'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairBackToBackTauIdEffZtoMuTauCombinedFitTightMuonIso'),
            title = cms.string('dPhi(Muon,Tau) > 160 deg.')
        ),
        cms.PSet(
            filter = cms.string('muTauPairZeroChargeTauIdEffZtoMuTauCombinedFitTightMuonIso'),
            title = cms.string('Charge(Muon+Tau) = 0'),
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
                'caloMEtHistManager',
                'pfMEtHistManager',
                'tauIdEffZtoMuTauHistManagerCombinedFit',
                'dataBinnerTauIdEffZtoMuTauCombinedFit',
                'tauIdEffDataBinner2regions',
                'tauIdEffBinGridHistManager2regions'
            ),
            replace = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauCombinedFit.muonSource = muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIsoCumulative',
                'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFit.diTauCandidateSource = muTauPairsTauIdEffZtoMuTauCombinedFitZeroChargeTightMuonIsoCumulative',
                'diTauCandidateEventActivityHistManagerTauIdEffZtoMuTauCombinedFit.diTauCandidateSource = muTauPairsTauIdEffZtoMuTauCombinedFitZeroChargeTightMuonIsoCumulative'
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
            filter = cms.string('muTauPairTauIdEffZtoMuTauCombinedFitTightMuonIso'),
            title = cms.string('dR(Muon-Tau) > 0.7')
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
                'tauIdEffDataBinnerComb3dWplusJets',
                'tauIdEffBinGridHistManagerComb3dWplusJets'
            ),
            replace = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauCombinedFit.muonSource = muonsForTauIdEffZtoMuTauCombinedFitTrkIPtightIsoCumulative',
                'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitWplusJets.diTauCandidateSource = muTauPairsTauIdEffZtoMuTauCombinedFitTightMuonIso'
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
            filter = cms.string('muTauPairTauIdEffZtoMuTauCombinedFitLooseMuonIso'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairBackToBackTauIdEffZtoMuTauCombinedFitLooseMuonIso'),
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
                'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD',
                'tauIdEffZtoMuTauHistManagerCombinedFitQCD',
                'dataBinnerTauIdEffZtoMuTauCombinedFit',
                'tauIdEffDataBinnerComb3dQCD',
                'tauIdEffBinGridHistManagerComb3dQCD'
            ),
            replace = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauCombinedFit.muonSource = muonsForTauIdEffZtoMuTauCombinedFitTrkIPlooseIsoCumulative',
                'diTauCandidateHistManagerTauIdEffZtoMuTauCombinedFitQCD.diTauCandidateSource = muTauPairsTauIdEffZtoMuTauCombinedFitBackToBackLooseMuonIsoCumulative'
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
   + analyzeEventsTauIdEffZtoMuTauCombinedFit
   + analyzeEventsTauIdEffZtoMuTauCombinedFitWplusJets + analyzeEventsTauIdEffZtoMuTauCombinedFitQCD
)
