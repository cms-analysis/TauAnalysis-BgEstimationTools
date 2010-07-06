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

muonsForTauIdEffZtoMuTauCombRelIso = cms.EDFilter("PATMuonSelector",
    ##cut = cms.string('(userIsolation("pat::TrackIso") + userIsolation("pat::EcalIso")) < (0.06*pt)'),
    #
    # CV: apply loose cut on (relative) muon isolation only,
    #     in order to have still some discrimination between Ztautau/WplusJets and QCD left
    #     for "generalized matrix method"                                                   
    #                                              
    cut = cms.string('(userIsolation("pat::TrackIso") + userIsolation("pat::EcalIso")) < (0.15*pt)'),
    filter = cms.bool(False)
)

muonsForTauIdEffZtoMuTauPionVetoRelIsolation = copy.deepcopy(selectedPatMuonsPionVeto)

muonsForTauIdEffZtoMuTauTrkRelIsolation = copy.deepcopy(selectedPatMuonsTrk)

muonsForTauIdEffZtoMuTauTrkIPrelIsolation = copy.deepcopy(selectedPatMuonsTrkIP)

muonSelConfiguratorTauIdEffZtoMuTau = objSelConfigurator(
    [ muonsForTauIdEffZtoMuTauCombRelIso,
      muonsForTauIdEffZtoMuTauPionVetoRelIsolation,
      muonsForTauIdEffZtoMuTauTrkRelIsolation,
      muonsForTauIdEffZtoMuTauTrkIPrelIsolation ],
    src = "selectedPatMuonsPt15Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuonsForTauIdEffZtoMuTau = muonSelConfiguratorTauIdEffZtoMuTau.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patPFTauSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *

tausForTauIdEffZtoMuTauMuonVeto = copy.deepcopy(selectedPatTausMuonVeto)

tauSelConfiguratorTauIdEffZtoMuTau = objSelConfigurator(
    [ tausForTauIdEffZtoMuTauMuonVeto ],
    src = "selectedPatTausForMuTauLeadTrkPtCumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectTausForTauIdEffZtoMuTau = tauSelConfiguratorTauIdEffZtoMuTau.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of pat::MEt objects
#--------------------------------------------------------------------------------

pfMEtForTauIdEffZtoMuTauPt20 = cms.EDFilter("PATMETSelector",
    src = cms.InputTag("patMETs"),                                 
    cut = cms.string('pt > 20.'),
    filter = cms.bool(False)
)

selectMEtForTauIdEffZtoMuTau = cms.Sequence(pfMEtForTauIdEffZtoMuTauPt20)

#--------------------------------------------------------------------------------  
# produce collection of muon + tau-jet combinations
#--------------------------------------------------------------------------------

from TauAnalysis.CandidateTools.resolutions_cfi import *

muTauPairsTauIdEffZtoMuTauAbsMuonIsolation = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedPatMuonsTrkIPcumulative'),
    srcLeg2 = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative'),
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

muTauPairsTauIdEffZtoMuTauMt1METabsMuonIsolation = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauAbsMuonIsolation'),                             
    cut = cms.string('mt1MET < 50.'),
    filter = cms.bool(False)
)


muTauPairsTauIdEffZtoMuTauPzetaDiffAbsMuonIsolation = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauMt1METabsMuonIsolation'),
    cut = cms.string('(pZeta - 1.5*pZetaVis) > -20.'),
    filter = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauNonBackToBackAbsMuonIsolation = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauPzetaDiffAbsMuonIsolation'),                                                 
    cut = cms.string('dPhi12 < 2.793'),
    ##cut = cms.string('dPhi12 < 3.15'),                                  
    filter = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauValidCollinearApproxAbsMuonIsolation = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauNonBackToBackAbsMuonIsolation'),                                       
    cut = cms.string('collinearApproxIsValid()'),
    filter = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(muTauPairsTauIdEffZtoMuTauAbsMuonIsolation)
muTauPairsTauIdEffZtoMuTauRelMuonIsolation.srcLeg1 = cms.InputTag('muonsForTauIdEffZtoMuTauTrkIPrelIsolationCumulative')

muTauPairsTauIdEffZtoMuTauMt1METrelMuonIsolation = copy.deepcopy(muTauPairsTauIdEffZtoMuTauMt1METabsMuonIsolation)
muTauPairsTauIdEffZtoMuTauMt1METrelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauRelMuonIsolation')

muTauPairsTauIdEffZtoMuTauPzetaDiffRelMuonIsolation = copy.deepcopy(muTauPairsTauIdEffZtoMuTauPzetaDiffAbsMuonIsolation)
muTauPairsTauIdEffZtoMuTauPzetaDiffRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauMt1METrelMuonIsolation')

muTauPairsTauIdEffZtoMuTauNonBackToBackRelMuonIsolation = copy.deepcopy(muTauPairsTauIdEffZtoMuTauNonBackToBackAbsMuonIsolation)
muTauPairsTauIdEffZtoMuTauNonBackToBackRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauPzetaDiffRelMuonIsolation')

muTauPairsTauIdEffZtoMuTauValidCollinearApproxRelMuonIsolation = copy.deepcopy(muTauPairsTauIdEffZtoMuTauValidCollinearApproxAbsMuonIsolation)
muTauPairsTauIdEffZtoMuTauValidCollinearApproxRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauNonBackToBackRelMuonIsolation')

selectMuTauPairsTauIdEffZtoMuTau = cms.Sequence(
    muTauPairsTauIdEffZtoMuTauAbsMuonIsolation
   + muTauPairsTauIdEffZtoMuTauMt1METabsMuonIsolation + muTauPairsTauIdEffZtoMuTauPzetaDiffAbsMuonIsolation
   + muTauPairsTauIdEffZtoMuTauNonBackToBackAbsMuonIsolation + muTauPairsTauIdEffZtoMuTauValidCollinearApproxAbsMuonIsolation
   + muTauPairsTauIdEffZtoMuTauRelMuonIsolation
   + muTauPairsTauIdEffZtoMuTauMt1METrelMuonIsolation + muTauPairsTauIdEffZtoMuTauPzetaDiffRelMuonIsolation 
   + muTauPairsTauIdEffZtoMuTauNonBackToBackRelMuonIsolation + muTauPairsTauIdEffZtoMuTauValidCollinearApproxRelMuonIsolation
)

#--------------------------------------------------------------------------------  
# produce boolean event selection flags
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.selectZtoMuTau_cff import *

cfgMuonCombRelIsoCutTauIdEffZtoMuTau = copy.deepcopy(cfgMuonTrkIsoCut)
cfgMuonCombRelIsoCutTauIdEffZtoMuTau.pluginName = cms.string('muonCombRelIsoCutTauIdEffZtoMuTau')
cfgMuonCombRelIsoCutTauIdEffZtoMuTau.src_cumulative = cms.InputTag('muonsForTauIdEffZtoMuTauCombRelIsoCumulative')
cfgMuonCombRelIsoCutTauIdEffZtoMuTau.systematics = cms.vstring()

cfgMuonAntiPionCutTauIdEffZtoMuTauRelIsolation = copy.deepcopy(cfgMuonAntiPionCut)
cfgMuonAntiPionCutTauIdEffZtoMuTauRelIsolation.pluginName = cms.string('muonAntiPionCutTauIdEffZtoMuTauRelIsolation')
cfgMuonAntiPionCutTauIdEffZtoMuTauRelIsolation.src_cumulative = cms.InputTag('muonsForTauIdEffZtoMuTauPionVetoRelIsolationCumulative')
cfgMuonAntiPionCutTauIdEffZtoMuTauRelIsolation.systematics = cms.vstring()

cfgMuonTrkIPcutTauIdEffZtoMuTauRelIsolation = copy.deepcopy(cfgMuonTrkIPcut)
cfgMuonTrkIPcutTauIdEffZtoMuTauRelIsolation.pluginName = cms.string('muonTrkIPcutTauIdEffZtoMuTauRelIsolation')
cfgMuonTrkIPcutTauIdEffZtoMuTauRelIsolation.src_cumulative = cms.InputTag('muonsForTauIdEffZtoMuTauTrkIPrelIsolationCumulative')
cfgMuonTrkIPcutTauIdEffZtoMuTauRelIsolation.systematics = cms.vstring()

cfgTauEtaCutTauIdEffZtoMuTau = copy.deepcopy(cfgTauEtaCut)
cfgTauEtaCutTauIdEffZtoMuTau.pluginName = cms.string('tauEtaCutTauIdEffZtoMuTau')
cfgTauEtaCutTauIdEffZtoMuTau.src_cumulative = cms.InputTag('selectedPatTausForMuTauEta21Cumulative')
cfgTauEtaCutTauIdEffZtoMuTau.systematics = cms.vstring()

cfgTauPtCutTauIdEffZtoMuTau = copy.deepcopy(cfgTauEtaCut)
cfgTauPtCutTauIdEffZtoMuTau.pluginName = cms.string('tauPtCutTauIdEffZtoMuTau')
cfgTauPtCutTauIdEffZtoMuTau.src_cumulative = cms.InputTag('selectedPatTausForMuTauPt20Cumulative')
cfgTauPtCutTauIdEffZtoMuTau.systematics = cms.vstring()

cfgTauLeadTrkCutTauIdEffZtoMuTau = copy.deepcopy(cfgTauLeadTrkCut)
cfgTauLeadTrkCutTauIdEffZtoMuTau.pluginName = cms.string('tauLeadTrkCutTauIdEffZtoMuTau')
cfgTauLeadTrkCutTauIdEffZtoMuTau.src_cumulative = cms.InputTag('selectedPatTausForMuTauLeadTrkCumulative')
cfgTauLeadTrkCutTauIdEffZtoMuTau.systematics = cms.vstring()

cfgTauLeadTrkPtCutTauIdEffZtoMuTau = copy.deepcopy(cfgTauLeadTrkPtCut)
cfgTauLeadTrkPtCutTauIdEffZtoMuTau.pluginName = cms.string('tauLeadTrkPtCutTauIdEffZtoMuTau')
cfgTauLeadTrkPtCutTauIdEffZtoMuTau.src_cumulative = cms.InputTag('selectedPatTausForMuTauLeadTrkPtCumulative')
cfgTauLeadTrkPtCutTauIdEffZtoMuTau.systematics = cms.vstring()

cfgTauMuonVetoTauIdEffZtoMuTau = copy.deepcopy(cfgTauMuonVeto)
cfgTauMuonVetoTauIdEffZtoMuTau.pluginName = cms.string('tauMuonVetoTauIdEffZtoMuTau')
cfgTauMuonVetoTauIdEffZtoMuTau.src_cumulative = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative')
cfgTauMuonVetoTauIdEffZtoMuTau.systematics = cms.vstring()

cfgPFMEtPt20TauIdEffZtoMuTau = cms.PSet(
    pluginName = cms.string('pfMEtPt20CutTauIdEffZtoMuTau'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('pfMEtForTauIdEffZtoMuTauPt20'),
    systematics = cms.vstring(),
    minNumber = cms.uint32(1)
)

cfgMuTauPairTauIdEffZtoMuTauAbsMuonIsolation = cms.PSet(
    pluginName = cms.string('muTauPairTauIdEffZtoMuTauAbsMuonIsolation'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauAbsMuonIsolation'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairMt1METtauIdEffZtoMuTauAbsMuonIsolation = cms.PSet(
    pluginName = cms.string('muTauPairMt1METtauIdEffZtoMuTauAbsMuonIsolation'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauMt1METabsMuonIsolation'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairPzetaDiffTauIdEffZtoMuTauAbsMuonIsolation = cms.PSet(
    pluginName = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauAbsMuonIsolation'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauPzetaDiffAbsMuonIsolation'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairNonBackToBackTauIdEffZtoMuTauAbsMuonIsolation = cms.PSet(
    pluginName = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTauAbsMuonIsolation'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauNonBackToBackAbsMuonIsolation'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauAbsMuonIsolation = cms.PSet(
    pluginName = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTauAbsMuonIsolation'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauValidCollinearApproxAbsMuonIsolation'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairTauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(cfgMuTauPairTauIdEffZtoMuTauAbsMuonIsolation)
cfgMuTauPairTauIdEffZtoMuTauRelMuonIsolation.pluginName = cms.string('muTauPairTauIdEffZtoMuTauRelMuonIsolation')
cfgMuTauPairTauIdEffZtoMuTauRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauRelMuonIsolation')

cfgMuTauPairMt1METtauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(cfgMuTauPairMt1METtauIdEffZtoMuTauAbsMuonIsolation)
cfgMuTauPairMt1METtauIdEffZtoMuTauRelMuonIsolation.pluginName = cms.string('muTauPairMt1METtauIdEffZtoMuTauRelMuonIsolation')
cfgMuTauPairMt1METtauIdEffZtoMuTauRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauMt1METrelMuonIsolation')

cfgMuTauPairPzetaDiffTauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(cfgMuTauPairPzetaDiffTauIdEffZtoMuTauAbsMuonIsolation)
cfgMuTauPairPzetaDiffTauIdEffZtoMuTauRelMuonIsolation.pluginName = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauRelMuonIsolation')
cfgMuTauPairPzetaDiffTauIdEffZtoMuTauRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauPzetaDiffRelMuonIsolation')

cfgMuTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(cfgMuTauPairNonBackToBackTauIdEffZtoMuTauAbsMuonIsolation)
cfgMuTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation.pluginName = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation')
cfgMuTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauNonBackToBackRelMuonIsolation')

cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauAbsMuonIsolation)
cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation.pluginName = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation')
cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauValidCollinearApproxRelMuonIsolation')

##uniqueTauCandidateCutTauIdEffZtoMuTau = cms.EDFilter("BoolEventSelFlagProducer",
##    pluginName = cms.string("uniqueTauCandidateCutTauIdEffZtoMuTau"),
##    pluginType = cms.string("PATCandViewMaxEventSelector"),
##    src = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative'),
##    maxNumber = cms.uint32(1)
##)

uniqueMuonCandidateCutTauIdEffZtoMuTau = cms.EDFilter("BoolEventSelFlagProducer",
    pluginName = cms.string("uniqueMuonCandidateCutTauIdEffZtoMuTau"),
    pluginType = cms.string("PATCandViewMaxEventSelector"),
    src = cms.InputTag('selectedPatMuonsGlobalCumulative'),
    maxNumber = cms.uint32(1)
)

evtSelConfiguratorTauIdEffZtoMuTau = eventSelFlagProdConfigurator(
    [ cfgMuonCombRelIsoCutTauIdEffZtoMuTau,
      cfgMuonAntiPionCutTauIdEffZtoMuTauRelIsolation,
      cfgMuonTrkIPcutTauIdEffZtoMuTauRelIsolation,
      cfgTauEtaCutTauIdEffZtoMuTau,
      cfgTauPtCutTauIdEffZtoMuTau,
      cfgTauLeadTrkCutTauIdEffZtoMuTau,
      cfgTauLeadTrkPtCutTauIdEffZtoMuTau,
      cfgTauMuonVetoTauIdEffZtoMuTau,
      cfgPFMEtPt20TauIdEffZtoMuTau,
      cfgMuTauPairTauIdEffZtoMuTauAbsMuonIsolation,
      cfgMuTauPairMt1METtauIdEffZtoMuTauAbsMuonIsolation,
      cfgMuTauPairPzetaDiffTauIdEffZtoMuTauAbsMuonIsolation,
      cfgMuTauPairNonBackToBackTauIdEffZtoMuTauAbsMuonIsolation,
      cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauAbsMuonIsolation,
      cfgMuTauPairTauIdEffZtoMuTauRelMuonIsolation,
      cfgMuTauPairMt1METtauIdEffZtoMuTauRelMuonIsolation,
      cfgMuTauPairPzetaDiffTauIdEffZtoMuTauRelMuonIsolation,
      cfgMuTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation,
      cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation,
      ##uniqueTauCandidateCutTauIdEffZtoMuTau,
      uniqueMuonCandidateCutTauIdEffZtoMuTau ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

selectEventsTauIdEffZtoMuTau = evtSelConfiguratorTauIdEffZtoMuTau.configure()

#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *

muonHistManagerTauIdEffZtoMuTau = copy.deepcopy(muonHistManager)
muonHistManagerTauIdEffZtoMuTau.pluginName = cms.string('muonHistManagerTauIdEffZtoMuTau')
muonHistManagerTauIdEffZtoMuTau.muonSource = cms.InputTag('selectedPatMuonsTrkIPcumulative')

tauHistManagerTauIdEffZtoMuTau = copy.deepcopy(tauHistManager)
tauHistManagerTauIdEffZtoMuTau.pluginName = cms.string('tauHistManagerTauIdEffZtoMuTau')
tauHistManagerTauIdEffZtoMuTau.tauSource = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative')

diTauCandidateHistManagerTauIdEffZtoMuTau = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTauIdEffZtoMuTau.pluginName = cms.string('diTauCandidateHistManagerTauIdEffZtoMuTau')
diTauCandidateHistManagerTauIdEffZtoMuTau.diTauCandidateSource = cms.InputTag('muTauPairsTauIdEffZtoMuTauValidCollinearApproxAbsMuonIsolation')
diTauCandidateHistManagerTauIdEffZtoMuTau.visMassHypothesisSource = cms.InputTag('')

from TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauHistManager_cfi import *

from TauAnalysis.BgEstimationTools.bgEstBinGridZtoMuTau_cfi import *

dataBinnerTauIdEffZtoMuTau = copy.deepcopy(dataBinner)
dataBinnerTauIdEffZtoMuTau.pluginName = cms.string('dataBinnerTauIdEffZtoMuTau')

tauIdEffBinning = cms.PSet(
    name = cms.string("tauIdEffBinning"),
    config = cms.VPSet(
        cms.PSet(
            extractor = cms.PSet(
                pluginType = cms.string("PATTauValExtractor"),
                src = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative'),
                value = cms.string("tauID('ewkTauId')")
            ),
            branchName = cms.string('ewkTauId'),
            binning = cms.PSet(
                boundaries = cms.vdouble(0.5),
                min = cms.double(-0.01),
                max = cms.double(1.01)
            )
        )
    )
)

analyzeEventsTauIdEffZtoMuTauAbsMuonIsolation = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('TauIdEffAnalyzerZtoMuTau_absMuonIsolation'), 
                            
    filters = cms.VPSet(
        genPhaseSpaceCut,
        evtSelTrigger,
        evtSelPrimaryEventVertex,
        evtSelPrimaryEventVertexQuality,
        evtSelPrimaryEventVertexPosition,
        evtSelGlobalMuon,
        evtSelMuonEta,
        evtSelMuonPt,
        evtSelMuonTrkIso,
        evtSelMuonEcalIso,
        evtSelMuonAntiPion,
        evtSelMuonTrkIP,
        cms.PSet(
            pluginName = cms.string('muonCombRelIsoCutTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonCombRelIsoCutTauIdEffZtoMuTau', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('muonAntiPionCutTauIdEffZtoMuTauRelIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonAntiPionCutTauIdEffZtoMuTauRelIsolation', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('muonTrkIPcutTauIdEffZtoMuTauRelIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonTrkIPcutTauIdEffZtoMuTauRelIsolation', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauEtaCutTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauEtaCutTauIdEffZtoMuTau', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauPtCutTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauPtCutTauIdEffZtoMuTau', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauLeadTrkCutTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauLeadTrkCutTauIdEffZtoMuTau', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauLeadTrkPtCutTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauLeadTrkPtCutTauIdEffZtoMuTau', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauMuonVetoTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauMuonVetoTauIdEffZtoMuTau', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('pfMEtPt20CutTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('pfMEtPt20CutTauIdEffZtoMuTau')
        ),                                                   
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauAbsMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairTauIdEffZtoMuTauAbsMuonIsolation')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairMt1METtauIdEffZtoMuTauAbsMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairMt1METtauIdEffZtoMuTauAbsMuonIsolation')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauAbsMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairPzetaDiffTauIdEffZtoMuTauAbsMuonIsolation')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTauAbsMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairNonBackToBackTauIdEffZtoMuTauAbsMuonIsolation')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTauAbsMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairValidCollinearApproxTauIdEffZtoMuTauAbsMuonIsolation')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauRelMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairTauIdEffZtoMuTauRelMuonIsolation')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairMt1METtauIdEffZtoMuTauRelMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairMt1METtauIdEffZtoMuTauRelMuonIsolation')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauRelMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairPzetaDiffTauIdEffZtoMuTauRelMuonIsolation')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation')
        ),
        evtSelDiMuPairZmumuHypothesisVeto,
        ##cms.PSet(
        ##    pluginName = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTau'),
        ##    pluginType = cms.string('BoolEventSelector'),
        ##    src = cms.InputTag('uniqueTauCandidateCutTauIdEffZtoMuTau')
        ##),
        cms.PSet(
            pluginName = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('uniqueMuonCandidateCutTauIdEffZtoMuTau')
        )
    ),
  
    analyzers = cms.VPSet(
        muonHistManagerTauIdEffZtoMuTau,
        tauHistManagerTauIdEffZtoMuTau,
        diTauCandidateHistManagerTauIdEffZtoMuTau,
        caloMEtHistManager,
        pfMEtHistManager,
        tauIdEffZtoMuTauHistManager,
        dataBinnerTauIdEffZtoMuTau,
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinner2regions'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinning,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResults2regions')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManager2regions'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinning,
            histManagers = cms.VPSet(
                muonHistManagerTauIdEffZtoMuTau,
                tauHistManagerTauIdEffZtoMuTau,
                diTauCandidateHistManagerTauIdEffZtoMuTau,
                tauIdEffZtoMuTauHistManager
            ),
            dqmDirectory_store = cms.string('tauIdEffHistograms2regions')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinnerGenMatrix3d'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_genMatrix3d,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResultsGenMatrix3d')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManagerGenMatrix3d'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_genMatrix3d,
            histManagers = cms.VPSet(
                tauHistManager
            ),
            dqmDirectory_store = cms.string('tauIdEffHistogramsGenMatrix3d')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinnerGenMatrix4d'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_genMatrix4d,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResultsGenMatrix4d')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManagerGenMatrix4d'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_genMatrix4d,
            histManagers = cms.VPSet(
                tauHistManager
            ),
            dqmDirectory_store = cms.string('tauIdEffHistogramsGenMatrix4d')
        ),
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
            filter = cms.string('tauEtaCutTauIdEffZtoMuTau'),
            title = cms.string('-2.1 < eta(Tau) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('tauPtCutTauIdEffZtoMuTau'),
            title = cms.string('Pt(Tau) > 20 GeV')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonTrkIso'),
            title = cms.string('Muon Track iso.')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonEcalIso'),
            title = cms.string('Muon ECAL iso.')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonAntiPion'),
            title = cms.string('Muon pi-Veto')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonTrkIP'),
            title = cms.string('Muon Track IP')
        ),
        cms.PSet(
            filter = cms.string('tauLeadTrkCutTauIdEffZtoMuTau'),
            title = cms.string('Tau lead. Track find.')
        ),
        cms.PSet(
            filter = cms.string('tauLeadTrkPtCutTauIdEffZtoMuTau'),
            title = cms.string('Tau lead. Track Pt')
        ),
        cms.PSet(
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTau'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('pfMEtPt20CutTauIdEffZtoMuTau'),
            title = cms.string('PFMEt Pt > 20 GeV')
        ),                                                            
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTauAbsMuonIsolation'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairMt1METtauIdEffZtoMuTauAbsMuonIsolation'),
            title = cms.string('M_{T}(Muon-MET) < 50 GeV')
        ),
        cms.PSet(
            filter = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauAbsMuonIsolation'),
            title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV')
        ),
        ##cms.PSet(
        ##    filter = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTauAbsMuonIsolation'),
        ##    title = cms.string('dPhi(Muon,Tau) < 160 deg.'),
        ##    saveRunEventNumbers = cms.vstring('passed_cumulative')
        ##),
        cms.PSet(
            filter = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTauAbsMuonIsolation'),
            title = cms.string('0 < x_1 < 1 && 0 < x_2 < 1')
        ),
        cms.PSet(
            filter = cms.string('evtSelDiMuPairZmumuHypothesisVeto'),
            title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ),
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTau'),
            title = cms.string('num. global Muons < 2')
        ),
        ##cms.PSet(
        ##    filter = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTau'),
        ##    title = cms.string('num. Tau-Jet Candidates < 2'),
        ##    saveRunEventNumbers = cms.vstring('')
        ##),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau',
                'tauHistManagerTauIdEffZtoMuTau',
                'diTauCandidateHistManagerTauIdEffZtoMuTau',
                'caloMEtHistManager',
                'pfMEtHistManager',
                'dataBinnerTauIdEffZtoMuTau',
                'tauIdEffZtoMuTauHistManager',
                'tauIdEffDataBinner2regions',
                'tauIdEffBinGridHistManager2regions'
            ),
            replace = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau.muonSource = selectedPatMuonsTrkIPcumulative',
                'diTauCandidateHistManagerTauIdEffZtoMuTau.diTauCandidateSource = muTauPairsTauIdEffZtoMuTauValidCollinearApproxAbsMuonIsolation',
                'pfMEtHistManager.metSource = pfMEtForTauIdEffZtoMuTauPt20'
            )
        )
    )
)

analyzeEventsTauIdEffZtoMuTauRelMuonIsolation = analyzeEventsTauIdEffZtoMuTauAbsMuonIsolation.clone(

    name = cms.string('TauIdEffAnalyzerZtoMuTau_relMuonIsolation'), 
    
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
            filter = cms.string('tauEtaCutTauIdEffZtoMuTau'),
            title = cms.string('-2.1 < eta(Tau) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('tauPtCutTauIdEffZtoMuTau'),
            title = cms.string('Pt(Tau) > 20 GeV')
        ),        
        cms.PSet(
            filter = cms.string('muonCombRelIsoCutTauIdEffZtoMuTau'),
            title = cms.string('Muon Track + ECAL relative iso.')
        ),
        cms.PSet(
            filter = cms.string('muonAntiPionCutTauIdEffZtoMuTauRelIsolation'),
            title = cms.string('Muon pi-Veto')
        ),
        cms.PSet(
            filter = cms.string('muonTrkIPcutTauIdEffZtoMuTauRelIsolation'),
            title = cms.string('Muon Track IP')
        ),
        cms.PSet(
            filter = cms.string('tauLeadTrkCutTauIdEffZtoMuTau'),
            title = cms.string('Tau lead. Track find.')
        ),
        cms.PSet(
            filter = cms.string('tauLeadTrkPtCutTauIdEffZtoMuTau'),
            title = cms.string('Tau lead. Track Pt')
        ),
        cms.PSet(
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTau'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('pfMEtPt20CutTauIdEffZtoMuTau'),
            title = cms.string('PFMEt Pt > 20 GeV')
        ), 
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTauRelMuonIsolation'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairMt1METtauIdEffZtoMuTauRelMuonIsolation'),
            title = cms.string('M_{T}(Muon-MET) < 50 GeV')
        ),
        cms.PSet(
            filter = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauRelMuonIsolation'),
            title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV')
        ),
        ##cms.PSet(
        ##    filter = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation'),
        ##    title = cms.string('dPhi(Muon,Tau) < 160 deg.')
        ##),
        cms.PSet(
            filter = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation'),
            title = cms.string('0 < x_1 < 1 && 0 < x_2 < 1')
        ),
        cms.PSet(
            filter = cms.string('evtSelDiMuPairZmumuHypothesisVeto'),
            title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ),        
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTau'),
            title = cms.string('num. global Muons < 2')
        ),
        ##cms.PSet(
        ##    filter = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTau'),
        ##    title = cms.string('num. Tau-Jet Candidates < 2')
        ##),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau',
                'tauHistManagerTauIdEffZtoMuTau',
                'diTauCandidateHistManagerTauIdEffZtoMuTau',
                'caloMEtHistManager',
                'pfMEtHistManager',
                'tauIdEffZtoMuTauHistManager',
                'dataBinnerTauIdEffZtoMuTau',
                'tauIdEffDataBinner2regions',
                'tauIdEffBinGridHistManager2regions'
            ),
            replace = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau.muonSource = muonsForTauIdEffZtoMuTauTrkIPrelIsolationCumulative',
                'diTauCandidateHistManagerTauIdEffZtoMuTau.diTauCandidateSource = muTauPairsTauIdEffZtoMuTauValidCollinearApproxRelMuonIsolation',
                'pfMEtHistManager.metSource = pfMEtForTauIdEffZtoMuTauPt20'
            )
        )
    )
)

analyzeEventsTauIdEffZtoMuTauGenMatrixRelMuonIsolation = analyzeEventsTauIdEffZtoMuTauRelMuonIsolation.clone(

    name = cms.string('TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation'), 
    
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
            filter = cms.string('tauEtaCutTauIdEffZtoMuTau'),
            title = cms.string('-2.1 < eta(Tau) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('tauPtCutTauIdEffZtoMuTau'),
            title = cms.string('Pt(Tau) > 20 GeV')
        ),        
        cms.PSet(
            filter = cms.string('muonCombRelIsoCutTauIdEffZtoMuTau'),
            title = cms.string('Muon Track + ECAL relative iso.')
        ),
        cms.PSet(
            filter = cms.string('muonAntiPionCutTauIdEffZtoMuTauRelIsolation'),
            title = cms.string('Muon pi-Veto')
        ),
        cms.PSet(
            filter = cms.string('muonTrkIPcutTauIdEffZtoMuTauRelIsolation'),
            title = cms.string('Muon Track IP')
        ),
        cms.PSet(
            filter = cms.string('tauLeadTrkCutTauIdEffZtoMuTau'),
            title = cms.string('Tau lead. Track find.')
        ),
        cms.PSet(
            filter = cms.string('tauLeadTrkPtCutTauIdEffZtoMuTau'),
            title = cms.string('Tau lead. Track Pt')
        ),
        cms.PSet(
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTau'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTauRelMuonIsolation'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTau'),
            title = cms.string('num. global Muons < 2')
        ),        
        ##cms.PSet(
        ##    filter = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTau'),
        ##    title = cms.string('num. Tau-Jet Candidates < 2')
        ##),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau',
                'tauHistManagerTauIdEffZtoMuTau',
                'diTauCandidateHistManagerTauIdEffZtoMuTau',
                'caloMEtHistManager',
                'pfMEtHistManager',
                'tauIdEffZtoMuTauHistManager',
                'dataBinnerTauIdEffZtoMuTau',
                'tauIdEffDataBinner2regions',
                'tauIdEffBinGridHistManager2regions',
                'tauIdEffDataBinnerGenMatrix3d',
                'tauIdEffBinGridHistManagerGenMatrix3d',
                'tauIdEffDataBinnerGenMatrix4d',
                'tauIdEffBinGridHistManagerGenMatrix4d'
            ),
            replace = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau.muonSource = muonsForTauIdEffZtoMuTauTrkIPrelIsolationCumulative',
                'diTauCandidateHistManagerTauIdEffZtoMuTau.diTauCandidateSource = muTauPairsTauIdEffZtoMuTauRelMuonIsolation',
                'pfMEtHistManager.metSource = patMETs'
            )
        )
    )
)    

#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstTauIdEffZtoMuTauAnalysisSequence = cms.Sequence(
    selectMuonsForTauIdEffZtoMuTau
   + selectTausForTauIdEffZtoMuTau
   + selectMEtForTauIdEffZtoMuTau
   + selectMuTauPairsTauIdEffZtoMuTau
   + selectEventsTauIdEffZtoMuTau
   + analyzeEventsTauIdEffZtoMuTauAbsMuonIsolation + analyzeEventsTauIdEffZtoMuTauRelMuonIsolation
   + analyzeEventsTauIdEffZtoMuTauGenMatrixRelMuonIsolation
)
