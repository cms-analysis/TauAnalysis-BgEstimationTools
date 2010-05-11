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
    cut = cms.string('(userIsolation("pat::TrackIso") + userIsolation("pat::EcalIso")) < (0.06*pt)'),
    filter = cms.bool(False)
)

muonsForTauIdEffZtoMuTauPionVetoRelIsolation = copy.deepcopy(selectedLayer1MuonsPionVeto)

muonsForTauIdEffZtoMuTauTrkRelIsolation = copy.deepcopy(selectedLayer1MuonsTrk)

muonsForTauIdEffZtoMuTauTrkIPrelIsolation = copy.deepcopy(selectedLayer1MuonsTrkIP)

muonSelConfiguratorTauIdEffZtoMuTau = objSelConfigurator(
    [ muonsForTauIdEffZtoMuTauCombRelIso,
      muonsForTauIdEffZtoMuTauPionVetoRelIsolation,
      muonsForTauIdEffZtoMuTauTrkRelIsolation,
      muonsForTauIdEffZtoMuTauTrkIPrelIsolation ],
    src = "selectedLayer1MuonsPt15Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuonsForTauIdEffZtoMuTau = muonSelConfiguratorTauIdEffZtoMuTau.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patPFTauSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *

tausForTauIdEffZtoMuTauMuonVeto = copy.deepcopy(selectedLayer1TausMuonVeto)

tauSelConfiguratorTauIdEffZtoMuTau = objSelConfigurator(
    [ tausForTauIdEffZtoMuTauMuonVeto ],
    ##src = "selectedLayer1TausPt20Cumulative",
    src = "selectedLayer1TausLeadTrkPtCumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectTausForTauIdEffZtoMuTau = tauSelConfiguratorTauIdEffZtoMuTau.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of muon + tau-jet combinations
#--------------------------------------------------------------------------------

muTauPairsTauIdEffZtoMuTau = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedLayer1MuonsTrkIPcumulative'),
    srcLeg2 = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('layer1METs'),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                        
    verbosity = cms.untracked.int32(0)
)
#
# Note: do not apply cuts on transverse mass and CDF Pzeta variable
#       in order not to bias muon Pt and eta distributions
#
muTauPairsTauIdEffZtoMuTauBackToBack = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTau'),                                                 
    cut = cms.string('dPhi12 > 2.967'),
    filter = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauNonBackToBack = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTau'),                                                 
    cut = cms.string('dPhi12 < 2.793'),
    filter = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauValidCollinearApprox = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauNonBackToBack'),                                       
    cut = cms.string('collinearApproxIsValid()'),
    filter = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(muTauPairsTauIdEffZtoMuTau)
muTauPairsTauIdEffZtoMuTauRelMuonIsolation.srcLeg1 = cms.InputTag('muonsForTauIdEffZtoMuTauTrkIPrelIsolationCumulative')

muTauPairsTauIdEffZtoMuTauBackToBackRelMuonIsolation = copy.deepcopy(muTauPairsTauIdEffZtoMuTauBackToBack)
muTauPairsTauIdEffZtoMuTauBackToBackRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauRelMuonIsolation')

muTauPairsTauIdEffZtoMuTauNonBackToBackRelMuonIsolation = copy.deepcopy(muTauPairsTauIdEffZtoMuTauNonBackToBack)
muTauPairsTauIdEffZtoMuTauNonBackToBackRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauRelMuonIsolation')

muTauPairsTauIdEffZtoMuTauValidCollinearApproxRelMuonIsolation = copy.deepcopy(muTauPairsTauIdEffZtoMuTauValidCollinearApprox)
muTauPairsTauIdEffZtoMuTauValidCollinearApproxRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauNonBackToBackRelMuonIsolation')

selectMuTauPairsTauIdEffZtoMuTau = cms.Sequence(
    muTauPairsTauIdEffZtoMuTau
   + muTauPairsTauIdEffZtoMuTauBackToBack + muTauPairsTauIdEffZtoMuTauNonBackToBack + muTauPairsTauIdEffZtoMuTauValidCollinearApprox
   + muTauPairsTauIdEffZtoMuTauRelMuonIsolation
   + muTauPairsTauIdEffZtoMuTauBackToBackRelMuonIsolation + muTauPairsTauIdEffZtoMuTauNonBackToBackRelMuonIsolation + muTauPairsTauIdEffZtoMuTauValidCollinearApproxRelMuonIsolation
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
cfgTauEtaCutTauIdEffZtoMuTau.src_cumulative = cms.InputTag('selectedLayer1TausEta21Cumulative')
cfgTauEtaCutTauIdEffZtoMuTau.systematics = cms.vstring()

cfgTauPtCutTauIdEffZtoMuTau = copy.deepcopy(cfgTauEtaCut)
cfgTauPtCutTauIdEffZtoMuTau.pluginName = cms.string('tauPtCutTauIdEffZtoMuTau')
cfgTauPtCutTauIdEffZtoMuTau.src_cumulative = cms.InputTag('selectedLayer1TausPt20Cumulative')
cfgTauPtCutTauIdEffZtoMuTau.systematics = cms.vstring()

cfgTauMuonVetoTauIdEffZtoMuTau = copy.deepcopy(cfgTauMuonVeto)
cfgTauMuonVetoTauIdEffZtoMuTau.pluginName = cms.string('tauMuonVetoTauIdEffZtoMuTau')
cfgTauMuonVetoTauIdEffZtoMuTau.src_cumulative = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative')
cfgTauMuonVetoTauIdEffZtoMuTau.systematics = cms.vstring()

cfgMuTauPairTauIdEffZtoMuTau = cms.PSet(
    pluginName = cms.string('muTauPairTauIdEffZtoMuTau'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTau'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairBackToBackTauIdEffZtoMuTau = cms.PSet(
    pluginName = cms.string('muTauPairBackToBackTauIdEffZtoMuTau'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauBackToBack'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairNonBackToBackTauIdEffZtoMuTau = cms.PSet(
    pluginName = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTau'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauNonBackToBack'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTau = cms.PSet(
    pluginName = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTau'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffZtoMuTauValidCollinearApprox'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairTauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(cfgMuTauPairTauIdEffZtoMuTau)
cfgMuTauPairTauIdEffZtoMuTauRelMuonIsolation.pluginName = cms.string('muTauPairTauIdEffZtoMuTauRelMuonIsolation')
cfgMuTauPairTauIdEffZtoMuTauRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauRelMuonIsolation')

cfgMuTauPairBackToBackTauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(cfgMuTauPairBackToBackTauIdEffZtoMuTau)
cfgMuTauPairBackToBackTauIdEffZtoMuTauRelMuonIsolation.pluginName = cms.string('muTauPairBackToBackTauIdEffZtoMuTauRelMuonIsolation')
cfgMuTauPairBackToBackTauIdEffZtoMuTauRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauBackToBackRelMuonIsolation')

cfgMuTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(cfgMuTauPairNonBackToBackTauIdEffZtoMuTau)
cfgMuTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation.pluginName = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation')
cfgMuTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauNonBackToBackRelMuonIsolation')

cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation = copy.deepcopy(cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTau)
cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation.pluginName = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation')
cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation.src = cms.InputTag('muTauPairsTauIdEffZtoMuTauValidCollinearApproxRelMuonIsolation')

uniqueTauCandidateCutTauIdEffZtoMuTau = cms.EDFilter("BoolEventSelFlagProducer",
    pluginName = cms.string("uniqueTauCandidateCutTauIdEffZtoMuTau"),
    pluginType = cms.string("PATCandViewMaxEventSelector"),
    src = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative'),
    maxNumber = cms.uint32(1)
)

uniqueMuonCandidateCutTauIdEffZtoMuTau = cms.EDFilter("BoolEventSelFlagProducer",
    pluginName = cms.string("uniqueMuonCandidateCutTauIdEffZtoMuTau"),
    pluginType = cms.string("PATCandViewMaxEventSelector"),
    src = cms.InputTag('selectedLayer1MuonsGlobalCumulative'),
    maxNumber = cms.uint32(1)
)

evtSelConfiguratorTauIdEffZtoMuTau = eventSelFlagProdConfigurator(
    [ cfgMuonCombRelIsoCutTauIdEffZtoMuTau,
      cfgMuonAntiPionCutTauIdEffZtoMuTauRelIsolation,
      cfgMuonTrkIPcutTauIdEffZtoMuTauRelIsolation,
      cfgTauEtaCutTauIdEffZtoMuTau,
      cfgTauPtCutTauIdEffZtoMuTau,
      cfgTauMuonVetoTauIdEffZtoMuTau,
      cfgMuTauPairTauIdEffZtoMuTau,
      cfgMuTauPairBackToBackTauIdEffZtoMuTau,
      cfgMuTauPairNonBackToBackTauIdEffZtoMuTau,
      cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTau,
      cfgMuTauPairTauIdEffZtoMuTauRelMuonIsolation,
      cfgMuTauPairBackToBackTauIdEffZtoMuTauRelMuonIsolation,
      cfgMuTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation,
      cfgMuTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation,
      uniqueTauCandidateCutTauIdEffZtoMuTau,
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
muonHistManagerTauIdEffZtoMuTau.muonSource = cms.InputTag('selectedLayer1MuonsTrkIPcumulative')

tauHistManagerTauIdEffZtoMuTau = copy.deepcopy(tauHistManager)
tauHistManagerTauIdEffZtoMuTau.pluginName = cms.string('tauHistManagerTauIdEffZtoMuTau')
tauHistManagerTauIdEffZtoMuTau.tauSource = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative')

diTauCandidateHistManagerTauIdEffZtoMuTau = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTauIdEffZtoMuTau.pluginName = cms.string('diTauCandidateHistManagerTauIdEffZtoMuTau')
##diTauCandidateHistManagerTauIdEffZtoMuTau.diTauCandidateSource = cms.InputTag('muTauPairsTauIdEffZtoMuTauBackToBack')
diTauCandidateHistManagerTauIdEffZtoMuTau.diTauCandidateSource = cms.InputTag('muTauPairsTauIdEffZtoMuTauValidCollinearApprox')
diTauCandidateHistManagerTauIdEffZtoMuTau.visMassHypothesisSource = cms.InputTag('')

from TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauHistManager_cfi import *

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

analyzeEventsTauIdEffZtoMuTau = cms.EDAnalyzer("GenericAnalyzer",
  
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
            src = cms.InputTag('muonCombRelIsoCutTauIdEffZtoMuTau', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('muonAntiPionCutTauIdEffZtoMuTauRelIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonAntiPionCutTauIdEffZtoMuTauRelIsolation', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('muonTrkIPcutTauIdEffZtoMuTauRelIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonTrkIPcutTauIdEffZtoMuTauRelIsolation', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('tauEtaCutTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauEtaCutTauIdEffZtoMuTau', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('tauPtCutTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauPtCutTauIdEffZtoMuTau', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('tauMuonVetoTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauMuonVetoTauIdEffZtoMuTau', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairTauIdEffZtoMuTau')
        ),
        ##cms.PSet(
        ##    pluginName = cms.string('muTauPairBackToBackTauIdEffZtoMuTau'),
        ##    pluginType = cms.string('BoolEventSelector'),
        ##    src = cms.InputTag('muTauPairBackToBackTauIdEffZtoMuTau')
        ##),
        cms.PSet(
            pluginName = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairBackToBackTauIdEffZtoMuTau')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairValidCollinearApproxTauIdEffZtoMuTau')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauRelMuonIsolation'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairTauIdEffZtoMuTauRelMuonIsolation')
        ),
        ##cms.PSet(
        ##    pluginName = cms.string('muTauPairBackToBackTauIdEffZtoMuTauRelMuonIsolation'),
        ##    pluginType = cms.string('BoolEventSelector'),
        ##    src = cms.InputTag('muTauPairBackToBackTauIdEffZtoMuTauRelMuonIsolation')
        ##),
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
        cms.PSet(
            pluginName = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('uniqueTauCandidateCutTauIdEffZtoMuTau')
        ),
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
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTau'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTau'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau',
                'tauHistManagerTauIdEffZtoMuTau',
                'diTauCandidateHistManagerTauIdEffZtoMuTau'
            ),
            replace = cms.vstring(
                'diTauCandidateHistManagerTauIdEffZtoMuTau.diTauCandidateSource = muTauPairsTauIdEffZtoMuTau'
            )
        ),
        ##cms.PSet(
        ##    filter = cms.string('muTauPairBackToBackTauIdEffZtoMuTau'),
        ##    title = cms.string('dPhi(Muon,Tau) > 170 deg.'),
        ##    saveRunEventNumbers = cms.vstring('')
        ##),
        cms.PSet(
            filter = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTau'),
            title = cms.string('dPhi(Muon,Tau) < 160 deg.'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau',
                'tauHistManagerTauIdEffZtoMuTau',
                'diTauCandidateHistManagerTauIdEffZtoMuTau'
            )
        ),
        cms.PSet(
            filter = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTau'),
            title = cms.string('0 < x_1 < 1 && 0 < x_2 < 1'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau',
                'tauHistManagerTauIdEffZtoMuTau',
                'diTauCandidateHistManagerTauIdEffZtoMuTau'
            )
        ),
        cms.PSet(
            filter = cms.string('evtSelDiMuPairZmumuHypothesisVeto'),
            title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau',
                'tauHistManagerTauIdEffZtoMuTau',
                'diTauCandidateHistManagerTauIdEffZtoMuTau'
            )
        ),
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTau'),
            title = cms.string('num. global Muons < 2')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTau',
                'tauHistManagerTauIdEffZtoMuTau',
                'diTauCandidateHistManagerTauIdEffZtoMuTau'
            )
        ),
        cms.PSet(
            filter = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTau'),
            title = cms.string('num. Tau-Jet Candidates < 2'),
            saveRunEventNumbers = cms.vstring('')
        ),
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
            )
        )
    )
)

analyzeEventsTauIdEffZtoMuTauRelMuonIsolation = analyzeEventsTauIdEffZtoMuTau.clone(

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
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTau'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTauRelMuonIsolation'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        ##cms.PSet(
        ##    filter = cms.string('muTauPairBackToBackTauIdEffZtoMuTauRelMuonIsolation'),
        ##    title = cms.string('dPhi(Muon,Tau) > 170 deg.'),
        ##    saveRunEventNumbers = cms.vstring('')
        ##),
        cms.PSet(
            filter = cms.string('muTauPairNonBackToBackTauIdEffZtoMuTauRelMuonIsolation'),
            title = cms.string('dPhi(Muon,Tau) < 160 deg.'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            filter = cms.string('muTauPairValidCollinearApproxTauIdEffZtoMuTauRelMuonIsolation'),
            title = cms.string('0 < x_1 < 1 && 0 < x_2 < 1'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            filter = cms.string('evtSelDiMuPairZmumuHypothesisVeto'),
            title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ),        
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTau'),
            title = cms.string('num. global Muons < 2')
        ),
        cms.PSet(
            filter = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTau'),
            title = cms.string('num. Tau-Jet Candidates < 2'),
            saveRunEventNumbers = cms.vstring('')
        ),
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
   + selectMuTauPairsTauIdEffZtoMuTau
   + selectEventsTauIdEffZtoMuTau
   + analyzeEventsTauIdEffZtoMuTau + analyzeEventsTauIdEffZtoMuTauRelMuonIsolation
)
