import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

#--------------------------------------------------------------------------------
# select loosely selected Z --> tau+ tau- --> muon + tau-jet event sample
# for measuring tau id. efficiency
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patPFTauSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *

tausForTauIdEffZtoMuTauMuonVeto = copy.deepcopy(selectedLayer1TausMuonVeto)

tauSelConfiguratorTauIdEffZtoMuTau = objSelConfigurator(
    [ tausForTauIdEffZtoMuTauMuonVeto ],
    src = "selectedLayer1TausPt20Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectTausTauIdEffZtoMuTau = tauSelConfiguratorTauIdEffZtoMuTau.configure(pyNameSpace = locals())

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

selectMuTauPairsTauIdEffZtoMuTau = cms.Sequence(muTauPairsTauIdEffZtoMuTau + muTauPairsTauIdEffZtoMuTauBackToBack)

#--------------------------------------------------------------------------------  
# produce boolean event selection flags
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.selectZtoMuTau_cff import *

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
    [ cfgTauEtaCutTauIdEffZtoMuTau,
      cfgTauPtCutTauIdEffZtoMuTau,
      cfgTauMuonVetoTauIdEffZtoMuTau,
      cfgMuTauPairTauIdEffZtoMuTau,
      cfgMuTauPairBackToBackTauIdEffZtoMuTau,
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
diTauCandidateHistManagerTauIdEffZtoMuTau.diTauCandidateSource = cms.InputTag('muTauPairsTauIdEffZtoMuTauBackToBack')
diTauCandidateHistManagerTauIdEffZtoMuTau.visMassHypothesisSource = cms.InputTag('')

from TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauHistManager_cfi import *

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
  
    name = cms.string('TauIdEffAnalyzerZtoMuTau'), 
                            
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
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairTauIdEffZtoMuTau')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairBackToBackTauIdEffZtoMuTau'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairBackToBackTauIdEffZtoMuTau')
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
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinner'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinning,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResults')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManager'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinning,
            histManagers = cms.VPSet(
                muonHistManager,
                tauHistManager,
                diTauCandidateHistManagerForMuTau,
                tauIdEffZtoMuTauHistManager
            ),
            dqmDirectory_store = cms.string('tauIdEffHistograms')
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

        # primary event vertex selection
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
            filter = cms.string('muTauPairBackToBackTauIdEffZtoMuTau'),
            title = cms.string('dPhi(Muon,Tau) > 170 deg.'),
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
                'tauIdEffDataBinner',
                'tauIdEffBinGridHistManager'
            )
        )
    )
)

#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstTauIdEffZtoMuTauAnalysisSequence = cms.Sequence(
    selectTausTauIdEffZtoMuTau
   + selectMuTauPairsTauIdEffZtoMuTau
   + selectEventsTauIdEffZtoMuTau
   + analyzeEventsTauIdEffZtoMuTau
)
