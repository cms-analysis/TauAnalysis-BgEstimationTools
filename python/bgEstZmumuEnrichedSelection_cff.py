import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

#--------------------------------------------------------------------------------
# select Z --> mu+ mu- background enriched event sample
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------  
# produce collection of muon + tau-jet combinations
#--------------------------------------------------------------------------------

muTauPairsBgEstZmumuEnriched = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedLayer1MuonsPionVetoCumulative'),
    srcLeg2 = cms.InputTag('selectedLayer1TausChargeCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('layer1METs'),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)
)

muTauPairsBgEstZmumuEnrichedDiscrAgainstMuonsFailed = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsBgEstZmumuEnriched'),                                                   
    cut = cms.string('leg2.tauID("againstMuon") < 0.5'),
    filter = cms.bool(False)
)

selectMuTauPairsBgEstZmumuEnriched = cms.Sequence(muTauPairsBgEstZmumuEnriched + muTauPairsBgEstZmumuEnrichedDiscrAgainstMuonsFailed)

#--------------------------------------------------------------------------------  
# produce boolean event selection flags
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.selectZtoMuTau_cff import *

cfgTauEtaCutBgEstZmumuEnriched = copy.deepcopy(cfgTauEtaCut)
cfgTauEtaCutBgEstZmumuEnriched.pluginName = cms.string('tauEtaCutBgEstZmumuEnriched')
cfgTauEtaCutBgEstZmumuEnriched.src_cumulative = cms.InputTag('selectedLayer1TausEta21Cumulative')
cfgTauEtaCutBgEstZmumuEnriched.systematics = cms.vstring()

cfgTauPtCutBgEstZmumuEnriched = copy.deepcopy(cfgTauEtaCut)
cfgTauPtCutBgEstZmumuEnriched.pluginName = cms.string('tauPtCutBgEstZmumuEnriched')
cfgTauPtCutBgEstZmumuEnriched.src_cumulative = cms.InputTag('selectedLayer1TausPt20Cumulative')
cfgTauPtCutBgEstZmumuEnriched.systematics = cms.vstring()

cfgTauLeadTrkCutBgEstZmumuEnriched = copy.deepcopy(cfgTauLeadTrkCut)
cfgTauLeadTrkCutBgEstZmumuEnriched.pluginName = cms.string('tauLeadTrkCutBgEstZmumuEnriched')
cfgTauLeadTrkCutBgEstZmumuEnriched.src_cumulative = cms.InputTag('selectedLayer1TausLeadTrkCumulative')
cfgTauLeadTrkCutBgEstZmumuEnriched.systematics = cms.vstring()

cfgTauLeadTrkPtCutBgEstZmumuEnriched = copy.deepcopy(cfgTauLeadTrkPtCut)
cfgTauLeadTrkPtCutBgEstZmumuEnriched.pluginName = cms.string('tauLeadTrkPtCutBgEstZmumuEnriched')
cfgTauLeadTrkPtCutBgEstZmumuEnriched.src_cumulative = cms.InputTag('selectedLayer1TausLeadTrkPtCumulative')
cfgTauLeadTrkPtCutBgEstZmumuEnriched.systematics = cms.vstring()

cfgTauTrkIsoCutBgEstZmumuEnriched = copy.deepcopy(cfgTauTrkIsoCut)
cfgTauTrkIsoCutBgEstZmumuEnriched.pluginName = cms.string('tauTrkIsoCutBgEstZmumuEnriched')
cfgTauTrkIsoCutBgEstZmumuEnriched.src_cumulative = cms.InputTag('selectedLayer1TausTrkIsoCumulative')
cfgTauTrkIsoCutBgEstZmumuEnriched.systematics = cms.vstring()

cfgTauEcalIsoCutBgEstZmumuEnriched = copy.deepcopy(cfgTauEcalIsoCut)
cfgTauEcalIsoCutBgEstZmumuEnriched.pluginName = cms.string('tauEcalIsoCutBgEstZmumuEnriched')
cfgTauEcalIsoCutBgEstZmumuEnriched.src_cumulative = cms.InputTag('selectedLayer1TausEcalIsoCumulative')
cfgTauEcalIsoCutBgEstZmumuEnriched.systematics = cms.vstring()

cfgTauProngCutBgEstZmumuEnriched = copy.deepcopy(cfgTauProngCut)
cfgTauProngCutBgEstZmumuEnriched.pluginName = cms.string('tauProngCutBgEstZmumuEnriched')
cfgTauProngCutBgEstZmumuEnriched.src_cumulative = cms.InputTag('selectedLayer1TausProngCumulative')
cfgTauProngCutBgEstZmumuEnriched.systematics = cms.vstring()

cfgTauChargeCutBgEstZmumuEnriched = copy.deepcopy(cfgTauChargeCut)
cfgTauChargeCutBgEstZmumuEnriched.pluginName = cms.string('tauChargeCutBgEstZmumuEnriched')
cfgTauChargeCutBgEstZmumuEnriched.src_cumulative = cms.InputTag('selectedLayer1TausChargeCumulative')
cfgTauChargeCutBgEstZmumuEnriched.systematics = cms.vstring()

cfgMuTauPairBgEstZmumuEnriched = cms.PSet(
    pluginName = cms.string('muTauPairBgEstZmumuEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsBgEstZmumuEnriched'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairDiscrAgainstMuonsFailedBgEstZmumuEnriched = cms.PSet(
    pluginName = cms.string('muTauPairDiscrAgainstMuonsFailedBgEstZmumuEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsBgEstZmumuEnrichedDiscrAgainstMuonsFailed'),
    minNumber = cms.uint32(1)
)

cfgDiMuonCutBgEstZmumuEnriched = cms.PSet(
    pluginName = cms.string('diMuonCutBgEstZmumuEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('selectedLayer1MuonsGlobalCumulative'),
    minNumber = cms.uint32(2)
)

evtSelConfiguratorBgEstZmumuEnriched = eventSelFlagProdConfigurator(
    [ cfgTauEtaCutBgEstZmumuEnriched,
      cfgTauPtCutBgEstZmumuEnriched,
      cfgTauLeadTrkCutBgEstZmumuEnriched,
      cfgTauLeadTrkPtCutBgEstZmumuEnriched,
      cfgTauTrkIsoCutBgEstZmumuEnriched,
      cfgTauEcalIsoCutBgEstZmumuEnriched,
      cfgTauProngCutBgEstZmumuEnriched,
      cfgTauChargeCutBgEstZmumuEnriched,
      cfgMuTauPairBgEstZmumuEnriched,
      cfgMuTauPairDiscrAgainstMuonsFailedBgEstZmumuEnriched,
      cfgDiMuonCutBgEstZmumuEnriched ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

selectEventsBgEstZmumuEnriched = evtSelConfiguratorBgEstZmumuEnriched.configure()

#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *

muonHistManagerBgEstZmumuEnriched = copy.deepcopy(muonHistManager)
muonHistManagerBgEstZmumuEnriched.pluginName = cms.string('muonHistManagerBgEstZmumuEnriched')
muonHistManagerBgEstZmumuEnriched.muonSource = cms.InputTag('selectedLayer1MuonsPionVetoCumulative')

tauHistManagerBgEstZmumuEnriched = copy.deepcopy(tauHistManager)
tauHistManagerBgEstZmumuEnriched.pluginName = cms.string('tauHistManagerBgEstZmumuEnriched')
tauHistManagerBgEstZmumuEnriched.tauSource = cms.InputTag('selectedLayer1TausChargeCumulative')

diTauCandidateHistManagerBgEstZmumuEnriched = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerBgEstZmumuEnriched.pluginName = cms.string('diTauCandidateHistManagerBgEstZmumuEnriched')
diTauCandidateHistManagerBgEstZmumuEnriched.diTauCandidateSource = cms.InputTag('muTauPairsBgEstZmumuEnriched')
diTauCandidateHistManagerBgEstZmumuEnriched.visMassHypothesisSource = cms.InputTag('')

analyzeEventsBgEstZmumuEnriched = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('BgEstTemplateAnalyzer_ZmumuEnriched'), 
                            
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
        cms.PSet(
            pluginName = cms.string('tauEtaCutBgEstZmumuEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauEtaCutBgEstZmumuEnriched', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('tauPtCutBgEstZmumuEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauPtCutBgEstZmumuEnriched', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('tauLeadTrkCutBgEstZmumuEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauLeadTrkCutBgEstZmumuEnriched', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('tauLeadTrkPtCutBgEstZmumuEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauLeadTrkPtCutBgEstZmumuEnriched', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('tauTrkIsoCutBgEstZmumuEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauTrkIsoCutBgEstZmumuEnriched', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('tauEcalIsoCutBgEstZmumuEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauEcalIsoCutBgEstZmumuEnriched', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('tauProngCutBgEstZmumuEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauProngCutBgEstZmumuEnriched', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('tauChargeCutBgEstZmumuEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauChargeCutBgEstZmumuEnriched', 'cumulative'),
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairBgEstZmumuEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairBgEstZmumuEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairDiscrAgainstMuonsFailedOrDiMuonCutBgEstZmumuEnriched'),
            pluginType = cms.string('OrEventSelector'),
            selectors = cms.VPSet(
                cms.PSet(
                     pluginName = cms.string('muTauPairDiscrAgainstMuonsFailedBgEstZmumuEnriched'),
                     pluginType = cms.string('BoolEventSelector'),
                     src = cms.InputTag('muTauPairDiscrAgainstMuonsFailedBgEstZmumuEnriched')
                ),
                cms.PSet(
                     pluginName = cms.string('diMuonCutBgEstZmumuEnriched'),
                     pluginType = cms.string('BoolEventSelector'),
                     src = cms.InputTag('diMuonCutBgEstZmumuEnriched')
                )
            )
        )
    ),
  
    analyzers = cms.VPSet(
        muonHistManagerBgEstZmumuEnriched,
        tauHistManagerBgEstZmumuEnriched,
        diTauCandidateHistManagerBgEstZmumuEnriched
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
            filter = cms.string('tauEtaCutBgEstZmumuEnriched'),
            title = cms.string('-2.1 < eta(Tau) < +2.1'),
        ),
        cms.PSet(
            filter = cms.string('tauPtCutBgEstZmumuEnriched'),
            title = cms.string('Pt(Tau) > 20 GeV'),
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
            filter = cms.string('tauLeadTrkCutBgEstZmumuEnriched'),
            title = cms.string('Tau lead. Track find.'),
        ),
        cms.PSet(
            filter = cms.string('tauLeadTrkPtCutBgEstZmumuEnriched'),
            title = cms.string('Tau lead. Track Pt'),
        ),
        cms.PSet(
            filter = cms.string('tauTrkIsoCutBgEstZmumuEnriched'),
            title = cms.string('Tau Track iso.')
        ),
        cms.PSet(
            filter = cms.string('tauEcalIsoCutBgEstZmumuEnriched'),
            title = cms.string('Tau ECAL iso.')
        ),
        cms.PSet(
             filter = cms.string('tauProngCutBgEstZmumuEnriched'),
             title = cms.string('Tau 1||3-Prong')
        ),
        cms.PSet(
             filter = cms.string('tauChargeCutBgEstZmumuEnriched'),
             title = cms.string('Charge(Tau) = +/-1')
        ),
        cms.PSet(
            filter = cms.string('muTauPairBgEstZmumuEnriched'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairDiscrAgainstMuonsFailedOrDiMuonCutBgEstZmumuEnriched'),
            title = cms.string('Tau mu-Veto failed || di-Muon')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerBgEstZmumuEnriched',
                'tauHistManagerBgEstZmumuEnriched',
                'diTauCandidateHistManagerBgEstZmumuEnriched'
            )
        )
    )
)

#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstZmumuEnrichedAnalysisSequence = cms.Sequence(
    selectMuTauPairsBgEstZmumuEnriched
   + selectEventsBgEstZmumuEnriched
   + analyzeEventsBgEstZmumuEnriched
)
