import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

#--------------------------------------------------------------------------------
# select W + jets background enriched event sample
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------  
# produce collection of pat::Muons
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patMuonSelection_cfi import *
#
# Note: muonPt cut improves WplusJets/QCD ratio by about a factor five,
#       but significantly shifts the muon + tau-jet visible invariant mass distribution towards higher values.
#       In order to supress QCD background contamination (on a statistical basis),
#       could extract W + jets template shape from difference in muon + tau-jet visible invariant mass distributions
#       of opposite sign - same sign muon and tau-jet combinations.
#      (SS/OS ratio is close to one for QCD background; significant charge asymmetry expected for W + jets background)
#    
muonsBgEstWplusJetsEnrichedPt = copy.deepcopy(selectedLayer1MuonsPt15)
muonsBgEstWplusJetsEnrichedPt.cut = cms.string('pt > 15.')

muonsBgEstWplusJetsEnrichedTrkIso = copy.deepcopy(selectedLayer1MuonsTrkIso)
muonsBgEstWplusJetsEnrichedTrkIso.sumPtMax = cms.double(1.)

muonsBgEstWplusJetsEnrichedEcalIso = copy.deepcopy(selectedLayer1MuonsEcalIso)
muonsBgEstWplusJetsEnrichedEcalIso.cut = cms.string('ecalIso < 1.')

muonSelConfiguratorBgEstWplusJetsEnriched = objSelConfigurator(
    [ muonsBgEstWplusJetsEnrichedPt,
      muonsBgEstWplusJetsEnrichedTrkIso,
      muonsBgEstWplusJetsEnrichedEcalIso ],
    src = "selectedLayer1MuonsEta21Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuonsBgEstWplusJetsEnriched = muonSelConfiguratorBgEstWplusJetsEnriched.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patPFTauSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *
#
# Note: probability for quark/gluon jets to pass tau track and ECAL isolation criteria
#       is higher for low Pt than for high Pt jets; the consequence is that muon + tau-jet visible invariant mass distribution
#       gets shifted towards higher values in case tau track and ECAL isolation criteria are not applied.
#       For this reason, either need to apply tau track and ECAL isolation criteria in selection of W + jets background enriched sample
#       or correct for template shape distortion by reweighting
#      (would gain a factor of about 2.5 in event statistics; reweighting of tauPt distribution not implemented yet, however)
#   
tausBgEstWplusJetsEnrichedTrkIso = copy.deepcopy(selectedLayer1TausTrkIso)
#tausBgEstWplusJetsEnrichedTrkIso.cut = cms.string('tauID("trackIsolation") > 0.5 | chargedHadronIso < 8.')
tausBgEstWplusJetsEnrichedTrkIso.cut = cms.string('tauID("trackIsolation") > 0.5')

tausBgEstWplusJetsEnrichedEcalIso = copy.deepcopy(selectedLayer1TausEcalIso)
#tausBgEstWplusJetsEnrichedEcalIso.cut = cms.string('tauID("ecalIsolation") > 0.5 | photonIso < 8.')
tausBgEstWplusJetsEnrichedEcalIso.cut = cms.string('tauID("ecalIsolation") > 0.5')

tausBgEstWplusJetsEnrichedMuonVeto = copy.deepcopy(selectedLayer1TausMuonVeto)

tauSelConfiguratorBgEstWplusJetsEnriched = objSelConfigurator(
    [ tausBgEstWplusJetsEnrichedTrkIso,
      tausBgEstWplusJetsEnrichedEcalIso,
      tausBgEstWplusJetsEnrichedMuonVeto ],
    src = "selectedLayer1TausForMuTauLeadTrkPtCumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectTausBgEstWplusJetsEnriched = tauSelConfiguratorBgEstWplusJetsEnriched.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of muon + tau-jet combinations
#--------------------------------------------------------------------------------

muTauPairsBgEstWplusJetsEnriched = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsBgEstWplusJetsEnrichedEcalIsoCumulative'),
    srcLeg2 = cms.InputTag('tausBgEstWplusJetsEnrichedMuonVetoCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('layer1METs'),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)
)

muTauPairsBgEstWplusJetsEnrichedMt1MET = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsBgEstWplusJetsEnriched'),                                                   
    cut = cms.string('mt1MET > 30.'),
    filter = cms.bool(False)
)

selectMuTauPairsBgEstWplusJetsEnriched = cms.Sequence(muTauPairsBgEstWplusJetsEnriched + muTauPairsBgEstWplusJetsEnrichedMt1MET)

#--------------------------------------------------------------------------------  
# produce collection of pat::Jets used for central jet veto
# (in order to reject QCD di-jet events)
#--------------------------------------------------------------------------------

jetsBgEstWplusJetsEnrichedAntiOverlapWithLeptonsVeto = cms.EDFilter("PATJetAntiOverlapSelector",
    src = cms.InputTag("selectedLayer1JetsEt20Cumulative"),                                                                  
    srcNotToBeFiltered = cms.VInputTag(
        "selectedLayer1ElectronsTrkIPcumulative",
        "muonsBgEstWplusJetsEnrichedEcalIsoCumulative",
        "tausBgEstWplusJetsEnrichedMuonVetoCumulative"
    ),
    dRmin = cms.double(0.7),
    filter = cms.bool(False)                                           
)

jetsBgEstWplusJetsEnrichedAlpha0point1 = cms.EDFilter("PATJetAlphaSelector",
    src = cms.InputTag("jetsBgEstWplusJetsEnrichedAntiOverlapWithLeptonsVeto"),
    alphaMin = cms.double(0.1),
    filter = cms.bool(False)
)

selectJetsBgEstWplusJetsEnriched = cms.Sequence(
    jetsBgEstWplusJetsEnrichedAntiOverlapWithLeptonsVeto + jetsBgEstWplusJetsEnrichedAlpha0point1
)

#--------------------------------------------------------------------------------  
# produce boolean event selection flags
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.selectZtoMuTau_cff import *

cfgMuonPtCutBgEstWplusJetsEnriched = copy.deepcopy(cfgMuonPtCut)
cfgMuonPtCutBgEstWplusJetsEnriched.pluginName = cms.string('muonPtCutBgEstWplusJetsEnriched')
cfgMuonPtCutBgEstWplusJetsEnriched.src_cumulative = cms.InputTag('muonsBgEstWplusJetsEnrichedPtCumulative')
cfgMuonPtCutBgEstWplusJetsEnriched.systematics = cms.vstring()

cfgMuonTrkIsoCutBgEstWplusJetsEnriched = copy.deepcopy(cfgMuonTrkIsoCut)
cfgMuonTrkIsoCutBgEstWplusJetsEnriched.pluginName = cms.string('muonTrkIsoCutBgEstWplusJetsEnriched')
cfgMuonTrkIsoCutBgEstWplusJetsEnriched.src_cumulative = cms.InputTag('muonsBgEstWplusJetsEnrichedTrkIsoCumulative')
cfgMuonTrkIsoCutBgEstWplusJetsEnriched.systematics = cms.vstring()

cfgMuonEcalIsoCutBgEstWplusJetsEnriched = copy.deepcopy(cfgMuonEcalIsoCut)
cfgMuonEcalIsoCutBgEstWplusJetsEnriched.pluginName = cms.string('muonEcalIsoCutBgEstWplusJetsEnriched')
cfgMuonEcalIsoCutBgEstWplusJetsEnriched.src_cumulative = cms.InputTag('muonsBgEstWplusJetsEnrichedEcalIsoCumulative')
cfgMuonEcalIsoCutBgEstWplusJetsEnriched.systematics = cms.vstring()

cfgTauTrkIsoCutBgEstWplusJetsEnriched = copy.deepcopy(cfgTauTrkIsoCut)
cfgTauTrkIsoCutBgEstWplusJetsEnriched.pluginName = cms.string('tauTrkIsoCutBgEstWplusJetsEnriched')
cfgTauTrkIsoCutBgEstWplusJetsEnriched.src_cumulative = cms.InputTag('tausBgEstWplusJetsEnrichedTrkIsoCumulative')
cfgTauTrkIsoCutBgEstWplusJetsEnriched.systematics = cms.vstring()

cfgTauEcalIsoCutBgEstWplusJetsEnriched = copy.deepcopy(cfgTauEcalIsoCut)
cfgTauEcalIsoCutBgEstWplusJetsEnriched.pluginName = cms.string('tauEcalIsoCutBgEstWplusJetsEnriched')
cfgTauEcalIsoCutBgEstWplusJetsEnriched.src_cumulative = cms.InputTag('tausBgEstWplusJetsEnrichedEcalIsoCumulative')
cfgTauEcalIsoCutBgEstWplusJetsEnriched.systematics = cms.vstring()

cfgTauMuonVetoBgEstWplusJetsEnriched = copy.deepcopy(cfgTauMuonVeto)
cfgTauMuonVetoBgEstWplusJetsEnriched.pluginName = cms.string('tauMuonVetoBgEstWplusJetsEnriched')
cfgTauMuonVetoBgEstWplusJetsEnriched.src_cumulative = cms.InputTag('tausBgEstWplusJetsEnrichedMuonVetoCumulative')
cfgTauMuonVetoBgEstWplusJetsEnriched.systematics = cms.vstring()

cfgMuTauPairBgEstWplusJetsEnriched = cms.PSet(
    pluginName = cms.string('muTauPairBgEstWplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsBgEstWplusJetsEnriched'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairMt1MEtBgEstWplusJetsEnriched= cms.PSet(
    pluginName = cms.string('muTauPairMt1MEtBgEstWplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsBgEstWplusJetsEnrichedMt1MET'),
    minNumber = cms.uint32(1)
)

cfgCentralJetVetoBgEstWplusJetsEnriched = cms.PSet(
    pluginName = cms.string('centralJetVetoBgEstWplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMaxEventSelector'),
    src = cms.InputTag('jetsBgEstWplusJetsEnrichedAlpha0point1'),
    maxNumber = cms.uint32(0)
)

cfgDiMuonVetoBgEstWplusJetsEnriched = cms.PSet(
    pluginName = cms.string('diMuonVetoBgEstWplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMaxEventSelector'),
    src = cms.InputTag('selectedLayer1MuonsGlobalIndividual'),
    maxNumber = cms.uint32(1)
)

evtSelConfiguratorBgEstWplusJetsEnriched = eventSelFlagProdConfigurator(
    [ cfgMuonPtCutBgEstWplusJetsEnriched,
      cfgMuonTrkIsoCutBgEstWplusJetsEnriched,
      cfgMuonEcalIsoCutBgEstWplusJetsEnriched,
      cfgTauTrkIsoCutBgEstWplusJetsEnriched,
      cfgTauEcalIsoCutBgEstWplusJetsEnriched,
      cfgTauMuonVetoBgEstWplusJetsEnriched,
      cfgMuTauPairBgEstWplusJetsEnriched,
      cfgMuTauPairMt1MEtBgEstWplusJetsEnriched,
      cfgCentralJetVetoBgEstWplusJetsEnriched,
      cfgDiMuonVetoBgEstWplusJetsEnriched ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

selectEventsBgEstWplusJetsEnriched = evtSelConfiguratorBgEstWplusJetsEnriched.configure()

#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *

muonHistManagerWplusJetsEnriched = copy.deepcopy(muonHistManager)
muonHistManagerWplusJetsEnriched.pluginName = cms.string('muonHistManagerWplusJetsEnriched')
muonHistManagerWplusJetsEnriched.muonSource = cms.InputTag('muonsBgEstWplusJetsEnrichedEcalIsoCumulative')

tauHistManagerWplusJetsEnriched = copy.deepcopy(tauHistManager)
tauHistManagerWplusJetsEnriched.pluginName = cms.string('tauHistManagerWplusJetsEnriched')
tauHistManagerWplusJetsEnriched.tauSource = cms.InputTag('tausBgEstWplusJetsEnrichedMuonVetoCumulative')

diTauCandidateHistManagerWplusJetsEnriched = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerWplusJetsEnriched.pluginName = cms.string('diTauCandidateHistManagerWplusJetsEnriched')
diTauCandidateHistManagerWplusJetsEnriched.diTauCandidateSource = cms.InputTag('muTauPairsBgEstWplusJetsEnrichedMt1MET')

analyzeEventsBgEstWplusJetsEnriched = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('BgEstTemplateAnalyzer_WplusJetsEnriched'), 
                            
    filters = cms.VPSet(
        genPhaseSpaceCut,
        evtSelTrigger,
        evtSelPrimaryEventVertex,
        evtSelPrimaryEventVertexQuality,
        evtSelPrimaryEventVertexPosition,
        evtSelGlobalMuon,
        evtSelMuonEta,
        cms.PSet(
            pluginName = cms.string('muonPtCutBgEstWplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonPtCutBgEstWplusJetsEnriched', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('muonTrkIsoCutBgEstWplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonTrkIsoCutBgEstWplusJetsEnriched', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('muonEcalIsoCutBgEstWplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonEcalIsoCutBgEstWplusJetsEnriched', 'cumulative')
        ),
        evtSelTauAntiOverlapWithMuonsVeto,
        evtSelTauEta,
        evtSelTauPt,
        evtSelTauLeadTrk,
        evtSelTauLeadTrkPt,
        cms.PSet(
            pluginName = cms.string('tauTrkIsoCutBgEstWplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauTrkIsoCutBgEstWplusJetsEnriched', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauEcalIsoCutBgEstWplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauEcalIsoCutBgEstWplusJetsEnriched', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauMuonVetoBgEstWplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauMuonVetoBgEstWplusJetsEnriched', 'cumulative')
        ),      
        cms.PSet(
            pluginName = cms.string('muTauPairBgEstWplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairBgEstWplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairMt1MEtBgEstWplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairMt1MEtBgEstWplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('centralJetVetoBgEstWplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('centralJetVetoBgEstWplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('diMuonVetoBgEstWplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diMuonVetoBgEstWplusJetsEnriched')
        )
    ),
  
    analyzers = cms.VPSet(
        muonHistManagerWplusJetsEnriched,
        tauHistManagerWplusJetsEnriched,
        diTauCandidateHistManagerWplusJetsEnriched
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
            filter = cms.string('muonPtCutBgEstWplusJetsEnriched'),
            title = cms.string('Pt(Muon) > 15 GeV')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauAntiOverlapWithMuonsVeto'),
            title = cms.string('Tau not overlapping w. Muon'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauEta'),
            title = cms.string('-2.1 < eta(Tau) < +2.1'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauPt'),
            title = cms.string('Pt(Tau) > 20 GeV'),
        ),
        cms.PSet(
            filter = cms.string('muonTrkIsoCutBgEstWplusJetsEnriched'),
            title = cms.string('Muon Track iso.')
        ),
        cms.PSet(
            filter = cms.string('muonEcalIsoCutBgEstWplusJetsEnriched'),
            title = cms.string('Muon ECAL iso.')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrk'),
            title = cms.string('Tau lead. Track find.'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrkPt'),
            title = cms.string('Tau lead. Track Pt'),
        ),
        cms.PSet(
            filter = cms.string('tauTrkIsoCutBgEstWplusJetsEnriched'),
            title = cms.string('Tau Track iso.')
        ),
        cms.PSet(
            filter = cms.string('tauEcalIsoCutBgEstWplusJetsEnriched'),
            title = cms.string('Tau ECAL iso.')
        ),
        cms.PSet(
            filter = cms.string('tauMuonVetoBgEstWplusJetsEnriched'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('muTauPairBgEstWplusJetsEnriched'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairMt1MEtBgEstWplusJetsEnriched'),
            title = cms.string('M_{T}(Muon-MET) > 30 GeV')
        ),
        cms.PSet(
            filter = cms.string('centralJetVetoBgEstWplusJetsEnriched'),
            title = cms.string('central Jet Veto')
        ),
        cms.PSet(
            filter = cms.string('diMuonVetoBgEstWplusJetsEnriched'),
            title = cms.string('di-Muon Veto')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerWplusJetsEnriched',
                'tauHistManagerWplusJetsEnriched',
                'diTauCandidateHistManagerWplusJetsEnriched'
            )
        )
    )
)

#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstWplusJetsEnrichedAnalysisSequence = cms.Sequence(
    selectMuonsBgEstWplusJetsEnriched
   + selectTausBgEstWplusJetsEnriched
   + selectMuTauPairsBgEstWplusJetsEnriched
   + selectJetsBgEstWplusJetsEnriched 
   + selectEventsBgEstWplusJetsEnriched
   + analyzeEventsBgEstWplusJetsEnriched
)