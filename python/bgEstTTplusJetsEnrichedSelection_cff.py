import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

#--------------------------------------------------------------------------------
# select ttbar + jets background enriched event sample
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------  
# produce collection of pat::Muons
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patMuonSelection_cfi import *

muonsBgEstTTplusJetsEnrichedTrkIso = copy.deepcopy(selectedLayer1MuonsTrkIso)
muonsBgEstTTplusJetsEnrichedTrkIso.sumPtMax = cms.double(2.)

muonsBgEstTTplusJetsEnrichedEcalIso = copy.deepcopy(selectedLayer1MuonsEcalIso)
muonsBgEstTTplusJetsEnrichedEcalIso.cut = cms.string('ecalIso < 2.')

muonSelConfiguratorBgEstTTplusJetsEnriched = objSelConfigurator(
    [ muonsBgEstTTplusJetsEnrichedTrkIso,
      muonsBgEstTTplusJetsEnrichedEcalIso ],
    src = "selectedLayer1MuonsPt15Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuonsBgEstTTplusJetsEnriched = muonSelConfiguratorBgEstTTplusJetsEnriched.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of muon + tau-jet combinations
#--------------------------------------------------------------------------------

muTauPairsBgEstTTplusJetsEnriched = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsBgEstTTplusJetsEnrichedEcalIsoCumulative'),
    srcLeg2 = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('layer1METs'),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)
)

muTauPairsBgEstTTplusJetsEnrichedZeroCharge = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsBgEstTTplusJetsEnriched'),                                                   
    cut = cms.string('abs(charge) < 0.5'),
    filter = cms.bool(False)
)

selectMuTauPairsBgEstTTplusJetsEnriched = cms.Sequence(muTauPairsBgEstTTplusJetsEnriched + muTauPairsBgEstTTplusJetsEnrichedZeroCharge)

#--------------------------------------------------------------------------------  
# produce collection of pat::Jets used for central jet veto
# (in order to reject QCD di-jet events)
#--------------------------------------------------------------------------------

jetsBgEstTTplusJetsEnrichedAntiOverlapWithLeptonsVeto = cms.EDFilter("PATJetAntiOverlapSelector",
    src = cms.InputTag("selectedLayer1JetsEt20Cumulative"),                                                                  
    srcNotToBeFiltered = cms.VInputTag(
        "selectedLayer1ElectronsTrkIPcumulative",
        "muonsBgEstTTplusJetsEnrichedEcalIsoCumulative",
        "selectedLayer1TausForMuTauMuonVetoCumulative"
    ),
    dRmin = cms.double(0.7),
    filter = cms.bool(False)                                           
)

jetsBgEstTTplusJetsEnrichedEt40 = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("jetsBgEstTTplusJetsEnrichedAntiOverlapWithLeptonsVeto"),
    cut = cms.string('et > 40.'),
    filter = cms.bool(False)
)

jetsBgEstTTplusJetsEnrichedEt40bTag = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("jetsBgEstTTplusJetsEnrichedEt40"),
    cut = cms.string('et > 40.'),
    filter = cms.bool(False)
)

jetsBgEstTTplusJetsEnrichedEt60 = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("jetsBgEstTTplusJetsEnrichedAntiOverlapWithLeptonsVeto"),
    cut = cms.string('et > 60.'),
    filter = cms.bool(False)
)

selectJetsBgEstTTplusJetsEnriched = cms.Sequence(
    jetsBgEstTTplusJetsEnrichedAntiOverlapWithLeptonsVeto
   + jetsBgEstTTplusJetsEnrichedEt40 + jetsBgEstTTplusJetsEnrichedEt40bTag
   + jetsBgEstTTplusJetsEnrichedEt60
)

#--------------------------------------------------------------------------------  
# produce boolean event selection flags
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.selectZtoMuTau_cff import *

cfgMuonTrkIsoCutBgEstTTplusJetsEnriched = copy.deepcopy(cfgMuonTrkIsoCut)
cfgMuonTrkIsoCutBgEstTTplusJetsEnriched.pluginName = cms.string('muonTrkIsoCutBgEstTTplusJetsEnriched')
cfgMuonTrkIsoCutBgEstTTplusJetsEnriched.src_cumulative = cms.InputTag('muonsBgEstTTplusJetsEnrichedTrkIsoCumulative')
cfgMuonTrkIsoCutBgEstTTplusJetsEnriched.systematics = cms.vstring()

cfgMuonEcalIsoCutBgEstTTplusJetsEnriched = copy.deepcopy(cfgMuonEcalIsoCut)
cfgMuonEcalIsoCutBgEstTTplusJetsEnriched.pluginName = cms.string('muonEcalIsoCutBgEstTTplusJetsEnriched')
cfgMuonEcalIsoCutBgEstTTplusJetsEnriched.src_cumulative = cms.InputTag('muonsBgEstTTplusJetsEnrichedEcalIsoCumulative')
cfgMuonEcalIsoCutBgEstTTplusJetsEnriched.systematics = cms.vstring()

cfgMuTauPairBgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('muTauPairBgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsBgEstTTplusJetsEnriched'),
    minNumber = cms.uint32(1)
)

cfgMuTauPairZeroChargeBgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('muTauPairZeroChargeBgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsBgEstTTplusJetsEnrichedZeroCharge'),
    minNumber = cms.uint32(1)
)

cfgJetsEt40BgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('jetsEt40BgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('jetsBgEstTTplusJetsEnrichedEt40'),
    minNumber = cms.uint32(2)
)

cfgJetEt40bTagBgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('jetEt40bTagBgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('jetsBgEstTTplusJetsEnrichedEt40bTag'),
    minNumber = cms.uint32(1)
)

cfgJetEt60BgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('jetEt60BgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('jetsBgEstTTplusJetsEnrichedEt60'),
    minNumber = cms.uint32(1)
)

cfgDiMuonVetoBgEstTTplusJetsEnriched = cms.PSet(
    pluginName = cms.string('diMuonVetoBgEstTTplusJetsEnriched'),
    pluginType = cms.string('PATCandViewMaxEventSelector'),
    src = cms.InputTag('selectedLayer1MuonsGlobalIndividual'),
    maxNumber = cms.uint32(1)
)

evtSelConfiguratorBgEstTTplusJetsEnriched = eventSelFlagProdConfigurator(
    [ cfgMuonTrkIsoCutBgEstTTplusJetsEnriched,
      cfgMuonEcalIsoCutBgEstTTplusJetsEnriched,
      cfgMuTauPairBgEstTTplusJetsEnriched,
      cfgMuTauPairZeroChargeBgEstTTplusJetsEnriched,
      cfgJetsEt40BgEstTTplusJetsEnriched,
      cfgJetEt40bTagBgEstTTplusJetsEnriched,
      cfgJetEt60BgEstTTplusJetsEnriched,
      cfgDiMuonVetoBgEstTTplusJetsEnriched ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

selectEventsBgEstTTplusJetsEnriched = evtSelConfiguratorBgEstTTplusJetsEnriched.configure()

#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *

muonHistManagerTTplusJetsEnriched = copy.deepcopy(muonHistManager)
muonHistManagerTTplusJetsEnriched.pluginName = cms.string('muonHistManagerTTplusJetsEnriched')
muonHistManagerTTplusJetsEnriched.muonSource = cms.InputTag('muonsBgEstTTplusJetsEnrichedEcalIsoCumulative')

tauHistManagerTTplusJetsEnriched = copy.deepcopy(tauHistManager)
tauHistManagerTTplusJetsEnriched.pluginName = cms.string('tauHistManagerTTplusJetsEnriched')
tauHistManagerTTplusJetsEnriched.tauSource = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative')

diTauCandidateHistManagerTTplusJetsEnriched = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTTplusJetsEnriched.pluginName = cms.string('diTauCandidateHistManagerTTplusJetsEnriched')
diTauCandidateHistManagerTTplusJetsEnriched.diTauCandidateSource = cms.InputTag('muTauPairsBgEstTTplusJetsEnrichedZeroCharge')

analyzeEventsBgEstTTplusJetsEnriched = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('BgEstTemplateAnalyzer_TTplusJetsEnriched'), 
                            
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
            pluginName = cms.string('muonTrkIsoCutBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonTrkIsoCutBgEstTTplusJetsEnriched', 'cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('muonEcalIsoCutBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muonEcalIsoCutBgEstTTplusJetsEnriched', 'cumulative')
        ),
        evtSelTauAntiOverlapWithMuonsVeto,
        evtSelTauEta,
        evtSelTauPt,
        evtSelTauLeadTrk,
        evtSelTauLeadTrkPt,
        evtSelTauTrkIso,
        evtSelTauEcalIso,
        evtSelTauProng,
        evtSelTauCharge,
        evtSelTauMuonVeto,
        cms.PSet(
            pluginName = cms.string('muTauPairBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairBgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairZeroChargeBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('muTauPairZeroChargeBgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('jetEt40BgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('jetsEt40BgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('jetEt40bTagBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('jetEt40bTagBgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('jetEt60BgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('jetEt60BgEstTTplusJetsEnriched')
        ),
        cms.PSet(
            pluginName = cms.string('diMuonVetoBgEstTTplusJetsEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diMuonVetoBgEstTTplusJetsEnriched')
        )
    ),
  
    analyzers = cms.VPSet(
        muonHistManagerTTplusJetsEnriched,
        tauHistManagerTTplusJetsEnriched,
        diTauCandidateHistManagerTTplusJetsEnriched
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
            title = cms.string('-2.1 < eta(Tau) < +2.1'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauPt'),
            title = cms.string('Pt(Tau) > 20 GeV'),
        ),
        cms.PSet(
            filter = cms.string('muonTrkIsoCutBgEstTTplusJetsEnriched'),
            title = cms.string('Muon Track iso.')
        ),
        cms.PSet(
            filter = cms.string('muonEcalIsoCutBgEstTTplusJetsEnriched'),
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
            filter = cms.string('evtSelTauTrkIso'),
            title = cms.string('Tau Track iso.'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauEcalIso'),
            title = cms.string('Tau ECAL iso.'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauProng'),
            title = cms.string('Tau 1||3-Prong'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauCharge'),
            title = cms.string('Charge(Tau) = +/-1'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauMuonVeto'),
            title = cms.string('Tau mu-Veto'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            filter = cms.string('muTauPairBgEstTTplusJetsEnriched'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairZeroChargeBgEstTTplusJetsEnriched'),
            title = cms.string('Charge(Muon+Tau) = 0')
        ),
        cms.PSet(
            filter = cms.string('jetEt40BgEstTTplusJetsEnriched'),
            title = cms.string('two E_{T} > 40 GeV Jets')
        ),
        cms.PSet(
            filter = cms.string('jetEt40bTagBgEstTTplusJetsEnriched'),
            title = cms.string('one E_{T} > 40 GeV Jet with b-Tag')
        ),
        cms.PSet(
            filter = cms.string('jetEt60BgEstTTplusJetsEnriched'),
            title = cms.string('one E_{T} > 60 GeV Jet')
        ),
        cms.PSet(
            filter = cms.string('diMuonVetoBgEstTTplusJetsEnriched'),
            title = cms.string('di-Muon Veto')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTTplusJetsEnriched',
                'tauHistManagerTTplusJetsEnriched',
                'diTauCandidateHistManagerTTplusJetsEnriched'
            )
        )
    )
)

#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstTTplusJetsEnrichedAnalysisSequence = cms.Sequence(
    selectMuonsBgEstTTplusJetsEnriched
   + selectMuTauPairsBgEstTTplusJetsEnriched
   + selectJetsBgEstTTplusJetsEnriched 
   + selectEventsBgEstTTplusJetsEnriched
   + analyzeEventsBgEstTTplusJetsEnriched
)