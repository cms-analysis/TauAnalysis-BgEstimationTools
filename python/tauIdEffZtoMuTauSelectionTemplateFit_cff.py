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

muonsForTauIdEffZtoMuTauTemplateFitCombRelIso = cms.EDFilter("PATMuonSelector",
    cut = cms.string('(userIsolation("pat::TrackIso") + userIsolation("pat::EcalIso")) < (0.06*pt)'),
    filter = cms.bool(False)
)

muonsForTauIdEffZtoMuTauTemplateFitPionVeto = copy.deepcopy(selectedPatMuonsPionVeto)

muonsForTauIdEffZtoMuTauTemplateFitTrk = copy.deepcopy(selectedPatMuonsTrk)

muonsForTauIdEffZtoMuTauTemplateFitTrkIP = copy.deepcopy(selectedPatMuonsTrkIP)

muonSelConfiguratorTauIdEffZtoMuTauTemplateFit = objSelConfigurator(
    [ muonsForTauIdEffZtoMuTauTemplateFitCombRelIso,
      muonsForTauIdEffZtoMuTauTemplateFitPionVeto,
      muonsForTauIdEffZtoMuTauTemplateFitTrk,
      muonsForTauIdEffZtoMuTauTemplateFitTrkIP ],
    src = "selectedPatMuonsPt15Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuonsForTauIdEffZtoMuTauTemplateFit = muonSelConfiguratorTauIdEffZtoMuTauTemplateFit.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patPFTauSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *

tausForTauIdEffZtoMuTauTemplateFitMuonVeto = copy.deepcopy(selectedPatTausMuonVeto)

tauSelConfiguratorTauIdEffZtoMuTauTemplateFit = objSelConfigurator(
    [ tausForTauIdEffZtoMuTauTemplateFitMuonVeto ],
    src = "selectedPatTausForMuTauLeadTrkPtCumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectTausForTauIdEffZtoMuTauTemplateFit = tauSelConfiguratorTauIdEffZtoMuTauTemplateFit.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of pat::MEt objects
#--------------------------------------------------------------------------------

pfMEtForTauIdEffZtoMuTauTemplateFitPt20 = cms.EDFilter("PATMETSelector",
    src = cms.InputTag("patMETs"),                                 
    cut = cms.string('pt > 20.'),
    filter = cms.bool(False)
)

selectMEtForTauIdEffZtoMuTauTemplateFit = cms.Sequence(pfMEtForTauIdEffZtoMuTauTemplateFitPt20)

#--------------------------------------------------------------------------------  
# produce collection of Z --> muon+ muon- hypotheses
#--------------------------------------------------------------------------------

# require muon candidates considered for Z --> mu+ mu- hypothesis
# to be reconstructed in muon system
# (with or without a track reconstructed in Pixel/SiStrip tracking detectors linked to it)
muonsLooseForZmumuHypothesesTauIdEffZtoMuTauTemplateFit = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string('isGlobalMuon | isStandAloneMuon | pt > 5.'),
    filter = cms.bool(False)
)

muonsTightForZmumuHypothesesTauIdEffZtoMuTauTemplateFit = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string('isGlobalMuon | isStandAloneMuon'),
    filter = cms.bool(False)
)

allDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsTightForZmumuHypothesesTauIdEffZtoMuTauTemplateFit'),
    srcLeg2 = cms.InputTag('muonsLooseForZmumuHypothesesTauIdEffZtoMuTauTemplateFit'),
    dRmin12 = cms.double(0.5),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                        
    verbosity = cms.untracked.int32(0)
)

selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag("allDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit"),                                   
    cut = cms.string('p4Vis.mass > 70. & p4Vis.mass < 110.'),
    filter = cms.bool(False)
)

produceDiMuPairsTauIdEffZtoMuTauTemplateFit = cms.Sequence(
    muonsLooseForZmumuHypothesesTauIdEffZtoMuTauTemplateFit * muonsTightForZmumuHypothesesTauIdEffZtoMuTauTemplateFit
   * allDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit
   * selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit
)

#--------------------------------------------------------------------------------  
# produce collection of muon + tau-jet combinations
#--------------------------------------------------------------------------------

from TauAnalysis.CandidateTools.resolutions_cfi import *

muTauPairsTauIdEffZtoMuTauTemplateFit = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsForTauIdEffZtoMuTauTemplateFitTrkIPcumulative'),
    srcLeg2 = cms.InputTag('tausForTauIdEffZtoMuTauTemplateFitMuonVetoCumulative'),
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

muTauPairsTauIdEffZtoMuTauTemplateFitMt1MET = cms.EDFilter("PATMuTauPairSelector",
    cut = cms.string('mt1MET < 40.'),
    filter = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauTemplateFitPzetaDiff = cms.EDFilter("PATMuTauPairSelector",
    cut = cms.string('(pZeta - 1.5*pZetaVis) > -20.'),
    filter = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauTemplateFitBackToBack = cms.EDFilter("PATMuTauPairSelector",
    cut = cms.string('dPhi12 > 2.793'),
    filter = cms.bool(False)
)

muTauPairConfiguratorTauIdEffZtoMuTauTemplateFit = objSelConfigurator(
    [ muTauPairsTauIdEffZtoMuTauTemplateFitMt1MET,
      muTauPairsTauIdEffZtoMuTauTemplateFitPzetaDiff,
      muTauPairsTauIdEffZtoMuTauTemplateFitBackToBack ],
    src = "muTauPairsTauIdEffZtoMuTauTemplateFit",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuTauPairsTauIdEffZtoMuTauTemplateFit = muTauPairConfiguratorTauIdEffZtoMuTauTemplateFit.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *

muonHistManagerTauIdEffZtoMuTauTemplateFit = copy.deepcopy(muonHistManager)
muonHistManagerTauIdEffZtoMuTauTemplateFit.pluginName = 'muonHistManagerTauIdEffZtoMuTauTemplateFit'
muonHistManagerTauIdEffZtoMuTauTemplateFit.muonSource = 'muonsForTauIdEffZtoMuTauTemplateFitTrkIPcumulative'

tauHistManagerTauIdEffZtoMuTauTemplateFit = copy.deepcopy(tauHistManager)
tauHistManagerTauIdEffZtoMuTauTemplateFit.pluginName = 'tauHistManagerTauIdEffZtoMuTauTemplateFit'
tauHistManagerTauIdEffZtoMuTauTemplateFit.tauSource = 'tausForTauIdEffZtoMuTauTemplateFitMuonVetoCumulative'

diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit.pluginName = 'diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit'
diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit.diTauCandidateSource = 'muTauPairsTauIdEffZtoMuTauTemplateFitBackToBackCumulative'
diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit.visMassHypothesisSource = cms.InputTag('')

from TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauHistManager_cfi import *
tauIdEffZtoMuTauHistManagerTemplateFit = copy.deepcopy(tauIdEffZtoMuTauHistManager)
tauIdEffZtoMuTauHistManagerTemplateFit.pluginName = 'tauIdEffZtoMuTauHistManagerTemplateFit'
tauIdEffZtoMuTauHistManagerTemplateFit.muonSource = 'muonsForTauIdEffZtoMuTauTemplateFitTrkIPcumulative'  
tauIdEffZtoMuTauHistManagerTemplateFit.tauSource = 'tausForTauIdEffZtoMuTauTemplateFitMuonVetoCumulative'
tauIdEffZtoMuTauHistManagerTemplateFit.diTauSource = 'muTauPairsTauIdEffZtoMuTauTemplateFitBackToBackCumulative'
tauIdEffZtoMuTauHistManagerTemplateFit.diTauChargeSignExtractor.src = 'muTauPairsTauIdEffZtoMuTauTemplateFitBackToBackCumulative'

from TauAnalysis.BgEstimationTools.bgEstBinGridZtoMuTau_cfi import *

dataBinnerTauIdEffZtoMuTauTemplateFit = copy.deepcopy(dataBinner)
dataBinnerTauIdEffZtoMuTauTemplateFit.pluginName = 'dataBinnerTauIdEffZtoMuTauTemplateFit'

binningTauIdEffZtoMuTauTemplateFit_ewkTauId = copy.deepcopy(binning_ewkTauId)
binningTauIdEffZtoMuTauTemplateFit_ewkTauId.extractor.src = 'tausForTauIdEffZtoMuTauTemplateFitMuonVetoCumulative'

tauIdEffBinningZtoMuTau_template1d = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_template1d"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauTemplateFit_ewkTauId
    )
)

analyzeEventsTauIdEffZtoMuTauTemplateFit = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('TauIdEffAnalyzerZtoMuTauTemplateFit'), 
                            
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
            pluginName = cms.string('muonCombRelIsoCutTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauTemplateFitCombRelIsoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muonAntiPionCutTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauTemplateFitPionVetoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muonTrkIPcutTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauTemplateFitTrkCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('tauEtaCutTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauTemplateFitTrkIPcumulative'),
            minNumber = cms.uint32(1)
        ),
        evtSelTauAntiOverlapWithMuonsVeto,
        evtSelTauEta,
        evtSelTauPt,
        evtSelTauLeadTrk,
        evtSelTauLeadTrkPt,
        cms.PSet(
            pluginName = cms.string('tauMuonVetoTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('tausForTauIdEffZtoMuTauTemplateFitMuonVetoCumulative'),
            minNumber = cms.uint32(1)
        ),    
        cms.PSet(
            pluginName = cms.string('pfMEtPt20CutTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('pfMEtForTauIdEffZtoMuTauTemplateFitPt20'),
            minNumber = cms.uint32(1)
        ),                                                   
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauTemplateFit'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairMt1METtauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauTemplateFitMt1METcumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauTemplateFitPzetaDiffCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairBackToBackTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauTemplateFitBackToBackCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('diMuPairZmumuHypothesisVeto'),
            pluginType = cms.string('PATCandViewMaxEventSelector'),
            src = cms.InputTag('selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit'),
            maxNumber = cms.uint32(0)
        ),
        ##cms.PSet(
        ##    pluginName = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTauTemplateFit'),
        ##    pluginType = cms.string('PATCandViewMaxEventSelector'),
        ##    src = cms.InputTag('tausForTauIdEffZtoMuTauTemplateFitMuonVeto'),
        ##    maxNumber = cms.uint32(0)
        ##),
        cms.PSet(
            pluginName = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMaxEventSelector'),
            src = cms.InputTag('muonsTightForZmumuHypothesesTauIdEffZtoMuTauTemplateFit'),
            maxNumber = cms.uint32(1)
        )
    ),
  
    analyzers = cms.VPSet(
        muonHistManagerTauIdEffZtoMuTauTemplateFit,
        tauHistManagerTauIdEffZtoMuTauTemplateFit,
        diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit,
        caloMEtHistManager,
        pfMEtHistManager,
        tauIdEffZtoMuTauHistManagerTemplateFit,
        dataBinnerTauIdEffZtoMuTauTemplateFit,
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinner2regions'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_template1d,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResults2regions')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManager2regions'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_template1d,
            histManagers = cms.VPSet(
                muonHistManagerTauIdEffZtoMuTauTemplateFit,
                tauHistManagerTauIdEffZtoMuTauTemplateFit,
                diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit,
                tauIdEffZtoMuTauHistManagerTemplateFit
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
            filter = cms.string('muonCombRelIsoCutTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('Muon Track + ECAL relative iso.')
        ),
        cms.PSet(
            filter = cms.string('muonAntiPionCutTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('Muon pi-Veto')
        ),
        cms.PSet(
            filter = cms.string('muonTrkIPcutTauIdEffZtoMuTauTemplateFit'),
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
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('pfMEtPt20CutTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('PFMEt Pt > 20 GeV')
        ), 
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairMt1METtauIdEffZtoMuTauTemplateFit'),
            title = cms.string('M_{T}(Muon-MET) < 50 GeV')
        ),
        cms.PSet(
            filter = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV')
        ),
        cms.PSet(
            filter = cms.string('muTauPairBackToBackTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('dPhi(Muon,Tau) > 160 deg.')
        ),
        cms.PSet(
            filter = cms.string('diMuPairZmumuHypothesisVeto'),
            title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ),        
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('num. global Muons < 2')
        ),
        ##cms.PSet(
        ##    filter = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTauTemplateFit'),
        ##    title = cms.string('num. Tau-Jet Candidates < 2')
        ##),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauTemplateFit',
                'tauHistManagerTauIdEffZtoMuTauTemplateFit',
                'diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit',
                'caloMEtHistManager',
                'pfMEtHistManager',
                'tauIdEffZtoMuTauHistManagerTemplateFit',
                'dataBinnerTauIdEffZtoMuTauTemplateFit',
                'tauIdEffDataBinner2regions',
                'tauIdEffBinGridHistManager2regions'
            ),
            replace = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauTemplateFit.muonSource = muonsForTauIdEffZtoMuTauTemplateFitTrkIPcumulative',
                'diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit.diTauCandidateSource = muTauPairsTauIdEffZtoMuTauTemplateFitBackToBackCumulative',
                'pfMEtHistManager.metSource = pfMEtForTauIdEffZtoMuTauTemplateFitPt20'
            )
        )
    )
)

#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstTauIdEffZtoMuTauTemplateFitAnalysisSequence = cms.Sequence(
    selectMuonsForTauIdEffZtoMuTauTemplateFit 
   + selectTausForTauIdEffZtoMuTauTemplateFit 
   + selectMEtForTauIdEffZtoMuTauTemplateFit
   + produceDiMuPairsTauIdEffZtoMuTauTemplateFit
   + muTauPairsTauIdEffZtoMuTauTemplateFit + selectMuTauPairsTauIdEffZtoMuTauTemplateFit 
   + analyzeEventsTauIdEffZtoMuTauTemplateFit 
)
