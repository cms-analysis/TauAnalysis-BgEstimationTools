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

muonsForTauIdEffZtoMuTauGenMatrixFitCombRelIso = cms.EDFilter("PATMuonSelector",
    #
    # CV: apply loose cut on (relative) muon isolation only,
    #     in order to have still some discrimination between Ztautau/WplusJets and QCD left
    #     for "generalized matrix method"                                                   
    #                                              
    cut = cms.string('(userIsolation("pat::TrackIso") + userIsolation("pat::EcalIso")) < (0.15*pt)'),
    filter = cms.bool(False)
)

muonsForTauIdEffZtoMuTauGenMatrixFitPionVeto = copy.deepcopy(selectedPatMuonsPionVeto)

muonsForTauIdEffZtoMuTauGenMatrixFitTrk = copy.deepcopy(selectedPatMuonsTrk)

muonsForTauIdEffZtoMuTauGenMatrixFitTrkIP = copy.deepcopy(selectedPatMuonsTrkIP)

muonSelConfiguratorTauIdEffZtoMuTauGenMatrixFit = objSelConfigurator(
    [ muonsForTauIdEffZtoMuTauGenMatrixFitCombRelIso,
      muonsForTauIdEffZtoMuTauGenMatrixFitPionVeto,
      muonsForTauIdEffZtoMuTauGenMatrixFitTrk,
      muonsForTauIdEffZtoMuTauGenMatrixFitTrkIP ],
    src = "selectedPatMuonsPt10Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuonsForTauIdEffZtoMuTauGenMatrixFit = muonSelConfiguratorTauIdEffZtoMuTauGenMatrixFit.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patPFTauSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *

tausForTauIdEffZtoMuTauGenMatrixFitMuonVeto = copy.deepcopy(selectedPatTausMuonVeto)

tauSelConfiguratorTauIdEffZtoMuTauGenMatrixFit = objSelConfigurator(
    [ tausForTauIdEffZtoMuTauGenMatrixFitMuonVeto ],
    src = "selectedPatTausForMuTauLeadTrkPtCumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectTausForTauIdEffZtoMuTauGenMatrixFit = tauSelConfiguratorTauIdEffZtoMuTauGenMatrixFit.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of Z --> muon+ muon- hypotheses
#--------------------------------------------------------------------------------

# require muon candidates considered for Z --> mu+ mu- hypothesis
# to be reconstructed in muon system
# (with or without a track reconstructed in Pixel/SiStrip tracking detectors linked to it)
muonsLooseForZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string('isGlobalMuon | isStandAloneMuon | pt > 5.'),
    filter = cms.bool(False)
)

muonsTightForZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string('isGlobalMuon | isStandAloneMuon'),
    filter = cms.bool(False)
)

allDiMuPairZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsTightForZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit'),
    srcLeg2 = cms.InputTag('muonsLooseForZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit'),
    dRmin12 = cms.double(0.5),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                        
    verbosity = cms.untracked.int32(0)
)

selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag("allDiMuPairZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit"),                                   
    cut = cms.string('p4Vis.mass > 70. & p4Vis.mass < 110.'),
    filter = cms.bool(False)
)

produceDiMuPairsTauIdEffZtoMuTauGenMatrixFit = cms.Sequence(
    muonsLooseForZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit * muonsTightForZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit
   * allDiMuPairZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit
   * selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit
)

#--------------------------------------------------------------------------------  
# produce collection of muon + tau-jet combinations
#--------------------------------------------------------------------------------

from TauAnalysis.CandidateTools.muTauPairProduction_cff import *
from TauAnalysis.CandidateTools.resolutions_cfi import *

muTauPairsTauIdEffZtoMuTauGenMatrixFit = allMuTauPairs.clone(
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsForTauIdEffZtoMuTauGenMatrixFitTrkIPcumulative'),
    srcLeg2 = cms.InputTag('tausForTauIdEffZtoMuTauGenMatrixFitMuonVetoCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('patMETs'),
    recoMode = cms.string(""),
    doSVreco = cms.bool(False)
)

#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *

muonHistManagerTauIdEffZtoMuTauGenMatrixFit = copy.deepcopy(muonHistManager)
muonHistManagerTauIdEffZtoMuTauGenMatrixFit.pluginName = 'muonHistManagerTauIdEffZtoMuTauGenMatrixFit'
muonHistManagerTauIdEffZtoMuTauGenMatrixFit.muonSource = 'muonsForTauIdEffZtoMuTauGenMatrixFitTrkIPcumulative'

tauHistManagerTauIdEffZtoMuTauGenMatrixFit = copy.deepcopy(tauHistManager)
tauHistManagerTauIdEffZtoMuTauGenMatrixFit.pluginName = 'tauHistManagerTauIdEffZtoMuTauGenMatrixFit'
tauHistManagerTauIdEffZtoMuTauGenMatrixFit.tauSource = 'tausForTauIdEffZtoMuTauGenMatrixFitMuonVetoCumulative'

diTauCandidateHistManagerTauIdEffZtoMuTauGenMatrixFit = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTauIdEffZtoMuTauGenMatrixFit.pluginName = 'diTauCandidateHistManagerTauIdEffZtoMuTauGenMatrixFit'
diTauCandidateHistManagerTauIdEffZtoMuTauGenMatrixFit.diTauCandidateSource = 'muTauPairsTauIdEffZtoMuTauGenMatrixFit'
diTauCandidateHistManagerTauIdEffZtoMuTauGenMatrixFit.visMassHypothesisSource = cms.InputTag('')

from TauAnalysis.BgEstimationTools.bgEstBinGridZtoMuTau_cfi import *

dataBinnerTauIdEffZtoMuTauGenMatrixFit = copy.deepcopy(dataBinner)
dataBinnerTauIdEffZtoMuTauGenMatrixFit.pluginName = 'dataBinnerTauIdEffZtoMuTauGenMatrixFit'

binningTauIdEffZtoMuTauGenMatrixFit_ewkTauId = copy.deepcopy(binning_ewkTauId)
binningTauIdEffZtoMuTauGenMatrixFit_ewkTauId.extractor.src = 'tausForTauIdEffZtoMuTauGenMatrixFitMuonVetoCumulative'

binningTauIdEffZtoMuTauGenMatrixFit_relMuonIso = copy.deepcopy(binning_relMuonIso)
binningTauIdEffZtoMuTauGenMatrixFit_relMuonIso.extractor.src = 'muonsForTauIdEffZtoMuTauGenMatrixFitTrkIPcumulative'

binningTauIdEffZtoMuTauGenMatrixFit_diTauMt1MET = copy.deepcopy(binning_diTauMt1MET)
binningTauIdEffZtoMuTauGenMatrixFit_diTauMt1MET.extractor.src = 'muTauPairsTauIdEffZtoMuTauGenMatrixFit'

binningTauIdEffZtoMuTauGenMatrixFit_diTauAbsCharge = copy.deepcopy(binning_diTauAbsCharge)
binningTauIdEffZtoMuTauGenMatrixFit_diTauAbsCharge.extractor.src = 'muTauPairsTauIdEffZtoMuTauGenMatrixFit'

tauIdEffBinningZtoMuTau_genMatrix1d = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_genMatrix1d"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauGenMatrixFit_ewkTauId
    )
)

tauIdEffBinningZtoMuTau_genMatrix3d = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_genMatrix3d"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauGenMatrixFit_ewkTauId,
        binningTauIdEffZtoMuTauGenMatrixFit_relMuonIso,
        binningTauIdEffZtoMuTauGenMatrixFit_diTauMt1MET
    )
)

tauIdEffBinningZtoMuTau_genMatrix4d = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_genMatrix4d"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauGenMatrixFit_ewkTauId,
        binningTauIdEffZtoMuTauGenMatrixFit_relMuonIso,
        binningTauIdEffZtoMuTauGenMatrixFit_diTauMt1MET,
        binningTauIdEffZtoMuTauGenMatrixFit_diTauAbsCharge
    )
)

analyzeEventsTauIdEffZtoMuTauGenMatrixFit = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('TauIdEffAnalyzerZtoMuTauGenMatrixFit'), 
                            
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
            pluginName = cms.string('muonCombRelIsoCutTauIdEffZtoMuTauGenMatrixFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauGenMatrixFitCombRelIsoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muonAntiPionCutTauIdEffZtoMuTauGenMatrixFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauGenMatrixFitPionVetoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muonTrkIPcutTauIdEffZtoMuTauGenMatrixFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauGenMatrixFitTrkCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('tauEtaCutTauIdEffZtoMuTauGenMatrixFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauGenMatrixFitTrkIPcumulative'),
            minNumber = cms.uint32(1)
        ),
        evtSelTauAntiOverlapWithMuonsVeto,
        evtSelTauEta,
        evtSelTauPt,
        evtSelTauLeadTrk,
        evtSelTauLeadTrkPt,
        cms.PSet(
            pluginName = cms.string('tauMuonVetoTauIdEffZtoMuTauGenMatrixFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('tausForTauIdEffZtoMuTauGenMatrixFitMuonVetoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauGenMatrixFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauGenMatrixFit'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('diMuPairZmumuHypothesisVeto'),
            pluginType = cms.string('PATCandViewMaxEventSelector'),
            src = cms.InputTag('selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit'),
            maxNumber = cms.uint32(0)
        ),
        ##cms.PSet(
        ##    pluginName = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTauGenMatrixFit'),
        ##    pluginType = cms.string('PATCandViewMaxEventSelector'),
        ##    src = cms.InputTag('tausForTauIdEffZtoMuTauGenMatrixFitMuonVeto'),
        ##    maxNumber = cms.uint32(0)
        ##),
        cms.PSet(
            pluginName = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTauGenMatrixFit'),
            pluginType = cms.string('PATCandViewMaxEventSelector'),
            src = cms.InputTag('muonsTightForZmumuHypothesesTauIdEffZtoMuTauGenMatrixFit'),
            maxNumber = cms.uint32(1)
        )
    ),
  
    analyzers = cms.VPSet(
        muonHistManagerTauIdEffZtoMuTauGenMatrixFit,
        tauHistManagerTauIdEffZtoMuTauGenMatrixFit,
        diTauCandidateHistManagerTauIdEffZtoMuTauGenMatrixFit,
        caloMEtHistManager,
        pfMEtHistManager,
        dataBinnerTauIdEffZtoMuTauGenMatrixFit,
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinner2regions'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_genMatrix1d,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResults2regions')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManager2regions'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_genMatrix1d,
            histManagers = cms.VPSet(
                muonHistManagerTauIdEffZtoMuTauGenMatrixFit,
                tauHistManagerTauIdEffZtoMuTauGenMatrixFit,
                diTauCandidateHistManagerTauIdEffZtoMuTauGenMatrixFit
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
            filter = cms.string('muonCombRelIsoCutTauIdEffZtoMuTauGenMatrixFit'),
            title = cms.string('Muon Track + ECAL relative iso.')
        ),
        cms.PSet(
            filter = cms.string('muonAntiPionCutTauIdEffZtoMuTauGenMatrixFit'),
            title = cms.string('Muon pi-Veto')
        ),
        cms.PSet(
            filter = cms.string('muonTrkIPcutTauIdEffZtoMuTauGenMatrixFit'),
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
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTauGenMatrixFit'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTauGenMatrixFit'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('diMuPairZmumuHypothesisVeto'),
            title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ),        
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTauGenMatrixFit'),
            title = cms.string('num. global Muons < 2')
        ),
        ##cms.PSet(
        ##    filter = cms.string('uniqueTauCandidateCutTauIdEffZtoMuTauGenMatrixFit'),
        ##    title = cms.string('num. Tau-Jet Candidates < 2')
        ##),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauGenMatrixFit',
                'tauHistManagerTauIdEffZtoMuTauGenMatrixFit',
                'diTauCandidateHistManagerTauIdEffZtoMuTauGenMatrixFit',
                'caloMEtHistManager',
                'pfMEtHistManager',
                'dataBinnerTauIdEffZtoMuTauGenMatrixFit',
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

bgEstTauIdEffZtoMuTauGenMatrixFitAnalysisSequence = cms.Sequence(
    selectMuonsForTauIdEffZtoMuTauGenMatrixFit 
   + selectTausForTauIdEffZtoMuTauGenMatrixFit
   + produceDiMuPairsTauIdEffZtoMuTauGenMatrixFit
   + muTauPairsTauIdEffZtoMuTauGenMatrixFit
   + analyzeEventsTauIdEffZtoMuTauGenMatrixFit 
)
