import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('runTauIdEffAnalysisZtoMuTau')

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'IDEAL_V12::All'

#--------------------------------------------------------------------------------
# import sequences for PAT-tuple production
process.load("TauAnalysis.Configuration.producePatTuple_cff")
process.load("TauAnalysis.Configuration.producePatTupleZtoMuTauSpecific_cff")

# import sequence for event selection
process.load("TauAnalysis.Configuration.selectZtoMuTau_cff")

# import configuration parameters for submission of jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)
from TauAnalysis.BgEstimationTools.recoSampleDefinitionsTauIdEffZtoMuTau_cfi import *
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.saveZtoMuTauPlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('plotsTauIdEffZtoMuTau_temp.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) 
    #input = cms.untracked.int32(1000)    
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/v/veelken/scratch0/CMSSW_2_2_10/src/TauAnalysis/Configuration/test/muTauSkim.root'
        #'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_Ztautau_part01.root'
        #'rfio:/castor/cern.ch/user/v/veelken/patTuples/tauIdEff/patTupleZtoMuTau_Ztautau_part01.root'
        #'file:/afs/cern.ch/user/v/veelken/scratch0/CMSSW_2_2_10/src/TauAnalysis/BgEstimationTools/test/patTuple.root'
    )
    #skipBadFiles = cms.untracked.bool(True) 
)

#--------------------------------------------------------------------------------
# import utility function for switching pat::Tau input
# to different reco::Tau collection stored on AOD
from PhysicsTools.PatAlgos.tools.tauTools import * 

# comment-out to take reco::CaloTaus instead of reco::PFTaus
# as input for pat::Tau production
#switchToCaloTau(process)

# comment-out to take shrinking dR = 5.0/Et(PFTau) signal cone
# instead of fixed dR = 0.07 signal cone reco::PFTaus
# as input for pat::Tau production
switchToPFTauShrinkingCone(process)
setattr(process.allLayer1Taus.tauIDSources, "ewkTauId", cms.InputTag('ewkTauId'))
#switchToPFTauFixedCone(process)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::METs
from TauAnalysis.Configuration.tools.metTools import *

# uncomment to add pfMET
# (set boolean parameter to true/false to enable/disable type-1 MET corrections)
addPFMet(process, correct = False)

# uncomment to replace caloMET by pfMET in all di-tau objects
process.load("TauAnalysis.CandidateTools.diTauPairProductionAllKinds_cff")
replaceMETforDiTaus(process, cms.InputTag('layer1METs'), cms.InputTag('layer1PFMETs'))
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for changing cut values
from TauAnalysis.Configuration.tools.changeCut import changeCut

# loosen cuts on muon Track and ECAL isolation
#changeCut(process, "selectedLayer1MuonsTrkIso", 8., attribute = "sumPtMax")
#changeCut(process, "selectedLayer1MuonsEcalIso", "ecalIso < 8.")

# disable tau id. & isolation cuts
changeCut(process, "selectedLayer1TausForMuTauLeadTrk", "tauID('leadingTrackFinding') > -1.")
changeCut(process, "selectedLayer1TausForMuTauLeadTrkPt", "tauID('leadingTrackPtCut') > -1.")
changeCut(process, "selectedLayer1TausForMuTauTrkIso", "tauID('trackIsolation') > -1.")
changeCut(process, "selectedLayer1TausForMuTauEcalIso", "tauID('ecalIsolation') > -1.")
changeCut(process, "selectedLayer1TausForMuTauProng", "signalTracks.size() > -1")
changeCut(process, "selectedLayer1TausForMuTauCharge", "abs(charge) > -1")
changeCut(process, "selectedLayer1TausForMuTauMuonVeto", "tauID('againstMuon') > -1.")

# disable cut on muon + tau-jet charge
changeCut(process, "selectedMuTauPairsZeroCharge", "charge > -1000.")

# require muon and tau-jet to be "back-to-back" (acoplanarity angle > 170 degrees)
# (in order to reduce contamination from quark/gluon jets in Z --> tau+ tau- signal events)
##changeCut(process, "selectedMuTauPairsAcoplanarity12", "dPhi12 > 2.967")
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define event selection criteria specific to tau id. efficiency measurement

process.uniqueTauCandidateCut = cms.EDFilter("BoolEventSelFlagProducer",
    pluginName = cms.string("uniqueTauCandidateCut"),
    pluginType = cms.string("PATCandViewMaxEventSelector"),
    src = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative'),
    maxNumber = cms.uint32(1)
)

process.selectZtoMuTauEvents._seq = process.selectZtoMuTauEvents._seq * process.uniqueTauCandidateCut

process.uniqueMuonCandidateCut = cms.EDFilter("BoolEventSelFlagProducer",
    pluginName = cms.string("uniqueMuonCandidateCut"),
    pluginType = cms.string("PATCandViewMaxEventSelector"),
    src = cms.InputTag('selectedLayer1MuonsGlobalCumulative'),
    maxNumber = cms.uint32(1)
)

process.selectZtoMuTauEvents._seq = process.selectZtoMuTauEvents._seq * process.uniqueMuonCandidateCut

process.selectedMuTauPairsBackToBack = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('selectedMuTauPairsPzetaDiffCumulative'),                                                             
    cut = cms.string('dPhi12 > 2.967'),
    filter = cms.bool(False)
)

process.diTauCandidateBackToBackCut = cms.EDFilter("BoolEventSelFlagProducer",
    pluginName = cms.string('diTauCandidateBackToBackCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('selectedMuTauPairsBackToBack'),
    minNumber = cms.uint32(1)
)

process.selectZtoMuTauEvents._seq = process.selectZtoMuTauEvents._seq * process.selectedMuTauPairsBackToBack * process.diTauCandidateBackToBackCut
#--------------------------------------------------------------------------------

# import config for event selection, event print-out and analysis sequence
# from Z --> muon + tau-jet cross-section analysis
from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *

muonHistManager.muonSource = cms.InputTag('selectedLayer1MuonsTrkIPcumulative')
tauHistManager.tauSource = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative')
diTauCandidateHistManagerForMuTau.pluginName = cms.string('diTauCandidateHistManagerForMuTau')
diTauCandidateHistManagerForMuTau.pluginType = cms.string('PATMuTauPairHistManager')
diTauCandidateHistManagerForMuTau.diTauCandidateSource = cms.InputTag('selectedMuTauPairsPzetaDiffCumulative')
diTauCandidateHistManagerForMuTau.visMassHypothesisSource = cms.InputTag('')

# import config for histogram manager specific to tau id. efficiency measurement
from TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauHistManager_cfi import *

tauIdEffBinning2regions = cms.PSet(
    name = cms.string("tauIdEffBinning2regions"),
    config = cms.VPSet(
        cms.PSet(
            extractor = cms.PSet(
                pluginType = cms.string("PATTauValExtractor"),
                src = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative'),
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

tauIdEffBinning4regions = cms.PSet(
    name = cms.string("tauIdEffBinning4regions"),
    config = cms.VPSet(
        cms.PSet(
            extractor = cms.PSet(
                pluginType = cms.string("PATTauValExtractor"),
                src = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative'),
                value = cms.string("tauID('ewkTauId')")
            ),
            branchName = cms.string('ewkTauId'),
            binning = cms.PSet(
                boundaries = cms.vdouble(0.5),
                min = cms.double(-0.01),
                max = cms.double(1.01)
            )
        ),
        cms.PSet(
            extractor = cms.PSet(
                pluginType = cms.string("PATMuTauPairChargeSignExtractor"),
                src = cms.InputTag('selectedMuTauPairsPzetaDiffCumulative')
            ),
            branchName = cms.string('diTauChargeSign'),
            binning = cms.PSet(
                boundaries = cms.vdouble(-0.5, +0.5),
                min = cms.double(-1.01),
                max = cms.double(+1.01)
            )
        )
    )
)

process.analyzeZtoMuTauEvents = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('TauIdEffAnalyzerZtoMuTau'), 
                            
    filters = cms.VPSet(
    
        #------------------------------------------------------------------------
        # selection criteria of Z --> muon + tau-jet cross-section analysis
        # (tau id. cuts and cut on muon + tau-jet charge disabled)
        #------------------------------------------------------------------------
    
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
        evtSelDiTauCandidateForMuTauAntiOverlapVeto,
        evtSelDiTauCandidateForMuTauZeroCharge,
        evtSelDiTauCandidateForMuTauAcoplanarity12,
        evtSelDiTauCandidateForMuTauMt1MET,
        evtSelDiTauCandidateForMuTauPzetaDiff,
        evtSelDiMuPairZmumuHypothesisVeto,

        #------------------------------------------------------------------------
        # selection criteria specific to tau id. efficiency measurement
        #------------------------------------------------------------------------

        # veto events with >= 2 muons
        cms.PSet(
            pluginName = cms.string('uniqueMuonCandidateCut'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('uniqueMuonCandidateCut')
        ),

        # veto events with >= 2 tau-jet candidates
        cms.PSet(
            pluginName = cms.string('uniqueTauCandidateCut'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('uniqueTauCandidateCut')
        ),

        # require muon and tau-jet to be back-to-back
        cms.PSet(
            pluginName = cms.string('diTauCandidateBackToBackCut'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateBackToBackCut')
        )
    ),
  
    analyzers = cms.VPSet(
        muonHistManager,
        tauHistManager,
        diTauCandidateHistManagerForMuTau,
        caloMEtHistManager,
        pfMEtHistManager,
        tauIdEffZtoMuTauHistManager,
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinner2regions'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinning2regions,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResults2regions')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManager2regions'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinning2regions,
            histManagers = cms.VPSet(
                muonHistManager,
                tauHistManager,
                diTauCandidateHistManagerForMuTau,
                tauIdEffZtoMuTauHistManager
            ),
            dqmDirectory_store = cms.string('tauIdEffHistograms2regions')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinner4regions'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinning4regions,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResults4regions')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManager4regions'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinning4regions,
            histManagers = cms.VPSet(
                muonHistManager,
                tauHistManager,
                diTauCandidateHistManagerForMuTau,
                tauIdEffZtoMuTauHistManager
            ),
            dqmDirectory_store = cms.string('tauIdEffHistograms4regions')
        )
    ),

    eventDumps = cms.VPSet(
        cms.PSet(
            pluginName = cms.string('muTauEventDump'),
            pluginType = cms.string('MuTauEventDump'),

            l1GtReadoutRecordSource = cms.InputTag('hltGtDigis::HLT'),
            l1GtObjectMapRecordSource = cms.InputTag('hltL1GtObjectMap::HLT'),
            l1BitsToPrint = cms.vstring('L1_SingleMu3', 'L1_SingleMu5', 'L1_SingleMu7', 'L1_SingleMu10', 'L1_SingleMu14'),
    
            hltResultsSource = cms.InputTag('TriggerResults::HLT'),
            hltPathsToPrint = cms.vstring('HLT_Mu15', 'HLT_IsoMu11'),
    
            genParticleSource = cms.InputTag('genParticles'),
            genJetSource = cms.InputTag('iterativeCone5GenJets'),
            genTauJetSource = cms.InputTag('tauGenJets'),
            
            electronSource = cms.InputTag('cleanLayer1Electrons'),
            muonSource = cms.InputTag('selectedLayer1MuonsTrkIPcumulative'),
            tauSource = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative'),
            diTauCandidateSource = cms.InputTag('selectedMuTauPairsPzetaDiffCumulative'),
            muTauZmumuHypothesisSource = cms.InputTag('muTauPairZmumuHypotheses'),
            diMuZmumuHypothesisSource = cms.InputTag('allDiMuPairZmumuHypotheses'),
            jetSource = cms.InputTag('selectedLayer1JetsEt20Cumulative'),
            caloMEtSource = cms.InputTag('layer1METs'),
            pfMEtSource = cms.InputTag('layer1PFMETs'),
            genMEtSource = cms.InputTag('genMETWithMu'),

            output = cms.string("std::cout"),
    
            triggerConditions = cms.vstring("uniqueTauCandidateCut: passed_cumulative")
        )
    ),
   
    analysisSequence = cms.VPSet(
    
        #------------------------------------------------------------------------
        # selection criteria of Z --> muon + tau-jet cross-section analysis
        # (tau id. cuts and cut on muon + tau-jet charge disabled)
        #------------------------------------------------------------------------
    
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

        # muon acceptance cuts
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

        # tau acceptance cuts
        cms.PSet(
            filter = cms.string('evtSelTauAntiOverlapWithMuonsVeto'),
            title = cms.string('Tau not overlapping w. Muon')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauEta'),
            title = cms.string('-2.1 < eta(Tau) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauPt'),
            title = cms.string('Pt(Tau) > 20 GeV')
        ),
    
        # selection of muon candidate (isolation & id.)
        # produced in muonic tau decay
        cms.PSet(
            analyzers = cms.vstring('muonHistManager')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonTrkIso'),
            title = cms.string('Muon Track iso.')
        ),
        cms.PSet(
            analyzers = cms.vstring('muonHistManager')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonEcalIso'),
            title = cms.string('Muon ECAL iso.')
        ),
        cms.PSet(
            analyzers = cms.vstring('muonHistManager')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonAntiPion'),
            title = cms.string('Muon pi-Veto')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonTrkIP'),
            title = cms.string('Muon Track IP')
        ),
    
        # selection of tau-jet candidate (id.)
        # produced in hadronic tau decay
        ##cms.PSet(
        ##    filter = cms.string('evtSelTauLeadTrk'),
        ##    title = cms.string('Tau lead. Track find.')
        ##),
        ##cms.PSet(
        ##    filter = cms.string('evtSelTauLeadTrkPt'),
        ##    title = cms.string('Tau lead. Track Pt')
        ##),
        ##cms.PSet(
        ##    filter = cms.string('evtSelTauTrkIso'),
        ##    title = cms.string('Tau Track iso.')
        ##),
        ##cms.PSet(
        ##    filter = cms.string('evtSelTauEcalIso'),
        ##    title = cms.string('Tau ECAL iso.')
        ##),
        ##cms.PSet(
        ##    filter = cms.string('evtSelTauProng'),
        ##    title = cms.string('Tau 1||3-Prong')
        ##),
        ##cms.PSet(
        ##    filter = cms.string('evtSelTauCharge'),
        ##    title = cms.string('Charge(Tau) = +/-1')
        ##),
        ##cms.PSet(
        ##    filter = cms.string('evtSelTauMuonVeto'),
        ##    title = cms.string('Tau mu-Veto')
        ##),

        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManager',
                'tauHistManager',
                'diTauCandidateHistManagerForMuTau'
            )
        ),

        # selection of muon + tau-jet combinations
        cms.PSet(
            filter = cms.string('evtSelDiTauCandidateForMuTauAntiOverlapVeto'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        ##cms.PSet(
        ##    filter = cms.string('evtSelDiTauCandidateForMuTauZeroCharge'),
        ##    title = cms.string('Charge(Muon+Tau) = 0')
        ##),
        cms.PSet(
            filter = cms.string('evtSelDiTauCandidateForMuTauAcoplanarity12'),
            title = cms.string('Acoplanarity(Muon+Tau)')
        ),
        ##cms.PSet(
        ##    filter = cms.string('evtSelDiTauCandidateForMuTauMt1MET'),
        ##    title = cms.string('M_{T}(Muon-MET) < 50 GeV')
        ##),
        ##cms.PSet(
        ##    filter = cms.string('evtSelDiTauCandidateForMuTauPzetaDiff'),
        ##    title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV')
        ##),
        ##cms.PSet(
        ##    filter = cms.string('evtSelDiMuPairZmumuHypothesisVeto'),
        ##    title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ##),

        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManager',
                'tauHistManager',
                'diTauCandidateHistManagerForMuTau'
            )
        ),

        #------------------------------------------------------------------------
        # selection criteria specific to tau id. efficiency measurement
        #------------------------------------------------------------------------

        # veto events with >= 2 muons
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCut'),
            title = cms.string('num. global Muons < 2')
        ),

        # veto events with >= 2 tau-jet candidates
        cms.PSet(
            filter = cms.string('uniqueTauCandidateCut'),
            title = cms.string('num. Tau-Jet Candidates < 2'),
            saveRunEventNumbers = cms.vstring('')
        ),

        # fill histograms & binning tables
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManager',
                'tauHistManager',
                'diTauCandidateHistManagerForMuTau',
                'caloMEtHistManager',
                'pfMEtHistManager',
                'tauIdEffZtoMuTauHistManager',
                'tauIdEffDataBinner2regions',
                'tauIdEffBinGridHistManager2regions',
                'tauIdEffDataBinner4regions',
                'tauIdEffBinGridHistManager4regions'
            )
        ),

        # require muon and tau-jet to be back-to-back
        # (cut added to end of selection criteria in order to be able to make plots of muon + tau-jet
        #  mass reconstructed via collinear approximation before acoplanarity cut is applied)
        cms.PSet(
            filter = cms.string('diTauCandidateBackToBackCut'),
            title = cms.string('dPhi(Muon,Tau) > 170 deg.'),
            saveRunEventNumbers = cms.vstring('')
        ),

        # fill histograms & binning tables
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManager',
                'tauHistManager',
                'diTauCandidateHistManagerForMuTau',
                'caloMEtHistManager',
                'pfMEtHistManager',
                'tauIdEffZtoMuTauHistManager',
                'tauIdEffDataBinner2regions',
                'tauIdEffBinGridHistManager2regions',
                'tauIdEffDataBinner4regions',
                'tauIdEffBinGridHistManager4regions'
            )
        )
    )
)

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system
#
#__process.source.fileNames = #inputFileNames#
#__process.maxEvents.input = cms.untracked.int32(#maxEvents#)
#__process.analyzeZtoMuTauEvents.filters[0] = copy.deepcopy(#genPhaseSpaceCut#)
#__process.saveZtoMuTauPlots.outputFileName = #plotsOutputFileName#
#
#--------------------------------------------------------------------------------

process.p = cms.Path(
   process.producePatTupleZtoMuTauSpecific
  + process.selectZtoMuTauEvents
  + process.analyzeZtoMuTauEvents
  + process.saveZtoMuTauPlots 
)

#--------------------------------------------------------------------------------
#
process.producePatTupleAll = cms.Sequence(process.producePatTuple + process.producePatTupleZtoMuTauSpecific)
#
# define "hook" for enabling/disabling production of PAT-tuple event content,
# depending on whether RECO/AOD or PAT-tuples are used as input for analysis
#
#__#patTupleProduction_line01#
#__#patTupleProduction_line02#
if not hasattr(process, 'batchMode'):
    process.p.replace(process.producePatTupleZtoMuTauSpecific, process.producePatTuple + process.producePatTupleZtoMuTauSpecific)
#--------------------------------------------------------------------------------

# print-out all python configuration parameter information
#print process.dumpPython()
