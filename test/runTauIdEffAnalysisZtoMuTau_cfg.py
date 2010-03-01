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
process.GlobalTag.globaltag = cms.string('MC_31X_V9::All')

#--------------------------------------------------------------------------------
# import sequences for PAT-tuple production
process.load("TauAnalysis.Configuration.producePatTuple_cff")
process.load("TauAnalysis.Configuration.producePatTupleZtoMuTauSpecific_cff")

# import sequence for event selection
process.load("TauAnalysis.Configuration.selectZtoMuTau_cff")

# import configuration parameters for submission of jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)
from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_cfi import *
from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_10TeV_cfi import *
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.saveTauIdEffZtoMuTauPlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('plotsTauIdEffZtoMuTau.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) 
    #input = cms.untracked.int32(1000)    
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_3_1_2/RelValZTT/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/A4DD1FAE-B178-DE11-B608-001D09F24EAC.root',
        #'/store/relval/CMSSW_3_1_2/RelValZTT/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/9408B54D-CB78-DE11-9AEB-001D09F2503C.root'
        'rfio:/castor/cern.ch/user/l/lusito/SkimOctober09/ZtautauSkimMT314_3/muTauSkim_1.root',
        'rfio:/castor/cern.ch/user/l/lusito/SkimOctober09/ZtautauSkimMT314_3/muTauSkim_2.root'
    )
    #skipBadFiles = cms.untracked.bool(True) 
)

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system
#
#__process.source.fileNames = #inputFileNames#
#__process.maxEvents.input = cms.untracked.int32(#maxEvents#)
#__setattr(process, "genPhaseSpaceCut", copy.deepcopy(#genPhaseSpaceCut#))
#__process.saveTauIdEffZtoMuTauPlots.outputFileName = #plotsOutputFileName#
#__#isBatchMode#
#
#--------------------------------------------------------------------------------

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
#switchToPFTauFixedCone(process)
setattr(process.allLayer1Taus.tauIDSources, "ewkTauId", cms.InputTag('ewkTauId'))
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::Jets
from PhysicsTools.PatAlgos.tools.jetTools import *

# uncomment to replace caloJets by pfJets
switchJetCollection(process, "iterativeCone5PFJets")
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::METs
from TauAnalysis.Configuration.tools.metTools import *

# uncomment to add pfMET
# set Boolean swich to true in order to apply type-1 corrections
addPFMet(process, correct = False)

# uncomment to replace caloMET by pfMET in all di-tau objects
process.load("TauAnalysis.CandidateTools.diTauPairProductionAllKinds_cff")
replaceMETforDiTaus(process, cms.InputTag('layer1METs'), cms.InputTag('layer1PFMETs'))
#--------------------------------------------------------------------------------

process.load('TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauSelection_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauWplusJetsEnrichedSelection_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauTTplusJetsEnrichedSelection_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauZmumuEnrichedSelection_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauQCDenrichedSelection_cff')

# set generator level phase-space selection
# (to avoid overlap of different  Monte Carlo samples in simulated phase-space)
if hasattr(process, "isBatchMode"):
    process.analyzeEventsTauIdEffZtoMuTau.filters[0] = getattr(process, "genPhaseSpaceCut")
    process.analyzeEventsBgEstWplusJetsEnriched.filters[0] = getattr(process, "genPhaseSpaceCut")
    process.analyzeEventsBgEstTTplusJetsEnriched.filters[0] = getattr(process, "genPhaseSpaceCut")
    process.analyzeEventsBgEstZmumuJetMisIdEnriched.filters[0] = getattr(process, "genPhaseSpaceCut")
    process.analyzeEventsBgEstZmumuMuonMisIdEnriched.filters[0] = getattr(process, "genPhaseSpaceCut")
    process.analyzeEventsBgEstQCDenriched.filters[0] = getattr(process, "genPhaseSpaceCut")

# produce event weight variable for correcting "bias"
# of muon |eta| distribution caused by cuts on muon track and ECAL isolation variables
# in QCD background events
process.kineEventReweightTauIdEffQCD = cms.EDProducer("ObjValProducer",
    config = cms.PSet(
        pluginType = cms.string("KineEventReweightExtractor"),
        weightLookupTable = cms.PSet(
            fileName = cms.string(
                'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_3_x/kineEventReweights/muonKineReweightsTauIdEffZtoMuTau.root'
            ),
            meName = cms.string('DQMData/tauIdEffKineEventReweights/QCDenrichedMuonCombIsoFactorized_data/MuonPtVsAbsEta')
        ),
        variables = cms.PSet(
            x = cms.PSet(
                pluginType = cms.string("PATMuTauPairValExtractor"),
                src = cms.InputTag('muTauPairsBgEstQCDenriched'),
                value = cms.string("abs(leg1.eta)"),
                indices = cms.vuint32(0)
            ),
            y = cms.PSet(
                pluginType = cms.string("PATMuTauPairValExtractor"),
                src = cms.InputTag('muTauPairsBgEstQCDenriched'),
                value = cms.string("leg1.pt"),
                indices = cms.vuint32(0)
            )
        )
    )
)

# add another analysis sequence for producing QCD templates
# in which the events are reweighted in order to correct for "bias" of muon Pt and eta distributions
# caused by cuts on (absolute) muon track and ECAL isolation
process.analyzeEventsBgEstQCDenriched_reweighted = copy.deepcopy(process.analyzeEventsBgEstQCDenriched)
process.analyzeEventsBgEstQCDenriched_reweighted.name = cms.string('BgEstTemplateAnalyzer_QCDenriched_reweighted')
setattr(process.analyzeEventsBgEstQCDenriched_reweighted, "eventWeightSource", cms.VInputTag("kineEventReweightTauIdEffQCD"))

process.p = cms.Path(
   process.producePatTupleZtoMuTauSpecific
  + process.selectZtoMuTauEvents
  + process.bgEstTauIdEffZtoMuTauAnalysisSequence
  + process.bgEstWplusJetsEnrichedAnalysisSequence
  + process.bgEstTTplusJetsEnrichedAnalysisSequence
  + process.bgEstZmumuEnrichedAnalysisSequence
  + process.bgEstQCDenrichedAnalysisSequence
  + process.kineEventReweightTauIdEffQCD + process.analyzeEventsBgEstQCDenriched_reweighted
  + process.saveTauIdEffZtoMuTauPlots 
)

#--------------------------------------------------------------------------------
# disable estimation of systematic uncertainties
from TauAnalysis.Configuration.tools.sysUncertaintyTools import disableSysUncertainties_runZtoMuTau
#
disableSysUncertainties_runZtoMuTau(process)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
process.producePatTupleAll = cms.Sequence(process.producePatTuple + process.producePatTupleZtoMuTauSpecific)
#
# define "hook" for enabling/disabling production of PAT-tuple event content,
# depending on whether RECO/AOD or PAT-tuples are used as input for analysis
#
#__#patTupleProduction#
if not hasattr(process, "isBatchMode"):
    process.p.replace(process.producePatTupleZtoMuTauSpecific, process.producePatTuple + process.producePatTupleZtoMuTauSpecific)
#--------------------------------------------------------------------------------

# print-out all python configuration parameter information
#print process.dumpPython()
