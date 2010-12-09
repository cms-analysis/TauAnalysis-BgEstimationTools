import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('runTauIdEffAnalysisZtoMuTau')

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
#process.MessageLogger.suppressInfo = cms.untracked.vstring()
process.MessageLogger.suppressWarning = cms.untracked.vstring("PATTriggerProducer",)
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START38_V14::All')

#--------------------------------------------------------------------------------
# import sequences for PAT-tuple production
process.load("TauAnalysis.Configuration.producePatTuple_cff")
process.load("TauAnalysis.Configuration.producePatTupleZtoMuTauSpecific_cff")

# import sequence for event selection
process.load("TauAnalysis.Configuration.selectZtoMuTau_cff")
process.load("TauAnalysis.RecoTools.filterDataQuality_cfi")

# import configuration parameters for submission of jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)
from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_cfi import *
from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_10TeV_cfi import *
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.saveZtoMuTauPlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('plotsTauIdEffZtoMuTau.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1) 
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0021/F405BC9A-525D-DF11-AB96-002618943811.root',
        #'/store/relval/CMSSW_3_6_1/RelValZTT/GEN-SIM-RECO/START36_V7-v1/0020/EE3E8F74-365D-DF11-AE3D-002618FDA211.root'
        'file:/data1/veelken/CMSSW_3_6_x/skims/Ztautau_1_1_sXK.root'
        #'file:/data1/veelken/CMSSW_3_8_x/skims/AHtoMuTau/selEvents_AHtoMuTau_woBtag_runs145762to148058_RECO.root'
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
#switchToPFTauShrinkingCone(process)
#switchToPFTauFixedCone(process)

# comment-out to take new HPS + TaNC combined tau id. algorithm
switchToPFTauHPSpTaNC(process)

# disable preselection on of pat::Taus
# (disabled also in TauAnalysis/RecoTools/python/patPFTauConfig_cfi.py ,
#  but re-enabled after switching tau collection)
process.cleanPatTaus.preselection = cms.string('')

# add "ewkTauId" flag
setattr(process.patTaus.tauIDSources, "ewkTauId", cms.InputTag('ewkTauId'))
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::Jets
from PhysicsTools.PatAlgos.tools.jetTools import *

# uncomment to replace caloJets by pfJets
switchJetCollection(process, jetCollection = cms.InputTag("ak5PFJets"), outputModule = '')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::METs
from TauAnalysis.Configuration.tools.metTools import *

# uncomment to add pfMET
# set Boolean swich to true in order to apply type-1 corrections
addPFMet(process, correct = False)

# uncomment to replace caloMET by pfMET in all di-tau objects
process.load("TauAnalysis.CandidateTools.diTauPairProductionAllKinds_cff")
replaceMETforDiTaus(process, cms.InputTag('patMETs'), cms.InputTag('patPFMETs'))
#--------------------------------------------------------------------------------

process.load('TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauSelectionTemplateFit_cff')
process.load('TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauSelectionGenMatrixFit_cff')
process.load('TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauSelectionCombinedFit_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauWplusJetsEnrichedSelection_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauTTplusJetsEnrichedSelection_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauZmumuEnrichedSelection_cff')
process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauQCDenrichedSelection_cff')

# produce event weight variable for correcting "bias"
# of muon |eta| distribution caused by cuts on muon track and ECAL isolation variables
# in QCD background events
##process.kineEventReweightTauIdEffQCD = cms.EDProducer("ObjValProducer",
##    config = cms.PSet(
##        pluginType = cms.string("KineEventReweightExtractor"),
##        weightLookupTable = cms.PSet(
##            fileName = cms.string(
##                'file:/afs/cern.ch/user/v/veelken/public/TauAnalysis/CMSSW_3_3_x/muonKineReweightsTauIdEffZtoMuTau.root'
##            ),
##            meName = cms.string('DQMData/tauIdEffKineEventReweights/QCDenrichedMuonCombIsoFactorized_data/MuonPtVsAbsEta')
##        ),
##        variables = cms.PSet(
##            x = cms.PSet(
##                pluginType = cms.string("PATMuTauPairValExtractor"),
##                src = cms.InputTag('muTauPairsBgEstQCDenriched'),
##                value = cms.string("abs(leg1.eta)"),
##                indices = cms.vuint32(0)
##            ),
##            y = cms.PSet(
##                pluginType = cms.string("PATMuTauPairValExtractor"),
##                src = cms.InputTag('muTauPairsBgEstQCDenriched'),
##                value = cms.string("leg1.pt"),
##                indices = cms.vuint32(0)
##            )
##        )
##    )
##)

# add another analysis sequence for producing QCD templates
# in which the events are reweighted in order to correct for "bias" of muon Pt and eta distributions
# caused by cuts on (absolute) muon track and ECAL isolation
##process.analyzeEventsBgEstQCDenriched_reweighted = copy.deepcopy(process.analyzeEventsBgEstQCDenriched)
##process.analyzeEventsBgEstQCDenriched_reweighted.name = cms.string('BgEstTemplateAnalyzer_QCDenriched_reweighted')
##setattr(process.analyzeEventsBgEstQCDenriched_reweighted, "eventWeightSource", cms.VInputTag("kineEventReweightTauIdEffQCD"))

# add analysis sequences for applying fake-rate weights
# to estimate contribution of W + jets and QCD backgrounds
# to SS 'control', tau id. passed region
##from TauAnalysis.BgEstimationTools.tools.fakeRateTools import enableFakeRates_runTauIdEffAnalysisZtoMuTau
##enableFakeRates_runTauIdEffAnalysisZtoMuTau(process)

#--------------------------------------------------------------------------------
# switch to pfMET in all analysis sequences 
if hasattr(process, "patPFMETs"):
    process.muTauPairsTauIdEffZtoMuTauTemplateFit.srcMET = cms.InputTag('patPFMETs')

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for changing cut values
from TauAnalysis.Configuration.tools.changeCut import changeCut

# change muon Pt threshold to 20 GeV
changeCut(process, "selectedPatMuonsPt10", "pt > 20.")
changeCut(process, "muonsBgEstWplusJetsEnrichedPt", "pt > 15.")

# disable cut on muon calo. + segment compatibility
# (check that muon calo. compatibility is not affected by pile-up before re-enabling this cut)
changeCut(process, "selectedPatMuonsPionVeto", -1000., attribute = "AntiPionCut")
changeCut(process, "selectedPatMuonsPionVetoLooseIsolation", -1000., attribute = "AntiPionCut")
changeCut(process, "muonsForTauIdEffZtoMuTauTemplateFitPionVeto", -1000., attribute = "AntiPionCut")
changeCut(process, "muonsForTauIdEffZtoMuTauCombinedFitPionVetoLooseIso", -1000., attribute = "AntiPionCut")


# change upper limit on tranverse impact parameter of muon track to 2mm
changeCut(process, "selectedPatMuonsTrkIP", 0.2, attribute = "IpMax")
changeCut(process, "muonsForTauIdEffZtoMuTauTemplateFitTrkIP", 0.2, attribute = "IpMax")
changeCut(process, "muonsForTauIdEffZtoMuTauCombinedFitTrkIPlooseIso", 0.2, attribute = "IpMax")

# change eta acceptance for tau-jets to |eta| < 2.3
changeCut(process, "selectedPatTausForMuTauEta21", "abs(eta) < 2.3")
changeCut(process, "selectedPatTausForMuTauEta21", "abs(eta) < 2.3")

# change cut on TaNC output in case using new HPS + TaNC combined tau id. algorithm
# CV: discriminators by Track/ECAL isolation do not exist for the combined tau id. algorithm
#     and need to be replaced by dummy cuts, in order to avoid run-time exceptions
changeCut(process, "selectedPatTausTrkIso", "tauID('byTaNCloose') > -1.")
changeCut(process, "selectedPatTausForMuTauTrkIso", "tauID('byTaNCloose') > -1.")
changeCut(process, "selectedPatTausEcalIso", "tauID('byTaNCloose') > -1.")
changeCut(process, "selectedPatTausForMuTauEcalIso", "tauID('byTaNCloose') > -1.")
changeCut(process, "selectedPatTausTaNCdiscr", "tauID('byTaNCmedium') > 0.5")
changeCut(process, "selectedPatTausForMuTauTaNCdiscr", "tauID('byTaNCmedium') > 0.5")
changeCut(process, "tausBgEstQCDenrichedTaNCdiscr", "tauID('byTaNCloose') > 0.5 & tauID('byTaNCmedium') < 0.5")
changeCut(process, "tausBgEstQCDenrichedTrkIso", "tauID('byTaNCloose') > -1.")
changeCut(process, "tausBgEstQCDenrichedEcalIso", "tauID('byTaNCloose') > -1.")
changeCut(process, "tausBgEstWplusJetsEnrichedTaNCdiscr", "tauID('byTaNCloose') > 0.5 & tauID('byTaNCmedium') < 0.5")
changeCut(process, "tausBgEstWplusJetsEnrichedTrkIso", "tauID('byTaNCloose') > -1.")
changeCut(process, "tausBgEstWplusJetsEnrichedEcalIso", "tauID('byTaNCloose') > -1.")
changeCut(process, "tausBgEstZmumuJetMisIdEnrichedTaNCdiscr", "tauID('byTaNCloose') > 0.5 & tauID('byTaNCmedium') < 0.5")
changeCut(process, "tausBgEstZmumuJetMisIdEnrichedTrkIso", "tauID('byTaNCloose') > -1.")
changeCut(process, "tausBgEstZmumuJetMisIdEnrichedEcalIso", "tauID('byTaNCloose') > -1.")

# change lower limit on separation required between muon and tau-jet to dR > 0.5
changeCut(process, "selectedMuTauPairsAntiOverlapVeto", "dR12 > 0.5")
changeCut(process, "selectedMuTauPairsAntiOverlapVetoLooseMuonIsolation", "dR12 > 0.5")
changeCut(process, "muTauPairsTauIdEffZtoMuTauTemplateFit", 0.5, attribute = "dRmin12")
changeCut(process, "muTauPairsTauIdEffZtoMuTauCombinedFit", 0.5, attribute = "dRmin12")
changeCut(process, "muTauPairsBgEstQCDenriched", 0.5, attribute = "dRmin12")
changeCut(process, "muTauPairsBgEstWplusJetsEnriched", 0.5, attribute = "dRmin12")
changeCut(process, "muTauPairsBgEstTTplusJetsEnriched", 0.5, attribute = "dRmin12") 
changeCut(process, "muTauPairsBgEstZmumuJetMisIdEnriched", 0.5, attribute = "dRmin12") 
changeCut(process, "muTauPairsBgEstZmumuMuonMisIdEnriched", 0.5, attribute = "dRmin12")

# change upper limit on muon + MET transverse mass to 40 GeV
changeCut(process, "selectedMuTauPairsMt1MET", "mt1MET < 40.")
changeCut(process, "selectedMuTauPairsMt1METlooseMuonIsolation", "mt1MET < 40.")
changeCut(process, "muTauPairsBgEstQCDenrichedMt1MET", "mt1MET < 40.")
changeCut(process, "muTauPairsBgEstZmumuJetMisIdEnrichedMt1MET", "mt1MET < 40.")

# disable cut on Pzeta variable
changeCut(process, "selectedMuTauPairsPzetaDiff", "(pZeta - 1.5*pZetaVis) > -1000.")
changeCut(process, "selectedMuTauPairsPzetaDiffLooseMuonIsolation", "(pZeta - 1.5*pZetaVis) > -1000.")
changeCut(process, "muTauPairsBgEstQCDenrichedPzetaDiff", "(pZeta - 1.5*pZetaVis) > -1000.")
changeCut(process, "muTauPairsBgEstZmumuJetMisIdEnrichedPzetaDiff", "(pZeta - 1.5*pZetaVis) > -1000.")

# disable b-tagging for now
#changeCut(process, "jetsBgEstTTplusJetsEnrichedEt40bTag", "bDiscriminator('trackCountingHighEffBJetTags') > -1000.")
#--------------------------------------------------------------------------------

process.p = cms.Path(
   process.producePatTupleZtoMuTauSpecific
  + process.selectZtoMuTauEvents
  + process.bgEstTauIdEffZtoMuTauTemplateFitAnalysisSequence
  + process.bgEstTauIdEffZtoMuTauGenMatrixFitAnalysisSequence
  + process.bgEstTauIdEffZtoMuTauCombinedFitAnalysisSequence
  + process.bgEstWplusJetsEnrichedAnalysisSequence
  + process.bgEstTTplusJetsEnrichedAnalysisSequence
  + process.bgEstZmumuEnrichedAnalysisSequence
  + process.bgEstQCDenrichedAnalysisSequence
  ##+ process.kineEventReweightTauIdEffQCD + process.analyzeEventsBgEstQCDenriched_reweighted
  + process.saveZtoMuTauPlots 
)

process.q = cms.Path(process.dataQualityFilters)

process.schedule = cms.Schedule(process.q, process.p)

#--------------------------------------------------------------------------------
# disable accessing generator level information
# if running on data
#from TauAnalysis.Configuration.tools.switchToData import switchToData
#switchToData(process)
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
