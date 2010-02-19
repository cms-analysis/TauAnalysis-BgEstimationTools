import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('prodKineEventReweightsTauIdEffZtoMuTau')

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

process.saveTemplatesZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('bgEstTemplatesZtoMuTau.root')
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
#__process.saveTemplatesZtoMuTau.outputFileName = #plotsOutputFileName#
#__#batchMode#
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

process.load('TauAnalysis.BgEstimationTools.bgEstQCDenrichedSelection_cff')

# set generator level phase-space selection
# (to avoid overlap of different  Monte Carlo samples in simulated phase-space)
if hasattr(process, "batchMode"):
    process.analyzeEventsBgEstQCDenriched.filters[0] = getattr(process, "genPhaseSpaceCut")

# define additional event selection criteria
# to increase purity of QCD enriched background sample
process.muonsTauIdEffQCDenrichedLooseTrkIso = copy.deepcopy(selectedLayer1MuonsTrkIso)
process.muonsTauIdEffQCDenrichedLooseTrkIso.sumPtMax = cms.double(8.)

process.muonsTauIdEffQCDenrichedLooseEcalIso = copy.deepcopy(selectedLayer1MuonsEcalIso)
process.muonsTauIdEffQCDenrichedLooseEcalIso.cut = cms.string('userIsolation(1) < 8.')

process.muonsTauIdEffQCDenrichedPionVeto = copy.deepcopy(selectedLayer1MuonsPionVeto)

muonSelConfiguratorTauIdEffQCDenriched = objSelConfigurator(
    [ muonsTauIdEffQCDenrichedLooseTrkIso,
      muonsTauIdEffQCDenrichedLooseEcalIso,
      muonsTauIdEffQCDenrichedPionVeto ],
    src = "selectedLayer1MuonsPt15Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

process.muonsTauIdEffQCDenrichedTightTrkIso = copy.deepcopy(selectedLayer1MuonsTrkIso)
process.muonsTauIdEffQCDenrichedTightTrkIso.src = cms.InputTag('muonsTauIdEffQCDenrichedPionVetoCumulative')
process.muonsTauIdEffQCDenrichedTightTrkIso.sumPtMax = cms.double(1.)

process.muonsTauIdEffQCDenrichedTightEcalIso = copy.deepcopy(selectedLayer1MuonsTrkIso)
process.muonsTauIdEffQCDenrichedTightEcalIso.src = cms.InputTag('muonsTauIdEffQCDenrichedPionVetoCumulative')
process.muonsTauIdEffQCDenrichedTightEcalIso.cut = cms.string('userIsolation(1) < 1.')

process.selectMuonsTauIdEffQCDenriched = muonSelConfiguratorTauIdEffQCDenriched.configure(pyNameSpace = locals())
process.selectMuonsTauIdEffQCDenriched.seq_ = process.selectMuonsTauIdEffQCDenriched.seq_
                                             * process.muonsTauIdEffQCDenrichedTightTrkIso
                                             * process.muonsTauIdEffQCDenrichedTightEcalIso

process.muTauPairsTauIdEffQCDenriched = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsTauIdEffQCDenrichedPionVetoCumulative'),
    srcLeg2 = cms.InputTag('tausBgEstQCDenrichedMuonVetoCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('layer1METs'),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)
)

process.muTauPairsTauIdEffQCDenrichedMt1MET = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsTauIdEffQCDenriched'),                                                   
    cut = cms.string('mt1MET < 30.'),
    filter = cms.bool(False)
)

process.muTauPairsTauIdEffQCDenrichedNonZeroCharge = cms.EDFilter("PATMuTauPairSelector",
    src = cms.InputTag('muTauPairsTauIdEffQCDenrichedMt1MET'),                                                   
    cut = cms.string('abs(charge) > 0.5'),
    filter = cms.bool(False)
)

process.selectMuTauPairsTauIdEffQCDenriched = cms.Sequence(
    muTauPairsTauIdEffQCDenriched + muTauPairsTauIdEffQCDenrichedMt1MET + muTauPairsTauIdEffQCDenrichedNonZeroCharge
)

process.metTauIdEffQCDenrichedPt20 = cms.EDFilter("PATMETSelector",
    src = cms.InputTag("layer1METs"),                                 
    cut = cms.string('pt > 20.'),
    filter = cms.bool(False)
)

# produce boolean event selection flags
cfgMuonAntiPionCutTauIdEffQCDenrichedLooseIso = copy.deepcopy(cfgMuonAntiPionCut)
cfgMuonAntiPionCutTauIdEffQCDenrichedLooseIso.pluginName = cms.string('muonAntiPionCutTauIdEffQCDenrichedLooseIso')
cfgMuonAntiPionCutTauIdEffQCDenrichedLooseIso.src_cumulative = cms.InputTag('muonsTauIdEffQCDenrichedPionVetoCumulative')
cfgMuonAntiPionCutTauIdEffQCDenrichedLooseIso.systematics = cms.vstring()

cfgMuonTightTrkIsoCutTauIdEffQCDenriched = copy.deepcopy(cfgMuonTrkIsoCut)
cfgMuonTightTrkIsoCutTauIdEffQCDenriched.pluginName = cms.string('muonTightTrkIsoCutTauIdEffQCDenriched')
cfgMuonTightTrkIsoCutTauIdEffQCDenriched.src_cumulative = cms.InputTag('muonsTauIdEffQCDenrichedTightTrkIsoCumulative')
cfgMuonTightTrkIsoCutTauIdEffQCDenriched.systematics = cms.vstring()

cfgMuonTightEcalIsoCutTauIdEffQCDenriched = copy.deepcopy(cfgMuonEcalIsoCut)
cfgMuonTightEcalIsoCutTauIdEffQCDenriched.pluginName = cms.string('muonTightEcalIsoCutTauIdEffQCDenriched')
cfgMuonTightEcalIsoCutTauIdEffQCDenriched.src_cumulative = cms.InputTag('muonsTauIdEffQCDenrichedTightEcalIsoCumulative')
cfgMuonTightEcalIsoCutTauIdEffQCDenriched.systematics = cms.vstring()

cfgMuTauPairNonZeroChargeTauIdEffQCDenriched = cms.PSet(
    pluginName = cms.string('muTauPairNonZeroChargeTauIdEffQCDenriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffQCDenrichedNonZeroCharge'),
    minNumber = cms.uint32(1)
)

cfgMEtPt20TauIdEffQCDenriched = cms.PSet(
    pluginName = cms.string('metPt20TauIdEffQCDenriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('metTauIdEffQCDenrichedPt20'),
    minNumber = cms.uint32(1)
)

evtSelConfiguratorTauIdEffQCDenriched = eventSelFlagProdConfigurator(
    [ cfgMuonAntiPionCutTauIdEffQCDenrichedLooseIso,
      cfgMuonTightTrkIsoCutTauIdEffQCDenriched,
      cfgMuonTightEcalIsoCutTauIdEffQCDenriched,
      cfgMuTauPairNonZeroChargeTauIdEffQCDenriched,
      cfgMEtPt20TauIdEffQCDenriched ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

selectEventsTauIdEffQCDenriched = evtSelConfiguratorTauIdEffQCDenriched.configure()

# filter events based on bolean event selection flags

process.filterEventsTauIdEffQCDenrichedLooseMuonIso = cms.EDFilter("MultiBoolEventSelFlagFilter",
    flags = cms.VInputTag(
        cms.InputTag('genPhaseSpaceEventInfo'),
        cms.InputTag('Trigger'),
        cms.InputTag('primaryEventVertex'),
        cms.InputTag('primaryEventVertexQuality'),
        cms.InputTag('primaryEventVertexPosition'),
        cms.InputTag('muonsTauIdEffQCDenrichedPionVeto', 'cumulative'),
        cms.InputTag('tauMuonVetoBgEstQCDenriched', 'cumulative'),
        cms.InputTag('muTauPairNonZeroChargeTauIdEffQCDenriched'),
        cms.InputTag('metPt20TauIdEffQCDenriched')
    )
)

process.plotEventsTauIdEffQCDenrichedLooseMuonIso = cms.EDAnalyzer("MuonAnalyzer",
    muonSource = cms.InputTag('muonsTauIdEffQCDenrichedPionVeto'),
    vertexSource = cms.InputTag('selectedPrimaryVertexPosition'),
    jetSource = cms.InputTag('selectedLayer1JetsEt20Cumulative'),
    genParticleSource = cms.InputTag('genParticles'),
    dqmDirectory_store = cms.string('TauIdEffQCDenrichedLooseMuonIso/MuonQuantities'),
    requireGenMuonMatch = cms.bool(False),
    skipPdgIdsGenParticleMatch = cms.vint32(12, 14, 16),                                                               
    normalization = cms.string("events")
)

process.analysisPathLooseMuonIso = cms.Path(
    process.filterEventsTauIdEffQCDenrichedLooseMuonIso
   + process.plotEventsTauIdEffQCDenrichedLooseMuonIso
)

process.filterEventsTauIdEffQCDenrichedTightMuonTrkIso = cms.EDFilter("MultiBoolEventSelFlagFilter",
    flags = cms.VInputTag(cms.InputTag('muonTightTrkIsoCutTauIdEffQCDenriched'))
)

process.plotEventsTauIdEffQCDenrichedTightMuonTrkIso = copy.deepcopy(process.plotEventsTauIdEffQCDenrichedLooseMuonIso)
process.plotEventsTauIdEffQCDenrichedTightMuonTrkIso.muonSource = cms.InputTag('muonsTauIdEffQCDenrichedTightTrkIso')

process.analysisPathTightMuonTrkIso = cms.Path(
    process.filterEventsTauIdEffQCDenrichedLooseMuonIso
   + process.filterEventsTauIdEffQCDenrichedTightMuonTrkIso
   + process.plotEventsTauIdEffQCDenrichedTightMuonTrkIso
)

process.filterEventsTauIdEffQCDenrichedTightMuonEcalIso = cms.EDFilter("MultiBoolEventSelFlagFilter",
    flags = cms.VInputTag(cms.InputTag('muonTightEcalIsoCutTauIdEffQCDenriched'))
)

process.plotEventsTauIdEffQCDenrichedTightMuonEcalIso = copy.deepcopy(process.plotEventsTauIdEffQCDenrichedLooseMuonIso)
process.plotEventsTauIdEffQCDenrichedTightMuonEcalIso.muonSource = cms.InputTag('muonsTauIdEffQCDenrichedTightEcalIso')

process.analysisPathTightMuonEcalIso = cms.Path(
    process.filterEventsTauIdEffQCDenrichedLooseMuonIso
   + process.filterEventsTauIdEffQCDenrichedTightMuonEcalIso
   + process.plotEventsTauIdEffQCDenrichedTightMuonEcalIso
)

process.p = cms.Path(
   process.producePatTupleZtoMuTauSpecific
  + process.selectZtoMuTauEvents
  + process.selectMuonsTauIdEffQCDenriched
  + process.selectMuTauPairsTauIdEffQCDenriched
  + process.selectMEtTauIdEffQCDenrichedPt20
)

process.q = cms.Path(
   process.saveTemplatesZtoMuTau
)

process.schedule = cms.Schedule(
    process.p,
    process.analysisPathLooseMuonIso,
    process.analysisPathTightMuonTrkIso,
    process.analysisPathTightMuonEcalIso,
    process.q
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
if not hasattr(process, "batchMode"):
    process.p.replace(process.producePatTupleZtoMuTauSpecific, process.producePatTuple + process.producePatTupleZtoMuTauSpecific)
#--------------------------------------------------------------------------------

# print-out all python configuration parameter information
#print process.dumpPython()
