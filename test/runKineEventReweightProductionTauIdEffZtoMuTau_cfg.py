import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('runKineEventReweightProductionTauIdEffZtoMuTau')

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

process.saveTauIdEffKineReweightHistogramsZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('tauIdEffKineReweightHistogramsZtoMuTau.root')
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
#__process.saveTauIdEffKineReweightHistogramsZtoMuTau.outputFileName = #plotsOutputFileName#
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

process.load('TauAnalysis.BgEstimationTools.bgEstZtoMuTauQCDenrichedSelection_cff')

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

# set generator level phase-space selection
# (to avoid overlap of different  Monte Carlo samples in simulated phase-space)
process.genPhaseSpace = cms.EDProducer("BoolEventSelFlagProducer",
    pluginName = cms.string('genPhaseSpace'),
    pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
    src = cms.InputTag('genPhaseSpaceEventInfo'),
    cut = cms.string('')
)

if hasattr(process, "isBatchMode"):
    process.analyzeEventsBgEstQCDenriched.filters[0] = getattr(process, "genPhaseSpaceCut")

    genPhaseSpaceCut_instance = copy.deepcopy(getattr(process, "genPhaseSpaceCut"))
    setattr(genPhaseSpaceCut_instance, "instanceName", cms.string(""))
    process.genPhaseSpace = cms.EDProducer("BoolEventSelFlagProducer",
        selectors = cms.VPSet(genPhaseSpaceCut_instance)
    )                              

# define additional event selection criteria
# to increase purity of QCD enriched background sample
process.muonsTauIdEffQCDenrichedLooseTrkIso = process.selectedLayer1MuonsTrkIso.clone(
    sumPtMax = cms.double(8.)
)    
process.muonsTauIdEffQCDenrichedLooseTrkIso.sumPtMax = cms.double(8.)

process.muonsTauIdEffQCDenrichedLooseEcalIso = process.selectedLayer1MuonsEcalIso.clone(
    cut = cms.string('userIsolation("pat::EcalIso") < 8.')
)

process.muonsTauIdEffQCDenrichedPionVeto = copy.deepcopy(process.selectedLayer1MuonsPionVeto)

muonSelConfiguratorTauIdEffQCDenriched = objSelConfigurator(
    [ process.muonsTauIdEffQCDenrichedLooseTrkIso,
      process.muonsTauIdEffQCDenrichedLooseEcalIso,
      process.muonsTauIdEffQCDenrichedPionVeto ],
    src = "selectedLayer1MuonsPt15Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

process.selectMuonsTauIdEffQCDenrichedLooseIso = muonSelConfiguratorTauIdEffQCDenriched.configure(process = process)

#--------------------------------------------------------------------------------
# now apply muon track and ECAL isolation cuts in order to determine bias
# of muon Pt and eta distributions caused by isolation cuts
#
# NOTE:
#     (1) need to separately apply either muon track isolation or muon ECAL isolation cuts, but **not** both,
#         as too much QCD events get cut otherwise and the purity of the QCD background enriched sample drops
#         to 50% (from 80% in case only the muon track isolation cut is applied)
#     (2) determine bias on muon Pt and eta distributions resulting from **combination**
#         of muon track and ECAL isolation cuts by assuming that biases are independent,
#         i.e.
#          binContent(muon Pt/eta shape after muon track && ECAL isolation) =
#            binContent(shape after track isolation) * binContent(shape after ECAL isolation)/binContent(shape before isolation)^2
#     (3) as only shapes are relevant to correct for the bias, one does **not** need to assume
#         that the probabilities for QCD background events
#         to pass the muon track and ECAL isolation cuts are independent (only that the biases are),
#         which is clearly not the case (the probability for a muon to pass the ECAL isolation cut
#         is clearly higher for a muon that passes the track isolation cut than for muon which fails the track isolation cut)
#--------------------------------------------------------------------------------

process.muonsTauIdEffQCDenrichedTightTrkIso = copy.deepcopy(process.selectedLayer1MuonsTrkIso)
process.muonsTauIdEffQCDenrichedTightTrkIso.src = cms.InputTag('muonsTauIdEffQCDenrichedPionVetoCumulative')
process.muonsTauIdEffQCDenrichedTightTrkIso.sumPtMax = cms.double(1.)

process.muonsTauIdEffQCDenrichedTightEcalIso = copy.deepcopy(process.selectedLayer1MuonsEcalIso)
process.muonsTauIdEffQCDenrichedTightEcalIso.src = cms.InputTag('muonsTauIdEffQCDenrichedPionVetoCumulative')
process.muonsTauIdEffQCDenrichedTightEcalIso.cut = cms.string('userIsolation("pat::EcalIso") < 1.')

process.muonsTauIdEffQCDenrichedTightCombIso = copy.deepcopy(process.selectedLayer1MuonsEcalIso)
process.muonsTauIdEffQCDenrichedTightCombIso.src = cms.InputTag('muonsTauIdEffQCDenrichedTightTrkIso')
process.muonsTauIdEffQCDenrichedTightCombIso.cut = cms.string('userIsolation("pat::EcalIso") < 1.')

process.selectMuonsTauIdEffQCDenriched = cms.Sequence(
    process.selectMuonsTauIdEffQCDenrichedLooseIso
   + process.muonsTauIdEffQCDenrichedTightTrkIso + process.muonsTauIdEffQCDenrichedTightEcalIso
   + process.muonsTauIdEffQCDenrichedTightCombIso
)

process.muTauPairsTauIdEffQCDenriched = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsTauIdEffQCDenrichedPionVetoCumulative'),
    srcLeg2 = cms.InputTag('tausBgEstQCDenrichedMuonVetoCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('layer1METs'),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                                   
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
    process.muTauPairsTauIdEffQCDenriched
   + process.muTauPairsTauIdEffQCDenrichedMt1MET + process.muTauPairsTauIdEffQCDenrichedNonZeroCharge
)

process.metTauIdEffQCDenrichedPt20 = cms.EDFilter("PATMETSelector",
    src = cms.InputTag("layer1METs"),                                 
    cut = cms.string('pt > 20.'),
    filter = cms.bool(False)
)

process.selectMEtTauIdEffQCDenriched = cms.Sequence(process.metTauIdEffQCDenrichedPt20)

# produce boolean event selection flags
cfgMuonAntiPionCutTauIdEffQCDenrichedLooseIso = copy.deepcopy(process.cfgMuonAntiPionCut)
cfgMuonAntiPionCutTauIdEffQCDenrichedLooseIso.pluginName = cms.string('muonAntiPionCutTauIdEffQCDenrichedLooseIso')
cfgMuonAntiPionCutTauIdEffQCDenrichedLooseIso.src_cumulative = cms.InputTag('muonsTauIdEffQCDenrichedPionVetoCumulative')
cfgMuonAntiPionCutTauIdEffQCDenrichedLooseIso.systematics = cms.vstring()

cfgMuonTightTrkIsoCutTauIdEffQCDenriched = copy.deepcopy(process.cfgMuonTrkIsoCut)
cfgMuonTightTrkIsoCutTauIdEffQCDenriched.pluginName = cms.string('muonTightTrkIsoCutTauIdEffQCDenriched')
cfgMuonTightTrkIsoCutTauIdEffQCDenriched.src_cumulative = cms.InputTag('muonsTauIdEffQCDenrichedTightTrkIso')
cfgMuonTightTrkIsoCutTauIdEffQCDenriched.systematics = cms.vstring()

cfgMuonTightEcalIsoCutTauIdEffQCDenriched = copy.deepcopy(process.cfgMuonEcalIsoCut)
cfgMuonTightEcalIsoCutTauIdEffQCDenriched.pluginName = cms.string('muonTightEcalIsoCutTauIdEffQCDenriched')
cfgMuonTightEcalIsoCutTauIdEffQCDenriched.src_cumulative = cms.InputTag('muonsTauIdEffQCDenrichedTightEcalIso')
cfgMuonTightEcalIsoCutTauIdEffQCDenriched.systematics = cms.vstring()

cfgMuonTightCombIsoCutTauIdEffQCDenriched = copy.deepcopy(process.cfgMuonEcalIsoCut)
cfgMuonTightCombIsoCutTauIdEffQCDenriched.pluginName = cms.string('muonTightCombIsoCutTauIdEffQCDenriched')
cfgMuonTightCombIsoCutTauIdEffQCDenriched.src_cumulative = cms.InputTag('muonsTauIdEffQCDenrichedTightCombIso')
cfgMuonTightCombIsoCutTauIdEffQCDenriched.systematics = cms.vstring()

cfgMuTauPairNonZeroChargeTauIdEffQCDenriched = cms.PSet(
    pluginName = cms.string('muTauPairNonZeroChargeTauIdEffQCDenriched'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('muTauPairsTauIdEffQCDenrichedNonZeroCharge'),
    minNumber = cms.uint32(1)
)

cfgMEtPt20TauIdEffQCDenriched = cms.PSet(
    pluginName = cms.string('metPt20TauIdEffQCDenriched'),
    pluginType = cms.string('PATCandViewMaxEventSelector'),
    src = cms.InputTag('metTauIdEffQCDenrichedPt20'),
    maxNumber = cms.uint32(0)
)

evtSelConfiguratorTauIdEffQCDenriched = eventSelFlagProdConfigurator(
    [ cfgMuonAntiPionCutTauIdEffQCDenrichedLooseIso,
      cfgMuonTightTrkIsoCutTauIdEffQCDenriched,
      cfgMuonTightEcalIsoCutTauIdEffQCDenriched,
      cfgMuonTightCombIsoCutTauIdEffQCDenriched,
      cfgMuTauPairNonZeroChargeTauIdEffQCDenriched,
      cfgMEtPt20TauIdEffQCDenriched ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

process.selectEventsTauIdEffQCDenriched = evtSelConfiguratorTauIdEffQCDenriched.configure(process = process)

# filter events based on bolean event selection flags

process.filterEventsTauIdEffQCDenrichedLooseMuonIso = cms.EDFilter("MultiBoolEventSelFlagFilter",
    flags = cms.VInputTag(
        cms.InputTag('genPhaseSpace'),
        cms.InputTag('Trigger'),
        cms.InputTag('primaryEventVertex'),
        cms.InputTag('primaryEventVertexQuality'),
        cms.InputTag('primaryEventVertexPosition'),
        cms.InputTag('muonAntiPionCutTauIdEffQCDenrichedLooseIso', 'cumulative'),
        cms.InputTag('tauMuonVetoBgEstQCDenriched', 'cumulative'),
        cms.InputTag('muTauPairNonZeroChargeTauIdEffQCDenriched'),
        cms.InputTag('metPt20TauIdEffQCDenriched')
    )
)

process.plotEventsTauIdEffQCDenrichedLooseMuonIso = cms.EDAnalyzer("TauIdEffZtoMuTauAnalyzer",
    muonSource = cms.InputTag('muonsTauIdEffQCDenrichedPionVetoCumulative'),
    tauSource = cms.InputTag('tausBgEstQCDenrichedMuonVetoCumulative'),
    diTauSource = cms.InputTag('muTauPairsTauIdEffQCDenrichedNonZeroCharge'),
    centralJetSource = cms.InputTag('selectedLayer1JetsEt20Cumulative'),       
    dqmDirectory_store = cms.string('TauIdEffQCDenrichedLooseMuonIso/TauIdEffSpecificQuantities'),
    tauIdDiscriminator = cms.string("ewkTauId"),
    diTauChargeSignExtractor = cms.PSet(
        pluginType = cms.string("PATMuTauPairChargeSignExtractor"),
        src = cms.InputTag('muTauPairsTauIdEffQCDenrichedNonZeroCharge')
    )
)

process.analysisPathLooseMuonIso = cms.Path(
    process.filterEventsTauIdEffQCDenrichedLooseMuonIso
   + process.plotEventsTauIdEffQCDenrichedLooseMuonIso
)

process.filterEventsTauIdEffQCDenrichedTightMuonTrkIso = cms.EDFilter("MultiBoolEventSelFlagFilter",
    flags = cms.VInputTag(cms.InputTag('muonTightTrkIsoCutTauIdEffQCDenriched', 'cumulative'))
)

process.plotEventsTauIdEffQCDenrichedTightMuonTrkIso = process.plotEventsTauIdEffQCDenrichedLooseMuonIso.clone(
    muonSource = cms.InputTag('muonsTauIdEffQCDenrichedTightTrkIso'),
    dqmDirectory_store = cms.string('TauIdEffQCDenrichedTightMuonTrkIso/TauIdEffSpecificQuantities')
)

process.analysisPathTightMuonTrkIso = cms.Path(
    process.filterEventsTauIdEffQCDenrichedLooseMuonIso
   + process.filterEventsTauIdEffQCDenrichedTightMuonTrkIso
   + process.plotEventsTauIdEffQCDenrichedTightMuonTrkIso
)

process.filterEventsTauIdEffQCDenrichedTightMuonEcalIso = cms.EDFilter("MultiBoolEventSelFlagFilter",
    flags = cms.VInputTag(cms.InputTag('muonTightEcalIsoCutTauIdEffQCDenriched', 'cumulative'))
)

process.plotEventsTauIdEffQCDenrichedTightMuonEcalIso = process.plotEventsTauIdEffQCDenrichedLooseMuonIso.clone(
    muonSource = cms.InputTag('muonsTauIdEffQCDenrichedTightEcalIso'),
    dqmDirectory_store = cms.string('TauIdEffQCDenrichedTightMuonEcalIso/TauIdEffSpecificQuantities')
)

process.analysisPathTightMuonEcalIso = cms.Path(
    process.filterEventsTauIdEffQCDenrichedLooseMuonIso
   + process.filterEventsTauIdEffQCDenrichedTightMuonEcalIso
   + process.plotEventsTauIdEffQCDenrichedTightMuonEcalIso
)

process.filterEventsTauIdEffQCDenrichedTightMuonCombIso = cms.EDFilter("MultiBoolEventSelFlagFilter",
    flags = cms.VInputTag(cms.InputTag('muonTightCombIsoCutTauIdEffQCDenriched', 'cumulative'))
)

process.plotEventsTauIdEffQCDenrichedTightMuonCombIso = process.plotEventsTauIdEffQCDenrichedLooseMuonIso.clone(
    muonSource = cms.InputTag('muonsTauIdEffQCDenrichedTightCombIso'),
    dqmDirectory_store = cms.string('TauIdEffQCDenrichedTightMuonCombIso/TauIdEffSpecificQuantities')
)

process.analysisPathTightMuonCombIso = cms.Path(
    process.filterEventsTauIdEffQCDenrichedLooseMuonIso
   + process.filterEventsTauIdEffQCDenrichedTightMuonCombIso
   + process.plotEventsTauIdEffQCDenrichedTightMuonCombIso
)

process.p = cms.Path(
   process.producePatTupleZtoMuTauSpecific
  + process.selectZtoMuTauEvents   
  + process.selectMuonsBgEstQCDenriched + process.selectMuonsTauIdEffQCDenriched
  + process.selectTausBgEstQCDenriched
  + process.selectMuTauPairsBgEstQCDenriched + process.selectMuTauPairsTauIdEffQCDenriched
  + process.selectMEtTauIdEffQCDenriched
  + process.genPhaseSpace + process.selectEventsBgEstQCDenriched + process.selectEventsTauIdEffQCDenriched
)

process.q = cms.Path(
   process.saveTauIdEffKineReweightHistogramsZtoMuTau
)

process.schedule = cms.Schedule(
    process.p,
    process.analysisPathLooseMuonIso,
    process.analysisPathTightMuonTrkIso,
    process.analysisPathTightMuonEcalIso,
    process.analysisPathTightMuonCombIso,
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
if not hasattr(process, "isBatchMode"):
    process.p.replace(process.producePatTupleZtoMuTauSpecific, process.producePatTuple + process.producePatTupleZtoMuTauSpecific)
#--------------------------------------------------------------------------------

# print-out all python configuration parameter information
print process.dumpPython()
