import FWCore.ParameterSet.Config as cms

process = cms.Process("skimZtoElecTau")

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
# import sequence for PAT-tuple production
process.load("TauAnalysis.Configuration.producePatTuple_cff")

# import sequence for event selection
process.load("TauAnalysis.Configuration.selectZtoElecTau_cff")

# import configuration parameters for submission of jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)
from TauAnalysis.Configuration.recoSampleDefinitionsZtoElecTau_cfi import *
from TauAnalysis.BgEstimationTools.bgEstSampleDefinitionsZtoElecTau_cfi import *

# import event-content definition of products to be stored in patTuple
from TauAnalysis.Configuration.patTupleEventContent_cff import *
from TauAnalysis.Skimming.EventContent_cff import *
#--------------------------------------------------------------------------------

# print event content 
process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/v/veelken/scratch0/CMSSW_2_2_10/src/TauAnalysis/Configuration/test/muTauSkim.root'
    )
)

#--------------------------------------------------------------------------------
# preselect events considered for "template" and "generalized matrix"
# background estimation methods
process.genPhaseSpaceFilter = cms.EDFilter("EventSelPluginFilter",
    selector = cms.PSet(
        pluginName = cms.string('genPhaseSpaceCut'),
        pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
        src = cms.InputTag('genPhaseSpaceEventInfo'),
        cut = cms.string('')
    )
)

process.electronsBgEstPreselection = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag('selectedLayer1ElectronsEcalIsoLooseIsolationCumulative'),                                        
    cut = cms.string('gsfTrack.isNonnull'),
    filter = cms.bool(False)
)

process.electronTrkCutBgEstPreselection = cms.EDFilter("BoolEventSelFlagProducer",
    pluginName = cms.string("electronTrkCutBgEstPreselection"),
    pluginType = cms.string("PATCandViewMinEventSelector"),
    src = cms.InputTag('electronsBgEstPreselection'),
    minNumber = cms.uint32(1)
)

process.tauProngCutBgEstPreselection = cms.EDFilter("BoolEventSelFlagProducer",
    pluginName = cms.string("tauProngCutBgEstPreselection"),
    pluginType = cms.string("PATCandViewMinEventSelector"),
    src = cms.InputTag('selectedLayer1TausProngCumulative'),
    minNumber = cms.uint32(1)                                                
)

process.elecTauPairsBgEstPreselection = cms.EDProducer("PATElecTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('electronsBgEstPreselection'),
    srcLeg2 = cms.InputTag('selectedLayer1TausProngCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('layer1METs'),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)
)

process.elecTauPairCutBgEstPreselection = cms.EDFilter("BoolEventSelFlagProducer",
    pluginName = cms.string("elecTauPairCutBgEstPreselection"),
    pluginType = cms.string("PATCandViewMinEventSelector"),
    src = cms.InputTag('elecTauPairsBgEstPreselection'),
    minNumber = cms.uint32(1)
)                                                                             

process.produceBoolEventSelFlags = cms.Sequence(
    process.electronsBgEstPreselection + process.electronTrkCutBgEstPreselection
   +process.tauProngCutBgEstPreselection
   +process.elecTauPairsBgEstPreselection + process.elecTauPairCutBgEstPreselection
)

process.selectEventsByBoolEventSelFlags = cms.EDFilter("MultiBoolEventSelFlagFilter",
    flags = cms.VInputTag(
        #cms.InputTag('Trigger'),
        cms.InputTag('primaryEventVertex'),
        cms.InputTag('primaryEventVertexQuality'),
        cms.InputTag('primaryEventVertexPosition'),
        cms.InputTag('electronTrkCutBgEstPreselection'),
        cms.InputTag('tauProngCutBgEstPreselection'),
        cms.InputTag('elecTauPairCutBgEstPreselection')
    )
)
#--------------------------------------------------------------------------------

process.saveBgEstSample = cms.OutputModule("PoolOutputModule",
    tauAnalysisEventContent,                                        
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('bgEstSkimPath')
    ),
    fileName = cms.untracked.string('bgEstSample.root')
)

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system
#
#__process.source.fileNames = #inputFileNames#
#__process.maxEvents.input = cms.untracked.int32(#maxEvents#)
#__process.genPhaseSpaceFilter.selector = copy.deepcopy(#genPhaseSpaceCut#)
#__process.saveBgEstSample.fileName = #bgEstSampleOutputFileName#
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

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.bgEstSkimPath = cms.Path(
    process.producePatTuple
#  * process.printEventContent    # uncomment to enable dump of event content after PAT-tuple production
   * process.selectZtoElecTauEvents
   * process.genPhaseSpaceFilter
   * process.produceBoolEventSelFlags
   * process.selectEventsByBoolEventSelFlags
)

process.o = cms.EndPath( process.saveBgEstSample )

# print-out all python configuration parameter information
#print process.dumpPython()
