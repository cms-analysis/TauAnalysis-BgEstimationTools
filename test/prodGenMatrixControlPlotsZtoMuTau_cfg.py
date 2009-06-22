import FWCore.ParameterSet.Config as cms

process = cms.Process("prodGenMatrixControlPlotsZtoMuTau")

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff')
process.GlobalTag.globaltag = 'IDEAL_V9::All'

#--------------------------------------------------------------------------------
# import sequence for PAT-tuple production
process.load("TauAnalysis.Configuration.producePatTuple_cff")

# import sequence for event selection
process.load("TauAnalysis.Configuration.selectZtoMuTau_cff")

# import configuration parameters for submission of jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)
from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_cfi import *
from TauAnalysis.BgEstimationTools.bgEstSampleDefinitionsZtoMuTau_cfi import *

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
        #'rfio:/castor/cern.ch/user/v/veelken/CMSSW_2_2_3/muTauSkim.root'
        'file:/afs/cern.ch/user/v/veelken/scratch0/CMSSW_2_2_10/src/TauAnalysis/Configuration/test/muTauSkim.root'
    )
)

process.DQMStore = cms.Service("DQMStore")

process.saveZtoMuTauPlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('plotsZtoMuTau.root')
)

#--------------------------------------------------------------------------------
# define event selection flags specific
# to production ofcontrol plots for generalized matrix method

process.muonEcalIsoCutLooseIsolation = cms.EDFilter("BoolEventSelFlagProducer",
    selectors = cms.VPSet(
        cms.PSet(
            pluginName = cms.string("muonEcalIsoCutLooseIsolation"),
            pluginType = cms.string("PATCandViewMinEventSelector"),
            src = cms.InputTag('selectedLayer1MuonsEcalIsoLooseIsolationCumulative'),
            minNumber = cms.uint32(1),
            instanceName = cms.string('prodNtupleZtoMuTau')
        )
    )
)

process.muTauPairsLooseSelection = cms.EDProducer("PATMuTauPairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('selectedLayer1MuonsEcalIsoLooseIsolationCumulative'),
    srcLeg2 = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('layer1METs'),
    recoMode = cms.string(""),
    verbosity = cms.untracked.int32(0)
)

process.muTauPairCutLooseSelection = cms.EDProducer("BoolEventSelFlagProducer",
    selectors = cms.VPSet(
        cms.PSet(
            pluginName = cms.string("muTauPairCutLooseSelection"),
            pluginType = cms.string("PATCandViewMinEventSelector"),
            src = cms.InputTag('muTauPairsLooseSelection'),
            minNumber = cms.uint32(1),
            instanceName = cms.string('prodNtupleZtoMuTau')
        )
    )
)                                                                             

process.produceBoolEventSelFlags = cms.Sequence( process.muonEcalIsoCutLooseIsolation
                                                +process.muTauPairsLooseSelection + process.muTauPairCutLooseSelection )
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define analysis sequence

from TauAnalysis.BgEstimationTools.bgEstBinGridZtoMuTau_cfi import *
from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *

process.prodGenMatrixControlPlotsZtoMuTau = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('zMuTauGenMatrixControlPlots'), 
                            
    filters = cms.VPSet(
        # generator level phase-space selection
        # (NOTE: to be used in case of Monte Carlo samples
        #        overlapping in simulated phase-space only !!)
        cms.PSet(
            pluginName = cms.string('genPhaseSpaceCut'),
            pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
            src = cms.InputTag('genPhaseSpaceEventInfo'),
            cut = cms.string('')
        ),
        
        # reconstruction level event selection
        cms.PSet(
            pluginName = cms.string('genMatrixEventSelection'),
            pluginType = cms.string('MultiBoolEventSelFlagSelector'),
            flags = cms.VInputTag(
                cms.InputTag('Trigger'),
                cms.InputTag('primaryEventVertex'),
                cms.InputTag('primaryEventVertexQuality'),
                cms.InputTag('primaryEventVertexPosition'),
                cms.InputTag('muonEcalIsoCutLooseIsolation', 'prodNtupleZtoMuTau'),
                cms.InputTag('tauMuonVeto', 'cumulative'),                                                        
                cms.InputTag('muTauPairCutLooseSelection', 'prodNtupleZtoMuTau')
            )
        )
    ),
  
    analyzers = cms.VPSet(
        cms.PSet(
            pluginName = cms.string('muTauDataBinner'),
            pluginType = cms.string('DataBinner'),
            binning = genMatrixBinningZtoMuTau,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('genMatrixBinningResults')
        ),
        cms.PSet(
            pluginName = cms.string('muTauBinGridHistManager'),
            pluginType = cms.string('BinGridHistManager'),
            binning = genMatrixBinningZtoMuTau,
            histManagers = cms.VPSet(
                muonHistManager,
                tauHistManager,
                diTauCandidateHistManagerForMuTau,
                metHistManager
            ),
            dqmDirectory_store = cms.string('genMatrixBinningHistograms')
        )
    ),

    eventDumps = cms.VPSet(),
   
    analysisSequence = cms.VPSet(
        # generator level phase-space selection
        # (NOTE: to be used in case of Monte Carlo samples
        #        overlapping in simulated phase-space only !!)
        cms.PSet(
            filter = cms.string('genPhaseSpaceCut'),
            title = cms.string('gen. Phase-Space'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            filter = cms.string('genMatrixEventSelection'),
            title = cms.string('event selection'),
            saveRunEventNumbers = cms.vstring('')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muTauDataBinner',
                'muTauBinGridHistManager'
            )
        )
    )
)
#--------------------------------------------------------------------------------

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
switchToPFTauFixedCone(process)
#--------------------------------------------------------------------------------

process.p = cms.Path( process.producePatTuple
#                    +process.printEventContent    # uncomment to enable dump of event content after PAT-tuple production
                     +process.selectZtoMuTauEvents
                     +process.produceBoolEventSelFlags
                     +process.prodGenMatrixControlPlotsZtoMuTau
                     +process.saveZtoMuTauPlots )

# print-out all python configuration parameter information
#print process.dumpPython()
