import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
#
# Determine tau identification efficiency
# by fitting signal and background contributions
# in different regions via "generalized matrix" method
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

from TauAnalysis.DQMTools.plotterStyleDefinitions_cfi import *

process = cms.Process('fitGenMatrixTauIdEffZtoMuTau')

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

process.loadTauIdEffZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    smSum = cms.PSet(
        inputFileNames = cms.vstring(
            'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/plots/ZtoMuTau_tauIdEff/7TeV/plotsTauIdEffZtoMuTau_all.root'
        ),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

dqmDirectory_genMatrix3d = 'TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResultsGenMatrix3d/'
dqmDirectory_genMatrix4d = 'TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResultsGenMatrix4d/'

process.dumpTauIdEffZtoMuTauBinningResults = cms.EDAnalyzer("DQMDumpBinningResults",
    binningService = cms.PSet(
        dqmDirectories = cms.PSet(
            Ztautau = cms.string('harvested/Ztautau' + '/' + dqmDirectory_genMatrix4d),
            Zmumu = cms.string('harvested/Zmumu' + '/' + dqmDirectory_genMatrix4d),
            WplusJets = cms.string('harvested/WplusJets' + '/' + dqmDirectory_genMatrix4d),
            QCD = cms.string('harvested/qcdSum' + '/' + dqmDirectory_genMatrix4d),
            TTplusJets = cms.string('harvested/TTplusJets' + '/' + dqmDirectory_genMatrix4d)
        ),
        pluginType = cms.string("DataBinningService")
    )
)

process.dumpTauIdEffZtoMuTauGenMatrixProbabilities = cms.EDAnalyzer("DQMDumpGenMatrixProbabilities",
    binningService = cms.PSet(                                           
        dqmDirectories = cms.PSet(
            Ztautau = cms.string('harvested/Ztautau' + '/' + dqmDirectory_genMatrix4d),
            Zmumu = cms.string('harvested/Zmumu' + '/' + dqmDirectory_genMatrix4d),
            WplusJets = cms.string('harvested/WplusJets' + '/' + dqmDirectory_genMatrix4d),
            QCD = cms.string('harvested/qcdSum' + '/' + dqmDirectory_genMatrix4d),
            TTplusJets = cms.string('harvested/TTplusJets' + '/' + dqmDirectory_genMatrix4d)
        ),
        pluginType = cms.string("DataBinningService")
    )    
)

process.fitTauIdEffZtoMuTau = cms.EDAnalyzer("GenMatrixFit",
    processes = cms.PSet(
        Ztautau = cms.PSet(
            norm = cms.PSet(
                initial = cms.double(750.)
            ),
            par1 = cms.PSet(
                initial = cms.double(0.50)
            ),
            par2 = cms.PSet(
                initial = cms.double(0.50)
            ),
            par3 = cms.PSet(
                initial = cms.double(0.75)
            ),
            par4 = cms.PSet(
                initial = cms.double(0.90)
            ),
            drawOptions = drawOption_Ztautau
        ),
        WplusJets = cms.PSet(
            norm = cms.PSet(
                initial = cms.double(1500.)
            ),
            par1 = cms.PSet(
                initial = cms.double(0.06)
            ),
            par2 = cms.PSet(
                initial = cms.double(0.50)
            ),
            par3 = cms.PSet(
                initial = cms.double(0.25)
            ),
            par4 = cms.PSet(
                initial = cms.double(0.60)
            ),
            drawOptions = drawOption_WplusJets
        ),
        qcdSum = cms.PSet(
            norm = cms.PSet(
                initial = cms.double(2500.)
            ),
            par1 = cms.PSet(
                initial = cms.double(0.04)
            ),
            par2 = cms.PSet(
                initial = cms.double(0.05)
            ),
            par3 = cms.PSet(
                initial = cms.double(0.75)
            ),
            par4 = cms.PSet(
                initial = cms.double(0.40)
            ),
            drawOptions = drawOption_QCD
        )
    ),

    # use "pseudo" data-samples consisting of all Monte Carlo processes for testing                      
    data = cms.PSet(
        binningService = cms.PSet(
            pluginType = cms.string("DataBinningService")
        ),
        dqmDirectory = cms.string('harvested/smSum' + '/' + dqmDirectory_genMatrix4d)
    ),

    fit = cms.PSet(
        verbosity = cms.PSet(
            printLevel = cms.int32(1),
            printWarnings = cms.bool(True)
        )
    ),                                         

    output = cms.PSet(
        fitResults = cms.PSet(
            dqmDirectory = cms.string('fitGenMatrixTauIdEffZtoMuTau/fitResults')
        ),
        controlPlots = cms.PSet(
            fileName = cms.string("fitGenMatrixTauIdEffZtoMuTau.png")
        )
    )
)                          

process.saveFitResultsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('fitGenMatrixTauIdEffZtoMuTau_results.root'),
    outputCommands = cms.vstring(
        'drop harvested/*',
        'keep template/*',
        'keep harvested/*/zMuTauAnalyzer/afterEvtSelDiMuPairZmumuHypothesisVeto/*',
        'keep fitGenMatrixTauIdEffZtoMuTau/*'
    )
)

process.p = cms.Path(
    process.loadTauIdEffZtoMuTau
   + process.dumpDQMStore
   + process.dumpTauIdEffZtoMuTauBinningResults
   + process.dumpTauIdEffZtoMuTauGenMatrixProbabilities
   + process.fitTauIdEffZtoMuTau
   + process.saveFitResultsTauIdEffZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
