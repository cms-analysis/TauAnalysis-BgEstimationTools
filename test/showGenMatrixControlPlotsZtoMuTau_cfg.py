import FWCore.ParameterSet.Config as cms

process = cms.Process('showGenMatrixControlPlotsZtoMuTau')

from TauAnalysis.Configuration.plotterProcessDefinitions_cfi import *

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

#--------------------------------------------------------------------------------
# Load histograms and binning results information from ROOT file
#--------------------------------------------------------------------------------

genMatrixZeePlusJets = copy.deepcopy(process_ZeePlusJets.config_dqmFileLoader)
genMatrixZeePlusJets.inputFileNames = cms.vstring('')

genMatrixZmumu = copy.deepcopy(process_Zmumu.config_dqmFileLoader)
genMatrixZmumu.inputFileNames = cms.vstring('')

genMatrixZmumuPlusJets = copy.deepcopy(process_ZmumuPlusJets.config_dqmFileLoader)
genMatrixZmumuPlusJets.inputFileNames = cms.vstring('')

genMatrixZtautau = copy.deepcopy(process_Ztautau.config_dqmFileLoader)
genMatrixZtautau.inputFileNames = cms.vstring('')

genMatrixZtautauPlusJets = copy.deepcopy(process_ZtautauPlusJets.config_dqmFileLoader)
genMatrixZtautauPlusJets.inputFileNames = cms.vstring('plotsZtoMuTau.root')

genMatrixWplusJets = copy.deepcopy(process_WplusJets.config_dqmFileLoader)
genMatrixWplusJets.inputFileNames = cms.vstring('')

genMatrixInclusivePPmuX = copy.deepcopy(process_InclusivePPmuX.config_dqmFileLoader)
genMatrixInclusivePPmuX.inputFileNames = cms.vstring('')

genMatrixPPmuXptGt20 = copy.deepcopy(process_PPmuXptGt20.config_dqmFileLoader)
genMatrixPPmuXptGt20.inputFileNames = cms.vstring('')

process.loadZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    #ZeePlusJets = genMatrixZeePlusJets,
    #Zmumu = genMatrixZmumu,
    #ZmumuPlusJets = genMatrixZmumuPlusJets,
    #Ztautau = genMatrixZtautau,
    ZtautauPlusJets = genMatrixZtautauPlusJets
    #WplusJets = genMatrixWplusJets,
    #InclusivePPmuX = genMatrixInclusivePPmuX,
    #PPmuXptGt20 = genMatrixPPmuXptGt20
)

#--------------------------------------------------------------------------------
# Print-out binning results information for Z --> mu + tau-jet channel
#--------------------------------------------------------------------------------

process.dumpZtoMuTau = cms.EDAnalyzer("DQMDumpBinningResults",
    binningService = cms.PSet(
        pluginType = cms.string("DataBinningService"),
        dqmDirectories = cms.PSet(
            #ZeePlusJets = cms.string('ZeePlusJets/zMuTauGenMatrixControlPlots/afterGenMatrixEventSelection/genMatrixBinningResults/'),
            #Zmumu = cms.string('Zmumu/zMuTauGenMatrixControlPlots/afterGenMatrixEventSelection/genMatrixBinningResults/'),
            #ZmumuPlusJets = cms.string('ZmumuPlusJets/zMuTauGenMatrixControlPlots/afterGenMatrixEventSelection/genMatrixBinningResults/'),
            #Ztautau = cms.string('Ztautau/zMuTauGenMatrixControlPlots/afterGenMatrixEventSelection/genMatrixBinningResults/'),
            ZtautauPlusJets = cms.string('ZtautauPlusJets/zMuTauGenMatrixControlPlots/afterGenMatrixEventSelection/genMatrixBinningResults/')
            #WplusJets = cms.string('WplusJets/zMuTauGenMatrixControlPlots/afterGenMatrixEventSelection/genMatrixBinningResults/'),
            #QCD = cms.string('qcdSum/zMuTauGenMatrixControlPlots/afterGenMatrixEventSelection/genMatrixBinningResults/')
       )
    )
)

process.showZtoMuTauPlots = cms.Sequence( process.loadZtoMuTau
                                         +process.dumpZtoMuTau )

process.p = cms.Path(process.showZtoMuTauPlots)

# print-out all python configuration parameter information
#print process.dumpPython()
