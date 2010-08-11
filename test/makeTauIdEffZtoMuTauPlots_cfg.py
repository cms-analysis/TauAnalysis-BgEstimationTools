import FWCore.ParameterSet.Config as cms

process = cms.Process('makeTauIdEffZtoMuTauPlots')

process.load("TauAnalysis.Configuration.dumpZtoMuTau_cff")
process.load("TauAnalysis.BgEstimationTools.plotTauIdEffZtoMuTau_cff")

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

# define directory from which .root files containing the histograms get loaded
process.loadTauIdEffZtoMuTau.inputFilePath = cms.string("rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/plots/ZtoMuTau_tauIdEff/7TeVii/")

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

process.dumpTauIdEffZtoMuTauTemplateFit = cms.EDAnalyzer("DQMDumpFilterStatisticsTables",
    dqmDirectories = cms.PSet(
        Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTauTemplateFit/FilterStatistics'),
        Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTauTemplateFit/FilterStatistics/'),
        WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTauTemplateFit/FilterStatistics/'),
        QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTauTemplateFit/FilterStatistics/'),
        TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTauTemplateFit/FilterStatistics')
    ),
    columnsSummaryTable = cms.vstring("Passed", "cumul. Efficiency", "margin. Efficiency", "indiv. Efficiency")
)

process.dumpTauIdEffZtoMuTauGenMatrixFit = cms.EDAnalyzer("DQMDumpFilterStatisticsTables",
    dqmDirectories = cms.PSet(
        Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTauGenMatrixFit/FilterStatistics'),
        Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTauGenMatrixFit/FilterStatistics/'),
        WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTauGenMatrixFit/FilterStatistics/'),
        QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTauGenMatrixFit/FilterStatistics/'),
        TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTauGenMatrixFit/FilterStatistics')
    ),
    columnsSummaryTable = cms.vstring("Passed", "cumul. Efficiency", "margin. Efficiency", "indiv. Efficiency")
)

process.dumpTauIdEffZtoMuTauCombinedFit = cms.EDAnalyzer("DQMDumpFilterStatisticsTables",
    dqmDirectories = cms.PSet(
        Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTauCombinedFit/FilterStatistics'),
        Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTauCombinedFit/FilterStatistics/'),
        WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTauCombinedFit/FilterStatistics/'),
        QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFit/FilterStatistics/'),
        TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTauCombinedFit/FilterStatistics')
    ),
    columnsSummaryTable = cms.vstring("Passed", "cumul. Efficiency", "margin. Efficiency", "indiv. Efficiency")
)

process.dumpTauIdEffZtoMuTau = cms.Sequence(
    process.dumpTauIdEffZtoMuTauTemplateFit
   * process.dumpTauIdEffZtoMuTauGenMatrixFit
   * process.dumpTauIdEffZtoMuTauCombinedFit
)

process.reloadTauIdEffZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(        
        inputFileNames = cms.vstring('plotsTauIdEffZtoMuTau_all.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)
 
process.p = cms.Path(
   # process.loadTauIdEffZtoMuTau
   #+ process.addTauIdEffZtoMuTau
    process.reloadTauIdEffZtoMuTau
   + process.dumpDQMStore 
   #+ process.saveTauIdEffZtoMuTau
   + process.dumpTauIdEffZtoMuTau
   + process.plotTauIdEffZtoMuTau + process.plotTauIdEffZtoMuTauShapes
)

# print-out all python configuration parameter information
#print process.dumpPython()
