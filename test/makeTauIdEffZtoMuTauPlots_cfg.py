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
process.loadTauIdEffZtoMuTau.inputFilePath = cms.string("rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/plots/ZtoMuTau_tauIdEff/7TeV/")

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

process.dumpTauIdEffZtoMuTau_absMuonIsolation = cms.EDAnalyzer("DQMDumpFilterStatisticsTables",
    dqmDirectories = cms.PSet(
        Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/FilterStatistics'),
        Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/FilterStatistics/'),
        WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/FilterStatistics/'),
        QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/FilterStatistics/'),
        TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/FilterStatistics')
    ),
    columnsSummaryTable = cms.vstring("Passed", "cumul. Efficiency", "margin. Efficiency", "indiv. Efficiency")
)

process.dumpTauIdEffZtoMuTau_relMuonIsolation = cms.EDAnalyzer("DQMDumpFilterStatisticsTables",
    dqmDirectories = cms.PSet(
        Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTau_relMuonIsolation/FilterStatistics'),
        Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTau_relMuonIsolation/FilterStatistics/'),
        WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTau_relMuonIsolation/FilterStatistics/'),
        QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTau_relMuonIsolation/FilterStatistics/'),
        TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTau_relMuonIsolation/FilterStatistics')
    ),
    columnsSummaryTable = cms.vstring("Passed", "cumul. Efficiency", "margin. Efficiency", "indiv. Efficiency")
)

process.dumpTauIdEffZtoMuTau_genMatrixRelMuonIsolation = cms.EDAnalyzer("DQMDumpFilterStatisticsTables",
    dqmDirectories = cms.PSet(
        Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/FilterStatistics'),
        Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/FilterStatistics/'),
        WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/FilterStatistics/'),
        QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/FilterStatistics/'),
        TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/FilterStatistics')
    ),
    columnsSummaryTable = cms.vstring("Passed", "cumul. Efficiency", "margin. Efficiency", "indiv. Efficiency")
)

process.dumpTauIdEffZtoMuTau = cms.Sequence(
    process.dumpTauIdEffZtoMuTau_absMuonIsolation
   * process.dumpTauIdEffZtoMuTau_relMuonIsolation
   * process.dumpTauIdEffZtoMuTau_genMatrixRelMuonIsolation
)

process.dumpTauIdEffZtoMuTauBinningResults2regions = cms.EDAnalyzer("DQMDumpBinningResults",
    binningService = cms.PSet(
        pluginType = cms.string("DataBinningService"),
        dqmDirectories = cms.PSet(
            Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResults2regions/'),
            Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResults2regions/'),
            WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResults2regions/'),
            QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResults2regions/'),
            TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTau_genMatrixRelMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResults2regions/'),
        )
    )
)

process.reloadTauIdEffZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(        
        inputFileNames = cms.vstring('plotsTauIdEffZtoMuTau_all.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)
 
process.p = cms.Path(
    process.loadTauIdEffZtoMuTau
   + process.addTauIdEffZtoMuTau
    #process.reloadTauIdEffZtoMuTau
   + process.dumpDQMStore 
   + process.saveTauIdEffZtoMuTau
   + process.dumpTauIdEffZtoMuTau + process.dumpTauIdEffZtoMuTauBinningResults2regions
   + process.plotTauIdEffZtoMuTau + process.plotTauIdEffZtoMuTauShapes
)

# print-out all python configuration parameter information
#print process.dumpPython()
