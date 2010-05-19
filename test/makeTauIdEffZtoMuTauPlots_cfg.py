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
process.loadTauIdEffZtoMuTau.inputFilePath = cms.string("rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_3_x/plots/ZtoMuTau_tauIdEff/7TeVii/")

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

process.dumpTauIdEffZtoMuTau = cms.EDAnalyzer("DQMDumpFilterStatisticsTables",
    dqmDirectories = cms.PSet(
        Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/FilterStatistics'),
        Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/FilterStatistics/'),
        WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/FilterStatistics/'),
        QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/FilterStatistics/'),
        TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/FilterStatistics')
    ),
    columnsSummaryTable = cms.vstring("Passed", "cumul. Efficiency", "margin. Efficiency", "indiv. Efficiency")
)

process.dumpTauIdEffZtoMuTauBinningResults2regions = cms.EDAnalyzer("DQMDumpBinningResults",
    binningService = cms.PSet(
        pluginType = cms.string("DataBinningService"),
        dqmDirectories = cms.PSet(
            Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResults2regions/'),
            Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResults2regions/'),
            WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResults2regions/'),
            QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResults2regions/'),
            TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau/tauIdEffBinningResults2regions/'),
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
