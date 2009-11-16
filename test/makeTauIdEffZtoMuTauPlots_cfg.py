import FWCore.ParameterSet.Config as cms

process = cms.Process('makeTauIdEffZtoMuTauPlots')

process.load("TauAnalysis.Configuration.dumpZtoMuTau_cff")
process.load("TauAnalysis.BgEstimationTools.plotTauIdEffZtoMuTau_cff")

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

process.dumpTauIdEffZtoMuTau = cms.EDAnalyzer("DQMDumpFilterStatisticsTables",
    dqmDirectories = cms.PSet(
        Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTau/FilterStatistics'),
        Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTau/FilterStatistics/'),
        WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTau/FilterStatistics/'),
        QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTau/FilterStatistics/'),
        TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTau/FilterStatistics')
    ),
    columnsSummaryTable = cms.vstring("Passed", "cumul. Efficiency", "margin. Efficiency", "indiv. Efficiency")
)

process.dumpTauIdEffZtoMuTau_binning = cms.EDAnalyzer("DQMDumpBinningResults",
    binningService = cms.PSet(
        pluginType = cms.string("DataBinningService"),
        dqmDirectories = cms.PSet(
            Ztautau = cms.string('harvested/Ztautau/TauIdEffAnalyzerZtoMuTau/afterUniqueTauCandidateCut/tauIdEffBinningResults4regions/'),
            Zmumu = cms.string('harvested/Zmumu/TauIdEffAnalyzerZtoMuTau/afterUniqueTauCandidateCut/tauIdEffBinningResults4regions/'),
            WplusJets = cms.string('harvested/WplusJets/TauIdEffAnalyzerZtoMuTau/afterUniqueTauCandidateCut/tauIdEffBinningResults4regions/'),
            QCD = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTau/afterUniqueTauCandidateCut/tauIdEffBinningResults4regions/'),
            TTplusJets = cms.string('harvested/TTplusJets/TauIdEffAnalyzerZtoMuTau/afterUniqueTauCandidateCut/tauIdEffBinningResults4regions/'),
        )
    )
)
 
process.makeTauIdEffZtoMuTauPlots = cms.Sequence(
    process.loadTauIdEffZtoMuTau
   + process.addTauIdEffZtoMuTau
   + process.saveTauIdEffZtoMuTau
   + process.dumpTauIdEffZtoMuTau + process.dumpTauIdEffZtoMuTau_binning
   + process.plotTauIdEffZtoMuTau + process.plotTauIdEffZtoMuTauShapes
)

process.p = cms.Path(process.makeTauIdEffZtoMuTauPlots)

# print-out all python configuration parameter information
#print process.dumpPython()
