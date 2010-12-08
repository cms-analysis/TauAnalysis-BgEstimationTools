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
##process.loadTauIdEffZtoMuTau.inputFilePath = cms.string("rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/plots/ZtoMuTau_tauIdEff/7TeV/")
##process.loadTauIdEffZtoMuTau.inputFilePath = cms.string("file:/data1/veelken/CMSSW_3_8_x/plots/ZtoMuTau_tauIdEff/old")
process.loadTauIdEffZtoMuTau.inputFilePath = cms.string("file:/data1/veelken/CMSSW_3_8_x/plots/ZtoMuTau_tauIdEff/2010Nov01")

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

process.dumpTauIdEffZtoMuTauBinningResults2dCombinedFit = cms.EDAnalyzer("DQMDumpBinningResults",
    binningService = cms.PSet(
        pluginType = cms.string("DataBinningService"),
        dqmDirectories = cms.PSet(
            Ztautau = cms.string(
              'harvested/Ztautau/TauIdEffAnalyzerZtoMuTauCombinedFit/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResults2d/'
            ),
            Zmumu = cms.string(
              'harvested/Zmumu/TauIdEffAnalyzerZtoMuTauCombinedFit/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResults2d/'
            ),
            WplusJets = cms.string(
              'harvested/WplusJets/TauIdEffAnalyzerZtoMuTauCombinedFit/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResults2d/'
            ),
            QCD = cms.string(
              'harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFit/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResults2d/'
            ),
            TTplusJets = cms.string(
              'harvested/TTplusJets/TauIdEffAnalyzerZtoMuTauCombinedFit/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResults2d/'
            ),
        )
    )
)

process.dumpTauIdEffZtoMuTauBinningResults2dCombinedFitQCD = cms.EDAnalyzer("DQMDumpBinningResults",
    binningService = cms.PSet(
        pluginType = cms.string("DataBinningService"),
        dqmDirectories = cms.PSet(
            Ztautau = cms.string(
              'harvested/Ztautau/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResultsComb3dQCD/'
            ),
            Zmumu = cms.string(
              'harvested/Zmumu/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResultsComb3dQCD/'
            ),
            WplusJets = cms.string(
              'harvested/WplusJets/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResultsComb3dQCD/'
            ),
            QCD = cms.string(
              'harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResultsComb3dQCD/'
            ),
            TTplusJets = cms.string(
              'harvested/TTplusJets/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResultsComb3dQCD/'
            ),
        )
    )
)

process.dumpTauIdEffZtoMuTauBinningResults2dCombinedFitWplusJets = cms.EDAnalyzer("DQMDumpBinningResults",
    binningService = cms.PSet(
        pluginType = cms.string("DataBinningService"),
        dqmDirectories = cms.PSet(
            Ztautau = cms.string(
              'harvested/Ztautau/TauIdEffAnalyzerZtoMuTauCombinedFitWplusJets/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResultsComb2dWplusJets/'
            ),
            Zmumu = cms.string(
              'harvested/Zmumu/TauIdEffAnalyzerZtoMuTauCombinedFitWplusJets/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResultsComb2dWplusJets/'
            ),
            WplusJets = cms.string(
              'harvested/WplusJets/TauIdEffAnalyzerZtoMuTauCombinedFitWplusJets/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResultsComb2dWplusJets/'
            ),
            QCD = cms.string(
              'harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFitWplusJets/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResultsComb2dWplusJets/'
            ),
            TTplusJets = cms.string(
              'harvested/TTplusJets/TauIdEffAnalyzerZtoMuTauCombinedFitWplusJets/' \
             + 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit/tauIdEffBinningResultsComb2dWplusJets/'
            ),
        )
    )
)

process.dumpTauIdEffZtoMuTauBinningResults2d = cms.Sequence(
    process.dumpTauIdEffZtoMuTauBinningResults2dCombinedFit
   * process.dumpTauIdEffZtoMuTauBinningResults2dCombinedFitQCD
   * process.dumpTauIdEffZtoMuTauBinningResults2dCombinedFitWplusJets
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
   ## process.reloadTauIdEffZtoMuTau
   + process.dumpDQMStore 
   + process.saveTauIdEffZtoMuTau
   + process.dumpTauIdEffZtoMuTau + process.dumpTauIdEffZtoMuTauBinningResults2d
   + process.plotTauIdEffZtoMuTauTemplateFit + process.plotTauIdEffZtoMuTauShapes
)

# print-out all python configuration parameter information
#print process.dumpPython()
