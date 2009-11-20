import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
#
# Plot histograms for tau id. efficiency measurement in Z --> mu + tau-jet events
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

from TauAnalysis.BgEstimationTools.plotTauIdEffZtoMuTau_processes_cfi import *
from TauAnalysis.BgEstimationTools.plotTauIdEffZtoMuTau_drawJobs_cfi import *
from TauAnalysis.Configuration.plotZtoMuTau_cff import plotZtoMuTau
from TauAnalysis.DQMTools.plotterStyleDefinitions_cfi import *

loadTauIdEffZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    Ztautau = copy.deepcopy(processZtoMuTau_ZtautauSum.config_dqmFileLoader),
    Zmumu = copy.deepcopy(processZtoMuTau_ZmumuSum.config_dqmFileLoader),
    WplusJets = copy.deepcopy(processZtoMuTau_WplusJetsSum.config_dqmFileLoader),
    InclusivePPmuX = copy.deepcopy(processZtoMuTau_InclusivePPmuX.config_dqmFileLoader),
    PPmuXptGt20 = copy.deepcopy(processZtoMuTau_PPmuXptGt20Sum.config_dqmFileLoader),
    TTplusJets = copy.deepcopy(processZtoMuTau_TTplusJetsSum.config_dqmFileLoader),
    inputFilePath = cms.string("rfio:/castor/cern.ch/user/v/veelken/plots/TauIdEffV/")
)

addTauIdEffZtoMuTau_qcdSum = cms.EDAnalyzer("DQMHistAdder",
    qcdSum = cms.PSet(
        dqmDirectories_input = cms.vstring(
            #'harvested/InclusivePPmuX/TauIdEffAnalyzerZtoMuTau',
            'harvested/PPmuXptGt20/TauIdEffAnalyzerZtoMuTau'
        ),
        dqmDirectory_output = cms.string('harvested/qcdSum/TauIdEffAnalyzerZtoMuTau')
    )                          
)

addTauIdEffZtoMuTau_smSum = cms.EDAnalyzer("DQMHistAdder",
    smSum = cms.PSet(
        dqmDirectories_input = cms.vstring(
            'harvested/Ztautau/TauIdEffAnalyzerZtoMuTau',
            'harvested/Zmumu/TauIdEffAnalyzerZtoMuTau',
            'harvested/WplusJets/TauIdEffAnalyzerZtoMuTau',
            'harvested/TTplusJets/TauIdEffAnalyzerZtoMuTau',
            'harvested/qcdSum/TauIdEffAnalyzerZtoMuTau'
        ),
        dqmDirectory_output = cms.string('harvested/smSum/TauIdEffAnalyzerZtoMuTau')
    )
)

addTauIdEffZtoMuTau = cms.Sequence(addTauIdEffZtoMuTau_qcdSum + addTauIdEffZtoMuTau_smSum)

plotTauIdEffZtoMuTau = copy.deepcopy(plotZtoMuTau)
plotTauIdEffZtoMuTau.drawOptionSets.default.Ztautau = copy.deepcopy(drawOption_red_separate)
plotTauIdEffZtoMuTau.drawOptionSets.default.Zmumu = copy.deepcopy(drawOption_darkBlue_separate)
plotTauIdEffZtoMuTau.drawOptionSets.default.WplusJets = copy.deepcopy(drawOption_lightBlue_separate)
plotTauIdEffZtoMuTau.drawOptionSets.default.TTplusJets = copy.deepcopy(drawOption_violett_separate)
plotTauIdEffZtoMuTau.drawOptionSets.default.qcdSum = copy.deepcopy(drawOption_orange_separate)
plotTauIdEffZtoMuTau.drawJobs = drawJobConfigurator_TauIdEffZtoMuTau.configure()
plotTauIdEffZtoMuTau.indOutputFileName = cms.string('plotTauIdEffZtoMuTau_#PLOT#.png')
    
plotTauIdEffZtoMuTauShapes = copy.deepcopy(plotZtoMuTau)
plotTauIdEffZtoMuTauShapes.drawOptionEntries = cms.PSet(
    region01 = copy.deepcopy(drawOption_orange_eff),
    region02 = copy.deepcopy(drawOption_red_eff),
    region05 = copy.deepcopy(drawOption_lightBlue_eff),
    region06 = copy.deepcopy(drawOption_green_eff)
)
#plotTauIdEffZtoMuTauShapes.legends.regular.posX = cms.double(0.50)
#plotTauIdEffZtoMuTauShapes.legends.regular.posX = cms.double(0.64)
#plotTauIdEffZtoMuTauShapes.legends.regular.sizeX = cms.double(0.39)
#plotTauIdEffZtoMuTauShapes.legends.regular.sizeY = cms.double(0.25)
plotTauIdEffZtoMuTauShapes.drawJobs = drawJobs_TauIdEffZtoMuTau_shapes
plotTauIdEffZtoMuTauShapes.indOutputFileName = cms.string('plotTauIdEffZtoMuTauShapes_#PLOT#.png')

saveTauIdEffZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('plotsTauIdEffZtoMuTau_all.root')
)


  
