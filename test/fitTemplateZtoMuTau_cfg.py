import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
#
# Apply "template" method for data-driven background estimation
# to Z --> mu + tau-jet channel
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_cfi import *
from TauAnalysis.BgEstimationTools.bgEstNtupleDefinitionsZtoMuTau_cfi import *
from TauAnalysis.DQMTools.plotterStyleDefinitions_cfi import *
from TauAnalysis.BgEstimationTools.templateHistDefinitions_cfi import *

process = cms.Process('fitTemplateZtoMuTau')

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

#--------------------------------------------------------------------------------
# produce template histograms of visible muon + tau-jet mass distribution
#--------------------------------------------------------------------------------
process.prodTemplateHistZtoMuTau_Ztautau_pure = copy.deepcopy(prodTemplateHist)
process.prodTemplateHistZtoMuTau_Ztautau_pure.fileNames = fileNames_Ztautau
process.prodTemplateHistZtoMuTau_Ztautau_pure.meName = cms.string("fitTemplateZtoMuTau/Ztautau_pure/diTauMvis12")
process.prodTemplateHistZtoMuTau_Ztautau_real = copy.deepcopy(process.prodTemplateHistZtoMuTau_Ztautau_pure)
process.prodTemplateHistZtoMuTau_Ztautau_real.fileNames = fileNames_pseudoData
process.prodTemplateHistZtoMuTau_Ztautau_real.treeSelection = cms.string("")
process.prodTemplateHistZtoMuTau_Ztautau_real.meName = cms.string("fitTemplateZtoMuTau/Ztautau_real/diTauMvis12")
process.prodTemplateHistZtoMuTau_Ztautau = cms.Sequence( process.prodTemplateHistZtoMuTau_Ztautau_pure
                                                        +process.prodTemplateHistZtoMuTau_Ztautau_real )

process.prodTemplateHistZtoMuTau_Zmumu_pure = copy.deepcopy(process.prodTemplateHistZtoMuTau_Ztautau_pure)
process.prodTemplateHistZtoMuTau_Zmumu_pure.fileNames = fileNames_Zmumu
process.prodTemplateHistZtoMuTau_Zmumu_pure.meName = cms.string("fitTemplateZtoMuTau/Zmumu_pure/diTauMvis12")
process.prodTemplateHistZtoMuTau_Zmumu_real = copy.deepcopy(process.prodTemplateHistZtoMuTau_Zmumu_pure)
process.prodTemplateHistZtoMuTau_Zmumu_real.fileNames = fileNames_pseudoData
process.prodTemplateHistZtoMuTau_Zmumu_real.treeSelection = cms.string("")
process.prodTemplateHistZtoMuTau_Zmumu_real.meName = cms.string("fitTemplateZtoMuTau/Zmumu_real/diTauMvis12")
process.prodTemplateHistZtoMuTau_Zmumu = cms.Sequence( process.prodTemplateHistZtoMuTau_Zmumu_pure
                                                      +process.prodTemplateHistZtoMuTau_Zmumu_real )

process.prodTemplateHistZtoMuTau_WplusJets_pure = copy.deepcopy(process.prodTemplateHistZtoMuTau_Ztautau_pure)
process.prodTemplateHistZtoMuTau_WplusJets_pure.fileNames = fileNames_WplusJets
process.prodTemplateHistZtoMuTau_WplusJets_pure.meName = cms.string("fitTemplateZtoMuTau/WplusJets_pure/diTauMvis12")
process.prodTemplateHistZtoMuTau_WplusJets_real = copy.deepcopy(process.prodTemplateHistZtoMuTau_WplusJets_pure)
process.prodTemplateHistZtoMuTau_WplusJets_real.fileNames = fileNames_pseudoData
process.prodTemplateHistZtoMuTau_WplusJets_real.treeSelection = cms.string("")
process.prodTemplateHistZtoMuTau_WplusJets_real.meName = cms.string("fitTemplateZtoMuTau/WplusJets_real/diTauMvis12")
process.prodTemplateHistZtoMuTau_WplusJets = cms.Sequence( process.prodTemplateHistZtoMuTau_WplusJets_pure
                                                          +process.prodTemplateHistZtoMuTau_WplusJets_real )

process.prodTemplateHistZtoMuTau_qcdSum_pure = copy.deepcopy(process.prodTemplateHistZtoMuTau_Ztautau_pure)
process.prodTemplateHistZtoMuTau_qcdSum_pure.fileNames = fileNames_qcdSum
process.prodTemplateHistZtoMuTau_qcdSum_pure.meName = cms.string("fitTemplateZtoMuTau/qcdSum_pure/diTauMvis12")
process.prodTemplateHistZtoMuTau_qcdSum_real = copy.deepcopy(process.prodTemplateHistZtoMuTau_qcdSum_pure)
process.prodTemplateHistZtoMuTau_qcdSum_real.fileNames = fileNames_pseudoData
process.prodTemplateHistZtoMuTau_qcdSum_real.treeSelection = cms.string("")
process.prodTemplateHistZtoMuTau_qcdSum_real.meName = cms.string("fitTemplateZtoMuTau/qcdSum_real/diTauMvis12")
process.prodTemplateHistZtoMuTau_qcdSum = cms.Sequence( process.prodTemplateHistZtoMuTau_qcdSum_pure
                                                       +process.prodTemplateHistZtoMuTau_qcdSum_real )

process.prodTemplateHistZtoMuTau_pseudoData = copy.deepcopy(process.prodTemplateHistZtoMuTau_Ztautau_pure)
process.prodTemplateHistZtoMuTau_pseudoData.fileNames = fileNames_pseudoData
process.prodTemplateHistZtoMuTau_pseudoData.meName = cms.string("fitTemplateZtoMuTau/pseudoData/diTauMvis12")

process.prodTemplateHistZtoMuTau = cms.Sequence( process.prodTemplateHistZtoMuTau_Ztautau
                                                +process.prodTemplateHistZtoMuTau_Zmumu
                                                +process.prodTemplateHistZtoMuTau_WplusJets
                                                +process.prodTemplateHistZtoMuTau_qcdSum
                                                +process.prodTemplateHistZtoMuTau_pseudoData )

#--------------------------------------------------------------------------------
# plot template histograms of "pure" Monte Carlo processes
# compared to the shapes determined by background "enriched" regions in (pseudo)Data
#--------------------------------------------------------------------------------
process.drawJob_Ztautau = copy.deepcopy(process.drawJobTemplateHist)
process.drawJob_Ztautau.plots.dqmMonitorElements = cms.vstring('fitTemplateZtoMuTau/#PROCESSDIR#/diTauMvis12')
process.drawJob_Ztautau.plots.processes = cms.vstring('Ztautau_pure', 'Ztautau_real')
process.drawJob_Ztautau.title = cms.string("M_{vis}(Mu + Tau)")

process.drawJob_Zmumu = copy.deepcopy(process.drawJob_Ztautau)
process.drawJob_Zmumu.plots.processes = cms.vstring('Zmumu_pure', 'Zmumu_real')

process.drawJob_WplusJets = copy.deepcopy(process.drawJob_Ztautau)
process.drawJob_WplusJets.plots.processes = cms.vstring('WplusJets_pure', 'WplusJets_real')

process.drawJob_qcdSum = copy.deepcopy(process.drawJob_Ztautau)
process.drawJob_qcdSum.plots.processes = cms.vstring('qcdSum_pure', 'qcdSum_real')

process.plotTemplateHistZtoMuTau = cms.EDAnalyzer("DQMHistPlotter",
    processes = cms.PSet(
        Ztautau_pure = cms.PSet(
            dqmDirectory = cms.string('Ztautau_pure'),
            legendEntry = cms.string('Z #rightarrow #tau^{+} #tau^{-}'),
            type = cms.string('smMC')
        ),
        Ztautau_real = cms.PSet(
            dqmDirectory = cms.string('Ztautau_real'),
            legendEntry = cms.string('(pseudo)Data, Z #rightarrow #tau^{+} #tau^{-} enriched'),
            type = cms.string('smMC')
        ),
        Zmumu_pure = cms.PSet(
            dqmDirectory = cms.string('Zmumu_pure'),
            legendEntry = cms.string('Z #rightarrow #mu^{+} #mu^{-}'),
            type = cms.string('smMC')
        ),
        Zmumu_real = cms.PSet(
            dqmDirectory = cms.string('Zmumu_real'),
            legendEntry = cms.string('(pseudo)Data, Z #rightarrow #mu^{+} #mu^{-} enriched'),
            type = cms.string('smMC')
        ),
        WplusJets_pure = cms.PSet(
            dqmDirectory = cms.string('WplusJets_pure'),
            legendEntry = cms.string('W + jets'),
            type = cms.string('smMC')
        ),
        WplusJets_real = cms.PSet(
            dqmDirectory = cms.string('WplusJets_real'),
            legendEntry = cms.string('(pseudo)Data, W + jets enriched'),
            type = cms.string('smMC')
        ),
        qcdSum_pure = cms.PSet(
            dqmDirectory = cms.string('qcdSum_pure'),
            legendEntry = cms.string('QCD'),
            type = cms.string('smMC')
        ),
        qcdSum_real = cms.PSet(
            dqmDirectory = cms.string('qcdSum_real'),
            legendEntry = cms.string('(pseudo)Data, QCD enriched'),
            type = cms.string('smMC')
        )
    ),

    xAxes = cms.PSet(
        M = copy.deepcopy(xAxis_mass)
    ),

    yAxes = cms.PSet(                         
        numEntries_linear = copy.deepcopy(yAxis_numEntries_linear),
        numEntries_log = copy.deepcopy(yAxis_numEntries_log)
    ),

    legends = cms.PSet(
        regular = cms.PSet(
            posX = cms.double(0.45),            
            posY = cms.double(0.69),             
            sizeX = cms.double(0.44),        
            sizeY = cms.double(0.20),            
            header = cms.string(''),          
            option = cms.string('brNDC'),       
            borderSize = cms.int32(0),          
            fillColor = cms.int32(0)             
        )
    ),

    labels = cms.PSet(
        mcNormScale = copy.deepcopy(label_mcNormScale)
    ),

    drawOptionSets = cms.PSet(
        default = cms.PSet(
            Ztautau_pure = copy.deepcopy(drawOption_red_separate),
            Ztautau_real = copy.deepcopy(drawOption_black_separate),
            Zmumu_pure = copy.deepcopy(drawOption_red_separate),
            Zmumu_real = copy.deepcopy(drawOption_black_separate),
            WplusJets_pure = copy.deepcopy(drawOption_red_separate),
            WplusJets_real = copy.deepcopy(drawOption_black_separate),
            qcdSum_pure = copy.deepcopy(drawOption_red_separate),
            qcdSum_real = copy.deepcopy(drawOption_black_separate)
        )
    ),

    drawJobs = cms.PSet(
        Ztautau = copy.deepcopy(process.drawJob_Ztautau),
        Zmumu = copy.deepcopy(process.drawJob_Zmumu),
        WplusJets = copy.deepcopy(process.drawJob_WplusJets),
        qcdSum = copy.deepcopy(process.drawJob_qcdSum)
    ),

    canvasSizeX = cms.int32(800),
    canvasSizeY = cms.int32(640),                         

    outputFilePath = cms.string('./plots/'),
    #outputFileName = cms.string('plotsTemplateHistZtoMuTau.ps')
    indOutputFileName = cms.string('plotTemplateHistZtoMuTau_#PLOT#.png')
)

process.saveTemplateHistZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('templatesZtoMuTau.root')
)

#--------------------------------------------------------------------------------
# fit template histograms to distribution of visible muon + tau-jet mass in (pseudo)Data,
# in order to determine normalization factors of individual background processes
#--------------------------------------------------------------------------------
process.bgEstFitZtoMuTau = cms.EDAnalyzer("TemplateBgEstFit",                                          
    processes = cms.PSet(
        Ztautau = cms.PSet(
            meName = process.prodTemplateHistZtoMuTau_Ztautau_pure.meName,
            drawOptions = drawOption_Ztautau
        ),
        Zmumu = cms.PSet(
            meName = process.prodTemplateHistZtoMuTau_Zmumu_pure.meName,
            drawOptions = drawOption_Zmumu
        ),
        WplusJets = cms.PSet(
            meName = process.prodTemplateHistZtoMuTau_WplusJets_pure.meName,
            drawOptions = drawOption_WplusJets
        ),
        qcdSum = cms.PSet(
            meName = process.prodTemplateHistZtoMuTau_qcdSum_pure.meName,
            drawOptions = drawOption_QCD
        )
    ),

    # use "pseudo" data-samples consisting of all Monte Carlo processes for testing                      
    data = cms.PSet(
        meName = process.prodTemplateHistZtoMuTau_pseudoData.meName
    ),

    fit = cms.PSet(
        variableName = process.prodTemplateHistZtoMuTau_Ztautau_pure.branchName,
        variableTitle = cms.string("M_{vis}^{#mu + #tau-jet}"),
        xMin = cms.double(0.),
        xMax = cms.double(150.)
    ),

    output = cms.PSet(
        controlPlots = cms.PSet(
            fileName = cms.string("fitTemplateZtoMuTau.png")
        )
    )                                      
)                          

process.fitZtoMuTau = cms.Sequence( process.bgEstFitZtoMuTau )

process.p = cms.Path( process.prodTemplateHistZtoMuTau
                     +process.plotTemplateHistZtoMuTau
                     +process.saveTemplateHistZtoMuTau
                     +process.fitZtoMuTau )

# print-out all python configuration parameter information
print process.dumpPython()


  
