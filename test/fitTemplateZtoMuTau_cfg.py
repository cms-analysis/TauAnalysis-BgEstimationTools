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
#
# NOTE:
#  1.) template histogram for Ztautau signal process
#       taken from Z --> mu mu events selected in (pseudo)data,
#       using MCEmbeddingTools
#  2.) distribution observed in (pseudo)data taken from plotsZtoMuTau_all.root file
#      produced by "official" analysis workflow documented at
#       https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideTauAnalysisZtoMuTau
#--------------------------------------------------------------------------------
#
# define event selection of background enriched samples
# from which template histograms are obtained
#
bgEstEventSelection_Zmumu = "muonTrackIso < 1. && muonEcalIso < 1. && tauDiscrAgainstMuons < 0.5"
bgEstEventSelection_WplusJets = "muonTrackIso < 1. && muonEcalIso < 1. && tauDiscrAgainstMuons > 0.5"
bgEstEventSelection_WplusJets += " && diTauMt1MET > 40. && numGlobalMuons < 2"
bgEstEventSelection_QCD = "muonTrackIso > 4. && muonEcalIso > 4. && tauDiscrAgainstMuons > 0.5"
bgEstEventSelection_QCD += " && numGlobalMuons < 2"

print("bgEstEventSelection_Zmumu = " + bgEstEventSelection_Zmumu)
print("bgEstEventSelection_WplusJets = " + bgEstEventSelection_WplusJets)
print("bgEstEventSelection_QCD = " + bgEstEventSelection_QCD)
#
# define observable to be used as template
# and histogram binning options
#
prodTemplateHistZtoMuTau = copy.deepcopy(prodTemplateHist)
prodTemplateHistZtoMuTau.branchName = cms.string('diTauMvis12')
prodTemplateHistZtoMuTau.numBinsX = cms.uint32(30)
prodTemplateHistZtoMuTau.xMin = cms.double(0.)
prodTemplateHistZtoMuTau.xMax = cms.double(150.)
#
# produce template histograms for Z --> mu mu (+ jets) background
#
process.prodTemplateHistZtoMuTau_Zmumu_data = copy.deepcopy(prodTemplateHistZtoMuTau)
process.prodTemplateHistZtoMuTau_Zmumu_data.fileNames = fileNames_pseudoData
process.prodTemplateHistZtoMuTau_Zmumu_data.treeSelection = cms.string(bgEstEventSelection_Zmumu)
process.prodTemplateHistZtoMuTau_Zmumu_data.meName = cms.string("fitTemplateZtoMuTau/Zmumu/data/diTauMvis12")
process.prodTemplateHistZtoMuTau_Zmumu_pure = copy.deepcopy(process.prodTemplateHistZtoMuTau_Zmumu_data)
process.prodTemplateHistZtoMuTau_Zmumu_pure.fileNames = fileNames_ZmumuPlusJets
process.prodTemplateHistZtoMuTau_Zmumu_pure.meName = cms.string("fitTemplateZtoMuTau/Zmumu/pure/diTauMvis12")
process.prodTemplateHistZtoMuTau_Zmumu_Ztautau = copy.deepcopy(process.prodTemplateHistZtoMuTau_Zmumu_data)
process.prodTemplateHistZtoMuTau_Zmumu_Ztautau.fileNames = fileNames_Ztautau
process.prodTemplateHistZtoMuTau_Zmumu_Ztautau.meName = cms.string("fitTemplateZtoMuTau/Zmumu/Ztautau/diTauMvis12")
process.prodTemplateHistZtoMuTau_Zmumu_WplusJets = copy.deepcopy(process.prodTemplateHistZtoMuTau_Zmumu_data)
process.prodTemplateHistZtoMuTau_Zmumu_WplusJets.fileNames = fileNames_WplusJets
process.prodTemplateHistZtoMuTau_Zmumu_WplusJets.meName = cms.string("fitTemplateZtoMuTau/Zmumu/WplusJets/diTauMvis12")
process.prodTemplateHistZtoMuTau_Zmumu_QCD = copy.deepcopy(process.prodTemplateHistZtoMuTau_Zmumu_data)
process.prodTemplateHistZtoMuTau_Zmumu_QCD.fileNames = fileNames_qcdSum
process.prodTemplateHistZtoMuTau_Zmumu_QCD.meName = cms.string("fitTemplateZtoMuTau/Zmumu/QCD/diTauMvis12")
process.prodTemplateHistZtoMuTau_Zmumu = cms.Sequence(
    process.prodTemplateHistZtoMuTau_Zmumu_data
   +process.prodTemplateHistZtoMuTau_Zmumu_pure
   +process.prodTemplateHistZtoMuTau_Zmumu_Ztautau
   +process.prodTemplateHistZtoMuTau_Zmumu_WplusJets
   +process.prodTemplateHistZtoMuTau_Zmumu_QCD
)
#
# produce template histograms for W + jets background
#
process.prodTemplateHistZtoMuTau_WplusJets_data = copy.deepcopy(prodTemplateHistZtoMuTau)
process.prodTemplateHistZtoMuTau_WplusJets_data.fileNames = fileNames_pseudoData
process.prodTemplateHistZtoMuTau_WplusJets_data.treeSelection = cms.string(bgEstEventSelection_WplusJets)
process.prodTemplateHistZtoMuTau_WplusJets_data.meName = cms.string("fitTemplateZtoMuTau/WplusJets/data/diTauMvis12")
process.prodTemplateHistZtoMuTau_WplusJets_pure = copy.deepcopy(process.prodTemplateHistZtoMuTau_WplusJets_data)
process.prodTemplateHistZtoMuTau_WplusJets_pure.fileNames = fileNames_WplusJets
process.prodTemplateHistZtoMuTau_WplusJets_pure.meName = cms.string("fitTemplateZtoMuTau/WplusJets/pure/diTauMvis12")
process.prodTemplateHistZtoMuTau_WplusJets_Ztautau = copy.deepcopy(process.prodTemplateHistZtoMuTau_WplusJets_data)
process.prodTemplateHistZtoMuTau_WplusJets_Ztautau.fileNames = fileNames_Ztautau
process.prodTemplateHistZtoMuTau_WplusJets_Ztautau.meName = cms.string("fitTemplateZtoMuTau/WplusJets/Ztautau/diTauMvis12")
process.prodTemplateHistZtoMuTau_WplusJets_Zmumu = copy.deepcopy(process.prodTemplateHistZtoMuTau_WplusJets_data)
process.prodTemplateHistZtoMuTau_WplusJets_Zmumu.fileNames = fileNames_ZmumuPlusJets
process.prodTemplateHistZtoMuTau_WplusJets_Zmumu.meName = cms.string("fitTemplateZtoMuTau/WplusJets/ZmumuPlusJets/diTauMvis12")
process.prodTemplateHistZtoMuTau_WplusJets_QCD = copy.deepcopy(process.prodTemplateHistZtoMuTau_WplusJets_data)
process.prodTemplateHistZtoMuTau_WplusJets_QCD.fileNames = fileNames_qcdSum
process.prodTemplateHistZtoMuTau_WplusJets_QCD.meName = cms.string("fitTemplateZtoMuTau/WplusJets/QCD/diTauMvis12")
process.prodTemplateHistZtoMuTau_WplusJets = cms.Sequence(
    process.prodTemplateHistZtoMuTau_WplusJets_data
   +process.prodTemplateHistZtoMuTau_WplusJets_pure
   +process.prodTemplateHistZtoMuTau_WplusJets_Ztautau
   +process.prodTemplateHistZtoMuTau_WplusJets_Zmumu
   +process.prodTemplateHistZtoMuTau_WplusJets_QCD
)
#
# produce template histograms for QCD background
#
process.prodTemplateHistZtoMuTau_qcdSum_data = copy.deepcopy(prodTemplateHistZtoMuTau)
process.prodTemplateHistZtoMuTau_qcdSum_data.fileNames = fileNames_pseudoData
process.prodTemplateHistZtoMuTau_qcdSum_data.treeSelection = cms.string(bgEstEventSelection_QCD)
process.prodTemplateHistZtoMuTau_qcdSum_data.meName = cms.string("fitTemplateZtoMuTau/qcdSum/data/diTauMvis12")
process.prodTemplateHistZtoMuTau_qcdSum_pure = copy.deepcopy(process.prodTemplateHistZtoMuTau_qcdSum_data)
process.prodTemplateHistZtoMuTau_qcdSum_pure.fileNames = fileNames_qcdSum
process.prodTemplateHistZtoMuTau_qcdSum_pure.meName = cms.string("fitTemplateZtoMuTau/qcdSum/pure/diTauMvis12")
process.prodTemplateHistZtoMuTau_qcdSum_Ztautau = copy.deepcopy(process.prodTemplateHistZtoMuTau_qcdSum_data)
process.prodTemplateHistZtoMuTau_qcdSum_Ztautau.fileNames = fileNames_Ztautau
process.prodTemplateHistZtoMuTau_qcdSum_Ztautau.meName = cms.string("fitTemplateZtoMuTau/qcdSum/Ztautau/diTauMvis12")
process.prodTemplateHistZtoMuTau_qcdSum_Zmumu = copy.deepcopy(process.prodTemplateHistZtoMuTau_qcdSum_data)
process.prodTemplateHistZtoMuTau_qcdSum_Zmumu.fileNames = fileNames_ZmumuPlusJets
process.prodTemplateHistZtoMuTau_qcdSum_Zmumu.meName = cms.string("fitTemplateZtoMuTau/qcdSum/Zmumu/diTauMvis12")
process.prodTemplateHistZtoMuTau_qcdSum_WplusJets = copy.deepcopy(process.prodTemplateHistZtoMuTau_qcdSum_data)
process.prodTemplateHistZtoMuTau_qcdSum_WplusJets.fileNames = fileNames_WplusJets
process.prodTemplateHistZtoMuTau_qcdSum_WplusJets.meName = cms.string("fitTemplateZtoMuTau/qcdSum/WplusJets/diTauMvis12")
process.prodTemplateHistZtoMuTau_qcdSum = cms.Sequence(
    process.prodTemplateHistZtoMuTau_qcdSum_data
   +process.prodTemplateHistZtoMuTau_qcdSum_pure
   +process.prodTemplateHistZtoMuTau_qcdSum_Ztautau
   +process.prodTemplateHistZtoMuTau_qcdSum_Zmumu 
   +process.prodTemplateHistZtoMuTau_qcdSum_WplusJets
)

process.prodTemplateHistZtoMuTau = cms.Sequence( process.prodTemplateHistZtoMuTau_Zmumu
                                                +process.prodTemplateHistZtoMuTau_WplusJets
                                                +process.prodTemplateHistZtoMuTau_qcdSum )

#--------------------------------------------------------------------------------
# load template histogram of visible muon + tau-jet mass distribution
# produced from by MCEmbeddingTools from Z --> mu mu events selected in (pseudo)data
#--------------------------------------------------------------------------------

process.loadTemplateHistZtoMuTau_Ztautau = cms.EDAnalyzer("DQMFileLoader",
    Ztautau = cms.PSet(
        inputFileNames = cms.vstring('rfio:/castor/cern.ch/user/v/veelken/bgEstTemplates/ZtoMuTau_from_selZmumu.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('Ztautau_from_selZmumu/pure')
    )
)

#--------------------------------------------------------------------------------
# load distribution of visible muon + tau-jet observed in (pseudo)data
#--------------------------------------------------------------------------------

process.loadHistZtoMuTau_data = cms.EDAnalyzer("DQMFileLoader",
    data = cms.PSet(
        inputFileNames = cms.vstring('rfio:/castor/cern.ch/user/v/veelken/bgEstPlots/ZtoMuTau/plotsZtoMuTau_all_fixedCone.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

#--------------------------------------------------------------------------------
# plot template histograms of "pure" Monte Carlo processes
# compared to the shapes determined by background "enriched" regions in (pseudo)Data
#--------------------------------------------------------------------------------
process.drawJob_Ztautau = copy.deepcopy(drawJobTemplateHist)
process.drawJob_Ztautau.plots.dqmMonitorElements = cms.vstring(
    'Ztautau_from_selZmumu/#PROCESSDIR#/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/DiTauCandidateQuantities/VisMass'
)
process.drawJob_Ztautau.plots.processes = cms.vstring('pure')
process.drawJob_Ztautau.title = cms.string("M_{vis}(Muon + Tau)")

process.drawJob_Zmumu = copy.deepcopy(process.drawJob_Ztautau)
process.drawJob_Zmumu.plots.dqmMonitorElements = cms.vstring('fitTemplateZtoMuTau/#PROCESSDIR#/diTauMvis12')
process.drawJob_Zmumu.plots.processes = cms.vstring('data', 'pure')

process.drawJob_WplusJets = copy.deepcopy(process.drawJob_Zmumu)
process.drawJob_WplusJets.plots.processes = cms.vstring('data', 'pure')

process.drawJob_qcdSum = copy.deepcopy(process.drawJob_Zmumu)
process.drawJob_qcdSum.plots.processes = cms.vstring('data', 'pure')

process.plotTemplateHistZtoMuTau = cms.EDAnalyzer("DQMHistPlotter",
    processes = cms.PSet(
        data = cms.PSet(
            dqmDirectory = cms.string('Ztautau_data_from_selZmumu'),
            legendEntry = cms.string('Process enriched in Data'),
            type = cms.string('smMC')
        ),
        pure = cms.PSet(
            dqmDirectory = cms.string('Ztautau_pure_from_selZmumu'),
            legendEntry = cms.string('pure Process'),
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
            data = copy.deepcopy(drawOption_black_separate),
            pure = copy.deepcopy(drawOption_blue_separate),
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
            meName = cms.string(
                'Ztautau_from_selZmumu/pure/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/DiTauCandidateQuantities/VisMass'
            ),
            drawOptions = drawOption_Ztautau
        ),
        Zmumu = cms.PSet(
            meName = process.prodTemplateHistZtoMuTau_Zmumu_data.meName,
            drawOptions = drawOption_Zmumu
        ),
        WplusJets = cms.PSet(
            meName = process.prodTemplateHistZtoMuTau_WplusJets_data.meName,
            drawOptions = drawOption_WplusJets
        ),
        qcdSum = cms.PSet(
            meName = process.prodTemplateHistZtoMuTau_qcdSum_data.meName,
            drawOptions = drawOption_QCD
        )
    ),

    # use "pseudo" data-samples consisting of all Monte Carlo processes for testing                      
    data = cms.PSet(
        meName = cms.string(
            'smSum/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/DiTauCandidateQuantities/VisMass'
        )
    ),

    fit = cms.PSet(
        variableName = prodTemplateHistZtoMuTau.branchName,
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
                     +process.loadTemplateHistZtoMuTau_Ztautau
                     +process.loadHistZtoMuTau_data
                     +process.plotTemplateHistZtoMuTau
                     +process.saveTemplateHistZtoMuTau
                     +process.fitZtoMuTau )

# print-out all python configuration parameter information
#print process.dumpPython()


  
