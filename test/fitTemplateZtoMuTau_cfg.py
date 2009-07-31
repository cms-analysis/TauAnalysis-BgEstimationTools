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
bgEstEventSelection_Zmumu = "(numSelDiTaus >= 1 && selMuonTrackIso_0 < 1. && selMuonEcalIso_0 < 1. && selTauDiscrAgainstMuons_0 < 0.5)"
bgEstEventSelection_Zmumu += " || (numSelDiTaus >= 2 && selMuonTrackIso_1 < 1. && selMuonEcalIso_1 < 1. && selTauDiscrAgainstMuons_1 < 0.5)"
bgEstEventSelection_WplusJets = "((numSelDiTaus >= 1 && selMuonTrackIso_0 < 1. && selMuonEcalIso_0 < 1. && selTauDiscrAgainstMuons_0 > 0.5"
bgEstEventSelection_WplusJets += " && selDiTauMt1MET_0 > 40. && selDiTauPzetaDiff_0 < -25.)"
bgEstEventSelection_WplusJets += " || (numSelDiTaus >= 2 && selMuonTrackIso_1 < 1. && selMuonEcalIso_1 < 1. && selTauDiscrAgainstMuons_1 > 0.5"
bgEstEventSelection_WplusJets += " && selDiTauMt1MET_1 > 40. && selDiTauPzetaDiff_1 < -25.))"
bgEstEventSelection_WplusJets += " && numGlobalMuons < 2"
bgEstEventSelection_QCD = "((numSelDiTaus >= 1 && selMuonTrackIso_0 > 4. && selMuonEcalIso_0 > 4. && selTauDiscrAgainstMuons_0 > 0.5)"
bgEstEventSelection_QCD += " || (numSelDiTaus >= 2 && selMuonTrackIso_1 > 4. && selMuonEcalIso_1 > 4. && selTauDiscrAgainstMuons_1 > 0.5))"
bgEstEventSelection_QCD += " && numGlobalMuons < 2"

print("bgEstEventSelection_Zmumu = " + bgEstEventSelection_Zmumu)
print("bgEstEventSelection_WplusJets = " + bgEstEventSelection_WplusJets)
print("bgEstEventSelection_QCD = " + bgEstEventSelection_QCD)
#
# define observable to be used as template
# and histogram binning options
#
# NOTE: binning needs to match that of the final analysis,
#       defined in TauAnalysis/Core/plugins/..HistManager.cc
# 
#
prodTemplateHistZtoMuTau = copy.deepcopy(prodTemplateHist)
prodTemplateHistZtoMuTau.branchName = cms.string('selDiTauMvis12_0')
prodTemplateHistZtoMuTau.numBinsX = cms.uint32(40)
prodTemplateHistZtoMuTau.xMin = cms.double(0.)
prodTemplateHistZtoMuTau.xMax = cms.double(200.)
prodTemplateHistZtoMuTau.norm = cms.double(1.)
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
process.prodTemplateHistZtoMuTau_QCD_data = copy.deepcopy(prodTemplateHistZtoMuTau)
process.prodTemplateHistZtoMuTau_QCD_data.fileNames = fileNames_pseudoData
process.prodTemplateHistZtoMuTau_QCD_data.treeSelection = cms.string(bgEstEventSelection_QCD)
process.prodTemplateHistZtoMuTau_QCD_data.meName = cms.string("fitTemplateZtoMuTau/QCD/data/diTauMvis12")
process.prodTemplateHistZtoMuTau_QCD_pure = copy.deepcopy(process.prodTemplateHistZtoMuTau_QCD_data)
process.prodTemplateHistZtoMuTau_QCD_pure.fileNames = fileNames_qcdSum
process.prodTemplateHistZtoMuTau_QCD_pure.meName = cms.string("fitTemplateZtoMuTau/QCD/pure/diTauMvis12")
process.prodTemplateHistZtoMuTau_QCD_Ztautau = copy.deepcopy(process.prodTemplateHistZtoMuTau_QCD_data)
process.prodTemplateHistZtoMuTau_QCD_Ztautau.fileNames = fileNames_Ztautau
process.prodTemplateHistZtoMuTau_QCD_Ztautau.meName = cms.string("fitTemplateZtoMuTau/QCD/Ztautau/diTauMvis12")
process.prodTemplateHistZtoMuTau_QCD_Zmumu = copy.deepcopy(process.prodTemplateHistZtoMuTau_QCD_data)
process.prodTemplateHistZtoMuTau_QCD_Zmumu.fileNames = fileNames_ZmumuPlusJets
process.prodTemplateHistZtoMuTau_QCD_Zmumu.meName = cms.string("fitTemplateZtoMuTau/QCD/Zmumu/diTauMvis12")
process.prodTemplateHistZtoMuTau_QCD_WplusJets = copy.deepcopy(process.prodTemplateHistZtoMuTau_QCD_data)
process.prodTemplateHistZtoMuTau_QCD_WplusJets.fileNames = fileNames_WplusJets
process.prodTemplateHistZtoMuTau_QCD_WplusJets.meName = cms.string("fitTemplateZtoMuTau/QCD/WplusJets/diTauMvis12")
process.prodTemplateHistZtoMuTau_QCD = cms.Sequence(
    process.prodTemplateHistZtoMuTau_QCD_data
   +process.prodTemplateHistZtoMuTau_QCD_pure
   +process.prodTemplateHistZtoMuTau_QCD_Ztautau
   +process.prodTemplateHistZtoMuTau_QCD_Zmumu 
   +process.prodTemplateHistZtoMuTau_QCD_WplusJets
)

process.prodTemplateHistZtoMuTau = cms.Sequence( process.prodTemplateHistZtoMuTau_Zmumu
                                                +process.prodTemplateHistZtoMuTau_WplusJets
                                                +process.prodTemplateHistZtoMuTau_QCD )

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
#
# define template histograms for Z --> tau tau signal
#
# NOTE: backgrounds contributing to Z --> mu mu sample from which
#       template histogram for Ztautau signal process is determined using MCEmbeddingTools
#       not included yet, so use "pure" sample as approximation for sample determined in "deta"
#       for the time being...
#
process.drawJob_Ztautau = copy.deepcopy(drawJobTemplateHist)
process.drawJob_Ztautau.plots[0].dqmMonitorElements = cms.vstring(
    'Ztautau_from_selZmumu/pure/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/DiTauCandidateQuantities/VisMass'
)
process.drawJob_Ztautau.plots[1].dqmMonitorElements = cms.vstring(
    'Ztautau_from_selZmumu/pure/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/DiTauCandidateQuantities/VisMass'
)
process.drawJob_Ztautau.plots[2].dqmMonitorElements = cms.vstring(
    'Ztautau/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/DiTauCandidateQuantities/VisMass'
)
process.drawJob_Ztautau.title = cms.string('M_{vis}^{#mu + #tau-jet} in Z #rightarrow #tau^{+} #tau^{-} Signal')
#
# define template histograms for Z --> mu mu background
#
process.drawJob_Zmumu = copy.deepcopy(drawJobTemplateHist)
process.drawJob_Zmumu.plots[0].dqmMonitorElements = cms.vstring(
    'fitTemplateZtoMuTau/Zmumu/data/diTauMvis12'
)
process.drawJob_Zmumu.plots[1].dqmMonitorElements = cms.vstring(
    'fitTemplateZtoMuTau/Zmumu/pure/diTauMvis12'
)
process.drawJob_Zmumu.plots[2].dqmMonitorElements = cms.vstring(
    'ZmumuPlusJets/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/DiTauCandidateQuantities/VisMass'
)
process.drawJob_Zmumu.title = cms.string('M_{vis}^{#mu + #tau-jet} in Z #rightarrow #mu^{+} #mu^{-} Background')
#
# define template histograms for W + jets background
#
process.drawJob_WplusJets = copy.deepcopy(drawJobTemplateHist)
process.drawJob_WplusJets.plots[0].dqmMonitorElements = cms.vstring(
    'fitTemplateZtoMuTau/WplusJets/data/diTauMvis12'
)
process.drawJob_WplusJets.plots[1].dqmMonitorElements = cms.vstring(
    'fitTemplateZtoMuTau/WplusJets/pure/diTauMvis12'
)
process.drawJob_WplusJets.plots[2].dqmMonitorElements = cms.vstring(
    'WplusJets/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/DiTauCandidateQuantities/VisMass'
)
process.drawJob_WplusJets.title = cms.string('M_{vis}^{#mu + #tau-jet} in W + jets Background')
#
# define template histograms for QCD background
#
process.drawJob_QCD = copy.deepcopy(drawJobTemplateHist)
process.drawJob_QCD.plots[0].dqmMonitorElements = cms.vstring(
    'fitTemplateZtoMuTau/QCD/data/diTauMvis12'
)
process.drawJob_QCD.plots[1].dqmMonitorElements = cms.vstring(
    'fitTemplateZtoMuTau/QCD/pure/diTauMvis12'
)
process.drawJob_QCD.plots[2].dqmMonitorElements = cms.vstring(
    'qcdSum/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/DiTauCandidateQuantities/VisMass'
)
process.drawJob_QCD.title = cms.string('M_{vis}^{#mu + #tau-jet} in QCD Background')

process.plotTemplateHistZtoMuTau = cms.EDAnalyzer("DQMHistPlotter",
    processes = cms.PSet(
        bgEstData = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist.plots[0].legendEntry,
            type = cms.string('smMC')
        ),
        bgEstPure = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist.plots[1].legendEntry,
            type = cms.string('smMC')
        ),
        finalEvtSel = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist.plots[2].legendEntry,
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

    drawOptionEntries = cms.PSet(
        bgEstData = copy.deepcopy(drawOption_black_points),
        bgEstPure = copy.deepcopy(drawOption_darkBlue_eff),
        finalEvtSel = copy.deepcopy(drawOption_red_eff)
    ),

    drawJobs = cms.PSet(
        Ztautau = copy.deepcopy(process.drawJob_Ztautau),
        Zmumu = copy.deepcopy(process.drawJob_Zmumu),
        WplusJets = copy.deepcopy(process.drawJob_WplusJets),
        QCD = copy.deepcopy(process.drawJob_QCD)
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
            meName = cms.string(process.drawJob_Ztautau.plots[0].dqmMonitorElements[0]),
            drawOptions = drawOption_Ztautau
        ),
        Zmumu = cms.PSet(
            meName = cms.string(process.drawJob_Zmumu.plots[0].dqmMonitorElements[0]),
            drawOptions = drawOption_Zmumu
        ),
        WplusJets = cms.PSet(
            meName = cms.string(process.drawJob_WplusJets.plots[0].dqmMonitorElements[0]),
            drawOptions = drawOption_WplusJets
        ),
        qcdSum = cms.PSet(
            meName = cms.string(process.drawJob_QCD.plots[0].dqmMonitorElements[0]),
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
            fileName = cms.string("./plots/fitTemplateZtoMuTau.png")
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


  
