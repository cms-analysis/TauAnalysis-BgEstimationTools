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
from TauAnalysis.BgEstimationTools.tools.prodTemplateHistConfigurator import prodTemplateHistConfigurator
from TauAnalysis.BgEstimationTools.tools.drawTemplateHistConfigurator import drawTemplateHistConfigurator

processName = 'fitTemplateZtoMuTau'
process = cms.Process(processName)

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

#--------------------------------------------------------------------------------
# define event selection of background enriched samples
# from which template histograms are obtained
#--------------------------------------------------------------------------------

bgEstEventSelection_Zmumu = (
    "(numSelDiTaus >= 1 && selMuonTrackIso_0 < 1. && selMuonEcalIso_0 < 1. && selTauDiscrAgainstMuons_0 < 0.5)"
    " || (numSelDiTaus >= 2 && selMuonTrackIso_1 < 1. && selMuonEcalIso_1 < 1. && selTauDiscrAgainstMuons_1 < 0.5)"
)

print("bgEstEventSelection_Zmumu = " + bgEstEventSelection_Zmumu)

bgEstEventSelection_WplusJets = (
    "((numSelDiTaus >= 1 && selMuonTrackIso_0 < 1. && selMuonEcalIso_0 < 1. && selTauDiscrAgainstMuons_0 > 0.5"
    " && selDiTauMt1MET_0 > 40. && selDiTauPzetaDiff_0 < -25.)"
    " || (numSelDiTaus >= 2 && selMuonTrackIso_1 < 1. && selMuonEcalIso_1 < 1. && selTauDiscrAgainstMuons_1 > 0.5"
    " && selDiTauMt1MET_1 > 40. && selDiTauPzetaDiff_1 < -25.))"
    " && numGlobalMuons < 2"
)

print("bgEstEventSelection_WplusJets = " + bgEstEventSelection_WplusJets)

bgEstEventSelection_QCD = (
    "((numSelDiTaus >= 1 && selMuonTrackIso_0 > 4. && selMuonEcalIso_0 > 4. && selTauDiscrAgainstMuons_0 > 0.5)"
    " || (numSelDiTaus >= 2 && selMuonTrackIso_1 > 4. && selMuonEcalIso_1 > 4. && selTauDiscrAgainstMuons_1 > 0.5))"
    " && numGlobalMuons < 2"
)

print("bgEstEventSelection_QCD = " + bgEstEventSelection_QCD)

#--------------------------------------------------------------------------------
# define directories in which histograms are stored in DQMStore
#--------------------------------------------------------------------------------

dqmDirectory_Ztautau_finalEvtSel = 'Ztautau/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/'
dqmDirectory_Ztautau_ZmumuTemplate = 'Ztautau_from_selZmumu/pure/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/'
dqmDirectory_Ztautau_systematics = processName + '/tautau/systematics/'

dqmDirectory_Zmumu_finalEvtSel = 'ZmumuPlusJets/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/'
dqmDirectory_Zmumu_bgEstEnriched_pure = processName + '/Zmumu/pure/'
dqmDirectory_Zmumu_bgEstEnriched_data = processName + '/Zmumu/data/'
dqmDirectory_Zmumu_systematics = processName + '/Zmumu/systematics/'

dqmDirectory_WplusJets_finalEvtSel = 'WplusJets/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/'
dqmDirectory_WplusJets_bgEstEnriched_pure = processName + '/WplusJets/pure/'
dqmDirectory_WplusJets_bgEstEnriched_data = processName + '/WplusJets/data/'
dqmDirectory_WplusJets_systematics = processName + '/WplusJets/systematics/'

dqmDirectory_QCD_finalEvtSel = 'qcdSum/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/'
dqmDirectory_QCD_bgEstEnriched_pure = processName + '/QCD/pure/'
dqmDirectory_QCD_bgEstEnriched_data = processName + '/QCD/data/'
dqmDirectory_QCD_systematics = processName + '/QCD/systematics/'

dqmDirectory_data_finalEvtSel = 'smSum/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/'

#--------------------------------------------------------------------------------
# define names of histograms used in fit
#--------------------------------------------------------------------------------

meName_diTauMvis12 = "VisMass"
meName_diTauMvis12_norm = "VisMassShape"

meName_diTauMt1MET = "Mt1MET"
meName_diTauMt1MET_norm = "Mt1METshape"

#--------------------------------------------------------------------------------
# define names of branches in Ntuple from which template histograms get produced
#--------------------------------------------------------------------------------

branchName_diTauMvis12 = "selDiTauMvis12_0"
branchName_diTauMt1MET = "selDiTauMt1MET_0"

#--------------------------------------------------------------------------------
# produce template histograms 
#--------------------------------------------------------------------------------

prodTemplateHistConfiguratorZtoMuTau = prodTemplateHistConfigurator("prodTemplateHistZtoMuTau", prodTemplateHist, dqmDirectory = processName)
prodTemplateHistConfiguratorZtoMuTau.addProcess("Ztautau", fileNames_Ztautau)
prodTemplateHistConfiguratorZtoMuTau.addProcess("Zmumu", fileNames_ZmumuPlusJets)
prodTemplateHistConfiguratorZtoMuTau.addSelection("Zmumu", bgEstEventSelection_Zmumu)
prodTemplateHistConfiguratorZtoMuTau.addProcess("WplusJets", fileNames_WplusJets)
prodTemplateHistConfiguratorZtoMuTau.addSelection("WplusJets", bgEstEventSelection_WplusJets)
prodTemplateHistConfiguratorZtoMuTau.addProcess("QCD", fileNames_qcdSum)
prodTemplateHistConfiguratorZtoMuTau.addSelection("QCD", bgEstEventSelection_QCD)
prodTemplateHistConfiguratorZtoMuTau.addProcess("data", fileNames_pseudoData)
prodTemplateHistConfiguratorZtoMuTau.addTemplate(meName_diTauMvis12_norm, branchName_diTauMvis12, 40, 0., 200.)
prodTemplateHistConfiguratorZtoMuTau.addTemplate(meName_diTauMt1MET_norm, branchName_diTauMt1MET, 40, 0., 200.)

process.prodTemplateHistZtoMuTau = prodTemplateHistConfiguratorZtoMuTau.configure(process)

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
# normalize to unit area distribution of visible muon + tau-jet mass
# produced from by MCEmbeddingTools from Z --> mu mu events selected in (pseudo)data
#--------------------------------------------------------------------------------

process.normalizeTemplateHistZtoMuTau_Ztautau = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
            meNameInput = cms.string(dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12),
            meNameOutput = cms.string(dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm)
        ),
        cms.PSet(
            meNameInput = cms.string(dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET),
            meNameOutput = cms.string(dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm)
        )
    ),   
    norm = cms.double(1.)
)

#--------------------------------------------------------------------------------
# load distribution of visible muon + tau-jet mass
# expected for different signal/background processes and observed in (pseudo)data
# in events passing final analysis selection criteria
#--------------------------------------------------------------------------------

process.loadAnalysisHistZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    data = cms.PSet(
        inputFileNames = cms.vstring('rfio:/castor/cern.ch/user/v/veelken/bgEstPlots/ZtoMuTau/plotsZtoMuTau_all_fixedCone.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

#--------------------------------------------------------------------------------
# normalize to unit area distribution of visible muon + tau-jet mass
# in simulated signal/background events passing final analysis selection criteria
#--------------------------------------------------------------------------------

process.normalizeAnalysisHistZtoMuTau_Ztautau = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_Ztautau_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12),
             meNameOutput = cms.string(dqmDirectory_Ztautau_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm)
        ),
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_Ztautau_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET),
             meNameOutput = cms.string(dqmDirectory_Ztautau_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm)
        )
    ),
    norm = cms.double(1.)
)

process.normalizeAnalysisHistZtoMuTau_Zmumu = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12),
             meNameOutput = cms.string(dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm)
        ),
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET),
             meNameOutput = cms.string(dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm)
        )
    ),
    norm = cms.double(1.)
)

process.normalizeAnalysisHistZtoMuTau_WplusJets = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12),
             meNameOutput = cms.string(dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm)
        ),
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET),
             meNameOutput = cms.string(dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm)
        )
    ),
    norm = cms.double(1.)
)

process.normalizeAnalysisHistZtoMuTau_QCD = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12),
             meNameOutput = cms.string(dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm)
        ),
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET),
             meNameOutput = cms.string(dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm)
        )
    ),
    norm = cms.double(1.)
)

process.normalizeAnalysisHistZtoMuTau = cms.Sequence(
    process.normalizeAnalysisHistZtoMuTau_Ztautau
   + process.normalizeAnalysisHistZtoMuTau_Zmumu
   + process.normalizeAnalysisHistZtoMuTau_WplusJets
   + process.normalizeAnalysisHistZtoMuTau_QCD
)

#--------------------------------------------------------------------------------
# plot template histograms of "pure" Monte Carlo processes
# compared to the shapes determined by background enriched regions in (pseudo)Data
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorZtoMuTau = drawTemplateHistConfigurator(
    template = drawJobTemplateHist
)

#--------------------------------------------------------------------------------
# define draw jobs for Z --> tau tau signal
#
# NOTE: backgrounds contributing to Z --> mu mu sample from which
#       template histogram for Ztautau signal process is determined using MCEmbeddingTools
#       not included yet, so use "pure" sample as approximation for sample determined in "deta"
#       for the time being...
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorZtoMuTau.add(
    meNames = [
        dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm,
        dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm,
        dqmDirectory_Ztautau_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm
    ],
    name = "Ztautau_diTauMvis12",
    title = "M_{vis}^{#mu + #tau-jet} in Z #rightarrow #tau^{+} #tau^{-} Signal" 
)

drawTemplateHistConfiguratorZtoMuTau.add(
    meNames = [
        dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm,
        dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm,
        dqmDirectory_Ztautau_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm
    ],
    name = "Ztautau_diTauMt1MET",
    title = "M_{T}^{#mu + MET} in Z #rightarrow #tau^{+} #tau^{-} Signal" 
)

#--------------------------------------------------------------------------------
# define draw jobs for Z --> mu mu background
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorZtoMuTau.add(
    meNames = [
        dqmDirectory_Zmumu_bgEstEnriched_data + meName_diTauMvis12_norm,
        dqmDirectory_Zmumu_bgEstEnriched_pure + meName_diTauMvis12_norm,
        dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm
    ],
    name = "Zmumu_diTauMvis12",
    title = "M_{vis}^{#mu + #tau-jet} in Z #rightarrow #mu^{+} #mu^{-} Background"
)

drawTemplateHistConfiguratorZtoMuTau.add(
    meNames = [
        dqmDirectory_Zmumu_bgEstEnriched_data + meName_diTauMt1MET_norm,
        dqmDirectory_Zmumu_bgEstEnriched_pure + meName_diTauMt1MET_norm,
        dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm
    ],
    name = "Zmumu_diTauMt1MET",
    title = "M_{T}^{#mu + MET} in Z #rightarrow #mu^{+} #mu^{-} Background" 
)

#--------------------------------------------------------------------------------
# define draw jobs for W + jets background
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorZtoMuTau.add(
    meNames = [
        dqmDirectory_WplusJets_bgEstEnriched_data + meName_diTauMvis12_norm,
        dqmDirectory_WplusJets_bgEstEnriched_pure + meName_diTauMvis12_norm,
        dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm
    ],
    name = "WplusJets_diTauMvis12",
    title = "M_{vis}^{#mu + #tau-jet} in W + jets Background"
)

drawTemplateHistConfiguratorZtoMuTau.add(
    meNames = [
        dqmDirectory_WplusJets_bgEstEnriched_data + meName_diTauMt1MET_norm,
        dqmDirectory_WplusJets_bgEstEnriched_pure + meName_diTauMt1MET_norm,
        dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm
    ],
    name = "WplusJets_diTauMt1MET",
    title = "M_{T}^{#mu + MET} in W + jets Background" 
)

#--------------------------------------------------------------------------------
# define draw jobs for QCD background
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorZtoMuTau.add(
    meNames = [
        dqmDirectory_QCD_bgEstEnriched_data + meName_diTauMvis12_norm,
        dqmDirectory_QCD_bgEstEnriched_pure + meName_diTauMvis12_norm,
        dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm
    ],
    name = "QCD_diTauMvis12",
    title = "M_{vis}^{#mu + #tau-jet} in QCD Background"
)

drawTemplateHistConfiguratorZtoMuTau.add(
    meNames = [
        dqmDirectory_QCD_bgEstEnriched_data + meName_diTauMt1MET_norm,
        dqmDirectory_QCD_bgEstEnriched_pure + meName_diTauMt1MET_norm,
        dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm
    ],
    name = "QCD_diTauMt1MET",
    title = "M_{T}^{#mu + MET} in QCD Background" 
)

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

    drawJobs = drawTemplateHistConfiguratorZtoMuTau.configure(),

    canvasSizeX = cms.int32(800),
    canvasSizeY = cms.int32(640),                         

    outputFilePath = cms.string('./plots/'),
    #outputFileName = cms.string('plotsTemplateHistZtoMuTau.ps')
    indOutputFileName = cms.string('plotTemplateHistZtoMuTau_#PLOT#.png')
)

#--------------------------------------------------------------------------------
# produce auxiliary histograms representing bias of visible muon + tau-jet mass distribution
# introduced by differences in event selection between final analysis and background enriched samples
#
# NOTE:
#  minuend    = contribution of (pure) signal/background process expected in final analysis
#               (estimated by Monte Carlo)
#  subtrahend = template histogram taken from background enriched sample,
#               including contributions from other signal/background processes
#               (determined by (pseudo)data)
#  difference = minuend - subtrahend
#
# --> in order to account for bias between distribution observed in final analysis
#     and the shapes of signal/background templates fitted to that distribution
#     one needs an **upward** fluctuation of the histogram representing the difference,
#     using a Gaussian of mean 0. and variance 1.
#
#--------------------------------------------------------------------------------

process.prodSysBiasHistZtoMuTau_Ztautau = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_Ztautau_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
            meNameSubtrahend = cms.string(dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
            meNameDifference = cms.string(dqmDirectory_Ztautau_systematics + "bias" + "/" + meName_diTauMvis12_norm)
        ),
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_Ztautau_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm),
            meNameSubtrahend = cms.string(dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm),
            meNameDifference = cms.string(dqmDirectory_Ztautau_systematics + "bias" + "/" + meName_diTauMt1MET_norm)
        )
    )    
)

process.prodSysBiasHistZtoMuTau_Zmumu = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
            meNameSubtrahend = cms.string(dqmDirectory_Zmumu_bgEstEnriched_data + meName_diTauMvis12_norm),
            meNameDifference = cms.string(dqmDirectory_Zmumu_systematics + "bias" + "/" + meName_diTauMvis12_norm)
        ),
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm),
            meNameSubtrahend = cms.string(dqmDirectory_Zmumu_bgEstEnriched_data + meName_diTauMt1MET_norm),
            meNameDifference = cms.string(dqmDirectory_Zmumu_systematics + "bias" + "/" + meName_diTauMt1MET_norm)
        )
    )                                                       
)

process.prodSysBiasHistZtoMuTau_WplusJets = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
            meNameSubtrahend = cms.string(dqmDirectory_WplusJets_bgEstEnriched_data + meName_diTauMvis12_norm),
            meNameDifference = cms.string(dqmDirectory_WplusJets_systematics + "bias" + "/" + meName_diTauMvis12_norm)
        ),
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm),
            meNameSubtrahend = cms.string(dqmDirectory_WplusJets_bgEstEnriched_data + meName_diTauMt1MET_norm),
            meNameDifference = cms.string(dqmDirectory_WplusJets_systematics + "bias" + "/" + meName_diTauMt1MET_norm)
        )
    )                                                       
)

process.prodSysBiasHistZtoMuTau_QCD = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
            meNameSubtrahend = cms.string(dqmDirectory_QCD_bgEstEnriched_data + meName_diTauMvis12_norm),
            meNameDifference = cms.string(dqmDirectory_QCD_systematics + "bias" + "/" + meName_diTauMvis12_norm)
        ),
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm),
            meNameSubtrahend = cms.string(dqmDirectory_QCD_bgEstEnriched_data + meName_diTauMt1MET_norm),
            meNameDifference = cms.string(dqmDirectory_QCD_systematics + "bias" + "/" + meName_diTauMt1MET_norm)
        )
    )                                                       
)

process.prodSysBiasHistZtoMuTau = cms.Sequence(
    process.prodSysBiasHistZtoMuTau_Ztautau
   + process.prodSysBiasHistZtoMuTau_Zmumu
   + process.prodSysBiasHistZtoMuTau_WplusJets
   + process.prodSysBiasHistZtoMuTau_QCD
)

#--------------------------------------------------------------------------------
# store all histograms into ROOT file
#--------------------------------------------------------------------------------

process.saveAllHistZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('templatesZtoMuTau.root')
)

#--------------------------------------------------------------------------------
# fit template histograms to distribution of visible muon + tau-jet mass in (pseudo)Data,
# in order to determine normalization factors of individual background processes
#--------------------------------------------------------------------------------
process.fitZtoMuTau = cms.EDAnalyzer("TemplateBgEstFit",                                          
    processes = cms.PSet(
        Ztautau = cms.PSet(
            meNames = cms.PSet(
                diTauMvis12 = cms.string(dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
                diTauMt1MET = cms.string(dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET_norm)
            ),    
            drawOptions = drawOption_Ztautau
        ),
        Zmumu = cms.PSet(
            meNames = cms.PSet(
                diTauMvis12 = cms.string(dqmDirectory_Zmumu_bgEstEnriched_data + meName_diTauMvis12_norm),
                diTauMt1MET = cms.string(dqmDirectory_Zmumu_bgEstEnriched_data + meName_diTauMt1MET_norm)
            ),    
            drawOptions = drawOption_Zmumu
        ),
        WplusJets = cms.PSet(
            meNames = cms.PSet(
                diTauMvis12 = cms.string(dqmDirectory_WplusJets_bgEstEnriched_data + meName_diTauMvis12_norm),
                diTauMt1MET = cms.string(dqmDirectory_WplusJets_bgEstEnriched_data + meName_diTauMt1MET_norm)
            ),  
            drawOptions = drawOption_WplusJets
        ),
        QCD = cms.PSet(
            meNames = cms.PSet(
                diTauMvis12 = cms.string(dqmDirectory_QCD_bgEstEnriched_data + meName_diTauMvis12_norm),
                diTauMt1MET = cms.string(dqmDirectory_QCD_bgEstEnriched_data + meName_diTauMt1MET_norm)
            ),  
            drawOptions = drawOption_QCD
        )
    ),

    # use "pseudo" data-samples consisting of all Monte Carlo processes for testing                      
    data = cms.PSet(
        meNames = cms.PSet(
            diTauMvis12 = cms.string(dqmDirectory_data_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12),
            diTauMt1MET = cms.string(dqmDirectory_data_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMt1MET)
        )
    ),

    fit = cms.PSet(
        mode = cms.string("Nd"),
        variables = cms.PSet(
            diTauMvis12 = cms.PSet(
               name = cms.string("diTauMvis12"),
               title = cms.string("M_{vis}^{#mu + #tau-jet}"),
               xMin = cms.double(0.),
               xMax = cms.double(150.)
            #),
            #diTauMt1MET = cms.PSet(
            #    name = cms.string("diTauMt1MET"),
            #    title = cms.string("M_{T}^{#mu + MET}"),
            #    xMin = cms.double(40.),
            #    xMax = cms.double(50.)
            )
        ),
        cutUnfittedRegion = cms.bool(False),
        #cutUnfittedRegion = cms.bool(True),
        verbosity = cms.PSet(
            printLevel = cms.int32(1),
            printWarnings = cms.bool(True)
        )
    ),

    estStatUncertainties = cms.PSet(
        numSamplings = cms.PSet(
            stat = cms.int32(1000)
        ),
        chi2redMax = cms.double(3.),
        verbosity = cms.PSet(
            printLevel = cms.int32(-1),
            printWarnings = cms.bool(False)
        )
    ),

    estSysUncertainties = cms.PSet(
        fluctuations = cms.PSet(
            bias = cms.PSet(
                meNames = cms.PSet(
                    Ztautau = cms.PSet(
                        diTauMvis12 = cms.string(dqmDirectory_Ztautau_systematics + "bias" + "/" + meName_diTauMvis12_norm),
                        diTauMt1MET = cms.string(dqmDirectory_Ztautau_systematics + "bias" + "/" + meName_diTauMt1MET_norm)
                    ),
                    Zmumu = cms.PSet(
                        diTauMvis12 = cms.string(dqmDirectory_Zmumu_systematics + "bias" + "/" + meName_diTauMvis12_norm),
                        diTauMt1MET = cms.string(dqmDirectory_Zmumu_systematics + "bias" + "/" + meName_diTauMt1MET_norm)
                    ),
                    WplusJets = cms.PSet(
                        diTauMvis12 = cms.string(dqmDirectory_WplusJets_systematics + "bias" + "/" + meName_diTauMvis12_norm),
                        diTauMt1MET = cms.string(dqmDirectory_WplusJets_systematics + "bias" + "/" + meName_diTauMt1MET_norm)
                    ),
                    QCD = cms.PSet(
                        diTauMvis12 = cms.string(dqmDirectory_QCD_systematics + "bias" + "/" + meName_diTauMvis12_norm),
                        diTauMt1MET = cms.string(dqmDirectory_QCD_systematics + "bias" + "/" + meName_diTauMt1MET_norm)
                    )
                ),
                pullRMS = cms.double(1.),
                pullMin = cms.double(0.),
                pullMax = cms.double(1.),
                mode = cms.string("coherent") # coherent/incoherent
            )
        ),       
        numSamplings = cms.PSet(
            stat = cms.int32(100),
            sys = cms.int32(100)
        ),
        chi2redMax = cms.double(3.),
        verbosity = cms.PSet(
            printLevel = cms.int32(-1),
            printWarnings = cms.bool(False)
        )
    ),                                     

    output = cms.PSet(
        controlPlots = cms.PSet(
            fileName = cms.string("./plots/fitTemplateZtoMuTau_#PLOT#.png")
        )
    )                                      
)                          

process.prodAllHistZtoMuTau = cms.Sequence(
    process.prodTemplateHistZtoMuTau
   + process.loadTemplateHistZtoMuTau_Ztautau
   + process.normalizeTemplateHistZtoMuTau_Ztautau 
   + process.loadAnalysisHistZtoMuTau
   + process.normalizeAnalysisHistZtoMuTau
   + process.plotTemplateHistZtoMuTau
   + process.prodSysBiasHistZtoMuTau
   + process.saveAllHistZtoMuTau
)

process.loadAllHistZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(
        inputFileNames = cms.vstring('templatesZtoMuTau.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

process.p = cms.Path(
    #process.prodAllHistZtoMuTau
    process.loadAllHistZtoMuTau
   + process.fitZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
