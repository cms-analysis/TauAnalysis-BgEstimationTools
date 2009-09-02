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
    "numDiTausZmumu >= 1 && muonTrackIsoZmumu_0 < 1. && muonEcalIsoZmumu_0 < 1. && tauDiscrAgainstMuonsZmumu_0 < 0.5"
    " && diTauAbsChargeZmumu_0 < 0.5"
)

print("bgEstEventSelection_Zmumu = " + bgEstEventSelection_Zmumu)

bgEstEventSelection_WplusJets = (
    "numDiTausWplusJets >= 1 && muonPtWplusJets_0 > 20. && muonTrackIsoWplusJets_0 < 2. && muonEcalIsoWplusJets_0 < 2."
    #"numDiTausWplusJets >= 1 && muonPtWplusJets_0 > 25. && muonTrackIsoWplusJets_0 < 1. && muonEcalIsoWplusJets_0 < 1."
    " && tauTrackIsoDiscrWplusJets_0 < 0.5 && tauTrackIsoWplusJets_0 > 2. && tauDiscrAgainstMuonsWplusJets_0 > 0.5"
    " && diTauMt1MEtWplusJets_0 > 30."
    " && numGlobalMuons < 2"
    " && numJetsAlpha0point1WplusJets < 1"
)

print("bgEstEventSelection_WplusJets = " + bgEstEventSelection_WplusJets)

bgEstEventSelection_TTplusJets = (
    "numDiTausTTplusJets >= 1 && muonTrackIsoTTplusJets_0 < 2. && muonEcalIsoTTplusJets_0 < 2."
    " && tauDiscrAgainstMuonsTTplusJets_0 > 0.5"
    " && diTauAbsChargeTTplusJets_0 < 0.5"
    " && numGlobalMuons < 2"
    " && ((numJetsEt40TTplusJets >= 1 && jetEt40bTaggingDiscrTrackCountingHighEffTTplusJets_0 > 4.5)"
    " || (numJetsEt40TTplusJets >= 2 && jetEt40bTaggingDiscrTrackCountingHighEffTTplusJets_1 > 4.5)"
    " || (numJetsEt40TTplusJets >= 3 && jetEt40bTaggingDiscrTrackCountingHighEffTTplusJets_2 > 4.5))"
    " && numJetsEt40TTplusJets >= 2 && numJetsEt60TTplusJets >= 1"
)

print("bgEstEventSelection_TTplusJets = " + bgEstEventSelection_TTplusJets)

bgEstEventSelection_QCD = (
    "numDiTausQCD >= 1 && muonTrackIsoQCD_0 > 4. && muonEcalIsoQCD_0 > 4."
    " && tauDiscrAgainstMuonsQCD_0 > 0.5"
    " && numGlobalMuons < 2"
)

print("bgEstEventSelection_QCD = " + bgEstEventSelection_QCD)

#--------------------------------------------------------------------------------
# define directories in which histograms are stored in DQMStore
#--------------------------------------------------------------------------------

dqmDirectory_Ztautau_finalEvtSel = 'harvested/Ztautau/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauPzetaDiff/'
dqmDirectory_Ztautau_ZmumuTemplate = 'Ztautau_from_selZmumu/pure/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauMt1MET/'
dqmDirectory_Ztautau_systematics = processName + '/Ztautau/systematics/'

dqmDirectory_Zmumu_finalEvtSel = 'harvested/ZmumuPlusJets/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauPzetaDiff/'
dqmDirectory_Zmumu_bgEstEnriched_pure = processName + '/Zmumu/pure/'
dqmDirectory_Zmumu_bgEstEnriched_data = processName + '/Zmumu/data/'
dqmDirectory_Zmumu_systematics = processName + '/Zmumu/systematics/'

dqmDirectory_WplusJets_finalEvtSel = 'harvested/WplusJets/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauPzetaDiff/'
dqmDirectory_WplusJets_bgEstEnriched_pure = processName + '/WplusJets/pure/'
dqmDirectory_WplusJets_bgEstEnriched_data = processName + '/WplusJets/data/'
dqmDirectory_WplusJets_systematics = processName + '/WplusJets/systematics/'

dqmDirectory_TTplusJets_finalEvtSel = 'harvested/TTplusJets/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauPzetaDiff/'
dqmDirectory_TTplusJets_bgEstEnriched_pure = processName + '/TTplusJets/pure/'
dqmDirectory_TTplusJets_bgEstEnriched_data = processName + '/TTplusJets/data/'
dqmDirectory_TTplusJets_systematics = processName + '/TTplusJets/systematics/'

dqmDirectory_QCD_finalEvtSel = 'harvested/qcdSum/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauPzetaDiff/'
dqmDirectory_QCD_bgEstEnriched_pure = processName + '/QCD/pure/'
dqmDirectory_QCD_bgEstEnriched_data = processName + '/QCD/data/'
dqmDirectory_QCD_systematics = processName + '/QCD/systematics/'

dqmDirectory_data_finalEvtSel = 'harvested/smSum/zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauPzetaDiff/'

#--------------------------------------------------------------------------------
# define names of histograms used in fit
#--------------------------------------------------------------------------------

meName_diTauMvis12 = "VisMass"
meName_diTauMvis12_norm = "VisMassShape"

#--------------------------------------------------------------------------------
# define names of branches in Ntuple from which template histograms get produced
#--------------------------------------------------------------------------------

branchName_diTauMvis12_Zmumu = "diTauMvis12Zmumu_0"
branchName_diTauMvis12_WplusJets = "diTauMvis12WplusJets_0"
branchName_diTauMvis12_TTplusJets = "diTauMvis12TTplusJets_0"
branchName_diTauMvis12_QCD = "diTauMvis12QCD_0"

#--------------------------------------------------------------------------------
# produce template histograms 
#--------------------------------------------------------------------------------

prodTemplateHistConfiguratorZmumuEnriched = prodTemplateHistConfigurator(
    "prodTemplateHistBgEstZmumuEnriched", prodTemplateHist, dqmDirectory = processName
)
prodTemplateHistConfiguratorZmumuEnriched.addProcess("Ztautau", fileNames_Ztautau)
prodTemplateHistConfiguratorZmumuEnriched.addProcess("Zmumu", fileNames_ZmumuPlusJets)
prodTemplateHistConfiguratorZmumuEnriched.addProcess("WplusJets", fileNames_WplusJets)
prodTemplateHistConfiguratorZmumuEnriched.addProcess("TTplusJets", fileNames_TTplusJets)
prodTemplateHistConfiguratorZmumuEnriched.addProcess("QCD", fileNames_qcdSum)
prodTemplateHistConfiguratorZmumuEnriched.addProcess("data", fileNames_pseudoData)
prodTemplateHistConfiguratorZmumuEnriched.addSelection("Zmumu", bgEstEventSelection_Zmumu)
prodTemplateHistConfiguratorZmumuEnriched.addTemplate(meName_diTauMvis12_norm, branchName_diTauMvis12_Zmumu, 40, 0., 200.)

process.prodTemplateHistBgEstZmumuEnriched = prodTemplateHistConfiguratorZmumuEnriched.configure(process)

prodTemplateHistConfiguratorWplusJetsEnriched = prodTemplateHistConfigurator(
    "prodTemplateHistBgEstWplusJetsEnriched", prodTemplateHist, dqmDirectory = processName
)
prodTemplateHistConfiguratorWplusJetsEnriched.addProcess("Ztautau", fileNames_Ztautau)
prodTemplateHistConfiguratorWplusJetsEnriched.addProcess("Zmumu", fileNames_ZmumuPlusJets)
prodTemplateHistConfiguratorWplusJetsEnriched.addProcess("WplusJets", fileNames_WplusJets)
prodTemplateHistConfiguratorWplusJetsEnriched.addProcess("TTplusJets", fileNames_TTplusJets)
prodTemplateHistConfiguratorWplusJetsEnriched.addProcess("QCD", fileNames_qcdSum)
prodTemplateHistConfiguratorWplusJetsEnriched.addProcess("data", fileNames_pseudoData)
prodTemplateHistConfiguratorWplusJetsEnriched.addSelection("WplusJets", bgEstEventSelection_WplusJets)
prodTemplateHistConfiguratorWplusJetsEnriched.addTemplate(meName_diTauMvis12_norm, branchName_diTauMvis12_WplusJets, 40, 0., 200.)

process.prodTemplateHistBgEstWplusJetsEnriched = prodTemplateHistConfiguratorWplusJetsEnriched.configure(process)

prodTemplateHistConfiguratorTTplusJetsEnriched = prodTemplateHistConfigurator(
    "prodTemplateHistBgEstTTplusJetsEnriched", prodTemplateHist, dqmDirectory = processName
)
prodTemplateHistConfiguratorTTplusJetsEnriched.addProcess("Ztautau", fileNames_Ztautau)
prodTemplateHistConfiguratorTTplusJetsEnriched.addProcess("Zmumu", fileNames_ZmumuPlusJets)
prodTemplateHistConfiguratorTTplusJetsEnriched.addProcess("WplusJets", fileNames_WplusJets)
prodTemplateHistConfiguratorTTplusJetsEnriched.addProcess("TTplusJets", fileNames_TTplusJets)
prodTemplateHistConfiguratorTTplusJetsEnriched.addProcess("QCD", fileNames_qcdSum)
prodTemplateHistConfiguratorTTplusJetsEnriched.addProcess("data", fileNames_pseudoData)
prodTemplateHistConfiguratorTTplusJetsEnriched.addSelection("TTplusJets", bgEstEventSelection_TTplusJets)
prodTemplateHistConfiguratorTTplusJetsEnriched.addTemplate(meName_diTauMvis12_norm, branchName_diTauMvis12_TTplusJets, 40, 0., 200.)

process.prodTemplateHistBgEstTTplusJetsEnriched = prodTemplateHistConfiguratorTTplusJetsEnriched.configure(process)

prodTemplateHistConfiguratorQCDenriched = prodTemplateHistConfigurator(
    "prodTemplateHistBgEstWplusJetsEnriched", prodTemplateHist, dqmDirectory = processName
)
prodTemplateHistConfiguratorQCDenriched.addProcess("Ztautau", fileNames_Ztautau)
prodTemplateHistConfiguratorQCDenriched.addProcess("Zmumu", fileNames_ZmumuPlusJets)
prodTemplateHistConfiguratorQCDenriched.addProcess("WplusJets", fileNames_WplusJets)
prodTemplateHistConfiguratorQCDenriched.addProcess("TTplusJets", fileNames_TTplusJets)
prodTemplateHistConfiguratorQCDenriched.addProcess("QCD", fileNames_qcdSum)
prodTemplateHistConfiguratorQCDenriched.addProcess("data", fileNames_pseudoData)
prodTemplateHistConfiguratorQCDenriched.addSelection("QCD", bgEstEventSelection_QCD)
prodTemplateHistConfiguratorQCDenriched.addTemplate(meName_diTauMvis12_norm, branchName_diTauMvis12_QCD, 40, 0., 200.)

process.prodTemplateHistBgEstQCDenriched = prodTemplateHistConfiguratorQCDenriched.configure(process)

process.prodTemplateHistZtoMuTau = cms.Sequence( process.prodTemplateHistBgEstZmumuEnriched
                                                * process.prodTemplateHistBgEstWplusJetsEnriched
                                                * process.prodTemplateHistBgEstTTplusJetsEnriched
                                                * process.prodTemplateHistBgEstQCDenriched )

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
        inputFileNames = cms.vstring('rfio:/castor/cern.ch/user/v/veelken/bgEstPlots/ZtoMuTau/plotsZtoMuTau_all_shrinkingCone.root'),
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
        )
    ),
    norm = cms.double(1.)
)

process.normalizeAnalysisHistZtoMuTau_Zmumu = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12),
             meNameOutput = cms.string(dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm)
        )
    ),
    norm = cms.double(1.)
)

process.normalizeAnalysisHistZtoMuTau_WplusJets = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12),
             meNameOutput = cms.string(dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm)
        )
    ),
    norm = cms.double(1.)
)

process.normalizeAnalysisHistZtoMuTau_TTplusJets = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_TTplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12),
             meNameOutput = cms.string(dqmDirectory_TTplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm)
        )
    ),
    norm = cms.double(1.)
)

process.normalizeAnalysisHistZtoMuTau_QCD = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
             meNameInput = cms.string(dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12),
             meNameOutput = cms.string(dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm)
        )
    ),
    norm = cms.double(1.)
)

process.normalizeAnalysisHistZtoMuTau = cms.Sequence(
    process.normalizeAnalysisHistZtoMuTau_Ztautau
   + process.normalizeAnalysisHistZtoMuTau_Zmumu
   + process.normalizeAnalysisHistZtoMuTau_WplusJets
   + process.normalizeAnalysisHistZtoMuTau_TTplusJets 
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

#--------------------------------------------------------------------------------
# define draw jobs for TTbar + jets background
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorZtoMuTau.add(
    meNames = [
        dqmDirectory_TTplusJets_bgEstEnriched_data + meName_diTauMvis12_norm,
        dqmDirectory_TTplusJets_bgEstEnriched_pure + meName_diTauMvis12_norm,
        dqmDirectory_TTplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm
    ],
    name = "TTplusJets_diTauMvis12",
    title = "M_{vis}^{#mu + #tau-jet} in t#bar{t} + jets Background"
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
        bgEstData = copy.deepcopy(drawOption_darkBlue_eff),
        bgEstPure = copy.deepcopy(drawOption_lightBlue_eff),
        finalEvtSel = copy.deepcopy(drawOption_red_eff)
    ),

    drawJobs = drawTemplateHistConfiguratorZtoMuTau.configure(),

    canvasSizeX = cms.int32(800),
    canvasSizeY = cms.int32(640),                         

    outputFilePath = cms.string('./plots/'),
    #outputFileName = cms.string('plotsTemplateHistZtoMuTau.ps')
    indOutputFileName = cms.string('plotTemplateHistZtoMuTau_#PLOT#.eps')
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
        )
    )    
)

process.prodSysBiasHistZtoMuTau_Zmumu = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_Zmumu_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
            meNameSubtrahend = cms.string(dqmDirectory_Zmumu_bgEstEnriched_data + meName_diTauMvis12_norm),
            meNameDifference = cms.string(dqmDirectory_Zmumu_systematics + "bias" + "/" + meName_diTauMvis12_norm)
        )
    )                                                       
)

process.prodSysBiasHistZtoMuTau_WplusJets = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_WplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
            meNameSubtrahend = cms.string(dqmDirectory_WplusJets_bgEstEnriched_data + meName_diTauMvis12_norm),
            meNameDifference = cms.string(dqmDirectory_WplusJets_systematics + "bias" + "/" + meName_diTauMvis12_norm)
        )
    )                                                       
)

process.prodSysBiasHistZtoMuTau_TTplusJets = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_TTplusJets_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
            meNameSubtrahend = cms.string(dqmDirectory_TTplusJets_bgEstEnriched_data + meName_diTauMvis12_norm),
            meNameDifference = cms.string(dqmDirectory_TTplusJets_systematics + "bias" + "/" + meName_diTauMvis12_norm)
        )
    )                                                       
)

process.prodSysBiasHistZtoMuTau_QCD = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meNameMinuend = cms.string(dqmDirectory_QCD_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
            meNameSubtrahend = cms.string(dqmDirectory_QCD_bgEstEnriched_data + meName_diTauMvis12_norm),
            meNameDifference = cms.string(dqmDirectory_QCD_systematics + "bias" + "/" + meName_diTauMvis12_norm)
        )
    )                                                       
)

process.prodSysBiasHistZtoMuTau = cms.Sequence(
    process.prodSysBiasHistZtoMuTau_Ztautau
   + process.prodSysBiasHistZtoMuTau_Zmumu
   + process.prodSysBiasHistZtoMuTau_WplusJets
   + process.prodSysBiasHistZtoMuTau_TTplusJets
   + process.prodSysBiasHistZtoMuTau_QCD
)

#--------------------------------------------------------------------------------
# store all histograms into ROOT file
#--------------------------------------------------------------------------------

process.saveAllHistZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('bgEstTemplatesZtoMuTau.root')
)

#--------------------------------------------------------------------------------
# fit template histograms to distribution of visible muon + tau-jet mass in (pseudo)Data,
# in order to determine normalization factors of individual background processes
#--------------------------------------------------------------------------------

diTauMvis12_smoothing = cms.PSet(
    pluginName = cms.string("Landau convoluted with Gaussian")
    pluginType = cms.string("TF1landauXgausWrapper"), # defaults to TF1Wrapper
    xMin = cms.double(0.),
    xMax = cms.double(200.),
    parameter = cms.PSet(
        par0 = cms.PSet( # width (scale) parameter of Landau density
            initial = cms.double(5.),
            min = cms.double(1.),
            max = cms.double(50.)
        ),
        par1 = cms.PSet( # most probable (MP, location) parameter of Landau density
            initial = cms.double(50.),
            min = cms.double(20.),
            max = cms.double(150.)
        ),
        par2 = cms.PSet( # total area (integral from -inf to +inf, normalization constant)
            initial = cms.double(1.),
            min = cms.double(0.1),
            max = cms.double(2.)
        ),
        par3 = cms.PSet( # width (sigma) of convoluted Gaussian function
            initial = cms.double(10.),
            min = cms.double(1.),
            max = cms.double(50.)
        )
    )
)

process.fitZtoMuTau = cms.EDAnalyzer("TemplateBgEstFit",                                          
    processes = cms.PSet(
        Ztautau = cms.PSet(
            templates = cms.PSet(
                diTauMvis12 = cms.PSet(
                    meName = cms.string(dqmDirectory_Ztautau_ZmumuTemplate + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12_norm),
                    smoothing = diTauMvis12_smoothing.clone(
                        pluginName = cms.string("diTauMvis12SmoothingZtautau")
                    )
                )
            ),    
            norm = cms.PSet(
                initial = cms.double(1000.)
            ),
            drawOptions = drawOption_Ztautau                
        ),
        Zmumu = cms.PSet(
            templates = cms.PSet(
                diTauMvis12 = cms.PSet(
                    meName = cms.string(dqmDirectory_Zmumu_bgEstEnriched_data + meName_diTauMvis12_norm),
                    smoothing = diTauMvis12_smoothing.clone(
                        pluginName = cms.string("diTauMvis12SmoothingZmumu")
                    )
                )
            ),    
            norm = cms.PSet(
                initial = cms.double(25.)
            ),
            drawOptions = drawOption_Zmumu
        ),
        WplusJets = cms.PSet(
            templates = cms.PSet(
                diTauMvis12 = cms.PSet(
                    meName = cms.string(dqmDirectory_WplusJets_bgEstEnriched_data + meName_diTauMvis12_norm),
                    smoothing = diTauMvis12_smoothing.clone(
                        pluginName = cms.string("diTauMvis12SmoothingWplusJets")
                    )
                )
            ),    
            norm = cms.PSet(
                initial = cms.double(500.)
            ),
            drawOptions = drawOption_WplusJets
        ),
        TTplusJets = cms.PSet(
            templates = cms.PSet(
                diTauMvis12 = cms.PSet(
                    meName = cms.string(dqmDirectory_TTplusJets_bgEstEnriched_data + meName_diTauMvis12_norm),
                    smoothing = diTauMvis12_smoothing.clone(
                        pluginName = cms.string("diTauMvis12SmoothingTTplusJets")
                    )
                )
            ),    
            norm = cms.PSet(
                initial = cms.double(100.)
            ),
            drawOptions = drawOption_TTplusJets
        ),
        QCD = cms.PSet(
            templates = cms.PSet(
                diTauMvis12 = cms.PSet(
                    meName = cms.string(dqmDirectory_QCD_bgEstEnriched_data + meName_diTauMvis12_norm),
                    smoothing = diTauMvis12_smoothing.clone(
                        pluginName = cms.string("diTauMvis12SmoothingQCD")
                    )
                )
            ),    
            norm = cms.PSet(
                initial = cms.double(100.)
            ),
            drawOptions = drawOption_QCD,
        )
    ),

    # use "pseudo" data-samples consisting of all Monte Carlo processes for testing                      
    data = cms.PSet(
        distributions = cms.PSet(
            diTauMvis12 = cms.PSet(
                meName = cms.string(dqmDirectory_data_finalEvtSel + "DiTauCandidateQuantities" + "/" + meName_diTauMvis12)
            )
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
            )
        ),
        # constrain normalization of W + jets, ttbar + jets and QCD backgrounds
        # to Monte Carlo expectation multiplied by "k-factors" determined
        # in background enriched samples
        constraints = cms.PSet(
            Zmumu = cms.PSet(
                norm = cms.double(1.*25.),
                uncertainty = cms.double(25.)
            ),
            WplusJets = cms.PSet(
                norm = cms.double(1.*500.),
                uncertainty = cms.double(250.)
            ),
            TTplusJets = cms.PSet(
                norm = cms.double(1.*100.),
                uncertainty = cms.double(100.)
            ),
            QCD = cms.PSet(
                norm = cms.double(1.*100.),
                uncertainty = cms.double(100.)
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
            stat = cms.int32(10000)
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
                        diTauMvis12 = cms.string(dqmDirectory_Ztautau_systematics + "bias" + "/" + meName_diTauMvis12_norm)
                    ),
                    Zmumu = cms.PSet(
                        diTauMvis12 = cms.string(dqmDirectory_Zmumu_systematics + "bias" + "/" + meName_diTauMvis12_norm)
                    ),
                    WplusJets = cms.PSet(
                        diTauMvis12 = cms.string(dqmDirectory_WplusJets_systematics + "bias" + "/" + meName_diTauMvis12_norm)
                    ),
                    TTplusJets = cms.PSet(
                        diTauMvis12 = cms.string(dqmDirectory_TTplusJets_systematics + "bias" + "/" + meName_diTauMvis12_norm)
                    ),
                    QCD = cms.PSet(
                        diTauMvis12 = cms.string(dqmDirectory_QCD_systematics + "bias" + "/" + meName_diTauMvis12_norm)
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
            fileName = cms.string("./plots/fitTemplateZtoMuTau_#PLOT#.eps")
        )
    )                                      
)                          

process.prodAllHistZtoMuTau = cms.Sequence(
    process.prodTemplateHistZtoMuTau
   + process.loadTemplateHistZtoMuTau_Ztautau
   + process.normalizeTemplateHistZtoMuTau_Ztautau
   + process.loadAnalysisHistZtoMuTau
   + process.normalizeAnalysisHistZtoMuTau
   + process.prodSysBiasHistZtoMuTau
  #+ process.saveAllHistZtoMuTau  
)

process.loadAllHistZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(
        inputFileNames = cms.vstring('bgEstTemplatesZtoMuTau.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

process.p = cms.Path(
    process.prodAllHistZtoMuTau
   #process.loadAllHistZtoMuTau
   + process.plotTemplateHistZtoMuTau
   + process.fitZtoMuTau
)

# print-out all python configuration parameter information
print process.dumpPython()


  
