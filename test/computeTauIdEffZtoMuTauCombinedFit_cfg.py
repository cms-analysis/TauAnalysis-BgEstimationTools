import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
#
# Determine tau id. efficiency by fitting number of Z --> mu + tau-jet candidate events
# observed in regions of tau id. discriminator passed/failed
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

from TauAnalysis.DQMTools.plotterStyleDefinitions_cfi import *
from TauAnalysis.BgEstimationTools.tools.drawTemplateHistConfigurator import drawTemplateHistConfigurator
from TauAnalysis.BgEstimationTools.templateHistDefinitions_cfi import drawJobTemplateHist
from TauAnalysis.Configuration.plotZtoMuTau_cff import plotZtoMuTau
from TauAnalysis.Configuration.plotZtoMuTau_drawJobs_cfi import plots_ZtoMuTau
from TauAnalysis.DQMTools.tools.drawJobConfigurator import *
from TauAnalysis.DQMTools.plotterStyleDefinitions_cfi import *

processName = 'computeTauIdEffZtoMuTauCombinedFit'
process = cms.Process(processName)

process.DQMStore = cms.Service("DQMStore")

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

#--------------------------------------------------------------------------------
# load template histograms and histograms to be fitted
#--------------------------------------------------------------------------------

process.loadTauIdEffZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(        
        inputFileNames = cms.vstring(
            'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/plots/ZtoMuTau_tauIdEff/7TeVii/plotsTauIdEffZtoMuTau_all.root'
        ),        
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

#--------------------------------------------------------------------------------
# define directories in which histograms are stored in DQMStore
#--------------------------------------------------------------------------------

dqmDirectory_Ztautau = 'harvested/Ztautau'
dqmDirectory_Ztautau_normalized = 'harvested/Ztautau_normalized'
dqmDirectory_Ztautau_scaled = 'harvested/Ztautau_scaled'
dqmDirectory_WplusJets = 'harvested/WplusJets'
dqmDirectory_WplusJets_normalized = 'harvested/WplusJets_normalized'
dqmDirectory_WplusJets_scaled = 'harvested/WplusJets_scaled'
dqmDirectory_QCD = 'harvested/qcdSum'
dqmDirectory_QCD_normalized = 'harvested/qcdSum_normalized'
dqmDirectory_QCD_scaled = 'harvested/qcdSum_scaled'
dqmDirectory_Zmumu = 'harvested/Zmumu'
dqmDirectory_TTplusJets = 'harvested/TTplusJets'

dqmDirectory_smSum = 'harvested/smSum'
dqmDirectory_knownBgSum = 'harvested/knownBgSum'
dqmDirectory_knownBgCorr = 'harvested/knownBgCorr'

dqmDirectory_Data = 'harvested/smSum'

processes = cms.vstring(
    'Ztautau',
    'WplusJets',
    'qcdSum',
    'Zmumu',
    'TTplusJets',
    'smSum'
)

dqmDirectory_processDir = 'harvested/#PROCESSDIR#'

tauIdPassed = '2'
tauIdFailed = '1'

dqmSubDirectory_analyzer = 'TauIdEffAnalyzerZtoMuTauCombinedFit'
dqmSubDirectory_afterCuts = 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmSubDirectory_results = dqmSubDirectory_analyzer + '/' + dqmSubDirectory_afterCuts
dqmSubDirectory_binning = dqmSubDirectory_results + '/' + 'tauIdEffBinningResults2regions'
dqmSubDirectory_histograms = dqmSubDirectory_results + '/' + 'tauIdEffHistograms2regions'
dqmSubDirectory_histograms_tauIdPassed = dqmSubDirectory_histograms + '/' + 'region0' + tauIdPassed
dqmSubDirectory_histograms_tauIdFailed = dqmSubDirectory_histograms + '/' + 'region0' + tauIdFailed

sidebandWplusJets_OS_tauIdPassed = '2'
sidebandWplusJets_OS_tauIdFailed = '1'
sidebandWplusJets_SS_tauIdPassed = '4'
sidebandWplusJets_SS_tauIdFailed = '3'

dqmSubDirectory_sidebandWplusJets_analyzer = 'TauIdEffAnalyzerZtoMuTauCombinedFitWplusJets'
dqmSubDirectory_sidebandWplusJets_afterCuts = 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmSubDirectory_sidebandWplusJets_results = dqmSubDirectory_sidebandWplusJets_analyzer + '/' + dqmSubDirectory_sidebandWplusJets_afterCuts
dqmSubDirectory_sidebandWplusJets_binning = dqmSubDirectory_sidebandWplusJets_results + '/' + 'tauIdEffBinningResultsComb2dWplusJets'
dqmSubDirectory_sidebandWplusJets_histograms = dqmSubDirectory_sidebandWplusJets_results + '/' + 'tauIdEffHistogramsComb2dWplusJets'
dqmSubDirectory_sidebandWplusJets_histograms_OS_tauIdPassed = \
  dqmSubDirectory_sidebandWplusJets_histograms + '/' + 'region0' + sidebandWplusJets_OS_tauIdPassed
dqmSubDirectory_sidebandWplusJets_histograms_OS_tauIdFailed = \
  dqmSubDirectory_sidebandWplusJets_histograms + '/' + 'region0' + sidebandWplusJets_OS_tauIdFailed
dqmSubDirectory_sidebandWplusJets_histograms_SS_tauIdPassed = \
  dqmSubDirectory_sidebandWplusJets_histograms + '/' + 'region0' + sidebandWplusJets_SS_tauIdPassed
dqmSubDirectory_sidebandWplusJets_histograms_SS_tauIdFailed = \
  dqmSubDirectory_sidebandWplusJets_histograms + '/' + 'region0' + sidebandWplusJets_SS_tauIdFailed

sidebandQCD_OS_tauIdPassed = '6'
sidebandQCD_OS_tauIdFailed = '5'
sidebandQCD_SS_tauIdPassed = '12'
sidebandQCD_SS_tauIdFailed = '11'

extrapolQCD_SS_tauIdPassed = '8'
extrapolQCD_SS_tauIdFailed = '7'
extrapolQCD_OS_tauIdPassed = '2'
extrapolQCD_OS_tauIdFailed = '1'

dqmSubDirectory_sidebandQCD_analyzer = 'TauIdEffAnalyzerZtoMuTauCombinedFitQCD'
dqmSubDirectory_sidebandQCD_afterCuts = 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmSubDirectory_sidebandQCD_results = dqmSubDirectory_sidebandQCD_analyzer + '/' + dqmSubDirectory_sidebandQCD_afterCuts
dqmSubDirectory_sidebandQCD_binning = dqmSubDirectory_sidebandQCD_results + '/' + 'tauIdEffBinningResultsComb3dQCD'
dqmSubDirectory_sidebandQCD_histograms = dqmSubDirectory_sidebandQCD_results + '/' + 'tauIdEffHistogramsComb3dQCD'
dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed = \
  dqmSubDirectory_sidebandQCD_histograms + '/' + 'region0' + sidebandQCD_OS_tauIdPassed
dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed = \
  dqmSubDirectory_sidebandQCD_histograms + '/' + 'region0' + sidebandQCD_OS_tauIdFailed
dqmSubDirectory_sidebandQCD_histograms_SS_tauIdPassed = \
  dqmSubDirectory_sidebandQCD_histograms + '/' + 'region' + sidebandQCD_SS_tauIdPassed
dqmSubDirectory_sidebandQCD_histograms_SS_tauIdFailed = \
  dqmSubDirectory_sidebandQCD_histograms + '/' + 'region' + sidebandQCD_SS_tauIdFailed

dqmSubDirectory_sidebandQCDfr_analyzer = 'TauIdEffAnalyzerZtoMuTauCombinedFitQCD_fr'
dqmSubDirectory_sidebandQCDfr_afterCuts = 'afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmSubDirectory_sidebandQCDfr_results = dqmSubDirectory_sidebandQCDfr_analyzer + '/' + dqmSubDirectory_sidebandQCDfr_afterCuts
dqmSubDirectory_sidebandQCDfr_binning = dqmSubDirectory_sidebandQCDfr_results + '/' + 'tauIdEffBinningResultsComb3dQCD'
dqmSubDirectory_sidebandQCDfr_histograms = dqmSubDirectory_sidebandQCDfr_results + '/' + 'tauIdEffHistogramsComb3dQCD'
dqmSubDirectory_sidebandQCDfr_histograms_SS_tauIdFailed = \
  dqmSubDirectory_sidebandQCDfr_histograms + '/' + 'region' + sidebandQCD_SS_tauIdFailed

fitRegion_tauIdPassed = extrapolQCD_SS_tauIdPassed
fitRegion_tauIdFailed = extrapolQCD_SS_tauIdFailed

dqmDirectory_WplusJets_template_tauIdPassed = \
  dqmDirectory_WplusJets + '/' + dqmSubDirectory_sidebandQCD_histograms + '/' + 'region0' + fitRegion_tauIdPassed
dqmDirectory_WplusJets_template_tauIdFailed = \
  dqmDirectory_WplusJets + '/' + dqmSubDirectory_sidebandQCD_histograms + '/' + 'region0' + fitRegion_tauIdFailed
#
# CV: take QCD template for tau id. passed region
#     from tau id. failed region, in order to increase event statistics
#
dqmDirectory_QCD_template_pure_tauIdPassed = \
  dqmDirectory_QCD + '/' + dqmSubDirectory_sidebandQCD_histograms + '/' + 'region0' + fitRegion_tauIdPassed
dqmDirectory_QCD_template_pure_tauIdFailed = \
  dqmDirectory_QCD + '/' + dqmSubDirectory_sidebandQCD_histograms + '/' + 'region0' + fitRegion_tauIdFailed
dqmDirectory_QCD_template_data_tauIdPassed = dqmDirectory_QCD + '/' + dqmSubDirectory_sidebandQCD_histograms + '/' + 'region11'
dqmDirectory_QCD_template_data_tauIdFailed = dqmDirectory_QCD + '/' + dqmSubDirectory_sidebandQCD_histograms + '/' + 'region11'

dqmDirectory_Ztautau_template_tauIdPassed = \
  dqmDirectory_Ztautau + '/' + dqmSubDirectory_sidebandQCD_histograms + '/' + 'region0' + fitRegion_tauIdPassed
dqmDirectory_Ztautau_template_tauIdFailed = \
  dqmDirectory_Ztautau + '/' + dqmSubDirectory_sidebandQCD_histograms + '/' + 'region0' + fitRegion_tauIdFailed

#--------------------------------------------------------------------------------
# define names of histograms used in fit
#--------------------------------------------------------------------------------

meName_Mt = 'DiTauCandidateQuantities/Mt1MET'

meOptionsSeparator = "#"
meOptionsNumWeighted = "".join([meOptionsSeparator, "a1", meOptionsSeparator, "s1"])
meOptionsErrWeighted = "".join([meOptionsSeparator, "a2", meOptionsSeparator, "s1"])

#--------------------------------------------------------------------------------
# print-out contributions of different processes
# to tau id. passed/tau id. failed regions and to sidebands
#--------------------------------------------------------------------------------

process.dumpTauIdEffZtoMuTauBinningResults2regions = cms.EDAnalyzer("DQMDumpMonitorElement",
    config = cms.VPSet(
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_binning + '/' + 'binContent_region' + tauIdPassed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_binning + '/' + 'binError_region' + tauIdPassed + meOptionsErrWeighted
           ),
           label = cms.string("tau id. passed region"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_binning + '/' + 'binContent_region' + tauIdFailed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_binning + '/' + 'binError_region' + tauIdFailed + meOptionsErrWeighted
           ),
           label = cms.string("tau id. failed region"),
           processes = processes
        )
    )
)

process.dumpTauIdEffZtoMuTauBinningResultsSidebandQCD = cms.EDAnalyzer("DQMDumpMonitorElement",
    config = cms.VPSet(
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binContent_region' + sidebandQCD_OS_tauIdPassed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binError_region' + sidebandQCD_OS_tauIdPassed + meOptionsErrWeighted
           ),
           label = cms.string("QCD enriched sideband, OS, tau id. passed"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binContent_region' + sidebandQCD_SS_tauIdPassed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binError_region' + sidebandQCD_SS_tauIdPassed + meOptionsErrWeighted
           ),
           label = cms.string("QCD enriched sideband, SS, tau id. passed"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binContent_region' + sidebandQCD_OS_tauIdFailed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binError_region' + sidebandQCD_OS_tauIdFailed + meOptionsErrWeighted
           ),
           label = cms.string("QCD enriched sideband, OS, tau id. failed"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binContent_region' + sidebandQCD_SS_tauIdFailed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binError_region' + sidebandQCD_SS_tauIdFailed + meOptionsErrWeighted
           ),
           label = cms.string("QCD enriched sideband, SS, tau id. failed"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binContent_region' + extrapolQCD_OS_tauIdPassed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binError_region' + extrapolQCD_OS_tauIdPassed + meOptionsErrWeighted
           ),
           label = cms.string("OS 'signal' region, tau id. passed"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binContent_region' + extrapolQCD_SS_tauIdPassed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binError_region' + extrapolQCD_SS_tauIdPassed + meOptionsErrWeighted
           ),
           label = cms.string("SS 'control' region, tau id. passed"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binContent_region' + extrapolQCD_OS_tauIdFailed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binError_region' + extrapolQCD_OS_tauIdFailed + meOptionsErrWeighted
           ),
           label = cms.string("OS control region, tau id. failed"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binContent_region' + extrapolQCD_SS_tauIdFailed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCD_binning + '/' \
            + 'binError_region' + extrapolQCD_SS_tauIdFailed + meOptionsErrWeighted
           ),
           label = cms.string("SS control region, tau id. failed"),
           processes = processes
        )
    )
)

process.dumpTauIdEffZtoMuTauBinningResultsSidebandWplusJets = cms.EDAnalyzer("DQMDumpMonitorElement",
    config = cms.VPSet(
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandWplusJets_binning + '/' \
            + 'binContent_region' + sidebandWplusJets_OS_tauIdPassed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandWplusJets_binning + '/' \
            + 'binError_region' + sidebandWplusJets_OS_tauIdPassed + meOptionsErrWeighted
           ),
           label = cms.string("W + jets enriched sideband, OS, tau id. passed"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandWplusJets_binning + '/' \
            + 'binContent_region' + sidebandWplusJets_SS_tauIdPassed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandWplusJets_binning + '/' \
            + 'binError_region' + sidebandWplusJets_SS_tauIdPassed + meOptionsErrWeighted
           ),
           label = cms.string("W + jets enriched sideband, SS, tau id. passed"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandWplusJets_binning + '/' \
            + 'binContent_region' + sidebandWplusJets_OS_tauIdFailed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandWplusJets_binning + '/' \
            + 'binError_region' + sidebandWplusJets_OS_tauIdFailed + meOptionsErrWeighted
           ),
           label = cms.string("W + jets enriched sideband, OS, tau id. failed"),
           processes = processes
        ),
        cms.PSet(
           meName = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandWplusJets_binning + '/' \
            + 'binContent_region' + sidebandWplusJets_SS_tauIdFailed + meOptionsNumWeighted
           ),
           meName_err = cms.string(
             dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandWplusJets_binning + '/' \
            + 'binError_region' + sidebandWplusJets_SS_tauIdFailed + meOptionsErrWeighted
           ),
           label = cms.string("W + jets enriched sideband, SS, tau id. failed"),
           processes = processes
        )
    )
)

process.dumpTauIdEffZtoMuTauBinningResults = cms.Sequence(
    process.dumpTauIdEffZtoMuTauBinningResults2regions
   * process.dumpTauIdEffZtoMuTauBinningResultsSidebandQCD
   * process.dumpTauIdEffZtoMuTauBinningResultsSidebandWplusJets
)

#--------------------------------------------------------------------------------
# make control plots comparing shape of template histograms
# obtained from background enriched regions in (pseudo)Data
# to "pure" background contributions in 'control' region
#--------------------------------------------------------------------------------

drawJobTemplateHist_Mt = copy.deepcopy(drawJobTemplateHist)
drawJobTemplateHist_Mt.plots[0].process = cms.string('templateData')
drawJobTemplateHist_Mt.plots[0].drawOptionEntry = cms.string('templateData')
drawJobTemplateHist_Mt.plots[0].legendEntry = cms.string("Sideband, Bg. enriched Data")
drawJobTemplateHist_Mt.plots[1].process = cms.string('templatePure')
drawJobTemplateHist_Mt.plots[1].drawOptionEntry = cms.string('templatePure')
drawJobTemplateHist_Mt.plots[1].legendEntry = cms.string("Sideband, pure Bg.")
drawJobTemplateHist_Mt.plots[2].process = cms.string('fittedDistr')
drawJobTemplateHist_Mt.plots[2].drawOptionEntry = cms.string('fitted')
drawJobTemplateHist_Mt.plots[2].legendEntry = cms.string("'control' Region, pure Bg.")
drawJobTemplateHist_Mt.title = cms.string('M_{T}')
drawJobTemplateHist_Mt.xAxis = cms.string('Mass')
drawJobTemplateHist_Mt.yAxis = cms.string('numEntries_log')

drawTemplateHistConfiguratorTauIdEffZtoMuTau_Mt = drawTemplateHistConfigurator(
    template = drawJobTemplateHist_Mt
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_Mt.add(
    meNames = [        
        dqmDirectory_QCD_template_data_tauIdFailed + '/' + meName_Mt,
        dqmDirectory_QCD_template_pure_tauIdFailed + '/' + meName_Mt,
        dqmDirectory_QCD + '/' + dqmSubDirectory_histograms_tauIdPassed + '/' + meName_Mt
    ],
    name = "QCD_Mt",
    title = "M_{T} in QCD Background"
)

process.plotTauIdEffZtoMuTauCombinedFit_Mt = cms.EDAnalyzer("DQMHistPlotter",
    processes = cms.PSet(
        templateData = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist_Mt.plots[0].legendEntry,
            type = cms.string('smMC')
        ),
        templatePure = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist_Mt.plots[1].legendEntry,
            type = cms.string('smMC')
        ),
        fitted = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist_Mt.plots[2].legendEntry,
            type = cms.string('smMC')
        )
    ),

    xAxes = cms.PSet(
        Mass = copy.deepcopy(xAxis_mass)                                      
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
        templateData = copy.deepcopy(drawOption_green_eff),
        templatePure = copy.deepcopy(drawOption_lightBlue_eff),
        fitted = copy.deepcopy(drawOption_red_eff)
    ),

    drawJobs = drawTemplateHistConfiguratorTauIdEffZtoMuTau_Mt.configure(),

    canvasSizeX = cms.int32(800),
    canvasSizeY = cms.int32(640),                         

    outputFilePath = cms.string('./plots/'),
    #outputFileName = cms.string('plotsTauIdEffZtoMuTauCombinedFit.ps')
    indOutputFileName = cms.string('plotTauIdEffZtoMuTauCombinedFit_#PLOT#.png')
)

process.plotTauIdEffZtoMuTauCombinedFit = cms.Sequence(
    process.plotTauIdEffZtoMuTauCombinedFit_Mt
)    

#--------------------------------------------------------------------------------
# compute contributions of W + jets and QCD backgrounds in tau id. passed region
# by scaling number of events observed in tau id. failed region by fake-rate
#--------------------------------------------------------------------------------

process.compFakeRateBgEstimatesTauIdEffZtoMuTauSideband_histogram = cms.EDAnalyzer("DQMHistErrorBandProducer",
    config = cms.VPSet(
        cms.PSet(
            dqmDirectories_inputVariance = cms.vstring(
                dqmDirectory_WplusJets + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_qcdMuEnriched'),
                dqmDirectory_WplusJets + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_qcdDiJetLeadJet'),
                dqmDirectory_WplusJets + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_qcdDiJetSecondLeadJet'),
                dqmDirectory_WplusJets + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_WplusJets')
            ),
            dqmDirectory_output = cms.string(dqmDirectory_WplusJets + '/' + dqmSubDirectory_sidebandQCDfr_histograms),
            method = cms.string("min_max")
        ),
        cms.PSet(
            dqmDirectories_inputVariance = cms.vstring(
                dqmDirectory_QCD + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_qcdMuEnriched'),
                dqmDirectory_QCD + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_qcdDiJetLeadJet'),
                dqmDirectory_QCD + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_qcdDiJetSecondLeadJet'),
                dqmDirectory_QCD + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_WplusJets')
            ),
            dqmDirectory_output = cms.string(dqmDirectory_QCD + '/' + dqmSubDirectory_sidebandQCDfr_histograms),
            method = cms.string("min_max")
        ),
        cms.PSet(
            dqmDirectories_inputVariance = cms.vstring(
                dqmDirectory_Data + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_qcdMuEnriched'),
                dqmDirectory_Data + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_qcdDiJetLeadJet'),
                dqmDirectory_Data + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_qcdDiJetSecondLeadJet'),
                dqmDirectory_Data + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_histograms).replace('QCD_fr', 'QCD_fr_WplusJets')
            ),
            dqmDirectory_output = cms.string(dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCDfr_histograms),
            method = cms.string("min_max")
        )
    )
)               

process.compFakeRateBgEstimatesTauIdEffZtoMuTauSideband_binning = cms.EDAnalyzer("DQMBinningErrorBandProducer",
    config = cms.VPSet(
        cms.PSet(
            dqmDirectories_inputVariance = cms.vstring(
                dqmDirectory_WplusJets + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_qcdMuEnriched'),
                dqmDirectory_WplusJets + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_qcdDiJetLeadJet'),
                dqmDirectory_WplusJets + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_qcdDiJetSecondLeadJet'),
                dqmDirectory_WplusJets + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_WplusJets')
            ),
            dqmDirectory_output = cms.string(dqmDirectory_WplusJets + '/' + dqmSubDirectory_sidebandQCDfr_binning),
            method = cms.string("min_max")
        ),
        cms.PSet(
            dqmDirectories_inputVariance = cms.vstring(
                dqmDirectory_QCD + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_qcdMuEnriched'),
                dqmDirectory_QCD + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_qcdDiJetLeadJet'),
                dqmDirectory_QCD + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_qcdDiJetSecondLeadJet'),
                dqmDirectory_QCD + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_WplusJets')    
            ),
            dqmDirectory_output = cms.string(dqmDirectory_QCD + '/' + dqmSubDirectory_sidebandQCDfr_binning),
            method = cms.string("min_max")
        ),
        cms.PSet(
            dqmDirectories_inputVariance = cms.vstring(
                dqmDirectory_Data + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_qcdMuEnriched'),
                dqmDirectory_Data + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_qcdDiJetLeadJet'),
                dqmDirectory_Data + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_qcdDiJetSecondLeadJet'),
                dqmDirectory_Data + '/' \
               + copy.deepcopy(dqmSubDirectory_sidebandQCDfr_binning).replace('QCD_fr', 'QCD_fr_WplusJets')
            ),
            dqmDirectory_output = cms.string(dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCDfr_binning),
            method = cms.string("min_max")
        )
    )
)

#--------------------------------------------------------------------------------
# print-out background contributions comparing number of events
# in tau id. failed region scaled by fake-rates
# to number of events in tau id. passed region
#--------------------------------------------------------------------------------

process.dumpTauIdEffZtoMuTauBinningResultsSidebandQCDfr = cms.EDAnalyzer("DQMDumpMonitorElement",
    config = cms.VPSet(
        cms.PSet(
            meName = cms.string(
              dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCDfr_binning + '/' \
             + 'binContent_region' + extrapolQCD_SS_tauIdFailed + meOptionsNumWeighted
            ),
            meName_err = cms.string(
              dqmDirectory_processDir + '/' + dqmSubDirectory_sidebandQCDfr_binning + '/' \
             + 'binError_region' + extrapolQCD_SS_tauIdFailed + meOptionsErrWeighted
            ),
            label = cms.string("SS 'control' region, tau id. passed; estimate obtained by fake-rate technique"),
            processes = cms.vstring(
                'WplusJets',
                'qcdSum'
            )
        )
    )
)

#--------------------------------------------------------------------------------
# make control plots comparing distributions
# in tau id. failed region scaled by fake-rates
# to distributions in tau id. passed region
#--------------------------------------------------------------------------------

drawJobTemplateHist_Mt_fr = drawJobTemplateHist_Mt.clone(
    plots = cms.VPSet(
        cms.PSet(
            dqmMonitorElements = cms.vstring(''),
            process = cms.string('tauIdPassed'),
            drawOptionEntry = cms.string('tauIdPassed'),
            legendEntry = cms.string('tau id. passed')
        ),
        cms.PSet(
            dqmMonitorElements = cms.vstring(''),
            process = cms.string('tauIdFailed_scaled'),
            drawOptionEntry = cms.string('tauIdFailed_scaled'),
            legendEntry = cms.string('tau id. failed, scaled by fake-rate')
        )
   )
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_Mt_fr = drawTemplateHistConfigurator(
    template = drawJobTemplateHist_Mt_fr
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_Mt_fr.add(
    meNames = [
        dqmDirectory_WplusJets + '/' + dqmSubDirectory_sidebandQCDfr_histograms_SS_tauIdFailed + '/' + meName_Mt,
        dqmDirectory_WplusJets + '/' + dqmSubDirectory_sidebandQCD_histograms_SS_tauIdPassed + '/' + meName_Mt
    ],
    name = "WplusJets_Mt",
    title = "M_{T} in W + jets Background"
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_Mt_fr.add(
    meNames = [        
        dqmDirectory_QCD + '/' + dqmSubDirectory_sidebandQCDfr_histograms_SS_tauIdFailed + '/' + meName_Mt,
        dqmDirectory_QCD + '/' + dqmSubDirectory_sidebandQCD_histograms_SS_tauIdPassed + '/' + meName_Mt
    ],
    name = "QCD_Mt",
    title = "M_{T} in QCD Background"
)

process.plotTauIdEffZtoMuTauBinningResultsSidebandQCDfr = process.plotTauIdEffZtoMuTauCombinedFit_Mt.clone(
    processes = cms.PSet(
        tauIdPassed = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist_Mt_fr.plots[0].legendEntry,
            type = cms.string('smMC')
        ),
        tauIdFailed_scaled = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist_Mt_fr.plots[1].legendEntry,
            type = cms.string('smMC')
        )
    ),

    drawOptionEntries = cms.PSet(
        tauIdPassed = copy.deepcopy(drawOption_green_eff),
        tauIdFailed_scaled = copy.deepcopy(drawOption_red_eff)
    ),

    drawJobs = drawTemplateHistConfiguratorTauIdEffZtoMuTau_Mt_fr.configure(),

    outputFilePath = cms.string('./plots/'),
    #outputFileName = cms.string('plotsTauIdEffZtoMuTauCombinedFit_fr.ps')
    indOutputFileName = cms.string('plotTauIdEffZtoMuTauCombinedFit_fr_#PLOT#.png')
)

process.compFakeRateBgEstimatesTauIdEffZtoMuTauSideband = cms.Sequence(
    process.compFakeRateBgEstimatesTauIdEffZtoMuTauSideband_histogram
   * process.compFakeRateBgEstimatesTauIdEffZtoMuTauSideband_binning
   * process.dumpTauIdEffZtoMuTauBinningResultsSidebandQCDfr
   * process.plotTauIdEffZtoMuTauBinningResultsSidebandQCDfr
)

#--------------------------------------------------------------------------------
# compute sum of "known" backgrounds
# (backgrounds the contributions of which are assumed
#  to be reliable predicted by Monte Carlo simulation)
#--------------------------------------------------------------------------------

process.addKnownBackgrounds = cms.EDAnalyzer("DQMHistAdder",
    knownBgSum = cms.PSet(
        dqmDirectories_input = cms.vstring(
            dqmDirectory_Zmumu + '/',
            dqmDirectory_TTplusJets + '/'
        ),
        dqmDirectory_output = cms.string(dqmDirectory_knownBgSum)
    )
)

#--------------------------------------------------------------------------------
# subtract sum of "known" backgrounds
# from QCD and W + jets sideband regions
#--------------------------------------------------------------------------------

process.compKnownBgCorrections_tauIdPassed = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meName_minuend = cms.string(
              dqmDirectory_smSum + '/' + dqmSubDirectory_sidebandQCD_histograms \
             + '/' + 'region0' + fitRegion_tauIdPassed + '/' + meName_Mt
            ),
            meName_subtrahend = cms.string(
              dqmDirectory_knownBgSum + '/' + dqmSubDirectory_sidebandQCD_histograms \
             + '/' + 'region0' + fitRegion_tauIdPassed + '/' + meName_Mt
            ),
            meName_difference = cms.string(
              dqmDirectory_knownBgCorr + '/' + dqmSubDirectory_sidebandQCD_histograms \
             + '/' + 'region0' + fitRegion_tauIdPassed + '/' + meName_Mt
            )
        )
    )
)

process.compKnownBgCorrections_tauIdFailed = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meName_minuend = cms.string(
              dqmDirectory_smSum + '/' + dqmSubDirectory_sidebandQCD_histograms \
             + '/' + 'region0' + fitRegion_tauIdFailed + '/' + meName_Mt
            ),
            meName_subtrahend = cms.string(
              dqmDirectory_knownBgSum + '/' + dqmSubDirectory_sidebandQCD_histograms \
             + '/' + 'region0' + fitRegion_tauIdFailed + '/' + meName_Mt
            ),
            meName_difference = cms.string(
              dqmDirectory_knownBgCorr + '/' + dqmSubDirectory_sidebandQCD_histograms \
             + '/' + 'region0' + fitRegion_tauIdFailed + '/' + meName_Mt
            )
        )
    )
)

process.compKnownBgCorrections = cms.Sequence(
    process.compKnownBgCorrections_tauIdPassed
   * process.compKnownBgCorrections_tauIdFailed
)    

#--------------------------------------------------------------------------------
# fit number of events observed in tau id. discriminator passed/failed regions
#--------------------------------------------------------------------------------

process.fitTauIdEffZtoMuTauSideband_tauIdPassed = cms.EDAnalyzer("TemplateHistFitter",                                          
    processes = cms.PSet(
        WplusJets = cms.PSet(
            templates = cms.PSet(
                Mt = cms.PSet(
                    meName = cms.string(dqmDirectory_WplusJets_template_tauIdPassed + '/' + meName_Mt),
                    fitSimultaneously = cms.bool(False)
                )                                                     
            ),    
            norm = cms.PSet(
                initial = cms.double(100.)
            ),
            drawOptions = drawOption_WplusJets_separate
        ),
        QCD = cms.PSet(
            templates = cms.PSet(
                Mt = cms.PSet(
                    meName = cms.string(dqmDirectory_QCD_template_data_tauIdPassed + '/' + meName_Mt),
                    fitSimultaneously = cms.bool(False)
                )
            ),    
            norm = cms.PSet(
                initial = cms.double(100.)
            ),
            drawOptions = drawOption_QCD_separate
        ),
        Ztautau = cms.PSet(
            templates = cms.PSet(
                Mt = cms.PSet(
                    meName = cms.string(dqmDirectory_Ztautau_template_tauIdPassed + '/' + meName_Mt),
                    fitSimultaneously = cms.bool(False)
                )                                                     
            ),    
            norm = cms.PSet(
                initial = cms.double(100.)
            ),
            drawOptions = drawOption_Ztautau_separate
        )
    ),

    # use "pseudo" data-samples consisting of all Monte Carlo processes for testing                      
    data = cms.PSet(
        distributions = cms.PSet(
            Mt = cms.PSet(
                meName = process.compKnownBgCorrections_tauIdPassed.config[0].meName_difference
            )                                                     
        )
    ),

    fit = cms.PSet(
        algorithm = cms.PSet(
            pluginName = cms.string("fitTauIdEffZtoMuTauSideband_tauIdPassed"),
            #pluginType = cms.string("TemplateFitAdapter_TFractionFitter")
            pluginType = cms.string("TemplateFitAdapter_RooFit")
        ),
        variables = cms.PSet(
            Mt = cms.PSet(
                name = cms.string("Mt"),
                title = cms.string("M_{T}"),
                min = cms.double(0.),
                max = cms.double(200.)
            )
        ),
        verbosity = cms.PSet(
            printLevel = cms.int32(1),
            printWarnings = cms.bool(True)
        )
    ),
    estStatUncertainties = cms.PSet(
        numSamplings = cms.int32(0),
        chi2redMax = cms.double(3.),
        verbosity = cms.PSet(
            printLevel = cms.int32(-1),
            printWarnings = cms.bool(False)
        )
    ),

    estSysUncertainties = cms.PSet(
        fluctuations = cms.PSet(),
        numSamplings = cms.int32(0),
        chi2redMax = cms.double(3.),
        verbosity = cms.PSet(
            printLevel = cms.int32(-1),
            printWarnings = cms.bool(False)
        )
    ),                                                                

    output = cms.PSet(
        controlPlots = cms.PSet(
            fileName = cms.string("./plots/fitTauIdEffZtoMuTauSideband_tauIdPassed_#PLOT#.png")
        ),
        fitResults = cms.PSet(
            dqmDirectory = cms.string("fitTauIdEffZtoMuTauSideband_tauIdPassed/fitResults")
        )
    )                                      
)                          

process.fitTauIdEffZtoMuTauSideband_tauIdFailed = copy.deepcopy(process.fitTauIdEffZtoMuTauSideband_tauIdPassed)
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.WplusJets.templates.Mt.meName = \
  dqmDirectory_WplusJets_template_tauIdFailed + '/' + meName_Mt
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.WplusJets.norm.initial = cms.double(2500.)
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.QCD.templates.Mt.meName = \
  dqmDirectory_QCD_template_data_tauIdFailed + '/' + meName_Mt
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.QCD.norm.initial = cms.double(1000.)
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.Ztautau.templates.Mt.meName = \
  dqmDirectory_Ztautau_template_tauIdFailed + '/' + meName_Mt
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.Ztautau.norm.initial = cms.double(100.)
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.data.distributions.Mt.meName = \
  process.compKnownBgCorrections_tauIdFailed.config[0].meName_difference                                         
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.output.controlPlots.fileName = \
  "./plots/fitTauIdEffZtoMuTauSideband_tauIdFailed_#PLOT#.png"
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.output.fitResults.dqmDirectory = \
  "fitTauIdEffZtoMuTauSideband_tauIdFailed/fitResults"

process.fitTauIdEffZtoMuTauSideband = cms.Sequence(
    process.addKnownBackgrounds
   * process.compKnownBgCorrections
   * process.fitTauIdEffZtoMuTauSideband_tauIdPassed
   * process.fitTauIdEffZtoMuTauSideband_tauIdFailed
)

#--------------------------------------------------------------------------------
# extrapolate contributions of QCD and W + jet backgrounds
# from sideband into signal region
#--------------------------------------------------------------------------------

meName_WplusJets_norm = 'WplusJets/norm/value' + meOptionsNumWeighted
meName_WplusJets_normErr = 'WplusJets/norm/error' + meOptionsErrWeighted

meName_QCD_norm = 'QCD/norm/value' + meOptionsNumWeighted
meName_QCD_normErr = 'QCD/norm/error' + meOptionsErrWeighted

meName_Ztautau_norm = 'Ztautau/norm/value' + meOptionsNumWeighted
meName_Ztautau_normErr = 'Ztautau/norm/error' + meOptionsErrWeighted

process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed = cms.EDAnalyzer("DQMHistScaler",
    config = cms.VPSet(
        cms.PSet(
            meName_input = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdPassed/fitResults' + '/' + meName_WplusJets_norm
            ),
            meName_inputErr = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdPassed/fitResults' + '/' + meName_WplusJets_normErr
            ),
            meName_numerator = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandWplusJets_binning \
             + '/' + 'binContent_region' + sidebandWplusJets_OS_tauIdPassed + meOptionsNumWeighted
            ),
            meName_numeratorErr = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandWplusJets_binning \
             + '/' + 'binError_region' + sidebandWplusJets_OS_tauIdPassed + meOptionsErrWeighted
            ),
            meName_denominator = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandWplusJets_binning \
             + '/' + 'binContent_region' + sidebandWplusJets_SS_tauIdPassed + meOptionsNumWeighted
            ),
            meName_denominatorErr = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandWplusJets_binning \
             + '/' + 'binError_region' + sidebandWplusJets_SS_tauIdPassed + meOptionsErrWeighted
            ),
            meType = cms.string("real"),
            meName_output = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/norm_WplusJets' + meOptionsNumWeighted
            ),
            meName_outputErr = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/normErr_WplusJets' + meOptionsErrWeighted
            )
        ),
        cms.PSet(
            meName_input = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdPassed/fitResults' + '/' + meName_QCD_norm
            ),
            meName_inputErr = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdPassed/fitResults' + '/' + meName_QCD_normErr
            ),
            meName_numerator = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binContent_region' + sidebandQCD_OS_tauIdPassed + meOptionsNumWeighted
            ),
            meName_numeratorErr = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binError_region' + sidebandQCD_OS_tauIdPassed + meOptionsErrWeighted
            ),
            meName_denominator = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binContent_region' + sidebandQCD_SS_tauIdPassed + meOptionsNumWeighted
            ),
            meName_denominatorErr = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binError_region' + sidebandQCD_SS_tauIdPassed + meOptionsErrWeighted
            ),
            meType = cms.string("real"),
            meName_output = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/norm_QCD' + meOptionsNumWeighted
            ),
            meName_outputErr = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/normErr_QCD' + meOptionsErrWeighted
            )
        )
    )                                                                             
)

process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed = cms.EDAnalyzer("DQMHistScaler",
    config = cms.VPSet(
        cms.PSet(
            meName_input = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdFailed/fitResults' + '/' + meName_WplusJets_norm
            ),
            meName_inputErr = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdFailed/fitResults' + '/' + meName_WplusJets_normErr
            ),
            meName_numerator = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandWplusJets_binning \
             + '/' + 'binContent_region' + sidebandWplusJets_OS_tauIdFailed + meOptionsNumWeighted
            ),
            meName_numeratorErr = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandWplusJets_binning \
             + '/' + 'binError_region' + sidebandWplusJets_OS_tauIdFailed + meOptionsErrWeighted
            ),
            meName_denominator = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandWplusJets_binning \
             + '/' + 'binContent_region' + sidebandWplusJets_SS_tauIdFailed + meOptionsNumWeighted
            ),
            meName_denominatorErr = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandWplusJets_binning \
             + '/' + 'binError_region' + sidebandWplusJets_SS_tauIdFailed + meOptionsErrWeighted
            ),
            meType = cms.string("real"),
            meName_output = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/norm_WplusJets' + meOptionsNumWeighted
            ),
            meName_outputErr = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/normErr_WplusJets' + meOptionsErrWeighted
            )
        ),
        cms.PSet(
            meName_input = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdFailed/fitResults' + '/' + meName_QCD_norm
            ),
            meName_inputErr = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdFailed/fitResults' + '/' + meName_QCD_normErr
            ),
            meName_numerator = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binContent_region' + sidebandQCD_OS_tauIdFailed + meOptionsNumWeighted
            ),
            meName_numeratorErr = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binError_region' + sidebandQCD_OS_tauIdFailed + meOptionsErrWeighted
            ),
            meName_denominator = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binContent_region' + sidebandQCD_SS_tauIdFailed + meOptionsNumWeighted
            ),
            meName_denominatorErr = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binError_region' + sidebandQCD_SS_tauIdFailed + meOptionsErrWeighted
            ),
            meType = cms.string("real"),
            meName_output = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/norm_QCD' + meOptionsNumWeighted
            ),
            meName_outputErr = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/normErr_QCD' + meOptionsErrWeighted
            )
        )
    )                                                                             
)

process.extrapolateTauIdEffZtoMuTauSideband = cms.Sequence(
    process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed
   * process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed
)    

#--------------------------------------------------------------------------------
# print-out contributions of W + jets and QCD backgrounds
# to 'signal' region in which tau id. efficiency is determined
#--------------------------------------------------------------------------------

process.dumpTauIdEffZtoMuTauBinningResultsSidebandQCDextrapol = cms.EDAnalyzer("DQMDumpMonitorElement",
    config = cms.VPSet(
        cms.PSet(
           meName = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[0].meName_output,
           meName_err = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[0].meName_outputErr,
           label = cms.string("W + jets background estimated contribution to OS 'signal' region, tau id. passed")
        ),
        cms.PSet(
           meName = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[1].meName_output,
           meName_err = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[1].meName_outputErr,
           label = cms.string("QCD background estimated contribution to OS 'signal' region, tau id. passed")
        ),
        cms.PSet(
           meName = process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[0].meName_output,
           meName_err = process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[0].meName_outputErr,
           label = cms.string("W + jets background estimated contribution to OS 'signal' region, tau id. failed")
        ),
        cms.PSet(
           meName = process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[1].meName_output,
           meName_err = process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[1].meName_outputErr,
           label = cms.string("QCD background estimated contribution to OS 'signal' region, tau id. failed")
        )
    )
)

#--------------------------------------------------------------------------------
# compute sum of "known" backgrounds plus
# QCD and W + jets background contributions extrapolated from sidebands
#--------------------------------------------------------------------------------

process.compBgSumSideband_tauIdPassed = cms.EDAnalyzer("DQMHistAdder",
    config = cms.VPSet(                                                               
        cms.PSet(
            input = cms.VPSet( 
                cms.PSet( 
                    meName = cms.string(
                      dqmDirectory_Zmumu + '/' + dqmSubDirectory_sidebandQCD_binning \
                     + '/' + 'binContent_region' + extrapolQCD_OS_tauIdPassed + meOptionsNumWeighted
                    ),
                    meName_err = cms.string(
                      dqmDirectory_Zmumu + '/' + dqmSubDirectory_sidebandQCD_binning \
                     + '/' + 'binError_region' + extrapolQCD_OS_tauIdPassed + meOptionsErrWeighted
                    )
                ),
                cms.PSet(
                    meName = cms.string(
                      dqmDirectory_TTplusJets + '/' + dqmSubDirectory_sidebandQCD_binning \
                     + '/' + 'binContent_region' + extrapolQCD_OS_tauIdPassed + meOptionsNumWeighted
                    ),
                    meName_err = cms.string(
                      dqmDirectory_TTplusJets + '/' + dqmSubDirectory_sidebandQCD_binning \
                     + '/' + 'binError_region' + extrapolQCD_OS_tauIdPassed + meOptionsErrWeighted
                    )
                ),
                cms.PSet(
                    meName = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[0].meName_output,
                    meName_err = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[0].meName_outputErr
                ),
                cms.PSet(
                    meName = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[1].meName_output,
                    meName_err = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[1].meName_outputErr
                )
            ),
            output = cms.PSet(
                meName = cms.string(
                  'extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/norm_smBgSum' + meOptionsNumWeighted
                ),
                meName_err = cms.string(
                  'extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/normErr_smBgSum' + meOptionsErrWeighted
                )
            )
        )
    )
)

process.compBgSumSideband_tauIdFailed = cms.EDAnalyzer("DQMHistAdder",
    config = cms.VPSet(                                                               
        cms.PSet(
            input = cms.VPSet(
                cms.PSet( 
                    meName = cms.string(
                      dqmDirectory_Zmumu + '/' + dqmSubDirectory_sidebandQCD_binning \
                     + '/' + 'binContent_region' + extrapolQCD_OS_tauIdFailed + meOptionsNumWeighted
                    ),
                    meName_err = cms.string(
                      dqmDirectory_Zmumu + '/' + dqmSubDirectory_sidebandQCD_binning \
                     + '/' + 'binError_region' + extrapolQCD_OS_tauIdFailed + meOptionsErrWeighted
                    )
                ),
                cms.PSet(
                    meName = cms.string(
                      dqmDirectory_TTplusJets + '/' + dqmSubDirectory_sidebandQCD_binning \
                     + '/' + 'binContent_region' + extrapolQCD_OS_tauIdFailed + meOptionsNumWeighted
                    ),
                    meName_err = cms.string(
                      dqmDirectory_TTplusJets + '/' + dqmSubDirectory_sidebandQCD_binning \
                     + '/' + 'binError_region' + extrapolQCD_OS_tauIdFailed + meOptionsErrWeighted
                    )
                ),
                cms.PSet(
                    meName = process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[0].meName_output,
                    meName_err = process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[0].meName_outputErr
                ),
                cms.PSet(
                    meName = process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[1].meName_output,
                    meName_err = process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[1].meName_outputErr
                )
            ),
            output = cms.PSet(
                meName = cms.string(
                  'extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/norm_smBgSum' + meOptionsNumWeighted
                ),
                meName_err = cms.string(
                  'extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/normErr_smBgSum' + meOptionsErrWeighted
                )
            )
        )
   )
)

process.compBgSumSideband = cms.Sequence(
    process.compBgSumSideband_tauIdPassed
   * process.compBgSumSideband_tauIdFailed
)

#--------------------------------------------------------------------------------
# subtract sum of "known" backgrounds plus
# QCD and W + jets background contributions extrapolated from sidebands
# from number of events observed in data
# (to obtain an estimate for the number of Ztautau events passing/failing tau id. criteria)
#--------------------------------------------------------------------------------

process.compBgCorrectionsSideband_tauIdPassed = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meName_minuend = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binContent_region' + extrapolQCD_OS_tauIdPassed + meOptionsNumWeighted
            ),
            meName_minuendErr = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binError_region' + extrapolQCD_OS_tauIdPassed + meOptionsErrWeighted
            ),
            meName_subtrahend = process.compBgSumSideband_tauIdPassed.config[0].output.meName,
            meName_subtrahendErr = process.compBgSumSideband_tauIdPassed.config[0].output.meName_err,
            meName_difference = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/norm_Ztautau' + meOptionsNumWeighted
            ),
            meName_differenceErr = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/normErr_Ztautau' + meOptionsErrWeighted
            )
        )
    )
)

process.compBgCorrectionsSideband_tauIdFailed = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meName_minuend = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binContent_region' + extrapolQCD_OS_tauIdFailed + meOptionsNumWeighted
            ),
            meName_minuendErr = cms.string(
              dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_binning \
             + '/' + 'binError_region' + extrapolQCD_OS_tauIdFailed + meOptionsErrWeighted
            ),
            meName_subtrahend = process.compBgSumSideband_tauIdFailed.config[0].output.meName,
            meName_subtrahendErr = process.compBgSumSideband_tauIdFailed.config[0].output.meName_err,
            meName_difference = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/norm_Ztautau' + meOptionsNumWeighted
            ),
            meName_differenceErr = cms.string(
              'extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/normErr_Ztautau' + meOptionsErrWeighted
            )
        )
    )
)

process.compBgCorrectionsSideband = cms.Sequence(
    process.compBgCorrectionsSideband_tauIdPassed
   * process.compBgCorrectionsSideband_tauIdFailed
)    

#--------------------------------------------------------------------------------
# compute tau id. efficiencies
# (using binomial errors for events statistics expected in data)
#--------------------------------------------------------------------------------

process.dumpBinErrorsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMBinErrorCalculator",
    config = cms.VPSet(
        cms.PSet(
            meName_passed = cms.string(
              dqmDirectory_Ztautau + '/' + dqmSubDirectory_binning + '/' + 'binContent_region' + tauIdPassed + meOptionsNumWeighted
            ),
            meName_passedErr = cms.string(
              dqmDirectory_Ztautau + '/' + dqmSubDirectory_binning + '/' + 'binError_region' + tauIdPassed + meOptionsErrWeighted
            ),
            meName_failed = cms.string(
              dqmDirectory_Ztautau + '/' + dqmSubDirectory_binning + '/' + 'binContent_region' + tauIdFailed + meOptionsNumWeighted
            ),
            meName_failedErr = cms.string(
              dqmDirectory_Ztautau + '/' + dqmSubDirectory_binning + '/' + 'binError_region' + tauIdFailed + meOptionsErrWeighted
            ),
            label = cms.string("Tau id., true")
        ),
        cms.PSet(
            meName_passed = process.compBgCorrectionsSideband_tauIdPassed.config[0].meName_difference,
            meName_passedErr = process.compBgCorrectionsSideband_tauIdPassed.config[0].meName_differenceErr,
            meName_failed = process.compBgCorrectionsSideband_tauIdFailed.config[0].meName_difference,
            meName_failedErr = process.compBgCorrectionsSideband_tauIdFailed.config[0].meName_differenceErr,
            label = cms.string("Tau id., fitted")
        )
    )
)

#--------------------------------------------------------------------------------
# make control plots of muon Pt and eta, tau-jet Pt and eta,
# transverse mass of muon + MEt and of visible invariant mass of muon + tau-jet
#--------------------------------------------------------------------------------

process.compNormalizedDistributions_tauIdPassed = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
            dqmDirectory_input = cms.string(
              dqmDirectory_Ztautau + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
            ),
            dqmDirectory_output = cms.string(
              dqmDirectory_Ztautau_normalized + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
            )
        ),
        cms.PSet(
            dqmDirectory_input = cms.string(
              dqmDirectory_WplusJets + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
            ),
            dqmDirectory_output = cms.string(
              dqmDirectory_WplusJets_normalized + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
            )
        ),
        cms.PSet(
            dqmDirectory_input = cms.string(
              dqmDirectory_QCD + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
            ),
            dqmDirectory_output = cms.string(
              dqmDirectory_QCD_normalized + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
            )
        )
    ),
    norm = cms.double(1.)                                                             
)

process.compScaledDistributions_tauIdPassed = cms.EDAnalyzer("DQMHistScaler",
    config = cms.VPSet(
        cms.PSet(
            dqmDirectory_input = process.compNormalizedDistributions_tauIdPassed.config[0].dqmDirectory_output,
            dqmDirectory_output = cms.string(
              dqmDirectory_Ztautau_scaled + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
            ),
            meName_scaleFactor = process.compBgCorrectionsSideband_tauIdPassed.config[0].meName_difference,
            meName_scaleFactorErr = process.compBgCorrectionsSideband_tauIdPassed.config[0].meName_differenceErr,
            meType = cms.string("real")
        ),
        cms.PSet(
            dqmDirectory_input = process.compNormalizedDistributions_tauIdPassed.config[1].dqmDirectory_output,
            dqmDirectory_output = cms.string(
              dqmDirectory_WplusJets_scaled + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
            ),
            meName_scaleFactor = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[0].meName_output,
            meName_scaleFactorErr = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[0].meName_outputErr,
            meType = cms.string("real")
        ),
        cms.PSet(
            dqmDirectory_input = process.compNormalizedDistributions_tauIdPassed.config[2].dqmDirectory_output,
            dqmDirectory_output = cms.string(
              dqmDirectory_QCD_scaled + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
            ),
            meName_scaleFactor = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[1].meName_output,
            meName_scaleFactorErr = process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed.config[1].meName_outputErr,
            meType = cms.string("real")
        )
    )
)

process.plotScaledDistributions_template = copy.deepcopy(plots_ZtoMuTau)
process.plotScaledDistributions_template.plots.processes = cms.vstring(
    'Zmumu',
    'WplusJets',
    'TTplusJets',
    'qcdSum',
    'Ztautau',
    'Data'
)

process.plotScaledDistributions_tauIdPassed = copy.deepcopy(plotZtoMuTau)
process.plotScaledDistributions_tauIdPassed.processes.Ztautau.dqmDirectory = \
  process.compScaledDistributions_tauIdPassed.config[0].dqmDirectory_output
process.plotScaledDistributions_tauIdPassed.processes.WplusJets.dqmDirectory = \
  process.compScaledDistributions_tauIdPassed.config[1].dqmDirectory_output
process.plotScaledDistributions_tauIdPassed.processes.qcdSum.dqmDirectory = \
  process.compScaledDistributions_tauIdPassed.config[2].dqmDirectory_output
process.plotScaledDistributions_tauIdPassed.processes.Zmumu.dqmDirectory = \
  dqmDirectory_Zmumu + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
process.plotScaledDistributions_tauIdPassed.processes.TTplusJets.dqmDirectory = \
  dqmDirectory_TTplusJets + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed
process.plotScaledDistributions_tauIdPassed.processes.Data = cms.PSet(
   dqmDirectory = cms.string(dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdPassed),
   legendEntry = cms.string('Data'),
   type = cms.string('Data') # 'Data' / 'smMC' / 'bsmMC' / 'smSumMC'
)
process.plotScaledDistributions_tauIdPassed.drawOptionSets.default.Data = drawOption_Data
drawJobConfigurator_TauIdEffZtoMuTau_tauIdPassed = drawJobConfigurator(
    template = process.plotScaledDistributions_template,
    dqmDirectory = '#PROCESSDIR#/'
)
drawJobConfigurator_TauIdEffZtoMuTau_tauIdPassed.add(
    plots = [
        drawJobConfigEntry(
            meName = 'MuonQuantities/Muon#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Muon",
            xAxis = '#PAR#',
            name = "muon"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau",
            xAxis = '#PAR#',
            name = "tau"
        ),
        drawJobConfigEntry(
            meName = 'TauIdEffSpecificQuantities/DiTauVisMassFromJetP4',
            title = "M_{vis}(Muon + Jet)",
            xAxis = 'Mass',
            name = "mVisiblePFJet"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/Mt1MET',
            title = "M_{T}(Muon + MET)",
            xAxis = 'Mt',
            name = "mtMuonMET"
        ),
    ]
)
process.plotScaledDistributions_tauIdPassed.drawJobs = drawJobConfigurator_TauIdEffZtoMuTau_tauIdPassed.configure()
process.plotScaledDistributions_tauIdPassed.indOutputFileName = cms.string('computeTauIdEffZtoMuTauCombinedFit_tauIdPassed_#PLOT#.png')

process.makeControlPlotsTauIdEffZtoMuTauSideband_tauIdPassed = cms.Sequence(
    process.compNormalizedDistributions_tauIdPassed
   * process.compScaledDistributions_tauIdPassed
   * process.plotScaledDistributions_tauIdPassed
)

process.compNormalizedDistributions_tauIdFailed = copy.deepcopy(process.compNormalizedDistributions_tauIdPassed)
process.compNormalizedDistributions_tauIdFailed.config[0].dqmDirectory_input = \
  dqmDirectory_Ztautau + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed
process.compNormalizedDistributions_tauIdFailed.config[0].dqmDirectory_output = \
  dqmDirectory_Ztautau_normalized + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed
process.compNormalizedDistributions_tauIdFailed.config[1].dqmDirectory_input = \
  dqmDirectory_WplusJets + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed
process.compNormalizedDistributions_tauIdFailed.config[1].dqmDirectory_output = \
  dqmDirectory_WplusJets_normalized + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed
process.compNormalizedDistributions_tauIdFailed.config[2].dqmDirectory_input = \
  dqmDirectory_QCD + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed
process.compNormalizedDistributions_tauIdFailed.config[2].dqmDirectory_output = \
  dqmDirectory_QCD_normalized + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed

process.compScaledDistributions_tauIdFailed = copy.deepcopy(process.compScaledDistributions_tauIdPassed)
process.compScaledDistributions_tauIdFailed.config[0].dqmDirectory_input = \
  process.compNormalizedDistributions_tauIdFailed.config[0].dqmDirectory_output
process.compScaledDistributions_tauIdFailed.config[0].dqmDirectory_output = \
  dqmDirectory_Ztautau_scaled + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed
process.compScaledDistributions_tauIdFailed.config[0].meName_scaleFactor = \
  process.compBgCorrectionsSideband_tauIdFailed.config[0].meName_difference
process.compScaledDistributions_tauIdFailed.config[0].meName_scaleFactorErr = \
  process.compBgCorrectionsSideband_tauIdFailed.config[0].meName_differenceErr
process.compScaledDistributions_tauIdFailed.config[1].dqmDirectory_input = \
  process.compNormalizedDistributions_tauIdFailed.config[1].dqmDirectory_output
process.compScaledDistributions_tauIdFailed.config[1].dqmDirectory_output = \
  dqmDirectory_WplusJets_scaled + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed
process.compScaledDistributions_tauIdFailed.config[1].meName_scaleFactor = \
  process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[0].meName_output
process.compScaledDistributions_tauIdFailed.config[1].meName_scaleFactorErr = \
  process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[0].meName_outputErr
process.compScaledDistributions_tauIdFailed.config[2].dqmDirectory_input = \
  process.compNormalizedDistributions_tauIdFailed.config[2].dqmDirectory_output
process.compScaledDistributions_tauIdFailed.config[2].dqmDirectory_output = \
  dqmDirectory_QCD_scaled + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed
process.compScaledDistributions_tauIdFailed.config[2].meName_scaleFactor = \
  process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[1].meName_output
process.compScaledDistributions_tauIdFailed.config[2].meName_scaleFactorErr = \
  process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed.config[1].meName_outputErr

process.plotScaledDistributions_tauIdFailed = copy.deepcopy(process.plotScaledDistributions_tauIdPassed)
process.plotScaledDistributions_tauIdFailed.processes.Ztautau.dqmDirectory = \
  process.compScaledDistributions_tauIdFailed.config[0].dqmDirectory_output
process.plotScaledDistributions_tauIdFailed.processes.WplusJets.dqmDirectory = \
  process.compScaledDistributions_tauIdFailed.config[1].dqmDirectory_output
process.plotScaledDistributions_tauIdFailed.processes.qcdSum.dqmDirectory = \
  process.compScaledDistributions_tauIdFailed.config[2].dqmDirectory_output
process.plotScaledDistributions_tauIdFailed.processes.Zmumu.dqmDirectory = \
  dqmDirectory_Zmumu + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed  
process.plotScaledDistributions_tauIdFailed.processes.TTplusJets.dqmDirectory = \
  dqmDirectory_TTplusJets + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed
process.plotScaledDistributions_tauIdFailed.processes.Data.dqmDirectory = \
  dqmDirectory_Data + '/' + dqmSubDirectory_sidebandQCD_histograms_OS_tauIdFailed
process.plotScaledDistributions_tauIdFailed.indOutputFileName = cms.string('computeTauIdEffZtoMuTauCombinedFit_tauIdFailed_#PLOT#.png')

process.makeControlPlotsTauIdEffZtoMuTauSideband_tauIdFailed = cms.Sequence(
    process.compNormalizedDistributions_tauIdFailed
   * process.compScaledDistributions_tauIdFailed
   * process.plotScaledDistributions_tauIdFailed
)

process.makeControlPlotsTauIdEffZtoMuTauSideband = cms.Sequence(
    process.makeControlPlotsTauIdEffZtoMuTauSideband_tauIdPassed
   * process.makeControlPlotsTauIdEffZtoMuTauSideband_tauIdFailed
)    

process.saveFitResultsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('fitTauIdEffZtoMuTau_results.root'),
    outputCommands = cms.vstring('drop harvested/*')
)

process.p = cms.Path(
    process.loadTauIdEffZtoMuTau
   #+ process.dumpDQMStore
   + process.dumpTauIdEffZtoMuTauBinningResults
   + process.plotTauIdEffZtoMuTauCombinedFit
   + process.compFakeRateBgEstimatesTauIdEffZtoMuTauSideband
   + process.fitTauIdEffZtoMuTauSideband
   + process.extrapolateTauIdEffZtoMuTauSideband
   + process.dumpTauIdEffZtoMuTauBinningResultsSidebandQCDextrapol 
   + process.compBgSumSideband
   + process.compBgCorrectionsSideband
   + process.dumpBinErrorsTauIdEffZtoMuTau
   + process.makeControlPlotsTauIdEffZtoMuTauSideband
   + process.saveFitResultsTauIdEffZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
