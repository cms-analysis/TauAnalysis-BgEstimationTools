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
from TauAnalysis.BgEstimationTools.plotTauIdEffZtoMuTau_cff import plotTauIdEffZtoMuTau, plots_TauIdEffZtoMuTau
from TauAnalysis.DQMTools.tools.drawJobConfigurator import *

processName = 'computeTauIdEffZtoMuTauCombinedFit'
process = cms.Process(processName)

process.DQMStore = cms.Service("DQMStore")

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
            'rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_6_x/plots/ZtoMuTau_tauIdEff/7TeV/plotsTauIdEffZtoMuTau_all.root'
        ),        
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

#--------------------------------------------------------------------------------
# define directories in which histograms are stored in DQMStore
#--------------------------------------------------------------------------------

dqmSubDirectory_control_tauIdPassed = 'tauIdEffHistograms2regions/region02'
dqmSubDirectory_control_tauIdFailed = 'tauIdEffHistograms2regions/region01'

dqmDirectory_WplusJets_sideband = \
  'harvested/WplusJets/TauIdEffAnalyzerZtoMuTauCombinedFitWplusJets/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_WplusJets_sidebandOS_tauIdPassed = dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffHistogramsComb3dWplusJets/region02'
dqmDirectory_WplusJets_sidebandOS_tauIdFailed = dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffHistogramsComb3dWplusJets/region01'
dqmDirectory_WplusJets_sidebandSS_tauIdPassed = dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffHistogramsComb3dWplusJets/region06'
dqmDirectory_WplusJets_sidebandSS_tauIdFailed = dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffHistogramsComb3dWplusJets/region05'
dqmDirectory_WplusJets_template = \
  'harvested/WplusJets/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_WplusJets_template_tauIdPassed = dqmDirectory_WplusJets_template + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
dqmDirectory_WplusJets_template_tauIdFailed = dqmDirectory_WplusJets_template + '/' + 'tauIdEffHistogramsComb3dQCD/region07'
dqmDirectory_WplusJets_control = \
  'harvested/WplusJets/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_WplusJets_control_tauIdPassed = dqmDirectory_WplusJets_control + '/' + 'tauIdEffHistograms2regions/region02'
dqmDirectory_WplusJets_control_tauIdFailed = dqmDirectory_WplusJets_control + '/' + 'tauIdEffHistograms2regions/region01'
dqmDirectory_WplusJets_controlNormalized = \
  'harvested/WplusJets_normalized/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_WplusJets_controlNormalized_tauIdPassed = \
  dqmDirectory_WplusJets_controlNormalized + '/' + dqmSubDirectory_control_tauIdPassed
dqmDirectory_WplusJets_controlNormalized_tauIdFailed = \
  dqmDirectory_WplusJets_controlNormalized + '/' + dqmSubDirectory_control_tauIdFailed
dqmDirectory_WplusJets_controlScaled = \
  'harvested/WplusJets_scaled/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_WplusJets_controlScaled_tauIdPassed = \
  dqmDirectory_WplusJets_controlScaled + '/' + dqmSubDirectory_control_tauIdPassed
dqmDirectory_WplusJets_controlScaled_tauIdFailed = \
  dqmDirectory_WplusJets_controlScaled + '/' + dqmSubDirectory_control_tauIdFailed

dqmDirectory_QCD_sideband = \
  'harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_QCD_sidebandOS_tauIdPassed = dqmDirectory_QCD_sideband + '/' + 'tauIdEffHistogramsComb3dQCD/region06'
dqmDirectory_QCD_sidebandOS_tauIdFailed = dqmDirectory_QCD_sideband + '/' + 'tauIdEffHistogramsComb3dQCD/region05'
dqmDirectory_QCD_sidebandSS_tauIdPassed = dqmDirectory_QCD_sideband + '/' + 'tauIdEffHistogramsComb3dQCD/region12'
dqmDirectory_QCD_sidebandSS_tauIdFailed = dqmDirectory_QCD_sideband + '/' + 'tauIdEffHistogramsComb3dQCD/region11'
dqmDirectory_QCD_template_pure = \
  'harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_QCD_template_pure_tauIdPassed = dqmDirectory_QCD_template_pure + '/' + 'tauIdEffHistogramsComb3dQCD/region07'
dqmDirectory_QCD_template_pure_tauIdFailed = dqmDirectory_QCD_template_pure + '/' + 'tauIdEffHistogramsComb3dQCD/region07'
dqmDirectory_QCD_template_data = \
  'harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_QCD_template_data_tauIdPassed = dqmDirectory_QCD_template_data + '/' + 'tauIdEffHistogramsComb3dQCD/region11'
dqmDirectory_QCD_template_data_tauIdFailed = dqmDirectory_QCD_template_data + '/' + 'tauIdEffHistogramsComb3dQCD/region11'
dqmDirectory_QCD_control = \
  'harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_QCD_control_tauIdPassed = dqmDirectory_QCD_control + '/' + 'tauIdEffHistograms2regions/region02'
dqmDirectory_QCD_control_tauIdFailed = dqmDirectory_QCD_control + '/' + 'tauIdEffHistograms2regions/region01'
dqmDirectory_QCD_controlNormalized = \
  'harvested/qcdSum_normalized/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_QCD_controlNormalized_tauIdPassed = \
  dqmDirectory_QCD_controlNormalized + '/' + dqmSubDirectory_control_tauIdPassed
dqmDirectory_QCD_controlNormalized_tauIdFailed = \
  dqmDirectory_QCD_controlNormalized + '/' + dqmSubDirectory_control_tauIdFailed
dqmDirectory_QCD_controlScaled = \
  'harvested/qcdSum_scaled/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_QCD_controlScaled_tauIdPassed = \
  dqmDirectory_QCD_controlScaled + '/' + dqmSubDirectory_control_tauIdPassed
dqmDirectory_QCD_controlScaled_tauIdFailed = \
  dqmDirectory_QCD_controlScaled + '/' + dqmSubDirectory_control_tauIdFailed

dqmDirectory_Ztautau_template = \
  'harvested/Ztautau/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_Ztautau_template_tauIdPassed = dqmDirectory_Ztautau_template + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
dqmDirectory_Ztautau_template_tauIdFailed = dqmDirectory_Ztautau_template + '/' + 'tauIdEffHistogramsComb3dQCD/region07'
dqmDirectory_Ztautau_control = \
  'harvested/Ztautau/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_Ztautau_control_tauIdPassed = dqmDirectory_Ztautau_control + '/' + 'tauIdEffHistograms2regions/region02'
dqmDirectory_Ztautau_control_tauIdFailed = dqmDirectory_Ztautau_control + '/' + 'tauIdEffHistograms2regions/region01'
dqmDirectory_Ztautau_controlNormalized = \
  'harvested/Ztautau_normalized/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_Ztautau_controlNormalized_tauIdPassed = \
  dqmDirectory_Ztautau_controlNormalized + '/' + dqmSubDirectory_control_tauIdPassed
dqmDirectory_Ztautau_controlNormalized_tauIdFailed = \
  dqmDirectory_Ztautau_controlNormalized + '/' + dqmSubDirectory_control_tauIdFailed
dqmDirectory_Ztautau_controlScaled = \
  'harvested/Ztautau_scaled/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_Ztautau_controlScaled_tauIdPassed = \
  dqmDirectory_Ztautau_controlScaled + '/' + dqmSubDirectory_control_tauIdPassed
dqmDirectory_Ztautau_controlScaled_tauIdFailed = \
  dqmDirectory_Ztautau_controlScaled + '/' + dqmSubDirectory_control_tauIdFailed

dqmDirectory_Zmumu_control = \
  'harvested/Zmumu/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_Zmumu_control_tauIdPassed = dqmDirectory_Zmumu_control + '/' + 'tauIdEffHistograms2regions/region02'
dqmDirectory_Zmumu_control_tauIdFailed = dqmDirectory_Zmumu_control + '/' + 'tauIdEffHistograms2regions/region01'

dqmDirectory_TTplusJets_control = \
  'harvested/TTplusJets/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_TTplusJets_control_tauIdPassed = dqmDirectory_TTplusJets_control + '/' + 'tauIdEffHistograms2regions/region02'
dqmDirectory_TTplusJets_control_tauIdFailed = dqmDirectory_TTplusJets_control + '/' + 'tauIdEffHistograms2regions/region01'

dqmDirectory_smSum_control = \
  'harvested/smSum_extrapolated/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_smSum_control_tauIdPassed = dqmDirectory_smSum_control + '/' + 'tauIdEffHistograms2regions/region02'
dqmDirectory_smSum_control_tauIdFailed = dqmDirectory_smSum_control + '/' + 'tauIdEffHistograms2regions/region01'

dqmDirectory_Data = \
  'harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'

#--------------------------------------------------------------------------------
# define names of histograms used in fit
#--------------------------------------------------------------------------------

meName_Mt = 'DiTauCandidateQuantities/Mt1MET'

#--------------------------------------------------------------------------------
# plot template histograms of "pure" Monte Carlo processes
# compared to the shapes determined by background enriched regions in (pseudo)Data
#--------------------------------------------------------------------------------

drawJobTemplateHist_Mt = copy.deepcopy(drawJobTemplateHist)
drawJobTemplateHist_Mt.plots[0].process = cms.string('templateData')
drawJobTemplateHist_Mt.plots[0].drawOptionEntry = cms.string('templateData')
drawJobTemplateHist_Mt.plots[0].legendEntry = cms.string('Sideband, Bg. enriched Data')
drawJobTemplateHist_Mt.plots[1].process = cms.string('templatePure')
drawJobTemplateHist_Mt.plots[1].drawOptionEntry = cms.string('templatePure')
drawJobTemplateHist_Mt.plots[1].legendEntry = cms.string('Sideband, pure Bg.')
drawJobTemplateHist_Mt.plots[2].process = cms.string('fittedDistr')
drawJobTemplateHist_Mt.plots[2].drawOptionEntry = cms.string('fitted')
drawJobTemplateHist_Mt.plots[2].legendEntry = cms.string('fitted, pure Bg.')
drawJobTemplateHist_Mt.title = cms.string('M_{T}')
drawJobTemplateHist_Mt.xAxis = cms.string('Mass')
drawJobTemplateHist_Mt.yAxis = cms.string('numEntries_log')

drawTemplateHistConfiguratorTauIdEffZtoMuTau_Mt = drawTemplateHistConfigurator(
    template = drawJobTemplateHist_Mt
)

#--------------------------------------------------------------------------------
# define draw jobs for QCD background
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorTauIdEffZtoMuTau_Mt.add(
    meNames = [        
        dqmDirectory_QCD_template_data_tauIdFailed + '/' + meName_Mt,
        dqmDirectory_QCD_template_pure_tauIdFailed + '/' + meName_Mt,
        'harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFit/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit' \
          + '/' + 'tauIdEffHistograms2regions/region02' + '/' + meName_Mt
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
# compute sum of "known" backgrounds
# (backgrounds the contributions of which are assumed
#  to be reliable predicted by Monte Carlo simulation)
#--------------------------------------------------------------------------------

process.addKnownBackgrounds = cms.EDAnalyzer("DQMHistAdder",
    knownBgSum = cms.PSet(
        dqmDirectories_input = cms.vstring(
            'harvested/Zmumu/',
            'harvested/TTplusJets/'
        ),
        dqmDirectory_output = cms.string('harvested/knownBgSum/')
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
              'harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region08' + '/' + meName_Mt
            ),
            meName_subtrahend = cms.string(
              'harvested/knownBgSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region08' + '/' + meName_Mt
            ),
            meName_difference = cms.string(
              'harvested/knownBgCorr/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region08' + '/' + meName_Mt
            )
        )
    )
)

process.compKnownBgCorrections_tauIdFailed = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(    
        cms.PSet(
            meName_minuend = cms.string(
              'harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region07' + '/' + meName_Mt
            ),
            meName_subtrahend = cms.string(
              'harvested/knownBgSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region07' + '/' + meName_Mt
            ),
            meName_difference = cms.string(
              'harvested/knownBgCorr/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region07' + '/' + meName_Mt
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
                    meName = cms.string(dqmDirectory_QCD_template_pure_tauIdPassed + '/' + meName_Mt),
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
                meName = cms.string(
                  'harvested/knownBgCorr/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit' \
                 + '/' + 'tauIdEffHistogramsComb3dQCD/region08' + '/' + meName_Mt
                )
            )                                                     
        )
    ),

    fit = cms.PSet(
        algorithm = cms.PSet(
            pluginName = cms.string("fitTauIdEffZtoMuTauSideband_tauIdPassed"),
            pluginType = cms.string("TemplateFitAdapter_TFractionFitter")
            #pluginType = cms.string("TemplateFitAdapter_RooFit")
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
  dqmDirectory_QCD_template_pure_tauIdFailed + '/' + meName_Mt
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.QCD.norm.initial = cms.double(1000.)
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.Ztautau.templates.Mt.meName = \
  dqmDirectory_Ztautau_template_tauIdFailed + '/' + meName_Mt
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.Ztautau.norm.initial = cms.double(100.)
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.data.distributions.Mt.meName = \
  'harvested/knownBgCorr/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit' \
 + '/' + 'tauIdEffHistogramsComb3dQCD/region07' + '/' + meName_Mt                                                                 
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

meOptionsSeparator = "#"
meOptionsNumWeighted = "".join([meOptionsSeparator, "a1", meOptionsSeparator, "s1"])
meOptionsErrWeighted = "".join([meOptionsSeparator, "a2", meOptionsSeparator, "s1"])

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
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binContent_region2' + meOptionsNumWeighted
            ),
            meName_numeratorErr = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binError_region2' + meOptionsErrWeighted
            ),
            meName_denominator = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binContent_region6' + meOptionsNumWeighted
            ),
            meName_denominatorErr = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binError_region6' + meOptionsErrWeighted
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
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binContent_region6' + meOptionsNumWeighted
            ),
            meName_numeratorErr = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binError_region6' + meOptionsErrWeighted
            ),
            meName_denominator = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binContent_region12' + meOptionsNumWeighted
            ),
            meName_denominatorErr = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binError_region12' + meOptionsErrWeighted
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
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binContent_region1' + meOptionsNumWeighted
            ),
            meName_numeratorErr = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binError_region1' + meOptionsErrWeighted
            ),
            meName_denominator = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binContent_region5' + meOptionsNumWeighted
            ),
            meName_denominatorErr = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binError_region5' + meOptionsErrWeighted
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
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binContent_region5' + meOptionsNumWeighted
            ),
            meName_numeratorErr = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binError_region5' + meOptionsErrWeighted
            ),
            meName_denominator = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binContent_region11' + meOptionsNumWeighted
            ),
            meName_denominatorErr = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binError_region11' + meOptionsErrWeighted
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
# compute sum of "known" backgrounds plus
# QCD and W + jets background contributions extrapolated from sidebands
#--------------------------------------------------------------------------------

process.compBgSumSideband_tauIdPassed = cms.EDAnalyzer("DQMHistAdder",
    config = cms.VPSet(                                                               
        cms.PSet(
            input = cms.VPSet( 
                cms.PSet( 
                    meName = cms.string(
                      dqmDirectory_Zmumu_control + '/' + 'tauIdEffBinningResults2regions/binContent_region2' + meOptionsNumWeighted
                    ),
                    meName_err = cms.string(
                      dqmDirectory_Zmumu_control + '/' + 'tauIdEffBinningResults2regions/binError_region2' + meOptionsErrWeighted
                    )
                ),
                cms.PSet(
                    meName = cms.string(
                      dqmDirectory_TTplusJets_control + '/' + 'tauIdEffBinningResults2regions/binContent_region2' + meOptionsNumWeighted
                    ),
                    meName_err = cms.string(
                      dqmDirectory_TTplusJets_control + '/' + 'tauIdEffBinningResults2regions/binError_region2' + meOptionsErrWeighted
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
                      dqmDirectory_Zmumu_control + '/' + 'tauIdEffBinningResults2regions/binContent_region1' + meOptionsNumWeighted
                    ),
                    meName_err = cms.string(
                      dqmDirectory_Zmumu_control + '/' + 'tauIdEffBinningResults2regions/binError_region1' + meOptionsErrWeighted
                    )
                ),
                cms.PSet(
                    meName = cms.string(
                      dqmDirectory_TTplusJets_control + '/' + 'tauIdEffBinningResults2regions/binContent_region1' + meOptionsNumWeighted
                    ),
                    meName_err = cms.string(
                      dqmDirectory_TTplusJets_control + '/' + 'tauIdEffBinningResults2regions/binError_region1' + meOptionsErrWeighted
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
              dqmDirectory_Data + '/' + 'tauIdEffBinningResults2regions/binContent_region2' + meOptionsNumWeighted
            ),
            meName_minuendErr = cms.string(
              dqmDirectory_Data + '/' + 'tauIdEffBinningResults2regions/binError_region2' + meOptionsErrWeighted
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
              dqmDirectory_Data + '/' + 'tauIdEffBinningResults2regions/binContent_region1' + meOptionsNumWeighted
            ),
            meName_minuendErr = cms.string(
              dqmDirectory_Data + '/' + 'tauIdEffBinningResults2regions/binError_region1' + meOptionsErrWeighted
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
              dqmDirectory_Ztautau_control + '/' + 'tauIdEffBinningResults2regions/binContent_region2' + meOptionsNumWeighted
            ),
            meName_passedErr = cms.string(
              dqmDirectory_Ztautau_control + '/' + 'tauIdEffBinningResults2regions/binError_region2' + meOptionsErrWeighted
            ),
            meName_failed = cms.string(
              dqmDirectory_Ztautau_control + '/' + 'tauIdEffBinningResults2regions/binContent_region1' + meOptionsNumWeighted
            ),
            meName_failedErr = cms.string(
              dqmDirectory_Ztautau_control + '/' + 'tauIdEffBinningResults2regions/binError_region1' + meOptionsErrWeighted
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
            dqmDirectory_input = cms.string(dqmDirectory_Ztautau_control_tauIdPassed),
            dqmDirectory_output = cms.string(dqmDirectory_Ztautau_controlNormalized_tauIdPassed),
        ),
        cms.PSet(
            dqmDirectory_input = cms.string(dqmDirectory_WplusJets_control_tauIdPassed),
            dqmDirectory_output = cms.string(dqmDirectory_WplusJets_controlNormalized_tauIdPassed),
        ),
        cms.PSet(
            dqmDirectory_input = cms.string(dqmDirectory_QCD_control_tauIdPassed),
            dqmDirectory_output = cms.string(dqmDirectory_QCD_controlNormalized_tauIdPassed),
        )
    ),
    norm = cms.double(1.)                                                             
)

process.compScaledDistributions_tauIdPassed = cms.EDAnalyzer("DQMHistScaler",
    config = cms.VPSet(
        cms.PSet(
            dqmDirectory_input = cms.string(dqmDirectory_Ztautau_controlNormalized_tauIdPassed),
            dqmDirectory_output = cms.string(dqmDirectory_Ztautau_controlScaled_tauIdPassed),
            meName_scaleFactor = cms.string(
              process.fitTauIdEffZtoMuTauSideband_tauIdPassed.output.fitResults.dqmDirectory.value() + '/' + meName_Ztautau_norm
            ),
            meName_scaleFactorErr = cms.string(
              process.fitTauIdEffZtoMuTauSideband_tauIdPassed.output.fitResults.dqmDirectory.value() + '/' + meName_Ztautau_normErr
            ),
            meType = cms.string("real")
        ),
        cms.PSet(
            dqmDirectory_input = cms.string(dqmDirectory_WplusJets_controlNormalized_tauIdPassed),
            dqmDirectory_output = cms.string(dqmDirectory_WplusJets_controlScaled_tauIdPassed),
            meName_scaleFactor = cms.string(
              process.fitTauIdEffZtoMuTauSideband_tauIdPassed.output.fitResults.dqmDirectory.value() + '/' + meName_WplusJets_norm
            ),
            meName_scaleFactorErr = cms.string(
              process.fitTauIdEffZtoMuTauSideband_tauIdPassed.output.fitResults.dqmDirectory.value() + '/' + meName_WplusJets_normErr
            ),
            meType = cms.string("real")
        ),
        cms.PSet(
            dqmDirectory_input = cms.string(dqmDirectory_QCD_controlNormalized_tauIdPassed),
            dqmDirectory_output = cms.string(dqmDirectory_QCD_controlScaled_tauIdPassed),
            meName_scaleFactor = cms.string(
              process.fitTauIdEffZtoMuTauSideband_tauIdPassed.output.fitResults.dqmDirectory.value() + '/' + meName_QCD_norm
            ),
            meName_scaleFactorErr = cms.string(
              process.fitTauIdEffZtoMuTauSideband_tauIdPassed.output.fitResults.dqmDirectory.value() + '/' + meName_QCD_normErr
            ),
            meType = cms.string("real")
        )
    )
)

##process.addScaledDistributions_tauIdPassed = cms.EDAnalyzer("DQMHistAdder",
##    smSum = cms.PSet(
##        dqmDirectories_input = cms.vstring(
##            dqmDirectory_Ztautau_controlScaled_tauIdPassed,
##            dqmDirectory_WplusJets_controlScaled_tauIdPassed,
##            dqmDirectory_QCD_controlScaled_tauIdPassed,
##            dqmDirectory_Zmumu_control_tauIdPassed,
##            dqmDirectory_TTplusJets_control_tauIdPassed
##        ),
##        dqmDirectory_output = cms.string(dqmDirectory_smSum_control_tauIdPassed)
##    )
##)

process.plotScaledDistributions_tauIdPassed = copy.deepcopy(plotTauIdEffZtoMuTau)
process.plotScaledDistributions_tauIdPassed.processes.Ztautau.dqmDirectory = dqmDirectory_Ztautau_controlScaled_tauIdPassed
process.plotScaledDistributions_tauIdPassed.processes.WplusJets.dqmDirectory = dqmDirectory_WplusJets_controlScaled_tauIdPassed
process.plotScaledDistributions_tauIdPassed.processes.qcdSum.dqmDirectory = dqmDirectory_QCD_controlScaled_tauIdPassed
process.plotScaledDistributions_tauIdPassed.processes.Zmumu.dqmDirectory = dqmDirectory_Zmumu_control_tauIdPassed
process.plotScaledDistributions_tauIdPassed.processes.TTplusJets.dqmDirectory = dqmDirectory_TTplusJets_control_tauIdPassed
drawJobConfigurator_TauIdEffZtoMuTau_tauIdPassed = drawJobConfigurator(
    template = plots_TauIdEffZtoMuTau,
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
   ##* process.addScaledDistributions_tauIdPassed
   * process.plotScaledDistributions_tauIdPassed
)

process.compNormalizedDistributions_tauIdFailed = copy.deepcopy(process.compNormalizedDistributions_tauIdPassed)
process.compNormalizedDistributions_tauIdFailed.config[0].dqmDirectory_input = dqmDirectory_Ztautau_control_tauIdFailed
process.compNormalizedDistributions_tauIdFailed.config[0].dqmDirectory_output = dqmDirectory_Ztautau_controlNormalized_tauIdFailed
process.compNormalizedDistributions_tauIdFailed.config[1].dqmDirectory_input = dqmDirectory_WplusJets_control_tauIdFailed
process.compNormalizedDistributions_tauIdFailed.config[1].dqmDirectory_output = dqmDirectory_WplusJets_controlNormalized_tauIdFailed
process.compNormalizedDistributions_tauIdFailed.config[2].dqmDirectory_input = dqmDirectory_QCD_control_tauIdFailed
process.compNormalizedDistributions_tauIdFailed.config[2].dqmDirectory_output = dqmDirectory_QCD_controlNormalized_tauIdFailed

process.compScaledDistributions_tauIdFailed = copy.deepcopy(process.compScaledDistributions_tauIdPassed)
process.compScaledDistributions_tauIdFailed.config[0].dqmDirectory_input = dqmDirectory_Ztautau_controlNormalized_tauIdFailed
process.compScaledDistributions_tauIdFailed.config[0].dqmDirectory_output = dqmDirectory_Ztautau_controlScaled_tauIdFailed
process.compScaledDistributions_tauIdFailed.config[0].meName_scaleFactor = \
  process.fitTauIdEffZtoMuTauSideband_tauIdFailed.output.fitResults.dqmDirectory.value() + '/' + meName_Ztautau_norm
process.compScaledDistributions_tauIdFailed.config[0].meName_scaleFactorErr = \
  process.fitTauIdEffZtoMuTauSideband_tauIdFailed.output.fitResults.dqmDirectory.value() + '/' + meName_Ztautau_normErr
process.compScaledDistributions_tauIdFailed.config[1].dqmDirectory_input = dqmDirectory_WplusJets_controlNormalized_tauIdFailed
process.compScaledDistributions_tauIdFailed.config[1].dqmDirectory_output = dqmDirectory_WplusJets_controlScaled_tauIdFailed
process.compScaledDistributions_tauIdFailed.config[1].meName_scaleFactor = \
  process.fitTauIdEffZtoMuTauSideband_tauIdFailed.output.fitResults.dqmDirectory.value() + '/' + meName_WplusJets_norm
process.compScaledDistributions_tauIdFailed.config[1].meName_scaleFactorErr = \
  process.fitTauIdEffZtoMuTauSideband_tauIdFailed.output.fitResults.dqmDirectory.value() + '/' + meName_WplusJets_normErr
process.compScaledDistributions_tauIdFailed.config[2].dqmDirectory_input = dqmDirectory_QCD_controlNormalized_tauIdFailed
process.compScaledDistributions_tauIdFailed.config[2].dqmDirectory_output = dqmDirectory_QCD_controlScaled_tauIdFailed
process.compScaledDistributions_tauIdFailed.config[2].meName_scaleFactor = \
  process.fitTauIdEffZtoMuTauSideband_tauIdFailed.output.fitResults.dqmDirectory.value() + '/' + meName_QCD_norm
process.compScaledDistributions_tauIdFailed.config[2].meName_scaleFactorErr = \
  process.fitTauIdEffZtoMuTauSideband_tauIdFailed.output.fitResults.dqmDirectory.value() + '/' + meName_QCD_normErr

##process.addScaledDistributions_tauIdFailed = cms.EDAnalyzer("DQMHistAdder",
##    smSum = cms.PSet(
##        dqmDirectories_input = cms.vstring(
##            dqmDirectory_Ztautau_controlScaled_tauIdFailed,
##            dqmDirectory_WplusJets_controlScaled_tauIdFailed,
##            dqmDirectory_QCD_controlScaled_tauIdFailed,
##            dqmDirectory_Zmumu_control_tauIdFailed,
##            dqmDirectory_TTplusJets_control_tauIdFailed
##        ),
##        dqmDirectory_output = cms.string(dqmDirectory_smSum_control_tauIdFailed)
##    )
##)

process.plotScaledDistributions_tauIdFailed = copy.deepcopy(process.plotScaledDistributions_tauIdPassed)
process.plotScaledDistributions_tauIdFailed.processes.Ztautau.dqmDirectory = dqmDirectory_Ztautau_controlScaled_tauIdFailed
process.plotScaledDistributions_tauIdFailed.processes.WplusJets.dqmDirectory = dqmDirectory_WplusJets_controlScaled_tauIdFailed
process.plotScaledDistributions_tauIdFailed.processes.qcdSum.dqmDirectory = dqmDirectory_QCD_controlScaled_tauIdFailed
process.plotScaledDistributions_tauIdFailed.processes.Zmumu.dqmDirectory = dqmDirectory_Zmumu_control_tauIdFailed
process.plotScaledDistributions_tauIdFailed.processes.TTplusJets.dqmDirectory = dqmDirectory_TTplusJets_control_tauIdFailed
process.plotScaledDistributions_tauIdFailed.indOutputFileName = cms.string('computeTauIdEffZtoMuTauCombinedFit_tauIdFailed_#PLOT#.png')

process.makeControlPlotsTauIdEffZtoMuTauSideband_tauIdFailed = cms.Sequence(
    process.compNormalizedDistributions_tauIdFailed
   * process.compScaledDistributions_tauIdFailed
   ##* process.addScaledDistributions_tauIdFailed
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

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

process.p = cms.Path(
    process.loadTauIdEffZtoMuTau
   #+ process.dumpDQMStore
   + process.plotTauIdEffZtoMuTauCombinedFit
   + process.fitTauIdEffZtoMuTauSideband
   #+ process.dumpDQMStore
   + process.extrapolateTauIdEffZtoMuTauSideband
   + process.compBgSumSideband
   + process.compBgCorrectionsSideband
   + process.dumpBinErrorsTauIdEffZtoMuTau
   + process.makeControlPlotsTauIdEffZtoMuTauSideband
   + process.saveFitResultsTauIdEffZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
