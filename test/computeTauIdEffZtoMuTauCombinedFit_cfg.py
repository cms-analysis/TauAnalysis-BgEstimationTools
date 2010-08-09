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

dqmDirectory_WplusJets_sideband = \
  'harvested/WplusJets/TauIdEffAnalyzerZtoMuTauCombinedFitWplusJets/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_WplusJets_sidebandOS_tauIdPassed = dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffHistogramsComb3dWplusJets/region02'
dqmDirectory_WplusJets_sidebandOS_tauIdFailed = dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffHistogramsComb3dWplusJets/region01'
dqmDirectory_WplusJets_sidebandSS_tauIdPassed = dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffHistogramsComb3dWplusJets/region06'
dqmDirectory_WplusJets_sidebandSS_tauIdFailed = dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffHistogramsComb3dWplusJets/region05'
dqmDirectory_WplusJets_template = \
  'harvested/WplusJets/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_WplusJets_template_tauIdPassed = dqmDirectory_WplusJets_template + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
dqmDirectory_WplusJets_template_tauIdFailed = dqmDirectory_WplusJets_template + '/' + 'tauIdEffHistogramsComb3dQCD/region08'

dqmDirectory_QCD_sideband = \
  'harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_QCD_sidebandOS_tauIdPassed = dqmDirectory_QCD_sideband + '/' + 'tauIdEffHistogramsComb3dQCD/region06'
dqmDirectory_QCD_sidebandOS_tauIdFailed = dqmDirectory_QCD_sideband + '/' + 'tauIdEffHistogramsComb3dQCD/region05'
dqmDirectory_QCD_sidebandSS_tauIdPassed = dqmDirectory_QCD_sideband + '/' + 'tauIdEffHistogramsComb3dQCD/region12'
dqmDirectory_QCD_sidebandSS_tauIdFailed = dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/region11'
dqmDirectory_QCD_template_pure = \
  'harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_QCD_template_pure_tauIdPassed = dqmDirectory_QCD_template + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
dqmDirectory_QCD_template_pure_tauIdFailed = dqmDirectory_QCD_template + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
dqmDirectory_QCD_template_data = \
  'harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_QCD_template_data_tauIdPassed = dqmDirectory_QCD_template + '/' + 'tauIdEffHistogramsComb3dQCD/region12'
dqmDirectory_QCD_template_data_tauIdFailed = dqmDirectory_QCD_template + '/' + 'tauIdEffHistogramsComb3dQCD/region12'

dqmDirectory_Data = 'harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
dqmDirectory_Data_tauIdPassed = dqmDirectory_Data_all + '/' + 'tauIdEffBinningResults2regions/region02'
dqmDirectory_Data_tauIdFailed = dqmDirectory_Data_all + '/' + 'tauIdEffBinningResults2regions/region01'

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
        'harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit' \
          + '/' + 'tauIdEffHistograms2regions/region02' + '/' + meName_Mt
    ],
    name = "QCD_Mt",
    title = "M_{T} in QCD Background"
)

plotTauIdEffZtoMuTau = cms.EDAnalyzer("DQMHistPlotter",
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

    drawJobs = cms.PSet(),

    canvasSizeX = cms.int32(800),
    canvasSizeY = cms.int32(640),                         

    outputFilePath = cms.string('./plots/'),
    #outputFileName = cms.string('plotsTauIdEffZtoMuTauCombinedFit.ps')
    indOutputFileName = cms.string('plotTauIdEffZtoMuTauCombinedFit_#PLOT#.png')
)

process.plotTauIdEffZtoMuTauCombinedFit_Mt = copy.deepcopy(plotTauIdEffZtoMuTau)
process.plotTauIdEffZtoMuTauCombinedFit_Mt.drawJobs = drawTemplateHistConfiguratorTauIdEffZtoMuTau_Mt.configure()

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

process.compKnownBgCorrections = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meNameMinuend = cms.string(
              'harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region07'
            ),
            meNameSubtrahend = cms.string(
              'harvested/knownBgSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region07'
            ),
            meNameDifference = cms.string(
              'harvested/knownBgCorr/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region07'
            )
        ),
        cms.PSet(
            meNameMinuend = cms.string(
              'harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
            ),
            meNameSubtrahend = cms.string(
              'harvested/knownBgSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
            ),
            meNameDifference = cms.string(
              'harvested/knownBgCorr/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
            )
        )
    )
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
        )
    ),

    # use "pseudo" data-samples consisting of all Monte Carlo processes for testing                      
    data = cms.PSet(
        distributions = cms.PSet(
            Mt = cms.PSet(
                meName = cms.string(
                  'harvested/knownBgCorr/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
                 + '/' + 'tauIdEffHistogramsComb3dQCD/region08' + '/' + meName_Mt
                )
            )                                                     
        )
    ),

    fit = cms.PSet(
        algorithm = cms.PSet(
            pluginName = cms.string("fitTauIdEffZtoMuTauAlgorithm_tauIdPassed"),
            pluginType = cms.string("TemplateFitAdapter_TFractionFitter")
            #pluginType = cms.string("TemplateFitAdapter_RooFit")
        ),
        variables = cms.PSet(
            Mt = cms.PSet(
                name = cms.string("Mt"),
                title = cms.string("M_{T}"),
                min = cms.double(0.),
                max = cms.double(1000.)
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
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.QCD.templates.diTauMass.meName = \
  dqmDirectory_QCD_template_pure_tauIdFailed + '/' + meName_Mt
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.processes.QCD.norm.initial = cms.double(1000.)
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.data.distributions.Mt.meName = \
  'harvested/knownBgCorr/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
 + '/' + 'tauIdEffHistogramsComb3dQCD/region08' + '/' + meName_Mt                                                                 
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.output.controlPlots.fileName = \
  "./plots/fitTauIdEffZtoMuTauSideband_tauIdFailed_#PLOT#.png"
process.fitTauIdEffZtoMuTauSideband_tauIdFailed.output.fitResults.dqmDirectory = \
  "fitTauIdEffZtoMuTauSideband_tauIdFailed/fitResults"

process.fitTauIdEffZtoMuTauSideband = cms.Sequence(
    process.fitTauIdEffZtoMuTauSideband_tauIdPassed
   * process.fitTauIdEffZtoMuTauSideband_tauIdFailed
)

#--------------------------------------------------------------------------------
# extrapolate contributions of QCD and W + jet backgrounds
# from sideband into signal region
#--------------------------------------------------------------------------------

meOptionsSeparator = "#"
meOptionsNumWeighted = "".join([meOptionsSeparator, "a1", meOptionsSeparator, "s1"])
meOptionsErrWeighted = "".join([meOptionsSeparator, "a2", meOptionsSeparator, "s1"])

process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed = cms.EDAnalyzer("DQMHistScaler",
    config = cms.VPSet(
        cms.PSet(
            meName_input = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdPassed/fitResults/norm_WplusJets'
            ),
            meName_inputErr = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdPassed/fitResults/normErr_WplusJets'
            ),
            meNameNumerator = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binContent_region2' + meOptionsNumWeighted
            ),
            meNameNumeratorErr = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binError_region2' + meOptionsErrWeighted
            ),
            meNameDenominator = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binContent_region6' + meOptionsNumWeighted
            ),
            meNameDenominatorErr = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binError_region6' + meOptionsErrWeighted
            ),
            meType = cms.string("real"),
            meName_output = cms.string('extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/norm_WplusJets'),
            meName_outputErr = cms.string('extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/normErr_WplusJets')
        ),
        cms.PSet(
            meName_input = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdPassed/fitResults/norm_QCD'
            ),
            meName_inputErr = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdPassed/fitResults/normErr_QCD'
            ),
            meNameNumerator = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binContent_region6' + meOptionsNumWeighted
            ),
            meNameNumeratorErr = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binError_region6' + meOptionsErrWeighted
            ),
            meNameDenominator = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binContent_region12' + meOptionsNumWeighted
            ),
            meNameDenominatorErr = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binError_region12' + meOptionsErrWeighted
            ),
            meType = cms.string("real"),
            meName_output = cms.string('extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/norm_QCD'),
            meName_outputErr = cms.string('extrapolateTauIdEffZtoMuTauSideband_tauIdPassed/normErr_QCD')
        )
    )                                                                             
)

process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed = cms.EDAnalyzer("DQMHistScaler",
    config = cms.VPSet(
        cms.PSet(
            meName_input = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdFailed/fitResults/norm_WplusJets'
            ),
            meName_inputErr = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdFailed/fitResults/normErr_WplusJets'
            ),
            meNameNumerator = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binContent_region1' + meOptionsNumWeighted
            ),
            meNameNumeratorErr = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binError_region1' + meOptionsErrWeighted
            ),
            meNameDenominator = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binContent_region5' + meOptionsNumWeighted
            ),
            meNameDenominatorErr = cms.string(
              dqmDirectory_WplusJets_sideband + '/' + 'tauIdEffBinningResultsComb3dWplusJets/binError_region5' + meOptionsErrWeighted
            ),
            meType = cms.string("real"),
            meName_output = cms.string('extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/norm_WplusJets'),
            meName_outputErr = cms.string('extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/normErr_WplusJets')
        ),
        cms.PSet(
            meName_input = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdFailed/fitResults/norm_QCD'
            ),
            meName_inputErr = cms.string(
              'fitTauIdEffZtoMuTauSideband_tauIdFailed/fitResults/normErr_QCD'
            ),
            meNameNumerator = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binContent_region5' + meOptionsNumWeighted
            ),
            meNameNumeratorErr = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binError_region5' + meOptionsErrWeighted
            ),
            meNameDenominator = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binContent_region11' + meOptionsNumWeighted
            ),
            meNameDenominatorErr = cms.string(
              dqmDirectory_QCD_sideband + '/' + 'tauIdEffBinningResultsComb3dQCD/binError_region11' + meOptionsErrWeighted
            ),
            meType = cms.string("real"),
            meName_output = cms.string('extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/norm_QCD'),
            meName_outputErr = cms.string('extrapolateTauIdEffZtoMuTauSideband_tauIdFailed/normErr_QCD')
        )
    )                                                                             
)

process.extrapolateTauIdEffZtoMuTauSideband = cms.Sequence(
    process.extrapolateTauIdEffZtoMuTauSideband_tauIdPassed
   * process.extrapolateTauIdEffZtoMuTauSideband_tauIdFailed
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

process.compKnownBgCorrections = cms.EDAnalyzer("DQMHistSubtractor",
    config = cms.VPSet(
        cms.PSet(
            meNameMinuend = cms.string(
              'harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region07'
            ),
            meNameSubtrahend = cms.string(
              'harvested/knownBgSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region07'
            ),
            meNameDifference = cms.string(
              'harvested/knownBgCorr/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region07'
            )
        ),
        cms.PSet(
            meNameMinuend = cms.string(
              'harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
            ),
            meNameSubtrahend = cms.string(
              'harvested/knownBgSum/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
            ),
            meNameDifference = cms.string(
              'harvested/knownBgCorr/TauIdEffAnalyzerZtoMuTauCombinedFitQCD/afterUniqueMuonCandidateCutTauIdEffZtoMuTauCombinedFit'
             + '/' + 'tauIdEffHistogramsComb3dQCD/region08'
            )
        )
    )
)

#--------------------------------------------------------------------------------
# compute tau id. efficiencies
# (using binomial errors for events statistics expected in data)
#--------------------------------------------------------------------------------

meName_norm = 'Ztautau/norm/value#a1#s1'

process.dumpBinErrorsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMBinErrorCalculator",
    config = cms.VPSet(
        cms.PSet(
            meName_passed = cms.string(dqmDirectory_Ztautau_all + '/' + 'tauIdEffBinningResults2regions/binContent_region2#a1#s1'),
            meName_failed = cms.string(dqmDirectory_Ztautau_all + '/' + 'tauIdEffBinningResults2regions/binContent_region1#a1#s1'),
            label = cms.string("Tau id., true")
        ),
        cms.PSet(
            meName_passed = cms.string(process.fitTauIdEffZtoMuTau_tauIdPassed.output.fitResults.dqmDirectory.value() + '/' + meName_norm),
            meName_failed = cms.string(process.fitTauIdEffZtoMuTau_tauIdFailed.output.fitResults.dqmDirectory.value() + '/' + meName_norm),
            label = cms.string("Tau id., fitted")
        )
    )
)

MAKE CONTROL PLOTS MUON pT, MUON ETA, VIS MASS. TAU PT, TAU ETA, TRANSMASS

process.saveFitResultsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('fitTauIdEffZtoMuTau_results.root'),
    outputCommands = cms.vstring('drop harvested/*')
)

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

process.p = cms.Path(
    process.loadTauIdEffZtoMuTau
  #+ process.dumpDQMStore
   + process.plotTemplateHistTauIdEffZtoMuTau
   + process.fitTauIdEffZtoMuTau
   + process.dumpBinErrorsTauIdEffZtoMuTau 
   + process.saveFitResultsTauIdEffZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
