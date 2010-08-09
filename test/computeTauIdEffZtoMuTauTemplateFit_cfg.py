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

processName = 'computeTauIdEffZtoMuTauTemplateFit'
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

region_tauIdPassed = 'region02'
region_tauIdFailed = 'region01'

dqmDirectory_Ztautau_all = 'harvested/Ztautau/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau'
##dqmDirectory_Ztautau_all = 'harvested/Ztautau/TauIdEffAnalyzerZtoMuTau_relMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau'
dqmDirectory_Ztautau_tauIdPassed = dqmDirectory_Ztautau_all + '/' + 'tauIdEffHistograms2regions' + '/' + region_tauIdPassed 
dqmDirectory_Ztautau_tauIdFailed = dqmDirectory_Ztautau_all + '/' + 'tauIdEffHistograms2regions' + '/' + region_tauIdFailed 
dqmDirectory_Ztautau_template = dqmDirectory_Ztautau_all

dqmDirectory_WplusJets_all = 'harvested/WplusJets/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau'
##dqmDirectory_WplusJets_all = 'harvested/WplusJets/TauIdEffAnalyzerZtoMuTau_relMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau'
dqmDirectory_WplusJets_tauIdPassed = dqmDirectory_WplusJets_all + '/' + 'tauIdEffHistograms2regions' + '/' + region_tauIdPassed
dqmDirectory_WplusJets_tauIdFailed = dqmDirectory_WplusJets_all + '/' + 'tauIdEffHistograms2regions' + '/' + region_tauIdFailed
#dqmDirectory_WplusJets_template = 'harvested/smSum/BgEstTemplateAnalyzer_WplusJetsEnriched/afterDiMuonVetoBgEstWplusJetsEnriched'
dqmDirectory_WplusJets_template = dqmDirectory_WplusJets_all 

dqmDirectory_QCD_all = 'harvested/qcdSum/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau'
##dqmDirectory_QCD_all = 'harvested/qcdSum/TauIdEffAnalyzerZtoMuTau_relMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau'
dqmDirectory_QCD_tauIdPassed = dqmDirectory_QCD_all + '/' + 'tauIdEffHistograms2regions' + '/' + region_tauIdPassed
dqmDirectory_QCD_tauIdFailed = dqmDirectory_QCD_all + '/' + 'tauIdEffHistograms2regions' + '/' + region_tauIdFailed 
dqmDirectory_QCD_template = 'harvested/smSum/BgEstTemplateAnalyzer_QCDenriched_reweighted/afterDiMuonVetoBgEstQCDenriched'
##dqmDirectory_QCD_template = 'harvested/smSum/BgEstTemplateAnalyzer_QCDenriched/afterDiMuonVetoBgEstQCDenriched'
#dqmDirectory_QCD_template = dqmDirectory_QCD_all

dqmDirectory_Data_all = 'harvested/smSum/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau'
##dqmDirectory_Data_all = 'harvested/smSum/TauIdEffAnalyzerZtoMuTau_relMuonIsolation/afterUniqueMuonCandidateCutTauIdEffZtoMuTau'
dqmDirectory_Data_tauIdPassed = dqmDirectory_Data_all + '/' + 'tauIdEffHistograms2regions' + '/' + region_tauIdPassed 
dqmDirectory_Data_tauIdFailed = dqmDirectory_Data_all + '/' + 'tauIdEffHistograms2regions' + '/' + region_tauIdFailed

#--------------------------------------------------------------------------------
# define names of histograms used in fit
#--------------------------------------------------------------------------------

meName_muonPt = 'TauIdEffSpecificQuantities/MuonPt'
meName_muonAbsEta = 'TauIdEffSpecificQuantities/MuonAbsEta'
meName_muonPtVsAbsEta = 'TauIdEffSpecificQuantities/MuonPtVsAbsEta'
meName_diTauMass = 'TauIdEffSpecificQuantities/DiTauCollinearApproxMassFromJetP4'

#--------------------------------------------------------------------------------
# plot template histograms of "pure" Monte Carlo processes
# compared to the shapes determined by background enriched regions in (pseudo)Data
#--------------------------------------------------------------------------------

drawJobTemplateHist_muonPt = copy.deepcopy(drawJobTemplateHist)
drawJobTemplateHist_muonPt.plots[0].process = cms.string('tauIdPassed')
drawJobTemplateHist_muonPt.plots[0].drawOptionEntry = cms.string('tauIdPassed')
drawJobTemplateHist_muonPt.plots[0].legendEntry = cms.string('Tau Id. passed')
drawJobTemplateHist_muonPt.plots[1].process = cms.string('tauIdFailed')
drawJobTemplateHist_muonPt.plots[1].drawOptionEntry = cms.string('tauIdFailed')
drawJobTemplateHist_muonPt.plots[1].legendEntry = cms.string('Tau Id. failed')
drawJobTemplateHist_muonPt.plots[2].process = cms.string('template')
drawJobTemplateHist_muonPt.plots[2].drawOptionEntry = cms.string('template')
drawJobTemplateHist_muonPt.plots[2].legendEntry = cms.string('Template')
drawJobTemplateHist_muonPt.title = cms.string('Muon P_{T}')
drawJobTemplateHist_muonPt.xAxis = cms.string('Pt')
drawJobTemplateHist_muonPt.yAxis = cms.string('numEntries_log')

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt = drawTemplateHistConfigurator(
    template = drawJobTemplateHist_muonPt
)

drawJobTemplateHist_muonAbsEta = copy.deepcopy(drawJobTemplateHist_muonPt)
drawJobTemplateHist_muonAbsEta.title = cms.string('Muon |#eta|')
drawJobTemplateHist_muonAbsEta.xAxis = cms.string('Eta')

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta = drawTemplateHistConfigurator(
    template = drawJobTemplateHist_muonAbsEta
)

drawJobTemplateHist_diTauMass = copy.deepcopy(drawJobTemplateHist_muonPt)
drawJobTemplateHist_diTauMass.title = cms.string('M(Muon + Tau)')
drawJobTemplateHist_diTauMass.xAxis = cms.string('Mass')

drawTemplateHistConfiguratorTauIdEffZtoMuTau_diTauMass = drawTemplateHistConfigurator(
    template = drawJobTemplateHist_diTauMass
)

#--------------------------------------------------------------------------------
# define draw jobs for Z --> tau+ tau- signal
#
# NOTE: backgrounds contributing to Z --> mu+ mu- sample from which
#       template histogram for Ztautau signal process is determined using MCEmbeddingTools
#       not included yet, so use "pure" sample as approximation for sample determined in "deta"
#       for the time being...
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.add(
    meNames = [
        dqmDirectory_Ztautau_tauIdPassed + '/' + meName_muonPt,
        dqmDirectory_Ztautau_tauIdFailed + '/' + meName_muonPt,
        dqmDirectory_Ztautau_all + '/' + meName_muonPt
    ],
    name = "Ztautau_muonPt",
    title = "P_{T}^{#mu} in Z #rightarrow #tau^{+} #tau^{-} Signal" 
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.add(
    meNames = [
        dqmDirectory_Ztautau_tauIdPassed + '/' + meName_muonAbsEta,
        dqmDirectory_Ztautau_tauIdFailed + '/' + meName_muonAbsEta,
        dqmDirectory_Ztautau_all + '/' + meName_muonAbsEta
    ],
    name = "Ztautau_muonAbsEta",
    title = "|#eta_{#mu}| in Z #rightarrow #tau^{+} #tau^{-} Signal" 
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_diTauMass.add(
    meNames = [
        dqmDirectory_Ztautau_tauIdPassed + '/' + meName_diTauMass,
        dqmDirectory_Ztautau_tauIdFailed + '/' + meName_diTauMass,
        dqmDirectory_Ztautau_all + '/' + meName_diTauMass
    ],
    name = "Ztautau_diTauMass",
    title = "M^{#mu #tau} in Z #rightarrow #tau^{+} #tau^{-} Signal" 
)

#--------------------------------------------------------------------------------
# define draw jobs for Z --> mu+ mu- background
#--------------------------------------------------------------------------------

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.add(
##    meNames = [
##        dqmDirectory_Zmumu_tauIdPassed + '/' + meName_muonPt,
##        dqmDirectory_Zmumu_tauIdFailed + '/' + meName_muonPt,
##        #dqmDirectory_Zmumu_all + '/' + meName_muonPt
##        dqmDirectory_Zmumu_template + '/' + meName_muonPt
##    ],
##    name = "Zmumu_muonPt",
##    title = "P_{T}^{#mu} in Z #rightarrow #mu^{+} #mu^{-} Background" 
##)

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.add(
##    meNames = [
##        dqmDirectory_Zmumu_tauIdPassed + '/' + meName_muonAbsEta,
##        dqmDirectory_Zmumu_tauIdFailed + '/' + meName_muonAbsEta,
##        #dqmDirectory_Zmumu_all + '/' + meName_muonAbsEta
##        dqmDirectory_Zmumu_template + '/' + meName_muonAbsEta
##    ],
##    name = "Zmumu_muonAbsEta",
##    title = "|#eta_{#mu}| in Z #rightarrow #mu^{+} #mu^{-} Background"
##)

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_diTauMass.add(
##    meNames = [
##        dqmDirectory_Zmumu_tauIdPassed + '/' + meName_diTauMass,
##        dqmDirectory_Zmumu_tauIdFailed + '/' + meName_diTauMass,
##        #dqmDirectory_Zmumu_all + '/' + meName_diTauMass
##        dqmDirectory_Zmumu_template + '/' + meName_diTauMass
##    ],
##    name = "Zmumu_diTauMass",
##    title = "M^{#mu #tau} in Z #rightarrow #mu^{+} #mu^{-} Background"
##)

#--------------------------------------------------------------------------------
# define draw jobs for W + jets background
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.add(
    meNames = [
        dqmDirectory_WplusJets_tauIdPassed + '/' + meName_muonPt,
        dqmDirectory_WplusJets_tauIdFailed + '/' + meName_muonPt,
        #dqmDirectory_WplusJets_all + '/' + meName_muonPt
        dqmDirectory_WplusJets_template + '/' + meName_muonPt
    ],
    name = "WplusJets_muonPt",
    title = "P_{T}^{#mu} in W + jets Background"
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.add(
    meNames = [
        dqmDirectory_WplusJets_tauIdPassed + '/' + meName_muonAbsEta,
        dqmDirectory_WplusJets_tauIdFailed + '/' + meName_muonAbsEta,
        #dqmDirectory_WplusJets_all + '/' + meName_muonAbsEta
        dqmDirectory_WplusJets_template + '/' + meName_muonAbsEta
    ],
    name = "WplusJets_muonAbsEta",
    title = "|#eta_{#mu}| in W + jets Background"
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_diTauMass.add(
   meNames = [
       dqmDirectory_WplusJets_tauIdPassed + '/' + meName_diTauMass,
       dqmDirectory_WplusJets_tauIdFailed + '/' + meName_diTauMass,
       #dqmDirectory_WplusJets_all + '/' + meName_diTauMass
       dqmDirectory_WplusJets_template + '/' + meName_diTauMass
   ],
   name = "WplusJets_diTauMass",
   title = "M^{#mu #tau} in W + jets Background"
)

#--------------------------------------------------------------------------------
# define draw jobs for TTbar + jets background
#--------------------------------------------------------------------------------

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.add(
##    meNames = [
##        dqmDirectory_TTplusJets_tauIdPassed + '/' + meName_muonPt,
##        dqmDirectory_TTplusJets_tauIdFailed + '/' + meName_muonPt,
##        #dqmDirectory_TTplusJets_all + '/' + meName_muonPt
##        dqmDirectory_TTplusJets_template + '/' + meName_muonPt
##    ],
##    name = "TTplusJets_muonPt",
##    title = "P_{T}^{#mu} in t#bar{t} + jets Background"
##)

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.add(
##    meNames = [
##        dqmDirectory_TTplusJets_tauIdPassed + '/' + meName_muonAbsEta,
##        dqmDirectory_TTplusJets_tauIdFailed + '/' + meName_muonAbsEta,
##        #dqmDirectory_TTplusJets_all + '/' + meName_muonAbsEta
##        dqmDirectory_TTplusJets_template + '/' + meName_muonAbsEta
##    ],
##    name = "TTplusJets_muonAbsEta",
##    title = "|#eta_{#mu}| in t#bar{t} + jets Background"
##)

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_diTauMass.add(
##   meNames = [
##       dqmDirectory_TTplusJets_tauIdPassed + '/' + meName_diTauMass,
##       dqmDirectory_TTplusJets_tauIdFailed + '/' + meName_diTauMass,
##       #dqmDirectory_TTplusJets_all + '/' + meName_diTauMass
##       dqmDirectory_TTplusJets_template + '/' + meName_diTauMass
##   ],
##   name = "TTplusJets_diTauMass",
##   title = "M^{#mu #tau} in t#bar{t} + jets Background"
##)

#--------------------------------------------------------------------------------
# define draw jobs for QCD background
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.add(
    meNames = [
        dqmDirectory_QCD_tauIdPassed + '/' + meName_muonPt,
        dqmDirectory_QCD_tauIdFailed + '/' + meName_muonPt,
        #dqmDirectory_QCD_all + '/' + meName_muonPt
        dqmDirectory_QCD_template + '/' + meName_muonPt
    ],
    name = "QCD_muonPt",
    title = "P_{T}^{#mu} in QCD Background"
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.add(
    meNames = [
        dqmDirectory_QCD_tauIdPassed + '/' + meName_muonAbsEta,
        dqmDirectory_QCD_tauIdFailed + '/' + meName_muonAbsEta,
        #dqmDirectory_QCD_all + '/' + meName_muonAbsEta
        dqmDirectory_QCD_template + '/' + meName_muonAbsEta
    ],
    name = "QCD_muonAbsEta",
    title = "|#eta_{#mu}| in QCD Background"
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_diTauMass.add(
  meNames = [
      dqmDirectory_QCD_tauIdPassed + '/' + meName_diTauMass,
      dqmDirectory_QCD_tauIdFailed + '/' + meName_diTauMass,
      #dqmDirectory_QCD_all + '/' + meName_diTauMass
      dqmDirectory_QCD_template + '/' + meName_diTauMass
  ],
  name = "QCD_diTauMass",
  title = "M^{#mu #tau} in QCD Background"
)

plotTemplateHistZtoMuTau = cms.EDAnalyzer("DQMHistPlotter",
    processes = cms.PSet(
        tauIdPassed = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist.plots[0].legendEntry,
            type = cms.string('smMC')
        ),
        tauIdFailed = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist.plots[1].legendEntry,
            type = cms.string('smMC')
        ),
        template = cms.PSet(
            dqmDirectory = cms.string(''),
            legendEntry = drawJobTemplateHist.plots[2].legendEntry,
            type = cms.string('smMC')
        )
    ),

    xAxes = cms.PSet(
        Pt = copy.deepcopy(xAxis_pt),
        Eta = copy.deepcopy(xAxis_eta),
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
        tauIdPassed = copy.deepcopy(drawOption_green_eff),
        tauIdFailed = copy.deepcopy(drawOption_lightBlue_eff),
        template = copy.deepcopy(drawOption_red_eff)
    ),

    drawJobs = cms.PSet(),

    canvasSizeX = cms.int32(800),
    canvasSizeY = cms.int32(640),                         

    outputFilePath = cms.string('./plots/'),
    #outputFileName = cms.string('plotsTemplateHistTauIdEffZtoMuTau.ps')
    indOutputFileName = cms.string('plotTemplateHistTauIdEffZtoMuTau_#PLOT#.png')
)

process.plotTemplateHistTauIdEffZtoMuTau_muonPt = copy.deepcopy(plotTemplateHistZtoMuTau)
process.plotTemplateHistTauIdEffZtoMuTau_muonPt.drawJobs = drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.configure()

process.plotTemplateHistTauIdEffZtoMuTau_muonAbsEta = copy.deepcopy(plotTemplateHistZtoMuTau)
process.plotTemplateHistTauIdEffZtoMuTau_muonAbsEta.drawJobs = drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.configure()

process.plotTemplateHistTauIdEffZtoMuTau_diTauMass = copy.deepcopy(plotTemplateHistZtoMuTau)
process.plotTemplateHistTauIdEffZtoMuTau_diTauMass.drawJobs = drawTemplateHistConfiguratorTauIdEffZtoMuTau_diTauMass.configure()

process.plotTemplateHistTauIdEffZtoMuTau = cms.Sequence(
    process.plotTemplateHistTauIdEffZtoMuTau_muonPt
   * process.plotTemplateHistTauIdEffZtoMuTau_muonAbsEta
   * process.plotTemplateHistTauIdEffZtoMuTau_diTauMass
)    

#--------------------------------------------------------------------------------
# fit number of events observed in tau id. discriminator passed/failed regions
#--------------------------------------------------------------------------------

process.fitTauIdEffZtoMuTau_tauIdPassed = cms.EDAnalyzer("TemplateHistFitter",                                          
    processes = cms.PSet(
        Ztautau = cms.PSet(
            templates = cms.PSet(
                ##muonPtVsAbsEta = cms.PSet(
                ##    meName = cms.string(dqmDirectory_Ztautau_tauIdPassed + '/' + meName_muonPtVsAbsEta),
                ##    fitSimultaneously = cms.bool(False)
                ##)
                diTauMass = cms.PSet(
                    meName = cms.string(dqmDirectory_Ztautau_tauIdPassed + '/' + meName_diTauMass),
                    fitSimultaneously = cms.bool(False)
                )                                                     
            ),    
            norm = cms.PSet(
                initial = cms.double(1000.)
            ),
            drawOptions = drawOption_Ztautau_separate                
        ),
        WplusJets = cms.PSet(
            templates = cms.PSet(
                ##muonPtVsAbsEta = cms.PSet(
                ##    #meName = cms.string(dqmDirectory_WplusJets_all + '/' + meName_muonPtVsAbsEta),
                ##    meName = cms.string(dqmDirectory_WplusJets_template + '/' + meName_muonPtVsAbsEta),
                ##    fitSimultaneously = cms.bool(False)
                ##)
                diTauMass = cms.PSet(
                    #meName = cms.string(dqmDirectory_WplusJets_all + '/' + meName_diTauMass),
                    meName = cms.string(dqmDirectory_WplusJets_template + '/' + meName_diTauMass),
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
                ##muonPtVsAbsEta = cms.PSet(
                ##    #meName = cms.string(dqmDirectory_QCD_all + '/' + meName_muonPtVsAbsEta),
                ##    meName = cms.string(dqmDirectory_QCD_template + '/' + meName_muonPtVsAbsEta),
                ##    fitSimultaneously = cms.bool(False)
                ##)
                diTauMass = cms.PSet(
                    #meName = cms.string(dqmDirectory_QCD_all + '/' + meName_diTauMass),
                    meName = cms.string(dqmDirectory_QCD_template + '/' + meName_diTauMass),
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
            ##muonPtVsAbsEta = cms.PSet(
            ##    meName = cms.string(dqmDirectory_Data_tauIdPassed + '/' + meName_muonPtVsAbsEta)
            ##)
            diTauMass = cms.PSet(
                meName = cms.string(dqmDirectory_Data_tauIdPassed + '/' + meName_diTauMass)
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
            ##muonPtVsAbsEta = cms.PSet(
            ##    name = cms.string("muonPtVsAbsEta"),
            ##    x = cms.PSet(
            ##        title = cms.string("|#eta_{#mu}|"),
            ##        min = cms.double(0.),
            ##        max = cms.double(2.1)
            ##    ),
            ##    y = cms.PSet(
            ##        title = cms.string("P_{T}^{#mu}"),
            ##        min = cms.double(15.),
            ##        max = cms.double(120.)
            ##    )
            ##)
            diTauMass = cms.PSet(
                name = cms.string("diTauMass"),
                title = cms.string("M^{#mu #tau}"),
                min = cms.double(0.),
                max = cms.double(250.)
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
            fileName = cms.string("./plots/fitTauIdEffZtoMuTau_tauIdPassed_#PLOT#.png")
        ),
        fitResults = cms.PSet(
            dqmDirectory = cms.string("fitTauIdEffZtoMuTau_tauIdPassed/fitResults")
        )
    )                                      
)                          

process.fitTauIdEffZtoMuTau_tauIdFailed = copy.deepcopy(process.fitTauIdEffZtoMuTau_tauIdPassed)
##process.fitTauIdEffZtoMuTau_tauIdFailed.processes.Ztautau.templates.muonPtVsAbsEta.meName = cms.string(dqmDirectory_Ztautau_tauIdFailed + '/' + meName_muonPtVsAbsEta)
process.fitTauIdEffZtoMuTau_tauIdFailed.processes.Ztautau.templates.diTauMass.meName = cms.string(dqmDirectory_Ztautau_tauIdFailed + '/' + meName_diTauMass)
process.fitTauIdEffZtoMuTau_tauIdFailed.processes.Ztautau.norm.initial = cms.double(500.)
#process.fitTauIdEffZtoMuTau_tauIdFailed.processes.WplusJets.templates.muonPtVsAbsEta.meName = cms.string(dqmDirectory_WplusJets_tauIdFailed + '/' + meName_muonPtVsAbsEta)
##process.fitTauIdEffZtoMuTau_tauIdFailed.processes.WplusJets.templates.muonPtVsAbsEta.meName = cms.string(dqmDirectory_WplusJets_template + '/' + meName_muonPtVsAbsEta)
process.fitTauIdEffZtoMuTau_tauIdFailed.processes.WplusJets.templates.diTauMass.meName = cms.string(dqmDirectory_WplusJets_template + '/' + meName_diTauMass)
process.fitTauIdEffZtoMuTau_tauIdFailed.processes.WplusJets.norm.initial = cms.double(2500.)
#process.fitTauIdEffZtoMuTau_tauIdFailed.processes.QCD.templates.muonPtVsAbsEta.meName = cms.string(dqmDirectory_QCD_tauIdFailed + '/' + meName_muonPtVsAbsEta)
##process.fitTauIdEffZtoMuTau_tauIdFailed.processes.QCD.templates.muonPtVsAbsEta.meName = cms.string(dqmDirectory_QCD_template + '/' + meName_muonPtVsAbsEta)
process.fitTauIdEffZtoMuTau_tauIdFailed.processes.QCD.templates.diTauMass.meName = cms.string(dqmDirectory_QCD_template + '/' + meName_diTauMass)
process.fitTauIdEffZtoMuTau_tauIdFailed.processes.QCD.norm.initial = cms.double(2500.)
##process.fitTauIdEffZtoMuTau_tauIdFailed.data.distributions.muonPtVsAbsEta.meName = cms.string(dqmDirectory_Data_tauIdFailed + '/' + meName_muonPtVsAbsEta)
process.fitTauIdEffZtoMuTau_tauIdFailed.data.distributions.diTauMass.meName = cms.string(dqmDirectory_Data_tauIdFailed + '/' + meName_diTauMass)
process.fitTauIdEffZtoMuTau_tauIdFailed.fit.algorithm.pluginName = cms.string("fitTauIdEffZtoMuTauAlgorithm_tauIdFailed")
process.fitTauIdEffZtoMuTau_tauIdFailed.output.controlPlots.fileName = cms.string("./plots/fitTauIdEffZtoMuTau_tauIdFailed_#PLOT#.png")
process.fitTauIdEffZtoMuTau_tauIdFailed.output.fitResults.dqmDirectory = cms.string("fitTauIdEffZtoMuTau_tauIdFailed/fitResults")

process.fitTauIdEffZtoMuTau = cms.Sequence(
    process.fitTauIdEffZtoMuTau_tauIdPassed
   * process.fitTauIdEffZtoMuTau_tauIdFailed
)

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


  