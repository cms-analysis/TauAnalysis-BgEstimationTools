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

from TauAnalysis.BgEstimationTools.bgEstNtupleDefinitionsZtoMuTau_cfi import *
from TauAnalysis.DQMTools.plotterStyleDefinitions_cfi import *
from TauAnalysis.BgEstimationTools.templateHistDefinitions_cfi import *
from TauAnalysis.BgEstimationTools.tools.prodTemplateHistConfigurator import makeTemplateHistProdSequence
from TauAnalysis.BgEstimationTools.tools.drawTemplateHistConfigurator import drawTemplateHistConfigurator
from TauAnalysis.BgEstimationTools.bgEstTemplateEvtSelZtoMuTau_cfi import *

processName = 'fitTauIdEffZtoMuTau'
process = cms.Process(processName)

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

process.loadTauIdEffZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(
        inputFileNames = cms.vstring('plotsTauIdEffZtoMuTau_all.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

#--------------------------------------------------------------------------------
# define directories in which histograms are stored in DQMStore
#--------------------------------------------------------------------------------

region_tauIdPassed = "region02"
region_tauIdFailed = "region01"

dqmDirectory_Ztautau_all = "harvested/Ztautau/TauIdEffAnalyzerZtoMuTau/afterDiTauCandidateBackToBackCut/"
dqmDirectory_Ztautau_tauIdPassed = dqmDirectory_Ztautau_all + "tauIdEffHistograms2regions/" + region_tauIdPassed + "/"
dqmDirectory_Ztautau_tauIdFailed = dqmDirectory_Ztautau_all + "tauIdEffHistograms2regions/" + region_tauIdFailed + "/"

dqmDirectory_WplusJets_all = "harvested/WplusJets/TauIdEffAnalyzerZtoMuTau/afterDiTauCandidateBackToBackCut/"
dqmDirectory_WplusJets_tauIdPassed = dqmDirectory_WplusJets_all + "tauIdEffHistograms2regions/" + region_tauIdPassed + "/"
dqmDirectory_WplusJets_tauIdFailed = dqmDirectory_WplusJets_all + "tauIdEffHistograms2regions/" + region_tauIdFailed + "/"

dqmDirectory_QCD_all = "harvested/qcdSum/TauIdEffAnalyzerZtoMuTau/afterDiTauCandidateBackToBackCut/"
dqmDirectory_QCD_tauIdPassed = dqmDirectory_QCD_all + "tauIdEffHistograms2regions/" + region_tauIdPassed + "/"
dqmDirectory_QCD_tauIdFailed = dqmDirectory_QCD_all + "tauIdEffHistograms2regions/" + region_tauIdFailed + "/"

dqmDirectory_Data_all = "harvested/smSum/TauIdEffAnalyzerZtoMuTau/afterDiTauCandidateBackToBackCut/"
dqmDirectory_Data_tauIdPassed = dqmDirectory_Data_all + "tauIdEffHistograms2regions/" + region_tauIdPassed + "/"
dqmDirectory_Data_tauIdFailed = dqmDirectory_Data_all + "tauIdEffHistograms2regions/" + region_tauIdFailed + "/"

#--------------------------------------------------------------------------------
# define names of histograms used in fit
#--------------------------------------------------------------------------------

dqmSubDirectory_muonPt  = 'TauIdEffSpecificQuantities/'
meName_muonPt = 'MuonPt'
dqmSubDirectory_muonAbsEta  = 'TauIdEffSpecificQuantities/'
meName_muonAbsEta = 'MuonAbsEta'
dqmSubDirectory_muonExtTrackIso = 'TauIdEffSpecificQuantities/'
meName_muonExtTrackIso = 'MuonExtTrkIsoPt'
meName_muonExtTrackIso_rebinned = 'MuonExtTrkIsoPt_rebinned'

#--------------------------------------------------------------------------------
# load template histograms
#--------------------------------------------------------------------------------

process.loadTemplateHistTauIdEffZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(
        inputFileNames = cms.vstring('fitTauIdEffZtoMuTau.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

#--------------------------------------------------------------------------------
# define names of branches in Ntuple from which template histograms get produced
#--------------------------------------------------------------------------------

branchNames_muonPt = dict()
branchNames_muonPt["Zmumu"] = "muonPtZmumu_0"
branchNames_muonPt["WplusJets"] = "muonPtWplusJets_0"
branchNames_muonPt["TTplusJets"] = "muonPtTTplusJets_0"
branchNames_muonPt["QCD"] = "muonPtQCD_0"

branchNames_muonAbsEta = dict()
branchNames_muonAbsEta["Zmumu"] = "muonAbsEtaZmumu_0"
branchNames_muonAbsEta["WplusJets"] = "muonAbsEtaWplusJets_0"
branchNames_muonAbsEta["TTplusJets"] = "muonAbsEtaTTplusJets_0"
branchNames_muonAbsEta["QCD"] = "muonAbsEtaQCD_0"

branchNames_muonExtTrackIso = dict()
branchNames_muonExtTrackIso["Zmumu"] = "muonExtTrackIsoZmumu_0"
branchNames_muonExtTrackIso["WplusJets"] = "muonExtTrackIsoWplusJets_0"
branchNames_muonExtTrackIso["TTplusJets"] = "muonExtTrackIsoTTplusJets_0"
branchNames_muonExtTrackIso["QCD"] = "muonExtTrackIsoQCD_0"

kineEventReweights = dict()
kineEventReweights["Zmumu"] = None
kineEventReweights["WplusJets"] = None
kineEventReweights["TTplusJets"] = None
kineEventReweights["QCD"] = None

#--------------------------------------------------------------------------------
# produce template histograms 
#--------------------------------------------------------------------------------

fileNames = dict()
fileNames["Ztautau"] = fileNamesZtoMuTau_Ztautau
fileNames["Zmumu"] = fileNamesZtoMuTau_ZmumuPlusJets
fileNames["WplusJets"] = fileNamesZtoMuTau_WplusJets
fileNames["TTplusJets"] = fileNamesZtoMuTau_TTplusJets
fileNames["QCD"] = fileNamesZtoMuTau_qcdSum
fileNames["data"] = fileNamesZtoMuTau_pseudoData

bgEstEventSelections = dict()
bgEstEventSelections["Zmumu"] = bgEstEventSelection_Zmumu
bgEstEventSelections["WplusJets"] = bgEstEventSelection_WplusJets 
bgEstEventSelections["TTplusJets"] = bgEstEventSelection_TTplusJets
bgEstEventSelections["QCD"] = bgEstEventSelection_QCD

print("bgEstEventSelection_Zmumu = " + bgEstEventSelections["Zmumu"])
print("bgEstEventSelection_WplusJets = " + bgEstEventSelections["WplusJets"])
print("bgEstEventSelection_TTplusJets = " + bgEstEventSelections["TTplusJets"])
print("bgEstEventSelection_QCD = " + bgEstEventSelections["QCD"])

process.prodTemplateHistTauIdEffZtoMuTau_muonPt = makeTemplateHistProdSequence(
    process, prodTemplateHist, fileNames, bgEstEventSelections, branchNames_muonPt, kineEventReweights,
    dqmDirectory = processName, meName = meName_muonPt, numBinsX = 9, xBins = [ 15., 20., 25., 30., 35., 40., 50., 60., 80., 120. ]
)

process.prodTemplateHistTauIdEffZtoMuTau_muonAbsEta = makeTemplateHistProdSequence(
    process, prodTemplateHist, fileNames, bgEstEventSelections, branchNames_muonAbsEta, kineEventReweights,
    dqmDirectory = processName, meName = meName_muonAbsEta, numBinsX = 14, xMin = 0., xMax = 2.1
)

process.prodTemplateHistTauIdEffZtoMuTau_muonExtTrackIso = makeTemplateHistProdSequence(
    process, prodTemplateHist, fileNames, bgEstEventSelections, branchNames_muonExtTrackIso, kineEventReweights,
    dqmDirectory = processName, meName = meName_muonExtTrackIso_rebinned, numBinsX = 8, xBins = [ 0., 0.5, 1., 1.5, 2., 3., 4., 6., 8. ]
)

process.prodTemplateHistTauIdEffZtoMuTau = cms.Sequence(
    process.prodTemplateHistTauIdEffZtoMuTau_muonPt
   * process.prodTemplateHistTauIdEffZtoMuTau_muonAbsEta
   * process.prodTemplateHistTauIdEffZtoMuTau_muonExtTrackIso
)    

#--------------------------------------------------------------------------------
# rebin template histogram for muon track isolation in (extended) dR = 1.0 cone
#--------------------------------------------------------------------------------

process.rebinTemplateHistTauIdEffZtoMuTau_muonExtTrackIso = cms.EDAnalyzer("DQMHistRebinner",
    config = cms.VPSet(
        cms.PSet(
            meName_original = cms.string(dqmDirectory_Ztautau_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_Ztautau_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_Ztautau_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_Ztautau_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_Ztautau_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_Ztautau_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        ##cms.PSet(
        ##    meName_original = cms.string(dqmDirectory_Zmumu_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
        ##    meName_rebinned = cms.string(dqmDirectory_Zmumu_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ##),
        ##cms.PSet(
        ##    meName_original = cms.string(dqmDirectory_Zmumu_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
        ##    meName_rebinned = cms.string(dqmDirectory_Zmumu_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ##),
        ##cms.PSet(
        ##    meName_original = cms.string(dqmDirectory_Zmumu_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
        ##    meName_rebinned = cms.string(dqmDirectory_Zmumu_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ##),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_WplusJets_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_WplusJets_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_WplusJets_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_WplusJets_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_WplusJets_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_WplusJets_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        ##cms.PSet(
        ##    meName_original = cms.string(dqmDirectory_TTplusJets_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
        ##    meName_rebinned = cms.string(dqmDirectory_TTplusJets_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ##),
        ##cms.PSet(
        ##    meName_original = cms.string(dqmDirectory_TTplusJets_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
        ##    meName_rebinned = cms.string(dqmDirectory_TTplusJets_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ##),
        ##cms.PSet(
        ##    meName_original = cms.string(dqmDirectory_TTplusJets_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
        ##    meName_rebinned = cms.string(dqmDirectory_TTplusJets_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ##),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_QCD_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_QCD_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_QCD_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_QCD_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_QCD_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_QCD_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_Data_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_Data_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_Data_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_Data_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        ),
        cms.PSet(
            meName_original = cms.string(dqmDirectory_Data_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso),
            meName_rebinned = cms.string(dqmDirectory_Data_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
        )
    ),
    numBinsX = cms.uint32(8),
    xBins = cms.vdouble( 0., 0.5, 1., 1.5, 2., 3., 4., 6., 8. )
)                                                                           

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

drawJobTemplateHist_muonExtTrackIso = copy.deepcopy(drawJobTemplateHist_muonPt)
drawJobTemplateHist_muonExtTrackIso.title = cms.string('Muon Track iso. P_{T}')
drawJobTemplateHist_muonExtTrackIso.xAxis = cms.string('Pt')

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonExtTrackIso = drawTemplateHistConfigurator(
    template = drawJobTemplateHist_muonExtTrackIso
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
        dqmDirectory_Ztautau_tauIdPassed + dqmSubDirectory_muonPt + meName_muonPt,
        dqmDirectory_Ztautau_tauIdFailed + dqmSubDirectory_muonPt + meName_muonPt,
        dqmDirectory_Ztautau_all + dqmSubDirectory_muonPt + meName_muonPt
    ],
    name = "Ztautau_muonPt",
    title = "P_{T}^{#mu} in Z #rightarrow #tau^{+} #tau^{-} Signal" 
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.add(
    meNames = [
        dqmDirectory_Ztautau_tauIdPassed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta,
        dqmDirectory_Ztautau_tauIdFailed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta,
        dqmDirectory_Ztautau_all + dqmSubDirectory_muonAbsEta + meName_muonAbsEta
    ],
    name = "Ztautau_muonAbsEta",
    title = "|#eta_{#mu}| in Z #rightarrow #tau^{+} #tau^{-} Signal" 
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonExtTrackIso.add(
    meNames = [
        dqmDirectory_Ztautau_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned,
        dqmDirectory_Ztautau_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned,
        dqmDirectory_Ztautau_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned
    ],
    name = "Ztautau_muonExtTrackIso",
    title = "|Track iso. P_{T}| in Z #rightarrow #tau^{+} #tau^{-} Signal" 
)

#--------------------------------------------------------------------------------
# define draw jobs for Z --> mu+ mu- background
#--------------------------------------------------------------------------------

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.add(
##    meNames = [
##        dqmDirectory_Zmumu_tauIdPassed + dqmSubDirectory_muonPt + meName_muonPt,
##        dqmDirectory_Zmumu_tauIdFailed + dqmSubDirectory_muonPt + meName_muonPt,
##        dqmDirectory_Zmumu_all + dqmSubDirectory_muonPt + meName_muonPt
##        #processName + '/Zmumu/data/' + meName_muonPt
##    ],
##    name = "Zmumu_muonPt",
##    title = "P_{T}^{#mu} in Z #rightarrow #mu^{+} #mu^{-} Background" 
##)

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.add(
##    meNames = [
##        dqmDirectory_Zmumu_tauIdPassed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta,
##        dqmDirectory_Zmumu_tauIdFailed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta,
##        dqmDirectory_Zmumu_all + dqmSubDirectory_muonAbsEta + meName_muonAbsEta
##        #processName + '/Zmumu/data/' + meName_muonAbsEta
##    ],
##    name = "Zmumu_muonAbsEta",
##    title = "|#eta_{#mu}| in Z #rightarrow #mu^{+} #mu^{-} Background"
##)

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonExtTrackIso.add(
##    meNames = [
##        dqmDirectory_Zmumu_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned,
##        dqmDirectory_Zmumu_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned,
##        dqmDirectory_Zmumu_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned
##        #processName + '/Zmumu/data/' + meName_muonExtTrackIso_rebinned
##    ],
##    name = "Zmumu_muonExtTrackIso",
##    title = "|Track iso. P_{T}| in Z #rightarrow #mu^{+} #mu^{-} Background"
##)

#--------------------------------------------------------------------------------
# define draw jobs for W + jets background
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.add(
    meNames = [
        dqmDirectory_WplusJets_tauIdPassed + dqmSubDirectory_muonPt + meName_muonPt,
        dqmDirectory_WplusJets_tauIdFailed + dqmSubDirectory_muonPt + meName_muonPt,
        dqmDirectory_WplusJets_all + dqmSubDirectory_muonPt + meName_muonPt
        #processName + '/WplusJets/data/' + meName_muonPt
        
    ],
    name = "WplusJets_muonPt",
    title = "P_{T}^{#mu} in W + jets Background"
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.add(
    meNames = [
        dqmDirectory_WplusJets_tauIdPassed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta,
        dqmDirectory_WplusJets_tauIdFailed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta,
        dqmDirectory_WplusJets_all + dqmSubDirectory_muonAbsEta + meName_muonAbsEta
        #processName + '/WplusJets/data/' + meName_muonAbsEta
    ],
    name = "WplusJets_muonAbsEta",
    title = "|#eta_{#mu}| in W + jets Background"
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonExtTrackIso.add(
    meNames = [
        dqmDirectory_WplusJets_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned,
        dqmDirectory_WplusJets_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned,
        dqmDirectory_WplusJets_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned
        #processName + '/WplusJets/data/' + meName_muonExtTrackIso_rebinned
    ],
    name = "WplusJets_muonExtTrackIso",
    title = "|Track iso. P_{T}| in W + jets Background"
)

#--------------------------------------------------------------------------------
# define draw jobs for TTbar + jets background
#--------------------------------------------------------------------------------

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.add(
##    meNames = [
##        dqmDirectory_TTplusJets_tauIdPassed + dqmSubDirectory_muonPt + meName_muonPt,
##        dqmDirectory_TTplusJets_tauIdFailed + dqmSubDirectory_muonPt + meName_muonPt,
##        dqmDirectory_TTplusJets_all + dqmSubDirectory_muonPt + meName_muonPt
##        #processName + '/TTplusJets/data/' + meName_muonPt
##    ],
##    name = "TTplusJets_muonPt",
##    title = "P_{T}^{#mu} in t#bar{t} + jets Background"
##)

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.add(
##    meNames = [
##        dqmDirectory_TTplusJets_tauIdPassed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta,
##        dqmDirectory_TTplusJets_tauIdFailed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta,
##        dqmDirectory_TTplusJets_all + dqmSubDirectory_muonAbsEta + meName_muonAbsEta
##        #processName + '/TTplusJets/data/' + meName_muonAbsEta
##    ],
##    name = "TTplusJets_muonAbsEta",
##    title = "|#eta_{#mu}| in t#bar{t} + jets Background"
##)

##drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonExtTrackIso.add(
##    meNames = [
##        dqmDirectory_TTplusJets_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned,
##        dqmDirectory_TTplusJets_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned,
##        dqmDirectory_TTplusJets_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned
##        #processName + '/TTplusJets/data/' + meName_muonExtTrackIso_rebinned
##    ],
##    name = "TTplusJets_muonExtTrackIso",
##    title = "|Track iso. P_{T}| in t#bar{t} + jets Background"
##)

#--------------------------------------------------------------------------------
# define draw jobs for QCD background
#--------------------------------------------------------------------------------

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.add(
    meNames = [
        dqmDirectory_QCD_tauIdPassed + dqmSubDirectory_muonPt + meName_muonPt,
        dqmDirectory_QCD_tauIdFailed + dqmSubDirectory_muonPt + meName_muonPt,
        dqmDirectory_QCD_all + dqmSubDirectory_muonPt + meName_muonPt
        #processName + '/QCD/data/' + meName_muonPt
    ],
    name = "QCD_muonPt",
    title = "P_{T}^{#mu} in QCD Background"
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.add(
    meNames = [
        dqmDirectory_QCD_tauIdPassed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta,
        dqmDirectory_QCD_tauIdFailed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta,
        dqmDirectory_QCD_all + dqmSubDirectory_muonAbsEta + meName_muonAbsEta
        #processName + '/QCD/data/' + meName_muonAbsEta
    ],
    name = "QCD_muonAbsEta",
    title = "|#eta_{#mu}| in QCD Background"
)

drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonExtTrackIso.add(
    meNames = [
        dqmDirectory_QCD_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned,
        dqmDirectory_QCD_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned,
        dqmDirectory_QCD_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned
        #processName + '/QCD/data/' + meName_muonExtTrackIso_rebinned
    ],
    name = "QCD_muonExtTrackIso",
    title = "|Track iso. P_{T}| in QCD Background"
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
        Eta = copy.deepcopy(xAxis_eta)
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
    indOutputFileName = cms.string('plotTemplateHistTauIdEffZtoMuTau_#PLOT#.eps')
)

process.plotTemplateHistTauIdEffZtoMuTau_muonPt = copy.deepcopy(plotTemplateHistZtoMuTau)
process.plotTemplateHistTauIdEffZtoMuTau_muonPt.drawJobs = drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonPt.configure()

process.plotTemplateHistTauIdEffZtoMuTau_muonAbsEta = copy.deepcopy(plotTemplateHistZtoMuTau)
process.plotTemplateHistTauIdEffZtoMuTau_muonAbsEta.drawJobs = drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonAbsEta.configure()

process.plotTemplateHistTauIdEffZtoMuTau_muonExtTrackIso = copy.deepcopy(plotTemplateHistZtoMuTau)
process.plotTemplateHistTauIdEffZtoMuTau_muonExtTrackIso.drawJobs = drawTemplateHistConfiguratorTauIdEffZtoMuTau_muonExtTrackIso.configure()

process.plotTemplateHistTauIdEffZtoMuTau = cms.Sequence(
    process.plotTemplateHistTauIdEffZtoMuTau_muonPt
   * process.plotTemplateHistTauIdEffZtoMuTau_muonAbsEta
   * process.plotTemplateHistTauIdEffZtoMuTau_muonExtTrackIso
)    

#--------------------------------------------------------------------------------
# fit number of events observed in tau id. discriminator passed/failed regions
#--------------------------------------------------------------------------------

process.fitTauIdEffZtoMuTau_tauIdPassed = cms.EDAnalyzer("TemplateBgEstFit",                                          
    processes = cms.PSet(
        Ztautau = cms.PSet(
            templates = cms.PSet(
                muonPt = cms.PSet(
                    meName = cms.string(dqmDirectory_Ztautau_all + dqmSubDirectory_muonPt + meName_muonPt),
                    fitSimultaneously = cms.bool(False)
                ),
                muonAbsEta = cms.PSet(
                    meName = cms.string(dqmDirectory_Ztautau_all + dqmSubDirectory_muonAbsEta + meName_muonAbsEta),
                    fitSimultaneously = cms.bool(False)
                ),
                muonExtTrackIso = cms.PSet(
                    meName = cms.string(dqmDirectory_Ztautau_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned),
                    fitSimultaneously = cms.bool(False)
                )
            ),    
            norm = cms.PSet(
                initial = cms.double(1000.)
            ),
            drawOptions = drawOption_Ztautau                
        ),
        WplusJets = cms.PSet(
            templates = cms.PSet(
                muonPt = cms.PSet(
                    meName = cms.string(dqmDirectory_WplusJets_all + dqmSubDirectory_muonPt + meName_muonPt),
                    fitSimultaneously = cms.bool(False)
                ),
                muonAbsEta = cms.PSet(
                    meName = cms.string(dqmDirectory_WplusJets_all + dqmSubDirectory_muonAbsEta + meName_muonAbsEta),
                    fitSimultaneously = cms.bool(False)
                ),
                muonExtTrackIso = cms.PSet(
                    meName = cms.string(dqmDirectory_WplusJets_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned),
                    fitSimultaneously = cms.bool(False)
                )
            ),    
            norm = cms.PSet(
                initial = cms.double(100.)
            ),
            drawOptions = drawOption_WplusJets
        ),
        QCD = cms.PSet(
            templates = cms.PSet(
                muonPt = cms.PSet(
                    meName = cms.string(dqmDirectory_QCD_all + dqmSubDirectory_muonPt + meName_muonPt),
                    fitSimultaneously = cms.bool(False)
                ),
                muonAbsEta = cms.PSet(
                    meName = cms.string(dqmDirectory_QCD_all + dqmSubDirectory_muonAbsEta + meName_muonAbsEta),
                    fitSimultaneously = cms.bool(False)
                ),
                muonExtTrackIso = cms.PSet(
                    meName = cms.string(dqmDirectory_QCD_all + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned),
                    fitSimultaneously = cms.bool(False)
                )
            ),    
            norm = cms.PSet(
                initial = cms.double(100.)
            ),
            drawOptions = drawOption_QCD
        )
    ),

    # use "pseudo" data-samples consisting of all Monte Carlo processes for testing                      
    data = cms.PSet(
        distributions = cms.PSet(
            muonPt = cms.PSet(
                meName = cms.string(dqmDirectory_Data_tauIdPassed + dqmSubDirectory_muonPt + meName_muonPt)
            ),
            muonAbsEta = cms.PSet(
                meName = cms.string(dqmDirectory_Data_tauIdPassed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta)
            ),
            muonExtTrackIso = cms.PSet(
                meName = cms.string(dqmDirectory_Data_tauIdPassed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
            )
        )
    ),

    fit = cms.PSet(
        variables = cms.PSet(
            muonPt = cms.PSet(
               name = cms.string("muonPt"),
               title = cms.string("P_{T}^{#mu}"),
               xMin = cms.double(15.),
               xMax = cms.double(120.)
            ),
            muonAbsEta = cms.PSet(
               name = cms.string("muonAbsEta"),
               title = cms.string("|#eta_{#mu}|"),
               xMin = cms.double(0.),
               xMax = cms.double(2.1)
            ),
            muonExtTrackIso = cms.PSet(
               name = cms.string("muonExtTrackIso"),
               title = cms.string("Track iso. P_{T}"),
               xMin = cms.double(0.),
               xMax = cms.double(8.)
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
            fileName = cms.string("./plots/fitTauIdEffZtoMuTau_tauIdPassed_#PLOT#.eps")
        ),
        fitResults = cms.PSet(
            dqmDirectory = cms.string("fitTauIdEffZtoMuTau/fitResults/")
        )
    )                                      
)                          

process.fitTauIdEffZtoMuTau_tauIdFailed = copy.deepcopy(process.fitTauIdEffZtoMuTau_tauIdPassed)
process.fitTauIdEffZtoMuTau_tauIdFailed.processes.Ztautau.norm.initial = cms.double(500.)
process.fitTauIdEffZtoMuTau_tauIdFailed.processes.WplusJets.norm.initial = cms.double(2500.)
process.fitTauIdEffZtoMuTau_tauIdFailed.processes.QCD.norm.initial = cms.double(2500.)
process.fitTauIdEffZtoMuTau_tauIdFailed.data.distributions.muonPt.meName = cms.string(dqmDirectory_Data_tauIdFailed + dqmSubDirectory_muonPt + meName_muonPt)
process.fitTauIdEffZtoMuTau_tauIdFailed.data.distributions.muonAbsEta.meName = cms.string(dqmDirectory_Data_tauIdFailed + dqmSubDirectory_muonAbsEta + meName_muonAbsEta)
process.fitTauIdEffZtoMuTau_tauIdFailed.data.distributions.muonExtTrackIso.meName = cms.string(dqmDirectory_Data_tauIdFailed + dqmSubDirectory_muonExtTrackIso + meName_muonExtTrackIso_rebinned)
process.fitTauIdEffZtoMuTau_tauIdFailed.output.controlPlots.fileName = cms.string("./plots/fitTauIdEffZtoMuTau_tauIdFailed_#PLOT#.eps")

process.fitTauIdEffZtoMuTau = cms.Sequence(
    process.fitTauIdEffZtoMuTau_tauIdPassed
   * process.fitTauIdEffZtoMuTau_tauIdFailed
)

process.saveAllHistTauIdEffZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('fitTauIdEffZtoMuTau.root'),
    drop = cms.vstring('')                       
)

process.p = cms.Path(
    process.loadTauIdEffZtoMuTau
   #+ process.loadTemplateHistTauIdEffZtoMuTau
   + process.prodTemplateHistTauIdEffZtoMuTau
   + process.rebinTemplateHistTauIdEffZtoMuTau_muonExtTrackIso
   + process.saveAllHistTauIdEffZtoMuTau
   + process.plotTemplateHistTauIdEffZtoMuTau
   + process.fitTauIdEffZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
