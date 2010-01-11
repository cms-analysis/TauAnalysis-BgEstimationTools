import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
#
# Compute reweighting parameters needed to correct for shape distortion
# of muon |eta| distribution in QCD background events
# in event sample from which tau iso. && id. efficiencies are obtained,
# caused by cuts on muon track (and ECAL) isolation variables
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

from TauAnalysis.BgEstimationTools.bgEstNtupleDefinitionsZtoMuTau_cfi import *
from TauAnalysis.BgEstimationTools.templateHistDefinitions_cfi import *
from TauAnalysis.BgEstimationTools.tools.prodTemplateHistConfigurator import makeTemplateHistProdSequence2d
from TauAnalysis.BgEstimationTools.tools.prodTemplateHistConfigurator import makeTemplateHistProdSequence1d

process = cms.Process('prodKineEventReweightsTauIdEffZtoMuTau')

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

##process.loadTauIdEffZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
##    all = cms.PSet(
##        inputFileNames = cms.vstring('/afs/cern.ch/user/v/veelken/scratch0/CMSSW_2_2_10/src/TauAnalysis/BgEstimationTools/test/plotsTauIdEffZtoMuTau_all.root'),
##        scaleFactor = cms.double(1.),
##        dqmDirectory_store = cms.string('')
##    )
##)

#--------------------------------------------------------------------------------
# produce histogram for kinematic reweighting
#--------------------------------------------------------------------------------

fileNames = dict()
fileNames["QCD"] = fileNamesZtoMuTau_qcdSum
fileNames["data"] = fileNamesZtoMuTau_pseudoData

bgEstEventSelection_QCDbeforeMuonTrkIso = (
    "numDiTausQCDnoMuonIso >= 1"
    " && tauDiscrAgainstMuonsQCDnoMuonIso_0 > 0.5"
    " && numGlobalMuons < 2"
    " && diTauMt1MEtQCDnoMuonIso_0 < 30."
    " && diTauAbsChargeQCDnoMuonIso_0 > 0.5"
    " && metPt_0 < 20."
)

bgEstEventSelection_QCDafterMuonTrkIso = (
    "numDiTausQCDnoMuonIso >= 1 && muonTrackIsoQCDnoMuonIso_0 < 1."
    " && tauDiscrAgainstMuonsQCDnoMuonIso_0 > 0.5"
    " && numGlobalMuons < 2"
    " && diTauMt1MEtQCDnoMuonIso_0 < 30."
    " && diTauAbsChargeQCDnoMuonIso_0 > 0.5"
    " && metPt_0 < 20."
)

bgEstEventSelections = dict()
bgEstEventSelections["QCDbeforeMuonTrkIso"] = bgEstEventSelection_QCDbeforeMuonTrkIso
bgEstEventSelections["QCDafterMuonTrkIso"] = bgEstEventSelection_QCDafterMuonTrkIso

print("bgEstEventSelection_QCDbeforeMuonTrkIso = " + bgEstEventSelections["QCDbeforeMuonTrkIso"])
print("bgEstEventSelection_QCDafterMuonTrkIso = " + bgEstEventSelections["QCDafterMuonTrkIso"])

branchNames_muonPt = dict()
branchNames_muonPt["QCDbeforeMuonTrkIso"] = "muonPtQCD_0"
branchNames_muonPt["QCDafterMuonTrkIso"] = branchNames_muonPt["QCDbeforeMuonTrkIso"]

branchNames_muonAbsEta = dict()
branchNames_muonAbsEta["QCDbeforeMuonTrkIso"] = "muonAbsEtaQCD_0"
branchNames_muonAbsEta["QCDafterMuonTrkIso"] = branchNames_muonAbsEta["QCDbeforeMuonTrkIso"]

kineEventReweights = dict()
kineEventReweights["QCDbeforeMuonTrkIso"] = None
kineEventReweights["QCDafterMuonTrkIso"] = None

binEdges_muonPt = [ 15., 20., 25., 30., 40., 60., 120. ]
binEdges_muonAbsEta = [ 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1 ]

meName_muonPtVsAbsEta_norm = 'MuonPtVsAbsEtaShape'
meName_muonAbsEta_norm = 'MuonEtaShape'

process.prodMuonKineReweightHistQCDenrichedPtVsAbsEta = makeTemplateHistProdSequence2d(
    process, prodTemplateHist, fileNames, bgEstEventSelections, kineEventReweights,
    dqmDirectory = "muonKineReweights", meName = meName_muonPtVsAbsEta_norm,
    branchNamesX = branchNames_muonAbsEta, numBinsX = 7, binEdgesX = binEdges_muonAbsEta,
    branchNamesY = branchNames_muonPt, numBinsY = 6, binEdgesY = binEdges_muonPt
)

process.prodMuonKineReweightHistQCDenrichedAbsEta = makeTemplateHistProdSequence1d(
    process, prodTemplateHist, fileNames, bgEstEventSelections, kineEventReweights,
    dqmDirectory = "muonKineReweights", meName = meName_muonAbsEta_norm,
    branchNames = branchNames_muonAbsEta, numBins = 60, min = -3., max = +3.
)

process.prodMuonKineReweightHistQCDenriched = cms.Sequence(
    process.prodMuonKineReweightHistQCDenrichedPtVsAbsEta
   + process.prodMuonKineReweightHistQCDenrichedAbsEta
)

process.dumpMuonKineReweightHistQCDenriched = cms.EDAnalyzer("DQMStoreDump")

meName = 'MuonQuantities/MuonPtVsAbsEta'
meName_numerator = 'TauIdEffAnalyzerZtoMuTau/afterEvtSelMuonTrkIso_beforeEvtSelMuonEcalIso/' + meName
meName_denominator = 'TauIdEffAnalyzerZtoMuTau/afterEvtSelTauPt_beforeEvtSelMuonTrkIso/' + meName

process.prodKineEventReweightsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMHistEffProducer",
    config = cms.VPSet(
        ##cms.PSet(
        ##    meName_numerator = cms.string('harvested/qcdSum/' + meName_numerator),
        ##    meName_denominator = cms.string('harvested/qcdSum/' + meName_denominator),
        ##    meName_efficiency = cms.string("muonKineReweights/QCDmc/muonPtVsAbsEta")
        ##),
        cms.PSet(
            meName_numerator = cms.string('muonKineReweights/QCDbeforeMuonTrkIso/' + 'QCD' + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string('muonKineReweights/QCDafterMuonTrkIso/' + 'QCD' + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string("muonKineReweights/QCDbgEnriched_pure/muonPtVsAbsEta")
        ),
        cms.PSet(
            meName_numerator = cms.string('muonKineReweights/QCDbeforeMuonTrkIso/' + 'data' + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string('muonKineReweights/QCDafterMuonTrkIso/' + 'data' + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string("muonKineReweights/QCDbgEnriched_data/muonPtVsAbsEta")
        ),
    )                                                     
)

process.saveKineEventReweightsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('muonKineReweightsTauIdEffZtoMuTau.root'),
    outputCommands = cms.vstring('drop harvested/*')
)

process.p = cms.Path(
   #process.loadTauIdEffZtoMuTau
    process.prodMuonKineReweightHistQCDenriched
  #+ process.dumpMuonKineReweightHistQCDenriched
   + process.prodKineEventReweightsTauIdEffZtoMuTau
   + process.saveKineEventReweightsTauIdEffZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
