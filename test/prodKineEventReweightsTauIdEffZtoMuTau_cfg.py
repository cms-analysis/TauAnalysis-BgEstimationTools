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

#--------------------------------------------------------------------------------
# produce histogram for kinematic reweighting
#--------------------------------------------------------------------------------

fileNames = dict()
fileNames["QCD"] = fileNamesZtoMuTau_qcdSum
fileNames["data"] = fileNamesZtoMuTau_pseudoData

bgEstEventSelection_QCDbeforeMuonIso = (
    "numDiTausQCDnoIso >= 1"
    " && tauDiscrAgainstMuonsQCDnoIso_0 > 0.5"
    " && numGlobalMuons < 2"
    " && diTauMt1MEtQCDnoIso_0 < 30."
    " && diTauAbsChargeQCDnoIso_0 > 0.5"
    " && metPt_0 < 20."
)

#--------------------------------------------------------------------------------
# CV: need to apply either muon track isolation or muon ECAL isolation cuts, but **not** both,
#     as too much QCD events get cut otherwise and the purity of the QCD background enriched sample drops
#     to 50% (from 80% in case only the muon track isolation cut is applied)
#--------------------------------------------------------------------------------

bgEstEventSelection_QCDafterMuonTrkIso = (
    "numDiTausQCDnoIso >= 1"
    " && muonTrackIsoQCDnoIso_0 < 1."
    " && tauDiscrAgainstMuonsQCDnoIso_0 > 0.5"
    " && numGlobalMuons < 2"
    " && diTauMt1MEtQCDnoIso_0 < 30."
    " && diTauAbsChargeQCDnoIso_0 > 0.5"
    " && metPt_0 < 20."
)

bgEstEventSelection_QCDafterMuonEcalIso = (
    "numDiTausQCDnoIso >= 1"
    " && muonEcalIsoQCDnoIso_0 < 1."
    " && tauDiscrAgainstMuonsQCDnoIso_0 > 0.5"
    " && numGlobalMuons < 2"
    " && diTauMt1MEtQCDnoIso_0 < 30."
    " && diTauAbsChargeQCDnoIso_0 > 0.5"
    " && metPt_0 < 20."
)

bgEstEventSelections = dict()
bgEstEventSelections["QCDbeforeMuonIso"] = bgEstEventSelection_QCDbeforeMuonIso
bgEstEventSelections["QCDafterMuonTrkIso"] = bgEstEventSelection_QCDafterMuonTrkIso
bgEstEventSelections["QCDafterMuonEcalIso"] = bgEstEventSelection_QCDafterMuonEcalIso

print("bgEstEventSelection_QCDbeforeMuonIso = " + bgEstEventSelections["QCDbeforeMuonIso"])
print("bgEstEventSelection_QCDafterMuonTrkIso = " + bgEstEventSelections["QCDafterMuonTrkIso"])
print("bgEstEventSelection_QCDafterMuonEcalIso = " + bgEstEventSelections["QCDafterMuonEcalIso"])

branchNames_muonPt = dict()
branchNames_muonPt["QCDbeforeMuonIso"] = "muonPtQCDnoIso_0"
branchNames_muonPt["QCDafterMuonTrkIso"] = branchNames_muonPt["QCDbeforeMuonIso"]
branchNames_muonPt["QCDafterMuonEcalIso"] = branchNames_muonPt["QCDbeforeMuonIso"]

branchNames_muonAbsEta = dict()
branchNames_muonAbsEta["QCDbeforeMuonIso"] = "muonAbsEtaQCDnoIso_0"
branchNames_muonAbsEta["QCDafterMuonTrkIso"] = branchNames_muonAbsEta["QCDbeforeMuonIso"]
branchNames_muonAbsEta["QCDafterMuonEcalIso"] = branchNames_muonAbsEta["QCDbeforeMuonIso"]

kineEventReweights = dict()
kineEventReweights["QCDbeforeMuonIso"] = None
kineEventReweights["QCDafterMuonTrkIso"] = None
kineEventReweights["QCDafterMuonEcalIso"] = None

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

process.prodKineEventReweightsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMHistEffProducer",
    config = cms.VPSet(
        cms.PSet(
            meName_numerator = cms.string('muonKineReweights/QCDafterMuonTrkIso/' + 'QCD' + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string('muonKineReweights/QCDbeforeMuonIso/' + 'QCD' + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string("muonKineReweights/QCDbgEnrichedTrkIso_pure/muonPtVsAbsEta")
        ),
        cms.PSet(
            meName_numerator = cms.string('muonKineReweights/QCDafterMuonTrkIso/' + 'data' + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string('muonKineReweights/QCDbeforeMuonIso/' + 'data' + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string("muonKineReweights/QCDbgEnrichedTrkIso_data/muonPtVsAbsEta")
        ),
        cms.PSet(
            meName_numerator = cms.string('muonKineReweights/QCDafterMuonEcalIso/' + 'QCD' + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string('muonKineReweights/QCDbeforeMuonIso/' + 'QCD' + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string("muonKineReweights/QCDbgEnrichedEcalIso_pure/muonPtVsAbsEta")
        ),
        cms.PSet(
            meName_numerator = cms.string('muonKineReweights/QCDafterMuonEcalIso/' + 'data' + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string('muonKineReweights/QCDbeforeMuonIso/' + 'data' + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string("muonKineReweights/QCDbgEnrichedEcalIso_data/muonPtVsAbsEta")
        ),
        cms.PSet(
            meNames_numerator = cms.vstring(
                 'muonKineReweights/QCDafterMuonTrkIso/' + 'QCD' + '/' + meName_muonPtVsAbsEta_norm,
                 'muonKineReweights/QCDafterMuonEcalIso/' + 'QCD' + '/' + meName_muonPtVsAbsEta_norm
            ),
            meNames_denominator = cms.vstring(
                 'muonKineReweights/QCDbeforeMuonIso/' + 'QCD' + '/' + meName_muonPtVsAbsEta_norm,
                 'muonKineReweights/QCDbeforeMuonIso/' + 'QCD' + '/' + meName_muonPtVsAbsEta_norm
            ),                 
            meName_efficiency = cms.string("muonKineReweights/QCDbgEnrichedCombIso_pure/muonPtVsAbsEta")
        ),
        cms.PSet(
            meNames_numerator = cms.vstring(
                 'muonKineReweights/QCDafterMuonTrkIso/' + 'data' + '/' + meName_muonPtVsAbsEta_norm,
                 'muonKineReweights/QCDafterMuonEcalIso/' + 'data' + '/' + meName_muonPtVsAbsEta_norm
            ),
            meNames_denominator = cms.vstring(
                 'muonKineReweights/QCDbeforeMuonIso/' + 'data' + '/' + meName_muonPtVsAbsEta_norm,
                 'muonKineReweights/QCDbeforeMuonIso/' + 'data' + '/' + meName_muonPtVsAbsEta_norm
            ),                 
            meName_efficiency = cms.string("muonKineReweights/QCDbgEnrichedCombIso_data/muonPtVsAbsEta")
        )
    )                                                     
)

process.saveKineEventReweightsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('muonKineReweightsTauIdEffZtoMuTau.root'),
    outputCommands = cms.vstring('drop harvested/*')
)

process.p = cms.Path(
    process.prodMuonKineReweightHistQCDenriched
  #+ process.dumpMuonKineReweightHistQCDenriched
   + process.prodKineEventReweightsTauIdEffZtoMuTau
   + process.saveKineEventReweightsTauIdEffZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
