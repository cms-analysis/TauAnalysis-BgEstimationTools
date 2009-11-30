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

process = cms.Process('prodKineEventReweightsTauIdEffZtoMuTau')

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

process.loadTauIdEffZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(
        inputFileNames = cms.vstring('/afs/cern.ch/user/v/veelken/scratch0/CMSSW_2_2_10/src/TauAnalysis/BgEstimationTools/test/plotsTauIdEffZtoMuTau_all.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

meName = 'MuonQuantities/MuonPtVsAbsEta'
meName_numerator = 'TauIdEffAnalyzerZtoMuTau/afterEvtSelMuonTrkIso_beforeEvtSelMuonEcalIso/' + meName
meName_denominator = 'TauIdEffAnalyzerZtoMuTau/afterEvtSelTauPt_beforeEvtSelMuonTrkIso/' + meName

process.prodKineEventReweightsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMHistEffProducer",
    config = cms.VPSet(
        cms.PSet(
            meName_numerator = cms.string('harvested/Ztautau/' + meName_numerator),
            meName_denominator = cms.string('harvested/Ztautau/' + meName_denominator),
            meName_efficiency = cms.string("bgEstKineEventReweights/Ztautau/muonPtVsAbsEta")
        ),
        cms.PSet(
            meName_numerator = cms.string('harvested/WplusJets/' + meName_numerator),
            meName_denominator = cms.string('harvested/WplusJets/' + meName_denominator),
            meName_efficiency = cms.string("bgEstKineEventReweights/WplusJets/muonPtVsAbsEta")
        ),
        cms.PSet(
            meName_numerator = cms.string('harvested/TTplusJets/' + meName_numerator),
            meName_denominator = cms.string('harvested/TTplusJets/' + meName_denominator),
            meName_efficiency = cms.string("bgEstKineEventReweights/TTplusJets/muonPtVsAbsEta")
        ),
        cms.PSet(
            meName_numerator = cms.string('harvested/qcdSum/' + meName_numerator),
            meName_denominator = cms.string('harvested/qcdSum/' + meName_denominator),
            meName_efficiency = cms.string("bgEstKineEventReweights/QCD/muonPtVsAbsEta")
        )
    )                                                     
)


process.saveKineEventReweightsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('bgEstKineEventReweightsTauIdEffZtoMuTau.root'),
    drop = cms.vstring('harvested')
)

process.p = cms.Path(
    process.loadTauIdEffZtoMuTau
   +process.prodKineEventReweightsTauIdEffZtoMuTau
   +process.saveKineEventReweightsTauIdEffZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
