import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
#
# Compute reweighting parameters needed to correct for shape distortion
# of visible muon + tau-jet mass distribution in background enriched samples
# caused by cuts on transverse mass and CDF Pzeta variable
# in final event selection criteria of Z --> mu + tau-jet channel
# (in case no reweighting is done to correct for the distortion,
#  the normalization factors determined by the "template" method
#  are systematically wrong (in particular for the W + jets background))
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

process = cms.Process('prodKineEventReweightsZtoMuTau')

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

process.loadZtoMuTau = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(
        inputFileNames = cms.vstring('rfio:/castor/cern.ch/user/v/veelken/bgEstPlots/ZtoMuTau/plotsZtoMuTau_all_shrinkingCone.root'),
        scaleFactor = cms.double(1.),
        dqmDirectory_store = cms.string('')
    )
)

meName = 'DiTauCandidateQuantities/DPhi12'
meName_numerator = 'zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauPzetaDiff/' + meName
meName_denominator = 'zMuTauAnalyzer/afterEvtSelDiTauCandidateForMuTauAcoplanarity12_beforeEvtSelDiTauCandidateForMuTauMt1MET/' + meName

process.prodKineEventReweightsZtoMuTau = cms.EDAnalyzer("DQMHistEffProducer",
    plots = cms.PSet(
        diTauDPhi12Reweight_Ztautau = cms.PSet(
            numerator = cms.string('harvested/Ztautau/' + meName_numerator),
            denominator = cms.string('harvested/Ztautau/' + meName_denominator),
            efficiency = cms.string("bgEstKineEventReweights/Ztautau/diTauDPhi12")
        ),
        diTauDPhi12Reweight_WplusJets = cms.PSet(
            numerator = cms.string('harvested/WplusJets/' + meName_numerator),
            denominator = cms.string('harvested/WplusJets/' + meName_denominator),
            efficiency = cms.string("bgEstKineEventReweights/WplusJets/diTauDPhi12")
        ),
        diTauDPhi12Reweight_TTplusJets = cms.PSet(
            numerator = cms.string('harvested/TTplusJets/' + meName_numerator),
            denominator = cms.string('harvested/TTplusJets/' + meName_denominator),
            efficiency = cms.string("bgEstKineEventReweights/TTplusJets/diTauDPhi12")
        ),
        diTauDPhi12Reweight_QCD = cms.PSet(
            numerator = cms.string('harvested/qcdSum/' + meName_numerator),
            denominator = cms.string('harvested/qcdSum/' + meName_denominator),
            efficiency = cms.string("bgEstKineEventReweights/QCD/diTauDPhi12")
        )
    )                                                     
)


process.saveKineEventReweightsZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('bgEstKineEventReweightsZtoMuTau.root'),
    drop = cms.vstring('harvested')
)

process.p = cms.Path(
    process.loadZtoMuTau
   +process.prodKineEventReweightsZtoMuTau
   +process.saveKineEventReweightsZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
