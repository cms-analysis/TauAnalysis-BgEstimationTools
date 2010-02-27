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

#--------------------------------------------------------------------------------
# load histograms for kinematic reweighting
#--------------------------------------------------------------------------------

process.load("TauAnalysis.Configuration.plotZtoMuTau_cff")

process.loadZtoMuTau.inputFilePath = cms.string("rfio:/castor/cern.ch/user/v/veelken/CMSSW_3_3_x/plots/ZtoMuTau_tauIdEffKineEventReweights/7TeV/")

process.addTauIdEffKineReweightHistogramsZtoMuTau_qcdSum = cms.EDAnalyzer("DQMHistAdder",
    qcdSum = cms.PSet(
        dqmDirectories_input = cms.vstring(
            'harvested/InclusivePPmuX/',
            'harvested/PPmuXptGt20/'
        ),
        dqmDirectory_output = cms.string('harvested/qcdSum/')
    )                          
)

process.addTauIdEffKineReweightHistogramsZtoMuTau_smSum = cms.EDAnalyzer("DQMHistAdder",
    smSum = cms.PSet(
        dqmDirectories_input = cms.vstring(
            'harvested/Ztautau/',
            'harvested/Zmumu/',
            'harvested/WplusJets/',
            'harvested/TTplusJets/',
            'harvested/qcdSum/'
        ),
        dqmDirectory_output = cms.string('harvested/smSum/')
    )
)

process.addTauIdEffKineReweightHistogramsZtoMuTau = cms.Sequence(
    process.addTauIdEffKineReweightHistogramsZtoMuTau_qcdSum + process.addTauIdEffKineReweightHistogramsZtoMuTau_smSum
)

process.dumpTauIdEffKineReweightHistogramsZtoMuTau = cms.EDAnalyzer("DQMStoreDump")

dqmDirectory_qcdSumLooseMuonIso = 'harvested/qcdSum/TauIdEffQCDenrichedLooseMuonIso/TauIdEffSpecificQuantities'
dqmDirectory_qcdSumTightMuonTrkIso = 'harvested/qcdSum/TauIdEffQCDenrichedTightMuonTrkIso/TauIdEffSpecificQuantities'
dqmDirectory_qcdSumTightMuonEcalIso = 'harvested/qcdSum/TauIdEffQCDenrichedTightMuonEcalIso/TauIdEffSpecificQuantities'
dqmDirectory_qcdSumTightMuonCombIso = 'harvested/qcdSum/TauIdEffQCDenrichedTightMuonCombIso/TauIdEffSpecificQuantities'

dqmDirectory_smSumLooseMuonIso = 'harvested/smSum/TauIdEffQCDenrichedLooseMuonIso/TauIdEffSpecificQuantities'
dqmDirectory_smSumTightMuonTrkIso = 'harvested/smSum/TauIdEffQCDenrichedTightMuonTrkIso/TauIdEffSpecificQuantities'
dqmDirectory_smSumTightMuonEcalIso = 'harvested/smSum/TauIdEffQCDenrichedTightMuonEcalIso/TauIdEffSpecificQuantities'
dqmDirectory_smSumTightMuonCombIso = 'harvested/smSum/TauIdEffQCDenrichedTightMuonCombIso/TauIdEffSpecificQuantities'

meName_muonPtVsAbsEta = 'MuonPtVsAbsEta'
meName_muonPtVsAbsEta_norm = meName_muonPtVsAbsEta + '_norm'

process.normalizeTauIdEffKineReweightHistogramsZtoMuTau = cms.EDAnalyzer("DQMHistNormalizer",
    config = cms.VPSet(
        cms.PSet(
            meNameInput = cms.string(dqmDirectory_qcdSumLooseMuonIso + '/' + meName_muonPtVsAbsEta),
            meNameOutput = cms.string(dqmDirectory_qcdSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm)
        ),
        cms.PSet(
            meNameInput = cms.string(dqmDirectory_qcdSumTightMuonTrkIso + '/' + meName_muonPtVsAbsEta),
            meNameOutput = cms.string(dqmDirectory_qcdSumTightMuonTrkIso + '/' + meName_muonPtVsAbsEta_norm)
        ),
        cms.PSet(
            meNameInput = cms.string(dqmDirectory_qcdSumTightMuonEcalIso + '/' + meName_muonPtVsAbsEta),
            meNameOutput = cms.string(dqmDirectory_qcdSumTightMuonEcalIso + '/' + meName_muonPtVsAbsEta_norm)
        ),
        cms.PSet(
            meNameInput = cms.string(dqmDirectory_qcdSumTightMuonCombIso + '/' + meName_muonPtVsAbsEta),
            meNameOutput = cms.string(dqmDirectory_qcdSumTightMuonCombIso + '/' + meName_muonPtVsAbsEta_norm)
        ),
        cms.PSet(
            meNameInput = cms.string(dqmDirectory_smSumLooseMuonIso + '/' + meName_muonPtVsAbsEta),
            meNameOutput = cms.string(dqmDirectory_smSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm)
        ),
        cms.PSet(
            meNameInput = cms.string(dqmDirectory_smSumTightMuonTrkIso + '/' + meName_muonPtVsAbsEta),
            meNameOutput = cms.string(dqmDirectory_smSumTightMuonTrkIso + '/' + meName_muonPtVsAbsEta_norm)
        ),
        cms.PSet(
            meNameInput = cms.string(dqmDirectory_smSumTightMuonEcalIso + '/' + meName_muonPtVsAbsEta),
            meNameOutput = cms.string(dqmDirectory_smSumTightMuonEcalIso + '/' + meName_muonPtVsAbsEta_norm)
        ),
        cms.PSet(
            meNameInput = cms.string(dqmDirectory_smSumTightMuonCombIso + '/' + meName_muonPtVsAbsEta),
            meNameOutput = cms.string(dqmDirectory_smSumTightMuonCombIso + '/' + meName_muonPtVsAbsEta_norm)
        )
    ),
    norm = cms.double(1.)                                                                     
)

dqmDirectory_effTightMuonTrkIso = 'tauIdEffKineEventReweights/QCDenrichedMuonTrkIso'
dqmDirectory_effTightMuonEcalIso = 'tauIdEffKineEventReweights/QCDenrichedMuonEcalIso'
dqmDirectory_effTightMuonCombIsoFactorized = 'tauIdEffKineEventReweights/QCDenrichedMuonCombIsoFactorized'
dqmDirectory_effTightMuonCombIso = 'tauIdEffKineEventReweights/QCDenrichedMuonCombIso'

process.prodKineEventReweightsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMHistEffProducer",
    config = cms.VPSet(
        cms.PSet(
            meName_numerator = cms.string(dqmDirectory_qcdSumTightMuonTrkIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string(dqmDirectory_qcdSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string(dqmDirectory_effTightMuonTrkIso + '_pure' + '/' + meName_muonPtVsAbsEta)
        ),
        cms.PSet(
            meName_numerator = cms.string(dqmDirectory_qcdSumTightMuonEcalIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string(dqmDirectory_qcdSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string(dqmDirectory_effTightMuonEcalIso + '_pure' + '/' + meName_muonPtVsAbsEta)
        ),
        cms.PSet(
            meNames_numerator = cms.vstring(
                dqmDirectory_qcdSumTightMuonTrkIso + '/' + meName_muonPtVsAbsEta_norm,
                dqmDirectory_qcdSumTightMuonEcalIso + '/' + meName_muonPtVsAbsEta_norm
            ),
            meNames_denominator = cms.vstring(
                dqmDirectory_qcdSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm,
                dqmDirectory_qcdSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm
            ),
            meName_efficiency = cms.string(dqmDirectory_effTightMuonCombIsoFactorized + '_pure' + '/' + meName_muonPtVsAbsEta)
        ),
        cms.PSet(
            meName_numerator = cms.string(dqmDirectory_qcdSumTightMuonCombIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string(dqmDirectory_qcdSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string(dqmDirectory_effTightMuonCombIso + '_pure' + '/' + meName_muonPtVsAbsEta)
        ),
        cms.PSet(
            meName_numerator = cms.string(dqmDirectory_smSumTightMuonTrkIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string(dqmDirectory_smSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string(dqmDirectory_effTightMuonTrkIso + '_data' + '/' + meName_muonPtVsAbsEta)
        ),
        cms.PSet(
            meName_numerator = cms.string(dqmDirectory_smSumTightMuonEcalIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string(dqmDirectory_smSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string(dqmDirectory_effTightMuonEcalIso + '_data' + '/' + meName_muonPtVsAbsEta)
        ),
        cms.PSet(
            meNames_numerator = cms.vstring(
                dqmDirectory_smSumTightMuonTrkIso + '/' + meName_muonPtVsAbsEta_norm,
                dqmDirectory_smSumTightMuonEcalIso + '/' + meName_muonPtVsAbsEta_norm
            ),
            meNames_denominator = cms.vstring(
                dqmDirectory_smSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm,
                dqmDirectory_smSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm
            ),
            meName_efficiency = cms.string(dqmDirectory_effTightMuonCombIsoFactorized + '_data' + '/' + meName_muonPtVsAbsEta)
        ),
        cms.PSet(
            meName_numerator = cms.string(dqmDirectory_smSumTightMuonCombIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_denominator = cms.string(dqmDirectory_smSumLooseMuonIso + '/' + meName_muonPtVsAbsEta_norm),
            meName_efficiency = cms.string(dqmDirectory_effTightMuonCombIso + '_data' + '/' + meName_muonPtVsAbsEta)
        )
    )                                                     
)

process.dumpKineEventReweightsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMDumpHistogram",
    meNames = cms.vstring(
        dqmDirectory_effTightMuonCombIsoFactorized + '_pure' + '/' + meName_muonPtVsAbsEta,
        dqmDirectory_effTightMuonCombIso + '_pure' + '/' + meName_muonPtVsAbsEta,
        dqmDirectory_effTightMuonCombIsoFactorized + '_data' + '/' + meName_muonPtVsAbsEta,
        dqmDirectory_effTightMuonCombIso + '_data' + '/' + meName_muonPtVsAbsEta
    )
)

process.saveKineEventReweightsTauIdEffZtoMuTau = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('muonKineReweightsTauIdEffZtoMuTau.root'),
    outputCommands = cms.vstring('drop harvested/*')
)

process.p = cms.Path(
     process.loadZtoMuTau
   + process.addTauIdEffKineReweightHistogramsZtoMuTau
   + process.dumpTauIdEffKineReweightHistogramsZtoMuTau
   + process.normalizeTauIdEffKineReweightHistogramsZtoMuTau
   + process.prodKineEventReweightsTauIdEffZtoMuTau
   + process.saveKineEventReweightsTauIdEffZtoMuTau
)

# print-out all python configuration parameter information
#print process.dumpPython()


  
