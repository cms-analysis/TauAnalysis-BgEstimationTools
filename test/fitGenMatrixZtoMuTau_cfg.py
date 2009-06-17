import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
#
# Apply "generalized matrix" method for data-driven background estimation
# to Z --> mu + tau-jet channel
#
# Author: Christian Veelken, UC Davis
#
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_cfi import *
from TauAnalysis.BgEstimationTools.bgEstNtupleDefinitionsZtoMuTau_cfi import *
from TauAnalysis.BgEstimationTools.bgEstBinGridZtoMuTau_cfi import *
from TauAnalysis.DQMTools.plotterStyleDefinitions_cfi import *

process = cms.Process('fitGenMatrixZtoMuTau')

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

process.compObjValCorrelationZtoMuTau_Ztautau = cms.EDAnalyzer("ObjValCorrelationAnalyzer",
    processName = cms.string("Ztautau"),                                
    #fileNames = cms.vstring(ntupleOutputDirectoryName.value() + 'bgEstNtuple.root'),
    fileNames = fileNames_Ztautau,
    treeName = cms.string("ntupleProducer/bgEstEvents"),
    branches = cms.vstring(
        'muonIso',
        'diTauAbsCharge',
        'diTauMt1MET'
    )
)

process.compObjValCorrelationZtoMuTau_Zmumu = copy.deepcopy(process.compObjValCorrelationZtoMuTau_Ztautau)
process.compObjValCorrelationZtoMuTau_Zmumu.processName = cms.string("Zmumu")
process.compObjValCorrelationZtoMuTau_Zmumu.fileNames = fileNames_Zmumu

process.compObjValCorrelationZtoMuTau_WplusJets = copy.deepcopy(process.compObjValCorrelationZtoMuTau_Ztautau)
process.compObjValCorrelationZtoMuTau_WplusJets.processName = cms.string("WplusJets")
process.compObjValCorrelationZtoMuTau_WplusJets.fileNames = fileNames_WplusJets

process.compObjValCorrelationZtoMuTau_qcdSum = copy.deepcopy(process.compObjValCorrelationZtoMuTau_Ztautau)
process.compObjValCorrelationZtoMuTau_qcdSum.processName = cms.string("qcdSum")
process.compObjValCorrelationZtoMuTau_qcdSum.fileNames = fileNames_WplusJets

process.compObjValCorrelationZtoMuTau = cms.Sequence( process.compObjValCorrelationZtoMuTau_Ztautau 
                                                     +process.compObjValCorrelationZtoMuTau_Zmumu
                                                     +process.compObjValCorrelationZtoMuTau_WplusJets
                                                     +process.compObjValCorrelationZtoMuTau_qcdSum )

process.bgEstFitZtoMuTau = cms.EDAnalyzer("GenMatrixBgEstFit",
    processes = cms.PSet(
        Ztautau = cms.PSet(
            fileNames = fileNames_Ztautau,
            #fixNorm = cms.bool(True),
            #fixP1 = cms.bool(True),
            #fixP2 = cms.bool(True),
            #fixP3 = cms.bool(True),
            drawOptions = drawOption_Ztautau
        ),
        #Zmumu = cms.PSet(
        #    fileNames = fileNames_Zmumu,
        #    fixNorm = cms.bool(True),
        #    fixP1 = cms.bool(True),
        #    fixP2 = cms.bool(True),
        #    fixP3 = cms.bool(True),
        #    drawOptions = drawOption_Zmumu
        #),
        WplusJets = cms.PSet(
            fileNames = fileNames_WplusJets,
            # no constraints
            drawOptions = drawOption_WplusJets
        #),
        #qcdSum = cms.PSet(
        #    fileNames = fileNames_qcdSum,
        #    # constrain OS/SS ratio to value of 1:1
        #    # (assumption of equal number of opposite-sign and like-sign muon + tau-jet pairs
        #    #  holds in first approximation, to be multiplied by correction factor
        #    #  determined by Monte Carlo later...)
        #    #P2 = cms.double(0.50),
        #    fixP2 = cms.bool(True),
        #    drawOptions = drawOption_QCD
        )
    ),

    # use "pseudo" data-samples consisting of all Monte Carlo processes for testing                      
    data = cms.vstring(fileNames_pseudoData),
                          
    treeName = cms.string("ntupleProducer/bgEstEvents"),
    treeSelection = cms.string("muonComp > 1. && muonTrkIP < 0.05"),                      
                          
    branches = genMatrixBinningZtoMuTau,
    #branchNameEventWeight = cms.string('eventWeight'),                                          

    output = cms.PSet(
        scaleFactors = cms.PSet(
            fileName = cms.string("bgEstFitProcessScaleFactors.py"),
            signalRegion = cms.vdouble(0., 0., 0.)
        ),
        controlPlots = cms.PSet(
            fileName = cms.string("fitGenMatrixZtoMuTau.png")
        )
    )
)                          

process.fitZtoMuTau = cms.Sequence( process.bgEstFitZtoMuTau )

process.p = cms.Path( process.compObjValCorrelationZtoMuTau
                     +process.fitZtoMuTau )

# print-out all python configuration parameter information
#print process.dumpPython()


  
