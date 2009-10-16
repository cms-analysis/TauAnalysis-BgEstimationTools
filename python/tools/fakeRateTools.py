import FWCore.ParameterSet.Config as cms
import copy

def reconfigDQMFileLoader(dqmFileLoaderConfig, plotDirectoryName_old, plotDirectoryName_new, dqmDirectory):

    # configure attributes of DQMFileLoader module
    # corresponding to different signal/background processes
    for processName in dir(dqmFileLoaderConfig):
        processConfigEntry = getattr(dqmFileLoaderConfig, processName)
        
        if isinstance(processConfigEntry, cms.PSet):
            if hasattr(processConfigEntry, "inputFileNames"):

                inputFileNames_old = getattr(processConfigEntry, "inputFileNames")
                inputFileNames_new = []
                
                for iInputFile in range(len(inputFileNames_old)):
                    directory_old = plotDirectoryName_old.value()
                    if not directory_old.endswith("/"):
                        directory_old += "/"

                    directory_new = plotDirectoryName_new.value()
                    if not directory_new.endswith("/"):
                        directory_new += "/"
                    
                    inputFileName_old = inputFileNames_old[iInputFile]
                    inputFileName_new = inputFileName_old.replace(directory_old, directory_new)

                    inputFileNames_new.append(inputFileName_new)

                setattr(processConfigEntry, "inputFileNames", cms.vstring(inputFileNames_new))

            if hasattr(processConfigEntry, "dqmDirectory_store"):
                if not dqmDirectory.endswith("/"):
                    dqmDirectory += "/"
                
                dqmDirectory_old = getattr(processConfigEntry, "dqmDirectory_store").value()
                dqmDirectory_new = dqmDirectory.replace("#PROCESSDIR#", dqmDirectory_old)

                setattr(processConfigEntry, "dqmDirectory_store", dqmDirectory_new)

#--------------------------------------------------------------------------------
# utility functions specific to application of fake-rate technique
# for data-driven background estimation to Z --> mu + tau-jet channel
#--------------------------------------------------------------------------------

def enableFakeRates_runZtoMuTau(process):

    # import utility function for changing cut values
    from TauAnalysis.Configuration.tools.changeCut import changeCut

    # disable cuts on tau id. discriminators
    changeCut(process, "selectedLayer1TausForMuTauLeadTrk", "tauID('leadingTrackFinding') > -1.")
    changeCut(process, "selectedLayer1TausForMuTauLeadTrkPt", "tauID('leadingTrackPtCut') > -1.")
    changeCut(process, "selectedLayer1TausForMuTauTrkIso", "tauID('trackIsolation') > -1.")
    changeCut(process, "selectedLayer1TausForMuTauEcalIso", "tauID('ecalIsolation') > -1.")
    changeCut(process, "selectedLayer1TausForMuTauProng", "signalTracks.size() > -1")
    changeCut(process, "selectedLayer1TausForMuTauCharge", "abs(charge) > -1")
    #changeCut(process, "selectedLayer1TausForMuTauMuonVeto", "tauID('againstMuon') > -1.")

    # add fake-rates to pat::Tau
    from TauAnalysis.RecoTools.patPFTauConfig_cfi import *
    setattr(allLayer1Taus.efficiencies, "bgEstFakeRateJetWeight", cms.InputTag("bgEstFakeRateJetWeights"))

    # weight events by fake-rate
    #
    # Note: special care is needed to avoid double-counting
    #       in case there is more than one (loosely selected) tau-jet candidate in the event
    #       when filling histograms that are sensitive to the tau-jet multiplicity
    #
    setattr(process.analyzeZtoMuTauEvents, "eventWeightSource", cms.VInputTag(cms.InputTag('bgEstFakeRateEventWeights')))

    if hasattr(process, "tauHistManager"):
        setattr(process.tauHistManager, "tauJetWeightSource", cms.vstring("bgEstFakeRateJetWeight"))
    if hasattr(process, "diTauCandidateZmumuHypothesisHistManagerForMuTau"):
        setattr(process.diTauCandidateZmumuHypothesisHistManagerForMuTau, "lepton2WeightSource", cms.vstring("bgEstFakeRateJetWeight"))
    if hasattr(process, "diTauCandidateHistManagerForMuTau"):
        setattr(process.diTauCandidateHistManagerForMuTau, "diTauLeg2WeightSource", cms.vstring("bgEstFakeRateJetWeight"))    
