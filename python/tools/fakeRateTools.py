import FWCore.ParameterSet.Config as cms
import copy

def reconfigDQMFileLoader(dqmFileLoaderConfig, dqmDirectory):

    # configure attributes of DQMFileLoader module
    # corresponding to different signal/background processes
    for processName in dir(dqmFileLoaderConfig):
        processConfigEntry = getattr(dqmFileLoaderConfig, processName)
        
        if isinstance(processConfigEntry, cms.PSet):

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

def enableFakeRates_runZtoMuTau(process, method = None):

    # check validity of method parameter
    if method is None:
        raise ValueError("Undefined method Parameter !!")
    else:
        if method != "simple" and method != "CDF":
            raise ValueError("Invalid method Parameter !!")

    # set method parameter in fakeRateWeight producer modules
    process.bgEstFakeRateJetWeightsForMuTau.method = method
    process.bgEstFakeRateEventWeightsForMuTau.method = method

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
    # require muon and loosely selected tau-jet candidate to have opposite charges
    #
    # NOTE: because tau-jet candidate may well have charge != +1||-1,
    #       cannot require sum of muon + tau-jet charges to be zero;
    #       instead, require that charge sum of muon + leading track within tau-jet equals zero
    #
    changeCut(process, "selectedMuTauPairsZeroCharge", "(leg1.charge + leg2.leadTrack.charge) = 0")

    # add fake-rates to pat::Tau
    from TauAnalysis.RecoTools.patPFTauConfig_cfi import *
    setattr(allLayer1Taus.efficiencies, "bgEstFakeRateJetWeight", cms.InputTag('bgEstFakeRateJetWeightsForMuTau'))

    # weight events by fake-rate
    #
    # Note: special care is needed to avoid double-counting
    #       in case there is more than one (loosely selected) tau-jet candidate in the event
    #       when filling histograms that are sensitive to the tau-jet multiplicity
    #
    setattr(process.analyzeZtoMuTauEvents, "eventWeightSource", cms.VInputTag(cms.InputTag('bgEstFakeRateEventWeightsForMuTau')))

    # check if factorization is enabled;
    # if so, apply fake-rate event weights to analysis paths without/with muon isolation
    if hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation"):
        setattr(process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation, "eventWeightSource", cms.VInputTag(cms.InputTag('bgEstFakeRateEventWeightsForMuTau')))
    if hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation"):
        setattr(process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation, "eventWeightSource", cms.VInputTag(cms.InputTag('bgEstFakeRateEventWeightsForMuTau')))

    if hasattr(process, "tauHistManager"):
        setattr(process.tauHistManager, "tauJetWeightSource", cms.vstring("bgEstFakeRateJetWeight"))
    if hasattr(process, "diTauCandidateZmumuHypothesisHistManagerForMuTau"):
        setattr(process.diTauCandidateZmumuHypothesisHistManagerForMuTau, "lepton2WeightSource", cms.vstring("bgEstFakeRateJetWeight"))
    if hasattr(process, "diTauCandidateHistManagerForMuTau"):
        setattr(process.diTauCandidateHistManagerForMuTau, "diTauLeg2WeightSource", cms.vstring("bgEstFakeRateJetWeight"))    
