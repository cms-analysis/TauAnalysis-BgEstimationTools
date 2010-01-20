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
# auxiliary functions needed for reconfiguration of analysis sequences
# of Z --> e + mu, Z --> e + tau-jet, Z --> mu + tau-jet, ... channels
#--------------------------------------------------------------------------------

def setAnalyzerParameter(genericAnalyzerModule, pluginName, parameterName, parameterValue):
    for analyzerPlugin in genericAnalyzerModule.analyzers:
        if hasattr(analyzerPlugin, "pluginName"):
            analyzerPluginName = getattr(analyzerPlugin, "pluginName")
            if analyzerPluginName == pluginName:
                setattr(analyzerPlugin, parameterName, parameterValue)

def addFakeRateAnalyzer(process, analyzer, bgEstFakeRateName, bgEstFakeRateAnalysisSequence):
    
    bgEstFakeRateAnalyzer = copy.deepcopy(analyzer)
    
    srcFakeRateEventWeight = cms.VInputTag(cms.InputTag("bgEstFakeRateJetWeights", bgEstFakeRateName))
    setattr(bgEstFakeRateAnalyzer, "eventWeightSource", srcFakeRateEventWeight)

    bgEstFakeRateAnalyzerName = analyzer.label() + "_" + bgEstFakeRateName
    bgEstFakeRateAnalyzer.name = cms.string(bgEstFakeRateAnalyzerName)

    setAnalyzerParameter(analyzer, "tauHistManager", "tauJetWeightSource", srcFakeRateEventWeight)
    setAnalyzerParameter(analyzer, "diTauCandidateZmumuHypothesisHistManagerForMuTau", "lepton2WeightSource", srcFakeRateEventWeight)
    setAnalyzerParameter(analyzer, "diTauCandidateHistManagerForMuTau", "diTauLeg2WeightSource", srcFakeRateEventWeight)

    setattr(process, bgEstFakeRateAnalyzerName, bgEstFakeRateAnalyzer)

    # add module to sequence
    if bgEstFakeRateAnalysisSequence == None:
        bgEstFakeRateAnalysisSequence = analyzer
    else:
        bgEstFakeRateAnalysisSequence *= getattr(process, bgEstFakeRateAnalyzerName)

    return bgEstFakeRateAnalysisSequence

#--------------------------------------------------------------------------------
# utility functions specific to application of fake-rate technique
# for data-driven background estimation to Z --> mu + tau-jet channel
#
# NOTE: in case factorization of the muon isolation efficiency
#       is used in order to improve the statistical precision of Monte Carlo estimates,
#       the function enableFakeRates_runZtoMuTau needs to be called **after**
#       the function enableFactorization_runZtoMuTau has been called
#
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
    bgEstFakeRates = [ [ "bgEstFakeRateJetWeightQCDmuEnriched",         "qcdMuEnriched"         ],
                       [ "bgEstFakeRateJetWeightQCDdiJetLeadJet",       "qcdDiJetLeadJet"       ],
                       [ "bgEstFakeRateJetWeightQCDdiJetSecondLeadJet", "qcdDiJetSecondLeadJet" ],
                       [ "bgEstFakeRateJetWeightWplusJets",             "WplusJets"             ] ]
    for bgEstFakeRate in bgEstFakeRates:
        setattr(allLayer1Taus.efficiencies, bgEstFakeRate[0], cms.InputTag("bgEstFakeRateJetWeights", bgEstFakeRate[1]))

    # weight events by fake-rate
    #
    # Note: special care is needed to avoid double-counting
    #       in case there is more than one (loosely selected) tau-jet candidate in the event
    #       when filling histograms that are sensitive to the tau-jet multiplicity
    #
    bgEstFakeRateAnalysisSequence = None
    for bgEstFakeRate in bgEstFakeRates:
        
         # check if factorization is enabled;
         # if so, apply fake-rate event weights to analysis paths without/with muon isolation,
         # else apply fake-rate event weights to "standard" analysis path
         if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and 
              hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
             bgEstFakeRateAnalysisSequence = addFakeRateAnalyzer(process, process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation,
                                                                 bgEstFakeRate[1], bgEstFakeRateAnalysisSequence)
             bgEstFakeRateAnalysisSequence = addFakeRateAnalyzer(process, process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation,
                                                                 bgEstFakeRate[1], bgEstFakeRateAnalysisSequence)
         else:
             bgEstFakeRateAnalysisSequence = addFakeRateAnalyzer(process, process.analyzeZtoMuTauEvents,
                                                                 bgEstFakeRate[1], bgEstFakeRateAnalysisSequence)

    setattr(process, "bgEstFakeRateAnalysisSequence", cms.Sequence(bgEstFakeRateAnalysisSequence))

    if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and
         hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
        process.analyzeZtoMuTauEventsAll = cms.Sequence(process.analyzeZtoMuTauEvents_factorized * process.bgEstFakeRateAnalysisSequence)
        process.p.replace(process.analyzeZtoMuTauEvents_factorized, process.analyzeZtoMuTauEventsAll)
    else:
        process.analyzeZtoMuTauEventsAll = cms.Sequence(process.analyzeZtoMuTauEvents * process.bgEstFakeRateAnalysisSequence)
        process.p.replace(process.analyzeZtoMuTauEvents, process.analyzeZtoMuTauEventsAll)
    
