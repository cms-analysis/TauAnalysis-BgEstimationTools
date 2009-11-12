import FWCore.ParameterSet.Config as cms
import copy

# import utility function for changing cut values
from TauAnalysis.Configuration.tools.changeCut import changeCut

# import utility function to enable factorization
from TauAnalysis.Configuration.factorizationTools import enableFactorization_makeZtoMuTauPlots

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

def addGenAnalyzer(process, genAnalyzerName, genAnalyzerModuleName, label,
                   pythonSequence = None, eventWeightSource = None, histManagerNames = None, frType = None):
    oldGenAnalyzer = getattr(process, genAnalyzerName)
    
    newGenAnalyzer = copy.deepcopy(oldGenAnalyzer)
    newGenAnalyzer.name = cms.string(genAnalyzerModuleName + "_" + label)

    if eventWeightSource is not None:
        setattr(newGenAnalyzer, "eventWeightSource", eventWeightSource)

    if histManagerNames is not None:
        for histManagerName in histManagerNames:            
            oldHistManagerName = histManagerName
            oldHistManager = getattr(process, oldHistManagerName)

            newHistManagerName = histManagerName + "_" + frType
            newHistManager = getattr(process, newHistManagerName)

            histManagers = newGenAnalyzer.analyzers;
            for iHistManager in range(len(histManagers)):
                histManager = histManagers[iHistManager]
                if histManager.pluginName.value() == oldHistManagerName:
                    histManagers[iHistManager] = newHistManager
            
            for analysisSequenceEntry in newGenAnalyzer.analysisSequence:
                if hasattr(analysisSequenceEntry, "analyzers"):
                    histManagerNames = getattr(analysisSequenceEntry, "analyzers")
                    for iHistManagerName in range(len(histManagerNames)):
                        if histManagerNames[iHistManagerName] == oldHistManagerName:
                            histManagerNames[iHistManagerName] = newHistManagerName
                if hasattr(analysisSequenceEntry, "replace"):
                    replaceStatements = getattr(analysisSequenceEntry, "replace")
                    for iReplaceStatement in range(len(replaceStatements)):
                        oldReplaceStatement = replaceStatements[iReplaceStatement]
                        newReplaceStatement = oldReplaceStatement.replace(oldHistManagerName, newHistManagerName)
                        replaceStatements[iReplaceStatement] = newReplaceStatement

    setattr(process, genAnalyzerName + "_" + label, newGenAnalyzer);

    if pythonSequence is not None:    
        pythonSequence._seq = pythonSequence._seq * newGenAnalyzer

    return newGenAnalyzer

def disableHistogramFilling(genAnalyzer):
    for iAnalysisSequenceEntry in range(len(genAnalyzer.analysisSequence)):

        # disable filling of histograms after all stages of the event selection
        # except after all cuts have been applied
        if iAnalysisSequenceEntry < (len(genAnalyzer.analysisSequence) - 1):
            if hasattr(genAnalyzer.analysisSequence[iAnalysisSequenceEntry], "analyzers"):
                setattr(genAnalyzer.analysisSequence[iAnalysisSequenceEntry], "analyzers", cms.vstring())
                setattr(genAnalyzer.analysisSequence[iAnalysisSequenceEntry], "replace", cms.vstring())

        # disable storing run and events numbers for events passing selection
        if hasattr(genAnalyzer.analysisSequence[iAnalysisSequenceEntry], "saveRunEventNumbers"):
            setattr(genAnalyzer.analysisSequence[iAnalysisSequenceEntry], "saveRunEventNumbers", cms.vstring())

#--------------------------------------------------------------------------------
# utility functions specific to application of fake-rate technique
# for data-driven background estimation to Z --> mu + tau-jet channel
#--------------------------------------------------------------------------------

def enableFakeRates_runZtoMuTau(process, frTypes = None, method = None):

    # check validity of frTypes parameter
    if frTypes is None or len(frTypes) == 0:
        raise ValueError("Undefined frTypes Parameter !!")

    # check validity of method parameter
    if method is None:
        raise ValueError("Undefined method Parameter !!")
    else:
        if method != "simple" and method != "CDF":
            raise ValueError("Invalid method Parameter !!")

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
    # NOTE:
    #  (1) because tau-jet candidate may well have charge != +1||-1,
    #      cannot require sum of muon + tau-jet charges to be zero;
    #      instead, require that charge sum of muon + leading track within tau-jet equals zero
    #  (2) this requirement implies that there is at least one track within the tau-jet
    #     (leading to a small underestimation of backgrounds)
    #
    changeCut(process, "selectedMuTauPairsZeroCharge", "leg2.leadTrack.isNonnull & (leg1.charge + leg2.leadTrack.charge) = 0")

    # set method parameter in fakeRateWeight producer modules
    process.bgEstFakeRateJetWeightsForMuTau.method = method
    process.bgEstFakeRateEventWeightsForMuTau.method = method

    # add fake-rates to pat::Tau
    for frType in dir(process.bgEstFakeRateJetWeightsForMuTau.frTypes):
        
        # check that "attribute" is not an internal attribute or method of cms.PSet
        isInternalAttribute = False
        for classAttribute in dir(cms.PSet):
            if frType == classAttribute:
                isInternalAttribute = True
        if frType.startswith("_"):
                isInternalAttribute = True             
        if not isInternalAttribute:
            frLabel = "bgEstFakeRateJetWeight" + "_" + frType
            frInputTag = cms.InputTag('bgEstFakeRateJetWeightsForMuTau', frType)
            setattr(process.allLayer1Taus.efficiencies, frLabel, frInputTag)

    # fill histograms only for events passing all event selection critera;
    # disable storing run and event numbers for events passing selection
    # (in order to save space in the .root files)
    if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and \
         hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
        disableHistogramFilling(process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation)
        disableHistogramFilling(process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation)
    else:
        disableHistogramFilling(process.analyzeZtoMuTauEvents)

    # duplicate analysis sequence:
    #  1.) fake-rate weights not applied
    #  2.) events weighted by fake-rate
    # for each type of fake-rate weights given as function argument
    #
    # Note: special care is needed to avoid double-counting
    #       in case there is more than one (loosely selected) tau-jet candidate in the event
    #       when filling histograms that are sensitive to the tau-jet multiplicity
    #
    # check if factorization is enabled;
    # if so, apply fake-rate event weights to analysis paths without/with muon isolation
    analyzeZtoMuTauFakeRateSequence = None
    if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and \
         hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
        seqEntry = addGenAnalyzer(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation",
                                  "zMuTauAnalyzer_factorizedWithoutMuonIsolation", "noWeights")
        analyzeZtoMuTauFakeRateSequence = cms.Sequence( seqEntry )
        addGenAnalyzer(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation",
                       "zMuTauAnalyzer_factorizedWithMuonIsolation", "noWeights", analyzeZtoMuTauFakeRateSequence)
    else:
        seqEntry = addGenAnalyzer(process, "analyzeZtoMuTauEvents", "zMuTauAnalyzer", "noWeights")
        analyzeZtoMuTauFakeRateSequence = cms.Sequence( seqEntry )
    for frType in frTypes:
        if hasattr(process, "tauHistManager"):
            tauHistManager = copy.deepcopy(process.tauHistManager)
            tauHistManager.pluginName = cms.string(process.tauHistManager.pluginName.value() + "_" + frType)
            setattr(tauHistManager, "tauJetWeightSource", cms.vstring("bgEstFakeRateJetWeight" + "_" + frType))
            setattr(process, "tauHistManager" + "_" + frType, tauHistManager)
        if hasattr(process, "diTauCandidateHistManagerForMuTau"):
            diTauCandidateHistManager = copy.deepcopy(process.diTauCandidateHistManagerForMuTau)
            diTauCandidateHistManager.pluginName = cms.string(process.diTauCandidateHistManagerForMuTau.pluginName.value() + "_" + frType)
            setattr(diTauCandidateHistManager, "diTauLeg2WeightSource", cms.vstring("bgEstFakeRateJetWeight" + "_" + frType)) 
            setattr(process, "diTauCandidateHistManagerForMuTau" + "_" + frType, diTauCandidateHistManager)
        if hasattr(process, "diTauCandidateZmumuHypothesisHistManagerForMuTau"):
            diTauCandidateZmumuHypothesisHistManager = copy.deepcopy(process.diTauCandidateZmumuHypothesisHistManagerForMuTau)
            diTauCandidateZmumuHypothesisHistManager.pluginName = cms.string(process.diTauCandidateZmumuHypothesisHistManagerForMuTau.pluginName.value() + "_" + frType)
            setattr(diTauCandidateZmumuHypothesisHistManager, "lepton2WeightSource", cms.vstring("bgEstFakeRateJetWeight" + "_" + frType))
            setattr(process, "diTauCandidateZmumuHypothesisHistManagerForMuTau" + "_" + frType, diTauCandidateZmumuHypothesisHistManager)

        frInputTag = cms.VInputTag( cms.InputTag('bgEstFakeRateEventWeightsForMuTau', frType) )
        histManagers = [ "tauHistManager", "diTauCandidateHistManagerForMuTau", "diTauCandidateZmumuHypothesisHistManagerForMuTau" ]
        if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and
             hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
            addGenAnalyzer(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation",
                           "zMuTauAnalyzer_factorizedWithoutMuonIsolation", "frWeights" + "_" + frType, analyzeZtoMuTauFakeRateSequence,
                           frInputTag, histManagers, frType)
            addGenAnalyzer(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation",
                           "zMuTauAnalyzer_factorizedWithMuonIsolation", "frWeights" + "_" + frType, analyzeZtoMuTauFakeRateSequence,
                           frInputTag, histManagers, frType)          
        else:
            addGenAnalyzer(process, "analyzeZtoMuTauEvents",
                           "zMuTauAnalyzer", "frWeights" + "_" + frType, analyzeZtoMuTauFakeRateSequence,
                           frInputTag, histManagers, frType)

    setattr(process, "analyzeZtoMuTauFakeRateSequence", analyzeZtoMuTauFakeRateSequence)

    if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and
         hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
        process.p.replace(process.analyzeZtoMuTauEvents_factorized, process.analyzeZtoMuTauFakeRateSequence)
    else:
        process.p.replace(process.analyzeZtoMuTauEvents, process.analyzeZtoMuTauFakeRateSequence)

def enableFakeRates_makeZtoMuTauPlots(process, frTypes = None)

    isFirstModule_addSequence = True
    isFirstModule_dumpSequence = True
    dumpSequence = None

    for frType in frTypes:
        enableFactorization_makeZtoMuTauPlots(process,
             dqmDirectoryIn_InclusivePPmuX = 'tauFakeRate/harvested/InclusivePPmuX/zMuTauAnalyzer' + '_' + frType,
             dqmDirectoryOut_InclusivePPmuX = 'tauFakeRate/harvested/InclusivePPmuX_factorized/zMuTauAnalyzer' + '_' + frType,
             dqmDirectoryIn_PPmuXptGt20 = 'tauFakeRate/harvested/PPmuXptGt20/zMuTauAnalyzer' + '_' + frType,
             dqmDirectoryOut_PPmuXptGt20 = 'tauFakeRate/harvested/PPmuXptGt20_factorized/zMuTauAnalyzer' + '_' + frType,
             moduleName_addZtoMuTau_qcdSum = "addBgEstFakeRateZtoMuTau_qcdSum_tauFakeRate",
             moduleName_addZtoMuTau = "addBgEstFakeRateZtoMuTau_tauFakeRate")

        for processName in [ "InclusivePPmuX", "PPmuXptGt20" ]:
            for moduleType in [ "plotsFactorizedTightEvtSel",
                                "filterStatFactorizedTightEvtSel",
                                "plotsFactorizedLooseEvtSel",
                                "filterStatFactorizedLooseEvtSel" ]:
                oldModuleName_scale = "dqmHistScaler" + "_" + moduleType + "_" + processName
                module_scale = getattr(process, oldModuleName_scale)
                newModuleName_scale = oldModuleName_scale + "_" + frType
                setattr(process, newModName_scale, module_scale)

                if isFirstModule_addSequence:
                     process.addZtoMuTau = cms.Sequence( module_scale )
                     isFirstModule_addSequence = False
                else:
                     process.addZtoMuTau._seq = process.addZtoMuTau._seq * module_scale

        oldModuleName_qcdSum = "addBgEstFakeRateZtoMuTau_qcdSum_tauFakeRate"
        module_qcdSum = getattr(process, oldModuleName_qcdSum)
        newModuleName_qcdSum = oldModuleName_qcdSum + "_" + frType
        setattr(process, newModuleName_qcdSum, module_qcdSum)
        process.addBgEstFakeRateZtoMuTau_tauFakeRate._seq = process.addBgEstFakeRateZtoMuTau_tauFakeRate._seq * module_qcdSum

        if hasattr(process, "addBgEstFakeRateZtoMuTau_smSum_tauFakeRate"):
            oldModuleName_smSum = "addBgEstFakeRateZtoMuTau_smSum_tauFakeRate"
            module_smSum = getattr(process, oldModuleName_smSum)
            newModuleName_smSum = oldModuleName_smSum + "_" + frType
            setattr(process, newModuleName_smSum, module_smSum)
            process.addBgEstFakeRateZtoMuTau_tauFakeRate._seq = process.addBgEstFakeRateZtoMuTau_tauFakeRate._seq * module_smSum

        oldModuleName_dump = "dumpZtoMuTau"
        oldModule_dump = getattr(process, oldModuleName_dump)
        newModule_dump = copy.deepcopy(oldModule_dump)
        for processName in dir(newModule_dump.dqmDirectories):            
            # check that "attribute" is not an internal attribute or method of cms.PSet
            isInternalAttribute = False
            for classAttribute in dir(cms.PSet):
                if frType == classAttribute:
                    isInternalAttribute = True
            if frType.startswith("_"):
                isInternalAttribute = True             
            if not isInternalAttribute:
                 oldDirectory = getattr(newModule_dump.dqmDirectories, processName)
                 newDirectory = oldDirectory.replace("/zMuTauAnalyzer/", "/zMuTauAnalyzer" + "_" + frType + "/")
                 setattr(newModule_dump.dqmDirectories, processName, newDirectory)
        newModuleName_dump = oldModuleName_dump + "_" + frType
        setattr(process, newModuleName_dump, module_dump)

        if isFirstModule_dumpSequence:
            dumpSequence = cms.Sequence( module_dump )
                isFirstModule_dumpSequence = False
            else:
                dumpSequence._seq = dumpSequence._seq * module_dump

    process.dumpZtoMuTau = dumpSequence
