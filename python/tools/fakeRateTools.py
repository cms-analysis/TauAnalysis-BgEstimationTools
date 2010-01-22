import FWCore.ParameterSet.Config as cms
import copy

# import utility function for changing cut values
from TauAnalysis.Configuration.tools.changeCut import changeCut

# import utility function to enable factorization
from TauAnalysis.Configuration.factorizationTools import enableFactorization_makeZtoMuTauPlots

#--------------------------------------------------------------------------------
# auxiliary functions needed for reconfiguration of analysis sequences
#--------------------------------------------------------------------------------

def setAnalyzerParameter(genAnalyzer, pluginName, parameterName, parameterValue):
    for analyzerPlugin in genAnalyzer.analyzers:
        if hasattr(analyzerPlugin, "pluginName"):
            analyzerPluginName = getattr(analyzerPlugin, "pluginName").value()
            if analyzerPluginName == pluginName:
                setattr(analyzerPlugin, parameterName, parameterValue)

def pruneAnalysisSequence(genAnalyzer):

    # disable filling of histograms after all stages of the event selection
    # except for the last occurence (after all cuts have been applied)

    lastEntry = {}
    for iAnalysisSequenceEntry in range(len(genAnalyzer.analysisSequence)):
        analysisSequenceEntry = genAnalyzer.analysisSequence[iAnalysisSequenceEntry]
        if hasattr(analysisSequenceEntry, "analyzers"):
            analyzerPlugins = getattr(analysisSequenceEntry, "analyzers")
            for analyzerPlugin in analyzerPlugins:
                analyzerPluginName = analyzerPlugin
                lastEntry[analyzerPluginName] = iAnalysisSequenceEntry

    prunedAnalysisSequence = []
    for iAnalysisSequenceEntry in range(len(genAnalyzer.analysisSequence)):
        analysisSequenceEntry = genAnalyzer.analysisSequence[iAnalysisSequenceEntry]
        if hasattr(analysisSequenceEntry, "analyzers"):
            # keep analyzer entry only in case it contains at least one histogram manager
            # not filled at a later stage of the event selection
            
            keepAnalyzers = []
            for lastEntryKey, lastEntryValue in lastEntry.items():
                if lastEntryValue == iAnalysisSequenceEntry:
                    keepAnalyzers.append(lastEntryKey)
            
            if len(keepAnalyzers) > 0:
                # keep analysis sequence entry, but fill only histograms
                # which are not filled at a later stage of the event selection 
                analysisSequenceEntry.analyzers = cms.vstring(keepAnalyzers)
                
                # in all cases, disable storing run and events numbers
                setattr(genAnalyzer.analysisSequence[iAnalysisSequenceEntry], "saveRunEventNumbers", cms.vstring())
                
                prunedAnalysisSequence.append(analysisSequenceEntry)
        else:
            # keep all filter entries
            prunedAnalysisSequence.append(analysisSequenceEntry)

    genAnalyzer.analysisSequence = cms.VPSet(prunedAnalysisSequence)

def addFakeRateAnalyzer(process, genAnalyzer, frType, bgEstFakeRateAnalysisSequence):

    bgEstFakeRateAnalyzer = copy.deepcopy(genAnalyzer)

    srcFakeRateEventWeight = cms.VInputTag(cms.InputTag("bgEstFakeRateEventWeights", frType))
    setattr(bgEstFakeRateAnalyzer, "eventWeightSource", srcFakeRateEventWeight)

    bgEstFakeRateAnalyzer.name = cms.string("".join([getattr(genAnalyzer, "name").value(), "_", frType]))

    srcFakeRateJetWeight = cms.vstring("".join(["bgEstFakeRateJetWeight", "_", frType]))
    setAnalyzerParameter(genAnalyzer, "tauHistManager", "tauJetWeightSource", srcFakeRateJetWeight)
    setAnalyzerParameter(genAnalyzer, "diTauCandidateZmumuHypothesisHistManagerForMuTau", "lepton2WeightSource", srcFakeRateJetWeight)
    setAnalyzerParameter(genAnalyzer, "diTauCandidateHistManagerForMuTau", "diTauLeg2WeightSource", srcFakeRateJetWeight)

    bgEstFakeRateAnalyzerName = "".join([genAnalyzer.label(), "_", frType])
    setattr(process, bgEstFakeRateAnalyzerName, bgEstFakeRateAnalyzer)
    bgEstFakeRateAnalyzer = getattr(process, bgEstFakeRateAnalyzerName)

    # add module to sequence
    if bgEstFakeRateAnalysisSequence is None:
        bgEstFakeRateAnalysisSequence = bgEstFakeRateAnalyzer
    else:
        bgEstFakeRateAnalysisSequence *= bgEstFakeRateAnalyzer

    return bgEstFakeRateAnalysisSequence

def getPSetAttributes(object):

    attributes = []
    
    for attribute in dir(object):
        
        # check that "attribute" is not an internal attribute or method of cms.PSet
        isInternalAttribute = False
        
        for classAttribute in dir(cms.PSet):
            if attribute == classAttribute:
                isInternalAttribute = True
        if attribute.startswith("_"):
            isInternalAttribute = True
            
        if not isInternalAttribute:
            attributes.append(attribute)

    return attributes

#--------------------------------------------------------------------------------
# auxiliary functions needed for reconfiguration of DQM file loader modules
#--------------------------------------------------------------------------------

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
#
# NOTE: in case factorization of the muon isolation efficiency
#       is used in order to improve the statistical precision of Monte Carlo estimates,
#       the function enableFakeRates_runZtoMuTau needs to be called **after**
#       the function enableFactorization_runZtoMuTau has been called
#
#--------------------------------------------------------------------------------

def enableFakeRates_runZtoMuTau(process):

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

    # get list of fake-rates types to be processed
    frTypes = getPSetAttributes(process.bgEstFakeRateJetWeights.frTypes)

    # add fake-rates to pat::Tau
    for frType in frTypes: 
        frLabel = "".join(["bgEstFakeRateJetWeight", "_", frType])
        frInputTag = cms.InputTag('bgEstFakeRateJetWeights', frType)
        setattr(process.allLayer1Taus.efficiencies, frLabel, frInputTag)

    # fill histograms only for events passing all event selection critera;
    # disable storing run and event numbers for events passing selection
    # (in order to save space in the .root files)
    if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and \
         hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
        pruneAnalysisSequence(process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation)
        pruneAnalysisSequence(process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation)
    else:
        pruneAnalysisSequence(process.analyzeZtoMuTauEvents)

    bgEstFakeRateAnalysisSequence = None  

    # duplicate analysis sequence:
    #  1.) tau id. discriminators not applied
    #  2.) events weighted by fake-rate
    # for each type of fake-rate weights given as function argument
    #
    # Note: special care is needed to avoid double-counting
    #       in case there is more than one (loosely selected) tau-jet candidate in the event
    #       when filling histograms that are sensitive to the tau-jet multiplicity
    #
    for frType in frTypes:
        # check if factorization is enabled;
        # if so, apply fake-rate event weights to analysis paths without/with muon isolation,
        # else apply fake-rate event weights to "standard" analysis path
        if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and 
             hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
            bgEstFakeRateAnalysisSequence = addFakeRateAnalyzer(process, process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation,
                                                                frType, bgEstFakeRateAnalysisSequence)
            bgEstFakeRateAnalysisSequence = addFakeRateAnalyzer(process, process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation,
                                                                frType, bgEstFakeRateAnalysisSequence)
        else:
            bgEstFakeRateAnalysisSequence = addFakeRateAnalyzer(process, process.analyzeZtoMuTauEvents,
                                                                frType, bgEstFakeRateAnalysisSequence)

    setattr(process, "bgEstFakeRateAnalysisSequence", cms.Sequence(bgEstFakeRateAnalysisSequence))

    if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and
         hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
        process.p.replace(process.analyzeZtoMuTauEvents_factorized, process.bgEstFakeRateAnalysisSequence)
    else:
        process.p.replace(process.analyzeZtoMuTauEvents, process.bgEstFakeRateAnalysisSequence)

def enableFakeRates_makeZtoMuTauPlots(process):

    # get list of fake-rates types to be processed
    process.load("TauAnalysis.RecoTools.fakeRateJetWeightProducer_cfi")
    frTypes = getPSetAttributes(process.bgEstFakeRateJetWeights.frTypes)

    addSequence = None
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
                moduleName_scale = "dqmHistScaler" + "_" + moduleType + "_" + processName
                module_scale = getattr(process, moduleName_scale)
                bgEstFakeRateModuleName_scale = moduleName_scale + "_" + frType
                setattr(process, bgEstFakeRateModName_scale, module_scale)

                if addSequence is None:
                    addSequence = module_scale
                else:
                    addSequence *= module_scale

        moduleName_qcdSum = "addBgEstFakeRateZtoMuTau_qcdSum_tauFakeRate"
        module_qcdSum = getattr(process, moduleName_qcdSum)
        bgEstFakeRateModuleName_qcdSum = moduleName_qcdSum + "_" + frType
        setattr(process, bgEstFakeRateModuleName_qcdSum, module_qcdSum)
        process.addBgEstFakeRateZtoMuTau_tauFakeRate._seq = process.addBgEstFakeRateZtoMuTau_tauFakeRate._seq * module_qcdSum

        if hasattr(process, "addBgEstFakeRateZtoMuTau_smSum_tauFakeRate"):
            moduleName_smSum = "addBgEstFakeRateZtoMuTau_smSum_tauFakeRate"
            module_smSum = getattr(process, moduleName_smSum)
            bgEstFakeRateModuleName_smSum = moduleName_smSum + "_" + frType
            setattr(process, bgEstFakeRateModuleName_smSum, module_smSum)
            process.addBgEstFakeRateZtoMuTau_tauFakeRate._seq = process.addBgEstFakeRateZtoMuTau_tauFakeRate._seq * module_smSum

        moduleName_dump = "dumpZtoMuTau"
        module_dump = getattr(process, moduleName_dump)
        bgEstFakeRateModule_dump = copy.deepcopy(module_dump)

        processNames = getPSetAttributes(bgEstFakeRateModule_dump.dqmDirectories)
        for processName in processNames:
            directory = getattr(bgEstFakeRateModule_dump.dqmDirectories, processName)
            bgEstFakeRateDirectory = directory.replace("/zMuTauAnalyzer/", "/zMuTauAnalyzer" + "_" + frType + "/")
            setattr(bgEstFakeRateModule_dump.dqmDirectories, processName, bgEstFakeRateDirectory)
        bgEstFakeRateModuleName_dump = moduleName_dump + "_" + frType
        setattr(process, bgEstFakeRateModuleName_dump, module_dump)

        if dumpSequence is None:
            dumpSequence = module_dump
        else:
            dumpSequence *= module_dump

    process.addZtoMuTau = cms.Sequence(addSequence)
    process.dumpZtoMuTau = cms.Sequence(dumpSequence)
