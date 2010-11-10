import FWCore.ParameterSet.Config as cms
import copy

# import utility function for changing cut values
from TauAnalysis.Configuration.tools.changeCut import changeCut

# import utility function to enable factorization
from TauAnalysis.Configuration.tools.factorizationTools import enableFactorization_makeZtoMuTauPlots_grid

# import configuration parameters of histogram manager
# for validation of tau id. efficiencies/fake-rates
from TauAnalysis.BgEstimationTools.tauIdEffValidationHistManager_cfi import tauIdEffValidationHistManager

# import utility function for add histogram manager to analysis sequence
from TauAnalysis.Configuration.tools.analysisSequenceTools import addAnalyzer

#--------------------------------------------------------------------------------
# auxiliary functions needed for reconfiguration of analysis sequences
#--------------------------------------------------------------------------------

def setAnalyzerParameter(genAnalyzerModule, pluginName, parameterName, parameterValue):
    for analyzerPlugin in genAnalyzerModule.analyzers:
        if hasattr(analyzerPlugin, "pluginName"):
            analyzerPluginName = getattr(analyzerPlugin, "pluginName").value()
            if analyzerPluginName == pluginName:
                setattr(analyzerPlugin, parameterName, parameterValue)

def pruneAnalysisSequence(genAnalyzerModule):

    # disable filling of histograms after all stages of the event selection
    # except for the last occurence (after all cuts have been applied)
    
    lastEntry = {}
    for iAnalysisSequenceEntry in range(len(genAnalyzerModule.analysisSequence)):
        analysisSequenceEntry = genAnalyzerModule.analysisSequence[iAnalysisSequenceEntry]
        if hasattr(analysisSequenceEntry, "analyzers"):
            analyzerPlugins = getattr(analysisSequenceEntry, "analyzers")
            for analyzerPlugin in analyzerPlugins:
                analyzerPluginName = analyzerPlugin
                lastEntry[analyzerPluginName] = iAnalysisSequenceEntry

    prunedAnalysisSequence = []
    for iAnalysisSequenceEntry in range(len(genAnalyzerModule.analysisSequence)):
        analysisSequenceEntry = genAnalyzerModule.analysisSequence[iAnalysisSequenceEntry]
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
                setattr(genAnalyzerModule.analysisSequence[iAnalysisSequenceEntry], "saveRunLumiSectionEventNumbers", cms.vstring())
                
                prunedAnalysisSequence.append(analysisSequenceEntry)
        else:
            # keep all filter entries,
            # but disable saving of run and event numbers of events passing filter
            setattr(analysisSequenceEntry, "saveRunLumiSectionEventNumbers", cms.vstring(''))
            prunedAnalysisSequence.append(analysisSequenceEntry)

    genAnalyzerModule.analysisSequence = cms.VPSet(prunedAnalysisSequence)

def disableEventDump(genAnalyzerModule):

    # disable event print-out

    disabledEventDump = copy.deepcopy(genAnalyzerModule.eventDumps[0])
    disabledEventDump.output = cms.string("std::cout")
    disabledEventDump.triggerConditions = cms.vstring()
    genAnalyzerModule.eventDumps[0] = disabledEventDump

def makeGenAnalyzerModule(process, genAnalyzerModule, label):

    analyzer = copy.deepcopy(genAnalyzerModule)

    analyzer.name = cms.string("".join([getattr(genAnalyzerModule, "name").value(), "_", label]))

    analyzerName = "".join([genAnalyzerModule.label(), "_", label])
    setattr(process, analyzerName, analyzer)
    analyzer = getattr(process, analyzerName)

    return analyzer

def addGenAnalyzerModule(process, genAnalyzerModule, label, analysisSequence):

    analyzer = makeGenAnalyzerModule(process, genAnalyzerModule, label)

    # add module to sequence
    if analysisSequence is None:
        analysisSequence = analyzer
    else:
        analysisSequence *= analyzer

    return analysisSequence

def addFakeRateGenAnalyzerModule(process, genAnalyzerModule, frType, bgEstFakeRateAnalysisSequence):

    bgEstFakeRateAnalyzer = makeGenAnalyzerModule(process, genAnalyzerModule, "fr" + "_" + frType)

    srcFakeRateEventWeight = cms.VInputTag(cms.InputTag("bgEstFakeRateEventWeights", frType))
    setattr(bgEstFakeRateAnalyzer, "eventWeightSource", srcFakeRateEventWeight)

    srcFakeRateJetWeight = cms.vstring("".join(["bgEstFakeRateJetWeight", "_", frType]))
    setAnalyzerParameter(bgEstFakeRateAnalyzer, "tauHistManager", "tauJetWeightSource", srcFakeRateJetWeight)
    setAnalyzerParameter(bgEstFakeRateAnalyzer, "diTauCandidateZmumuHypothesisHistManagerForMuTau", "lepton2WeightSource", srcFakeRateJetWeight)
    setAnalyzerParameter(bgEstFakeRateAnalyzer, "diTauCandidateHistManagerForMuTau", "diTauLeg2WeightSource", srcFakeRateJetWeight)

    # add module to sequence
    if bgEstFakeRateAnalysisSequence is None:
        bgEstFakeRateAnalysisSequence = bgEstFakeRateAnalyzer
    else:
        bgEstFakeRateAnalysisSequence *= bgEstFakeRateAnalyzer

    return bgEstFakeRateAnalysisSequence

def makeDataBinningDumpSequence(process, dqmDirectory, processSubDirectories, frSubDirectories, moduleLabel):

    dataBinningDumpAnalysisSequence = None

    for processName, processSubDirectory in processSubDirectories.items():

        module = cms.EDAnalyzer("DQMDumpBinningResults",
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService"),
                dqmDirectories = cms.PSet()
            )
        )

        for frType, frSubDirectory in frSubDirectories.items():
            
            dqmDirectory_i = dqmDirectory
            dqmDirectory_i = dqmDirectory_i.replace("#PROCESSDIR#", processSubDirectory.value())
            dqmDirectory_i = dqmDirectory_i.replace("#FAKERATEDIR#", frSubDirectory)

            setattr(module.binningService.dqmDirectories, frType, cms.string(dqmDirectory_i))

        moduleName = "".join(["dumpDataBinningBgEstFakeRateZtoMuTau", "_", processName, "_", moduleLabel])
        setattr(process, moduleName, module)

        module = getattr(process, moduleName)

        if dataBinningDumpAnalysisSequence is None:
            dataBinningDumpAnalysisSequence = module
        else:
            dataBinningDumpAnalysisSequence *= module

    return dataBinningDumpAnalysisSequence

def makeFilterStatTableDumpSequence(process, dqmDirectory, processSubDirectories, frSubDirectories, moduleLabel):

    filterStatTableDumpAnalysisSequence = None

    for processName, processSubDirectory in processSubDirectories.items():

        module = cms.EDAnalyzer("DQMDumpFilterStatisticsTables",
            dqmDirectories = cms.PSet(),
            #columnsSummaryTable = cms.vstring("Passed", "cumul. Efficiency", "margin. Efficiency"),
            columnsSummaryTable = cms.vstring("Passed"),
            printSummaryTableOnly = cms.bool(True)                    
        )

        for frType, frSubDirectory in frSubDirectories.items():
            
            dqmDirectory_i = dqmDirectory
            dqmDirectory_i = dqmDirectory_i.replace("#PROCESSDIR#", processSubDirectory.value())
            dqmDirectory_i = dqmDirectory_i.replace("#FAKERATEDIR#", frSubDirectory)

            setattr(module.dqmDirectories, frType, cms.string(dqmDirectory_i))

        moduleName = "".join(["dumpFilterStatTableBgEstFakeRateZtoMuTau", "_", processName, "_", moduleLabel])
        setattr(process, moduleName, module)

        module = getattr(process, moduleName)

        if filterStatTableDumpAnalysisSequence is None:
            filterStatTableDumpAnalysisSequence = module
        else:
            filterStatTableDumpAnalysisSequence *= module

    return filterStatTableDumpAnalysisSequence

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
# auxiliary functions needed for adding fake-rate producer modules
# to PAT production sequence (**not** channel specific)
#--------------------------------------------------------------------------------

def configureFakeRateWeightProduction(process, method = None, preselPFTauJetSource = 'shrinkingConePFTauProducer',
                                      patTauJetCut = "tauID('againstElectron') > 0.5 & tauID('againstMuon') > 0.5"):

    # check validity of method parameter
    if method is None:
        raise ValueError("Undefined method Parameter !!")
    else:
        if method != "simple" and method != "CDF":
            raise ValueError("Invalid method Parameter !!")
    
    # compute fake-rate weights
    #
    # NOTE: jet weights are computed for all (shrinking signal cone) reco::PFTaus,
    #       but only those tau-jet candidates passing preselection on PAT level
    #       must enter event weight computation !!
    #    
    process.load("RecoTauTag.TauAnalysisTools.PFTauEfficiencyAssociator_cfi")
    process.producePrePat._seq = process.producePrePat._seq * process.associateTauFakeRates
    
    process.load("TauAnalysis.BgEstimationTools.fakeRateJetWeightProducer_cfi")    
    process.bgEstFakeRateJetWeights.preselTauJetSource = cms.InputTag(preselPFTauJetSource)
    process.bgEstFakeRateJetWeights.method = method
    process.producePrePat._seq = process.producePrePat._seq * process.bgEstFakeRateJetWeights

    process.tausForFakeRateEventWeights = cms.EDFilter("PATTauSelector",
        src = cms.InputTag('selectedPatTausForMuTauLeadTrkPtCumulative'),               
        cut = cms.string(patTauJetCut),
        filter = cms.bool(False)
    )
    
    process.load("TauAnalysis.BgEstimationTools.fakeRateEventWeightProducer_cfi")
    process.bgEstFakeRateEventWeights.preselTauJetSource = cms.InputTag('tausForFakeRateEventWeights')
    process.bgEstFakeRateEventWeights.method = method
    process.produceFakeRateEventWeights = cms.Sequence(process.tausForFakeRateEventWeights + process.bgEstFakeRateEventWeights)
    process.producePatTupleZtoMuTauSpecific._seq = process.producePatTupleZtoMuTauSpecific._seq * process.produceFakeRateEventWeights

    # add fake-rates to pat::Tau
    frTypes = getPSetAttributes(process.bgEstFakeRateJetWeights.frTypes)
    for frType in frTypes: 
        frLabel = "".join(["bgEstFakeRateJetWeight", "_", frType])
        frInputTag = cms.InputTag('bgEstFakeRateJetWeights', frType)
        setattr(process.patTaus.efficiencies, frLabel, frInputTag)
    process.patTaus.addEfficiencies = cms.bool(True)    

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

def enableFakeRates_runZtoMuTau(process, method = None):

    # check validity of method parameter
    if method is None:
        raise ValueError("Undefined method Parameter !!")
    else:
        if method != "simple" and method != "CDF":
            raise ValueError("Invalid method Parameter !!")

    # compute fake-rate weights
    configureFakeRateWeightProduction(process, method = method)
        
    # disable cuts on tau id. discriminators
    #
    # NOTE: tau lead. track finding and lead. track Pt discriminators
    #       must **not** be disabled, as these discriminators are already applied at the skimming stage !!
    #       Instead, need to apply TauAnalysis specific efficiency/fake-rate values,
    #       which represent the probability for a tau-jet candidate
    #       passing the lead. track finding and lead. track Pt discriminators
    #       to pass the track && ECAL isolation, 1||3 tracks in signal cone and charge = +/- 1 requirements as well
    #
    #changeCut(process, "selectedPatTausForMuTauLeadTrk", "tauID('leadingTrackFinding') > -1.")
    #changeCut(process, "selectedPatTausForMuTauLeadTrkPt", "tauID('leadingTrackPtCut') > -1.")
    changeCut(process, "selectedPatTausForMuTauTaNCdiscr", "tauID('byTaNCfrQuarterPercent') > -1.e+3")
    changeCut(process, "selectedPatTausForMuTauTrkIso", "tauID('trackIsolation') > -1.")
    changeCut(process, "selectedPatTausForMuTauEcalIso", "tauID('ecalIsolation') > -1.")    
    changeCut(process, "selectedPatTausForMuTauProng", "signalPFChargedHadrCands.size() > -1")
    changeCut(process, "selectedPatTausForMuTauCharge", "abs(charge) > -1")
    #changeCut(process, "selectedPatTausForMuTauMuonVeto", "tauID('againstMuon') > -1.")
    # require muon and loosely selected tau-jet candidate to have opposite charges
    #
    # NOTE:
    #  (1) because tau-jet candidate may well have charge != +1||-1,
    #      cannot require sum of muon + tau-jet charges to be zero;
    #      instead, require that charge sum of muon + leading track within tau-jet equals zero
    #  (2) this requirement implies that there is at least one track within the tau-jet
    #     (leading to a small underestimation of backgrounds)
    #
    changeCut(process, "selectedMuTauPairsZeroCharge", "leg2.leadPFChargedHadrCand.isNonnull & (leg1.charge + leg2.leadPFChargedHadrCand.charge) = 0")

    # fill histograms only for events passing all event selection critera;
    # disable storing run and event numbers for events passing selection
    # (in order to save space in the .root files)
    if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and \
         hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
        pruneAnalysisSequence(process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation)
        pruneAnalysisSequence(process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation)
    else:
        pruneAnalysisSequence(process.analyzeZtoMuTauEvents)

    # disable event print-out for all analysis sequences
    if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and \
         hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
        disableEventDump(process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation)
        disableEventDump(process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation)
    else:
        disableEventDump(process.analyzeZtoMuTauEvents)

    # enable checking of fake-rates and tau id. efficiencies
    # with event weights in tau-jet histogram manager
    setattr(process.tauHistManager, "checkWeightConsistency", cms.bool(True))

    # get list of fake-rates types to be processed
    frTypes = getPSetAttributes(process.bgEstFakeRateJetWeights.frTypes)

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
            bgEstFakeRateAnalysisSequence = \
              addFakeRateGenAnalyzerModule(process, process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation,
                                           frType, bgEstFakeRateAnalysisSequence)            
            bgEstFakeRateAnalysisSequence = \
              addFakeRateGenAnalyzerModule(process, process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation,
                                           frType, bgEstFakeRateAnalysisSequence)
        else:
            bgEstFakeRateAnalysisSequence = \
              addFakeRateGenAnalyzerModule(process, process.analyzeZtoMuTauEvents,
                                           frType, bgEstFakeRateAnalysisSequence)

    # add analysis sequence:
    #  1.) with tau id. discriminators not applied
    #  2.) events **not** weighted by fake-rate
    # (for the purpose of making control plots for the data sample from which contributions 
    #  of individual background processes are estimated via the fake-rate technique)
    if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and 
         hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):        
        bgEstFakeRateAnalysisSequence = \
          addGenAnalyzerModule(process, process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation,
                               "frUnweighted", bgEstFakeRateAnalysisSequence)
        addAnalyzer(process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation_frUnweighted,
                    tauIdEffValidationHistManager, "evtSelTauLeadTrkPt",
                    "tauIdEffValidationHistManager.tauSource = selectedPatTausForMuTauLeadTrkPtCumulative")
        addAnalyzer(process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation_frUnweighted,
                    tauIdEffValidationHistManager, "evtSelDiMuPairZmumuHypothesisVeto",
                    "tauIdEffValidationHistManager.tauSource = selectedPatTausForMuTauMuonVetoCumulative")
        bgEstFakeRateAnalysisSequence = \
          addGenAnalyzerModule(process, process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation,
                               "frUnweighted", bgEstFakeRateAnalysisSequence)
        addAnalyzer(process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation_frUnweighted,
                    tauIdEffValidationHistManager, "evtSelTauLeadTrkPt",
                    "tauIdEffValidationHistManager.tauSource = selectedPatTausForMuTauLeadTrkPtCumulative")
        addAnalyzer(process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation_frUnweighted,
                    tauIdEffValidationHistManager, "evtSelDiMuPairZmumuHypothesisVeto",
                    "tauIdEffValidationHistManager.tauSource = selectedPatTausForMuTauMuonVetoCumulative")
    else:
        bgEstFakeRateAnalysisSequence = \
          addGenAnalyzerModule(process, process.analyzeZtoMuTauEvents,
                               "frUnweighted", bgEstFakeRateAnalysisSequence)
        addAnalyzer(process.analyzeZtoMuTauEvents_frUnweighted,
                    tauIdEffValidationHistManager, "evtSelTauLeadTrkPt",
                    "tauIdEffValidationHistManager.tauSource = selectedPatTausForMuTauLeadTrkPtCumulative")
        addAnalyzer(process.analyzeZtoMuTauEvents_frUnweighted,
                    tauIdEffValidationHistManager, "evtSelDiMuPairZmumuHypothesisVeto",
                    "tauIdEffValidationHistManager.tauSource = selectedPatTausForMuTauMuonVetoCumulative")

    # if method is "simple", add one more analysis sequence:
    #  1.) with tau id. discriminators not applied
    #  2.) events weighted by tau id. efficiency
    # (for the purpose of checking the tau id. efficiency values
    #  which are used by the "CDF" method)
    if method == "simple":
        
        tauIdEfficiency = cms.PSet(
            tauJetDiscriminators = cms.VPSet(
                cms.PSet(
                    tauJetIdEffSource = cms.InputTag("shrinkingConeZTTEffSimAssociator", "effByStandardChainZTTsim"),
                    qcdJetFakeRateSource = cms.InputTag("shrinkingConeZTTEffSimAssociator", "effByStandardChainZTTsim"),
                    tauJetDiscrSource = cms.InputTag("ewkTauId")
                )
            )
        )
        
        setattr(process.bgEstFakeRateJetWeights.frTypes, "tauIdEfficiency", tauIdEfficiency)
        setattr( process.bgEstFakeRateEventWeights.frTypes, "tauIdEfficiency", tauIdEfficiency)
        frLabel = "".join(["bgEstFakeRateJetWeight", "_", "tauIdEfficiency"])
        frInputTag = cms.InputTag('bgEstFakeRateJetWeights', "tauIdEfficiency")
        setattr(process.patTaus.efficiencies, frLabel, frInputTag)
        
        if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and 
             hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
            bgEstFakeRateAnalysisSequence = \
              addFakeRateGenAnalyzerModule(process, process.analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation,
                                           "tauIdEfficiency", bgEstFakeRateAnalysisSequence)
            bgEstFakeRateAnalysisSequence = \
              addFakeRateGenAnalyzerModule(process, process.analyzeZtoMuTauEvents_factorizedWithMuonIsolation,
                                           "tauIdEfficiency", bgEstFakeRateAnalysisSequence)
        else:
            bgEstFakeRateAnalysisSequence = \
              addFakeRateGenAnalyzerModule(process, process.analyzeZtoMuTauEvents,
                                           "tauIdEfficiency", bgEstFakeRateAnalysisSequence)

    setattr(process, "bgEstFakeRateAnalysisSequence", cms.Sequence(bgEstFakeRateAnalysisSequence))

    if ( hasattr(process, "analyzeZtoMuTauEvents_factorizedWithoutMuonIsolation") and
         hasattr(process, "analyzeZtoMuTauEvents_factorizedWithMuonIsolation") ):
        process.p.replace(process.analyzeZtoMuTauEvents_factorized, process.bgEstFakeRateAnalysisSequence)
    else:
        process.p.replace(process.analyzeZtoMuTauEvents, process.bgEstFakeRateAnalysisSequence)

def enableFakeRates_makeZtoMuTauPlots(process, enableFactorization = True):

    # get list of fake-rates types to be processed
    process.load("TauAnalysis.BgEstimationTools.fakeRateJetWeightProducer_cfi")
    frTypes = getPSetAttributes(process.bgEstFakeRateJetWeights.frTypes)
    frTypes.append("frUnweighted")

    seq_isFirstModule = True

    for frType in frTypes:

        mod_addZtoMuTau_qcdSum = copy.deepcopy(process.addBgEstFakeRateZtoMuTau_qcdSum_tauFakeRate)
        modInputDir_addZtoMuTau_qcdSum = cms.vstring(
            "".join(['tauFakeRate/harvested/InclusivePPmuX/zMuTauAnalyzer', '_fr_', frType]),
            "".join(['tauFakeRate/harvested/PPmuXptGt20/zMuTauAnalyzer', '_fr_', frType])
        )
        setattr(mod_addZtoMuTau_qcdSum.qcdSum, "dqmDirectories_input", modInputDir_addZtoMuTau_qcdSum)
        modOutputDir_addZtoMuTau_qcdSum = cms.string("".join(['tauFakeRate/harvested/qcdSum/zMuTauAnalyzer', '_fr_', frType]))
        setattr(mod_addZtoMuTau_qcdSum.qcdSum, "dqmDirectory_output", modOutputDir_addZtoMuTau_qcdSum)
        modName_addZtoMuTau_qcdSum = "".join(["addBgEstFakeRateZtoMuTau_qcdSum_tauFakeRate", "_", frType])
        setattr(process, modName_addZtoMuTau_qcdSum, mod_addZtoMuTau_qcdSum)

        seq_addZtoMuTau = cms.Sequence(getattr(process, modName_addZtoMuTau_qcdSum))

        modName_addZtoMuTau_smBgSum = "undefined"
        if hasattr(process, "addBgEstFakeRateZtoMuTau_smBgSum_tauFakeRate"):
            mod_addZtoMuTau_smBgSum = copy.deepcopy(process.addBgEstFakeRateZtoMuTau_smBgSum_tauFakeRate)
            modInputDir_addZtoMuTau_smBgSum = cms.vstring(
                "".join(['tauFakeRate/harvested/Zmumu/zMuTauAnalyzer', '_fr_', frType]),
                "".join(['tauFakeRate/harvested/WplusJets/zMuTauAnalyzer', '_fr_', frType]),
                "".join(['tauFakeRate/harvested/TTplusJets/zMuTauAnalyzer', '_fr_', frType]),
                "".join(['tauFakeRate/harvested/qcdSum/zMuTauAnalyzer', '_fr_', frType])
            )
            setattr(mod_addZtoMuTau_smBgSum.smBgSum, "dqmDirectories_input", modInputDir_addZtoMuTau_smBgSum)
            modOutputDir_addZtoMuTau_smBgSum = cms.string("".join(['tauFakeRate/harvested/smBgSum/zMuTauAnalyzer', '_fr_', frType]))
            setattr(mod_addZtoMuTau_smBgSum.smBgSum, "dqmDirectory_output", modOutputDir_addZtoMuTau_smBgSum)
            modName_addZtoMuTau_smBgSum = "".join(["addBgEstFakeRateZtoMuTau_smBgSum_tauFakeRate", "_", frType])
            setattr(process, modName_addZtoMuTau_smBgSum, mod_addZtoMuTau_smBgSum)

            seq_addZtoMuTau._seq = seq_addZtoMuTau._seq * getattr(process, modName_addZtoMuTau_smBgSum)
             
        modName_addZtoMuTau_smSum = "undefined"
        if hasattr(process, "addBgEstFakeRateZtoMuTau_smSum_tauFakeRate"):
            mod_addZtoMuTau_smSum = copy.deepcopy(process.addBgEstFakeRateZtoMuTau_smSum_tauFakeRate)
            modInputDir_addZtoMuTau_smSum = cms.vstring(
                "".join(['tauFakeRate/harvested/Ztautau/zMuTauAnalyzer', '_fr_', frType]),
                "".join(['tauFakeRate/harvested/smBgSum/zMuTauAnalyzer', '_fr_', frType])
            )
            setattr(mod_addZtoMuTau_smSum.smSum, "dqmDirectories_input", modInputDir_addZtoMuTau_smSum)
            modOutputDir_addZtoMuTau_smSum = cms.string("".join(['tauFakeRate/harvested/smSum/zMuTauAnalyzer', '_fr_', frType]))
            setattr(mod_addZtoMuTau_smSum.smSum, "dqmDirectory_output", modOutputDir_addZtoMuTau_smSum)
            modName_addZtoMuTau_smSum = "".join(["addBgEstFakeRateZtoMuTau_smSum_tauFakeRate", "_", frType])
            setattr(process, modName_addZtoMuTau_smSum, mod_addZtoMuTau_smSum)

            seq_addZtoMuTau._seq = seq_addZtoMuTau._seq * getattr(process, modName_addZtoMuTau_smSum)
            
        seqName_addZtoMuTau = "".join(["addBgEstFakeRateZtoMuTau_tauFakeRate", "_", frType])
        setattr(process, seqName_addZtoMuTau, seq_addZtoMuTau)

        ##if enableFactorization:
        ##    enableFactorization_makeZtoMuTauPlots_grid(process,
        ##      dqmDirectoryIn_InclusivePPmuX = \
        ##        "".join(['tauFakeRate/harvested/InclusivePPmuX/zMuTauAnalyzer', '_fr_', frType]),
        ##      dqmDirectoryOut_InclusivePPmuX = \
        ##        "".join(['tauFakeRate/harvested/InclusivePPmuX_factorized/zMuTauAnalyzer', '_fr_', frType]),
        ##      dqmDirectoryIn_PPmuXptGt20 = \
        ##        "".join(['tauFakeRate/harvested/PPmuXptGt20/zMuTauAnalyzer', '_fr_', frType]),
        ##      dqmDirectoryOut_PPmuXptGt20 = \
        ##        "".join(['tauFakeRate/harvested/PPmuXptGt20_factorized/zMuTauAnalyzer', '_fr_', frType]),
        ##      modName_addZtoMuTau_qcdSum = modName_addZtoMuTau_qcdSum,
        ##      modName_addZtoMuTau_smBgSum = modName_addZtoMuTau_smBgSum,         
        ##      modName_addZtoMuTau_smSum = modName_addZtoMuTau_smSum,
        ##      seqName_addZtoMuTau = seqName_addZtoMuTau,
        ##      pyObjectLabel = frType)

        if seq_isFirstModule:
            setattr(process, "addBgEstFakeRateZtoMuTau_tauFakeRate", cms.Sequence(getattr(process, seqName_addZtoMuTau)))
            seq_isFirstModule = False
        else:
            process.addBgEstFakeRateZtoMuTau_tauFakeRate._seq = \
              process.addBgEstFakeRateZtoMuTau_tauFakeRate._seq * getattr(process, seqName_addZtoMuTau)
        
#--------------------------------------------------------------------------------
# utility functions specific to application of fake-rate weights
# to tau id. efficiency measurement analysis
#--------------------------------------------------------------------------------

def enableFakeRates_runTauIdEffAnalysisZtoMuTau(process):

    method = "simple"

    # compute fake-rate weights
    configureFakeRateWeightProduction(process, method = method)

    # enable checking of fake-rates and tau id. efficiencies
    # with event weights in tau-jet histogram manager
    setattr(process.tauHistManager, "checkWeightConsistency", cms.bool(True))

    # get list of fake-rates types to be processed
    frTypes = getPSetAttributes(process.bgEstFakeRateJetWeights.frTypes)

    fakeRateAnalysisSequence = None  

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
        fakeRateAnalysisSequence = \
          addFakeRateGenAnalyzerModule(process, process.analyzeEventsTauIdEffZtoMuTauCombinedFit,
                                       frType, fakeRateAnalysisSequence)

    setattr(process, "fakeRateAnalysisSequence", cms.Sequence(fakeRateAnalysisSequence))

    process.bgEstTauIdEffZtoMuTauCombinedFitAnalysisSequence._seq = \
      process.bgEstTauIdEffZtoMuTauCombinedFitAnalysisSequence._seq * process.fakeRateAnalysisSequence

def enableFakeRates_makeTauIdEffZtoMuTauPlots(process):

    # get list of fake-rates types to be processed
    process.load("TauAnalysis.BgEstimationTools.fakeRateJetWeightProducer_cfi")
    frTypes = getPSetAttributes(process.bgEstFakeRateJetWeights.frTypes)

    for frType in frTypes:

        mod_addZtoMuTau_qcdSum = copy.deepcopy(process.addTauIdEffZtoMuTau_qcdSum)
        modInputDir_addZtoMuTau_qcdSum = cms.vstring(
            "".join(['harvested/InclusivePPmuX/TauIdEffAnalyzerZtoMuTauCombinedFit', '_fr_', frType]),
            "".join(['harvested/PPmuXptGt20/TauIdEffAnalyzerZtoMuTauCombinedFit', '_fr_', frType])
        )
        setattr(mod_addZtoMuTau_qcdSum.qcdSum, "dqmDirectories_input", modInputDir_addZtoMuTau_qcdSum)
        modOutputDir_addZtoMuTau_qcdSum = \
          cms.string("".join(['tauFakeRate/harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFit', '_fr_', frType]))
        setattr(mod_addZtoMuTau_qcdSum.qcdSum, "dqmDirectory_output", modOutputDir_addZtoMuTau_qcdSum)
        modName_addZtoMuTau_qcdSum = "".join(["addTauIdEffZtoMuTau_qcdSum", "_", frType])
        setattr(process, modName_addZtoMuTau_qcdSum, mod_addZtoMuTau_qcdSum)
        
        process.addTauIdEffZtoMuTau._seq = process.addTauIdEffZtoMuTau._seq * mod_addZtoMuTau_qcdSum
             
        mod_addZtoMuTau_smSum = copy.deepcopy(process.addTauIdEffZtoMuTau_smSum)
        modInputDir_addZtoMuTau_smSum = cms.vstring(
            "".join(['harvested/Ztautau/TauIdEffAnalyzerZtoMuTauCombinedFit', '_fr_', frType]),
            "".join(['harvested/Zmumu/TauIdEffAnalyzerZtoMuTauCombinedFit', '_fr_', frType]),
            "".join(['harvested/WplusJets/TauIdEffAnalyzerZtoMuTauCombinedFit', '_fr_', frType]),
            "".join(['harvested/TTplusJets/TauIdEffAnalyzerZtoMuTauCombinedFit', '_fr_', frType]),
            "".join(['harvested/qcdSum/TauIdEffAnalyzerZtoMuTauCombinedFit', '_fr_', frType])
        )
        setattr(mod_addZtoMuTau_smSum.smSum, "dqmDirectories_input", modInputDir_addZtoMuTau_smSum)
        modOutputDir_addZtoMuTau_smSum = \
          cms.string("".join(['tauFakeRate/harvested/smSum/TauIdEffAnalyzerZtoMuTauCombinedFit', '_fr_', frType]))
        setattr(mod_addZtoMuTau_smSum.smSum, "dqmDirectory_output", modOutputDir_addZtoMuTau_smSum)
        modName_addZtoMuTau_smSum = "".join(["addTauIdEffZtoMuTau_smSum", "_", frType])
        setattr(process, modName_addZtoMuTau_smSum, mod_addZtoMuTau_smSum)
        
        process.addTauIdEffZtoMuTau._seq = process.addTauIdEffZtoMuTau._seq * mod_addZtoMuTau_smSum
