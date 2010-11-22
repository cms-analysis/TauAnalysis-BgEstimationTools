import FWCore.ParameterSet.Config as cms

tauIdEffValidationHistManager = cms.PSet(
    pluginName = cms.string('tauIdEffValidationHistManager'),
    pluginType = cms.string('TauIdEffValidationHistManager'),

    ##tauSource = cms.InputTag('selectedPatTausLeadTrkPtCumulative'),
    tauSource = cms.InputTag('selectedPatTausForMuTauLeadTrkPtCumulative'),

    dqmDirectory_store = cms.string('TauIdEffValidation'),

    frTypes = cms.vstring(
        #"WJetssim",
        #"MuEnrichedQCDsim",
        #"DiJetHighPtsim",
        #"DiJetSecondPtsim"
        "DiJetHighPtdata",
        "DiJetSecondPtdata"
    ),

    effTypes = cms.vstring(
        "ZTTsim"
    ),

    # define numerator histograms to be filled by:
    #  o applying tau id. discriminators/cuts
    #  o not applying tau id. discriminators/cuts,
    #    but weighting tau-jet candidates by efficiency/fake-rate values instead
    # (NOTE: denominator histograms are automatically added by TauIdEffValidationHistManager)
    numerators = cms.VPSet(
        cms.PSet(
            cutEffName = cms.string("ByEWKTauID"),
            cuts = cms.vstring(
                "tauID('byTaNCfrQuarterPercent') > 0.5",
                "abs(charge) == 1",
            )
        ),
        #cms.PSet(
            #cutEffName = cms.string("ByTrackIsolationSeq"),
            #cuts = cms.vstring(
                #"tauID('trackIsolation') > 0.5"
            #)
        #),
        #cms.PSet(
            #cutEffName = cms.string("ByEcalIsolationSeq"),
            #cuts = cms.vstring(
                #"tauID('trackIsolation') > 0.5",
                #"tauID('ecalIsolation') > 0.5"
            #)
        #),
        #cms.PSet(
            #cutEffName = cms.string("ByNTracksSeq"),
            #cuts = cms.vstring(
                #"tauID('trackIsolation') > 0.5",
                #"tauID('ecalIsolation') > 0.5",
                #"signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3"
            #)
        #),
        #cms.PSet(
            #cutEffName = cms.string("ByChargeSeq"),
            #cuts = cms.vstring(
                #"tauID('trackIsolation') > 0.5",
                #"tauID('ecalIsolation') > 0.5",
                #"signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3",
                #"abs(charge) > 0.5 & abs(charge) < 1.5"
            #)
        #),
        #cms.PSet(
            #cutEffName = cms.string("ByStandardChain"),
            #cuts = cms.vstring(
                #"tauID('trackIsolation') > 0.5",
                #"tauID('ecalIsolation') > 0.5",
                #"signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3",
                #"abs(charge) > 0.5 & abs(charge) < 1.5"
            #)
        #)
    )
)
