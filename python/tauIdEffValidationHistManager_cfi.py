import FWCore.ParameterSet.Config as cms

tauIdEffValidationHistManager = cms.PSet(    
    pluginName = cms.string('tauIdEffValidationHistManager'),
    pluginType = cms.string('TauIdEffValidationHistManager'),

    tauSource = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative'),

    dqmDirectory_store = cms.string('TauIdEffValidation'),

    frTypes = cms.vstring(
        "WJetssim",
        "MuEnrichedQCDsim",
        "DiJetHighPtSim",
        "DiJetSecondPtsim"
    ),

    effTypes = cms.vstring(
        "ZTTsim"
    ),

    frValidation = cms.VPSet(
        cms.PSet(
            cutEffName = cms.string("ByStandardChain"),
            cuts = cms.vstring(
                "tauID('trackIsolation') > 0.5",
                "tauID('ecalIsolation') > 0.5",
                "signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3",
                "abs(charge) > 0.5 & abs(charge) < 1.5"
            )
        )
    )
)
