import FWCore.ParameterSet.Config as cms

tauIdEffZtoMuTauHistManager = cms.PSet(    
    pluginName = cms.string('tauIdEffZtoMuTauHistManager'),
    pluginType = cms.string('TauIdEffZtoMuTauHistManager'),

    muonSource = cms.InputTag('selectedLayer1MuonsTrkIPcumulative'),  
    tauSource = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative'),
    diTauSource = cms.InputTag('muTauPairsTauIdEffZtoMuTauValidCollinearApproxAbsMuonIsolation'),
    centralJetSource = cms.InputTag('selectedLayer1JetsAntiOverlapWithLeptonsVetoCumulative'),

    dqmDirectory_store = cms.string('TauIdEffSpecificQuantities'),

    tauIdDiscriminator = cms.string("ewkTauId"),

    diTauChargeSignExtractor = cms.PSet(
        pluginType = cms.string("PATMuTauPairChargeSignExtractor"),
        src = cms.InputTag('muTauPairsTauIdEffZtoMuTauValidCollinearApproxAbsMuonIsolation')
    )
)
