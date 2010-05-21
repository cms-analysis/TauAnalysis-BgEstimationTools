import FWCore.ParameterSet.Config as cms

tauIdEffZtoMuTauHistManager = cms.PSet(    
    pluginName = cms.string('tauIdEffZtoMuTauHistManager'),
    pluginType = cms.string('TauIdEffZtoMuTauHistManager'),

    muonSource = cms.InputTag('selectedPatMuonsTrkIPcumulative'),  
    tauSource = cms.InputTag('tausForTauIdEffZtoMuTauMuonVetoCumulative'),
    diTauSource = cms.InputTag('muTauPairsTauIdEffZtoMuTauValidCollinearApproxAbsMuonIsolation'),
    centralJetSource = cms.InputTag('selectedPatJetsAntiOverlapWithLeptonsVetoCumulative'),

    dqmDirectory_store = cms.string('TauIdEffSpecificQuantities'),

    tauIdDiscriminator = cms.string("ewkTauId"),

    diTauChargeSignExtractor = cms.PSet(
        pluginType = cms.string("PATMuTauPairChargeSignExtractor"),
        src = cms.InputTag('muTauPairsTauIdEffZtoMuTauValidCollinearApproxAbsMuonIsolation')
    )
)
