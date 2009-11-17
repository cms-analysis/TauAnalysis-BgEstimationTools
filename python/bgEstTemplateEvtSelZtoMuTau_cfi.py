import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
# define event selection of background enriched samples
# from which template histograms are obtained
#--------------------------------------------------------------------------------

bgEstEventSelection_Zmumu = (
    "numDiTausZmumu >= 1 && muonTrackIsoZmumu_0 < 1. && muonEcalIsoZmumu_0 < 1."
#   " && tauDiscrAgainstMuonsZmumu_0 > 0.5 && numGlobalMuons >= 2 && numJetsAlpha0point1Zmumu < 1"
    " && (tauDiscrAgainstMuonsZmumu_0 < 0.5 || numGlobalMuons >= 2)"
#   " && diTauAbsChargeZmumu_0 < 0.5"
)

bgEstEventSelection_WplusJets = (
#   "numDiTausWplusJets >= 1 && muonPtWplusJets_0 > 25. && muonTrackIsoWplusJets_0 < 1. && muonEcalIsoWplusJets_0 < 1."
#
# Note: muonPt cut improves WplusJets/QCD ratio by about a factor five,
#       but significantly shifts the muon + tau-jet visible invariant mass distribution towards higher values.
#       In order to supress QCD background contamination (on a statistical basis),
#       could extract W + jets template shape from difference in muon + tau-jet visible invariant mass distributions
#       of opposite sign - same sign muon and tau-jet combinations.
#      (SS/OS ratio is close to one for QCD background; significant charge asymmetry expected for W + jets background)
#    
    "numDiTausWplusJets >= 1 && muonTrackIsoWplusJets_0 < 1. && muonEcalIsoWplusJets_0 < 1."
#   " && tauTrackIsoDiscrWplusJets_0 < 0.5 && tauTrackIsoWplusJets_0 > 2. && tauDiscrAgainstMuonsWplusJets_0 > 0.5"
#
# Note: probability for quark/gluon jets to pass tau track and ECAL isolation criteria
#       is higher for low Pt than for high Pt jets; the consequence is that muon + tau-jet visible invariant mass distribution
#       gets shifted towards higher values in case tau track and ECAL isolation criteria are not applied.
#       For this reason, either need to apply tau track and ECAL isolation criteria in selection of W + jets background enriched sample
#       or correct for template shape distortion by reweighting
#      (would gain a factor of about 2.5 in event statistics; reweighting of tauPt distribution not implemented yet, however)
#    
    " && tauTrackIsoDiscrWplusJets_0 > 0.5 && tauEcalIsoDiscrWplusJets_0 > 0.5 && tauDiscrAgainstMuonsWplusJets_0 > 0.5"
    " && diTauMt1MEtWplusJets_0 > 30."
    " && numGlobalMuons < 2"
    " && numJetsAlpha0point1WplusJets < 1"
)

bgEstEventSelection_TTplusJets = (
    "numDiTausTTplusJets >= 1 && muonTrackIsoTTplusJets_0 < 2. && muonEcalIsoTTplusJets_0 < 2."
    " && tauDiscrAgainstMuonsTTplusJets_0 > 0.5"
    " && diTauAbsChargeTTplusJets_0 < 0.5"
    " && numGlobalMuons < 2"
    " && ((numJetsEt40TTplusJets >= 1 && jetEt40bTaggingDiscrTrackCountingHighEffTTplusJets_0 > 4.5)"
    " || (numJetsEt40TTplusJets >= 2 && jetEt40bTaggingDiscrTrackCountingHighEffTTplusJets_1 > 4.5)"
    " || (numJetsEt40TTplusJets >= 3 && jetEt40bTaggingDiscrTrackCountingHighEffTTplusJets_2 > 4.5))"
    " && numJetsEt40TTplusJets >= 2 && numJetsEt60TTplusJets >= 1"
)

bgEstEventSelection_QCD = (
    "numDiTausQCD >= 1 && muonTrackIsoQCD_0 > 4. && muonEcalIsoQCD_0 > 4."
    " && tauDiscrAgainstMuonsQCD_0 > 0.5"
    " && numGlobalMuons < 2"
)

