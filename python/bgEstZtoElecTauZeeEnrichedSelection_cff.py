import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

#--------------------------------------------------------------------------------
# select Z --> e+ e- background enriched event sample
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patLeptonSelection_cff import *

#--------------------------------------------------------------------------------  
# produce collection of pat::Electrons
#--------------------------------------------------------------------------------

# require electron candidate to pass the eidRobustTight electron id. criteria
electronsBgEstZeeEnrichedId = copy.deepcopy(selectedPatElectronsForElecTauId)

# require electron candidate to not be within eta-crack
# between Barrel and Encap ECAL calorimeter
electronsBgEstZeeEnrichedAntiCrackCut = copy.deepcopy(selectedPatElectronsForElecTauAntiCrackCut)

# require electron candidate to be within geometric acceptance of electron trigger
electronsBgEstZeeEnrichedEta = copy.deepcopy(selectedPatElectronsForElecTauEta21)

# require electron candidate to have transverse momentum above threshold
electronsBgEstZeeEnrichedPt = copy.deepcopy(selectedPatElectronsForElecTauPt15)
electronsBgEstZeeEnrichedPt.cut = cms.string('pt > 30.')

# require electron candidate to be isolated
# with respect to tracks (of Pt >~ 0.3 GeV)
electronsBgEstZeeEnrichedTrkIso = copy.deepcopy(selectedPatElectronsForElecTauTrkIso)

# require electron candidate to be isolated
# with respect to energy deposits in ECAL
# (not associated to electron candidate)
electronsBgEstZeeEnrichedEcalIso = copy.deepcopy(selectedPatElectronsForElecTauEcalIso)

# require electron candidate to be linked to (GSF) track
#electronsBgEstZeeEnrichedTrk = copy.deepcopy(selectedPatElectronsTrk)

# require track of electron candidate to have small transverse impact parameter
# (in order to veto electrons resulting from b-quark decays)
#electronsBgEstZeeEnrichedTrkIP = copy.deepcopy(selectedPatElectronsForElecTauTrkIP)

# require electron to not be from a photon conversion
electronsBgEstZeeEnrichedConversionVeto = copy.deepcopy(selectedPatElectronsForElecTauConversionVeto)


electronSelConfiguratorBgEstZeeEnriched = objSelConfigurator(
    [ electronsBgEstZeeEnrichedId,
      electronsBgEstZeeEnrichedAntiCrackCut,
      electronsBgEstZeeEnrichedEta,
      electronsBgEstZeeEnrichedPt,
      electronsBgEstZeeEnrichedTrkIso,
      electronsBgEstZeeEnrichedEcalIso,
#      electronsBgEstZeeEnrichedTrkIP,
      electronsBgEstZeeEnrichedConversionVeto],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = False
)


selectElectronsBgEstZeeEnriched = electronSelConfiguratorBgEstZeeEnriched.configure(pyNameSpace = locals())

   
#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

#selectTausBgEstZeeEnriched = copy.deepcopy(selectPatTausForElecTau)

# require tau candidate not to overlap with selected electrons
# (in order to avoid double-counting one and the same physical particle
#  as electron and as tau candidate)
tausBgEstZeeEnrichedAntiOverlapWithElectronsVeto = cms.EDFilter("PATTauAntiOverlapSelector",
    srcNotToBeFiltered = cms.VInputTag("electronsBgEstZeeEnrichedPtCumulative"),
    dRmin = cms.double(0.3),
    filter = cms.bool(False)                                           
)

# require tau candidate to be within geometric acceptance of Pixel + SiTracker detectors
tausBgEstZeeEnrichedEta = copy.deepcopy(selectedPatTausForElecTauEta21)

# require tau candidate to have transverse energy above threshold
tausBgEstZeeEnrichedPt = copy.deepcopy(selectedPatTausForElecTauPt20)

# require tau candidate to have a leading track
# (track of Pt > 1. GeV within matching cone of size dR = 0.2 around jet-axis)
tausBgEstZeeEnrichedLeadTrk = copy.deepcopy(selectedPatTausForElecTauLeadTrk)

# require leading track of tau candidate to have Pt > 5. GeV
tausBgEstZeeEnrichedLeadTrkPt = copy.deepcopy(selectedPatTausForElecTauLeadTrkPt)

# require tau candidate to pass TaNC discriminator
tausBgEstZeeEnrichedTaNCdiscr = copy.deepcopy(selectedPatTausForElecTauTaNCdiscr)

# require tau candidate to have no tracks of Pt > 1. GeV
# in isolation cone of size dR = 0.8, surrounding signal cone of size dR = 5./Et
tausBgEstZeeEnrichedTrkIso = copy.deepcopy(selectedPatTausForElecTauTrkIso)

# require tau candidate to be isolated
# with respect to energy deposits in ECAL
tausBgEstZeeEnrichedEcalIso = copy.deepcopy(selectedPatTausForElecTauEcalIso)

# require tau candidate to have either one or three tracks within signal cone
tausBgEstZeeEnrichedProng = copy.deepcopy(selectedPatTausForElecTauProng)

# require tau candidate to have charge either +1 or -1
# (computed as sum of charges of tracks within signal cone)
tausBgEstZeeEnrichedCharge = copy.deepcopy(selectedPatTausForElecTauCharge)

# require tau candidate to pass electron veto
tausBgEstZeeEnrichedElectronVeto = copy.deepcopy(selectedPatTausForElecTauElectronVeto)

# require tau candidate not to be in ECAL barrel/endcap crack
tausBgEstZeeEnrichedEcalCrackVeto = copy.deepcopy(selectedPatTausForElecTauEcalCrackVeto)

# require tau candidate to pass muon veto
tausBgEstZeeEnrichedMuonVeto = copy.deepcopy(selectedPatTausForElecTauMuonVeto)



tauSelConfiguratorBgEstZeeEnriched = objSelConfigurator(
    [ tausBgEstZeeEnrichedAntiOverlapWithElectronsVeto,
      tausBgEstZeeEnrichedEta,
      tausBgEstZeeEnrichedPt,
      tausBgEstZeeEnrichedLeadTrk,
      tausBgEstZeeEnrichedLeadTrkPt,
      tausBgEstZeeEnrichedTaNCdiscr,
      tausBgEstZeeEnrichedTrkIso,
      tausBgEstZeeEnrichedEcalIso,
      tausBgEstZeeEnrichedProng,
      tausBgEstZeeEnrichedCharge,
      tausBgEstZeeEnrichedElectronVeto,
      tausBgEstZeeEnrichedEcalCrackVeto,
      tausBgEstZeeEnrichedMuonVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectTausBgEstZeeEnriched = tauSelConfiguratorBgEstZeeEnriched.configure(pyNameSpace = locals())


#--------------------------------------------------------------------------------  
# produce collection of electron + tau-jet combinations
#--------------------------------------------------------------------------------


from TauAnalysis.CandidateTools.elecTauPairProduction_cff import *

elecTauPairsBgEstZeeEnriched = copy.deepcopy(allElecTauPairs)
elecTauPairsBgEstZeeEnriched.srcLeg1 = cms.InputTag('electronsBgEstZeeEnrichedConversionVetoCumulative')
elecTauPairsBgEstZeeEnriched.srcLeg2 = cms.InputTag('tausBgEstZeeEnrichedMuonVetoCumulative')

produceElecTauPairsBgEstZeeEnriched = cms.Sequence(elecTauPairsBgEstZeeEnriched)


from TauAnalysis.CandidateTools.elecTauPairSelection_cfi import *

elecTauPairsBgEstZeeEnrichedAntiOverlapVeto = copy.deepcopy(selectedElecTauPairsAntiOverlapVeto)
elecTauPairsBgEstZeeEnrichedZeroCharge = copy.deepcopy(selectedElecTauPairsZeroCharge)
elecTauPairsBgEstZeeEnrichedAcoplanarity12 = copy.deepcopy(selectedElecTauPairsAcoplanarity12)
elecTauPairsBgEstZeeEnrichedMt1MET = copy.deepcopy(selectedElecTauPairsMt1MET)
elecTauPairsBgEstZeeEnrichedPzetaDiff = copy.deepcopy(selectedElecTauPairsPzetaDiff)



elecTauPairSelConfiguratorBgEstZeeEnriched = objSelConfigurator(
    [ elecTauPairsBgEstZeeEnrichedAntiOverlapVeto,
      elecTauPairsBgEstZeeEnrichedZeroCharge,
      elecTauPairsBgEstZeeEnrichedAcoplanarity12,
      #elecTauPairsBgEstZeeEnrichedMt1MET,
      elecTauPairsBgEstZeeEnrichedPzetaDiff ],
    src = "elecTauPairsBgEstZeeEnriched",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectElecTauPairsBgEstZeeEnriched = elecTauPairSelConfiguratorBgEstZeeEnriched.configure(pyNameSpace = locals())



#--------------------------------------------------------------------------------  
# produce boolean event selection flags
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.selectZtoElecTau_cff import *

# electron cuts

cfgElectronIdCutBgEstZeeEnriched = copy.deepcopy(cfgElectronIdCut)
cfgElectronIdCutBgEstZeeEnriched.pluginName = cms.string('electronIdCutBgEstZeeEnriched')
cfgElectronIdCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedIdCumulative')
cfgElectronIdCutBgEstZeeEnriched.systematics = cms.vstring()

cfgElectronAntiCrackCutBgEstZeeEnriched = copy.deepcopy(cfgElectronAntiCrackCut)
cfgElectronAntiCrackCutBgEstZeeEnriched.pluginName = cms.string('electronAntiCrackCutBgEstZeeEnriched')
cfgElectronAntiCrackCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedAntiCrackCutCumulative')
cfgElectronAntiCrackCutBgEstZeeEnriched.systematics = cms.vstring()

cfgElectronEtaCutBgEstZeeEnriched = copy.deepcopy(cfgElectronEtaCut)
cfgElectronEtaCutBgEstZeeEnriched.pluginName = cms.string('electronEtaCutBgEstZeeEnriched')
cfgElectronEtaCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedEtaCumulative')
cfgElectronEtaCutBgEstZeeEnriched.systematics = cms.vstring()

cfgElectronPtCutBgEstZeeEnriched = copy.deepcopy(cfgElectronPtCut)
cfgElectronPtCutBgEstZeeEnriched.pluginName = cms.string('electronPtCutBgEstZeeEnriched')
cfgElectronPtCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedPtCumulative')
cfgElectronPtCutBgEstZeeEnriched.systematics = cms.vstring()

cfgElectronTrkIsoCutBgEstZeeEnriched = copy.deepcopy(cfgElectronTrkIsoCut)
cfgElectronTrkIsoCutBgEstZeeEnriched.pluginName = cms.string('electronTrkIsoCutBgEstZeeEnriched')
cfgElectronTrkIsoCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedTrkIsoCumulative')
cfgElectronTrkIsoCutBgEstZeeEnriched.systematics = cms.vstring()

cfgElectronEcalIsoCutBgEstZeeEnriched = copy.deepcopy(cfgElectronEcalIsoCut)
cfgElectronEcalIsoCutBgEstZeeEnriched.pluginName = cms.string('electronEcalIsoCutBgEstZeeEnriched')
cfgElectronEcalIsoCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedEcalIsoCumulative')
cfgElectronEcalIsoCutBgEstZeeEnriched.systematics = cms.vstring()

## cfgElectronTrkIPcutBgEstZeeEnriched = copy.deepcopy(cfgElectronTrkIPcut)
## cfgElectronTrkIPcutBgEstZeeEnriched.pluginName = cms.string('electronTrkIPcutBgEstZeeEnriched')
## cfgElectronTrkIPcutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedTrkIPCumulative')
## cfgElectronTrkIPcutBgEstZeeEnriched.systematics = cms.vstring()

cfgElectronConversionVetoBgEstZeeEnriched = copy.deepcopy(cfgElectronConversionVeto)
cfgElectronConversionVetoBgEstZeeEnriched.pluginName = cms.string('electronConversionVetoBgEstZeeEnriched')
cfgElectronConversionVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedConversionVetoCumulative')
cfgElectronConversionVetoBgEstZeeEnriched.systematics = cms.vstring()

# tau cuts

cfgTauAntiOverlapWithElectronsVetoBgEstZeeEnriched = copy.deepcopy(cfgTauAntiOverlapWithElectronsVeto)
cfgTauAntiOverlapWithElectronsVetoBgEstZeeEnriched.pluginName = cms.string('tauAntiOverlapWithElectronsVetoBgEstZeeEnriched')
cfgTauAntiOverlapWithElectronsVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedAntiOverlapWithElectronsVetoCumulative')
cfgTauAntiOverlapWithElectronsVetoBgEstZeeEnriched.systematics = cms.vstring()

cfgTauEtaCutBgEstZeeEnriched = copy.deepcopy(cfgTauEtaCut)
cfgTauEtaCutBgEstZeeEnriched.pluginName = cms.string('tauEtaCutBgEstZeeEnriched')
cfgTauEtaCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedEtaCumulative')
cfgTauEtaCutBgEstZeeEnriched.systematics = cms.vstring()

cfgTauPtCutBgEstZeeEnriched = copy.deepcopy(cfgTauPtCut)
cfgTauPtCutBgEstZeeEnriched.pluginName = cms.string('tauPtCutBgEstZeeEnriched')
cfgTauPtCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedPtCumulative')
cfgTauPtCutBgEstZeeEnriched.systematics = cms.vstring()

cfgTauLeadTrkCutBgEstZeeEnriched = copy.deepcopy(cfgTauLeadTrkCut)
cfgTauLeadTrkCutBgEstZeeEnriched.pluginName = cms.string('tauLeadTrkCutBgEstZeeEnriched')
cfgTauLeadTrkCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedLeadTrkCumulative')
cfgTauLeadTrkCutBgEstZeeEnriched.systematics = cms.vstring()

cfgTauLeadTrkPtCutBgEstZeeEnriched = copy.deepcopy(cfgTauLeadTrkPtCut)
cfgTauLeadTrkPtCutBgEstZeeEnriched.pluginName = cms.string('tauLeadTrkPtCutBgEstZeeEnriched')
cfgTauLeadTrkPtCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedLeadTrkPtCumulative')
cfgTauLeadTrkPtCutBgEstZeeEnriched.systematics = cms.vstring()

cfgTauTaNCdiscrCutBgEstZeeEnriched = copy.deepcopy(cfgTauTaNCdiscrCut)
cfgTauTaNCdiscrCutBgEstZeeEnriched.pluginName = cms.string('tauTaNCdiscrCutBgEstZeeEnriched')
cfgTauTaNCdiscrCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedTaNCdiscrCumulative')
cfgTauTaNCdiscrCutBgEstZeeEnriched.systematics = cms.vstring()

cfgTauTrkIsoCutBgEstZeeEnriched = copy.deepcopy(cfgTauTrkIsoCut)
cfgTauTrkIsoCutBgEstZeeEnriched.pluginName = cms.string('tauTrkIsoCutBgEstZeeEnriched')
cfgTauTrkIsoCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedTrkIsoCumulative')
cfgTauTrkIsoCutBgEstZeeEnriched.systematics = cms.vstring()

cfgTauEcalIsoCutBgEstZeeEnriched = copy.deepcopy(cfgTauEcalIsoCut)
cfgTauEcalIsoCutBgEstZeeEnriched.pluginName = cms.string('tauEcalIsoCutBgEstZeeEnriched')
cfgTauEcalIsoCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedEcalIsoCumulative')
cfgTauEcalIsoCutBgEstZeeEnriched.systematics = cms.vstring()

cfgTauProngCutBgEstZeeEnriched = copy.deepcopy(cfgTauProngCut)
cfgTauProngCutBgEstZeeEnriched.pluginName = cms.string('tauProngCutBgEstZeeEnriched')
cfgTauProngCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedProngCumulative')
cfgTauProngCutBgEstZeeEnriched.systematics = cms.vstring()

cfgTauChargeCutBgEstZeeEnriched = copy.deepcopy(cfgTauChargeCut)
cfgTauChargeCutBgEstZeeEnriched.pluginName = cms.string('tauChargeCutBgEstZeeEnriched')
cfgTauChargeCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedChargeCumulative')
cfgTauChargeCutBgEstZeeEnriched.systematics = cms.vstring()

cfgTauElectronVetoBgEstZeeEnriched = copy.deepcopy(cfgTauElectronVeto)
cfgTauElectronVetoBgEstZeeEnriched.pluginName = cms.string('tauElectronVetoBgEstZeeEnriched')
cfgTauElectronVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedElectronVetoCumulative')
cfgTauElectronVetoBgEstZeeEnriched.systematics = cms.vstring()

cfgTauEcalCrackVetoBgEstZeeEnriched = copy.deepcopy(cfgTauEcalCrackVeto)
cfgTauEcalCrackVetoBgEstZeeEnriched.pluginName = cms.string('tauEcalCrackVetoBgEstZeeEnriched')
cfgTauEcalCrackVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedEcalCrackVetoCumulative')
cfgTauEcalCrackVetoBgEstZeeEnriched.systematics = cms.vstring()

cfgTauMuonVetoBgEstZeeEnriched = copy.deepcopy(cfgTauMuonVeto)
cfgTauMuonVetoBgEstZeeEnriched.pluginName = cms.string('tauMuonVetoBgEstZeeEnriched')
cfgTauMuonVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedMuonVetoCumulative')
cfgTauMuonVetoBgEstZeeEnriched.systematics = cms.vstring()

# di-tau candidate selection

cfgDiTauCandidateAntiOverlapVetoBgEstZeeEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauAntiOverlapVeto)
cfgDiTauCandidateAntiOverlapVetoBgEstZeeEnriched.pluginName = cms.string('diTauCandidateAntiOverlapVetoBgEstZeeEnriched')
cfgDiTauCandidateAntiOverlapVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZeeEnrichedAntiOverlapVetoCumulative')
cfgDiTauCandidateAntiOverlapVetoBgEstZeeEnriched.systematics = cms.vstring()

cfgDiTauCandidateForElecTauZeroChargeCutBgEstZeeEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauZeroChargeCut)
cfgDiTauCandidateForElecTauZeroChargeCutBgEstZeeEnriched.pluginName = cms.string('diTauCandidateZeroChargeCutBgEstZeeEnriched')
cfgDiTauCandidateForElecTauZeroChargeCutBgEstZeeEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZeeEnrichedZeroChargeCumulative')
cfgDiTauCandidateForElecTauZeroChargeCutBgEstZeeEnriched.systematics = cms.vstring()

cfgDiTauCandidateForElecTauAcoplanarity12CutBgEstZeeEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauAcoplanarity12Cut)
cfgDiTauCandidateForElecTauAcoplanarity12CutBgEstZeeEnriched.pluginName = cms.string('diTauCandidateAcoplanarity12CutBgEstZeeEnriched')
cfgDiTauCandidateForElecTauAcoplanarity12CutBgEstZeeEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZeeEnrichedAcoplanarity12Cumulative')
cfgDiTauCandidateForElecTauAcoplanarity12CutBgEstZeeEnriched.systematics = cms.vstring()

## cfgDiTauCandidateForElecTauMt1METCutBgEstZeeEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauMt1METCut)
## cfgDiTauCandidateForElecTauMt1METCutBgEstZeeEnriched.pluginName = cms.string('diTauCandidateMt1METCutBgEstZeeEnriched')
## cfgDiTauCandidateForElecTauMt1METCutBgEstZeeEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZeeEnrichedMt1METCumulative')
## cfgDiTauCandidateForElecTauMt1METCutBgEstZeeEnriched.systematics = cms.vstring()

cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZeeEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauPzetaDiffCut)
cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZeeEnriched.pluginName = cms.string('diTauCandidatePzetaDiffCutBgEstZeeEnriched')
cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZeeEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZeeEnrichedPzetaDiffCumulative')
cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZeeEnriched.systematics = cms.vstring()

cfgElecTauPairZeeHypothesisVetoBgEstZeeEnriched = copy.deepcopy(cfgElecTauPairZeeHypothesisVeto)
cfgElecTauPairZeeHypothesisVetoBgEstZeeEnriched.pluginName = cms.string('diTauCandidateZeeHypothesisVetoBgEstZeeEnriched')
cfgElecTauPairZeeHypothesisVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZeeEnrichedZeeHypothesisVetoCumulative')
cfgElecTauPairZeeHypothesisVetoBgEstZeeEnriched.systematics = cms.vstring()


evtSelConfiguratorBgEstZeeEnriched = eventSelFlagProdConfigurator(
    [ cfgElectronIdCutBgEstZeeEnriched,
      cfgElectronAntiCrackCutBgEstZeeEnriched,
      cfgElectronEtaCutBgEstZeeEnriched,
      cfgElectronPtCutBgEstZeeEnriched,
      cfgElectronTrkIsoCutBgEstZeeEnriched,
      cfgElectronEcalIsoCutBgEstZeeEnriched,
      cfgElectronConversionVetoBgEstZeeEnriched,
      #cfgElectronTrkIPcutBgEstZeeEnriched,
      cfgTauAntiOverlapWithElectronsVetoBgEstZeeEnriched,
      cfgTauEtaCutBgEstZeeEnriched,
      cfgTauPtCutBgEstZeeEnriched,
      cfgTauLeadTrkCutBgEstZeeEnriched,
      cfgTauLeadTrkPtCutBgEstZeeEnriched,
      cfgTauTaNCdiscrCutBgEstZeeEnriched,
      cfgTauTrkIsoCutBgEstZeeEnriched,
      cfgTauEcalIsoCutBgEstZeeEnriched,
      cfgTauProngCutBgEstZeeEnriched,
      cfgTauChargeCutBgEstZeeEnriched,
      cfgTauElectronVetoBgEstZeeEnriched,
      cfgTauEcalCrackVetoBgEstZeeEnriched,
      cfgTauMuonVetoBgEstZeeEnriched,
      cfgDiTauCandidateAntiOverlapVetoBgEstZeeEnriched,
      cfgDiTauCandidateForElecTauZeroChargeCutBgEstZeeEnriched,
      cfgDiTauCandidateForElecTauAcoplanarity12CutBgEstZeeEnriched,
      #cfgDiTauCandidateForElecTauMt1METCutBgEstZeeEnriched,
      cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZeeEnriched,
      cfgElecTauPairZeeHypothesisVetoBgEstZeeEnriched ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

selectEventsBgEstZeeEnriched = evtSelConfiguratorBgEstZeeEnriched.configure()


#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoElecTau_cff import *


diTauCandidateHistManagerBgEstZeeEnriched = copy.deepcopy(diTauCandidateHistManagerForElecTau)
diTauCandidateHistManagerBgEstZeeEnriched.pluginName = cms.string('diTauCandidateHistManagerBgEstZeeEnriched')
diTauCandidateHistManagerBgEstZeeEnriched.diTauCandidateSource = cms.InputTag('diTauCandidateZeeHypothesisVetoBgEstZeeEnriched')
diTauCandidateHistManagerBgEstZeeEnriched.visMassHypothesisSource = cms.InputTag('')


analyzeEventsBgEstZeeEnriched = cms.EDAnalyzer("GenericAnalyzer",

    name = cms.string('BgEstTemplateAnalyzer_ZeeEnriched'), 

    filters = cms.VPSet(
        genPhaseSpaceCut,
        evtSelTrigger,
        evtSelPrimaryEventVertex,
        evtSelPrimaryEventVertexQuality,
        evtSelPrimaryEventVertexPosition,

        #start electron cuts
        
        cms.PSet(
            pluginName = cms.string('electronIdCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronIdCutBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('electronAntiCrackCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronAntiCrackCutBgEstZeeEnriched','cumulative')
        ),        
        cms.PSet(
            pluginName = cms.string('electronEtaCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronEtaCutBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('electronPtCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronPtCutBgEstZeeEnriched','cumulative')
        ),        
        cms.PSet(
            pluginName = cms.string('electronTrkIsoCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronTrkIsoCutBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('electronEcalIsoCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronEcalIsoCutBgEstZeeEnriched','cumulative')
        ),
##         cms.PSet(
##             pluginName = cms.string('electronTrkIPcutBgEstZeeEnriched'),
##             pluginType = cms.string('BoolEventSelector'),
##             src = cms.InputTag('electronTrkIPcutBgEstZeeEnriched','cumulative')
##         ),
        cms.PSet(
            pluginName = cms.string('electronConversionVetoBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronConversionVetoBgEstZeeEnriched','cumulative')
        ),

        #start tau cuts
        
        cms.PSet(
            pluginName = cms.string('tauAntiOverlapWithElectronsVetoBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauAntiOverlapWithElectronsVetoBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauEtaCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauEtaCutBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauPtCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauPtCutBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauLeadTrkCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauLeadTrkCutBgEstZeeEnriched','cumulative')
        ),        
        cms.PSet(
            pluginName = cms.string('tauLeadTrkPtCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauLeadTrkPtCutBgEstZeeEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauTaNCdiscrCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauTaNCdiscrCutBgEstZeeEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauTrkIsoCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauTrkIsoCutBgEstZeeEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauEcalIsoCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauEcalIsoCutBgEstZeeEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauProngCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauProngCutBgEstZeeEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauChargeCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauChargeCutBgEstZeeEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauElectronVetoBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauElectronVetoBgEstZeeEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauEcalCrackVetoBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauEcalCrackVetoBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauMuonVetoBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauMuonVetoBgEstZeeEnriched','cumulative')
        ),

        #start ditau cuts
        
        cms.PSet(
            pluginName = cms.string('diTauCandidateAntiOverlapVetoBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateAntiOverlapVetoBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('diTauCandidateZeroChargeCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateZeroChargeCutBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('diTauCandidateAcoplanarity12CutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateAcoplanarity12CutBgEstZeeEnriched','cumulative')
        ),
##         cms.PSet(
##             pluginName = cms.string('diTauCandidateMt1METCutBgEstZeeEnriched'),
##             pluginType = cms.string('BoolEventSelector'),
##             src = cms.InputTag('diTauCandidateMt1METCutBgEstZeeEnriched','cumulative')
##         ),
        cms.PSet(
            pluginName = cms.string('diTauCandidatePzetaDiffCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidatePzetaDiffCutBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('diTauCandidateZeeHypothesisVetoBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateZeeHypothesisVetoBgEstZeeEnriched','cumulative')
        ),
        
        
    ),
  
    analyzers = cms.VPSet(
        diTauCandidateHistManagerBgEstZeeEnriched
    ),

    eventDumps = cms.VPSet(),
   
    analysisSequence = cms.VPSet(
    
        # generator level phase-space selection
        # (NOTE: (1) to be used in case of Monte Carlo samples
        #            overlapping in simulated phase-space only !!
        #        (2) genPhaseSpaceCut needs to be **always** the first entry in the list of cuts
        #           - otherwise the script submitToBatch.csh for submission of cmsRun jobs
        #            to the CERN batch system will not work !!)
        cms.PSet(
            filter = cms.string('genPhaseSpaceCut'),
            title = cms.string('gen. Phase-Space')
        ),
        cms.PSet(
            filter = cms.string('evtSelTrigger'),
            title = cms.string('Trigger')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertex'),
            title = cms.string('Vertex exists')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertexQuality'),
            title = cms.string('p(chi2Vertex) > 0.01')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertexPosition'),
            title = cms.string('-25 < zVertex < +25 cm')
        ),
        cms.PSet(
            filter = cms.string('electronIdCutBgEstZeeEnriched'),
            title = cms.string('Electron ID'),
        ),
        cms.PSet(
            filter = cms.string('electronAntiCrackCutBgEstZeeEnriched'),
            title = cms.string('Electron crack-Veto'),
        ),
        cms.PSet(
            filter = cms.string('electronEtaCutBgEstZeeEnriched'),
            title = cms.string('-2.1 < eta(Electron) < +2.1'),
        ),
        cms.PSet(
            filter = cms.string('electronPtCutBgEstZeeEnriched'),
            title = cms.string('Pt(Electron) > 30 GeV'),
        ),
        cms.PSet(
            filter = cms.string('tauAntiOverlapWithElectronsVetoBgEstZeeEnriched'),
            title = cms.string('Tau not overlapping with Elec.'),
        ),
        cms.PSet(
            filter = cms.string('tauEtaCutBgEstZeeEnriched'),
            title = cms.string('-2.1 < eta(Tau) < +2.1'),
        ),
        cms.PSet(
            filter = cms.string('tauPtCutBgEstZeeEnriched'),
            title = cms.string('Pt(Tau) > 20 GeV'),
        ),
        cms.PSet(
            filter = cms.string('electronTrkIsoCutBgEstZeeEnriched'),
            title = cms.string('Electron Track iso.'),
        ),
        cms.PSet(
            filter = cms.string('electronEcalIsoCutBgEstZeeEnriched'),
            title = cms.string('Electron ECAL iso.'),
        ),
        cms.PSet(
            filter = cms.string('electronConversionVetoBgEstZeeEnriched'),
            title = cms.string('Electron Track conv. veto'),
        ),
##         cms.PSet(
##             filter = cms.string('electronTrkIPcutBgEstZeeEnriched'),
##             title = cms.string('Electron Track IP'),
##         ),
        cms.PSet(
            filter = cms.string('tauLeadTrkCutBgEstZeeEnriched'),
            title = cms.string('Tau lead. Track find.'),
        ),
        cms.PSet(
            filter = cms.string('tauLeadTrkPtCutBgEstZeeEnriched'),
            title = cms.string('Tau lead. Track Pt'),
        ),
        cms.PSet(
            filter = cms.string('tauTaNCdiscrCutBgEstZeeEnriched'),
            title = cms.string('Tau TaNC at 0.5%'),
        ),
        cms.PSet(
            filter = cms.string('tauTrkIsoCutBgEstZeeEnriched'),
            title = cms.string('Tau Track iso.'),
        ),
        cms.PSet(
            filter = cms.string('tauEcalIsoCutBgEstZeeEnriched'),
            title = cms.string('Tau ECAL iso.'),
        ),
        cms.PSet(
            filter = cms.string('tauProngCutBgEstZeeEnriched'),
            title = cms.string('Tau 1||3-Prong'),
        ),
        cms.PSet(
            filter = cms.string('tauChargeCutBgEstZeeEnriched'),
            title = cms.string('Charge(Tau) = +/-1'),
        ),
        cms.PSet(
            filter = cms.string('tauElectronVetoBgEstZeeEnriched'),
            title = cms.string('Tau e-Veto'),
        ),
        cms.PSet(
            filter = cms.string('tauEcalCrackVetoBgEstZeeEnriched'),
            title = cms.string('Tau ECAL crack-Veto'),
        ),
        cms.PSet(
            filter = cms.string('tauMuonVetoBgEstZeeEnriched'),
            title = cms.string('Tau mu-Veto'),
        ),
        cms.PSet(
            filter = cms.string('diTauCandidateAntiOverlapVetoBgEstZeeEnriched'),
            title = cms.string('dR(Electron-Tau) > 0.7'),
        ),
        cms.PSet(
            filter = cms.string('diTauCandidateZeroChargeCutBgEstZeeEnriched'),
            title = cms.string('Charge(Electron+Tau) = 0'),
        ),
        cms.PSet(
            filter = cms.string('diTauCandidateAcoplanarity12CutBgEstZeeEnriched'),
            title = cms.string('Acoplanarity(Electron+Tau)'),
        ),
##         cms.PSet(
##             filter = cms.string('diTauCandidateMt1METCutBgEstZeeEnriched'),
##             title = cms.string('M_{T}(Electron-MET) < 50 GeV'),
##         ),
        cms.PSet(
            filter = cms.string('diTauCandidatePzetaDiffCutBgEstZeeEnriched'),
            title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV'),
        ),
        cms.PSet(
            filter = cms.string('diTauCandidateZeeHypothesisVetoBgEstZeeEnriched'),
            title = cms.string('not 85 < M_{vis} (Electron-Tau) < 100 GeV'),
        ),




        
        cms.PSet(
            analyzers = cms.vstring(
                'diTauCandidateHistManagerBgEstZeeEnriched',
            )
        )

    )

)




#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstZeeEnrichedAnalysisSequence = cms.Sequence(

    selectElectronsBgEstZeeEnriched
    + selectTausBgEstZeeEnriched
    + produceElecTauPairsBgEstZeeEnriched
    + selectElecTauPairsBgEstZeeEnriched
    + selectEventsBgEstZeeEnriched 
    + analyzeEventsBgEstZeeEnriched

)
