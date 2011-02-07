
import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

#--------------------------------------------------------------------------------
# select Z --> e+ e- background enriched event sample
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------  
# produce collection of good vertices
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.recoVertexSelection_cff import *
selectPrimaryVertexBgEst = cms.Sequence(selectPrimaryVertex)


#--------------------------------------------------------------------------------  
# produce collection of pat::Electrons
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patLeptonSelection_cff import *

# require electron candidate to pass the eidRobustTight electron id. criteria
electronsBgEstZtautauEnrichedId = copy.deepcopy(selectedPatElectronsForElecTauId)

# require electron candidate to not be within eta-crack
# between Barrel and Encap ECAL calorimeter
electronsBgEstZtautauEnrichedAntiCrackCut = copy.deepcopy(selectedPatElectronsForElecTauAntiCrackCut)

# require electron candidate to be within geometric acceptance of electron trigger
electronsBgEstZtautauEnrichedEta = copy.deepcopy(selectedPatElectronsForElecTauEta)

# require electron candidate to have transverse momentum above threshold
electronsBgEstZtautauEnrichedPt = copy.deepcopy(selectedPatElectronsForElecTauPt)

# require electron candidate to be isolated w.r.t combined relative PF candidates
electronsBgEstZtautauEnrichedIso = copy.deepcopy(selectedPatElectronsForElecTauIso)

# require electron to not be from a photon conversion
electronsBgEstZtautauEnrichedConversionVeto = copy.deepcopy(selectedPatElectronsForElecTauConversionVeto)

# require track of electron candidate to have small transverse impact parameter
# (in order to veto electrons resulting from b-quark decays)
electronsBgEstZtautauEnrichedTrkIP = copy.deepcopy(selectedPatElectronsForElecTauTrkIP)


electronSelConfiguratorBgEstZtautauEnriched = objSelConfigurator(
    [ electronsBgEstZtautauEnrichedId,
      electronsBgEstZtautauEnrichedAntiCrackCut,
      electronsBgEstZtautauEnrichedEta,
      electronsBgEstZtautauEnrichedPt,
      electronsBgEstZtautauEnrichedIso,
      electronsBgEstZtautauEnrichedConversionVeto,
      electronsBgEstZtautauEnrichedTrkIP ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = False
)


selectElectronsBgEstZtautauEnriched = electronSelConfiguratorBgEstZtautauEnriched.configure(pyNameSpace = locals())

   
#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

#selectTausBgEstZtautauEnriched = copy.deepcopy(selectPatTausForElecTau)

# require tau candidate not to overlap with selected electrons
# (in order to avoid double-counting one and the same physical particle
#  as electron and as tau candidate)
tausBgEstZtautauEnrichedAntiOverlapWithElectronsVeto = cms.EDFilter("PATTauAntiOverlapSelector",
    srcNotToBeFiltered = cms.VInputTag("electronsBgEstZtautauEnrichedPtCumulative"),
    dRmin = cms.double(0.3),
    filter = cms.bool(True)                     
)

# require tau candidate to be within geometric acceptance of Pixel + SiTracker detectors
tausBgEstZtautauEnrichedEta = copy.deepcopy(selectedPatTausForElecTauEta)

# require tau candidate to have transverse energy above threshold
tausBgEstZtautauEnrichedPt = copy.deepcopy(selectedPatTausForElecTauPt)

# require tau candidate to have a leading track
# (track of Pt > 1. GeV within matching cone of size dR = 0.2 around jet-axis)
tausBgEstZtautauEnrichedLeadTrk = copy.deepcopy(selectedPatTausForElecTauLeadTrk)

# require leading track of tau candidate to have Pt > 5. GeV
tausBgEstZtautauEnrichedLeadTrkPt = copy.deepcopy(selectedPatTausForElecTauLeadTrkPt)

# require tau candidate to pass TaNC discriminator
tausBgEstZtautauEnrichedTaNCdiscr = copy.deepcopy(selectedPatTausForElecTauTaNCdiscr)

# require tau candidate to have either one or three tracks within signal cone
tausBgEstZtautauEnrichedProng = copy.deepcopy(selectedPatTausForElecTauProng)

# require tau candidate to have charge either +1 or -1
# (computed as sum of charges of tracks within signal cone)
tausBgEstZtautauEnrichedCharge = copy.deepcopy(selectedPatTausForElecTauCharge)

# require tau candidate to pass electron veto
tausBgEstZtautauEnrichedElectronVeto = copy.deepcopy(selectedPatTausForElecTauElectronVeto)

# require tau candidate not to be in ECAL barrel/endcap crack
tausBgEstZtautauEnrichedEcalCrackVeto = copy.deepcopy(selectedPatTausForElecTauEcalCrackVeto)

# require tau candidate to pass muon veto
tausBgEstZtautauEnrichedMuonVeto = copy.deepcopy(selectedPatTausForElecTauMuonVeto)



tauSelConfiguratorBgEstZtautauEnriched = objSelConfigurator(
    [ tausBgEstZtautauEnrichedAntiOverlapWithElectronsVeto,
      tausBgEstZtautauEnrichedEta,
      tausBgEstZtautauEnrichedPt,
      tausBgEstZtautauEnrichedLeadTrk,
      tausBgEstZtautauEnrichedLeadTrkPt,
      tausBgEstZtautauEnrichedTaNCdiscr,
      tausBgEstZtautauEnrichedProng,
      tausBgEstZtautauEnrichedCharge,
      tausBgEstZtautauEnrichedElectronVeto,
      tausBgEstZtautauEnrichedEcalCrackVeto,
      tausBgEstZtautauEnrichedMuonVeto ],
    src = "cleanPatTaus",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectTausBgEstZtautauEnriched = tauSelConfiguratorBgEstZtautauEnriched.configure(pyNameSpace = locals())


#--------------------------------------------------------------------------------  
# produce collection of electron + tau-jet combinations
#--------------------------------------------------------------------------------

### production

from TauAnalysis.CandidateTools.elecTauPairProduction_cff import *

elecTauPairsBgEstZtautauEnriched = copy.deepcopy(allElecTauPairs)
elecTauPairsBgEstZtautauEnriched.srcLeg1 = cms.InputTag('electronsBgEstZtautauEnrichedTrkIPcumulative')
elecTauPairsBgEstZtautauEnriched.srcLeg2 = cms.InputTag('tausBgEstZtautauEnrichedMuonVetoCumulative')

produceElecTauPairsBgEstZtautauEnriched = cms.Sequence(elecTauPairsBgEstZtautauEnriched)

### selection

from TauAnalysis.CandidateTools.elecTauPairSelection_cfi import *

elecTauPairsBgEstZtautauEnrichedAntiOverlapVeto = copy.deepcopy(selectedElecTauPairsAntiOverlapVeto)
elecTauPairsBgEstZtautauEnrichedMt1MET = copy.deepcopy(selectedElecTauPairsMt1MET)
elecTauPairsBgEstZtautauEnrichedPzetaDiff = copy.deepcopy(selectedElecTauPairsPzetaDiff)
elecTauPairsBgEstZtautauEnrichedZeroCharge = copy.deepcopy(selectedElecTauPairsZeroCharge)

elecTauPairSelConfiguratorBgEstZtautauEnriched = objSelConfigurator(
    [ elecTauPairsBgEstZtautauEnrichedAntiOverlapVeto,
      elecTauPairsBgEstZtautauEnrichedMt1MET,
      elecTauPairsBgEstZtautauEnrichedPzetaDiff,
      elecTauPairsBgEstZtautauEnrichedZeroCharge,
      ],
    src = "elecTauPairsBgEstZtautauEnriched",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectElecTauPairsBgEstZtautauEnriched = elecTauPairSelConfiguratorBgEstZtautauEnriched.configure(pyNameSpace = locals())

####### anti Zee Cut #######

# Two methods of Z->ee rejection:

# 1)  Z->ee reconstruction veto
# 2)  Opposite-sign, isolated second lepton veto

#from TauAnalysis.RecoTools.elecTauPairZeeHypothesis_cff import *
from TauAnalysis.RecoTools.diElecPairZeeHypothesis_cff import *

### production

#elecTauPairBgEstZtautauEnrichedZeeHypotheses = copy.deepcopy(elecTauPairZeeHypotheses)
#elecTauPairBgEstZtautauEnrichedZeeHypotheses.diCandidatePairSource = cms.InputTag('elecTauPairsBgEstZtautauEnrichedPzetaDiffCumulative')
allDiElecPairBgEstZtautauEnrichedZeeHypothesesByLooseIsolation = allDiElecPairZeeHypothesesByLooseIsolation.clone(
	srcLeg1 = cms.InputTag("selectedPatElectronsBgEstZtautauEnrichedTrkIPcumulative")
)

### selection

#selectedElecTauPairBgEstZtautauEnrichedZeeHypotheses = copy.deepcopy(selectedElecTauPairZeeHypotheses)
#selectedElecTauPairBgEstZtautauEnrichedZeeHypotheses.src = cms.InputTag('elecTauPairBgEstZtautauEnrichedZeeHypotheses')
selectedDiElecPairBgEstZtautauEnrichedZeeHypothesesByLooseIsolation = selectedDiElecPairZeeHypothesesByLooseIsolation.clone(
	src = cms.InputTag("allDiElecPairBgEstZtautauEnrichedZeeHypothesesByLooseIsolation")
)

#produceElecTauPairZeeHypothesesBgEstZtautauEnriched = cms.Sequence(elecTauPairBgEstZtautauEnrichedZeeHypotheses*selectedElecTauPairBgEstZtautauEnrichedZeeHypotheses)
produceDiElecPairsBgEstZtautauEnriched = cms.Sequence(
	selectedPatElectronsForZeeHypotheses * 
	allDiElecPairBgEstZtautauEnrichedZeeHypothesesByLooseIsolation *
	selectedDiElecPairBgEstZtautauEnrichedZeeHypothesesByLooseIsolation
)

#--------------------------------------------------------------------------------  
# produce boolean event selection flags
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.selectZtoElecTau_cff import *

# trigger selection
cfgTriggerBgEst = copy.deepcopy(cfgTrigger)

# primary event vertex selection
cfgPrimaryEventVertexBgEst = copy.deepcopy(cfgPrimaryEventVertex)
cfgPrimaryEventVertexQualityBgEst = copy.deepcopy(cfgPrimaryEventVertexQuality)
cfgPrimaryEventVertexPositionBgEst = copy.deepcopy(cfgPrimaryEventVertexPosition)



# electron cuts

cfgElectronIdCutBgEstZtautauEnriched = copy.deepcopy(cfgElectronIdCut)
cfgElectronIdCutBgEstZtautauEnriched.pluginName = cms.string('electronIdCutBgEstZtautauEnriched')
cfgElectronIdCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('electronsBgEstZtautauEnrichedIdCumulative')

cfgElectronAntiCrackCutBgEstZtautauEnriched = copy.deepcopy(cfgElectronAntiCrackCut)
cfgElectronAntiCrackCutBgEstZtautauEnriched.pluginName = cms.string('electronAntiCrackCutBgEstZtautauEnriched')
cfgElectronAntiCrackCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('electronsBgEstZtautauEnrichedAntiCrackCutCumulative')

cfgElectronEtaCutBgEstZtautauEnriched = copy.deepcopy(cfgElectronEtaCut)
cfgElectronEtaCutBgEstZtautauEnriched.pluginName = cms.string('electronEtaCutBgEstZtautauEnriched')
cfgElectronEtaCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('electronsBgEstZtautauEnrichedEtaCumulative')

cfgElectronPtCutBgEstZtautauEnriched = copy.deepcopy(cfgElectronPtCut)
cfgElectronPtCutBgEstZtautauEnriched.pluginName = cms.string('electronPtCutBgEstZtautauEnriched')
cfgElectronPtCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('electronsBgEstZtautauEnrichedPtCumulative')

cfgElectronIsoCutBgEstZtautauEnriched = copy.deepcopy(cfgElectronIsoCut)
cfgElectronIsoCutBgEstZtautauEnriched.pluginName = cms.string('electronIsoCutBgEstZtautauEnriched')
cfgElectronIsoCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('electronsBgEstZtautauEnrichedIsoCumulative')

cfgElectronConversionVetoBgEstZtautauEnriched = copy.deepcopy(cfgElectronConversionVeto)
cfgElectronConversionVetoBgEstZtautauEnriched.pluginName = cms.string('electronConversionVetoBgEstZtautauEnriched')
cfgElectronConversionVetoBgEstZtautauEnriched.src_cumulative = cms.InputTag('electronsBgEstZtautauEnrichedConversionVetoCumulative')

cfgElectronTrkIPcutBgEstZtautauEnriched = copy.deepcopy(cfgElectronTrkIPcut)
cfgElectronTrkIPcutBgEstZtautauEnriched.pluginName = cms.string('electronTrkIPcutBgEstZtautauEnriched')
cfgElectronTrkIPcutBgEstZtautauEnriched.src_cumulative = cms.InputTag('electronsBgEstZtautauEnrichedTrkIPcumulative')

# tau cuts

cfgTauAntiOverlapWithElectronsVetoBgEstZtautauEnriched = copy.deepcopy(cfgTauAntiOverlapWithElectronsVeto)
cfgTauAntiOverlapWithElectronsVetoBgEstZtautauEnriched.pluginName = cms.string('tauAntiOverlapWithElectronsVetoBgEstZtautauEnriched')
cfgTauAntiOverlapWithElectronsVetoBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedAntiOverlapWithElectronsVetoCumulative')

cfgTauEtaCutBgEstZtautauEnriched = copy.deepcopy(cfgTauEtaCut)
cfgTauEtaCutBgEstZtautauEnriched.pluginName = cms.string('tauEtaCutBgEstZtautauEnriched')
cfgTauEtaCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedEtaCumulative')

cfgTauPtCutBgEstZtautauEnriched = copy.deepcopy(cfgTauPtCut)
cfgTauPtCutBgEstZtautauEnriched.pluginName = cms.string('tauPtCutBgEstZtautauEnriched')
cfgTauPtCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedPtCumulative')

cfgTauLeadTrkCutBgEstZtautauEnriched = copy.deepcopy(cfgTauLeadTrkCut)
cfgTauLeadTrkCutBgEstZtautauEnriched.pluginName = cms.string('tauLeadTrkCutBgEstZtautauEnriched')
cfgTauLeadTrkCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedLeadTrkCumulative')

cfgTauLeadTrkPtCutBgEstZtautauEnriched = copy.deepcopy(cfgTauLeadTrkPtCut)
cfgTauLeadTrkPtCutBgEstZtautauEnriched.pluginName = cms.string('tauLeadTrkPtCutBgEstZtautauEnriched')
cfgTauLeadTrkPtCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedLeadTrkPtCumulative')

cfgTauTaNCdiscrCutBgEstZtautauEnriched = copy.deepcopy(cfgTauTaNCdiscrCut)
cfgTauTaNCdiscrCutBgEstZtautauEnriched.pluginName = cms.string('tauTaNCdiscrCutBgEstZtautauEnriched')
cfgTauTaNCdiscrCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedTaNCdiscrCumulative')

cfgTauProngCutBgEstZtautauEnriched = copy.deepcopy(cfgTauProngCut)
cfgTauProngCutBgEstZtautauEnriched.pluginName = cms.string('tauProngCutBgEstZtautauEnriched')
cfgTauProngCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedProngCumulative')

cfgTauChargeCutBgEstZtautauEnriched = copy.deepcopy(cfgTauChargeCut)
cfgTauChargeCutBgEstZtautauEnriched.pluginName = cms.string('tauChargeCutBgEstZtautauEnriched')
cfgTauChargeCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedChargeCumulative')

cfgTauElectronVetoBgEstZtautauEnriched = copy.deepcopy(cfgTauElectronVeto)
cfgTauElectronVetoBgEstZtautauEnriched.pluginName = cms.string('tauElectronVetoBgEstZtautauEnriched')
cfgTauElectronVetoBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedElectronVetoCumulative')

cfgTauEcalCrackVetoBgEstZtautauEnriched = copy.deepcopy(cfgTauEcalCrackVeto)
cfgTauEcalCrackVetoBgEstZtautauEnriched.pluginName = cms.string('tauEcalCrackVetoBgEstZtautauEnriched')
cfgTauEcalCrackVetoBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedEcalCrackVetoCumulative')

cfgTauMuonVetoBgEstZtautauEnriched = copy.deepcopy(cfgTauMuonVeto)
cfgTauMuonVetoBgEstZtautauEnriched.pluginName = cms.string('tauMuonVetoBgEstZtautauEnriched')
cfgTauMuonVetoBgEstZtautauEnriched.src_cumulative = cms.InputTag('tausBgEstZtautauEnrichedMuonVetoCumulative')

# di-tau candidate selection

cfgDiTauCandidateAntiOverlapVetoBgEstZtautauEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauAntiOverlapVeto)
cfgDiTauCandidateAntiOverlapVetoBgEstZtautauEnriched.pluginName = cms.string('diTauCandidateAntiOverlapVetoBgEstZtautauEnriched')
cfgDiTauCandidateAntiOverlapVetoBgEstZtautauEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZtautauEnrichedAntiOverlapVetoCumulative')

cfgDiTauCandidateForElecTauMt1METCutBgEstZtautauEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauMt1METCut)
cfgDiTauCandidateForElecTauMt1METCutBgEstZtautauEnriched.pluginName = cms.string('diTauCandidateMt1METCutBgEstZtautauEnriched')
cfgDiTauCandidateForElecTauMt1METCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZtautauEnrichedMt1METcumulative')

cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZtautauEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauPzetaDiffCut)
cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZtautauEnriched.pluginName = cms.string('diTauCandidatePzetaDiffCutBgEstZtautauEnriched')
cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZtautauEnrichedPzetaDiffCumulative')

cfgDiTauCandidateForElecTauZeroChargeCutBgEstZtautauEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauZeroChargeCut)
cfgDiTauCandidateForElecTauZeroChargeCutBgEstZtautauEnriched.pluginName = cms.string('diTauCandidateZeroChargeCutBgEstZtautauEnriched')
cfgDiTauCandidateForElecTauZeroChargeCutBgEstZtautauEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZtautauEnrichedZeroChargeCumulative')

#cfgElecTauPairZeeHypothesisVetoBgEstZtautauEnriched = copy.deepcopy(cfgElecTauPairZeeHypothesisVeto)
#cfgElecTauPairZeeHypothesisVetoBgEstZtautauEnriched.pluginName = cms.string('diTauCandidateZeeHypothesisVetoBgEstZtautauEnriched')
#cfgElecTauPairZeeHypothesisVetoBgEstZtautauEnriched.src_cumulative = cms.InputTag('selectedElecTauPairBgEstZtautauEnrichedZeeHypotheses')

cfgDiElecPairZeeHypothesisVetoByLooseIsolationBgEstZtautauEnriched = cfgDiElecPairZeeHypothesisVetoByLooseIsolation.clone(
	pluginName = cms.string('diElecPairZeeHypothesisVetoByLooseIsolationBgEstZtautauEnriched'),
	src = cms.InputTag('selectedDiElecPairBgEstZtautauEnrichedZeeHypothesesByLooseIsolation'),
	maxNumber = cms.uint32(0)
)

evtSelConfiguratorBgEstZtautauEnriched = eventSelFlagProdConfigurator(
    [ cfgTriggerBgEst,
      cfgPrimaryEventVertexBgEst,
      cfgPrimaryEventVertexQualityBgEst,
      cfgPrimaryEventVertexPositionBgEst,
      cfgElectronIdCutBgEstZtautauEnriched,
      cfgElectronAntiCrackCutBgEstZtautauEnriched,
      cfgElectronEtaCutBgEstZtautauEnriched,
      cfgElectronPtCutBgEstZtautauEnriched,
      cfgElectronIsoCutBgEstZtautauEnriched,
      cfgElectronConversionVetoBgEstZtautauEnriched,
      cfgElectronTrkIPcutBgEstZtautauEnriched,
      cfgTauAntiOverlapWithElectronsVetoBgEstZtautauEnriched,
      cfgTauEtaCutBgEstZtautauEnriched,
      cfgTauPtCutBgEstZtautauEnriched,
      cfgTauLeadTrkCutBgEstZtautauEnriched,
      cfgTauLeadTrkPtCutBgEstZtautauEnriched,
      cfgTauTaNCdiscrCutBgEstZtautauEnriched,
      cfgTauProngCutBgEstZtautauEnriched,
      cfgTauChargeCutBgEstZtautauEnriched,
      cfgTauElectronVetoBgEstZtautauEnriched,
      cfgTauEcalCrackVetoBgEstZtautauEnriched,
      cfgTauMuonVetoBgEstZtautauEnriched,
      cfgDiTauCandidateAntiOverlapVetoBgEstZtautauEnriched,
      cfgDiTauCandidateForElecTauMt1METCutBgEstZtautauEnriched,
      cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZtautauEnriched,
      cfgDiTauCandidateForElecTauZeroChargeCutBgEstZtautauEnriched,
	  #cfgElecTauPairZeeHypothesisVetoBgEstZtautauEnriched
	  cfgDiElecPairZeeHypothesisVetoByLooseIsolationBgEstZtautauEnriched
      ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

selectEventsBgEstZtautauEnriched = evtSelConfiguratorBgEstZtautauEnriched.configure()


#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoElecTau_cff import *


diTauCandidateHistManagerForElecTauBgEstZtautauEnriched = copy.deepcopy(diTauCandidateHistManagerForElecTau)
diTauCandidateHistManagerForElecTauBgEstZtautauEnriched.pluginName = cms.string('diTauCandidateHistManagerForElecTauBgEstZtautauEnriched')
diTauCandidateHistManagerForElecTauBgEstZtautauEnriched.diTauCandidateSource = cms.InputTag('elecTauPairsBgEstZtautauEnrichedZeroChargeCumulative')
#diTauCandidateHistManagerForElecTauBgEstZtautauEnriched.ZllHypothesisSource = cms.InputTag('selectedElecTauPairBgEstZtautauEnrichedZeeHypotheses')



analyzeEventsBgEstZtautauEnriched = cms.EDAnalyzer("GenericAnalyzer",

    name = cms.string('BgEstTemplateAnalyzer_ZtautauEnriched'), 

    filters = cms.VPSet(
        evtSelGenPhaseSpace,
        evtSelTrigger,
        evtSelPrimaryEventVertex,
        evtSelPrimaryEventVertexQuality,
        evtSelPrimaryEventVertexPosition,

        #start electron cuts
        
        cms.PSet(
            pluginName = cms.string('electronIdCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronIdCutBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('electronAntiCrackCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronAntiCrackCutBgEstZtautauEnriched','cumulative')
        ),        
        cms.PSet(
            pluginName = cms.string('electronEtaCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronEtaCutBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('electronPtCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronPtCutBgEstZtautauEnriched','cumulative')
        ),        
        cms.PSet(
            pluginName = cms.string('electronIsoCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronIsoCutBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('electronConversionVetoBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronConversionVetoBgEstZtautauEnriched','cumulative')
        ),        
        cms.PSet(
            pluginName = cms.string('electronTrkIPcutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronTrkIPcutBgEstZtautauEnriched','cumulative')
        ),


        #start tau cuts
        
        cms.PSet(
            pluginName = cms.string('tauAntiOverlapWithElectronsVetoBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauAntiOverlapWithElectronsVetoBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauEtaCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauEtaCutBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauPtCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauPtCutBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauLeadTrkCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauLeadTrkCutBgEstZtautauEnriched','cumulative')
        ),        
        cms.PSet(
            pluginName = cms.string('tauLeadTrkPtCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauLeadTrkPtCutBgEstZtautauEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauTaNCdiscrCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauTaNCdiscrCutBgEstZtautauEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauProngCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauProngCutBgEstZtautauEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauChargeCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauChargeCutBgEstZtautauEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauElectronVetoBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauElectronVetoBgEstZtautauEnriched','cumulative')
        ),  
        cms.PSet(
            pluginName = cms.string('tauEcalCrackVetoBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauEcalCrackVetoBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('tauMuonVetoBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('tauMuonVetoBgEstZtautauEnriched','cumulative')
        ),

        #start ditau cuts
        
        cms.PSet(
            pluginName = cms.string('diTauCandidateAntiOverlapVetoBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateAntiOverlapVetoBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('diTauCandidateMt1METCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateMt1METCutBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('diTauCandidatePzetaDiffCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidatePzetaDiffCutBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('diTauCandidateZeroChargeCutBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateZeroChargeCutBgEstZtautauEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('diElecPairZeeHypothesisVetoByLooseIsolationBgEstZtautauEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diElecPairZeeHypothesisVetoByLooseIsolationBgEstZtautauEnriched','cumulative')
        )
		#cms.PSet(
		#    pluginName = cms.string('diTauCandidateZeeHypothesisVetoBgEstZtautauEnriched'),
		#    pluginType = cms.string('BoolEventSelector'),
		#    src = cms.InputTag('diTauCandidateZeeHypothesisVetoBgEstZtautauEnriched','cumulative')
		#)
        
        
    ),
  
    analyzers = cms.VPSet(
         diTauCandidateHistManagerForElecTauBgEstZtautauEnriched
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
            filter = cms.string('electronIdCutBgEstZtautauEnriched'),
            title = cms.string('Electron ID'),
        ),
        cms.PSet(
            filter = cms.string('electronAntiCrackCutBgEstZtautauEnriched'),
            title = cms.string('Electron crack-Veto'),
        ),
        cms.PSet(
            filter = cms.string('electronEtaCutBgEstZtautauEnriched'),
            title = cms.string('-2.1 < eta(Electron) < +2.1'),
        ),
        cms.PSet(
            filter = cms.string('electronPtCutBgEstZtautauEnriched'),
            title = cms.string('Pt(Electron) > 15 GeV'),
        ),
        cms.PSet(
            filter = cms.string('tauAntiOverlapWithElectronsVetoBgEstZtautauEnriched'),
            title = cms.string('Tau not overlapping with Elec.'),
        ),
        cms.PSet(
            filter = cms.string('tauEtaCutBgEstZtautauEnriched'),
            title = cms.string('-2.1 < eta(Tau) < +2.1'),
        ),
        cms.PSet(
            filter = cms.string('tauPtCutBgEstZtautauEnriched'),
            title = cms.string('Pt(Tau) > 20 GeV'),
        ),
        cms.PSet(
            filter = cms.string('electronIsoCutBgEstZtautauEnriched'),
            title = cms.string('Electron ECAL iso.'),
        ),
        cms.PSet(
            filter = cms.string('electronConversionVetoBgEstZtautauEnriched'),
            title = cms.string('Electron Track conv. veto'),
        ),
        cms.PSet(
            filter = cms.string('electronTrkIPcutBgEstZtautauEnriched'),
            title = cms.string('Electron Track IP'),
        ),
        cms.PSet(
            filter = cms.string('tauLeadTrkPtCutBgEstZtautauEnriched'),
            title = cms.string('Tau lead. Track Pt'),
        ),
        cms.PSet(
            filter = cms.string('tauTaNCdiscrCutBgEstZtautauEnriched'),
            title = cms.string('Tau TaNC at 0.5%'),
        ),
        cms.PSet(
            filter = cms.string('tauProngCutBgEstZtautauEnriched'),
            title = cms.string('Tau 1||3-Prong'),
        ),
        cms.PSet(
            filter = cms.string('tauChargeCutBgEstZtautauEnriched'),
            title = cms.string('Charge(Tau) = +/-1'),
        ),
        cms.PSet(
            filter = cms.string('tauElectronVetoBgEstZtautauEnriched'),
            title = cms.string('Tau e-Veto'),
        ),
        cms.PSet(
            filter = cms.string('tauEcalCrackVetoBgEstZtautauEnriched'),
            title = cms.string('Tau ECAL crack-Veto'),
        ),
        cms.PSet(
            filter = cms.string('tauMuonVetoBgEstZtautauEnriched'),
            title = cms.string('Tau mu-Veto'),
        ),
        cms.PSet(
            filter = cms.string('diTauCandidateAntiOverlapVetoBgEstZtautauEnriched'),
            title = cms.string('dR(Electron-Tau) > 0.7'),
        ),
        cms.PSet(
            filter = cms.string('diTauCandidateMt1METCutBgEstZtautauEnriched'),
            title = cms.string('M_{T}(Electron-MET) < 50 GeV'),
        ),
        cms.PSet(
            filter = cms.string('diTauCandidatePzetaDiffCutBgEstZtautauEnriched'),
            title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV'),
        ),
        cms.PSet(
            filter = cms.string('diTauCandidateZeroChargeCutBgEstZtautauEnriched'),
            title = cms.string('Charge(Electron+Tau) = 0'),
        ),
		cms.PSet(
		    filter = cms.string('diElecPairZeeHypothesisVetoByLooseIsolationBgEstZtautauEnriched'),
		    title = cms.string('charge(elec + looseIsoElec) != 0'),
		),
		#cms.PSet(
		#    filter = cms.string('diTauCandidateZeeHypothesisVetoBgEstZtautauEnriched'),
		#    title = cms.string('not 85 < M_{vis} (Electron-Tau) < 100 GeV'),
		#),




        
        cms.PSet(
            analyzers = cms.vstring(
                'diTauCandidateHistManagerForElecTauBgEstZtautauEnriched'
            )
        )

    )

)

saveEvents = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_output.root')
)




#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstZtautauEnrichedAnalysisSequence = cms.Sequence(

    selectPrimaryVertexBgEst
    + selectElectronsBgEstZtautauEnriched
    + selectTausBgEstZtautauEnriched
    + produceElecTauPairsBgEstZtautauEnriched
    + selectElecTauPairsBgEstZtautauEnriched
	#+ produceElecTauPairZeeHypothesesBgEstZtautauEnriched
	+ produceDiElecPairsBgEstZtautauEnriched
    + selectEventsBgEstZtautauEnriched 
    + analyzeEventsBgEstZtautauEnriched

)
