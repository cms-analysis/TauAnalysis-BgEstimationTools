import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

#--------------------------------------------------------------------------------
# select Z --> e+ e- background enriched event sample
#--------------------------------------------------------------------------------

from TauAnalysis.BgEstimationTools.bgEstZtoElecTauZtautauEnrichedSelection_cff import *

#--------------------------------------------------------------------------------  
# produce collection of pat::Electrons
#--------------------------------------------------------------------------------

# require electron candidate to pass the eidRobustTight electron id. criteria
electronsBgEstZeeEnrichedId = copy.deepcopy(electronsBgEstZtautauEnrichedId)

# require electron candidate to not be within eta-crack
# between Barrel and Encap ECAL calorimeter
electronsBgEstZeeEnrichedAntiCrackCut = copy.deepcopy(electronsBgEstZtautauEnrichedAntiCrackCut)

# require electron candidate to be within geometric acceptance of electron trigger
electronsBgEstZeeEnrichedEta = copy.deepcopy(electronsBgEstZtautauEnrichedEta)

# require electron candidate to have transverse momentum above threshold
electronsBgEstZeeEnrichedPt = copy.deepcopy(electronsBgEstZtautauEnrichedPt)

electronsBgEstZeeEnrichedIso = copy.deepcopy(electronsBgEstZtautauEnrichedIso)

# require electron to not be from a photon conversion
electronsBgEstZeeEnrichedConversionVeto = copy.deepcopy(electronsBgEstZtautauEnrichedConversionVeto)

electronsBgEstZeeEnrichedTrkIP = copy.deepcopy(electronsBgEstZtautauEnrichedTrkIP)


electronSelConfiguratorBgEstZeeEnriched = objSelConfigurator(
    [ electronsBgEstZeeEnrichedId,
      electronsBgEstZeeEnrichedAntiCrackCut,
      electronsBgEstZeeEnrichedEta,
      electronsBgEstZeeEnrichedPt,
      electronsBgEstZeeEnrichedIso,
      electronsBgEstZeeEnrichedConversionVeto,
      electronsBgEstZeeEnrichedTrkIP
      ],
    src = "cleanPatElectrons",
    pyModuleName = __name__,
    doSelIndividual = False
)


selectElectronsBgEstZeeEnriched = electronSelConfiguratorBgEstZeeEnriched.configure(pyNameSpace = locals())

   
#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

# require tau candidate not to overlap with selected electrons
# (in order to avoid double-counting one and the same physical particle
#  as electron and as tau candidate)
tausBgEstZeeEnrichedAntiOverlapWithElectronsVeto = copy.deepcopy(tausBgEstZtautauEnrichedAntiOverlapWithElectronsVeto)
tausBgEstZeeEnrichedAntiOverlapWithElectronsVeto.srcNotToBeFiltered = cms.VInputTag("electronsBgEstZeeEnrichedPtCumulative")

# require tau candidate to be within geometric acceptance of Pixel + SiTracker detectors
tausBgEstZeeEnrichedEta = copy.deepcopy(tausBgEstZtautauEnrichedEta)

# require tau candidate to have transverse energy above threshold
tausBgEstZeeEnrichedPt = copy.deepcopy(tausBgEstZtautauEnrichedPt)

# require leading track of tau candidate to have Pt > 5. GeV
tausBgEstZeeEnrichedLeadTrkPt = copy.deepcopy(tausBgEstZtautauEnrichedLeadTrkPt)

# require tau candidate to pass TaNC discriminator
tausBgEstZeeEnrichedTaNCdiscr = copy.deepcopy(tausBgEstZtautauEnrichedTaNCdiscr)
#tausBgEstZeeEnrichedTaNCdiscr.cut = cms.string('tauID("byHPSloose") > -1')
tausBgEstZeeEnrichedTaNCdiscr.cut = cms.string('tauID("byTaNCfrHalfPercent") > -1')


# require tau candidate to have either one or three tracks within signal cone
tausBgEstZeeEnrichedProng = copy.deepcopy(tausBgEstZtautauEnrichedProng)
tausBgEstZeeEnrichedProng.cut = cms.string("signalPFChargedHadrCands.size() != -1")

# require tau candidate to have charge either +1 or -1
# (computed as sum of charges of tracks within signal cone)
tausBgEstZeeEnrichedCharge = copy.deepcopy(tausBgEstZtautauEnrichedCharge)
tausBgEstZeeEnrichedCharge.cut = cms.string('abs(charge) > -1')

# require tau candidate to fail electron veto  !!!!! inverted cut
tausBgEstZeeEnrichedElectronVeto = copy.deepcopy(tausBgEstZtautauEnrichedElectronVeto)
#tausBgEstZeeEnrichedElectronVeto.cut = cms.string('!leadPFCand().isNonnull() | leadPFCand().mva_e_pi() > -0.1 | hcalTotOverPLead() < 0.1')
#tausBgEstZeeEnrichedElectronVeto.cut = cms.string('leadPFCand().mva_e_pi() > -0.1 | hcalTotOverPLead() < 0.1')
tausBgEstZeeEnrichedElectronVeto.cut = cms.string('leadPFCand().mva_e_pi() > -0.1')


# require tau candidate not to be in ECAL barrel/endcap crack
tausBgEstZeeEnrichedEcalCrackVeto = copy.deepcopy(tausBgEstZtautauEnrichedEcalCrackVeto)


# require tau candidate to pass muon veto
tausBgEstZeeEnrichedMuonVeto = copy.deepcopy(tausBgEstZtautauEnrichedMuonVeto)



tauSelConfiguratorBgEstZeeEnriched = objSelConfigurator(
    [ tausBgEstZeeEnrichedAntiOverlapWithElectronsVeto,
      tausBgEstZeeEnrichedEta,
      tausBgEstZeeEnrichedPt,
      tausBgEstZeeEnrichedLeadTrkPt,
      tausBgEstZeeEnrichedTaNCdiscr,
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

### production


elecTauPairsBgEstZeeEnriched = copy.deepcopy(elecTauPairsBgEstZtautauEnriched)
elecTauPairsBgEstZeeEnriched.srcLeg1 = cms.InputTag('electronsBgEstZeeEnrichedTrkIPcumulative')
elecTauPairsBgEstZeeEnriched.srcLeg2 = cms.InputTag('tausBgEstZeeEnrichedMuonVetoCumulative')

produceElecTauPairsBgEstZeeEnriched = cms.Sequence(elecTauPairsBgEstZeeEnriched)

### selection

elecTauPairsBgEstZeeEnrichedAntiOverlapVeto = copy.deepcopy(elecTauPairsBgEstZtautauEnrichedAntiOverlapVeto)
elecTauPairsBgEstZeeEnrichedMt1MET = copy.deepcopy(elecTauPairsBgEstZtautauEnrichedMt1MET)
elecTauPairsBgEstZeeEnrichedMt1MET.cut = cms.string('mt1MET < 500.')
elecTauPairsBgEstZeeEnrichedPzetaDiff = copy.deepcopy(elecTauPairsBgEstZtautauEnrichedPzetaDiff)
elecTauPairsBgEstZeeEnrichedPzetaDiff.cut = cms.string('(pZeta - 1.5*pZetaVis) > -2000.')
elecTauPairsBgEstZeeEnrichedZeroCharge = copy.deepcopy(elecTauPairsBgEstZtautauEnrichedZeroCharge)

elecTauPairSelConfiguratorBgEstZeeEnriched = objSelConfigurator(
    [ elecTauPairsBgEstZeeEnrichedAntiOverlapVeto,
      elecTauPairsBgEstZeeEnrichedMt1MET,
      elecTauPairsBgEstZeeEnrichedPzetaDiff,
      elecTauPairsBgEstZeeEnrichedZeroCharge
      ],
    src = "elecTauPairsBgEstZeeEnriched",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectElecTauPairsBgEstZeeEnriched = elecTauPairSelConfiguratorBgEstZeeEnriched.configure(pyNameSpace = locals())

####### anti Zee Cut #######



### production

allDiElecPairBgEstZeeEnrichedZeeHypothesesByLooseIsolation = copy.deepcopy(allDiElecPairBgEstZtautauEnrichedZeeHypothesesByLooseIsolation)
allDiElecPairBgEstZeeEnrichedZeeHypothesesByLooseIsolation.srcLeg1 = cms.InputTag("electronsBgEstZeeEnrichedTrkIPcumulative")

### selection

selectedDiElecPairBgEstZeeEnrichedZeeHypothesesByLooseIsolation = copy.deepcopy(selectedDiElecPairBgEstZtautauEnrichedZeeHypothesesByLooseIsolation)
selectedDiElecPairBgEstZeeEnrichedZeeHypothesesByLooseIsolation.src = cms.InputTag("allDiElecPairBgEstZeeEnrichedZeeHypothesesByLooseIsolation")

produceElecTauPairZeeHypothesesBgEstZeeEnriched = cms.Sequence(
	selectedPatElectronsForZeeHypotheses * 
	allDiElecPairBgEstZeeEnrichedZeeHypothesesByLooseIsolation *
	selectedDiElecPairBgEstZeeEnrichedZeeHypothesesByLooseIsolation
)


#--------------------------------------------------------------------------------  
# produce boolean event selection flags
#--------------------------------------------------------------------------------


# electron cuts

cfgElectronIdCutBgEstZeeEnriched = copy.deepcopy(cfgElectronIdCutBgEstZtautauEnriched)
cfgElectronIdCutBgEstZeeEnriched.pluginName = cms.string('electronIdCutBgEstZeeEnriched')
cfgElectronIdCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedIdCumulative')

cfgElectronAntiCrackCutBgEstZeeEnriched = copy.deepcopy(cfgElectronAntiCrackCutBgEstZtautauEnriched)
cfgElectronAntiCrackCutBgEstZeeEnriched.pluginName = cms.string('electronAntiCrackCutBgEstZeeEnriched')
cfgElectronAntiCrackCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedAntiCrackCutCumulative')

cfgElectronEtaCutBgEstZeeEnriched = copy.deepcopy(cfgElectronEtaCutBgEstZtautauEnriched)
cfgElectronEtaCutBgEstZeeEnriched.pluginName = cms.string('electronEtaCutBgEstZeeEnriched')
cfgElectronEtaCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedEtaCumulative')

cfgElectronPtCutBgEstZeeEnriched = copy.deepcopy(cfgElectronPtCutBgEstZtautauEnriched)
cfgElectronPtCutBgEstZeeEnriched.pluginName = cms.string('electronPtCutBgEstZeeEnriched')
cfgElectronPtCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedPtCumulative')

cfgElectronIsoCutBgEstZeeEnriched = copy.deepcopy(cfgElectronIsoCutBgEstZtautauEnriched)
cfgElectronIsoCutBgEstZeeEnriched.pluginName = cms.string('electronIsoCutBgEstZeeEnriched')
cfgElectronIsoCutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedIsoCumulative')

cfgElectronConversionVetoBgEstZeeEnriched = copy.deepcopy(cfgElectronConversionVetoBgEstZtautauEnriched)
cfgElectronConversionVetoBgEstZeeEnriched.pluginName = cms.string('electronConversionVetoBgEstZeeEnriched')
cfgElectronConversionVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedConversionVetoCumulative')

cfgElectronTrkIPcutBgEstZeeEnriched = copy.deepcopy(cfgElectronTrkIPcutBgEstZtautauEnriched)
cfgElectronTrkIPcutBgEstZeeEnriched.pluginName = cms.string('electronTrkIPcutBgEstZeeEnriched')
cfgElectronTrkIPcutBgEstZeeEnriched.src_cumulative = cms.InputTag('electronsBgEstZeeEnrichedTrkIPcumulative')


# tau cuts

cfgTauAntiOverlapWithElectronsVetoBgEstZeeEnriched = copy.deepcopy(cfgTauAntiOverlapWithElectronsVetoBgEstZtautauEnriched)
cfgTauAntiOverlapWithElectronsVetoBgEstZeeEnriched.pluginName = cms.string('tauAntiOverlapWithElectronsVetoBgEstZeeEnriched')
cfgTauAntiOverlapWithElectronsVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedAntiOverlapWithElectronsVetoCumulative')

cfgTauEtaCutBgEstZeeEnriched = copy.deepcopy(cfgTauEtaCutBgEstZtautauEnriched)
cfgTauEtaCutBgEstZeeEnriched.pluginName = cms.string('tauEtaCutBgEstZeeEnriched')
cfgTauEtaCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedEtaCumulative')

cfgTauPtCutBgEstZeeEnriched = copy.deepcopy(cfgTauPtCutBgEstZtautauEnriched)
cfgTauPtCutBgEstZeeEnriched.pluginName = cms.string('tauPtCutBgEstZeeEnriched')
cfgTauPtCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedPtCumulative')

## cfgTauLeadTrkCutBgEstZeeEnriched = copy.deepcopy(cfgTauLeadTrkCutBgEstZtautauEnriched)
## cfgTauLeadTrkCutBgEstZeeEnriched.pluginName = cms.string('tauLeadTrkCutBgEstZeeEnriched')
## cfgTauLeadTrkCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedLeadTrkCumulative')

cfgTauLeadTrkPtCutBgEstZeeEnriched = copy.deepcopy(cfgTauLeadTrkPtCutBgEstZtautauEnriched)
cfgTauLeadTrkPtCutBgEstZeeEnriched.pluginName = cms.string('tauLeadTrkPtCutBgEstZeeEnriched')
cfgTauLeadTrkPtCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedLeadTrkPtCumulative')

cfgTauTaNCdiscrCutBgEstZeeEnriched = copy.deepcopy(cfgTauTaNCdiscrCutBgEstZtautauEnriched)
cfgTauTaNCdiscrCutBgEstZeeEnriched.pluginName = cms.string('tauTaNCdiscrCutBgEstZeeEnriched')
cfgTauTaNCdiscrCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedTaNCdiscrCumulative')

cfgTauProngCutBgEstZeeEnriched = copy.deepcopy(cfgTauProngCutBgEstZtautauEnriched)
cfgTauProngCutBgEstZeeEnriched.pluginName = cms.string('tauProngCutBgEstZeeEnriched')
cfgTauProngCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedProngCumulative')

cfgTauChargeCutBgEstZeeEnriched = copy.deepcopy(cfgTauChargeCutBgEstZtautauEnriched)
cfgTauChargeCutBgEstZeeEnriched.pluginName = cms.string('tauChargeCutBgEstZeeEnriched')
cfgTauChargeCutBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedChargeCumulative')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! inverted cut
cfgTauElectronVetoBgEstZeeEnriched = copy.deepcopy(cfgTauElectronVetoBgEstZtautauEnriched)
cfgTauElectronVetoBgEstZeeEnriched.pluginName = cms.string('tauElectronVetoBgEstZeeEnriched')
cfgTauElectronVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedElectronVetoCumulative')

cfgTauEcalCrackVetoBgEstZeeEnriched = copy.deepcopy(cfgTauEcalCrackVetoBgEstZtautauEnriched)
cfgTauEcalCrackVetoBgEstZeeEnriched.pluginName = cms.string('tauEcalCrackVetoBgEstZeeEnriched')
cfgTauEcalCrackVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedEcalCrackVetoCumulative')

cfgTauMuonVetoBgEstZeeEnriched = copy.deepcopy(cfgTauMuonVetoBgEstZtautauEnriched)
cfgTauMuonVetoBgEstZeeEnriched.pluginName = cms.string('tauMuonVetoBgEstZeeEnriched')
cfgTauMuonVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('tausBgEstZeeEnrichedMuonVetoCumulative')

# di-tau candidate selection

cfgDiTauCandidateAntiOverlapVetoBgEstZeeEnriched = copy.deepcopy(cfgDiTauCandidateAntiOverlapVetoBgEstZtautauEnriched)
cfgDiTauCandidateAntiOverlapVetoBgEstZeeEnriched.pluginName = cms.string('diTauCandidateAntiOverlapVetoBgEstZeeEnriched')
cfgDiTauCandidateAntiOverlapVetoBgEstZeeEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZeeEnrichedAntiOverlapVetoCumulative')

cfgDiTauCandidateForElecTauMt1METCutBgEstZeeEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauMt1METCutBgEstZtautauEnriched)
cfgDiTauCandidateForElecTauMt1METCutBgEstZeeEnriched.pluginName = cms.string('diTauCandidateMt1METCutBgEstZeeEnriched')
cfgDiTauCandidateForElecTauMt1METCutBgEstZeeEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZeeEnrichedMt1METcumulative')

cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZeeEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZtautauEnriched)
cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZeeEnriched.pluginName = cms.string('diTauCandidatePzetaDiffCutBgEstZeeEnriched')
cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZeeEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZeeEnrichedPzetaDiffCumulative')

cfgDiTauCandidateForElecTauZeroChargeCutBgEstZeeEnriched = copy.deepcopy(cfgDiTauCandidateForElecTauZeroChargeCutBgEstZtautauEnriched)
cfgDiTauCandidateForElecTauZeroChargeCutBgEstZeeEnriched.pluginName = cms.string('diTauCandidateZeroChargeCutBgEstZeeEnriched')
cfgDiTauCandidateForElecTauZeroChargeCutBgEstZeeEnriched.src_cumulative = cms.InputTag('elecTauPairsBgEstZeeEnrichedZeroChargeCumulative')

cfgElecTauPairZeeHypothesisVetoBgEstZeeEnriched = copy.deepcopy(cfgDiElecPairZeeHypothesisVetoByLooseIsolation)
cfgElecTauPairZeeHypothesisVetoBgEstZeeEnriched.pluginName = cms.string('diTauCandidateZeeHypothesisVetoBgEstZeeEnriched')
cfgElecTauPairZeeHypothesisVetoBgEstZeeEnriched.src = cms.InputTag('selectedDiElecPairBgEstZeeEnrichedZeeHypothesesByLooseIsolation')


evtSelConfiguratorBgEstZeeEnriched = eventSelFlagProdConfigurator(
    [ cfgElectronIdCutBgEstZeeEnriched,
      cfgElectronAntiCrackCutBgEstZeeEnriched,
      cfgElectronEtaCutBgEstZeeEnriched,
      cfgElectronPtCutBgEstZeeEnriched,
      cfgElectronIsoCutBgEstZeeEnriched,
      cfgElectronConversionVetoBgEstZeeEnriched,
      cfgElectronTrkIPcutBgEstZeeEnriched,
      cfgTauAntiOverlapWithElectronsVetoBgEstZeeEnriched,
      cfgTauEtaCutBgEstZeeEnriched,
      cfgTauPtCutBgEstZeeEnriched,
      cfgTauLeadTrkPtCutBgEstZeeEnriched,
      cfgTauTaNCdiscrCutBgEstZeeEnriched,
      cfgTauProngCutBgEstZeeEnriched,
      cfgTauChargeCutBgEstZeeEnriched,
      cfgTauElectronVetoBgEstZeeEnriched,
      cfgTauEcalCrackVetoBgEstZeeEnriched,
      cfgTauMuonVetoBgEstZeeEnriched,
      cfgDiTauCandidateAntiOverlapVetoBgEstZeeEnriched,
      cfgDiTauCandidateForElecTauMt1METCutBgEstZeeEnriched,
      cfgDiTauCandidateForElecTauPzetaDiffCutBgEstZeeEnriched,
      cfgDiTauCandidateForElecTauZeroChargeCutBgEstZeeEnriched,
      cfgElecTauPairZeeHypothesisVetoBgEstZeeEnriched
      ],
    boolEventSelFlagProducer = "BoolEventSelFlagProducer",
    pyModuleName = __name__
)

selectEventsBgEstZeeEnriched = evtSelConfiguratorBgEstZeeEnriched.configure()


#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoElecTau_cff import *


diTauCandidateHistManagerForElecTauBgEstZeeEnriched = copy.deepcopy(diTauCandidateHistManagerForElecTauBgEstZtautauEnriched)
diTauCandidateHistManagerForElecTauBgEstZeeEnriched.pluginName = cms.string('diTauCandidateHistManagerForElecTauBgEstZeeEnriched')
diTauCandidateHistManagerForElecTauBgEstZeeEnriched.diTauCandidateSource = cms.InputTag('elecTauPairsBgEstZeeEnriched')

electronHistManagerForElecTauBgEstZeeEnriched = copy.deepcopy(electronHistManagerForElecTauBgEstZtautauEnriched)
electronHistManagerForElecTauBgEstZeeEnriched.pluginName = cms.string('electronHistManagerForElecTauBgEstZeeEnriched')

tauHistManagerForElecTauBgEstZeeEnriched = copy.deepcopy(tauHistManagerForElecTauBgEstZtautauEnriched)
tauHistManagerForElecTauBgEstZeeEnriched.pluginName = cms.string('tauHistManagerForElecTauBgEstZeeEnriched')
tauHistManagerForElecTauBgEstZeeEnriched.jetSource = cms.InputTag('selectedPatJets')



analyzeEventsBgEstZeeEnriched = cms.EDAnalyzer("GenericAnalyzer",

    name = cms.string('BgEstTemplateAnalyzer_ZeeEnriched'), 

    filters = cms.VPSet(
#        evtSelGenPhaseSpace,
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
            pluginName = cms.string('electronIsoCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronIsoCutBgEstZeeEnriched','cumulative')
        ),        
        cms.PSet(
            pluginName = cms.string('electronConversionVetoBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronConversionVetoBgEstZeeEnriched','cumulative')
        ),        
        cms.PSet(
            pluginName = cms.string('electronTrkIPcutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('electronTrkIPcutBgEstZeeEnriched','cumulative')
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
            pluginName = cms.string('diTauCandidateMt1METCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateMt1METCutBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('diTauCandidatePzetaDiffCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidatePzetaDiffCutBgEstZeeEnriched','cumulative')
        ),
        cms.PSet(
            pluginName = cms.string('diTauCandidateZeroChargeCutBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateZeroChargeCutBgEstZeeEnriched','cumulative')
        ),        
        cms.PSet(
            pluginName = cms.string('diTauCandidateZeeHypothesisVetoBgEstZeeEnriched'),
            pluginType = cms.string('BoolEventSelector'),
            src = cms.InputTag('diTauCandidateZeeHypothesisVetoBgEstZeeEnriched')
        )
    ),
  
    analyzers = cms.VPSet(
         diTauCandidateHistManagerForElecTauBgEstZeeEnriched,
         electronHistManagerForElecTauBgEstZeeEnriched,
         tauHistManagerForElecTauBgEstZeeEnriched,
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
            title = cms.string('Pt(Electron) > 15 GeV'),
        ),
        cms.PSet(
            filter = cms.string('tauAntiOverlapWithElectronsVetoBgEstZeeEnriched'),
            title = cms.string('Tau not overlapping with Elec.'),
        ),
        cms.PSet(
            filter = cms.string('tauEtaCutBgEstZeeEnriched'),
            title = cms.string('-2.3 < eta(Tau) < +2.3'),
        ),
        cms.PSet(
            filter = cms.string('tauPtCutBgEstZeeEnriched'),
            title = cms.string('Pt(Tau) > 18 GeV'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    ),
            replace = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched.electronSource = electronsBgEstZeeEnrichedPtCumulative',
                                  'tauHistManagerForElecTauBgEstZeeEnriched.tauSource = tausBgEstZeeEnrichedPtCumulative'
                                  )
        ),           
        cms.PSet(
            filter = cms.string('electronIsoCutBgEstZeeEnriched'),
            title = cms.string('Electron Isolation'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    ),
            replace = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched.electronSource = electronsBgEstZeeEnrichedIsoCumulative')
        ),          
        cms.PSet(
            filter = cms.string('electronConversionVetoBgEstZeeEnriched'),
            title = cms.string('Electron Track conversion veto'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    ),
            replace = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched.electronSource = electronsBgEstZeeEnrichedConversionVetoCumulative')
        ),          
        cms.PSet(
            filter = cms.string('electronTrkIPcutBgEstZeeEnriched'),
            title = cms.string('Electron Track IP'),
        ),        
        cms.PSet(
            filter = cms.string('tauLeadTrkPtCutBgEstZeeEnriched'),
            title = cms.string('Tau lead. Track Pt'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    ),
            replace = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched.electronSource = electronsBgEstZeeEnrichedTrkIPcumulative',
                                  'tauHistManagerForElecTauBgEstZeeEnriched.tauSource = tausBgEstZeeEnrichedLeadTrkPtCumulative')
        ),           
##         cms.PSet(
##             filter = cms.string('tauTaNCdiscrCutBgEstZeeEnriched'),
##             title = cms.string('Tau TaNC by HPS Loose (off)'),
##         ),
        cms.PSet(
            filter = cms.string('tauTaNCdiscrCutBgEstZeeEnriched'),
            title = cms.string('Tau TaNC by 0.5% (off)'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    ),
            replace = cms.vstring('tauHistManagerForElecTauBgEstZeeEnriched.tauSource = tausBgEstZeeEnrichedTaNCdiscrCumulative')
        ),         
        cms.PSet(
            filter = cms.string('tauProngCutBgEstZeeEnriched'),
            title = cms.string('Tau 1||3-Prong (off)'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    ),
            replace = cms.vstring('tauHistManagerForElecTauBgEstZeeEnriched.tauSource = tausBgEstZeeEnrichedProngCumulative')
        ),            
        cms.PSet(
            filter = cms.string('tauChargeCutBgEstZeeEnriched'),
            title = cms.string('Charge(Tau) = +/-1 (off)'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    ),
            replace = cms.vstring('tauHistManagerForElecTauBgEstZeeEnriched.tauSource = tausBgEstZeeEnrichedChargeCumulative')
        ),          
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!! inverted cut
        cms.PSet(
            filter = cms.string('tauElectronVetoBgEstZeeEnriched'),
            title = cms.string('Tau e-Veto (inverted)'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    ),
            replace = cms.vstring('tauHistManagerForElecTauBgEstZeeEnriched.tauSource = tausBgEstZeeEnrichedElectronVetoCumulative')
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
            title = cms.string('dR(Electron-Tau) > 0.5'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    'diTauCandidateHistManagerForElecTauBgEstZeeEnriched',
                                    ),
            replace = cms.vstring('tauHistManagerForElecTauBgEstZeeEnriched.tauSource = tausBgEstZeeEnrichedMuonVetoCumulative',
                                  'diTauCandidateHistManagerForElecTauBgEstZeeEnriched.diTauCandidate = elecTauPairsBgEstZeeEnrichedAntiOverlapVetoCumulative'
                                  )
        ),         
        cms.PSet(
            filter = cms.string('diTauCandidateMt1METCutBgEstZeeEnriched'),
            title = cms.string('M_{T}(Electron-MET) < 40 GeV (off)'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    'diTauCandidateHistManagerForElecTauBgEstZeeEnriched',
                                    ),
        ),        
        cms.PSet(
            filter = cms.string('diTauCandidatePzetaDiffCutBgEstZeeEnriched'),
            title = cms.string('Pzeta-1.5*Pzeta(vis) > -20 GeV (off)'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    'diTauCandidateHistManagerForElecTauBgEstZeeEnriched',
                                    ),
        ),            
        cms.PSet(
            filter = cms.string('diTauCandidateZeroChargeCutBgEstZeeEnriched'),
            title = cms.string('Charge(Electron+Tau) = 0'),
        ),
        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    'diTauCandidateHistManagerForElecTauBgEstZeeEnriched',
                                    ),
        ),          
        cms.PSet(
            filter = cms.string('diTauCandidateZeeHypothesisVetoBgEstZeeEnriched'),
            title = cms.string('no 2nd OS, loosely-isolated electron')
         ),

        cms.PSet(
            analyzers = cms.vstring('electronHistManagerForElecTauBgEstZeeEnriched',
                                    'tauHistManagerForElecTauBgEstZeeEnriched',
                                    'diTauCandidateHistManagerForElecTauBgEstZeeEnriched',
                                    ),
        ),     

    )

)

saveEvents = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test_output.root')
)




#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstZeeEnrichedAnalysisSequence = cms.Sequence(


      selectElectronsBgEstZeeEnriched
    + selectTausBgEstZeeEnriched
    + produceElecTauPairsBgEstZeeEnriched
    + selectElecTauPairsBgEstZeeEnriched
    + produceElecTauPairZeeHypothesesBgEstZeeEnriched
    + selectEventsBgEstZeeEnriched 
    + analyzeEventsBgEstZeeEnriched

)
