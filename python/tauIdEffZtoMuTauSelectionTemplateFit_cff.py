import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *
from TauAnalysis.RecoTools.tools.eventSelFlagProdConfigurator import *

#--------------------------------------------------------------------------------
# select loosely selected Z --> tau+ tau- --> muon + tau-jet event sample
# for measuring tau id. efficiency
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------  
# produce collection of pat::Muons
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patMuonSelection_cfi import *

muonsForTauIdEffZtoMuTauTemplateFitPFRelIso = cms.EDFilter("PATMuonPFIsolationSelector",
    pfCandidateSource = cms.InputTag('pfNoPileUp'),
    chargedHadronIso = cms.PSet(
        ptMin = cms.double(0.5),        
        dRvetoCone = cms.double(-1.),
        dRisoCone = cms.double(0.6)
    ),
    neutralHadronIso = cms.PSet(
        ptMin = cms.double(1000.),        
        dRvetoCone = cms.double(0.08),        
        dRisoCone = cms.double(0.6)
    ),
    photonIso = cms.PSet(
        ptMin = cms.double(0.5),        
        dPhiVeto = cms.double(-1.),  # asymmetric Eta x Phi veto region 
        dEtaVeto = cms.double(-1.),  # to account for photon conversions in electron isolation case        
        dRvetoCone = cms.double(-1.),
        dRisoCone = cms.double(0.6)
    ),
    sumPtMax = cms.double(0.10),
    sumPtMethod = cms.string("relative"), # either "relative" or "absolute"
    filter = cms.bool(False)                                                        
)

muonsForTauIdEffZtoMuTauTemplateFitTrk = copy.deepcopy(selectedPatMuonsTrk)

muonsForTauIdEffZtoMuTauTemplateFitTrkIP = copy.deepcopy(selectedPatMuonsTrkIP)

muonSelConfiguratorTauIdEffZtoMuTauTemplateFit = objSelConfigurator(
    [ muonsForTauIdEffZtoMuTauTemplateFitPFRelIso,
      muonsForTauIdEffZtoMuTauTemplateFitTrk,
      muonsForTauIdEffZtoMuTauTemplateFitTrkIP ],
    src = "selectedPatMuonsPt10Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuonsForTauIdEffZtoMuTauTemplateFit = muonSelConfiguratorTauIdEffZtoMuTauTemplateFit.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of pat::Taus
#--------------------------------------------------------------------------------

from TauAnalysis.RecoTools.patPFTauSelection_cfi import *
from TauAnalysis.RecoTools.patPFTauSelectionForMuTau_cfi import *

tausForTauIdEffZtoMuTauTemplateFitMuonVeto = copy.deepcopy(selectedPatTausForMuTauMuonVeto)

tausForTauIdEffZtoMuTauTemplateFitCaloMuonVeto = copy.deepcopy(selectedPatTausForMuTauCaloMuonVeto)

tausForTauIdEffZtoMuTauTemplateFitPFRelIsoLoose = cms.EDFilter("PATTauPFIsolationSelector",
    pfCandidateSource = cms.InputTag('pfNoPileUp'),
    chargedHadronIso = cms.PSet(
        ptMin = cms.double(1.0),        
        dRvetoCone = cms.double(-1.),
        dRisoCone = cms.double(0.6)
    ),
    neutralHadronIso = cms.PSet(
        ptMin = cms.double(1000.),        
        dRvetoCone = cms.double(0.08),        
        dRisoCone = cms.double(0.6)
    ),
    photonIso = cms.PSet(
        ptMin = cms.double(1.5),        
        dPhiVeto = cms.double(-1.),  # asymmetric Eta x Phi veto region 
        dEtaVeto = cms.double(-1.),  # to account for photon conversions in electron isolation case        
        dRvetoCone = cms.double(-1.),
        dRisoCone = cms.double(0.6)
    ),
    sumPtMax = cms.double(2.5),
    sumPtMethod = cms.string("absolute"), # either "relative" or "absolute"
    filter = cms.bool(False)                                                        
)

tauSelConfiguratorTauIdEffZtoMuTauTemplateFit = objSelConfigurator(
    [ tausForTauIdEffZtoMuTauTemplateFitMuonVeto,
      tausForTauIdEffZtoMuTauTemplateFitCaloMuonVeto,
      tausForTauIdEffZtoMuTauTemplateFitPFRelIsoLoose ],
    src = "selectedPatTausForMuTauLeadTrkPtCumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectTausForTauIdEffZtoMuTauTemplateFit = tauSelConfiguratorTauIdEffZtoMuTauTemplateFit.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# produce collection of Z --> muon+ muon- hypotheses
#--------------------------------------------------------------------------------

# require muon candidates considered for Z --> mu+ mu- hypothesis
# to be reconstructed in muon system
# (with or without a track reconstructed in Pixel/SiStrip tracking detectors linked to it)
muonsLooseForZmumuHypothesesTauIdEffZtoMuTauTemplateFit = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string('isGlobalMuon | isStandAloneMuon | pt > 5.'),
    filter = cms.bool(False)
)

muonsTightForZmumuHypothesesTauIdEffZtoMuTauTemplateFit = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string('isGlobalMuon | isStandAloneMuon'),
    filter = cms.bool(False)
)

allDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit = cms.EDProducer("DiCandidatePairProducer",
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsTightForZmumuHypothesesTauIdEffZtoMuTauTemplateFit'),
    srcLeg2 = cms.InputTag('muonsLooseForZmumuHypothesesTauIdEffZtoMuTauTemplateFit'),
    dRmin12 = cms.double(0.5),
    srcMET = cms.InputTag(''),
    recoMode = cms.string(""),
    scaleFuncImprovedCollinearApprox = cms.string('1'),                                        
    verbosity = cms.untracked.int32(0)
)

selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit = cms.EDFilter("DiCandidatePairSelector",
    src = cms.InputTag("allDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit"),                                   
    cut = cms.string('p4Vis.mass > 70. & p4Vis.mass < 110.'),
    filter = cms.bool(False)
)

produceDiMuPairsTauIdEffZtoMuTauTemplateFit = cms.Sequence(
    muonsLooseForZmumuHypothesesTauIdEffZtoMuTauTemplateFit * muonsTightForZmumuHypothesesTauIdEffZtoMuTauTemplateFit
   * allDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit
   * selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit
)

#--------------------------------------------------------------------------------  
# produce collection of muon + tau-jet combinations
#--------------------------------------------------------------------------------

from TauAnalysis.CandidateTools.muTauPairProduction_cff import *
from TauAnalysis.CandidateTools.resolutions_cfi import *

muTauPairsTauIdEffZtoMuTauTemplateFit = allMuTauPairs.clone(
    useLeadingTausOnly = cms.bool(False),
    srcLeg1 = cms.InputTag('muonsForTauIdEffZtoMuTauTemplateFitTrkIPcumulative'),
    srcLeg2 = cms.InputTag('tausForTauIdEffZtoMuTauTemplateFitPFRelIsoLooseCumulative'),
    dRmin12 = cms.double(0.7),
    srcMET = cms.InputTag('patMETs'),
    recoMode = cms.string(""),
    doSVreco = cms.bool(False)
)

produceMuTauPairsTauIdEffZtoMuTauTemplateFit = cms.Sequence(muTauPairsTauIdEffZtoMuTauTemplateFit)

muTauPairsTauIdEffZtoMuTauTemplateFitMt1MET = cms.EDFilter("PATMuTauPairSelector",
    cut = cms.string('mt1MET < 40.'),
    filter = cms.bool(False)
)

muTauPairsTauIdEffZtoMuTauTemplateFitPzetaDiff = cms.EDFilter("PATMuTauPairSelector",
    cut = cms.string('(pZeta - 1.5*pZetaVis) > -20.'),
    filter = cms.bool(False)
)

muTauPairConfiguratorTauIdEffZtoMuTauTemplateFit = objSelConfigurator(
    [ muTauPairsTauIdEffZtoMuTauTemplateFitMt1MET,
      muTauPairsTauIdEffZtoMuTauTemplateFitPzetaDiff ],
    src = "muTauPairsTauIdEffZtoMuTauTemplateFit",
    pyModuleName = __name__,
    doSelIndividual = False
)

selectMuTauPairsTauIdEffZtoMuTauTemplateFit = muTauPairConfiguratorTauIdEffZtoMuTauTemplateFit.configure(pyNameSpace = locals())

#--------------------------------------------------------------------------------  
# apply event selection criteria; fill histograms
#--------------------------------------------------------------------------------

from TauAnalysis.Configuration.analyzeZtoMuTau_cfi import *
from TauAnalysis.BgEstimationTools.selectZtoMuTauEventVertex_cff import *

muonHistManagerTauIdEffZtoMuTauTemplateFit = copy.deepcopy(muonHistManager)
muonHistManagerTauIdEffZtoMuTauTemplateFit.pluginName = 'muonHistManagerTauIdEffZtoMuTauTemplateFit'
muonHistManagerTauIdEffZtoMuTauTemplateFit.muonSource = muTauPairsTauIdEffZtoMuTauTemplateFit.srcLeg1

tauHistManagerTauIdEffZtoMuTauTemplateFit = copy.deepcopy(tauHistManager)
tauHistManagerTauIdEffZtoMuTauTemplateFit.pluginName = 'tauHistManagerTauIdEffZtoMuTauTemplateFit'
tauHistManagerTauIdEffZtoMuTauTemplateFit.tauSource = muTauPairsTauIdEffZtoMuTauTemplateFit.srcLeg2

from TauAnalysis.Core.leptonPFIsolationHistManager_cfi import *
tauPFIsolationHistManagerZtoMuTauTemplateFit = copy.deepcopy(tauPFIsolationHistManager)
tauPFIsolationHistManagerZtoMuTauTemplateFit.pluginName = 'tauPFIsolationHistManagerZtoMuTauTemplateFit'
tauPFIsolationHistManagerZtoMuTauTemplateFit.leptonSource = muTauPairsTauIdEffZtoMuTauTemplateFit.srcLeg2
del tauPFIsolationHistManagerZtoMuTauTemplateFit.genLeptonMatch

diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit = copy.deepcopy(diTauCandidateHistManagerForMuTau)
diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit.pluginName = 'diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit'
diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit.diTauCandidateSource = 'muTauPairsTauIdEffZtoMuTauTemplateFitPzetaDiffCumulative'
diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit.visMassHypothesisSource = cms.InputTag('')

from TauAnalysis.BgEstimationTools.tauIdEffZtoMuTauHistManager_cfi import *
tauIdEffZtoMuTauHistManagerTemplateFit = copy.deepcopy(tauIdEffZtoMuTauHistManager)
tauIdEffZtoMuTauHistManagerTemplateFit.pluginName = 'tauIdEffZtoMuTauHistManagerTemplateFit'
tauIdEffZtoMuTauHistManagerTemplateFit.muonSource = muTauPairsTauIdEffZtoMuTauTemplateFit.srcLeg1
tauIdEffZtoMuTauHistManagerTemplateFit.tauSource = muTauPairsTauIdEffZtoMuTauTemplateFit.srcLeg2
tauIdEffZtoMuTauHistManagerTemplateFit.diTauSource = diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit.diTauCandidateSource
tauIdEffZtoMuTauHistManagerTemplateFit.diTauChargeSignExtractor.src = tauIdEffZtoMuTauHistManagerTemplateFit.diTauSource

from TauAnalysis.Core.caloMEtHistManager_cfi import *
caloMEtHistManagerTemplateFit = copy.deepcopy(caloMEtHistManager)
caloMEtHistManagerTemplateFit.pluginName = 'caloMEtHistManagerTemplateFit'
caloMEtHistManagerTemplateFit.leg1Source = muTauPairsTauIdEffZtoMuTauTemplateFit.srcLeg1
caloMEtHistManagerTemplateFit.leg2Source = muTauPairsTauIdEffZtoMuTauTemplateFit.srcLeg2
from TauAnalysis.Core.pfMEtHistManager_cfi import *
pfMEtHistManagerTemplateFit = copy.deepcopy(pfMEtHistManager)
pfMEtHistManagerTemplateFit.pluginName = 'pfMEtHistManagerTemplateFit'
pfMEtHistManagerTemplateFit.leg1Source = muTauPairsTauIdEffZtoMuTauTemplateFit.srcLeg1
pfMEtHistManagerTemplateFit.leg2Source = muTauPairsTauIdEffZtoMuTauTemplateFit.srcLeg2

from TauAnalysis.BgEstimationTools.bgEstBinGridZtoMuTau_cfi import *

dataBinnerTauIdEffZtoMuTauTemplateFit = copy.deepcopy(dataBinner)
dataBinnerTauIdEffZtoMuTauTemplateFit.pluginName = 'dataBinnerTauIdEffZtoMuTauTemplateFit'

binningTauIdEffZtoMuTauTemplateFit_ewkTauId = copy.deepcopy(binning_ewkTauId)
binningTauIdEffZtoMuTauTemplateFit_ewkTauId.extractor.src = muTauPairsTauIdEffZtoMuTauTemplateFit.srcLeg2

binningTauIdEffZtoMuTauTemplateFit_diTauAbsCharge = copy.deepcopy(binning_diTauAbsCharge)
binningTauIdEffZtoMuTauTemplateFit_diTauAbsCharge.extractor.src = \
  diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit.diTauCandidateSource

tauIdEffBinningZtoMuTau_template1d = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_template1d"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauTemplateFit_ewkTauId
    )
)

tauIdEffBinningZtoMuTau_template2d = cms.PSet(
    name = cms.string("tauIdEffBinningZtoMuTau_template2d"),
    config = cms.VPSet(
        binningTauIdEffZtoMuTauTemplateFit_ewkTauId,
        binningTauIdEffZtoMuTauTemplateFit_diTauAbsCharge
    )
)

analyzeEventsTauIdEffZtoMuTauTemplateFit = cms.EDAnalyzer("GenericAnalyzer",
  
    name = cms.string('TauIdEffAnalyzerZtoMuTauTemplateFit'), 
                            
    filters = cms.VPSet(
        evtSelGenPhaseSpace,
        evtSelTrigger,
        evtSelPrimaryEventVertex,
        evtSelPrimaryEventVertexQuality,
        evtSelPrimaryEventVertexPosition,
        evtSelGlobalMuon,
        evtSelMuonEta,
        evtSelMuonPt,
        cms.PSet(
            pluginName = cms.string('muonPFRelIsoCutTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauTemplateFitPFRelIsoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muonTrkIPcutTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauTemplateFitTrkCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('tauEtaCutTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muonsForTauIdEffZtoMuTauTemplateFitTrkIPcumulative'),
            minNumber = cms.uint32(1)
        ),
        evtSelTauAntiOverlapWithMuonsVeto,
        evtSelTauEta,
        evtSelTauPt,
        evtSelTauLeadTrk,
        evtSelTauLeadTrkPt,
        cms.PSet(
            pluginName = cms.string('tauMuonVetoTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('tausForTauIdEffZtoMuTauTemplateFitCaloMuonVetoCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('tauPFRelIsoLooseTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('tausForTauIdEffZtoMuTauTemplateFitPFRelIsoLooseCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauTemplateFit'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairMt1METtauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauTemplateFitMt1METcumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMinEventSelector'),
            src = cms.InputTag('muTauPairsTauIdEffZtoMuTauTemplateFitPzetaDiffCumulative'),
            minNumber = cms.uint32(1)
        ),
        cms.PSet(
            pluginName = cms.string('diMuPairZmumuHypothesisVeto'),
            pluginType = cms.string('PATCandViewMaxEventSelector'),
            src = cms.InputTag('selectedDiMuPairZmumuHypothesesTauIdEffZtoMuTauTemplateFit'),
            maxNumber = cms.uint32(0)
        ),
        cms.PSet(
            pluginName = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTauTemplateFit'),
            pluginType = cms.string('PATCandViewMaxEventSelector'),
            src = cms.InputTag('muonsTightForZmumuHypothesesTauIdEffZtoMuTauTemplateFit'),
            maxNumber = cms.uint32(1)
        )
    ),
  
    analyzers = cms.VPSet(
        muonHistManagerTauIdEffZtoMuTauTemplateFit,
        tauHistManagerTauIdEffZtoMuTauTemplateFit,
        tauPFIsolationHistManagerZtoMuTauTemplateFit,
        diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit,
        caloMEtHistManagerTemplateFit,
        pfMEtHistManagerTemplateFit,        
        tauIdEffZtoMuTauHistManagerTemplateFit,
        dataBinnerTauIdEffZtoMuTauTemplateFit,
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinner2regions'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_template1d,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResults2regions')
        ),       
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManager2regions'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_template1d,
            histManagers = cms.VPSet(
                muonHistManagerTauIdEffZtoMuTauTemplateFit,
                tauHistManagerTauIdEffZtoMuTauTemplateFit,
                tauPFIsolationHistManagerZtoMuTauTemplateFit,
                diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit,
                tauIdEffZtoMuTauHistManagerTemplateFit
            ),
            dqmDirectory_store = cms.string('tauIdEffHistograms2regions')
        ),
        cms.PSet(
            pluginName = cms.string('tauIdEffDataBinner4regions'),
            pluginType = cms.string('DataBinner'),
            binning = tauIdEffBinningZtoMuTau_template2d,
            binningService = cms.PSet(
                pluginType = cms.string("DataBinningService")
            ),
            dqmDirectory_store = cms.string('tauIdEffBinningResults4regions')
        ),       
        cms.PSet(
            pluginName = cms.string('tauIdEffBinGridHistManager4regions'),
            pluginType = cms.string('BinGridHistManager'),
            binning = tauIdEffBinningZtoMuTau_template2d,
            histManagers = cms.VPSet(
                muonHistManagerTauIdEffZtoMuTauTemplateFit,
                tauHistManagerTauIdEffZtoMuTauTemplateFit,
                tauPFIsolationHistManagerZtoMuTauTemplateFit,
                diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit,
                tauIdEffZtoMuTauHistManagerTemplateFit
            ),
            dqmDirectory_store = cms.string('tauIdEffHistograms4regions')
        )
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
            title = cms.string('Muon Trigger')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertex'),
            title = cms.string('Vertex')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertexQuality'),
            title = cms.string('Vertex quality')
        ),
        cms.PSet(
            filter = cms.string('evtSelPrimaryEventVertexPosition'),
            title = cms.string('-24 < zVertex < +24 cm')
        ),
        cms.PSet(
            filter = cms.string('evtSelGlobalMuon'),
            title = cms.string('global Muon')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonEta'),
            title = cms.string('-2.1 < eta(Muon) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('evtSelMuonPt'),
            title = cms.string('Pt(Muon) > 15 GeV')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauAntiOverlapWithMuonsVeto'),
            title = cms.string('Tau not overlapping w. Muon'),
        ),
        cms.PSet(
            filter = cms.string('evtSelTauEta'),
            title = cms.string('-2.1 < eta(Tau) < +2.1')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauPt'),
            title = cms.string('Pt(Tau) > 20 GeV')
        ),        
        cms.PSet(
            filter = cms.string('muonPFRelIsoCutTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('Muon iso.')
        ),
        cms.PSet(
            filter = cms.string('muonTrkIPcutTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('Muon Track IP')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrk'),
            title = cms.string('Tau lead. Track find.')
        ),
        cms.PSet(
            filter = cms.string('evtSelTauLeadTrkPt'),
            title = cms.string('Tau lead. Track Pt')
        ),
        cms.PSet(
            filter = cms.string('tauMuonVetoTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('Tau mu-Veto')
        ),
        cms.PSet(
            filter = cms.string('tauPFRelIsoLooseTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('Tau loose iso.')
        ),
        cms.PSet(
            filter = cms.string('muTauPairTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('dR(Muon-Tau) > 0.7')
        ),
        cms.PSet(
            filter = cms.string('muTauPairMt1METtauIdEffZtoMuTauTemplateFit'),
            title = cms.string('M_{T}(Muon-MET) < 40 GeV')
        ),
        cms.PSet(
            filter = cms.string('muTauPairPzetaDiffTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('P_{#zeta} - 1.5*P_{#zeta}^{vis} > -20 GeV')
        ),
        cms.PSet(
            filter = cms.string('diMuPairZmumuHypothesisVeto'),
            title = cms.string('not 80 < M (Muon-Muon) < 100 GeV')
        ),        
        cms.PSet(
            filter = cms.string('uniqueMuonCandidateCutTauIdEffZtoMuTauTemplateFit'),
            title = cms.string('num. global || stand-alone Muons < 2')
        ),
        cms.PSet(
            analyzers = cms.vstring(
                'muonHistManagerTauIdEffZtoMuTauTemplateFit',
                'tauHistManagerTauIdEffZtoMuTauTemplateFit',
                'tauPFIsolationHistManagerZtoMuTauTemplateFit',
                'diTauCandidateHistManagerTauIdEffZtoMuTauTemplateFit',
                'caloMEtHistManagerTemplateFit',
                'pfMEtHistManagerTemplateFit',
                'tauIdEffZtoMuTauHistManagerTemplateFit',
                'dataBinnerTauIdEffZtoMuTauTemplateFit',
                'tauIdEffDataBinner2regions',
                'tauIdEffBinGridHistManager2regions',
                'tauIdEffDataBinner4regions',
                'tauIdEffBinGridHistManager4regions'
            )
        )
    )
)

analysisSequenceTauIdEffZtoMuTauTemplateFit = cms.Sequence(analyzeEventsTauIdEffZtoMuTauTemplateFit)

#--------------------------------------------------------------------------------  
# define (final) analysis sequence
#--------------------------------------------------------------------------------

bgEstTauIdEffZtoMuTauTemplateFitAnalysisSequence = cms.Sequence(
    selectMuonsForTauIdEffZtoMuTauTemplateFit 
   + selectTausForTauIdEffZtoMuTauTemplateFit 
   + produceDiMuPairsTauIdEffZtoMuTauTemplateFit
   + produceMuTauPairsTauIdEffZtoMuTauTemplateFit + selectMuTauPairsTauIdEffZtoMuTauTemplateFit 
   + analysisSequenceTauIdEffZtoMuTauTemplateFit
)
