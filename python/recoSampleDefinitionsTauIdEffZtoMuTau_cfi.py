import FWCore.ParameterSet.Config as cms
import copy

# define configuration parameters for submission of Z --> mu + tau-jet jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)

intLumiZtoMuTau_Data = float(200.)

#--------------------------------------------------------------------------------
# Z --> tau+ tau- sample generated with Pythia + Tauola (all decay modes)
#  integrated luminosity = 1135 pb^-1
# (corrected by scale factor of 1.28 for missing files)
#
intLumiZtoMuTau_Ztautau = float(1135.4)
corrFactorZtoMuTau_Ztautau = float(1.28)

fileNamesZtoMuTau_Ztautau_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_Ztautau_part01.root'
)

fileNamesZtoMuTau_Ztautau_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_Ztautau_part02.root'
)

genPhaseSpaceCutZtoMuTau_Ztautau = cms.PSet(
    pluginName = cms.string('genPhaseSpaceCut'),
    pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
    src = cms.InputTag('genPhaseSpaceEventInfo'),
    cut = cms.string('')
)

plotsOutputFileNameZtoMuTau_Ztautau = cms.string('plotsZtoMuTau_Ztautau_partXX.root')
patTupleOutputFileNameZtoMuTau_Ztautau = cms.untracked.string('patTupleZtoMuTau_Ztautau_partXX.root')
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# Z --> mu+ mu- sample generated with Pythia
#  integrated luminosity = 633 pb^-1
# (corrected by scale factor of 1.49 for missing files)
#
intLumiZtoMuTau_Zmumu = float(633.)
corrFactorZtoMuTau_Zmumu = float(1.49)

fileNamesZtoMuTau_Zmumu_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_Zmumu_part01.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_Zmumu_part02.root'
)

fileNamesZtoMuTau_Zmumu_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_Zmumu_part03.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_Zmumu_part04.root'
)

fileNamesZtoMuTau_Zmumu_part03 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_Zmumu_part05.root'
)

genPhaseSpaceCutZtoMuTau_Zmumu = cms.PSet(
    pluginName = cms.string('genPhaseSpaceCut'),
    pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
    src = cms.InputTag('genPhaseSpaceEventInfo'),
    cut = cms.string('')
)

plotsOutputFileNameZtoMuTau_Zmumu = cms.string('plotsZtoMuTau_Zmumu_partXX.root')
patTupleOutputFileNameZtoMuTau_Zmumu = cms.untracked.string('patTupleZtoMuTau_Zmumu_partXX.root')
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# Z + jets sample generated with Madgraph
#  integrated luminosity = 341 pb^-1
# (to be corrected for missing files)
#
# (NOTE: for Monte Carlo samples generated by Madgraph,
#        the filter efficiency is already included in the cross-sections
#        listed at https://twiki.cern.ch/twiki/bin/view/CMS/ProductionSummer2008 !!)
#
intLumiZtoMuTau_ZplusJets = float(341.)
corrFactorZtoMuTau_ZplusJets = float(1.16)

fileNamesZtoMuTau_ZeePlusJets = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_ZeePlusJets_part01.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_ZeePlusJets_part02.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_ZeePlusJets_part03.root'
)

fileNamesZtoMuTau_ZmumuPlusJets = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_ZmumuPlusJets_part01.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_ZmumuPlusJets_part02.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_ZmumuPlusJets_part03.root'
)

fileNamesZtoMuTau_ZtautauPlusJets = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_ZtautauPlusJets_part01.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_ZtautauPlusJets_part02.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_ZtautauPlusJets_part03.root'
)

corrFactorZtoMuTau_ZeePlusJets = corrFactorZtoMuTau_ZplusJets
intLumiZtoMuTau_ZeePlusJets = intLumiZtoMuTau_ZplusJets

plotsOutputFileNameZtoMuTau_ZeePlusJets = cms.string('plotsZtoMuTau_ZeePlusJets.root')
patTupleOutputFileNameZtoMuTau_ZeePlusJets = cms.untracked.string('patTupleZtoMuTau_ZeePlusJets.root')

genPhaseSpaceCutZtoMuTau_ZeePlusJets = cms.PSet(
    pluginName = cms.string('genPhaseSpaceCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('genElectronsFromZs'),
    minNumber = cms.uint32(2)
)

corrFactorZtoMuTau_ZmumuPlusJets = corrFactorZtoMuTau_ZplusJets
intLumiZtoMuTau_ZmumuPlusJets = intLumiZtoMuTau_ZplusJets

plotsOutputFileNameZtoMuTau_ZmumuPlusJets = cms.string('plotsZtoMuTau_ZmumuPlusJets.root')
patTupleOutputFileNameZtoMuTau_ZmumuPlusJets = cms.untracked.string('patTupleZtoMuTau_ZmumuPlusJets.root')

genPhaseSpaceCutZtoMuTau_ZmumuPlusJets = cms.PSet(
    pluginName = cms.string('genPhaseSpaceCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('genMuonsFromZs'),
    minNumber = cms.uint32(2)
)

corrFactorZtoMuTau_ZtautauPlusJets = corrFactorZtoMuTau_ZplusJets
intLumiZtoMuTau_ZtautauPlusJets = intLumiZtoMuTau_ZplusJets

plotsOutputFileNameZtoMuTau_ZtautauPlusJets = cms.string('plotsZtoMuTau_ZtautauPlusJets.root')
patTupleOutputFileNameZtoMuTau_ZtautauPlusJets = cms.untracked.string('patTupleZtoMuTau_ZtautauPlusJets.root')

genPhaseSpaceCutZtoMuTau_ZtautauPlusJets = cms.PSet(
    pluginName = cms.string('genPhaseSpaceCut'),
    pluginType = cms.string('PATCandViewMinEventSelector'),
    src = cms.InputTag('genTausFromZs'),
    minNumber = cms.uint32(2)
)
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# W + jets sample generated with Madgraph
#  integrated luminosity = 242 pb^-1
# (corrected by scale factor of 1.36 for missing files)
#
# (NOTE: for Monte Carlo samples generated by Madgraph,
#        the filter efficiency is already included in the cross-sections
#        listed at https://twiki.cern.ch/twiki/bin/view/CMS/ProductionSummer2008 !!)
#
intLumiZtoMuTau_WplusJets = float(242.)
corrFactorZtoMuTau_WplusJets = float(1.36)

fileNamesZtoMuTau_WplusJets_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part01.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part02.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part03.root'
)

fileNamesZtoMuTau_WplusJets_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part04.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part05.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part06.root'
)

fileNamesZtoMuTau_WplusJets_part03 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part07.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part08.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part09.root'
)

fileNamesZtoMuTau_WplusJets_part04 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part10.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_WplusJets_part11.root'
)

genPhaseSpaceCutZtoMuTau_WplusJets = cms.PSet(
    pluginName = cms.string('genPhaseSpaceCut'),
    pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
    src = cms.InputTag('genPhaseSpaceEventInfo'),
    cut = cms.string('')
)

plotsOutputFileNameZtoMuTau_WplusJets = cms.string('plotsZtoMuTau_WplusJets_partXX.root')
patTupleOutputFileNameZtoMuTau_WplusJets = cms.untracked.string('patTupleZtoMuTau_WplusJets_partXX.root')
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# muon enriched QCD sample generated with Pythia (no cut on Pt(hat))
#  integrated luminosity = 0.044 pb^-1
# (corrected by scale factor of 1.22 for missing files)
#
intLumiZtoMuTau_InclusivePPmuX = float(0.044)
corrFactorZtoMuTau_InclusivePPmuX = float(1.22)

fileNamesZtoMuTau_InclusivePPmuX = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_InclusivePPmuX_part01.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_InclusivePPmuX_part02.root'
)

genPhaseSpaceCutZtoMuTau_InclusivePPmuX = cms.PSet(
    pluginName = cms.string('genPhaseSpaceCut'),
    pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
    src = cms.InputTag('genPhaseSpaceEventInfo'),
    cut = cms.string('ptHat < 20. | leadingGenMuon.pt < 15.')
)

plotsOutputFileNameZtoMuTau_InclusivePPmuX = cms.string('plotsZtoMuTau_InclusivePPmuX.root')
patTupleOutputFileNameZtoMuTau_InclusivePPmuX = cms.untracked.string('patTupleZtoMuTau_InclusivePPmuX.root')
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# muon enriched QCD sample generated with Pythia (Pt(hat)> 20 GeV && PtMuon > 15 GeV)
#  integrated luminosity = 42 pb^-1
# (corrected by estimated scale factor of 1.1 for missing files)
#
intLumiZtoMuTau_PPmuXptGt20 = float(42.0)
corrFactorZtoMuTau_PPmuXptGt20 = float(1.1)

fileNamesZtoMuTau_PPmuXptGt20_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part01.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part02.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part03.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part04.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part05.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part06.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part07.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part08.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part09.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part10.root'
)

fileNamesZtoMuTau_PPmuXptGt20_part02 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part11.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part12.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part13.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part14.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part15.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part16.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part17.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part18.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part19.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part20.root'
)

fileNamesZtoMuTau_PPmuXptGt20_part03 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part21.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part22.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part23.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part24.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part25.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part26.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part27.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part28.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part29.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part30.root'
)

fileNamesZtoMuTau_PPmuXptGt20_part04 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part31.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part32.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part33.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part34.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part35.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part36.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part37.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part38.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part39.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part40.root'
)

fileNamesZtoMuTau_PPmuXptGt20_part05 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part41.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part42.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part43.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part44.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part45.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part46.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part47.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part48.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part49.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part50.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_PPmuXptGt20_part51.root'
)

genPhaseSpaceCutZtoMuTau_PPmuXptGt20 = cms.PSet(
    pluginName = cms.string('genPhaseSpaceCut'),
    pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
    src = cms.InputTag('genPhaseSpaceEventInfo'),
    cut = cms.string('ptHat > 20.')
)

plotsOutputFileNameZtoMuTau_PPmuXptGt20 = cms.string('plotsZtoMuTau_PPmuXptGt20_partXX.root')
patTupleOutputFileNameZtoMuTau_PPmuXptGt20 = cms.untracked.string('patTupleZtoMuTau_PPmuXptGt20_partXX.root')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# TTbar sample generated with Madgraph
#
intLumiZtoMuTau_TTplusJets = float(2986)
corrFactorZtoMuTau_TTplusJets = float(1.0)

fileNamesZtoMuTau_TTplusJets_part01 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part01.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part02.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part03.root'
)

fileNamesZtoMuTau_TTplusJets_part02 = cms.untracked.vstring(    
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part04.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part05.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part06.root'
)

fileNamesZtoMuTau_TTplusJets_part03 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part07.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part08.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part09.root'
)

fileNamesZtoMuTau_TTplusJets_part04 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part10.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part11.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part12.root'
)

fileNamesZtoMuTau_TTplusJets_part05 = cms.untracked.vstring(
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part13.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part14.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part15.root',
    'rfio:/castor/cern.ch/user/v/veelken/tauIdEff/selEvents_ZtoMuTau_TTplusJets_part16.root'
)

genPhaseSpaceCutZtoMuTau_TTplusJets = cms.PSet(
    pluginName = cms.string('genPhaseSpaceCut'),
    pluginType = cms.string('GenPhaseSpaceEventInfoSelector'),
    src = cms.InputTag('genPhaseSpaceEventInfo'),
    cut = cms.string('')
)

plotsOutputFileNameZtoMuTau_TTplusJets = cms.string('plotsZtoMuTau_TTplusJets_partXX.root')
patTupleOutputFileNameZtoMuTau_TTplusJets = cms.untracked.string('patTupleZtoMuTau_TTplusJets_partXX.root')
#--------------------------------------------------------------------------------
