import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.Configuration.plotterProcessDefinitions_cfi import *
from TauAnalysis.BgEstimationTools.recoSampleDefinitionsTauIdEffZtoMuTau_cfi import *

#--------------------------------------------------------------------------------
# define for Z --> mu + tau-jet analysis names of .root files containing histograms
#--------------------------------------------------------------------------------

processZtoMuTau_Ztautau = copy.deepcopy(process_Ztautau)
processZtoMuTau_Ztautau.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_Ztautau_part01.root',
    'plotsZtoMuTau_Ztautau_part02.root'
)
processZtoMuTau_Ztautau.config_dqmFileLoader.scaleFactor = cms.double(corrFactorZtoMuTau_Ztautau*intLumiZtoMuTau_Data/intLumiZtoMuTau_Ztautau)

processZtoMuTau_ZtautauSum = copy.deepcopy(process_Ztautau)
processZtoMuTau_ZtautauSum.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_ZtautauSum.root'
)
processZtoMuTau_ZtautauSum.config_dqmFileLoader.dqmDirectory_store = cms.string('harvested')
processZtoMuTau_ZtautauSum.config_dqmFileLoader.scaleFactor = cms.double(1.)

#--------------------------------------------------------------------------------

processZtoMuTau_Zmumu = copy.deepcopy(process_Zmumu)
processZtoMuTau_Zmumu.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_Zmumu_part01.root',
    'plotsZtoMuTau_Zmumu_part02.root',
    'plotsZtoMuTau_Zmumu_part03.root'
)
processZtoMuTau_Zmumu.config_dqmFileLoader.scaleFactor = cms.double(corrFactorZtoMuTau_Zmumu*intLumiZtoMuTau_Data/intLumiZtoMuTau_Zmumu)

processZtoMuTau_ZmumuSum = copy.deepcopy(process_Zmumu)
processZtoMuTau_ZmumuSum.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_ZmumuSum.root'
)
processZtoMuTau_ZmumuSum.config_dqmFileLoader.dqmDirectory_store = cms.string('harvested')
processZtoMuTau_ZmumuSum.config_dqmFileLoader.scaleFactor = cms.double(1.)

#--------------------------------------------------------------------------------

processZtoMuTau_ZeePlusJets = copy.deepcopy(process_ZeePlusJets)
processZtoMuTau_ZeePlusJets.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_ZeePlusJets.root'
)
processZtoMuTau_ZeePlusJets.config_dqmFileLoader.dqmDirectory_store = cms.string('harvested/ZeePlusJets')
processZtoMuTau_ZeePlusJets.config_dqmFileLoader.scaleFactor = cms.double(corrFactorZtoMuTau_ZeePlusJets*intLumiZtoMuTau_Data/intLumiZtoMuTau_ZeePlusJets)

processZtoMuTau_ZmumuPlusJets = copy.deepcopy(process_ZmumuPlusJets)
processZtoMuTau_ZmumuPlusJets.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_ZmumuPlusJets.root'
)
processZtoMuTau_ZmumuPlusJets.config_dqmFileLoader.dqmDirectory_store = cms.string('harvested/ZmumuPlusJets')
processZtoMuTau_ZmumuPlusJets.config_dqmFileLoader.scaleFactor = cms.double(corrFactorZtoMuTau_ZmumuPlusJets*intLumiZtoMuTau_Data/intLumiZtoMuTau_ZmumuPlusJets)

processZtoMuTau_ZtautauPlusJets = copy.deepcopy(process_ZtautauPlusJets)
processZtoMuTau_ZtautauPlusJets.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_ZtautauPlusJets.root'
)
processZtoMuTau_ZtautauPlusJets.config_dqmFileLoader.dqmDirectory_store = cms.string('harvested/ZtautauPlusJets')
processZtoMuTau_ZtautauPlusJets.config_dqmFileLoader.scaleFactor = cms.double(corrFactorZtoMuTau_ZtautauPlusJets*intLumiZtoMuTau_Data/intLumiZtoMuTau_ZtautauPlusJets)

#--------------------------------------------------------------------------------

processZtoMuTau_WplusJets = copy.deepcopy(process_WplusJets)
processZtoMuTau_WplusJets.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_WplusJets_part01.root',
    'plotsZtoMuTau_WplusJets_part02.root',
    'plotsZtoMuTau_WplusJets_part03.root',
    'plotsZtoMuTau_WplusJets_part04.root'
)

processZtoMuTau_WplusJets.config_dqmFileLoader.scaleFactor = cms.double(corrFactorZtoMuTau_WplusJets*intLumiZtoMuTau_Data/intLumiZtoMuTau_WplusJets)

processZtoMuTau_WplusJetsSum = copy.deepcopy(process_WplusJets)
processZtoMuTau_WplusJetsSum.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_WplusJetsSum.root'
)
processZtoMuTau_WplusJetsSum.config_dqmFileLoader.dqmDirectory_store = cms.string('harvested')
processZtoMuTau_WplusJetsSum.config_dqmFileLoader.scaleFactor = cms.double(1.)

#--------------------------------------------------------------------------------

processZtoMuTau_InclusivePPmuX = copy.deepcopy(process_InclusivePPmuX)
processZtoMuTau_InclusivePPmuX.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_InclusivePPmuX.root'
)
processZtoMuTau_InclusivePPmuX.config_dqmFileLoader.dqmDirectory_store = cms.string('harvested/InclusivePPmuX')
processZtoMuTau_InclusivePPmuX.config_dqmFileLoader.scaleFactor = cms.double(corrFactorZtoMuTau_InclusivePPmuX*intLumiZtoMuTau_Data/intLumiZtoMuTau_InclusivePPmuX)

#--------------------------------------------------------------------------------

processZtoMuTau_PPmuXptGt20 = copy.deepcopy(process_PPmuXptGt20)
processZtoMuTau_PPmuXptGt20.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_PPmuXptGt20_part01.root',
    'plotsZtoMuTau_PPmuXptGt20_part02.root',
    'plotsZtoMuTau_PPmuXptGt20_part03.root',
    'plotsZtoMuTau_PPmuXptGt20_part04.root',
    'plotsZtoMuTau_PPmuXptGt20_part05.root'
)
processZtoMuTau_PPmuXptGt20.config_dqmFileLoader.scaleFactor = cms.double(corrFactorZtoMuTau_PPmuXptGt20*intLumiZtoMuTau_Data/intLumiZtoMuTau_PPmuXptGt20)

processZtoMuTau_PPmuXptGt20Sum = copy.deepcopy(process_PPmuXptGt20)
processZtoMuTau_PPmuXptGt20Sum.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_PPmuXptGt20Sum.root'
)
processZtoMuTau_PPmuXptGt20Sum.config_dqmFileLoader.dqmDirectory_store = cms.string('harvested')
processZtoMuTau_PPmuXptGt20Sum.config_dqmFileLoader.scaleFactor = cms.double(1.)

#--------------------------------------------------------------------------------

processZtoMuTau_TTplusJets = copy.deepcopy(process_TTplusJets)
processZtoMuTau_TTplusJets.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_TTplusJets_part01.root',
    'plotsZtoMuTau_TTplusJets_part02.root',
    'plotsZtoMuTau_TTplusJets_part03.root',
    'plotsZtoMuTau_TTplusJets_part04.root',
    'plotsZtoMuTau_TTplusJets_part05.root'
)
processZtoMuTau_TTplusJets.config_dqmFileLoader.scaleFactor = cms.double(corrFactorZtoMuTau_TTplusJets*intLumiZtoMuTau_Data/intLumiZtoMuTau_TTplusJets)

processZtoMuTau_TTplusJetsSum = copy.deepcopy(process_TTplusJets)
processZtoMuTau_TTplusJetsSum.config_dqmFileLoader.inputFileNames = cms.vstring(
    'plotsZtoMuTau_TTplusJetsSum.root'
)
processZtoMuTau_TTplusJetsSum.config_dqmFileLoader.dqmDirectory_store = cms.string('harvested')
processZtoMuTau_TTplusJetsSum.config_dqmFileLoader.scaleFactor = cms.double(1.)
