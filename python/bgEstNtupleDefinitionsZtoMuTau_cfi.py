import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
# define Ntuple files used for data-driven background estimation methods
# in Z --> mu + tau-jet channel
#--------------------------------------------------------------------------------

bgEstNtupleDirectoryName = cms.string("rfio:/castor/cern.ch/user/v/veelken/bgEstNtuples/ZtoMuTauII/")

# Z --> tau+ tau+ sample
fileNamesZtoMuTau_Ztautau = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Ztautau.root'
)

# Z --> mu+ mu- sample
fileNamesZtoMuTau_Zmumu = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part01a.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part01b.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part02a.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part02b.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part03a.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part03b.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part04a.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part04b.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part05a.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part05b.root'
)

#--------------------------------------------------------------------------------

# Z --> tau+ tau- + jets sample
fileNamesZtoMuTau_ZtautauPlusJets = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZtautauPlusJets.root'
)

# Z --> mu+ mu- + jets sample
fileNamesZtoMuTau_ZmumuPlusJets = cms.vstring(
    #bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZmumuPlusJets_part01a.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZmumuPlusJets_part01b.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZmumuPlusJets_part02a.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZmumuPlusJets_part02b.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZmumuPlusJets_part03a.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZmumuPlusJets_part03b.root'    
)

# Z --> e+ e- + jets sample
fileNamesZtoMuTau_ZeePlusJets = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZeePlusJets.root'
)

#--------------------------------------------------------------------------------

# W + jets sample
fileNamesZtoMuTau_WplusJets = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_WplusJets_part01.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_WplusJets_part02.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_WplusJets_part03.root'
)

#--------------------------------------------------------------------------------

# pp --> mu X QCD sample
fileNamesZtoMuTau_qcdSum = cms.vstring(
    #bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_InclusivePPmuX.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part01.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part02.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part03.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part04.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part05.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part06.root',
    #bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part07.root',
    #bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part08.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part09.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part10.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part11.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part12.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part13.root'
    #bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part14.root',
    #bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part15.root',
    #bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part16.root',
    #bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part17.root'
)

#--------------------------------------------------------------------------------

# ttbar + jets sample
fileNamesZtoMuTau_TTplusJets = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_TTplusJets_part01.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_TTplusJets_part02.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_TTplusJets_part03.root'
)

#--------------------------------------------------------------------------------

# (pseudo)Data sample = sum of all signal + background Monte Carlo samples
fileNamesZtoMuTau_pseudoData = cms.vstring()
fileNamesZtoMuTau_pseudoData.extend(fileNamesZtoMuTau_Ztautau)
fileNamesZtoMuTau_pseudoData.extend(fileNamesZtoMuTau_ZmumuPlusJets)
fileNamesZtoMuTau_pseudoData.extend(fileNamesZtoMuTau_WplusJets)
fileNamesZtoMuTau_pseudoData.extend(fileNamesZtoMuTau_TTplusJets)
fileNamesZtoMuTau_pseudoData.extend(fileNamesZtoMuTau_qcdSum)
