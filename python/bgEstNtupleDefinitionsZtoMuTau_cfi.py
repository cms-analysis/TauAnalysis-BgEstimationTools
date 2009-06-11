import FWCore.ParameterSet.Config as cms
import copy

# define Ntuple files used for data-driven background estimation methods
# in Z --> mu + tau-jet channel

ntupleOutputDirectoryName = cms.string('./')

# Z --> tau tau sample
fileNames_Ztautau = cms.vstring(
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_Ztautau_part1.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_Ztautau_part2.root'
)

# Z --> mu mu sample
fileNames_Zmumu = cms.vstring(
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part1.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part2.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu_part3.root'
)

#--------------------------------------------------------------------------------

# Z --> tau tau + jets sample
fileNames_ZtautauPlusJets = cms.vstring(
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZtautauPlusJets_part1.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZtautauPlusJets_part2.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZtautauPlusJets_part3.root'
)

# Z --> mu mu + jets sample
fileNames_ZmumuPlusJets = cms.vstring(
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZmumuPlusJets_part1.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZmumuPlusJets_part2.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZmumuPlusJets_part3.root'
)

# Z --> e e + jets sample
fileNames_ZeePlusJets = cms.vstring(
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZeePlusJets_part1.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZeePlusJets_part2.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZeePlusJets_part3.root'
)

#--------------------------------------------------------------------------------

# W + jets sample
fileNames_WplusJets = cms.vstring(
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_WplusJets_part1.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_WplusJets_part2.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_WplusJets_part3.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_WplusJets_part4.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_WplusJets_part5.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_WplusJets_part6.root'
)

#--------------------------------------------------------------------------------

# pp --> mu X QCD sample
fileNames_qcdSum = cms.vstring(
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_InclusivePPmuX.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part1.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part2.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part3.root',
    #ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part4.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part5.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part6.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part7.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part8.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part9.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part10.root',
    #ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part11.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part12.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part13.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part14.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part16.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part18.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part19.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part21.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part27.root',
    #ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part28.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part30.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part31.root',
    ntupleOutputDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part35.root'
)

#--------------------------------------------------------------------------------

# (pseudo)Data sample = sum of all signal + background Monte Carlo samples
fileNames_pseudoData = cms.vstring()
fileNames_pseudoData.extend(fileNames_Ztautau)
#fileNames_pseudoData.extend(fileNames_Zmumu)
fileNames_pseudoData.extend(fileNames_WplusJets)
#fileNames_pseudoData.extend(fileNames_qcdSum)
