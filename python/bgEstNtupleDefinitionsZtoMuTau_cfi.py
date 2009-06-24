import FWCore.ParameterSet.Config as cms
import copy

# define Ntuple files used for data-driven background estimation methods
# in Z --> mu + tau-jet channel

from TauAnalysis.BgEstimationTools.bgEstSampleDefinitionsZtoMuTau_cfi import *

# Z --> tau tau sample
fileNames_Ztautau = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Ztautau.root'
)

# Z --> mu mu sample
fileNames_Zmumu = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_Zmumu.root'
)

#--------------------------------------------------------------------------------

# Z --> tau tau + jets sample
fileNames_ZtautauPlusJets = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZtautauPlusJets.root'
)

# Z --> mu mu + jets sample
fileNames_ZmumuPlusJets = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZmumuPlusJets.root'
)

# Z --> e e + jets sample
fileNames_ZeePlusJets = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_ZeePlusJets.root'
)

#--------------------------------------------------------------------------------

# W + jets sample
fileNames_WplusJets = cms.vstring(
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_WplusJets.root'
)

#--------------------------------------------------------------------------------

# pp --> mu X QCD sample
fileNames_qcdSum = cms.vstring(
    #bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_InclusivePPmuX.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part01.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part02.root',
    bgEstNtupleDirectoryName.value() + 'bgEstNtupleZtoMuTau_PPmuXptGt20_part03.root'
)

#--------------------------------------------------------------------------------

# (pseudo)Data sample = sum of all signal + background Monte Carlo samples
fileNames_pseudoData = cms.vstring()
fileNames_pseudoData.extend(fileNames_Ztautau)
#fileNames_pseudoData.extend(fileNames_Zmumu)
fileNames_pseudoData.extend(fileNames_WplusJets)
#fileNames_pseudoData.extend(fileNames_qcdSum)
