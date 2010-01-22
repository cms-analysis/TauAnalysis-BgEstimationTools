import FWCore.ParameterSet.Config as cms
import copy

# define configuration parameters for submission of Z --> mu + tau-jet jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)

bgEstSampleDirectoryName = cms.string("rfio:/castor/cern.ch/user/v/veelken/bgEstSkim/ZtoMuTauII/")

# definitions for producing skimmed samples
bgEstSampleOutputFileNameZtoMuTau_Ztautau = cms.untracked.string('bgEstSampleZtoMuTau_Ztautau_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesZtoMuTau_Ztautau = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Ztautau_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Ztautau_part02.root'
)

# definitions for producing ntuples
bgEstNtupleOutputFileNameZtoMuTau_Ztautau = cms.string('bgEstNtupleZtoMuTau_Ztautau.root')

#--------------------------------------------------------------------------------

# definitions for producing skimmed samples
bgEstSampleOutputFileNameZtoMuTau_Zmumu = cms.untracked.string('bgEstSampleZtoMuTau_Zmumu_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesZtoMuTau_Zmumu_part01 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Zmumu_part01.root'
)
bgEstSampleFileNamesZtoMuTau_Zmumu_part01a = bgEstSampleFileNamesZtoMuTau_Zmumu_part01
bgEstSampleFileNamesZtoMuTau_Zmumu_part01b = bgEstSampleFileNamesZtoMuTau_Zmumu_part01

bgEstSampleFileNamesZtoMuTau_Zmumu_part02 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Zmumu_part02.root'
)
bgEstSampleFileNamesZtoMuTau_Zmumu_part02a = bgEstSampleFileNamesZtoMuTau_Zmumu_part02
bgEstSampleFileNamesZtoMuTau_Zmumu_part02b = bgEstSampleFileNamesZtoMuTau_Zmumu_part02

bgEstSampleFileNamesZtoMuTau_Zmumu_part03 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Zmumu_part03.root'
)
bgEstSampleFileNamesZtoMuTau_Zmumu_part03a = bgEstSampleFileNamesZtoMuTau_Zmumu_part03
bgEstSampleFileNamesZtoMuTau_Zmumu_part03b = bgEstSampleFileNamesZtoMuTau_Zmumu_part03

bgEstSampleFileNamesZtoMuTau_Zmumu_part04 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Zmumu_part04.root'
)
bgEstSampleFileNamesZtoMuTau_Zmumu_part04a = bgEstSampleFileNamesZtoMuTau_Zmumu_part04
bgEstSampleFileNamesZtoMuTau_Zmumu_part04b = bgEstSampleFileNamesZtoMuTau_Zmumu_part04

bgEstSampleFileNamesZtoMuTau_Zmumu_part05 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Zmumu_part05.root'
)
bgEstSampleFileNamesZtoMuTau_Zmumu_part05a = bgEstSampleFileNamesZtoMuTau_Zmumu_part05
bgEstSampleFileNamesZtoMuTau_Zmumu_part05b = bgEstSampleFileNamesZtoMuTau_Zmumu_part05

# definitions for producing ntuples
bgEstNtupleOutputFileNameZtoMuTau_Zmumu = cms.string('bgEstNtupleZtoMuTau_Zmumu_partXX.root')

#--------------------------------------------------------------------------------

# definitions for producing skimmed samples
bgEstSampleOutputFileNameZtoMuTau_ZtautauPlusJets = cms.untracked.string('bgEstSampleZtoMuTau_ZtautauPlusJets_partXX.root')

bgEstSampleOutputFileNameZtoMuTau_ZmumuPlusJets = cms.untracked.string('bgEstSampleZtoMuTau_ZmumuPlusJets_partXX.root')

bgEstSampleOutputFileNameZtoMuTau_ZeePlusJets = cms.untracked.string('bgEstSampleZtoMuTau_ZeePlusJets_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesZtoMuTau_ZtautauPlusJets = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZtautauPlusJets_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZtautauPlusJets_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZtautauPlusJets_part03.root'
)

bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part01 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZmumuPlusJets_part01.root'
)
bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part01a = bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part01
bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part01b = bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part01

bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part02 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZmumuPlusJets_part02.root'
)
bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part02a = bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part02
bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part02b = bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part02

bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part03 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZmumuPlusJets_part03.root'
)
bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part03a = bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part03
bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part03b = bgEstSampleFileNamesZtoMuTau_ZmumuPlusJets_part03

bgEstSampleFileNamesZtoMuTau_ZeePlusJets = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZeePlusJets_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZeePlusJets_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZeePlusJets_part03.root'
)

# definitions for producing ntuples
bgEstNtupleOutputFileNameZtoMuTau_ZtautauPlusJets = cms.string('bgEstNtupleZtoMuTau_ZtautauPlusJets.root')

bgEstNtupleOutputFileNameZtoMuTau_ZmumuPlusJets = cms.string('bgEstNtupleZtoMuTau_ZmumuPlusJets_partXX.root')

bgEstNtupleOutputFileNameZtoMuTau_ZeePlusJets = cms.string('bgEstNtupleZtoMuTau_ZeePlusJets.root')

#--------------------------------------------------------------------------------

# definitions for producing skimmed samples
bgEstSampleOutputFileNameZtoMuTau_WplusJets = cms.untracked.string('bgEstSampleZtoMuTau_WplusJets_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesZtoMuTau_WplusJets_part01 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part03.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part04.root'
)

bgEstSampleFileNamesZtoMuTau_WplusJets_part02 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part05.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part06.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part07.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part08.root'
)

bgEstSampleFileNamesZtoMuTau_WplusJets_part03 = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part09.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part10.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part11.root'
)

# definitions for producing ntuples
bgEstNtupleOutputFileNameZtoMuTau_WplusJets = cms.string('bgEstNtupleZtoMuTau_WplusJets_partXX.root')

#--------------------------------------------------------------------------------

# definitions for producing skimmed samples
bgEstSampleOutputFileNameZtoMuTau_InclusivePPmuX = cms.untracked.string('bgEstSampleZtoMuTau_InclusivePPmuX.root')

bgEstSampleOutputFileNameZtoMuTau_PPmuXptGt20 = cms.untracked.string('bgEstSampleZtoMuTau_PPmuXptGt20_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part01  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part03.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part02  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part04.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part05.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part06.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part03  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part07.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part08.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part09.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part04  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part10.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part11.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part12.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part05  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part13.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part14.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part15.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part06  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part16.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part17.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part18.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part07  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part19.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part20.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part21.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part08  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part22.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part23.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part24.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part09  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part25.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part26.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part27.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part10  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part28.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part29.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part30.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part11  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part31.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part32.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part33.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part12  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part34.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part35.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part36.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part13  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part37.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part38.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part39.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part14  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part40.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part41.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part42.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part15  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part43.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part44.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part45.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part16  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part46.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part47.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part48.root'
)

bgEstSampleFileNamesZtoMuTau_PPmuXptGt20_part17  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part49.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part50.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part51.root'
)

# definitions for producing ntuples
bgEstNtupleOutputFileNameZtoMuTau_InclusivePPmuX = cms.string('bgEstNtupleZtoMuTau_InclusivePPmuX.root')

bgEstNtupleOutputFileNameZtoMuTau_PPmuXptGt20 = cms.string('bgEstNtupleZtoMuTau_PPmuXptGt20_partXX.root')

#--------------------------------------------------------------------------------

# definitions for producing skimmed samples
bgEstSampleOutputFileNameZtoMuTau_TTplusJets = cms.untracked.string('bgEstSampleZtoMuTau_TTplusJets_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesZtoMuTau_TTplusJets_part01  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part03.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part04.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part05.root'
)

bgEstSampleFileNamesZtoMuTau_TTplusJets_part02  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part06.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part07.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part08.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part09.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part10.root'
)

bgEstSampleFileNamesZtoMuTau_TTplusJets_part03  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part11.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part12.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part13.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part14.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part15.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_TTplusJets_part16.root'
)

# definitions for producing ntuples
bgEstNtupleOutputFileNameZtoMuTau_TTplusJets = cms.string('bgEstNtupleZtoMuTau_TTplusJets_partXX.root')

