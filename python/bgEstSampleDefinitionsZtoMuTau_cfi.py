import FWCore.ParameterSet.Config as cms
import copy

# define configuration parameters for submission of Z --> mu + tau-jet jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)

bgEstSampleDirectoryName = cms.string("rfio:/castor/cern.ch/user/v/veelken/bgEstSkim/ZtoMuTau/")
bgEstNtupleDirectoryName = cms.string("rfio:/castor/cern.ch/user/v/veelken/bgEstNtuples/ZtoMuTau/")

# definitions for producing skimmed samples
bgEstSampleOutputFileNameZtautau = cms.string('bgEstSampleZtoMuTau_Ztautau_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesZtautau = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Ztautau_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Ztautau_part02.root'
)

# definitions for producing ntuples
bgEstNtupleOutputFileNameZtautau = cms.string('bgEstNtupleZtoMuTau_Ztautau.root')

#--------------------------------------------------------------------------------

# definitions for producing skimmed samples
bgEstSampleOutputFileNameZmumu = cms.string('bgEstSampleZtoMuTau_Zmumu_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesZmumu = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Zmumu_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Zmumu_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_Zmumu_part03.root'
)

# definitions for producing ntuples
bgEstNtupleOutputFileNameZmumu = cms.string('bgEstNtupleZtoMuTau_Zmumu.root')

#--------------------------------------------------------------------------------

# definitions for producing skimmed samples
bgEstSampleOutputFileNameZtautauPlusJets = cms.string('bgEstSampleZtoMuTau_ZtautauPlusJets_partXX.root')

bgEstSampleOutputFileNameZmumuPlusJets = cms.string('bgEstSampleZtoMuTau_ZmumuPlusJets_partXX.root')

bgEstSampleOutputFileNameZeePlusJets = cms.string('bgEstSampleZtoMuTau_ZeePlusJets_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesZtautauPlusJets = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZtautauPlusJets_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZtautauPlusJets_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZtautauPlusJets_part03.root'
)

bgEstSampleFileNamesZmumuPlusJets = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZmumuPlusJets_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZmumuPlusJets_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZmumuPlusJets_part03.root'
)

bgEstSampleFileNamesZeePlusJets = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZeePlusJets_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZeePlusJets_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_ZeePlusJets_part03.root'
)

# definitions for producing ntuples
bgEstNtupleOutputFileNameZtautauPlusJets = cms.string('bgEstNtupleZtoMuTau_ZtautauPlusJets.root')

bgEstNtupleOutputFileNameZmumuPlusJets = cms.string('bgEstNtupleZtoMuTau_ZmumuPlusJets.root')

bgEstNtupleOutputFileNameZeePlusJets = cms.string('bgEstNtupleZtoMuTau_ZeePlusJets.root')

#--------------------------------------------------------------------------------

# definitions for producing skimmed samples
bgEstSampleOutputFileNameWplusJets = cms.string('bgEstSampleZtoMuTau_WplusJets_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesWplusJets = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part03.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part04.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part05.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_WplusJets_part06.root'
)

# definitions for producing ntuples
bgEstNtupleOutputFileNameWplusJets = cms.string('bgEstNtupleZtoMuTau_WplusJets.root')

#--------------------------------------------------------------------------------

# definitions for producing skimmed samples
bgEstSampleOutputFileNameInclusivePPmuX = cms.string('bgEstSampleZtoMuTau_InclusivePPmuX.root')

bgEstSampleOutputFileNamePPmuXptGt20 = cms.string('bgEstSampleZtoMuTau_PPmuXptGt20_partXX.root')

# definitions for reading skimmed samples
bgEstSampleFileNamesPPmuXptGt20_part01  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part01.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part02.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part03.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part04.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part05.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part06.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part07.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part08.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part09.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part10.root'
)

bgEstSampleFileNamesPPmuXptGt20_part02  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part11.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part12.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part13.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part14.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part16.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part18.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part19.root'
)

bgEstSampleFileNamesPPmuXptGt20_part03  = cms.untracked.vstring(
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part21.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part27.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part28.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part30.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part31.root',
    bgEstSampleDirectoryName.value() + 'bgEstSampleZtoMuTau_PPmuXptGt20_part35.root'
)

# definitions for producing ntuples
bgEstNtupleOutputFileNameInclusivePPmuX = cms.string('bgEstNtupleZtoMuTau_InclusivePPmuX.root')

bgEstNtupleOutputFileNamePPmuXptGt20 = cms.string('bgEstNtupleZtoMuTau_PPmuXptGt20_partXX.root')

