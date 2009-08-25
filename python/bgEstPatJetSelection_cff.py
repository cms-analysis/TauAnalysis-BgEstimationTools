import FWCore.ParameterSet.Config as cms

from TauAnalysis.CandidateTools.tools.objSelConfigurator import *

#--------------------------------------------------------------------------------  
# produce collections of pat::Jets passing selection criteria
# specific to data-driven background estimation methods
#--------------------------------------------------------------------------------

# see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
# on how to use the cut-string parser

# select jets with Et > 40 GeV
bgEstSelectedLayer1JetsEt40 = cms.EDFilter("PATJetSelector",
    cut = cms.string('et > 40.'),
    filter = cms.bool(False)
)

# select jets with Et > 60 GeV
bgEstSelectedLayer1JetsEt60 = cms.EDFilter("PATJetSelector",
    cut = cms.string('et > 60.'),
    filter = cms.bool(False)
)

bgEstPatJetSelConfiguratorEt = objSelConfigurator(
    [ bgEstSelectedLayer1JetsEt40,
      bgEstSelectedLayer1JetsEt60 ],
    src = "selectedLayer1JetsEt20Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

bgEstSelectLayer1JetsEt = bgEstPatJetSelConfiguratorEt.configure(namespace = locals())

# select jets with quantity alpha > 0.1
# (defined as ratio of sum of charged particle transverse momenta 
#  to sum of charged plus neutral particle transverse momenta)
bgEstSelectedLayer1JetsAlpha0_1 = cms.EDFilter("PATJetAlphaSelector",
    alphaMin = cms.double(0.1),
    filter = cms.bool(False)
)

# select jets with quantity alpha > 0.3
bgEstSelectedLayer1JetsAlpha0_3 = cms.EDFilter("PATJetAlphaSelector",
    alphaMin = cms.double(0.3),
    filter = cms.bool(False)
)

bgEstPatJetSelConfiguratorAlpha = objSelConfigurator(
    [ bgEstSelectedLayer1JetsAlpha0_1,
      bgEstSelectedLayer1JetsAlpha0_3 ],
    src = "selectedLayer1JetsEt20Cumulative",
    pyModuleName = __name__,
    doSelIndividual = False
)

bgEstSelectLayer1JetsAlpha = bgEstPatJetSelConfiguratorAlpha.configure(namespace = locals())

bgEstSelectLayer1Jets = cms.Sequence(bgEstSelectLayer1JetsEt * bgEstSelectLayer1JetsAlpha)
