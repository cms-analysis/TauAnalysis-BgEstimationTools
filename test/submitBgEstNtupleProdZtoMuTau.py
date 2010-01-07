#!/usr/bin/env python

from TauAnalysis.Configuration.submitToBatch import submitToBatch
from TauAnalysis.BgEstimationTools.makeReplacementsBgEstNtupleProd import makeReplacementsBgEstNtupleProd

# name of the directory (either on afs area or castor)
# to which all .root files produced by the cmsRun job will be copied
outputFilePath = "/castor/cern.ch/user/v/veelken/bgEstNtuples/ZtoMuTau/"

# small cmsRun job for testing purposes...
#submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "Ztautau",
#              replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = 100; skipEvents = 0",
#              job = "bgEstNtupleProd", queue = "1nh", outputFilePath = outputFilePath)

#--------------------------------------------------------------------------------
#
# Monte Carlo samples from Summer'08 production
# reprocessed with CMSSW_2_2_3, skimmed by Letizia and Monica
#
# NOTE: The jobs get submitted to the '1nd' queue,
#       which allows for an execution time of the cmsRun jobs of up to 24 hours
#       (the queues are {'1nh' (1 hour), '1nd' (24 hours) and '1nw' (1 week execution time limit);
#        see https://twiki.cern.ch/twiki/bin/view/CMS/CMSUKCMSSWBatch for details about the CERN batch system)           
#
#--------------------------------------------------------------------------------

# Z --> tau tau jobs
submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "Ztautau",
              replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = -1; skipEvents = 0",
              job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)

# Z --> mu mu jobs
for i in range(5):
    submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "Zmumu_part%(i)02da" % {"i" : (i + 1)},
                  replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = 25000; skipEvents = 0",
                  job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)
    submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "Zmumu_part%(i)02db" % {"i" : (i + 1)},
                  replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = -1; skipEvents = 25000",
                  job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)

# pp --> mu X QCD jobs
submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "InclusivePPmuX",
              replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = -1; skipEvents = 0",
              job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)

for i in range(5):
    submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "PPmuXptGt20_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = -1; skipEvents = 0",
                  job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)

# W/Z + jets jobs
for i in range(3):
    submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "WplusJets_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = -1; skipEvents = 0",
                  job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)

submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "ZeePlusJets",
              replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = -1; skipEvents = 0",
              job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)
for i in range(3):
    submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "ZmumuPlusJets_part%(i)02da" % {"i" : (i + 1)},
                  replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = 25000; skipEvents = 0",
                  job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)
    submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "ZmumuPlusJets_part%(i)02db" % {"i" : (i + 1)},
                  replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = -1; skipEvents = 25000",
                  job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)
submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "ZtautauPlusJets",
              replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = -1; skipEvents = 0",
              job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)

# ttbar + jets  jobs
for i in range(3):
    submitToBatch(configFile = "prodNtupleZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "TTplusJets_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsBgEstNtupleProd, replacements = "maxEvents = -1; skipEvents = 0",
                  job = "bgEstNtupleProd", queue = "1nd", outputFilePath = outputFilePath)
