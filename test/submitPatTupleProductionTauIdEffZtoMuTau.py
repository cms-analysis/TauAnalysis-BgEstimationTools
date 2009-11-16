#!/usr/bin/env python

from TauAnalysis.Configuration.submitToBatch import submitToBatch
from TauAnalysis.Configuration.makeReplacementsPatTupleProduction import makeReplacementsPatTupleProduction

# name of the directory (either on afs area or castor)
# to which all .root files produced by the cmsRun job will be copied
outputFilePath = "/castor/cern.ch/user/v/veelken/patTuples/tauIdEff/"

# small cmsRun job for testing purposes... 
#submitToBatch(configFile = "prodPatTupleTauIdEff_cfg.py", channel = "ZtoMuTau", sample = "Ztautau_part01",
#              replFunction = makeReplacementsPatTupleProduction, replacements = "maxEvents = 100",
#              job = "patTupleProduction", queue = "1nh", outputDirectory = outputDirectory)

#--------------------------------------------------------------------------------
#
# Monte Carlo samples from Summer'08 production
# reprocessed with CMSSW_2_2_3, first skimmed by Letizia and Monica,
# skimmed a second time by Christian
#
# NOTE: The jobs get submitted to the '1nd' queue,
#       which allows for an execution time of the cmsRun jobs of up to 24 hours
#       (the queues are {'1nh' (1 hour), '1nd' (24 hours) and '1nw' (1 week execution time limit);
#        see https://twiki.cern.ch/twiki/bin/view/CMS/CMSUKCMSSWBatch for details about the CERN batch system)           
#
#--------------------------------------------------------------------------------

# Z --> tau tau jobs
for i in range(2):
    submitToBatch(configFile = "prodPatTupleTauIdEff_cfg.py", channel = "ZtoMuTau",
                  sample = "Ztautau_part%(i)02d" % {"i" : (i + 1)},
		  replFunction = makeReplacementsPatTupleProduction, replacements = "maxEvents = -1",
                  job = "patTupleProduction", queue = "1nd", outputFilePath = outputFilePath)

# Z --> mu mu jobs
for i in range(3):
    submitToBatch(configFile = "prodPatTupleTauIdEff_cfg.py", channel = "ZtoMuTau",
                  sample = "Zmumu_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsPatTupleProduction, replacements = "maxEvents = -1",
                  job = "patTupleProduction", queue = "1nw", outputFilePath = outputFilePath)

# pp --> mu X QCD jobs
submitToBatch(configFile = "prodPatTupleTauIdEff_cfg.py", channel = "ZtoMuTau",
              sample = "InclusivePPmuX",
	      replFunction = makeReplacementsPatTupleProduction, replacements = "maxEvents = -1",
              job = "patTupleProduction", queue = "1nd", outputFilePath = outputFilePath)

for i in range(5):
    submitToBatch(configFile = "prodPatTupleTauIdEff_cfg.py", channel = "ZtoMuTau",
                  sample = "PPmuXptGt20_part%(i)02d" % {"i" : (i + 1)},
		  replFunction = makeReplacementsPatTupleProduction, replacements = "maxEvents = -1",
                  job = "patTupleProduction", queue = "1nd", outputFilePath = outputFilePath)

# W/Z + jets jobs
for i in range(4):
    submitToBatch(configFile = "prodPatTupleTauIdEff_cfg.py", channel = "ZtoMuTau",
                  sample = "WplusJets_part%(i)02d" % {"i" : (i + 1)},
		  replFunction = makeReplacementsPatTupleProduction, replacements = "maxEvents = -1",
                  job = "patTupleProduction", queue = "1nd", outputFilePath = outputFilePath)

submitToBatch(configFile = "prodPatTupleTauIdEff_cfg.py", channel = "ZtoMuTau",
              sample = "ZeePlusJets",
	      replFunction = makeReplacementsPatTupleProduction, replacements = "maxEvents = -1",
              job = "patTupleProduction", queue = "1nd", outputFilePath = outputFilePath)
submitToBatch(configFile = "prodPatTupleTauIdEff_cfg.py", channel = "ZtoMuTau",
              sample = "ZmumuPlusJets",
	      replFunction = makeReplacementsPatTupleProduction, replacements = "maxEvents = -1",
              job = "patTupleProduction", queue = "1nw", outputFilePath = outputFilePath)
submitToBatch(configFile = "prodPatTupleTauIdEff_cfg.py", channel = "ZtoMuTau",
              sample = "ZtautauPlusJets",
	      replFunction = makeReplacementsPatTupleProduction, replacements = "maxEvents = -1",
              job = "patTupleProduction", queue = "1nd", outputFilePath = outputFilePath)

# ttbar + jets  jobs
for i in range(5):
    submitToBatch(configFile = "prodPatTupleTauIdEff_cfg.py", channel = "ZtoMuTau",
                  sample = "TTplusJets_part%(i)02d" % {"i" : (i + 1)},
		  replFunction = makeReplacementsPatTupleProduction, replacements = "maxEvents = -1",
                  job = "patTupleProduction", queue = "1nd", outputFilePath = outputFilePath)


