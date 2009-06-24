#!/usr/bin/env python

from TauAnalysis.Configuration.submitToBatch import submitToBatch
from TauAnalysis.BgEstimationTools.makeReplacementsBgEstGenMatrixControlPlots import makeReplacementsBgEstGenMatrixControlPlots

# name of the directory (either on afs area or castor)
# to which all .root files produced by the cmsRun job will be copied
outputDirectory = "/castor/cern.ch/user/v/veelken/bgEstPlots/ZtoMuTau/"

# small cmsRun job for testing purposes...
#submitToBatch(configFile = "prodGenMatrixControlPlotsZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "Ztautau",
#              replFunction = makeReplacementsBgEstGenMatrixControlPlots, replacements = "maxEvents = 100",
#              job = "bgEstGenMatrixControlPlots", queue = "1nh", outputDirectory = outputDirectory)

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
submitToBatch(configFile = "prodGenMatrixControlPlotsZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "Ztautau",
              replFunction = makeReplacementsBgEstGenMatrixControlPlots, replacements = "maxEvents = -1",
              job = "bgEstGenMatrixControlPlots", queue = "1nd", outputDirectory = outputDirectory)

# Z --> mu mu jobs
submitToBatch(configFile = "prodGenMatrixControlPlotsZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "Zmumu",
              replFunction = makeReplacementsBgEstGenMatrixControlPlots, replacements = "maxEvents = -1",
              job = "bgEstGenMatrixControlPlots", queue = "1nd", outputDirectory = outputDirectory)

# pp --> mu X QCD jobs
submitToBatch(configFile = "prodGenMatrixControlPlotsZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "InclusivePPmuX",
              replFunction = makeReplacementsBgEstGenMatrixControlPlots, replacements = "maxEvents = -1",
              job = "bgEstGenMatrixControlPlots", queue = "1nd", outputDirectory = outputDirectory)

for i in range(3):
    submitToBatch(configFile = "prodGenMatrixControlPlotsZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "PPmuXptGt20_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsBgEstGenMatrixControlPlots, replacements = "maxEvents = -1",
                  job = "bgEstGenMatrixControlPlots", queue = "1nd", outputDirectory = outputDirectory)

# W/Z + jets jobs
submitToBatch(configFile = "prodGenMatrixControlPlotsZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "WplusJets",
              replFunction = makeReplacementsBgEstGenMatrixControlPlots, replacements = "maxEvents = -1",
              job = "bgEstGenMatrixControlPlots", queue = "1nd", outputDirectory = outputDirectory)

submitToBatch(configFile = "prodGenMatrixControlPlotsZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "ZeePlusJets",
              replFunction = makeReplacementsBgEstGenMatrixControlPlots, replacements = "maxEvents = -1",
              job = "bgEstGenMatrixControlPlots", queue = "1nd", outputDirectory = outputDirectory)
submitToBatch(configFile = "prodGenMatrixControlPlotsZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "ZmumuPlusJets",
              replFunction = makeReplacementsBgEstGenMatrixControlPlots, replacements = "maxEvents = -1",
              job = "bgEstGenMatrixControlPlots", queue = "1nd", outputDirectory = outputDirectory)
submitToBatch(configFile = "prodGenMatrixControlPlotsZtoMuTau_cfg.py", channel = "ZtoMuTau", sample = "ZtautauPlusJets",
              replFunction = makeReplacementsBgEstGenMatrixControlPlots, replacements = "maxEvents = -1",
              job = "bgEstGenMatrixControlPlots", queue = "1nd", outputDirectory = outputDirectory)

