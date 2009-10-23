#!/usr/bin/env python

from TauAnalysis.Configuration.submitToBatch import submitToBatch
from TauAnalysis.Configuration.makeReplacementsHarvesting import makeReplacementsHarvesting

# name of the directory (either on afs area or castor)
# to which all .root files produced by the cmsRun job will be copied
outputFilePath = "/castor/cern.ch/user/v/veelken/bgEstPlots/ZtoMuTau/"

inputFilePath = "rfio:" + outputFilePath

#--------------------------------------------------------------------------------
#
# Add histograms, numbers in FilterStatisticsTables and run + event numbers
# stored as DQM MonitorElements in different ROOT files
#
# NOTE: The jobs get submitted to the '1nh' queue,
#       which allows for an execution time of the cmsRun jobs of up to 1 hour
#       (the queues are {'1nh' (1 hour), '1nd' (24 hours) and '1nw' (1 week execution time limit);
#        see https://twiki.cern.ch/twiki/bin/view/CMS/CMSUKCMSSWBatch for details about the CERN batch system)           
#
#--------------------------------------------------------------------------------
# harvest Z --> tau tau 
submitToBatch(configFile = "../../Configuration/test/harvestZtoMuTauPlots_cfg.py", channel = "ZtoMuTau",
              sample = "Ztautau",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "frHarvesting", queue = "1nh", outputFilePath = outputFilePath)
submitToBatch(configFile = "../../Configuration/test/harvestZtoMuTauPlots_cfg.py", channel = "ZtoMuTau",
              sample = "Ztautau_from_selZmumu",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "frHarvesting", queue = "1nh", outputFilePath = outputFilePath)

# harvest Z --> mu mu
submitToBatch(configFile = "../../Configuration/test/harvestZtoMuTauPlots_cfg.py", channel = "ZtoMuTau",
              sample = "Zmumu",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "frHarvesting", queue = "1nh", outputFilePath = outputFilePath)

# harvest PPmuXptGt20
for i in range(2):
    submitToBatch(configFile = "../../Configuration/test/harvestZtoMuTauPlots_cfg.py", channel = "ZtoMuTau",
                  sample = "PPmuXptGt20_part%(i)02d" % {"i" : (i + 1)},
                  replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
                  job = "frHarvesting", queue = "1nh", outputFilePath = outputFilePath)

# harvest W/Z + jets
submitToBatch(configFile = "../../Configuration/test/harvestZtoMuTauPlots_cfg.py", channel = "ZtoMuTau",
              sample = "WplusJets",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "frHarvesting", queue = "1nh", outputFilePath = outputFilePath)

submitToBatch(configFile = "../../Configuration/test/harvestZtoMuTauPlots_cfg.py", channel = "ZtoMuTau",
              sample = "ZtautauPlusJets",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "frHarvesting", queue = "1nh", outputFilePath = outputFilePath)
submitToBatch(configFile = "../../Configuration/test/harvestZtoMuTauPlots_cfg.py", channel = "ZtoMuTau",
              sample = "ZmumuPlusJets",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "frHarvesting", queue = "1nh", outputFilePath = outputFilePath)
submitToBatch(configFile = "../../Configuration/test/harvestZtoMuTauPlots_cfg.py", channel = "ZtoMuTau",
              sample = "ZeePlusJets",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "frHarvesting", queue = "1nh", outputFilePath = outputFilePath)

# harvest TTplusJets
submitToBatch(configFile = "../../Configuration/test/harvestZtoMuTauPlots_cfg.py", channel = "ZtoMuTau",
              sample = "TTplusJets",
              replFunction = makeReplacementsHarvesting, replacements = "inputFilePath = " + inputFilePath,
              job = "frHarvesting", queue = "1nh", outputFilePath = outputFilePath)
