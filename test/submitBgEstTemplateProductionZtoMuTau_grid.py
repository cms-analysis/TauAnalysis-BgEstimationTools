#!/usr/bin/env python

from TauAnalysis.Configuration.recoSampleDefinitionsZtoMuTau_7TeV_grid_cfi import recoSampleDefinitionsZtoMuTau_7TeV
from TauAnalysis.Configuration.submitAnalysisToGrid import submitAnalysisToGrid
from TauAnalysis.Configuration.userRegistry import getAnalysisFilePath, getJobId

channel = 'ZtoMuTau_bgEstTemplate'
configFile = 'runBgEstTemplateProductionZtoMuTau_cfg.py'
outputFilePath = getOutputFilePath(channel)
jobId = getJobId(channel)

samplesToAnalyze = [
    # modify in case you want to submit crab jobs for some of the samples only...
]

# Submit analysis jobs to the grid;
# disable estimation of systematic uncertainties for all samples except Ztautau
submitAnalysisToGrid(configFile = configFile, channel = channel,
                     samples = recoSampleDefinitionsZtoMuTau_7TeV, outputFilePath = outputFilePath, jobId = jobId,
                     samplesToAnalyze = samplesToAnalyze, samplesToSkip = ['Ztautau'], disableSysUncertainties = True)
submitAnalysisToGrid(configFile = configFile, channel = channel,
                     samples = recoSampleDefinitionsZtoMuTau_7TeV, outputFilePath = outputFilePath, jobId = jobId,
                     samplesToAnalyze = ['Ztautau'], samplesToSkip = None, disableSysUncertainties = False)
