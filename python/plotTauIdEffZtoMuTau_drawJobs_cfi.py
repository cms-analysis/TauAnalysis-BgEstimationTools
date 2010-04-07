import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.DQMTools.tools.drawJobConfigurator import *

# define template for all kinds of plots
# (specific to tau id. efficiency measurement in Z --> mu + tau-jet events)
plots_TauIdEffZtoMuTau = cms.PSet(
    plots = cms.PSet(  
        dqmMonitorElements = cms.vstring(''),
        processes = cms.vstring(
            'Zmumu',
            'WplusJets',
            'TTplusJets',
            'qcdSum',
            'Ztautau'
        )
    ),
    xAxis = cms.string('unlabeled'),
    yAxis = cms.string('numEntries_linear'),
    #yAxis = cms.string('numEntries_log'),
    legend = cms.string('regular'),
    labels = cms.vstring('mcNormScale'),                   
    drawOptionSet = cms.string('default'),
    ##stack = cms.vstring(
    ##    'Zmumu',
    ##    'WplusJets',
    ##    'TTplusJets',
    ##    'qcdSum',
    ##    'Ztautau'
    ##)
)

drawJobConfigurator_TauIdEffZtoMuTau = drawJobConfigurator(
    template = plots_TauIdEffZtoMuTau,
    dqmDirectory = '#PROCESSDIR#/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/'
)

#--------------------------------------------------------------------------------
# define distributions to be plotted
# for events passing all event selection criteria
#--------------------------------------------------------------------------------

drawJobConfigurator_TauIdEffZtoMuTau.add(
    afterCut = "uniqueTauCandidateCutTauIdEffZtoMuTau",
    plots = [
        drawJobConfigEntry(
            meName = 'MuonQuantities/Muon#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Muon (all Events passing Selection)",
            xAxis = '#PAR#',
            name = "muon"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonTrkIsoPt',
            title = "Muon Track Isolation (all Events passing Selection)",
            xAxis = 'Pt',
            name = "muonTrackIsoPt"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonEcalIsoPt',
            title = "Muon ECAL Isolation (all Events passing Selection)",
            xAxis = 'Pt',
            name = "muonEcalIsoPt"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonIsoSumPt',
            title = "Muon Track+ ECAL Isolation (all Events passing Selection)",
            xAxis = 'Pt',
            name = "muonIsoSumPt"
        ),
        drawJobConfigEntry(
            meName = 'MuonQuantities/MuonMatchingGenParticlePdgId',
            title = "PdgId of gen. Particle matching Muon (all Events passing Selection)",
            xAxis = 'PdgId',
            name = "pdgIdGenParticleMatchingMuon"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/Tau#PAR#',
            PAR = [ 'Pt', 'Eta', 'Phi' ],
            title = "Tau (all Events passing Selection)",
            xAxis = '#PAR#',
            name = "tau"
        ),
        drawJobConfigEntry(
            meName = 'TauQuantities/TauMatchingGenParticlePdgId',
            title = "PdgId of gen. Particle matching Tau (all Events passing Selection)",
            xAxis = 'PdgId',
            name = "pdgIdGenParticleMatchingTau"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/VisMass',
            title = "M_{vis}(Muon + Tau) (final Event sample)",
            xAxis = 'Mass',
            name = "mVisible"
        ),
        drawJobConfigEntry(
            meName = 'DiTauCandidateQuantities/CollinearApproxMass',
            title = "M(Muon + Tau), collinear Approx. (final Event sample)",
            xAxis = 'Mass',
            name = "mCollApprox"
        )
    ]
)

#--------------------------------------------------------------------------------
# define distributions to be plotted
# for separately for events in the two different regions:
#  o tau id. passed
#  o tau id. failed#
# ("tau id." = leading track finding && leading track Pt cut && track isolation && ECAL isolation 
#             && 1||3 tracks in signal cone && charge +1||-1)
#--------------------------------------------------------------------------------

for iRegion in [ 1, 2]:
    
    dqmSubDirectory_region = 'tauIdEffHistograms2regions/region' + "%(i)02d" % {"i" : iRegion} + '/'
    
    title_region = { 1: "Tau id. failed",
                     2: "Tau id. passed" }[iRegion]
    
    name_region = "region%(i)02d" % {"i" : iRegion}
            
    drawJobConfigurator_TauIdEffZtoMuTau.add(
        afterCut = "uniqueTauCandidateCutTauIdEffZtoMuTau",
        plots = [        
            ##drawJobConfigEntry(
            ##    meName = dqmSubDirectory_region + 'MuonQuantities/Muon#PAR#',
            ##    PAR = [ 'Pt', 'Eta', 'Phi' ],
            ##    title = "Muon (" + title_region + ")",
            ##    xAxis = '#PAR#',
            ##    name = name_region + "_muon"
            ##),
            ##drawJobConfigEntry(
            ##    meName = dqmSubDirectory_region + 'MuonQuantities/MuonTrkIsoPt',
            ##    title = "Muon Track Isolation (" + title_region + ")",
            ##    xAxis = 'Pt',
            ##    name = name_region + "_muonTrackIsoPt"
            ##),
            drawJobConfigEntry(
                meName = dqmSubDirectory_region + 'TauIdEffSpecificQuantities/MuonPt',
                title = "Muon P_{T} (" + title_region + ")",
                xAxis = 'Pt',
                name = name_region + "_muonPt"
            ),
            drawJobConfigEntry(
                meName = dqmSubDirectory_region + 'TauIdEffSpecificQuantities/MuonAbsEta',
                title = "Muon |#eta| (" + title_region + ")",
                xAxis = 'Eta',
                name = name_region + "_muonAbsEta"
            ),
            drawJobConfigEntry(
                meName = dqmSubDirectory_region + 'TauIdEffSpecificQuantities/MuonExtTrkIsoPt',
                title = "Muon Track Isolation in ext. #DeltaR = 1.0 Cone (" + title_region + ")",
                xAxis = 'Pt',
                name = name_region + "_muonExtTrackIso"
            ),
            drawJobConfigEntry(
                meName = dqmSubDirectory_region + 'MuonQuantities/MuonEcalIsoPt',
                title = "Muon ECAL Isolation (" + title_region + ")",
                xAxis = 'Pt',
                name = name_region + "_muonEcalIsoPt"
            ),
            drawJobConfigEntry(
                meName = dqmSubDirectory_region + 'MuonQuantities/MuonIsoSumPt',
                title = "Muon Track + ECAL Isolation (" + title_region + ")",
                xAxis = 'Pt',
                name = name_region + "_muonIsoSumPt"
            ),
            drawJobConfigEntry(
                meName = dqmSubDirectory_region + 'MuonQuantities/MuonMatchingGenParticlePdgId',
                title = "PdgId of gen. Particle matching Muon (" + title_region + ")",
                xAxis = 'PdgId',
                name = name_region + "_pdgIdGenParticleMatchingMuon"
            ),
            drawJobConfigEntry(
                meName = dqmSubDirectory_region + 'TauQuantities/Tau#PAR#',
                PAR = [ 'Pt', 'Eta', 'Phi' ],
                title = "Tau (" + title_region + ")",
                xAxis = '#PAR#',
                name = name_region + "_tau"
            ),
            drawJobConfigEntry(
                meName = dqmSubDirectory_region + 'TauQuantities/TauMatchingGenParticlePdgId',
                title = "PdgId of gen. Particle matching Tau (" + title_region + ")",
                xAxis = 'PdgId',
                name = name_region + "_pdgIdGenParticleMatchingTau"
            ),
            drawJobConfigEntry(
                meName = 'DiTauCandidateQuantities/VisMass',
                title = "M_{vis}(Muon + Tau) (final Event sample)",
                xAxis = 'Mass',
                name = name_region + "_mVisible"
            ),
            drawJobConfigEntry(
                meName = 'DiTauCandidateQuantities/CollinearApproxMass',
                title = "M(Muon + Tau), collinear Approx. (final Event sample)",
                xAxis = 'Mass',
                name = name_region + "_mCollApprox"
            )
        ]
    )


#--------------------------------------------------------------------------------
# plot shapes of the muon track, ECAL and combined isolation distributions
# compared for the two different regions
#--------------------------------------------------------------------------------

plots_TauIdEffZtoMuTau_shapes = cms.PSet(
    plots = cms.VPSet(
        cms.PSet(
            dqmMonitorElements = cms.vstring(),
            process = cms.string(''),
            drawOptionEntry = cms.string('region01'),
            legendEntry = cms.string('')
        ),
        cms.PSet(
            dqmMonitorElements = cms.vstring(),
            process = cms.string(''),
            drawOptionEntry = cms.string('region02'),
            legendEntry = cms.string('')
        ),
        cms.PSet(
            dqmMonitorElements = cms.vstring(),
            process = cms.string(''),
            drawOptionEntry = cms.string('region05'),
            legendEntry = cms.string('')
        ),
        cms.PSet(
            dqmMonitorElements = cms.vstring(),
            process = cms.string(''),
            drawOptionEntry = cms.string('region06'),
            legendEntry = cms.string('')
        )
    ),    
    title = cms.string(''),
    xAxis = cms.string('unlabeled'),
    yAxis = cms.string('numEntries_linear'),
    #yAxis = cms.string('numEntries_log'),        
    legend = cms.string('regular'),
    labels = cms.vstring()
)

drawJobs_TauIdEffZtoMuTau_shapes = cms.PSet()

def addDrawJob(drawJobs, meName, xAxis, label):
    for process in [ 'Zmumu', 'WplusJets', 'qcdSum', 'Ztautau' ]:
        drawJob_process = copy.deepcopy(plots_TauIdEffZtoMuTau_shapes)

        for iRegion in [ 1, 2 ]:
        
            index = { 1: 0,
                      2: 1 }[iRegion]
        
            dqmDirectory = 'harvested/' + process + '/TauIdEffAnalyzerZtoMuTau_absMuonIsolation/afterUniqueTauCandidateCutTauIdEffZtoMuTau/'    
            dqmSubDirectory_region = 'tauIdEffHistograms2regions/region' + "%(i)02d" % {"i" : iRegion} + '/'
            meName_full = dqmDirectory + dqmSubDirectory_region + meName
            drawJob_process.plots[index].dqmMonitorElements = cms.vstring(meName_full)
            
            drawJob_process.plots[index].process = cms.string(process)
            
            title_region = { 1: "Tau id. failed",
                             2: "Tau id. passed" }[iRegion]
            title = title_region
            drawJob_process.plots[index].legendEntry = cms.string(title)
            
        drawJob_process.title = cms.string(label + "_" + process)
        drawJob_process.xAxis = cms.string(xAxis)
        setattr(drawJob_process, "norm", cms.double(1.))
            
        setattr(drawJobs, label + "_" + process, drawJob_process)

addDrawJob(drawJobs_TauIdEffZtoMuTau_shapes, 'TauIdEffSpecificQuantities/MuonPt', "Pt", "muonPt")
addDrawJob(drawJobs_TauIdEffZtoMuTau_shapes, 'TauIdEffSpecificQuantities/MuonAbsEta', "Eta", "muonAbsEta")
addDrawJob(drawJobs_TauIdEffZtoMuTau_shapes, 'TauIdEffSpecificQuantities/MuonExtTrkIsoPt', "Pt", "muonExtTrackIso")
addDrawJob(drawJobs_TauIdEffZtoMuTau_shapes, 'DiTauCandidateQuantities/VisMass', "Mass", "mVisible")
addDrawJob(drawJobs_TauIdEffZtoMuTau_shapes, 'DiTauCandidateQuantities/CollinearApproxMass', "Mass", "mCollApprox")
addDrawJob(drawJobs_TauIdEffZtoMuTau_shapes, 'TauQuantities/TauPt', "Pt", "tauPt")
addDrawJob(drawJobs_TauIdEffZtoMuTau_shapes, 'TauQuantities/TauEta', "Eta", "tauEta")

