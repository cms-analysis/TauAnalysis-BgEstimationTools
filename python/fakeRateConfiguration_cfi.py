import FWCore.ParameterSet.Config as cms
import os
import shutil
from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants

def frproducer_name(tautype, fake_rate):
    return "tauFakeRates%s%s" % (tautype, fake_rate)

def pateff_name(fake_rate, fake_rate_config):
    return ((fake_rate_config['is_eff'] and 'eff' or 'fr') +
            fake_rate_config['cut'] + fake_rate)

from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants

from RecoTauTag.RecoTau.PFRecoTauDiscriminationByLeadingTrackPtCut_cfi import \
        pfRecoTauDiscriminationByLeadingTrackPtCut

fake_rates = {
    'shrinkingCone' : {
        'producer_name' : 'hpsTancTaus',
        'fake_rates' : {
            # WARNING DUPLICATED DIJET HIGH PT FAKE RATES TO TEST SOFTWARE
            # FIXME FIXME FIXME
            'ppMuXData' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/public/fakeRates/'
                'fakerate_qcdDiJet1st/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'DiJetHighPtdata' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/public/fakeRates/'
                'fakerate_qcdDiJet1st/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'DiJetSecondPtdata' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/public/fakeRates/'
                'fakerate_qcdDiJet2nd/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'WplusJetsdata' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/public/fakeRates/'
                'fakerate_wjets/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'ZTTsim' : {
                'is_eff' : True,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/public/fakeRates/'
                'fakerate_zttEff/fakeRate.db',
                'tag' : 'FakeRate',
            },
        },
    },
    'hpsTancTaus' : {
        'producer_name' : 'hpsTancTaus',
        'fake_rates' : {
            # MC fake rates
            'ppMuXSim' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/scratch0/fakerates/fakerate_ewkTauId_PPmuX_mcPU156bx/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'DiJetHighPtSim' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/scratch0/fakerates/fakerate_ewkTauId_QCDdiJet1st_mc/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'DiJetSecondPtSim' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/scratch0/fakerates/fakerate_ewkTauId_QCDdiJet2nd_mc/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'WplusJetsSim' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/scratch0/fakerates/fakerate_ewkTauId_WplusJets_mcPU156bx/fakeRate.db',
                'tag' : 'FakeRate',
            },
            # Data fake rates
            'ppMuXData' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/scratch0/fakerates/fakerate_ewkTauId_PPmuX_data/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'DiJetHighPtdata' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/scratch0/fakerates/fakerate_ewkTauId_QCDdiJet1st_data/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'DiJetSecondPtdata' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/scratch0/fakerates/fakerate_ewkTauId_QCDdiJet2nd_data/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'WplusJetsdata' : {
                'is_eff' : False,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/scratch0/fakerates/fakerate_ewkTauId_WplusJets_data/fakeRate.db',
                'tag' : 'FakeRate',
            },
            'ZTTsim' : {
                'is_eff' : True,
                'cut' : 'ByEWKTauID',
                'db' : '/afs/cern.ch/user/f/friis/scratch0/fakerates/fakerate_ewkTauId_Ztautau_mcPU156bx/fakeRate.db',
                'tag' : 'FakeRate',
            },
        },
    }
}

#PRODUCER = 'shrinkingCone'
PRODUCER = 'hpsTancTaus'
fake_rates_to_add = ['DiJetHighPtdata', 'DiJetSecondPtdata',
                     'ZTTsim', 'WplusJetsdata', 'ppMuXData',
                     # MC fake rates
                     'DiJetHighPtSim', 'DiJetSecondPtSim',
                     'WplusJetsSim', 'ppMuXSim',
                    ]

data_directory = os.path.join(os.environ['CMSSW_BASE'], 'src', 'TauAnalysis',
                              'BgEstimationTools', 'data')

# File-in-path version
data_directory_fip = os.path.join('TauAnalysis', 'BgEstimationTools', 'data')

def setupFakeRates(process, patProducer):
    print "Warning: Deleting any existing efficiencies!"
    patProducer.efficiencies = cms.PSet()
    patProducer.addEfficiencies = cms.bool(True)
    producer = PRODUCER
    for fake_rate in fake_rates_to_add:
        print "Making DB loader: %s fake rates for %s" % (fake_rate, producer)
        fake_rate_config = fake_rates[producer]['fake_rates'][fake_rate]
        # Copy the db to the working area
        db_source_file = fake_rate_config['db']
        local_db_file_name = producer + fake_rate + os.path.basename(db_source_file)
        local_db_copy_name = os.path.join(
            data_directory, local_db_file_name)
        if not os.path.exists(data_directory):
            os.makedirs(data_directory)
        if not os.path.exists(local_db_copy_name):
            print "Copying fakerate db:", db_source_file, \
                    "->", local_db_copy_name
            shutil.copy(db_source_file, local_db_copy_name)

        local_db_file_in_path = os.path.join(
            data_directory_fip, local_db_file_name)
        loader = cms.ESSource(
            "PoolDBESSource",
            CondDBSetup,
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(cms.PSet(
                record = cms.string('TauTagMVAComputerRcd'),
                tag = cms.string(fake_rate_config['tag'])
            )),
            connect = cms.string('sqlite_fip:%s' % local_db_file_in_path),
            BlobStreamerName = cms.untracked.string(
                'TBufferBlobStreamingService'),
            # Give it a unique name in the EventSetup
            appendToDataLabel = cms.string("%s_%s" % (producer, fake_rate)),
        )
        setattr(process, fake_rate + "LoadDB", loader)

    # Ok, databases are loaded. Now build the sequences

    # So we can make "empty" sequences
    process.dummyForFakeRate = cms.EDProducer("DummyModule")
    # Sequence that produces PFTauDiscriminators containing the fake rates
    if not hasattr(process, "tauFakeRates"):
        setattr(process, "tauFakeRates", cms.Sequence(process.dummyForFakeRate))
    for fake_rate in fake_rates_to_add:
        fake_rate_config = fake_rates[producer]['fake_rates'][fake_rate]
        prediscriminator = pfRecoTauDiscriminationByLeadingTrackPtCut.clone(
            PFTauProducer = fake_rates[producer]['producer_name'],
            Prediscriminants = noPrediscriminants,
        )
        prediscriminator_name = "computeFakeRates%s%sLeadTrk" % (producer, fake_rate)
        setattr(process, prediscriminator_name, prediscriminator)
        process.tauFakeRates += prediscriminator
        # Build our discriminator that computes the fake rate
        discriminator = cms.EDProducer(
            "RecoTauMVADiscriminator",
            PFTauProducer = cms.InputTag(fake_rates[producer]['producer_name']),
            Prediscriminants = cms.PSet(
                BooleanOperator = cms.string("and"),
                leadTrack = cms.PSet(
                    Producer = cms.InputTag(prediscriminator_name),
                    cut = cms.double(0.5),
                )
            ),
            # Point it to the correct DB instance
            dbLabel = cms.string("%s_%s" % (producer, fake_rate)),
            # We don't specify for individual decay modes
            mvas = cms.VPSet(),
            # The name in the MVA computer container.  Set in
            # trainTauFakeRate_cfg
            defaultMVA = cms.string("train"),
            remapOutput = cms.bool(False),
        )
        discriminator_name = "computeFakeRates%s%s" % (producer, fake_rate)
        setattr(process, discriminator_name, discriminator)
        process.tauFakeRates += discriminator
        # Convert the fake rate to the pat::LookupTableRecord format
        converter = cms.EDProducer(
            "PFTauDiscriminatorToPatEfficiencies",
            discSrc = cms.InputTag(discriminator_name),
            tauSrc = cms.InputTag(fake_rates[producer]['producer_name']),
        )
        converter_name = frproducer_name(producer, fake_rate)
        setattr(process, converter_name, converter)
        process.tauFakeRates += converter
        # Build the efficiency name.  The convention (defined in the code of
        # TauIdEffValidationHistManager) is (fr/eff) + CutName + Sample +
        # (data/sim)
        fake_rate_config = fake_rates[producer]['fake_rates'][fake_rate]
        pat_eff_name = pateff_name(fake_rate, fake_rate_config)
        # Finally, add this to our efficiency collection
        setattr(patProducer.efficiencies,
                pat_eff_name, cms.InputTag(converter_name))
        fake_rate_config['pateffname'] = pat_eff_name
        print "Adding efficiency:", pat_eff_name, "to pat::Tau"
