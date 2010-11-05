import FWCore.ParameterSet.Config as cms
import os
import shutil
from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants

fake_rates = {
    'shrinkingCone' : {
        'producer_name' : 'shrinkingConePFTauProducer',
        'fake_rates' : {
            'qcdDijet' : {
                'db' : '/afs/cern.ch/user/v/veelken/public/fakeRate.db',
                'tag' : 'FakeRate',
            }
        },
    }
}

fake_rates_to_add = [('shrinkingCone', 'qcdDijet')]

data_directory = os.path.join(os.environ['CMSSW_BASE'], 'src', 'TauAnalysis',
                              'TauIdEfficiency', 'data')
# File-in-path version
data_directory_fip = os.path.join('TauAnalysis', 'TauIdEfficiency', 'data')

def setupFakeRates(process, patProducer):
    print "Warning: Deleting any existing efficiencies!"
    patProducer.efficiencies = cms.PSet()
    patProducer.addEfficiencies = cms.bool(True)
    for producer, fake_rate in fake_rates_to_add:
        print "Making DB loader: %s fake rates for %s" % (fake_rate, producer)
        fake_rate_config = fake_rates[producer]['fake_rates'][fake_rate]
        # Copy the db to the working area
        db_source_file = fake_rate_config['db']
        local_db_copy_name = os.path.join(
            data_directory, os.path.basename(db_source_file))
        if not os.path.exists(local_db_copy_name):
            print "Copying fakerate db:", db_source_file, \
                    "->", local_db_copy_name
            shutil.copy(db_source_file, local_db_copy_name)

        local_db_file_in_path = os.path.join(data_directory_fip,
                                             os.path.basename(db_source_file))
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
    for producer, fake_rate in fake_rates_to_add:
        # Build our discriminator that computes the fake rate
        discriminator = cms.EDProducer(
            "RecoTauMVADiscriminator",
            PFTauProducer = cms.InputTag(fake_rate_config['producer_name']),
            Prediscriminants = noPrediscriminants,
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
            tauSrc = cms.InputTag(fake_rate_config['producer_name']),
        )
        converter_name = "tauFakeRates%s%s" % (producer, fake_rate)
        setattr(process, converter_name, converter)
        process.tauFakeRates += converter
        # Finally, add this to our efficiency collection
        setattr(patProducer.efficiencies,
                fake_rate, cms.InputTag(converter_name))
