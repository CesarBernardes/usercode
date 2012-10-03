import FWCore.ParameterSet.Config as cms

process = cms.Process('ANALYZEDATA')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/RawToDigi_Data_cff')
process.load('RecoLocalTracker.Configuration.RecoLocalTracker_cff')
#process.load('RecoParticleFlow.Configuration.RecoParticleFlow_cff')
#process.load('RecoLocalCalo.CastorReco.CastorSimpleReconstructor_cfi')
process.load('Configuration/StandardSequences/Reconstruction_Data_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.load('Configuration/EventContent/EventContentHeavyIons_cff')
process.load('Validation.TrackerRecHits.trackerRecHitsValidation_cff')
process.load('Validation.TrackerRecHits.SiPixelRecHitsValid_cfi')
process.load('Validation.TrackerRecHits.SiPixelRecHitsValid_cfi')
process.load('Validation.TrackerRecHitAnalyzer.trackerrechitanalyzer_cfi')

process.maxEvents = cms.untracked.PSet(
  #  input = cms.untracked.int32(10)
   input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PoolSource",
    fileNames  = cms.untracked.vstring(
      #/ExpressPhysics/PARun2012-Express-v1/FEVT
      '/store/hidata/data/PARun2012/PAPhysics/RECO/PromptReco-v2/000/202/792/FCF7DAD0-1DFF-E111-9EDB-001D09F24EE3.root'
    )
)

# Other statements
process.GlobalTag.globaltag = 'GR_P_V42::All'



################################
#
#  ANALYZER
#
################################

process.TFileService=cms.Service("TFileService",fileName=cms.string('histos.root'))


process.all = cms.Path(process.siPixelRecHits*process.siStripMatchedRecHits*process.TrackerRecHitAnalyzer)

