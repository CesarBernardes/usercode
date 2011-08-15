import FWCore.ParameterSet.Config as cms

SiStripRecHitsAnalyzer = cms.EDAnalyzer('SiStripRecHitsToDetId',
     RecHitCollections = cms.VInputTag( cms.InputTag('siStripMatchedRecHits','rphiRecHit'), 
                                        cms.InputTag('siStripMatchedRecHits','stereoRecHit')
                                       ),

     RestoredAPVCollections = cms.InputTag('siStripVRDigis','VirginRaw')
)                                     
