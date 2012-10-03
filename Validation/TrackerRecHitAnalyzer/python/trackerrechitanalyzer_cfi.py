import FWCore.ParameterSet.Config as cms

TrackerRecHitAnalyzer = cms.EDAnalyzer('TrackerRecHitAnalyzer',
     StripRecHitCollections = cms.VInputTag( cms.InputTag('siStripMatchedRecHits','rphiRecHit'), 
                                        cms.InputTag('siStripMatchedRecHits','stereoRecHit')
                                          ),
     PixelRecHitCollections =  cms.InputTag('siPixelRecHits')
                                
)
