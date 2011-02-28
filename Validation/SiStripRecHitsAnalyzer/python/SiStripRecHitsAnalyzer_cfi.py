import FWCore.ParameterSet.Config as cms

SiStripRecHitsAnalyzer = cms.EDAnalyzer('SiStripRecHitsAnalyzer',
     RecHitCollections = cms.VInputTag( cms.InputTag('siStripMatchedRecHits','rphiRecHit'), 
                                        cms.InputTag('siStripMatchedRecHits','stereoRecHit')
                                          )

)
