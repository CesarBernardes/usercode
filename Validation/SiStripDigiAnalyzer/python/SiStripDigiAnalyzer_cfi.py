import FWCore.ParameterSet.Config as cms

from RecoLocalTracker.SiStripZeroSuppression.DefaultAlgorithms_cff import *
SiStripDigiAnalyzer = cms.EDAnalyzer("SiStripDigiAnalyzer",
    Algorithms = DefaultAlgorithms,
    src = cms.InputTag('siStripZeroSuppression','VirginRaw'),
	srcProcessedRawDigi = cms.InputTag('siStripZeroSuppression','VirginRaw'),
	srcAPVCM = cms.InputTag('siStripZeroSuppression','APVCM'),
    doComparisonRaw = cms.bool(True)
)
