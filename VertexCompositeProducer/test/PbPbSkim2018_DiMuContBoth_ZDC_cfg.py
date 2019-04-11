import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('ANASKIM',eras.Run2_2018_pp_on_AA)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# Define the input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/hidata/HIRun2018A/HIDoubleMuon/AOD/PromptReco-v2/000/327/148/00000/37F475FE-F3F1-0749-8455-C9DB3177598B.root'),
   inputCommands=cms.untracked.vstring('keep *', 'drop *_hiEvtPlane_*_*')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(200))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string('103X_dataRun2_v6')

# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
#process.hiCentrality.srcZDChits = cms.InputTag("QWzdcreco")
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)

# Add PbPb event plane
process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
process.load("RecoHI.HiEvtPlaneAlgos.hiEvtPlaneFlat_cfi")
process.hiEvtPlane.trackTag = cms.InputTag("generalTracks")
process.hiEvtPlane.vertexTag = cms.InputTag("offlinePrimaryVerticesRecovery")
process.hiEvtPlane.loadDB = cms.bool(True)
process.hiEvtPlane.useNtrk = cms.untracked.bool(False)
process.hiEvtPlane.caloCentRef = cms.double(-1)
process.hiEvtPlane.caloCentRefWidth = cms.double(-1)
process.hiEvtPlaneFlat.caloCentRef = cms.double(-1)
process.hiEvtPlaneFlat.caloCentRefWidth = cms.double(-1)
process.hiEvtPlaneFlat.vertexTag = cms.InputTag("offlinePrimaryVerticesRecovery")
process.hiEvtPlaneFlat.useNtrk = cms.untracked.bool(False)
process.CondDB.connect = "sqlite_file:HeavyIonRPRcd_PbPb2018_offline.db"
process.PoolDBESSource = cms.ESSource("PoolDBESSource",
                                       process.CondDB,
                                       toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
                                                                  tag = cms.string('HeavyIonRPRcd_PbPb2018_offline')
                                                                  )
                                                         )
                                      )
process.es_prefer_flatparms = cms.ESPrefer('PoolDBESSource','')
process.evtplane_seq = cms.Sequence(process.hiEvtPlane * process.hiEvtPlaneFlat)

# Add the VertexComposite producer
process.load("VertexCompositeAnalysis.VertexCompositeProducer.generalDiMuCandidates_cff")
process.generalMuMuMassMin7CandidatesWrongSign = process.generalMuMuMassMin7Candidates.clone(isWrongSign = cms.bool(True))
from VertexCompositeAnalysis.VertexCompositeProducer.PATAlgos_cff import doPATMuons
doPATMuons(process, False)

# Add muon event selection
process.goodMuon = cms.EDFilter("MuonSelector",
            src = cms.InputTag("muons"),
            cut = process.generalMuMuMassMin7Candidates.muonSelection,
            )
process.twoGoodMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodMuon"), minNumber = cms.uint32(2))
process.goodDimuon = cms.EDProducer("CandViewShallowCloneCombiner",
            cut = process.generalMuMuMassMin7Candidates.candidateSelection,
            checkCharge = cms.bool(False),
            decay = cms.string('goodMuon@+ goodMuon@-')
            )
process.oneGoodDimuon = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("goodDimuon"), minNumber = cms.uint32(1))
process.dimuonEvtSel = cms.Sequence(process.goodMuon * process.twoGoodMuons * process.goodDimuon * process.oneGoodDimuon)

# Add trigger selection
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilter.andOr = cms.bool(True)
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = [
    # Double muon triggers
    'HLT_HIL1DoubleMuOpen_OS_Centrality_40_100_v*', # Peripheral OS dimuons
    'HLT_HIL1DoubleMuOpen_Centrality_50_100_v*', # Peripheral dimuons
    'HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v*', # Bottomonia
    'HLT_HIL1DoubleMu10_v*', # Z bosons
    # Single muon triggers
    'HLT_HIL1MuOpen_Centrality_80_100_v*', # Peripheral muons
    'HLT_HIL3Mu12_v*', # Electroweak bosons
    ]

# Add PbPb collision event selection
process.load('VertexCompositeAnalysis.VertexCompositeProducer.collisionEventSelection_cff')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.clusterCompatibilityFilter_cfi')
process.load('VertexCompositeAnalysis.VertexCompositeProducer.hfCoincFilter_cff')
process.load("RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesRecovery_cfi")
process.colEvtSel = cms.Sequence(process.offlinePrimaryVerticesRecovery * process.hfCoincFilter2Th4 * process.primaryVertexFilter * process.clusterCompatibilityFilter)

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.hltFilter *
    process.dimuonEvtSel *
    process.colEvtSel
)
process.eventFilter_HM_step = cms.Path( process.eventFilter_HM )

#-----------------------------------------
# CMSSW/Hcal Related Module import
#-----------------------------------------
process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")

#set digi and analyzer
rawTag = ''
process.hcalDigis.InputLabel = rawTag

if rawTag == '':
    process.digis = cms.Sequence()
else:
    process.digis = cms.Sequence(process.hcalDigis)

# ZDC info
process.load('QWAna.QWZDC2018RecHit.QWZDC2018Producer_cfi')
process.load('QWAna.QWZDC2018RecHit.QWZDC2018RecHit_cfi')

process.zdc_step = cms.Path(
    process.eventFilter_HM *
    process.digis *
    process.zdcdigi *
    process.QWzdcreco
)

# Define the analysis steps
process.pcentandep_step = cms.Path(process.eventFilter_HM * process.cent_seq * process.evtplane_seq)
process.dimurereco_step = cms.Path(process.eventFilter_HM * process.patMuonSequence * process.generalMuMuMassMin7Candidates)
process.dimurerecowrongsign_step = cms.Path(process.eventFilter_HM * process.patMuonSequence * process.generalMuMuMassMin7CandidatesWrongSign)

# Define the output
process.load("VertexCompositeAnalysis.VertexCompositeProducer.ppanalysisSkimContentJPsi_cff")
process.output_HM = cms.OutputModule("PoolOutputModule",
    outputCommands = process.analysisSkimContent.outputCommands,
    fileName = cms.untracked.string('PbPb_DiMuCont.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('eventFilter_HM_step')),
)
process.output_HM_step = cms.EndPath(process.output_HM)

# Define the process schedule
process.schedule = cms.Schedule(
    process.eventFilter_HM_step,
    process.zdc_step,
    process.pcentandep_step,
    process.dimurereco_step,
    process.dimurerecowrongsign_step,
    process.output_HM_step
)

# Add the event selection filters
process.Flag_colEvtSel = cms.Path(process.colEvtSel)
process.Flag_hfCoincFilter2Th4 = cms.Path(process.hfCoincFilter2Th4)
process.Flag_primaryVertexFilter = cms.Path(process.primaryVertexFilter)
process.Flag_clusterCompatibilityFilter = cms.Path(process.clusterCompatibilityFilter)
eventFilterPaths = [ process.Flag_colEvtSel , process.Flag_hfCoincFilter2Th4 , process.Flag_primaryVertexFilter , process.Flag_clusterCompatibilityFilter ]
for P in eventFilterPaths:
    process.schedule.insert(0, P)

# Add recovery for offline primary vertex
from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"
