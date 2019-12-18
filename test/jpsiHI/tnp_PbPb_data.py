import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/hidata/HIRun2018A/HISingleMuon/AOD/04Apr2019-v1/270003/21A8A05E-B3C5-4445-A76E-1433602ED7FF.root'),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.GlobalTag.globaltag = cms.string('103X_dataRun2_Prompt_fixEcalADCToGeV_v2')

## PbPb centrality bin producer
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1031x02_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])

## ==== FILTERS ====
### PbPb Event Selection
process.load("MuonAnalysis.TagAndProbe.OfflinePrimaryVerticesRecovery_cfi")
process.load('MuonAnalysis.TagAndProbe.collisionEventSelection_cff')
### Trigger selection
process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
process.triggerResultsFilter.triggerConditions = cms.vstring( 'HLT_HIL3Mu3_NHitQ10*' , 'HLT_HIL3Mu5_NHitQ10*' , 'HLT_HIL3Mu7_NHitQ10*',
                                                              'HLT_HIL3Mu12_v*' , 'HLT_HIL3Mu15_v*' , 'HLT_HIL3Mu20_v*' )
process.triggerResultsFilter.hltResults = cms.InputTag("TriggerResults","","HLT")
process.triggerResultsFilter.l1tResults = cms.InputTag("") # keep empty!
process.triggerResultsFilter.throw = False
### Filter sequence
process.fastFilter = cms.Sequence(process.triggerResultsFilter + process.offlinePrimaryVerticesRecovery + process.collisionEventSelectionAODv2)

##    __  __
##   |  \/  |_   _  ___  _ __  ___
##   | |\/| | | | |/ _ \| '_ \/ __|
##   | |  | | |_| | (_) | | | \__ \
##   |_|  |_|\__,_|\___/|_| |_|___/
##
## ==== Merge CaloMuons and Tracks into the collection of reco::Muons  ====

InAcceptance_TRACK      = 'abs(eta)<2.4 && pt>=0.5'
InAcceptance_2018_Loose = '((abs(eta)<0.3 && pt>=3.4) || (0.3<=abs(eta) && abs(eta)<1.1 && pt>=3.3) || (1.1<=abs(eta) && abs(eta)<1.4 && pt>=7.7-4.0*abs(eta)) || (1.4<=abs(eta) && abs(eta)<1.55 && pt>=2.1) || (1.55<=abs(eta) && abs(eta)<2.2 && pt>=4.25-1.39*abs(eta)) || (2.2<=abs(eta) && abs(eta)<2.4 && pt>=1.2))'

from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",
    mergeTracks = cms.bool(True),
    mergeCaloMuons = cms.bool(False), # AOD
    muons     = cms.InputTag("muons"),
    caloMuons = cms.InputTag("calomuons"),
    tracks    = cms.InputTag("generalTracks"),
    minCaloCompatibility = calomuons.minCaloCompatibility,
    ## Apply some minimal pt cut
    muonsCut     = cms.string(InAcceptance_TRACK + " && track.isNonnull"),
    caloMuonsCut = cms.string(InAcceptance_TRACK),
    tracksCut    = cms.string(InAcceptance_TRACK),
)
process.twoMuons = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("mergedMuons"), minNumber = cms.uint32(2))

## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")
## For trigger muons
switchOffAmbiguityResolution(process)
## For L3 muons
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.patTrigger.collections.remove("hltIterL3MuonCandidates")
process.patTrigger.collections.append("hltIterL3MuonCandidatesPPOnAA")
process.muonMatchHLTL3.matchedCuts = cms.string('coll("hltIterL3MuonCandidatesPPOnAA")')

## ==== Tag and probe variables
from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")
## PbPb centrality variables
from MuonAnalysis.TagAndProbe.heavyIon_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.heavyIon_modules_cff")
## Flags
### Track Id
Track_Id = "track.isNonnull && track.quality('highPurity') && (track.ptError/track.pt < 0.1) && (track.normalizedChi2/track.hitPattern.trackerLayersWithMeasurement<0.18) && (track.numberOfValidHits >= 11)"
### Muon Id
HybridSoftIdReco_2018 = "isGlobalMuon && isTrackerMuon && innerTrack.hitPattern.trackerLayersWithMeasurement > 5 && innerTrack.hitPattern.pixelLayersWithMeasurement > 0"
HybridSoftId_2018 = HybridSoftIdReco_2018 + " && abs(dB('PV2D')) < 0.3 && abs(dB('PVDZ')) < 20.0"
### Trigger
TriggerTagFlags = cms.PSet(
    # Single Muon Trigger Paths
    HLT_HIL3Mu3_NHitQ10 = cms.string("!triggerObjectMatchesByPath('HLT_HIL3Mu3_NHitQ10*_v*',1,0).empty()"),
    HLT_HIL3Mu5_NHitQ10 = cms.string("!triggerObjectMatchesByPath('HLT_HIL3Mu5_NHitQ10*_v*',1,0).empty()"),
    HLT_HIL3Mu7_NHitQ10 = cms.string("!triggerObjectMatchesByPath('HLT_HIL3Mu7_NHitQ10*_v*',1,0).empty()"),
    HLT_HIL3Mu12 = cms.string("!triggerObjectMatchesByPath('HLT_HIL3Mu12_v*',1,0).empty()"),
    HLT_HIL3Mu15 = cms.string("!triggerObjectMatchesByPath('HLT_HIL3Mu15_v*',1,0).empty()"),
    HLT_HIL3Mu20 = cms.string("!triggerObjectMatchesByPath('HLT_HIL3Mu20_v*',1,0).empty()"),
)

## ==== Tag muons
process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(InAcceptance_2018_Loose+" && "+HybridSoftId_2018+" && (!triggerObjectMatchesByPath('HLT_HIL3Mu3_NHitQ10*_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_HIL3Mu5_NHitQ10*_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_HIL3Mu7_NHitQ10*_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_HIL3Mu12_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_HIL3Mu15_v*',1,0).empty() || !triggerObjectMatchesByPath('HLT_HIL3Mu20_v*',1,0).empty())"),
)
process.oneTag = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))
process.pseudoTag = cms.EDFilter("MuonSelector",
    src = cms.InputTag("mergedMuons"),
    cut = cms.string(InAcceptance_2018_Loose+" && "+HybridSoftIdReco_2018),
)
process.onePseudoTag = process.oneTag.clone(src = cms.InputTag("pseudoTag"))

#--------------------------------------------------------------------
##    ____                   _____               _      ____            _
##   | __ )  __ _ _ __ ___  |_   _| __ __ _  ___| | __ |  _ \ _ __ ___ | |__   ___  ___
##   |  _ \ / _` | '__/ _ \   | || '__/ _` |/ __| |/ / | |_) | '__/ _ \| '_ \ / _ \/ __|
##   | |_) | (_| | | |  __/   | || | | (_| | (__|   <  |  __/| | | (_) | |_) |  __/\__ \
##   |____/ \__,_|_|  \___|   |_||_|  \__,_|\___|_|\_\ |_|   |_|  \___/|_.__/ \___||___/
##
##

process.probeTracks = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string('innerTrack.isNonnull && ' + InAcceptance_TRACK),
)
process.oneProbeTrack = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("probeTracks"), minNumber = cms.uint32(1))
process.pseudoProbeTracks = cms.EDFilter("MuonSelector",
    src = cms.InputTag("mergedMuons"),
    cut = process.probeTracks.cut
)
process.onePseudoProbeTrack = process.oneProbeTrack.clone(src = cms.InputTag("pseudoProbeTracks"))

process.tpPairsTrk = cms.EDProducer("CandViewShallowCloneCombiner",
    cut = cms.string('2.5 < mass < 3.5'),
    decay = cms.string('tagMuons@+ probeTracks@-')
)
process.onePairTrk = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairsTrk"), minNumber = cms.uint32(1))
process.pseudoPairsTrk = process.tpPairsTrk.clone(
    decay = cms.string('pseudoTag@+ pseudoProbeTracks@-')
)
process.onePseudoPairTrk = process.onePairTrk.clone(src = cms.InputTag("pseudoPairsTrk"))

process.fastPseudoTnPTrk = cms.Sequence(process.mergedMuons + process.twoMuons + process.pseudoTag + process.onePseudoTag + process.pseudoProbeTracks + process.onePseudoProbeTrack + process.pseudoPairsTrk + process.onePseudoPairTrk)

process.tpTreeTrk = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # choice of tag and probe pairs, and arbitration
    tagProbePairs = cms.InputTag("tpPairsTrk"),
    arbitration   = cms.string("OneProbe"),
    # probe variables: all useful ones
    variables = cms.PSet(
      KinematicVariables,
      TrackQualityVariables,
      GlobalTrackQualityVariables,
      StaOnlyVariables,
      dxyPVdzmin = cms.InputTag("muonDxyPVdzminTrk","dxyPVdzmin"),
      dzPV = cms.InputTag("muonDxyPVdzminTrk","dzPV"),
      dxyPV = cms.InputTag("muonDxyPVdzminTrk","dxyPV"),
      dzError = cms.string("? innerTrack.isNull() ? -1 : innerTrack.dzError"),
      dxyError = cms.string("? innerTrack.isNull() ? -1 : innerTrack.dxyError"),
      trackOrigAlgo = cms.string("? innerTrack.isNull() ? -1 : innerTrack.originalAlgo"),
      trackAlgo = cms.string("? innerTrack.isNull() ? -1 : innerTrack.algo"),
    ),
    flags = cms.PSet(
      TrackQualityFlags,
      MuonIDFlags,
      isSTA = cms.string("isStandAloneMuon"),
      TrackID = cms.string(Track_Id),
      InAcceptance_2018_Loose = cms.string(InAcceptance_2018_Loose),
    ),
    tagVariables = cms.PSet(
      KinematicVariables,
      MuonIDVariables,
      TrackQualityVariables,
      GlobalTrackQualityVariables,
      CentralityVariables,
      nVertices = cms.InputTag("nverticesModule"),
      dxyPVdzmin = cms.InputTag("muonDxyPVdzminTags","dxyPVdzmin"),
      dzPV = cms.InputTag("muonDxyPVdzminTags","dzPV"),
      dxyPV = cms.InputTag("muonDxyPVdzminTags","dxyPV"),
    ),
    tagFlags = cms.PSet(
      TriggerTagFlags,
    ),
    pairVariables = cms.PSet(
      pt = cms.string("pt"),
      y = cms.string("rapidity"),
      absy = cms.string("abs(rapidity)"),
      deltaR = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"),
      dz = cms.string("daughter(0).vz - daughter(1).vz"),
      probeMultiplicity = cms.InputTag("probeMultiplicity"),
      probeMultiplicity_TMGM = cms.InputTag("probeMultiplicityTMGM"),
    ),
    pairFlags = cms.PSet(
      BestJPsi = cms.InputTag("bestPairByJpsiMass"),
    ),
    isMC = cms.bool(False),
    addRunLumiInfo = cms.bool(True),
)

process.nverticesModule = cms.EDProducer("VertexMultiplicityCounter",
    probes = cms.InputTag("tagMuons"),
    objects = cms.InputTag("offlinePrimaryVertices"),
    objectSelection = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2")
)
process.bestPairByJpsiMass = cms.EDProducer("BestPairByMass",
    pairs = cms.InputTag("tpPairsTrk"),
    mass  = cms.double(3.096916)
)
process.probeMultiplicity = cms.EDProducer("ProbeMulteplicityProducer",
    pairs = cms.InputTag("tpPairsTrk")
)
process.probeMultiplicityTMGM = process.probeMultiplicity.clone(
    probeCut = cms.string("isTrackerMuon || isGlobalMuon")
)
process.muonDxyPVdzminTrk = cms.EDProducer("MuonDxyPVdzmin",
    probes = cms.InputTag("probeTracks")
)
process.muonDxyPVdzminTags = process.muonDxyPVdzmin.clone(
    probes = cms.InputTag("tagMuons")
)

process.tnpSimpleSequenceTrk = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeTracks +
    process.oneProbeTrack +
    process.tpPairsTrk  +
    process.onePairTrk +
    process.nverticesModule +
    process.muonDxyPVdzminTrk + process.muonDxyPVdzminTags +
    process.probeMultiplicity + process.probeMultiplicityTMGM +
    process.bestPairByJpsiMass +
    process.centralityInfo +
    process.centralityBinInfo +
    process.tpTreeTrk
)

process.tagAndProbeTrk = cms.Path(
    process.fastFilter
    * process.fastPseudoTnPTrk
    * process.centralityBin
    * process.mergedMuons * process.patMuonsWithTriggerSequence
    * process.tnpSimpleSequenceTrk
)


##    _____               _    _
##   |_   _| __ __ _  ___| | _(_)_ __   __ _
##     | || '__/ _` |/ __| |/ / | '_ \ / _` |
##     | || | | (_| | (__|   <| | | | | (_| |
##     |_||_|  \__,_|\___|_|\_\_|_| |_|\__, |
##                                     |___/

## Then make another collection for standalone muons, using standalone track to define the 4-momentum
process.muonsSta = cms.EDProducer("RedefineMuonP4FromTrack",
    src   = cms.InputTag("muons"),
    track = cms.string("outer"),
)
## Match to trigger, to measure the efficiency of HLT tracking
from PhysicsTools.PatAlgos.tools.helpers import *
process.patMuonsWithTriggerSequenceSta = cloneProcessingSnippet(process, process.patMuonsWithTriggerSequence, "Sta")
process.patMuonsWithTriggerSequenceSta.replace(process.patTriggerFullSta, process.patTriggerFull)
process.patTriggerSta.src = cms.InputTag("patTriggerFull")
process.muonMatchHLTL3Sta.maxDeltaR = 0.5
massSearchReplaceAnyInputTag(process.patMuonsWithTriggerSequenceSta, "mergedMuons", "muonsSta")

## Define probes and T&P pairs
process.probeMuonsSta = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTriggerSta"),
    cut = cms.string("outerTrack.isNonnull"), # no real cut now
)
process.oneProbeSta = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("probeMuonsSta"), minNumber = cms.uint32(1))
process.pseudoProbeSta = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muonsSta"),
    cut = process.probeMuonsSta.cut,
)
process.onePseudoProbeSta = process.oneProbeSta.clone(src = cms.InputTag("pseudoProbeSta"))

process.tpPairsSta = process.tpPairsTrk.clone(
    cut = cms.string('1. < mass < 6.'),
    decay = cms.string('tagMuons@+ probeMuonsSta@-')
)
process.onePairSta = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairsSta"), minNumber = cms.uint32(1))
process.pseudoPairsSta = process.tpPairsSta.clone(
    decay = cms.string('pseudoTag@+ pseudoProbeSta@-')
)
process.onePseudoPairSta = process.onePairSta.clone(src = cms.InputTag("pseudoPairsSta"))

process.fastPseudoTnPSta = cms.Sequence(process.mergedMuons + process.pseudoTag + process.onePseudoTag + process.muonsSta + process.pseudoProbeSta + process.onePseudoProbeSta + process.pseudoPairsSta + process.onePseudoPairSta)

process.staToTkMatch.maxDeltaR = 0.3
process.staToTkMatch.maxDeltaPtRel = 2.
process.staToTkMatchNoJPsi.maxDeltaR = 0.3
process.staToTkMatchNoJPsi.maxDeltaPtRel = 2.
process.staToTkMatchNoBestJPsi.maxDeltaR = 0.3
process.staToTkMatchNoBestJPsi.maxDeltaPtRel = 2.

process.load("MuonAnalysis.TagAndProbe.tracking_reco_info_cff")

process.tpTreeSta = process.tpTreeTrk.clone(
    tagProbePairs = "tpPairsSta",
    arbitration   = "OneProbe",
    variables = cms.PSet(
        KinematicVariables,
        StaOnlyVariables,
        ## track matching variables
        tk_deltaR = cms.InputTag("staToTkMatch","deltaR"),
        tk_deltaEta = cms.InputTag("staToTkMatch","deltaEta"),
        tk_deltaR_NoJPsi = cms.InputTag("staToTkMatchNoJPsi","deltaR"),
        tk_deltaEta_NoJPsi = cms.InputTag("staToTkMatchNoJPsi","deltaEta"),
        tk_deltaR_NoBestJPsi = cms.InputTag("staToTkMatchNoBestJPsi","deltaR"),
        tk_deltaEta_NoBestJPsi = cms.InputTag("staToTkMatchNoBestJPsi","deltaEta"),
    ),
    flags = cms.PSet(
        outerValidHits = cms.string("outerTrack.numberOfValidHits > 0"),
        TM  = cms.string("isTrackerMuon"),
        Glb = cms.string("isGlobalMuon"),
        PF  = cms.string("isPFMuon"),
        tk  = cms.string("track.isNonnull"),
        isSTA = cms.string("isStandAloneMuon"),
        StaTkSameCharge = cms.string("outerTrack.isNonnull && innerTrack.isNonnull && (outerTrack.charge == innerTrack.charge)"),
        InAcceptance_2018_Loose = cms.string(InAcceptance_2018_Loose),
    ),
    tagVariables = cms.PSet(
        KinematicVariables,
        MuonIDVariables,
        TrackQualityVariables,
        GlobalTrackQualityVariables,
        CentralityVariables,
        nVertices = cms.InputTag("nverticesModule"),
        dxyPVdzmin = cms.InputTag("muonDxyPVdzminTags","dxyPVdzmin"),
        dzPV = cms.InputTag("muonDxyPVdzminTags","dzPV"),
        dxyPV = cms.InputTag("muonDxyPVdzminTags","dxyPV"),
    ),
    tagFlags = cms.PSet(
        TriggerTagFlags,
    ),
    pairVariables = cms.PSet(
        dz = cms.string("daughter(0).vz - daughter(1).vz"),
        pt = cms.string("pt"),
        rapidity = cms.string("rapidity"),
        absy = cms.string("abs(rapidity)"),
        deltaR = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)"),
    ),
    pairFlags = cms.PSet(),
)

process.tnpSimpleSequenceSta = cms.Sequence(
    process.tagMuons +
    process.oneTag     +
    process.probeMuonsSta   +
    process.oneProbeSta     +
    process.tpPairsSta      +
    process.onePairSta      +
    process.nverticesModule +
    process.muonDxyPVdzminTags +
    process.staToTkMatchSequenceJPsi +
    process.centralityInfo +
    process.centralityBinInfo +
    process.tpTreeSta
)

if True: # turn on for tracking efficiency from RECO/AOD + earlyGeneralTracks
    process.tracksNoMuonSeeded = cms.EDFilter("TrackSelector",
      src = cms.InputTag("generalTracks"),
      cut = cms.string(" || ".join("isAlgoInMask('%s')" % a for a in [
                    'initialStep', 'lowPtTripletStep', 'pixelPairStep', 'detachedTripletStep',
                    'mixedTripletStep', 'pixelLessStep', 'tobTecStep', 'jetCoreRegionalStep' ] ) )
    )
    process.pCutTracks0 = process.pCutTracks.clone(src = 'tracksNoMuonSeeded')
    process.tkTracks0 = process.tkTracks.clone(src = 'pCutTracks0')
    process.tkTracksNoJPsi0 = process.tkTracksNoJPsi.clone(src = 'tkTracks0')
    process.tkTracksNoBestJPsi0 = process.tkTracksNoBestJPsi.clone(src = 'tkTracks0')
    process.preTkMatchSequenceJPsi.replace(
            process.tkTracksNoJPsi, process.tkTracksNoJPsi +
            process.tracksNoMuonSeeded + process.pCutTracks0 + process.tkTracks0 + process.tkTracksNoJPsi0 +process.tkTracksNoBestJPsi0
    )
    process.staToTkMatch0 = process.staToTkMatch.clone(matched = 'tkTracks0')
    process.staToTkMatchNoJPsi0 = process.staToTkMatchNoJPsi.clone(matched = 'tkTracksNoJPsi0')
    process.staToTkMatchNoBestJPsi0 = process.staToTkMatchNoBestJPsi.clone(matched = 'tkTracksNoJPsi0')
    process.staToTkMatchSequenceJPsi.replace( process.staToTkMatch, process.staToTkMatch + process.staToTkMatch0 )
    process.staToTkMatchSequenceJPsi.replace( process.staToTkMatchNoJPsi, process.staToTkMatchNoJPsi + process.staToTkMatchNoJPsi0 )
    process.staToTkMatchSequenceJPsi.replace( process.staToTkMatchNoBestJPsi, process.staToTkMatchNoBestJPsi + process.staToTkMatchNoBestJPsi0 )
    process.tpTreeSta.variables.tk0_deltaR     = cms.InputTag("staToTkMatch0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta   = cms.InputTag("staToTkMatch0","deltaEta")
    process.tpTreeSta.variables.tk0_deltaR_NoJPsi   = cms.InputTag("staToTkMatchNoJPsi0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta_NoJPsi = cms.InputTag("staToTkMatchNoJPsi0","deltaEta")
    process.tpTreeSta.variables.tk0_deltaR_NoBestJPsi   = cms.InputTag("staToTkMatchNoBestJPsi0","deltaR")
    process.tpTreeSta.variables.tk0_deltaEta_NoBestJPsi = cms.InputTag("staToTkMatchNoBestJPsi0","deltaEta")

process.tagAndProbeSta = cms.Path(
    process.fastFilter
    * process.fastPseudoTnPSta
    * process.centralityBin
    * process.mergedMuons * process.patMuonsWithTriggerSequence
    * process.muonsSta * process.patMuonsWithTriggerSequenceSta
    * process.tnpSimpleSequenceSta
)

process.schedule = cms.Schedule(
   process.tagAndProbeTrk,
   process.tagAndProbeSta,
)

process.RandomNumberGeneratorService.tkTracksNoJPsi = cms.PSet( initialSeed = cms.untracked.uint32(81) )
process.RandomNumberGeneratorService.tkTracksNoBestJPsi = cms.PSet( initialSeed = cms.untracked.uint32(81) )
process.RandomNumberGeneratorService.tkTracksNoJPsi0 = cms.PSet( initialSeed = cms.untracked.uint32(81) )
process.RandomNumberGeneratorService.tkTracksNoBestJPsi0 = cms.PSet( initialSeed = cms.untracked.uint32(81) )

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnpJpsi_Data_PbPb.root"))

from HLTrigger.Configuration.CustomConfigs import MassReplaceInputTag
process = MassReplaceInputTag(process,"offlinePrimaryVertices","offlinePrimaryVerticesRecovery")
process.offlinePrimaryVerticesRecovery.oldVertexLabel = "offlinePrimaryVertices"
