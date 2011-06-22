import FWCore.ParameterSet.Config as cms

process = cms.Process('RAWSkim')

# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContentHeavyIons_cff')


process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.232.2.6.2.2 $'),
    annotation = cms.untracked.string('test nevts:2'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
#process.MessageLogger.cerr.FwkReport.reportEvery = 100

#process.Timing = cms.Service("Timing")

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'file:ROOTFiles/00DFBAEF-5741-E011-B023-0025901D6486.root'
    'file:ROOTFiles/FE9995AE-D8FF-DF11-AB4D-003048F024DE_Core_RAW.root',
    'file:ROOTFiles/84E1770D-D4FF-DF11-BDF9-001D09F2AF96_Core_RAW.root',
    'file:ROOTFiles/86FAD73C-B6FF-DF11-8E9D-0030487C608C_Core_RAW.root',
	'file:ROOTFiles/90A4BCED-D0FF-DF11-853B-0030487C2B86_Core_RAW.root',
	'file:ROOTFiles/F0DD858D-DDFF-DF11-95F6-0030487C90EE_Core_RAW.root'

    )
)

process.options = cms.untracked.PSet(

)

############ Filters ###########


process.load('L1Trigger.Skimmer.l1Filter_cfi')
process.DoubleMuOpen = process.l1Filter.clone(
   algorithms = cms.vstring('L1_DoubleMuOpen_BptxAND')
)

process.SingleMu3 = process.l1Filter.clone(
   algorithms = cms.vstring('L1_SingleMu3_BptxAND')
)

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")

### MinBias SD
process.hltMBHI = process.hltHighLevel.clone(HLTPaths = ['HLT_HIMinBiasHfOrBSC_Core'])
process.filterSdMBHI = cms.Path(process.hltMBHI)

### JetHI SD
process.hltJetHI = process.hltHighLevel.clone(HLTPaths = ['HLT_HIJet35U_Core'])
process.filterSdJetHI = cms.Path(process.hltJetHI)

### PhotonHI SD
process.hltPhotonHI = process.hltHighLevel.clone(HLTPaths = ['HLT_HIPhoton15_Cleaned_Core'])
process.filterSdPhotonHI = cms.Path(process.hltPhotonHI)

### MuHIL2DoubleMu3 SD
#process.hltL2DoubleMu3HI = process.hltHighLevel.clone(HLTPaths = ['HLT_HIL2DoubleMu3_Core'])
#process.filterSdL2DoubleMu3HI = cms.Path(process.RawToDigi*process.hltMBHI*process.hltL2DoubleMu3HI)

### MuHIL1_DoubleMuOpen SD
process.filterSdL1DoubleMuOpenHI = cms.Path(process.RawToDigi*process.hltMBHI*process.DoubleMuOpen)

### MuHIL1_SingleMu3 SD
process.filterSdL1SingleMu3HI = cms.Path(process.RawToDigi*process.hltMBHI*process.SingleMu3)

############ Output Modules ##########
myEvContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *_*_*_RAWSkim')
    )

### MinBias SD
process.outputSdMBHI = cms.OutputModule("PoolOutputModule",
										splitLevel = cms.untracked.int32(0),
                                         SelectEvents = cms.untracked.PSet(
   										 	SelectEvents = cms.vstring('filterSdMBHI')),                               
                                         dataset = cms.untracked.PSet(
    										dataTier = cms.untracked.string('RAW'),
    										filterName = cms.untracked.string('SD_MBHI')),
                                         outputCommands = process.RAWEventContent.outputCommands,
                                           	fileName = cms.untracked.string('SD_MBHI.root')
                                         )
process.outputSdMBHI.outputCommands.extend(myEvContent.outputCommands)

### JetHI SD
process.outputSdJetHI = cms.OutputModule("PoolOutputModule",
										splitLevel = cms.untracked.int32(0),
                                         SelectEvents = cms.untracked.PSet(
   										 	SelectEvents = cms.vstring('filterSdJetHI')),                               
                                         dataset = cms.untracked.PSet(
    										dataTier = cms.untracked.string('RAW'),
    										filterName = cms.untracked.string('SD_JetHI')),
                                         outputCommands = process.RAWEventContent.outputCommands,
                                           	fileName = cms.untracked.string('SD_Jet35HI.root')
                                         )
                                         
process.outputSdJetHI.outputCommands.extend(myEvContent.outputCommands)

### PhotonHI SD
process.outputSdPhotonHI = cms.OutputModule("PoolOutputModule",
										splitLevel = cms.untracked.int32(0),
                                         SelectEvents = cms.untracked.PSet(
   										 	SelectEvents = cms.vstring('filterSdPhotonHI')),                               
                                         dataset = cms.untracked.PSet(
    										dataTier = cms.untracked.string('RAW'),
    										filterName = cms.untracked.string('SD_PhotonHI')),
                                         outputCommands = process.RAWEventContent.outputCommands,
                                           	fileName = cms.untracked.string('SD_Photon15HI.root')
                                         )

process.outputSdPhotonHI.outputCommands.extend(myEvContent.outputCommands)



### MuHIL1SingleMu3 SD
process.outputSdL1SingleMu3HI = cms.OutputModule("PoolOutputModule",
										splitLevel = cms.untracked.int32(0),
                                         SelectEvents = cms.untracked.PSet(
   										 	SelectEvents = cms.vstring('filterSdL1SingleMu3HI')),                               
                                         dataset = cms.untracked.PSet(
    										dataTier = cms.untracked.string('RAW'),
    										filterName = cms.untracked.string('SD_L1SingleMu3HI')),
                                         outputCommands = process.RAWEventContent.outputCommands,
                                           	fileName = cms.untracked.string('SD_L1SingleMu3HI.root')
                                         )
                                         
process.outputSdL1SingleMu3HI.outputCommands.extend(myEvContent.outputCommands)                                         

### MuHIL1_DoubleMuOpen SD
process.outputSdL1DoubleMuOpenHI = cms.OutputModule("PoolOutputModule",
										splitLevel = cms.untracked.int32(0),
                                         SelectEvents = cms.untracked.PSet(
   										 	SelectEvents = cms.vstring('filterSdL1DoubleMuOpenHI')),                               
                                         dataset = cms.untracked.PSet(
    										dataTier = cms.untracked.string('RAW'),
    										filterName = cms.untracked.string('SD_L1DoubleMuOpenHI')),
                                         outputCommands = process.RAWEventContent.outputCommands,
                                            fileName = cms.untracked.string('SD_L1DoubleMuOpenHI.root')
                                         )

process.outputSdL1DoubleMuOpenHI.outputCommands.extend(myEvContent.outputCommands)

# Output definition

# process.RECOoutput = cms.OutputModule("PoolOutputModule",
#     splitLevel = cms.untracked.int32(0),
# #    outputCommands = process.RECOEventContent.outputCommands,
#     fileName = cms.untracked.string('test_RAW2DIGI_RECO_DQM_ALCA_REPACK.root'),
#     dataset = cms.untracked.PSet(
#         filterName = cms.untracked.string(''),
#         dataTier = cms.untracked.string('RECO')
#     )
# )
# 
# process.REPACKRAWoutput = cms.OutputModule("PoolOutputModule",
#     splitLevel = cms.untracked.int32(0),
#     outputCommands = process.REPACKRAWEventContent.outputCommands,
#     fileName = cms.untracked.string('test_RAW2DIGI_RECO_DQM_ALCA_REPACK_inREPACKRAW.root'),
#     dataset = cms.untracked.PSet(
#         filterName = cms.untracked.string(''),
#         dataTier = cms.untracked.string('RAW')
#     )
# )



# Other statements
process.GlobalTag.globaltag = 'GR_R_39X_V6B::All'


# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)

#process.reconstruction_step = cms.Path(process.reconstructionHeavyIons)

#process.digi2repack_step = cms.Path(process.DigiToSplitRawRepack)

#process.endjob_step = cms.EndPath(process.endOfProcess)

#process.RECOoutput_step = cms.EndPath(process.RECOoutput)

#process.REPACKRAWoutput_step = cms.EndPath(process.REPACKRAWoutput)


process.MB_output_step = cms.EndPath(process.outputSdMBHI)
process.Jet_output_step = cms.EndPath(process.outputSdJetHI)
process.Photon_output_step = cms.EndPath(process.outputSdPhotonHI)
process.SingleMu3_output_step = cms.EndPath(process.outputSdL1SingleMu3HI)
process.DoubleMuOpen_output_step = cms.EndPath(process.outputSdL1DoubleMuOpenHI)

#process.schedule = cms.Schedule(process.raw2digi_step,  process.MB_output_step, process.Jet_output_step, process.Photon_output_step, process.DoubleMu3_output_step, process.DoubleMuOpen_output_step)
#process.sequence = cms.Sequence(process.RawToDigi*process.outputSdMBHI*process.outputSdJetHI*process.outputSdPhotonHI*process.RECOoutput)