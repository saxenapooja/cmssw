import FWCore.ParameterSet.Config as cms
process = cms.Process("DumpTwinMuxRaw")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = '90X_upgrade2017_realistic_v20'
process.GlobalTag.globaltag = '90X_dataRun2_HLT_v1'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.threshold  = cms.untracked.string('WARNING')

process.source = cms.Source("PoolSource",
#fileNames = cms.untracked.vstring( *('file:/afs/cern.ch/work/p/pooja/work_area/outer-hardon/emulator/CMSSW_9_2_0_patch5/src/EventFilter/L1TXRawToDigi/data/singleMuon_run296174.root',)) 
#fileNames = cms.untracked.vstring( *('file:/afs/cern.ch/work/p/pooja/work_area/outer-hardon/emulator/CMSSW_9_2_0_patch5/src/EventFilter/L1TXRawToDigi/data/singleMuon_run2017B_297114.root',))
fileNames = cms.untracked.vstring( *('file:/afs/cern.ch/work/p/pooja/work_area/outer-hardon/emulator/CMSSW_9_2_0_patch5/src/EventFilter/L1TXRawToDigi/data/singleMuon_run2017B_297505.root',))
#fileNames = cms.untracked.vstring( *('file:/afs/cern.ch/work/p/pooja/work_area/outer-hardon/emulator/CMSSW_9_2_0_patch5/src/EventFilter/L1TXRawToDigi/data/SingleMuon_run2017C_299649.root',))
)

process.twinMuxStage2Digis = cms.EDProducer("L1TTwinMuxRawToDigi",
                                       #DTTM7_FED_Source = cms.InputTag("source"),
                                       DTTM7_FED_Source = cms.InputTag("rawDataCollector"),
                                       #DTTM7_FED_Source = cms.InputTag("NewEventStreamFileReader"),
#                                       feds = cms.untracked.vint32(     1394 ),
#                                       wheels = cms.untracked.vint32(   2 ),
   feds        = cms.untracked.vint32( 1395,           1391,           1390,           1393,           1394           ),
   wheels      = cms.untracked.vint32( -2,             -1,             0,              +1,             +2             ),
   amcsecmap   = cms.untracked.vint64( 0x123456789ABC, 0x123456789ABC, 0x123456789ABC, 0x123456789ABC, 0x123456789ABC ),
   HOFirstFED  = cms.untracked.int32(700),
   firstSample = cms.int32(0),
   lastSample  = cms.int32(9),
   InputLabel  = cms.InputTag("rawDataCollector"),
#   HOFEDs      = cms.untracked.vint32(8),
#   inputHOLUTs = cms.FileInPath('EventFilter/L1TXRawToDigi/data/unpacker-LUTS-v8-Dick.csv'),
#   inputHOLUTs = cms.FileInPath('EventFilter/L1TXRawToDigi/data/unpacker-LUTS-v7-HO.csv'),
#   inputHOLUTs = cms.FileInPath('EventFilter/L1TXRawToDigi/data/HCALmapHO_l.csv'),
   inputHOLUTs = cms.FileInPath('EventFilter/L1TXRawToDigi/data/v1/HCALmapHO_K.csv'),
   ElectronicsMap = cms.string(""),
#   inputHOLUTs = cms.FileInPath('EventFilter/L1TXRawToDigi/data/unpacker-LUTS-v5.csv'),
#   feds     = cms.untracked.vint32(         1394           ),
#   wheels   = cms.untracked.vint32(             +2             ),
#   amcsecmap= cms.untracked.vint64(  0x123456789ABC ),
# 
                                      # AMC / Sector mapping example
                                       # amcsecmap                     ( 0xF9FFFF3FFFFF ),
                                       # AmcId                         (   123456789... )
                                       # Sector                        (   -9----3----- )
                                       #amcsecmap = cms.untracked.vint64( 0xFFFFFFFF9FFF ),
#                                       amcsecmap = cms.untracked.vint64( 0x123456789ABC ),
                                       debug = cms.untracked.bool(True),
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.RawToDigi = cms.Sequence(process.twinMuxStage2Digis)
process.p = cms.Path(process.RawToDigi)


process.save = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('TwinMux_2017B_run297505_SOI_bcn.root'),
                                fileName = cms.untracked.string('TwinMux_2017B_GoodDataMarkerOn.root'),
#                               fileName = cms.untracked.string('test.root'),
#                                outputCommands = cms.untracked.vstring('drop *', 
#                                    'keep L1MuDTChambPhContainer_*_*_*', 
#                                    'keep L1MuDTChambThContainer_*_*_*'),
)

#process.print1 = cms.OutputModule("AsciiOutputModule")
process.out = cms.EndPath(process.save)
