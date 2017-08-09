import FWCore.ParameterSet.Config as cms

simHOTriggerPrimitiveDigis = cms.EDProducer("HOTrigPrimDigiProducer",
    weights = cms.vdouble(1.0, 1.0), ##hardware algo        
    latency = cms.int32(1),
    numberOfSamples = cms.int32(4),
    numberOfPresamples = cms.int32(2),
    muonBitThreshold_high = cms.uint32(1023), # For HO muon bit
    muonBitThreshold_low = cms.uint32(40),    # For HO muon bit
    MinSignalThreshold = cms.uint32(0), # For HF PMT veto

    # Input digi label (_must_ be without zero-suppression!)

    inputLabel = cms.VInputTag(cms.InputTag('simHcalUnsuppressedDigis')),

    InputTagFEDRaw = cms.InputTag("rawDataCollector"),
    FrontEndFormatError = cms.bool(False), # Front End Format Error, for real data only
    PeakFinderAlgorithm = cms.int32(2)
)
