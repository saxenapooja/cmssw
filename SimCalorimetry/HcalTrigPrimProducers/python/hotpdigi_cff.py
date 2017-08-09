# The following comments couldn't be translated into the new config version:

# if read_Ascii_LUTs is true then read Ascii LUTs via "inputLUTs" below

import FWCore.ParameterSet.Config as cms

from SimCalorimetry.HcalTrigPrimProducers.hotpdigi_cfi import *
from CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi import *
HcalTPGCoderULUT = cms.ESProducer("HcalTPGCoderULUT",
#    read_Ascii_LUTs = cms.bool(False),
    read_Ascii_LUTs = cms.bool(True),
    read_XML_LUTs = cms.bool(False),
    read_FG_LUTs = cms.bool(False),
    LUTGenerationMode = cms.bool(True),
    MaskBit = cms.int32(0x8000),
#    inputLUTs = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/inputLUTcoder_physics.dat'),
    inputLUTs = cms.FileInPath('SimCalorimetry/HcalTrigPrimProducers/data/HO_ped9_inputLUTcoderDec.txt'),
    FGLUTs = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/HBHE_FG_LUT.dat'),
    RCalibFile = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/RecHit-TPG-calib.dat')
)

#import Geometry.HcalEventSetup.hcalTopologyConstants_cfi as hcalTopologyConstants_cfi
#HcalTPGCoderULUT.hcalTopologyConstants = cms.PSet(hcalTopologyConstants_cfi.hcalTopologyConstants)

HcalTrigTowerGeometryESProducer = cms.ESProducer("HcalTrigTowerGeometryESProducer")
