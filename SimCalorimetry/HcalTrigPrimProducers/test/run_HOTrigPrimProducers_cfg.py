import FWCore.ParameterSet.Config as cms
process = cms.Process('test')

#process.source = cms.Source(
#    'PoolSource',
#    fileNames = cms.untracked.vstring("file:hcalDigis_74X_2015D_1000Evnt.root")
#    fileNames = cms.untracked.vstring("file:hcalDigis_74X_2015D.root")
#    )

process.source = cms.Source("PoolSource",
#fileNames = cms.untracked.vstring( *('file:/nfs/dust/cms/user/pooja/scratch/hadron_outer/emulation_HO/HcalTPCode/CMSSW_9_2_0/src/SimCalorimetry/HcalTrigPrimProducers/data/hcalDigis_74X_2017A.root',))
#fileNames = cms.untracked.vstring("file:/nfs/dust/cms/user/pooja/scratch/hadron_outer/emulation_HO/HcalTPCode/CMSSW_9_2_0/src/HOTpDigiAnalyzer/HOTpDigiAnalyzer/data/hcalDigis_74X_2017B.root")
fileNames = cms.untracked.vstring("file:/nfs/dust/cms/user/pooja/scratch/hadron_outer/emulation_HO/HcalTPCode/CMSSW_9_2_0/src/SimCalorimetry/HcalTrigPrimProducers/data/hcalDigis_92X_2017B_297505.root")
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5))

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hotpdigi_cff')

process.GlobalTag.globaltag = '92X_dataRun2_Prompt_v4'

process.emulDigis = process.simHOTriggerPrimitiveDigis.clone()
process.emulDigis.inputLabel = cms.VInputTag('hcalDigis')

process.Out = cms.OutputModule("PoolOutputModule",
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p")),
#                               fileName = cms.untracked.string ("HO_EmulDigis_2017B_2975050.root")
                               fileName = cms.untracked.string ("test.root")
)

process.p = cms.Path(process.emulDigis)
process.e = cms.EndPath(process.Out)

