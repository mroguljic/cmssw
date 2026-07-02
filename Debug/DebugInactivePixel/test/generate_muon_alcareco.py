import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2026_cff import Run3_2026
from RecoVertex.BeamSpotProducer.BeamSpot_cff import *

process = cms.Process("ALCARECO")
process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "160X_dataRun3_Prompt_v1", '')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "/store/data/Run2026D/Muon0/RAW-RECO/LogError-PromptReco-v1/000/403/461/00000/a04b8e7e-6ee6-4392-ad06-036c1d5afcc5.root"
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)

process.load(
    "Calibration.TkAlCaRecoProducers.ALCARECOSiPixelCalSingleMuonTight_cff"
)
process.load(
    "Calibration.TkAlCaRecoProducers.ALCARECOSiPixelCalSingleMuonTight_Output_cff"
)


process.pathALCARECOSiPixelCalSingleMuonTight = cms.Path(
    process.seqALCARECOSiPixelCalSingleMuonTight
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    process.OutALCARECOSiPixelCalSingleMuonTight,
    fileName=cms.untracked.string("alcareco.root"),
)

process.endpathALCARECOSiPixelCalSingleMuonTight = cms.EndPath(
    process.out
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)