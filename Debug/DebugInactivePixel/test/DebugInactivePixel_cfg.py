import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2026_cff import Run3_2026

process = cms.Process('DBG', Run3_2026)

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


from FWCore.ParameterSet.VarParsing import VarParsing

# Allow overriding dataTier from the cmsRun command line:
# Example: `cmsRun DebugInactivePixel_cfg.py dataTier=ALCARECO`
options = VarParsing('analysis')
options.register('dataTier', 'RECO', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Data tier: RECO or ALCARECO')
options.parseArguments()
dataTier = options.dataTier

if dataTier == "RECO":
    # dataset=/ZeroBias/Run2026D-LogError-PromptReco-v1/RAW-RECO run=403461
    inputFileName = '/store/data/Run2026D/ZeroBias/RAW-RECO/LogError-PromptReco-v1/000/403/461/00000/27795d9c-f1ff-474c-b53d-710bb653337a.root'
    inputTrackCollection = 'generalTracks'
    rawErrCollcetion = 'siPixelDigis'
    maxEvents = 5
elif dataTier == "ALCARECO":
    # dataset=/Muon0/Run2026D-SiPixelCalSingleMuonTight-PromptReco-v1/ALCARECO run=403461 
    inputFileName = '/store/data/Run2026D/Muon0/ALCARECO/SiPixelCalSingleMuonTight-PromptReco-v1/000/403/461/00000/3238a69a-bc16-438f-91fa-19629a794937.root'
    inputTrackCollection = 'ALCARECOSiPixelCalSingleMuonTight'
    rawErrCollcetion = 'siPixelDigis'  # Not present in ALCARECO
    maxEvents = 5000
else:
    raise ValueError("Invalid dataTier: {}. Should be RECO / ALCARECO".format(dataTier))

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(maxEvents)
)

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(inputFileName)
)

process.DebugInactivePixel = cms.EDAnalyzer('DebugInactivePixel',
    tracksLabel   = cms.InputTag(inputTrackCollection),
    refitTracksLabel = cms.InputTag('TrackRefitter'),
    rawErrLabel   = cms.InputTag(rawErrCollcetion)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000
process.options.numberOfThreads = 1
process.options.numberOfStreams = 0

process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

if dataTier == "ALCARECO":
    process.MeasurementTrackerEvent.pixelClusterProducer = 'ALCARECOSiPixelCalSingleMuonTight'
    process.MeasurementTrackerEvent.stripClusterProducer = 'ALCARECOSiPixelCalSingleMuonTight'
    process.TrackRefitter.src = 'ALCARECOSiPixelCalSingleMuonTight'

process.TrackRefitter.src = inputTrackCollection
process.TrackRefitter.TrajectoryInEvent = True
process.TrackRefitter.TTRHBuilder = 'WithTrackAngle'

process.p = cms.Path(process.offlineBeamSpot*process.MeasurementTrackerEvent*process.TrackRefitter*process.DebugInactivePixel)
