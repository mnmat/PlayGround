import FWCore.ParameterSet.Config as cms

eta = "17_27"
nevents = "100"
energy = sys.argv[1]
idx = sys.argv[2]
mb = "mb_ngun"
cap = "zpos"

input_dir = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.2/InputTracks/200_PU/"
output_dir = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.2/InputTracks/200_PU/"


infile_ = "file:"+ input_dir + "/" + cap + "/n" + nevents +  "/Eta_" + eta + "/singlemuon_flatEGun_hgcalCenter/step3/step3_singlemuon_e" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"
outfile_ = output_dir + "/" + cap + "/n" + nevents +  "/Eta_" + eta + "/singlemuon_flatEGun_hgcalCenter/step3/neighborhoodanalyzer_singlemuon_e" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"



from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("Demo",Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(infile_))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfile_),
                                   closeFileFast = cms.untracked.bool(True)
                               )


process.demo = cms.EDAnalyzer('NeighborhoodAnalyzer',
   caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
   simclusters = cms.InputTag("mix","MergedCaloTruth"),
   hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
   hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
   hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
   pCaloHitsEE = cms.untracked.InputTag('g4SimHits','HGCHitsEE'),
   pCaloHitsFH = cms.untracked.InputTag('g4SimHits','HGCHitsHEfront'),
   pCaloHitsBH = cms.untracked.InputTag('g4SimHits','HGCHitsHEback'),
   associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation"),
   simVertices = cms.InputTag("g4SimHits"),
)

process.p = cms.Path(process.demo)
process.add_(cms.Service("AdaptorConfig", native=cms.untracked.vstring("root")))

