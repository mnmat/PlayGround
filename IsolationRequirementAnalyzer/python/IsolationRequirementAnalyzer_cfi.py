import FWCore.ParameterSet.Config as cms

isolationRequirementAnalyzer = cms.EDAnalyzer('IsolationRequirementAnalyzer',
   caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
   Tracksters = cms.InputTag("ticlTrackstersMerge","","RECO"),
   hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
   hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
   hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
   KFHits = cms.InputTag("ticlTrackstersKalmanFilter","KFHits","RECO"),
   #PropHits = cms.InputTag("ticlTrackstersKalmanFilter","KFHits","RECO"),
   #PropHits = cms.InputTag("ticlTrackstersStandalonePropagator","KFHits","RECO"),
   hgcalLayerClusters = cms.InputTag("hgcalMergeLayerClusters", "", "RECO"),
   tracks    = cms.untracked.InputTag('generalTracks'),
   lcMask = cms.InputTag("ticlTrackstersCLUE3DHigh",""),
   associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation"),
   simVertices = cms.InputTag("g4SimHits"),
   #trackPtMin = cms.double(0.3)
)
