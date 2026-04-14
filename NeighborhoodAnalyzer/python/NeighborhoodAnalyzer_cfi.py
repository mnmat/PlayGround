import FWCore.ParameterSet.Config as cms

NeighborhoodAnalyzer = cms.EDAnalyzer('NeighborhoodAnalyzer',
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