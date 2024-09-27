import FWCore.ParameterSet.Config as cms

pcalohitAnalyzer = cms.EDAnalyzer('PCaloHitAnalyzer',
    tracks    = cms.untracked.InputTag('generalTracks'),
    #tracks    = cms.untracked.InputTag('standAloneMuons'), 
    caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
    pSimHits  = cms.untracked.InputTag('g4SimHits','MuonCSCHits'),
    pCaloHitsEE = cms.untracked.InputTag('g4SimHits','HGCHitsEE'),
    pCaloHitsFH = cms.untracked.InputTag('g4SimHits','HGCHitsHEfront'),
    pCaloHitsBH = cms.untracked.InputTag('g4SimHits','HGCHitsHEback'),
    hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
    hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    )
