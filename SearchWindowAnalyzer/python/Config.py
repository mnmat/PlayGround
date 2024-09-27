import FWCore.ParameterSet.Config as cms
from RecoHGCal.TICL.ticlLayerTileProducer_cfi import ticlLayerTileProducer

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('PROD',Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Specify Eta

import sys
import os

eta = sys.argv[1]
energy = sys.argv[2]
nevents = sys.argv[3]
idx = sys.argv[4]
input_dir = sys.argv[5]
output_dir = sys.argv[6]

eta = eta.replace(".","")

mb = "mb_ngun"
cap = "zpos"


process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:'+ input_dir + "/" + cap + "/n" + nevents +  "/Eta_" + eta + "/singlemuon_flatEGun_hgcalCenter/step3/step3_singlemuon_e" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"))

output_dir = output_dir + "/" + cap + "/n" + nevents +  "/Eta_" + eta + "/singlemuon_flatEGun_hgcalCenter/step3/"
outfile_ ="file:" + output_dir + "searchWindowAnalyzer_singlemuon_e" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"            

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfile_),
                                   closeFileFast = cms.untracked.bool(True)
                               )
                               
# Define Process
process.ticlRecHitTile = ticlLayerTileProducer.clone(
    isLC = False)

process.ticlLayerTile = ticlLayerTileProducer.clone(
    isLC = True)

#process.ticlRecHitTile = cms.Task(ticlRecHitTile)

process.demo = cms.EDAnalyzer('SearchWindowAnalyzer',
    recHitTiles = cms.InputTag('ticlRecHitTile'),
    layerTiles = cms.InputTag('ticlLayerTile'),
    layerClusters = cms.InputTag("hgcalMergeLayerClusters"),
    hgcalRecHitsEE =  cms.InputTag("HGCalRecHit", "HGCEERecHits"),
    hgcalRecHitsFH =  cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
    hgcalRecHitsBH =  cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
    caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
    simVertices = cms.InputTag("g4SimHits"),
    lcMask = cms.InputTag("ticlTrackstersCLUE3DHigh",""),
    KFHits = cms.InputTag("ticlTrackstersKalmanFilter","KFHits","RECO"))

process.p = cms.Path(process.ticlRecHitTile*process.ticlLayerTile*process.demo)