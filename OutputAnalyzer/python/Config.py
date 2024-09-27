import FWCore.ParameterSet.Config as cms
import argparse

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('PROD',Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Specify Eta

import sys
import os

eta = str(27)
energy = str(100)
nevents = str(100)
idx = str(1)
input_dir = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/Output/0_PU/"

mb = "mb_ngun"
cap = "zpos"

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:'+ input_dir + "/" + cap + "/n" + nevents +  "/Eta_" + eta + "/singlemuon_flatEGun_hgcalCenter/step3/step3_singlemuon_e" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"))


#outfile_ = 'file:/eos/home-m/mmatthew/Data/deleteme.root'
#fname = '/eos/home-m/mmatthew/Data/Analyzer/UpdatorStudies/'+propagator+'/' + cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'
#fname = '/eos/home-m/mmatthew/Data/KF/MaterialBudget/Radlen/0_25/'+ cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'


process.demo = cms.EDAnalyzer('OutputAnalyzer',
   tracks    = cms.untracked.InputTag('ticlTrackstersKalmanFilter','HGCALTracks'),
   tracksters = cms.untracked.InputTag('ticlTrackstersKalmanFilter'),
   HGCEEInput =  cms.InputTag("HGCalRecHit", "HGCEERecHits"),
   HGCFHInput =  cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
   HGCBHInput =  cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
)


process.p = cms.Path(process.demo)

