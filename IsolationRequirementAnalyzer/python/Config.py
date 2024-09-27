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

eta = sys.argv[1]
energy = sys.argv[2]
nevents = sys.argv[3]
idx = sys.argv[4]
#input_dir = sys.argv[6]
input_dir = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/samples/200_PU/"

eta = eta.replace(".","")

cap = "zpos"

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:'+ input_dir + "/" + cap + "/n" + nevents +  "/Eta_" + eta + "/singlemuon_flatEGun_hgcalCenter/step3/step3_singlemuon_e" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"))


#outfile_ = 'file:/eos/home-m/mmatthew/Data/deleteme.root'
#fname = '/eos/home-m/mmatthew/Data/Analyzer/UpdatorStudies/'+propagator+'/' + cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'
#fname = '/eos/home-m/mmatthew/Data/KF/MaterialBudget/Radlen/0_25/'+ cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'
output_dir = input_dir + "/" + cap + "/n" + nevents +  "/Eta_" + eta + "/singlemuon_flatEGun_hgcalCenter/step3/"
outfile_ ="file:" + output_dir + "isolationRequirements_singlemuon_e" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"            

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfile_),
                                   closeFileFast = cms.untracked.bool(True)
                               )

from PlayGround.IsolationRequirementAnalyzer.IsolationRequirementAnalyzer_cfi import isolationRequirementAnalyzer
process.demo = isolationRequirementAnalyzer.clone()
process.p = cms.Path(process.demo)

process.add_(cms.Service("AdaptorConfig", native=cms.untracked.vstring("root")))

