import FWCore.ParameterSet.Config as cms
import argparse

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process('PROD',Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("TrackingTools.TrackRefitter.TracksToTrajectories_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools.TrackFitters.TrackFitters_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T25', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Specify Eta

import sys
import os

eta = str(17)
energy = str(100)
nevents = str(10)
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

process.hgctracker = cms.ESProducer("DetHGCTrackerESProducer",
    radlen = cms.vdouble(1.53770787, 0.71064359, 1.45345887, 0.57315113, 1.02882455,
        0.92384098, 1.15461784, 0.72404336, 0.9948446 , 1.0107427 ,
        1.05947235, 0.91730167, 1.20028302, 0.6703572 , 0.98144224,
        1.01024202, 1.08452792, 0.86282149, 1.53923452, 0.99185102,
        1.67874405, 0.70709974, 1.63824099, 0.97162878, 1.74571227,
        0.69011827, 2.92834302, 3.01147101, 3.0583451 , 3.12601533,
        2.85205937, 2.95217992, 3.14263578, 3.07471756, 3.05502943,
        2.82345623, 3.0230636 , 4.29398744, 3.9234094 , 4.27748842,
        3.91229994, 4.23728221, 4.02845205, 4.21537293, 4.32452121,
        3.83363941, 4.32332509),
    xi = cms.vdouble(
        0.00264665, 0.00050171, 0.00081145, 0.0003883 , 0.00049233,
        0.00066116, 0.00059059, 0.00050953, 0.0004874 , 0.00069975,
        0.00051153, 0.00065396, 0.00062511, 0.00046753, 0.00048085,
        0.00070629, 0.00052536, 0.00061608, 0.00068004, 0.0006793 ,
        0.00079585, 0.0005112 , 0.00073527, 0.00067903, 0.00081588,
        0.00049695, 0.00292419, 0.0029901 , 0.00304991, 0.00309913,
        0.00283049, 0.00295079, 0.00310096, 0.00307401, 0.00304154,
        0.00278516, 0.00302603, 0.00425792, 0.00389764, 0.00425274,
        0.0038517 , 0.00426734, 0.00399329, 0.00420363, 0.00429176,
        0.00378711, 0.004295657)
    )


process.refitter = cms.EDAnalyzer('HGCDetLayerGeometryAnalyzer',
   seeds =  cms.untracked.InputTag('generalTracks'),
   tracks    = cms.untracked.InputTag('ticlTrackstersKalmanFilter','HGCALTracks'),
)

process.p = cms.Path(process.refitter)

process.add_(cms.Service("AdaptorConfig", native=cms.untracked.vstring("root")))