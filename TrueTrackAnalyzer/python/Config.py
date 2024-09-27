import FWCore.ParameterSet.Config as cms
import sys

# Settings
# eta = sys.argv[2]
# energy = sys.argv[3]
# nevents = sys.argv[4]
# idx = sys.argv[5]
# input_dir = sys.argv[6]
# eta = eta.replace(".","")

eta = "23"
energy = "100"
nevents = "100"
idx = "1"
mb = "mb_ngun"
cap = "zpos"
input_dir = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/debugging/TrackAnalysis/200_PU/"
output_dir = "/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/debugging/TrackAnalysis/200_PU/"


infile_ = "file:"+ input_dir + "/" + cap + "/n" + nevents +  "/Eta_" + eta + "/singlemuon_flatEGun_hgcalCenter/step3/step3_singlemuon_e" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"
outfile_ = output_dir + "/" + cap + "/n" + nevents +  "/Eta_" + eta + "/singlemuon_flatEGun_hgcalCenter/step3/trackAnalyzer" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"
#outfile_ ="file:" + output_dir + "ntuplizer_singlemuon_e" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"            


mb = "mb_ngun"
cap = "zpos"

# Define Process
process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:'+ input_dir + "/" + cap + "/n" + nevents +  "/Eta_" + eta + "/singlemuon_flatEGun_hgcalCenter/step3/step3_singlemuon_e" + energy + "GeV_eta" + eta +"_" + cap +"_events" + nevents + "_nopu_" +idx +".root"))



process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfile_),
                                   closeFileFast = cms.untracked.bool(True)
                               )

process.load("Validation.RecoMuon.associators_cff")
#from Validation.RecoMuon.associators_cff import *

process.demo = cms.EDAnalyzer('TrueTrackAnalyzer',
   tracks    = cms.untracked.InputTag('generalTracks'),
   muons = cms.untracked.InputTag('globalMuons'),
   standalone = cms.untracked.InputTag('standAloneMuons'),
   associators = cms.untracked.VInputTag(["trackingParticleRecoTrackAsssociation","tpToGlbMuonAssociation","tpToStaMuonAssociation"]),
   recoMuons = cms.untracked.InputTag('muons'),
   eta = cms.string(eta),
   energy = cms.string(energy), 
   muonassociators = cms.untracked.VInputTag("")
                              )

process.p = cms.Path(process.demo)
process.add_(cms.Service("AdaptorConfig", native=cms.untracked.vstring("root")))
