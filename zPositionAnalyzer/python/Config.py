import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("Demo",Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_13_2_0_pre3/KFv0.1/samples/0_PU/zpos/n1/Eta_16/singlemuon_flatEGun_hgcalCenter/step1/step1_singlemuon_e100GeV_eta16_zpos_events1_nopu_1.root'
            'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/samples/0_PU/zpos/n10/Eta_17/singlemuon_flatEGun_hgcalCenter/step1/step1_singlemuon_e100GeV_eta17_zpos_events10_nopu_1.root'
    )
)

process.demo = cms.EDAnalyzer('zPositionAnalyzer')

process.p = cms.Path(process.demo)