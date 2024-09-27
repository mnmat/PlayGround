import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("Demo",Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_13_2_0_pre3/KFv0.1/MCTruth/0_PU/zpos/n100/Eta_17/singlemuon_flatEGun_hgcalCenter/step3/step3_singlemuon_e10GeV_eta17_zpos_events100_nopu_1.root'
                                                       )                       )

outfile_ ="file:/afs/cern.ch/work/m/mmatthew/private/PatternRecognitionByKalmanFilter/CMSSW_13_2_0_pre3/src/Playground/PCaloHitAnalyzer/pcalohitAnalyzer.root"            

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfile_),
                                   closeFileFast = cms.untracked.bool(True))

from Playground.PCaloHitAnalyzer.PCaloHitAnalyzer_cfi import pcalohitAnalyzer
process.demo = pcalohitAnalyzer.clone()

process.p = cms.Path(process.demo)
