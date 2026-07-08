
import FWCore.ParameterSet.Config as cms
from Analyzers.Ntuplizer.Ntuplizer_cfi import ntuplizer

def customiseTICLForKFRefit(process):
    rescale = 5

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

    process.hgcdetlayergeometry = cms.ESProducer("HGCDetLayerGeometryESProducer")
    process.dummyGeometry = cms.ESProducer("DummyGeometryESProducer")

    # from SimG4Core.Application.g4SimHits_cfi import g4SimHits as _g4SimHits

    # process.load("TrackPropagation.Geant4e.Geant4ePropagator_cfi")
    # process.geopro = cms.EDProducer("GeometryProducer",
    #     GeoFromDD4hep = cms.bool(False),
    #     UseMagneticField = cms.bool(True),
    #     UseSensitiveDetectors = cms.bool(False),
    #     MagneticField =  _g4SimHits.MagneticField.clone()
    # )

    # from Configuration.ProcessModifiers.dd4hep_cff import dd4hep
    # dd4hep.toModify(process.geopro, GeoFromDD4hep = True )

    process.kfrefitter = cms.ESProducer( "KFTrajectoryFitterESProducer",
        appendToDataLabel = cms.string( "" ),
        ComponentName = cms.string( "testReFitter" ),
        RecoGeometry = cms.string( "DummyGeometry" ),
        #Propagator = cms.string("Geant4ePropagator")
    )

    process.kfsmoother = cms.ESProducer("KFTrajectorySmootherESProducer",
        ComponentName = cms.string("testSmoother"),
        RecoGeometry = cms.string("DummyGeometry"),
        errorRescaling = cms.double(rescale),
        #Propagator = cms.string("Geant4ePropagator")
    )

    #process.testTrackAnalyzer = cms.EDAnalyzer('TestTrackAnalyzer',
    #    tracks = cms.untracked.InputTag('generalTracks'))

    process.refitter = cms.EDAnalyzer('HGCRefitAnalyzer',
        seeds =  cms.untracked.InputTag('generalTracks'),
        tracks    = cms.untracked.InputTag('ticlTrackstersKalmanFilter','HGCALTracks'),
        tracksters = cms.untracked.InputTag('ticlTrackstersKalmanFilter'),
        HGCEEInput =  cms.InputTag("HGCalRecHit", "HGCEERecHits"),
        HGCFHInput =  cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
        HGCBHInput =  cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
        Fitter = cms.ESInputTag('','testReFitter'),
        Smoother = cms.ESInputTag('','testSmoother'),
    )

    #process.p = cms.Path(process.dummyGeometryAnalyzer*process.testTrackAnalyzer)
    #process.p = cms.Path(process.dummyGeometryAnalyzer)
    #process.p = cms.Path(process.refitter)

    process.kfhitproducer = cms.EDProducer("KFHitProducer",
        seeds =  cms.untracked.InputTag('generalTracks'),
        tracks    = cms.untracked.InputTag('ticlTrackstersKalmanFilter','HGCALTracks'),
        tracksters = cms.untracked.InputTag('ticlTrackstersKalmanFilter'),
        HGCEEInput =  cms.InputTag("HGCalRecHit", "HGCEERecHits"),
        HGCFHInput =  cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
        HGCBHInput =  cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
        Fitter = cms.ESInputTag('','testReFitter'),
        Smoother = cms.ESInputTag('','testSmoother'),
        KFHits = cms.InputTag("ticlTrackstersKalmanFilter","KFHits","RECO"),
    )

    process.prophitproducer = cms.EDProducer("KFHitProducer",
        seeds =  cms.untracked.InputTag('generalTracks'),
        tracks    = cms.untracked.InputTag('ticlTrackstersStandalonePropagator','HGCALTracks'),
        tracksters = cms.untracked.InputTag('ticlTrackstersKalmanFilter'),
        HGCEEInput =  cms.InputTag("HGCalRecHit", "HGCEERecHits"),
        HGCFHInput =  cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
        HGCBHInput =  cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
        Fitter = cms.ESInputTag('','testReFitter'),
        Smoother = cms.ESInputTag('','testSmoother'),
        KFHits = cms.InputTag("ticlTrackstersStandalonePropagator","KFHits","RECO"),
    )

    process.kfhitanalyzer = cms.EDAnalyzer("KFHitAnalyzer",
        kfhits = cms.untracked.InputTag('kfhitproducer','RefitUpdated')
    )

    process.ntuplizerRefitUpdated = ntuplizer.clone(
        KFHits = cms.InputTag("kfhitproducer","RefitUpdated"),
        PropHits = cms.InputTag("prophitproducer","RefitUpdated"),
        associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation"),
        tracks    = cms.untracked.InputTag('generalTracks'),
    )

    process.ntuplizerRefitForwardPredicted = ntuplizer.clone(
        KFHits = cms.InputTag("kfhitproducer","RefitForwardPredicted"),
        PropHits = cms.InputTag("prophitproducer","RefitForwardPredicted"),
        associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation"),
        tracks    = cms.untracked.InputTag('generalTracks'),
    )

    process.ntuplizerSmoothedUpdated = ntuplizer.clone(
        KFHits = cms.InputTag("kfhitproducer","SmoothedUpdated"),
        PropHits = cms.InputTag("prophitproducer","SmoothedUpdated"),
        associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation"),
        tracks    = cms.untracked.InputTag('generalTracks'),
    )

    process.ntuplizerSmoothedForwardPredicted = ntuplizer.clone(
        KFHits = cms.InputTag("kfhitproducer","SmoothedForwardPredicted"),
        PropHits = cms.InputTag("prophitproducer","SmoothedForwardPredicted"),
        associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation"),
        tracks    = cms.untracked.InputTag('generalTracks'),
    )

    process.ntuplizerSmoothedBackwardPredicted = ntuplizer.clone(
        KFHits = cms.InputTag("kfhitproducer","SmoothedBackwardPredicted"),
        PropHits = cms.InputTag("prophitproducer","SmoothedBackwardPredicted"),
        associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation"),
        tracks    = cms.untracked.InputTag('generalTracks'),
    )

    process.ntuplizer_refit = cms.Sequence(process.ntuplizerRefitForwardPredicted+process.ntuplizerRefitUpdated+process.ntuplizerSmoothedUpdated+process.ntuplizerSmoothedForwardPredicted+process.ntuplizerSmoothedBackwardPredicted )

    from Analyzers.EfficiencyAnalyzerDemo.EfficiencyAnalyzerDemo_cfi import efficiencyAnalyzerDemo

    process.efficiencyAnalyzerRefitUpdated = efficiencyAnalyzerDemo.clone(
        KFHits = cms.InputTag("kfhitproducer","RefitUpdated"),
        tracks    = cms.untracked.InputTag('generalTracks'),
        associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation")
    )

    process.efficiencyAnalyzerRefitForwardPredicted = efficiencyAnalyzerDemo.clone(
        KFHits = cms.InputTag("kfhitproducer","RefitForwardPredicted"),
        tracks    = cms.untracked.InputTag('generalTracks'),
        associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation")
    )

    process.efficiencyAnalyzerSmoothedForwardPredicted = efficiencyAnalyzerDemo.clone(
        KFHits = cms.InputTag("kfhitproducer","SmoothedForwardPredicted"),
        tracks    = cms.untracked.InputTag('generalTracks'),
        associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation")
    )

    process.efficiencyAnalyzerSmoothedBackwardPredicted = efficiencyAnalyzerDemo.clone(
        KFHits = cms.InputTag("kfhitproducer","SmoothedBackwardPredicted"),
        tracks    = cms.untracked.InputTag('generalTracks'),
        associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation")
    )

    process.efficiencyAnalyzerSmoothedUpdated = efficiencyAnalyzerDemo.clone(
        KFHits = cms.InputTag("kfhitproducer","SmoothedUpdated"),
        tracks    = cms.untracked.InputTag('generalTracks'),
        associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation")
    )

    process.efficiencyAnalyzer_refit = cms.Sequence(process.efficiencyAnalyzerRefitForwardPredicted+process.efficiencyAnalyzerRefitUpdated+process.efficiencyAnalyzerSmoothedUpdated+process.efficiencyAnalyzerSmoothedForwardPredicted+process.efficiencyAnalyzerSmoothedBackwardPredicted )

    from PlayGround.PCaloHitAnalyzer.customiseTICLForPCaloHitAnalyzer import customiseTICLForPCaloHitAnalyzer
    process = customiseTICLForPCaloHitAnalyzer(process)

    process.residuals = cms.EDAnalyzer('HGCResidualAnalyzer',
        seeds =  cms.untracked.InputTag('generalTracks'),
        tracks    = cms.untracked.InputTag('ticlTrackstersKalmanFilter','HGCALTracks'),
        tracksters = cms.untracked.InputTag('ticlTrackstersKalmanFilter'),
        HGCEEInput =  cms.InputTag("HGCalRecHit", "HGCEERecHits"),
        HGCFHInput =  cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
        HGCBHInput =  cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
        Fitter = cms.ESInputTag('','testReFitter'),
        Smoother = cms.ESInputTag('','testSmoother')
    )

    process.refit_path = cms.Sequence(process.kfhitproducer*process.prophitproducer*process.ntuplizer_refit*process.efficiencyAnalyzer_refit)
    #process.add_(cms.Service("AdaptorConfig", native=cms.untracked.vstring("root")))

    return process
