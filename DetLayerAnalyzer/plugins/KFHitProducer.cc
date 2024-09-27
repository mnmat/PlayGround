// -*- C++ -*-
//
// Package:    PlayGround/KFHitProducer
// Class:      KFHitProducer
//
/**\class KFHitProducer KFHitProducer.cc PlayGround/KFHitProducer/plugins/KFHitProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Mon, 05 Aug 2024 15:18:18 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HGCalReco/interface/HGCTrackingRecHit.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformerForGlobalCosmicMuons.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformerForCosmicMuons.h"
#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"


#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "PlayGround/DetLayerAnalyzer/plugins/HGCTransientTrackingRecHit.h"
#include "PlayGround/DetLayerAnalyzer/plugins/DetHGCTracker.h"
#include "DataFormats/HGCalReco/interface/KFHit.h"

//
// class declaration
//

using reco::TrackCollection;

class KFHitProducer : public edm::stream::EDProducer<> {
public:
  explicit KFHitProducer(const edm::ParameterSet&);
  ~KFHitProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;
  void mergeRecHitCollections(std::vector<const HGCRecHit*>& recHitCollection,
                              const HGCRecHitCollection& rechitsEE,
                              const HGCRecHitCollection& rechitsFH,
                              const HGCRecHitCollection& rechitsBH) const;


  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------

  edm::EDGetTokenT<TrackCollection> seedToken_;
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::ESGetToken<TrajectoryFitter, TrajectoryFitter::Record> theFitterToken_;
  edm::ESGetToken<TrajectorySmoother, TrajectoryFitter::Record> theSmootherToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfieldtoken_;
  edm::ESGetToken<DetHGCTracker,CaloGeometryRecord> DetHGCTrackerToken_;

  std::vector<const HGCRecHit*> recHitCollection;
  hgcal::RecHitTools rhtools_;
  const DetHGCTracker* hgcTracker_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
KFHitProducer::KFHitProducer(const edm::ParameterSet& iConfig)
  : seedToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("seeds"))),
    tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
    tracksterToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getUntrackedParameter<edm::InputTag>("tracksters"))),
    hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCEEInput"))),
    hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCFHInput"))),
    hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCBHInput"))),
    theFitterToken_(esConsumes<TrajectoryFitter,TrajectoryFitter::Record>(iConfig.getParameter<edm::ESInputTag>("Fitter"))),
    theSmootherToken_(esConsumes<TrajectorySmoother,TrajectoryFitter::Record>(iConfig.getParameter<edm::ESInputTag>("Smoother"))),
    caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
    bfieldtoken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
    DetHGCTrackerToken_(esConsumes<DetHGCTracker, CaloGeometryRecord>()) {

  // produces<std::vector<KFHit>>("Refit_updated");
  // produces<std::vector<KFHit>>("Refit_forward_predicted");

  // produces<std::vector<KFHit>>("Smoothed_updated");
  // produces<std::vector<KFHit>>("Smoothed_forward_predicted");
  // produces<std::vector<KFHit>>("Refit_backward_predicted"); 

  produces<std::vector<KFHit>>("RefitUpdated");
  produces<std::vector<KFHit>>("RefitForwardPredicted");

  produces<std::vector<KFHit>>("SmoothedUpdated");
  produces<std::vector<KFHit>>("SmoothedForwardPredicted");
  produces<std::vector<KFHit>>("SmoothedBackwardPredicted"); 

  //register your products
  /* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
  */
  //now do what ever other initialization is needed
}

KFHitProducer::~KFHitProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

void KFHitProducer::mergeRecHitCollections(std::vector<const HGCRecHit*>& recHitCollection,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH) const {
  
  recHitCollection.clear();
  for (const auto& hit : rechitsEE) {
    recHitCollection.push_back(&hit);
  }

  for (const auto& hit : rechitsFH) {
    recHitCollection.push_back(&hit);
  }

  for (const auto& hit : rechitsBH) {
    recHitCollection.push_back(&hit);
  }
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void KFHitProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  auto result_refit_updated = std::make_unique<std::vector<KFHit>>();
  auto result_refit_forward_predicted = std::make_unique<std::vector<KFHit>>();

  auto result_smoothed_updated = std::make_unique<std::vector<KFHit>>();
  auto result_smoothed_forward_predicted = std::make_unique<std::vector<KFHit>>();
  auto result_smoothed_backward_predicted = std::make_unique<std::vector<KFHit>>();

  edm::Handle<HGCRecHitCollection> ee_hits;
  edm::Handle<HGCRecHitCollection> fh_hits;
  edm::Handle<HGCRecHitCollection> bh_hits;

  iEvent.getByToken(hgcalRecHitsEEToken_, ee_hits);
  iEvent.getByToken(hgcalRecHitsFHToken_, fh_hits);
  iEvent.getByToken(hgcalRecHitsBHToken_, bh_hits);
  mergeRecHitCollections(recHitCollection, *ee_hits, *fh_hits, *bh_hits);

  edm::ESHandle<TrajectoryFitter> theFitter;
  theFitter = iSetup.getHandle(theFitterToken_);
  auto fitter = theFitter.product();

  edm::ESHandle<TrajectorySmoother> theSmoother;
  theSmoother = iSetup.getHandle(theSmootherToken_);
  auto smoother = theSmoother.product();

  const CaloGeometry* geom = &iSetup.getData(caloGeomToken_);
  rhtools_.setGeometry(*geom);

  hgcTracker_ = &iSetup.getData(DetHGCTrackerToken_);

  edm::ESHandle<MagneticField> bfield_;
  bfield_ = iSetup.getHandle(bfieldtoken_);
  auto magneticField = bfield_.product();

  // ----------
  // Tracks
  // ----------

  auto seedingTrack = iEvent.get(seedToken_);
  for (const auto& track : iEvent.get(tracksToken_)) {

    // TrackParameters
    //printTrackParameters(track);
    // Requires TrackExtra
    //printTrackExtraParameters(track);

    // Build TransientTrack (link with magnetic field )
    reco::TransientTrack ttrack(track, magneticField);

    // Build transient RecHits
    //std::cout << "Build Transient RecHits" << std::endl;
    //TransientTrackingRecHit::ConstRecHitContainer recHitsforReFit =  getTransientRecHits(ttrack);

      /*
      auto recHits = track.recHits();
      PropagationDirection propagationDirection = alongMomentum;
     
      auto hit = (*recHits.begin());
      const HGCDiskLayer* layerdisk = hgcTracker_->firstDisk(zside,propagationDirection);
      auto geomdet = layerdisk->second;
      auto hgchit = HGCTransientTrackingRecHit::specificBuild(geomdet, hit);
      */



    TransientTrackingRecHit::ConstRecHitContainer recHitsforReFit;
    auto recHits = track.recHits();

    PropagationDirection propagationDirection = alongMomentum;
    int zside = +1;
    auto geomdets = hgcTracker_->getDiskLayers(zside);
    for (auto& hit : recHits){
      if (hit->rawId() == 0) continue;
      //auto hit = *it;
      int layer = rhtools_.getLayerWithOffset(hit->rawId());
      std::cout << geomdets.size() << ", " << hit->rawId() << ", " << layer-1 << std::endl;
      auto geomdet = geomdets[layer-1]->second;
      auto hgchit = HGCTransientTrackingRecHit::specificBuild(geomdet, hit);
      recHitsforReFit.push_back(hgchit);
      std::cout << rhtools_.getPosition(hit->rawId()).x() << ", " << rhtools_.getPosition(hit->rawId()).y() << ", " << rhtools_.getPosition(hit->rawId()).z() << std::endl;
    }

    // Build first TSOS
    //TSOS firstTSOS = getFirstTSOS();
    //

    FreeTrajectoryState fts = trajectoryStateTransform::innerFreeState(track,magneticField);
    auto surface = hgcTracker_->firstDisk(zside,propagationDirection)->second->surface();
    TrajectoryStateOnSurface* firstTSOS = new TrajectoryStateOnSurface(fts,surface);
    //firstTSOS = TrajectoryStateOnSurface();

    /*
    std::cout << "Create FreeTrajectoryState" << std::endl;
    [[maybe_unused]] FreeTrajectoryState fts = trajectoryStateTransform::outerFreeState(seedingTrack[0],magneticField);
    std::cout << "Done creating FreeTrajectoryState" << std::endl;
    auto surface = hgcTracker_->firstDisk(zside,propagationDirection)->second->surface();
    std::cout << "Make first TSOS" << std::endl;
    TrajectoryStateOnSurface* firstTSOS = new TrajectoryStateOnSurface(fts,surface);
    std::cout << "done making first TSOS" << std::endl;    

    */

    std::cout << firstTSOS->globalPosition().x() << ", " << firstTSOS->globalPosition().y() << ", " << firstTSOS->globalPosition().z() << std::endl;

    // Make TrajectorySeeds
    TrajectorySeed seed({}, {}, propagationDirection);

    // Can I create a trajectoryMeasurement from a DetLayer?

    //auto& pred = firstTSOS;
    //TransientTrackingRecHit* ihit = recHitsforReFit.back();
    //HGCDetLayer detlayer = HGCDetLayer();
    //TrajectoryMeasurement TM = TrajectoryMeasurement(pred, ihit, 0, detlayer);

    // Can I create a DetLayer from HGCDiskGeomDet

    // Fit Trajectories
    std::vector<Trajectory> trajectories = fitter->fit(seed, recHitsforReFit, *firstTSOS);
    std::cout << "Size of trajectories: " << trajectories.size() << std::endl;

    if (trajectories.size() == 0){
      std::cout << "No trajectory created!!!!" << std::endl;
      continue;
    } 
    int trackId = 0;
    for (auto& tm : trajectories[0].measurements()){
      auto updated_tsos = tm.updatedState();
      auto forward_predicted_tsos = tm.forwardPredictedState();
      auto detid = tm.recHit()->geographicalId();
      int layer = rhtools_.getLayerWithOffset(detid.rawId());

      KFHit* refit_updated_hit = new KFHit(updated_tsos,detid,track,trackId,layer);
      KFHit* refit_forward_predicted_hit = new KFHit(forward_predicted_tsos,detid,track,trackId,layer);

      (*result_refit_updated).push_back(*refit_updated_hit);
      (*result_refit_forward_predicted).push_back(*refit_forward_predicted_hit);

      //std::cout << tm.updatedState().globalPosition().x() << ", " << tm.updatedState().globalPosition().y() << ", " << tm.updatedState().globalPosition().z() << ", " << std::endl;
      //std::cout << tm.predictedState().globalPosition().x() << ", " << tm.predictedState().globalPosition().y() << ", " << tm.predictedState().globalPosition().z() << ", " << std::endl;
      //std::cout << tm.forwardPredictedState().globalPosition().x() << ", " << tm.forwardPredictedState().globalPosition().y() << ", " << tm.forwardPredictedState().globalPosition().z() << ", " << std::endl;
      //std::cout << tm.backwardPredictedState().globalPosition().x() << ", " << tm.backwardPredictedState().globalPosition().y() << ", " << tm.backwardPredictedState().globalPosition().z() << ", " << std::endl;
      //std::cout << (tm.recHit()->geographicalId().rawId()) << std::endl;
    }

    std::cout << "Smoothe Trajectory" << std::endl;
    auto traj = smoother->trajectories(trajectories[0]);
    std::cout << "Size of trajectories: " << traj.size() << std::endl;
    for (auto& tm : traj[0].measurements()){
      auto updated_tsos = tm.updatedState();
      auto forward_predicted_tsos = tm.forwardPredictedState();
      auto backward_predicted_tsos = tm.backwardPredictedState();
      auto detid = tm.recHit()->geographicalId();
      int layer = rhtools_.getLayerWithOffset(detid.rawId());

      KFHit* smoothed_updated_hit = new KFHit(updated_tsos,detid,track,trackId,layer);
      KFHit* smoothed_forward_predicted_hit = new KFHit(forward_predicted_tsos,detid,track,trackId,layer);
      KFHit* smoothed_backward_predicted_hit = new KFHit(backward_predicted_tsos,detid,track,trackId,layer);

      (*result_smoothed_updated).push_back(*smoothed_updated_hit);
      (*result_smoothed_forward_predicted).push_back(*smoothed_forward_predicted_hit);
      (*result_smoothed_backward_predicted).push_back(*smoothed_backward_predicted_hit);

      // std::cout << tm.updatedState().globalPosition().x() << ", " << tm.updatedState().globalPosition().y() << ", " << tm.updatedState().globalPosition().z() << ", " << std::endl;
      // std::cout << tm.predictedState().globalPosition().x() << ", " << tm.predictedState().globalPosition().y() << ", " << tm.predictedState().globalPosition().z() << ", " << std::endl;
      // std::cout << tm.forwardPredictedState().globalPosition().x() << ", " << tm.forwardPredictedState().globalPosition().y() << ", " << tm.forwardPredictedState().globalPosition().z() << ", " << std::endl;
      // std::cout << tm.backwardPredictedState().globalPosition().x() << ", " << tm.backwardPredictedState().globalPosition().y() << ", " << tm.backwardPredictedState().globalPosition().z() << ", " << std::endl;
    }
  }   
  iEvent.put(std::move(result_refit_updated),"RefitUpdated");
  iEvent.put(std::move(result_refit_forward_predicted),"RefitForwardPredicted");
  iEvent.put(std::move(result_smoothed_updated),"SmoothedUpdated");
  iEvent.put(std::move(result_smoothed_forward_predicted),"SmoothedForwardPredicted");
  iEvent.put(std::move(result_smoothed_backward_predicted),"SmoothedBackwardPredicted");

  // iEvent.put(std::move(result_refit_updated));
  // iEvent.put(std::move(result_refit_forward_predicted));
  // iEvent.put(std::move(result_smoothed_updated));
  // iEvent.put(std::move(result_smoothed_forward_predicted));
  // iEvent.put(std::move(result_smoothed_backward_predicted));
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void KFHitProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void KFHitProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
KFHitProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
KFHitProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
KFHitProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
KFHitProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void KFHitProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(KFHitProducer);
