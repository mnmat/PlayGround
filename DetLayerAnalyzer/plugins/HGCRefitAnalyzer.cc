// -*- C++ -*-
//
// Package:    PlayGround/HGCResidualAnalyzer
// Class:      HGCResidualAnalyzer
//
/**\class HGCResidualAnalyzer HGCResidualAnalyzer.cc PlayGround/HGCResidualAnalyzer/plugins/HGCResidualAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Tue, 12 Dec 2023 08:30:06 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
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

#include "TrajectoryToResiduals.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//ROOT includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include <TVector.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include <algorithm>
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveStats.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



using reco::TrackCollection;

class HGCResidualAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit HGCResidualAnalyzer(const edm::ParameterSet&);
  ~HGCResidualAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  void printTrackParameters(const reco::Track& track) const;
  void printTrackExtraParameters(const reco::Track& track) const;
  void mergeRecHitCollections(std::vector<const HGCRecHit*>& recHitCollection,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH) const;
  TransientTrackingRecHit::ConstRecHitContainer getTransientRecHits(
    const reco::TransientTrack& track) const; 

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
  
  int eventnr = 0;
  
  // Output
  TTree* residualTree;
  std::vector<int> res_evt;
  std::vector<int> res_layer;
  std::vector<int> res_detid;
  std::vector<float> res_diff_x;
  std::vector<float> res_diff_y;
  std::vector<float> res_pull_x;
  std::vector<float> res_pull_y;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
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
HGCResidualAnalyzer::HGCResidualAnalyzer(const edm::ParameterSet& iConfig)
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
      DetHGCTrackerToken_(esConsumes<DetHGCTracker, CaloGeometryRecord>()){

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed


  usesResource("TFileService");
  edm::Service<TFileService> file;

  residualTree = file->make<TTree>("Residuals","Residuals");

  residualTree->Branch("evt", &res_evt);
  residualTree->Branch("layer", &res_layer);
  residualTree->Branch("detid", &res_detid);
  residualTree->Branch("diff_x", &res_diff_x);
  residualTree->Branch("diff_y", &res_diff_y);
  residualTree->Branch("pull_x", &res_pull_x);
  residualTree->Branch("pull_y", &res_pull_y);
}

HGCResidualAnalyzer::~HGCResidualAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

void HGCResidualAnalyzer::printTrackParameters(const reco::Track& track) const{
  std::cout << "Chi2: " << track.chi2() << std::endl;
  std::cout << "isTimeOk: " << track.isTimeOk() << std::endl;
  std::cout << "ndof: " << track.ndof() << std::endl;
  std::cout << "normalizedChi2: " << track.normalizedChi2() << std::endl;
  std::cout << "charge: " << track.charge() << std::endl;
  std::cout << "qoverp: " << track.qoverp() << std::endl;
  std::cout << "theta: " << track.theta() << std::endl;
  std::cout << "lambda: " << track.lambda() << std::endl;
  std::cout << "dxy: " << track.dxy() << std::endl;
  std::cout << "d0: " << track.d0() << std::endl;
  std::cout << "dsz: " << track.dsz() << std::endl;
  std::cout << "dz: " << track.dz() << std::endl;
  std::cout << "p2: " << track.p2() << std::endl;
  std::cout << "p: " << track.p() << std::endl;
  std::cout << "pt2: " << track.pt2() << std::endl;
  std::cout << "pt: " << track.pt() << std::endl;
  std::cout << "px: " << track.px() << std::endl;
  std::cout << "py: " << track.py() << std::endl;
  std::cout << "pz: " << track.pz() << std::endl;
  std::cout << "phi: " << track.phi() << std::endl;
  std::cout << "eta: " << track.eta() << std::endl;
  std::cout << "vx: " << track.vx() << std::endl;
  std::cout << "vy: " << track.vy() << std::endl;
  std::cout << "vz: " << track.vz() << std::endl;
  std::cout << "momentum: " << track.momentum() << std::endl;
  std::cout << "referencePoint: " << track.referencePoint() << std::endl;
  std::cout << "t0: " << track.t0() << std::endl;
  std::cout << "beta: " << track.beta() << std::endl;
  std::cout << "vertex: " << track.vertex() << std::endl;
  std::cout << "parameters: " << track.parameters() << std::endl;
  std::cout << "covariance: " << track.covariance() << std::endl;
  // std::cout << "hitpattern: " << track.hitPattern() << std::endl;
  std::cout << "numberOfValidHits: " << track.numberOfValidHits() << std::endl;
  std::cout << "numberOfLostHits: " << track.numberOfLostHits() << std::endl;
  std::cout << "missingInnerHits: " << track.missingInnerHits() << std::endl;
  std::cout << "missingOuterHits: " << track.missingOuterHits() << std::endl;
  std::cout << "validFraction: " << track.validFraction() << std::endl;
  std::cout << "algoName: " << track.algoName() << std::endl;
  std::cout << "Found: " << track.found() << std::endl;
  std::cout << "Lost: " << track.lost() << std::endl;
}

void HGCResidualAnalyzer::printTrackExtraParameters(const reco::Track& track) const{
  std::cout << "OuterOk: " << track.outerOk() << std::endl;
  std::cout << "InnerOk: " << track.innerOk() << std::endl;
  std::cout << "InnerPos: " << track.innerPosition() << std::endl;
  std::cout << "OuterPos: " << track.outerPosition() << std::endl;
  std::cout << "InnerMom: " << track.innerMomentum() << std::endl;
  std::cout << "OuterMom: " << track.outerMomentum() << std::endl;    
  std::cout << "InnerStateCov: " << track.innerStateCovariance() << std::endl;
  std::cout << "OuterStateCov: " << track.outerStateCovariance() << std::endl;
  std::cout << "InnerDetId: " << track.innerDetId() << std::endl;
  std::cout << "OuterDetid: " << track.outerDetId() << std::endl;
  std::cout << "Size of RecHits: " << track.recHitsSize() << std::endl;
  std::cout << "OuterEta: " << track.outerEta() << std::endl;
  std::cout << "OuterTheta: " << track.outerTheta() << std::endl;  
  std::cout << "OuterPhi: " << track.outerPhi() << std::endl; 
  // Measurements
  for (auto& h : track.recHits()){
    std::cout << h->rawId() << std::endl;
  }
}

void HGCResidualAnalyzer::mergeRecHitCollections(std::vector<const HGCRecHit*>& recHitCollection,
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

TransientTrackingRecHit::ConstRecHitContainer HGCResidualAnalyzer::getTransientRecHits(
    const reco::TransientTrack& track) const {
  TransientTrackingRecHit::ConstRecHitContainer tkHits;
  //TransientTrackingRecHit::ConstRecHitContainer staHits;


  /*
  for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
    std::cout << (&**hit)->surface()->zSpan().first << std::endl;
    tkHits.push_back(std::shared_ptr<const TrackingRecHit>((&**hit)->clone()));
  }
  return tkHits;



  for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
    if ((*hit)->isValid()) {
      if ((*hit)->geographicalId().det() == DetId::Tracker && TrackerKeep((*hit)->geographicalId())) {
        tkHits.push_back(theTrackerRecHitBuilder->build(&**hit));
      } else if ((*hit)->geographicalId().det() == DetId::Muon && MuonKeep((*hit)->geographicalId())) {
        if ((*hit)->geographicalId().subdetId() == 3 && !theRPCInTheFit) {
          LogTrace("Reco|TrackingTools|TrackTransformer") << "RPC Rec Hit discarged";
          continue;
        }
        staHits.push_back(theMuonRecHitBuilder->build(&**hit));
      }
    }
  }

  if (staHits.empty())
    return staHits;

  copy(staHits.begin(), staHits.end(), back_inserter(tkHits));

  for (TransientTrackingRecHit::ConstRecHitContainer::const_iterator hit = tkHits.begin(); hit != tkHits.end(); ++hit) {
    DetId hitId = (*hit)->geographicalId();
    GlobalPoint glbpoint = trackingGeometry()->idToDet(hitId)->position();

    if (hitId.det() == DetId::Tracker) {
      if (hitId.subdetId() == StripSubdetector::TIB)
        LogTrace("TrackFitters") << glbpoint << " I am TIB " << tTopo_->tibLayer(hitId);
      else if (hitId.subdetId() == StripSubdetector::TOB)
        LogTrace("TrackFitters") << glbpoint << " I am TOB " << tTopo_->tobLayer(hitId);
      else if (hitId.subdetId() == StripSubdetector::TEC)
        LogTrace("TrackFitters") << glbpoint << " I am TEC " << tTopo_->tecWheel(hitId);
      else if (hitId.subdetId() == StripSubdetector::TID)
        LogTrace("TrackFitters") << glbpoint << " I am TID " << tTopo_->tidWheel(hitId);
      else if (hitId.subdetId() == (int)PixelSubdetector::PixelBarrel)
        LogTrace("TrackFitters") << glbpoint << " I am PixBar " << tTopo_->pxbLayer(hitId);
      else if (hitId.subdetId() == (int)PixelSubdetector::PixelEndcap)
        LogTrace("TrackFitters") << glbpoint << " I am PixFwd " << tTopo_->pxfDisk(hitId);
      else
        LogTrace("TrackFitters") << " UNKNOWN TRACKER HIT TYPE ";
    } else if (hitId.det() == DetId::Muon) {
      if (hitId.subdetId() == MuonSubdetId::DT)
        LogTrace("TrackFitters") << glbpoint << " I am DT " << DTWireId(hitId);
      else if (hitId.subdetId() == MuonSubdetId::CSC)
        LogTrace("TrackFitters") << glbpoint << " I am CSC " << CSCDetId(hitId);
      else if (hitId.subdetId() == MuonSubdetId::RPC)
        LogTrace("TrackFitters") << glbpoint << " I am RPC " << RPCDetId(hitId);
      else
        LogTrace("TrackFitters") << " UNKNOWN MUON HIT TYPE ";
    } else
      LogTrace("TrackFitters") << " UNKNOWN HIT TYPE ";
  }

  */

  return tkHits;
}

// ------------ method called for each event  ------------
void HGCResidualAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::cout << "Hello" << std::endl;

  using namespace edm;

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


  std::cout << "---------------------" << std::endl;
  std::cout << "Tracks" << std::endl;
  std::cout << "---------------------" << std::endl;

  auto seedingTrack = iEvent.get(seedToken_);
  for (const auto& track : iEvent.get(tracksToken_)) {

    // TrackParameters
    //printTrackParameters(track);
    // Requires TrackExtra
    //printTrackExtraParameters(track);

    // Build TransientTrack (link with magnetic field )
    std::cout << "Build Transient Track" << std::endl;
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



    std::cout << "Build Transient RecHits for Refit" << std::endl;
    TransientTrackingRecHit::ConstRecHitContainer recHitsforReFit;
    auto recHits = track.recHits();

    PropagationDirection propagationDirection = alongMomentum;
    int zside = +1;
    auto geomdets = hgcTracker_->getDiskLayers(zside);
    std::cout << "Start Loop" << std::endl;
    for (auto& hit : recHits){
      if (hit->rawId() == 0) continue;
      //auto hit = *it;
      std::cout << "Get Layer" << std::endl;
      int layer = rhtools_.getLayerWithOffset(hit->rawId());
      std::cout << "Get GeomDet" << std::endl;
      std::cout << geomdets.size() << ", " << hit->rawId() << ", " << layer-1 << std::endl;
      auto geomdet = geomdets[layer-1]->second;
      std::cout << "Build TransientTrackingrecHit" << std::endl;
      auto hgchit = HGCTransientTrackingRecHit::specificBuild(geomdet, hit);
      std::cout << "Pushback" << std::endl;
      recHitsforReFit.push_back(hgchit);
      std::cout << rhtools_.getPosition(hit->rawId()).x() << ", " << rhtools_.getPosition(hit->rawId()).y() << ", " << rhtools_.getPosition(hit->rawId()).z() << std::endl;
      std::cout << "Done with one iteration" << std::endl;
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

    if (trajectories.size() == 0){
      std::cout << "No trajectories found!!!" << std::endl;
      continue;
    }
    for (auto& tm : trajectories[0].measurements()){
      tm.layer();
      std::cout << tm.updatedState().globalPosition().x() << ", " << tm.updatedState().globalPosition().y() << ", " << tm.updatedState().globalPosition().z() << ", " << std::endl;
      std::cout << tm.predictedState().globalPosition().x() << ", " << tm.predictedState().globalPosition().y() << ", " << tm.predictedState().globalPosition().z() << ", " << std::endl;
      std::cout << tm.forwardPredictedState().globalPosition().x() << ", " << tm.forwardPredictedState().globalPosition().y() << ", " << tm.forwardPredictedState().globalPosition().z() << ", " << std::endl;
      //std::cout << tm.backwardPredictedState().globalPosition().x() << ", " << tm.backwardPredictedState().globalPosition().y() << ", " << tm.backwardPredictedState().globalPosition().z() << ", " << std::endl;

      std::cout << (tm.recHit()->geographicalId().rawId()) << std::endl;
    }

    std::cout << "Smoothe Trajectory" << std::endl;
    auto traj = smoother->trajectories(trajectories[0]);
    std::cout << "Size of trajectories: " << traj.size() << std::endl;
    for (auto& tm : traj[0].measurements()){
      tm.layer();
      std::cout << tm.updatedState().globalPosition().x() << ", " << tm.updatedState().globalPosition().y() << ", " << tm.updatedState().globalPosition().z() << ", " << std::endl;
      std::cout << tm.predictedState().globalPosition().x() << ", " << tm.predictedState().globalPosition().y() << ", " << tm.predictedState().globalPosition().z() << ", " << std::endl;
      std::cout << tm.forwardPredictedState().globalPosition().x() << ", " << tm.forwardPredictedState().globalPosition().y() << ", " << tm.forwardPredictedState().globalPosition().z() << ", " << std::endl;
      std::cout << tm.backwardPredictedState().globalPosition().x() << ", " << tm.backwardPredictedState().globalPosition().y() << ", " << tm.backwardPredictedState().globalPosition().z() << ", " << std::endl;

      std::cout << (tm.recHit()->geographicalId().rawId()) << std::endl;
    }


    // Residuals

    auto res = trajectoryToResiduals(traj[0]);
    auto measurements = traj[0].measurements();
    
    for (int i = 0; i<int(measurements.size());i++){

      auto rhit = measurements[i].recHit();
      
      int detid = rhit->rawId();
      int layer = rhtools_.getLayerWithOffset(detid);

      res_evt.push_back(eventnr);
      res_layer.push_back(layer);
      res_detid.push_back(detid);
      res_diff_x.push_back(res.residualX(i));
      res_diff_y.push_back(res.residualY(i));
      res_pull_x.push_back(res.pullX(i));
      res_pull_y.push_back(res.pullY(i));
    }

    residualTree->Fill();
    res_evt.clear();
    res_layer.clear();
    res_detid.clear();
    res_diff_x.clear();
    res_diff_y.clear();
    res_pull_x.clear();
    res_pull_y.clear();

    eventnr++;

    std::cout << "Done Fitting" << std::endl;
    //std::cout << "Number of Trajectories: " << trajectories.size() << std::endl;

    //for (auto& traj: trajectories){
    //  auto tms = traj.measurements();
    //  for (auto& tm: tms){
    //    auto tsos = tm.updatedState();
    //    std::cout << tsos.globalPosition().x() << ", " << tsos.globalPosition().y() << ", " << tsos.globalPosition().z() << std::endl;
    //  }
    //}
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void HGCResidualAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void HGCResidualAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCResidualAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCResidualAnalyzer);


