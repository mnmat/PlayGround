// -*- C++ -*-
//
// Package:    Analyzer/IsolationRequirementAnalyzer
// Class:      IsolationRequirementAnalyzer
//
/**\class IsolationRequirementAnalyzer IsolationRequirementAnalyzer.cc Analyzer/IsolationRequirementAnalyzer/plugins/IsolationRequirementAnalyzer.cc

 Description: [one line class summary]

 Implementation:f
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Thu, 30 Jun 2022 10:49:35 GMT
//
//

// system include files
#include <memory>
#include <numeric>
#include <sstream>
#include <any>
#include <iomanip>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <fstream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TProfile.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/KFHit.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/Math/interface/Vector3D.h"

#include "Validation/HGCalValidation/interface/CaloParticleSelector.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelectorBase.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelectorBase.h"

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
using namespace ticl;


enum class Dtype{
  si120,
  si200,
  si300,
  sc,
  undefined
};


class IsolationRequirementAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:

  explicit IsolationRequirementAnalyzer(const edm::ParameterSet&);
  ~IsolationRequirementAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  bool inSensorCell(DetId detid_, GlobalPoint point);
  void clear_arrays();
  virtual void fillHitMap(std::map<DetId, std::pair<const HGCRecHit*,float>>& hitMap, 
      const HGCRecHitCollection& rechitsEE, 
      const HGCRecHitCollection& rechitsFH,
      const HGCRecHitCollection& rechitsBH,
      const std::vector<float>& lcmask) const;
  std::vector<int> matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_);
  LocalError calculateLocalError(DetId detid_, const HGCalDDDConstants* ddd);
  Dtype getDtype(const DetId& detid_) const;
  std::string writeDtype(Dtype& dtype) const;


  hgcal::RecHitTools recHitTools_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::EDGetTokenT<std::vector<KFHit>> KFHitsToken_;
  edm::EDGetTokenT<std::vector<KFHit>> PropHitsToken_;

  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
  edm::EDGetTokenT<std::vector<float>> lcMaskToken_;
  std::vector<edm::InputTag> associators;
  std::vector<edm::EDGetTokenT<reco::SimToRecoCollection>> associatormapStRs;
  std::vector<edm::EDGetTokenT<reco::RecoToSimCollection>> associatormapRtSs;
  edm::EDGetTokenT<std::vector<SimVertex>> simVerticesToken_;

  std::vector<std::string> detectors, objects, positions, hittypes;


  // variables
  
  int eventnr =0;
  std::shared_ptr<hgcal::RecHitTools> recHitTools;

  TTree* signalTree;
  TTree* neighborTree;

  // KF

  std::vector<float> sig_x;
  std::vector<float> sig_y;
  std::vector<float> sig_z;
  std::vector<float> sig_e;
  std::vector<int> sig_detid;
  std::vector<int> sig_layer;
  std::vector<int> sig_evt;
  std::vector<int> sig_inLC;
  std::vector<int> sig_isolated;
  std::vector<int> sig_totalNeighbors;
  std::vector<int> sig_hitNeighbors;
  std::vector<float> sig_efrac;
  std::vector<float> sig_time;

  std::vector<float> nbor_x;
  std::vector<float> nbor_y;
  std::vector<float> nbor_z;
  std::vector<float> nbor_e;
  std::vector<int> nbor_detid;
  std::vector<int> nbor_layer;
  std::vector<int> nbor_evt;
  std::vector<float> nbor_time;
  std::vector<std::string> nbor_dtype;

  //std::vector<float> sim_mask;


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
IsolationRequirementAnalyzer::IsolationRequirementAnalyzer(const edm::ParameterSet& iConfig) :
      caloParticlesToken_(consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticles"))), 
      hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
      hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
      hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
      //abs_failToken_(consumes<float>(iConfig.getParameter<edm::InputTag>("abs_fail"))),
      KFHitsToken_(consumes<std::vector<KFHit>>(iConfig.getParameter<edm::InputTag>("KFHits"))),
      hgcalLayerClustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("hgcalLayerClusters"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      tracksToken_(consumes<edm::View<reco::Track>>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
      lcMaskToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("lcMask"))),
      associators(iConfig.getUntrackedParameter<std::vector<edm::InputTag>>("associators")),
      simVerticesToken_(consumes<std::vector<SimVertex>>(iConfig.getParameter<edm::InputTag>("simVertices"))){


  associatormapRtSs.push_back(consumes<reco::RecoToSimCollection>(associators[0]));

  detectors = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};
  hittypes = {"Simhits","Rechits","KF"};
  objects = {"Simhits", "Rechits"};        
  positions = {"KF","Prop"};
  //recHitTools_.reset(new hgcal::RecHitTools());
  //now do what ever initialization is needed

  usesResource("TFileService");
  edm::Service<TFileService> file;
  signalTree = file->make<TTree>("KFHit","KFHit");
  neighborTree = file->make<TTree>("Neighbors","Neighbors");


  // KF

  signalTree->Branch("x", &sig_x);
  signalTree->Branch("y", &sig_y);
  signalTree->Branch("z", &sig_z);
  signalTree->Branch("e", &sig_e);
  signalTree->Branch("layer", &sig_layer);
  signalTree->Branch("evt", &sig_evt);
  signalTree->Branch("detid", &sig_detid);
  signalTree->Branch("inLC", &sig_inLC);
  signalTree->Branch("isolated", &sig_isolated);
  signalTree->Branch("totalNeighbors", &sig_totalNeighbors);
  signalTree->Branch("hitNeighbors", &sig_hitNeighbors);
  signalTree->Branch("efrac", &sig_efrac);
  signalTree->Branch("time", &sig_time);


  neighborTree->Branch("x", &nbor_x);
  neighborTree->Branch("y", &nbor_y);
  neighborTree->Branch("z", &nbor_z);
  neighborTree->Branch("e", &nbor_e);
  neighborTree->Branch("layer", &nbor_layer);
  neighborTree->Branch("evt", &nbor_evt);
  neighborTree->Branch("detid", &nbor_detid);
  neighborTree->Branch("time", &nbor_time);
  neighborTree->Branch("dtype", &nbor_dtype);

    std::vector<float> nbor_time;
  std::vector<std::string> nbor_dtype;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

IsolationRequirementAnalyzer::~IsolationRequirementAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------

LocalError IsolationRequirementAnalyzer::calculateLocalError(DetId id, const HGCalDDDConstants* ddd){
  if(recHitTools_.isSilicon(id)){
    float A;
    if(recHitTools_.getSiThickness(id) < 200) A = 1.18; // TODO: replace with non-hardcoded value; hardcoded value from TDR
    else  A = 0.52; // TODO: replace with non-hardcoded value; hardcoded value from TDR
    float a = sqrt(2*A/(3*sqrt(3)));
    double varx = pow(a,4)*5*sqrt(3)/(16*A); // x
    double vary = pow(a,4)*5*sqrt(3)/(16*A); // y 
    return LocalError(varx, 0, vary);
  }
  else{
    const GlobalPoint &pos = recHitTools_.getPosition(id);
    double r = sqrt(pos.x()*pos.x() + pos.y()*pos.y());
    auto radiusLayer = ddd->getRadiusLayer(recHitTools_.getLayer(id));
    int idx = static_cast<int>(std::lower_bound(radiusLayer.begin(), radiusLayer.end(),r)-radiusLayer.begin());
    float rmax = radiusLayer[idx];
    float rmin = radiusLayer[idx-1];

    double phi = recHitTools_.getPhi(id) + M_PI; // radians [0, 2pi]
    double dphi = recHitTools_.getScintDEtaDPhi(id).second; // radians
    double phimin = phi - 0.5*dphi;
    double phimax = phi + 0.5*dphi;

    double A = (rmax*rmax - rmin*rmin)*M_PI*dphi/(2*M_PI);

    double ex2 = 1/(8*A) * (pow(rmax,4) - pow(rmin,4)) * (-phimin - sin(phimin)*cos(phimin) + phimax + sin(phimax)*cos(phimax));
    double ex = 1/(3*A) * (pow(rmax,3) - pow(rmin,3)) * (sin(phimax) - sin(phimin));
    double varx = ex2 - ex*ex;

    double ey2 = 1/(8*A) * (pow(rmax,4) - pow(rmin,4)) * (-phimin + sin(phimin)*cos(phimin) + phimax - sin(phimax)*cos(phimax));
    double ey = 1/(3*A) * (pow(rmax,3) - pow(rmin,3)) * (cos(phimin) - cos(phimax));
    double vary = ey2 - ey*ey;

    double varxy = 1/(16*A)*(pow(rmax,4)-pow(rmin,4))*(cos(2*phimin)-cos(2*phimax)) - ex*ey;
    return LocalError(varx, varxy, vary);
  }
} 

bool IsolationRequirementAnalyzer::inSensorCell(DetId detid_, GlobalPoint point){
  // inSensorCell checks both hexagonal cells and the trapezoid based on the corners stored in the hgcalGeometry.
  // The corners are stored in a counter clockwise fashion.
  // The hexagonal boundaries are checked using the winding number algorithm.
  // The trapezoid boundaries are checked using simple polar coordinates.

  // Get Coordinates
  auto x = point.x();
  auto y = point.y();
  
  // GetCorners
  auto hgcalgeometry = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_));
  auto corners = hgcalgeometry->getCorners(detid_);

  /*
  std::cout << "v1=["<<corners[0].x()<<","<<corners[0].y()<<","<<corners[0].z()<<"]" << std::endl;
  std::cout << "v2=["<<corners[1].x()<<","<<corners[1].y()<<","<<corners[1].z()<<"]" << std::endl;
  std::cout << "v3=["<<corners[2].x()<<","<<corners[2].y()<<","<<corners[2].z()<<"]" << std::endl;
  std::cout << "v4=["<<corners[3].x()<<","<<corners[3].y()<<","<<corners[3].z()<<"]" << std::endl;
  std::cout << "v5=["<<corners[4].x()<<","<<corners[4].y()<<","<<corners[4].z()<<"]" << std::endl;
  std::cout << "v6=["<<corners[5].x()<<","<<corners[5].y()<<","<<corners[5].z()<<"]" << std::endl;
  std::cout << "v7=["<<corners[6].x()<<","<<corners[6].y()<<","<<corners[6].z()<<"]" << std::endl;
  std::cout << "p=["<<point.x()<<","<<point.y()<<","<<point.z()<<"]" << std::endl;
  */
  bool inside = false;
  if (recHitTools_.isSilicon(detid_)){
    // Check Hexagon
    int cornersSize = 7;
    int j=0;
    for (int i = 1; i < cornersSize; i++) {
      if (((corners[i].y() >+ point.y()) != (corners[j].y() >= point.y())) &&
          (x <= (corners[j].x() - corners[i].x()) * (y - corners[i].y()) / (corners[j].y() - corners[i].y()) + corners[i].x())) {
          inside = !inside;
      }
      j++;
    }
  }
  else {
    // Check Scintillator
    float phi_min = std::atan2(corners[0].y(),corners[0].x()); //radians
    float phi_max = std::atan2(corners[1].y(),corners[1].x()); //radians
    float r_min = std::sqrt(corners[3].x()*corners[3].x() + corners[3].y()*corners[3].y());
    float r_max = std::sqrt(corners[0].x()*corners[0].x() + corners[0].y()*corners[0].y());;

    auto r = std::sqrt(point.x()*point.x()+point.y()*point.y());
    if (point.phi()>=phi_min && point.phi()<=phi_max){
      if (r>=r_min && r<=r_max){
        inside = !inside;
      }
    }
  }
  return inside;
}
void IsolationRequirementAnalyzer::clear_arrays(){

  // KF
  sig_x.clear();
  sig_y.clear();
  sig_z.clear();
  sig_e.clear();
  sig_detid.clear();
  sig_layer.clear();
  sig_evt.clear();
  sig_inLC.clear();
  sig_isolated.clear();
  sig_totalNeighbors.clear();
  sig_hitNeighbors.clear();
  sig_efrac.clear();
  sig_time.clear();


  nbor_x.clear();
  nbor_y.clear();
  nbor_z.clear();
  nbor_e.clear();
  nbor_detid.clear();
  nbor_layer.clear();
  nbor_evt.clear();
  nbor_time.clear();
  nbor_dtype.clear();

}



void IsolationRequirementAnalyzer::fillHitMap(std::map<DetId, std::pair<const HGCRecHit*, float>>& hitMap,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH,
                                const std::vector<float>& lcmask) const {
  hitMap.clear();
  int counter = 0;
  for (const auto& hit : rechitsEE) {
    hitMap.emplace(hit.detid(), std::make_pair(&hit, lcmask[counter]));
    counter++;
  }

  for (const auto& hit : rechitsFH) {
    hitMap.emplace(hit.detid(), std::make_pair(&hit, lcmask[counter]));
    counter++;
  }

  for (const auto& hit : rechitsBH) {
    hitMap.emplace(hit.detid(), std::make_pair(&hit, lcmask[counter]));
    counter++;
  }
} // end of EfficiencyStudies::fillHitMap


std::vector<int> IsolationRequirementAnalyzer::matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_) {
  std::vector<int> matchedIdxs; matchedIdxs.clear();
  for (unsigned int i0=0; i0<rechitdetid_.size(); ++i0) {
    if (detid_ == rechitdetid_[i0]) { matchedIdxs.push_back(i0); }
  }
  return matchedIdxs;
} // end of matchRecHit2CPRecHits

Dtype IsolationRequirementAnalyzer::getDtype(const DetId& detid_) const{
  Dtype dtype;
  if((detid_.det()==8)||(detid_.det()==9)){
    int thickness = recHitTools_.getSiThickness(detid_); 
    if (thickness == 120) dtype = Dtype::si120;
    else if (thickness == 200) dtype = Dtype::si200;
    else if (thickness == 300) dtype = Dtype::si300;
  }
  else if (detid_.det()==10){
    dtype = Dtype::sc;
  }
  else{
    dtype = Dtype::undefined;
  }


  return dtype;
}

std::string IsolationRequirementAnalyzer::writeDtype(Dtype& dtype) const{
  switch (dtype){
    case Dtype::si120:
      return "Si 120";
    case Dtype::si200:
      return "Si 200";
    case Dtype::si300:
      return "Si 300";
    case Dtype::sc:
      return "Sc";
    case Dtype::undefined:
    default:
      return "Undefined";
  }
}

void IsolationRequirementAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //std::cout << "Analyze" << std::endl;

  edm::Handle<std::vector<CaloParticle>> CaloParticles;
  iEvent.getByToken(caloParticlesToken_, CaloParticles);
  const CaloParticleCollection& cps = *CaloParticles;
  
  edm::Handle<HGCRecHitCollection> recHitHandleEE;  
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);

  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);

  edm::Handle<HGCRecHitCollection> recHitHandleBH;
  iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);

  edm::Handle<reco::CaloClusterCollection> layerClusterHandle;
  iEvent.getByToken(hgcalLayerClustersToken_, layerClusterHandle);
  const reco::CaloClusterCollection &lcs = *layerClusterHandle;

  edm::Handle<std::vector<float>> lcMaskHandle;
  iEvent.getByToken(lcMaskToken_, lcMaskHandle);

  std::map<DetId, std::pair<const HGCRecHit*, float>> hitMap;
  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH, *lcMaskHandle);

  edm::Handle<std::vector<KFHit>> KFHitsHandle;
  iEvent.getByToken(KFHitsToken_, KFHitsHandle);
  const std::vector<KFHit> &kfhits = *KFHitsHandle;

  edm::Handle<edm::View<reco::Track>> tracks_h;
  iEvent.getByToken(tracksToken_,tracks_h);
  const edm::View<reco::Track> & tkx = *(tracks_h.product()); 

  // Match tracks to simtrack
  
  reco::RecoToSimCollection const* recSimCollP = nullptr;
  edm::Handle<reco::RecoToSimCollection> recotosimCollectionH;
  iEvent.getByToken(associatormapRtSs[0], recotosimCollectionH);
  recSimCollP = recotosimCollectionH.product();
  reco::RecoToSimCollection const& recSimColl = *recSimCollP;

  edm::Handle<std::vector<SimVertex>> simVerticesHandle;
  iEvent.getByToken(simVerticesToken_, simVerticesHandle);
  std::vector<SimVertex> const& simVertices = *simVerticesHandle;

  TrackingParticleSelector tpSelector_ = TrackingParticleSelector(0, 100, 1.5, 3,120,280,0,true,false,false,false,{13});
  std::vector<int> signalIdx;

  // Create LC map of RecHits
  std::vector<DetId> lcs_detids;
  for (const auto& it_lc : lcs) {
    const std::vector<std::pair<DetId, float>> &hf = it_lc.hitsAndFractions();
    // loop over the rechits of this specific layer cluster
    for (unsigned int j = 0; j < hf.size(); j++) { 
      DetId detid_ = hf[j].first;  
      lcs_detids.push_back(detid_);
    }
  }

  for (edm::View<reco::Track>::size_type i = 0; i < tkx.size(); ++i) {
    RefToBase<reco::Track> track(tracks_h, i);  
        
    if(recSimColl.find(track) == recSimColl.end()) continue;

    std::vector<std::pair<TrackingParticleRef, double> > tps;
    tps = recSimColl[track];
    for (auto tp:tps){
      if(!tpSelector_(*(tp.first))) continue;
      signalIdx.push_back(i); 
    }
  }

  if (signalIdx.size()==0){
    std::cout << "No Signal found!!!!!!" << std::endl;
    return;
  }


    // init vars
  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  recHitTools_.setGeometry(geom);

  const CaloSubdetectorGeometry *geomHB = geom.getSubdetectorGeometry(DetId::Detector(10), ForwardSubdetector::ForwardEmpty);
  auto subgeomHB = static_cast<const HGCalGeometry*>(geomHB);
  auto topoHB = subgeomHB->topology();
  
  const CaloSubdetectorGeometry *geomHF = geom.getSubdetectorGeometry(DetId::Detector(9), ForwardSubdetector::ForwardEmpty);
  auto subgeomHF = static_cast<const HGCalGeometry*>(geomHF);
  auto topoHF = subgeomHF->topology();

  const CaloSubdetectorGeometry *geomEE = geom.getSubdetectorGeometry(DetId::Detector(8), ForwardSubdetector::ForwardEmpty);
  auto subgeomEE = static_cast<const HGCalGeometry*>(geomEE);
  auto topoEE = subgeomEE->topology();

  // const CaloSubdetectorGeometry *subGeom = geom.getSubdetectorGeometry(DetId::Detector(10), ForwardSubdetector::ForwardEmpty);
  //auto geomEE = static_cast<const HGCalGeometry*>(geomEE);
  //const HGCalDDDConstants* ddd = &(geomEE->topology().dddConstants());
  //auto radiusLayer = ddd->rangeRLayer(8, true);
  //std::cout << radiusLayer.first << ", " << radiusLayer.second << std::endl;


  // Define CP selector
  CaloParticleSelector cpSelector_ = CaloParticleSelector(0, 100, 1.5, 3,120,280,0,1000000,true,false,false,false,false,{13},-3.2,3.2);

  // Loop over Caloparticles 
  std::vector<DetId> tmprechits_; tmprechits_.clear();
  int obj_id = 0;
  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    const CaloParticle& cp = ((it_cp)); 

    //Select only Signal
    if (!cpSelector_(cp,simVertices)) continue;

    auto pid = cp.particleId();
    const SimClusterRefVector& simclusters = cp.simClusters();

    for (const auto& it_simc : simclusters){
      const SimCluster& simc = (*(it_simc));
      const auto& sc_hae = simc.hits_and_energies();

      for (const auto& it_sc_hae : sc_hae){

        DetId detid_ = (it_sc_hae.first);
        std::map<DetId,std::pair<const HGCRecHit *, float>>::const_iterator itcheck = hitMap.find(detid_);
        if (itcheck == hitMap.end()) continue;

        unsigned int layer_ = recHitTools_.getLayerWithOffset(detid_);
        float time = (*hitMap.find(detid_)).second.first->time();
        //float time = (*itcheck).second.first->time();
        

        // Get Neighbors
        std::vector<DetId> neighbors;
        if (detid_.det() == 8){
          neighbors = topoEE.neighbors(detid_);
        } else if (detid_.det() == 9){
          neighbors =  topoHF.neighbors(detid_);
        } else if (detid_.det() == 10){
          neighbors =  topoHB.neighbors(detid_);
        }
        int totalNeighbors = neighbors.size();

        // Check if in LC
        int inLC;
        auto hit_check = std::find(lcs_detids.begin(),lcs_detids.end(),detid_);
        if (hit_check != lcs_detids.end()){
          inLC = 1;
        } else {
          inLC = 0;
        }

        // Check if neighbors are RecHits
        int hitNeighbors = 0;
        float energy = (it_sc_hae.second);

        for (auto& neighbor: neighbors){
          auto candidate = hitMap.find(neighbor);
          if (candidate != hitMap.end()){
            hitNeighbors++;
            energy += (*candidate).second.first->energy();

            Dtype dtype = getDtype((*candidate).first);

            nbor_x.push_back(recHitTools_.getPosition(detid_).x());
            nbor_y.push_back(recHitTools_.getPosition(detid_).y());
            nbor_z.push_back(recHitTools_.getPosition(detid_).z());
            nbor_e.push_back(it_sc_hae.second);
            nbor_detid.push_back(detid_);
            nbor_layer.push_back(layer_);
            nbor_evt.push_back(eventnr);
            nbor_time.push_back(time);
            nbor_dtype.push_back(writeDtype(dtype));
            //std::cout << (*candidate).second.first->time() << ", ";
            //auto hit_check = std::find(lcs_detids.begin(),lcs_detids.end(),(*candidate).first);
            //if (hit_check != lcs_detids.end()){
            //  std::cout << 1 << ", ";
            //} else {
            //  std::cout << 0 << ", ";
            //}
          }
        }
        float efrac = it_sc_hae.second/energy;
        int isIsolated = 1;
        if (hitNeighbors!= 0) isIsolated = 0;

        //std::cout << layer_ << ", " << detid_.det() << ", " << neighbors.size() << ", " << hitNeighbors << ", " << efrac <<std::endl;
        // Calculate the energy fractions in the area
        sig_x.push_back(recHitTools_.getPosition(detid_).x());
        sig_y.push_back(recHitTools_.getPosition(detid_).y());
        sig_z.push_back(recHitTools_.getPosition(detid_).z());
        sig_e.push_back(it_sc_hae.second);
        sig_detid.push_back(detid_);
        sig_layer.push_back(layer_);
        sig_evt.push_back(eventnr);
        sig_inLC.push_back(inLC);
        sig_isolated.push_back(isIsolated);
        sig_totalNeighbors.push_back(totalNeighbors);
        sig_hitNeighbors.push_back(hitNeighbors);
        sig_efrac.push_back(efrac);
        sig_time.push_back(time);
      }
    }
    obj_id++;
  }
  
  signalTree->Fill();
  neighborTree->Fill();

  // RecHit

  sig_x.clear();
  sig_y.clear();
  sig_z.clear();
  sig_e.clear();
  sig_detid.clear();
  sig_layer.clear();
  sig_evt.clear();
  sig_inLC.clear();
  sig_isolated.clear();
  sig_totalNeighbors.clear();
  sig_hitNeighbors.clear();
  sig_efrac.clear();
  sig_time.clear(); 

  nbor_x.clear();
  nbor_y.clear();
  nbor_z.clear();
  nbor_e.clear();
  nbor_detid.clear();
  nbor_layer.clear();
  nbor_evt.clear();
  nbor_time.clear();
  nbor_dtype.clear();

  //clear_arrays();
  eventnr=eventnr+1;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void IsolationRequirementAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void IsolationRequirementAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void IsolationRequirementAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsolationRequirementAnalyzer);
