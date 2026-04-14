// -*- C++ -*-
//
// Package:    PlayGround/NeighborhoodAnalyzer
// Class:      NeighborhoodAnalyzer
//
/**\class NeighborhoodAnalyzer NeighborhoodAnalyzer.cc PlayGround/NeighborhoodAnalyzer/plugins/NeighborhoodAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Fri, 30 Jan 2026 21:37:34 GMT
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

#include <Eigen/Dense>

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

class NeighborhoodAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit NeighborhoodAnalyzer(const edm::ParameterSet&);
  ~NeighborhoodAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  std::vector<DetId> getNeighbors(DetId detid_,HGCalTopology& topoEE, HGCalTopology& topoFH, HGCalTopology& topoBH) const;
  std::pair<int,double> getNeighborHitsAndEnergies(std::vector<DetId> neighbors, std::map<DetId, std::pair<const HGCRecHit*, float>>& hitMap) const;
  virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap, 
    const HGCRecHitCollection& rechitsEE, 
    const HGCRecHitCollection& rechitsFH,
    const HGCRecHitCollection& rechitsBH) const;
  void fillCPMap(std::unordered_multimap<uint32_t, const CaloParticle*>& cpMap,
      const CaloParticleCollection& cps) const;
  void fillSCMap(std::unordered_multimap<uint32_t, const SimCluster*>& scMap,
    const SimClusterCollection& scs) const;
  std::string getDtype(int detid_) const;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<SimCluster>> simclustersToken_;

  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::EDGetTokenT<std::vector<PCaloHit>> pCaloHitEEToken_;
  edm::EDGetTokenT<std::vector<PCaloHit>> pCaloHitFHToken_;
  edm::EDGetTokenT<std::vector<PCaloHit>> pCaloHitBHToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  std::vector<edm::InputTag> associators;
  std::vector<edm::EDGetTokenT<reco::SimToRecoCollection>> associatormapStRs;
  std::vector<edm::EDGetTokenT<reco::RecoToSimCollection>> associatormapRtSs;
  edm::EDGetTokenT<std::vector<SimVertex>> simVerticesToken_;

  hgcal::RecHitTools recHitTools_;

  TTree* h_tree;
  TTree* n_tree;
  TTree* s_tree;

  std::vector<int> h_tree_evt;
  std::vector<int> h_tree_cp;
  std::vector<int> h_tree_cp_hits;
  std::vector<int> h_tree_hits;

  std::vector<int> n_tree_detid;
  std::vector<int> n_tree_evt;
  std::vector<int> n_tree_layer;
  std::vector<float> n_tree_eta;
  std::vector<float> n_tree_phi;
  std::vector<float> n_tree_energy;
  std::vector<std::string> n_tree_dtype;
  std::vector<int> n_tree_ndetid;
  std::vector<int> n_tree_npdgid;
  std::vector<float> n_tree_nenergy;
  std::vector<std::string> n_tree_ndtype;

  std::vector<int> s_tree_detid;
  std::vector<int> s_tree_evt;
  std::vector<int> s_tree_layer;
  std::vector<float> s_tree_eta;
  std::vector<float> s_tree_phi;
  std::vector<std::string> s_tree_dtype;
  std::vector<float> s_tree_energy_total;
  std::vector<int> s_tree_pdgid;
  std::vector<float> s_tree_energy_pdgid;


  std::map<DetId, const HGCRecHit*> hitMap;
  std::unordered_multimap<uint32_t, const CaloParticle*> cpMap;
  std::unordered_multimap<uint32_t, const SimCluster*> scMap;

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
NeighborhoodAnalyzer::NeighborhoodAnalyzer(const edm::ParameterSet& iConfig)
: caloParticlesToken_(consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticles"))), 
simclustersToken_(consumes<std::vector<SimCluster> >(iConfig.getParameter<edm::InputTag>("simclusters"))), 
hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
pCaloHitEEToken_(consumes<std::vector<PCaloHit>>(iConfig.getUntrackedParameter<edm::InputTag>("pCaloHitsEE"))),
pCaloHitFHToken_(consumes<std::vector<PCaloHit>>(iConfig.getUntrackedParameter<edm::InputTag>("pCaloHitsFH"))),
pCaloHitBHToken_(consumes<std::vector<PCaloHit>>(iConfig.getUntrackedParameter<edm::InputTag>("pCaloHitsBH"))),
//abs_failToken_(consumes<float>(iConfig.getParameter<edm::InputTag>("abs_fail"))),
caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
associators(iConfig.getUntrackedParameter<std::vector<edm::InputTag>>("associators")),
simVerticesToken_(consumes<std::vector<SimVertex>>(iConfig.getParameter<edm::InputTag>("simVertices"))){

associatormapRtSs.push_back(consumes<reco::RecoToSimCollection>(associators[0]));
usesResource("TFileService");
edm::Service<TFileService> file;
h_tree = file->make<TTree>("h_tree","h_tree");
n_tree = file->make<TTree>("n_tree","n_tree"); // Replace this with an input tree value?
s_tree = file->make<TTree>("s_tree","s_tree"); // Replace this with an input tree value?

h_tree->Branch("evt", &h_tree_evt);
h_tree->Branch("cp", &h_tree_cp);
h_tree->Branch("cp_hits",&h_tree_cp_hits);
h_tree->Branch("hits",&h_tree_hits);

n_tree->Branch("detid", &n_tree_detid);
n_tree->Branch("evt",&n_tree_evt);
n_tree->Branch("layer", &n_tree_layer);
n_tree->Branch("eta", &n_tree_eta);
n_tree->Branch("phi", &n_tree_phi);
n_tree->Branch("energy", &n_tree_energy);
n_tree->Branch("dtype",&n_tree_dtype);
n_tree->Branch("ndetid", &n_tree_ndetid);    
n_tree->Branch("npdgid", &n_tree_npdgid);    
n_tree->Branch("nenergy", &n_tree_nenergy);
n_tree->Branch("ndtype",&n_tree_ndtype);


s_tree->Branch("detid", &s_tree_detid);
s_tree->Branch("evt",&s_tree_evt);
s_tree->Branch("layer", &s_tree_layer);
s_tree->Branch("eta", &s_tree_eta);
s_tree->Branch("phi", &s_tree_phi);
s_tree->Branch("energy_total", &s_tree_energy_total);
s_tree->Branch("dtype",&s_tree_dtype);
s_tree->Branch("pdgid",&s_tree_pdgid);
s_tree->Branch("pdgid_energy",&s_tree_energy_pdgid);

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

NeighborhoodAnalyzer::~NeighborhoodAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

void NeighborhoodAnalyzer::fillCPMap(std::unordered_multimap<uint32_t, const CaloParticle*>& cpMap,
    const CaloParticleCollection& cps) const {

  cpMap.clear();
  for (const auto& cp: cps){
    const SimClusterRefVector& simclusters = cp.simClusters();

    for (const auto& it_simc : simclusters){
      const SimCluster& simc = (*(it_simc));
      const auto& sc_hae = simc.hits_and_energies();

      for (const auto& it_sc_hae : sc_hae){
        DetId detid_ = (it_sc_hae.first);
        cpMap.emplace(detid_.rawId(),&cp);
      }
    }
  } 
}

void NeighborhoodAnalyzer::fillSCMap(std::unordered_multimap<uint32_t, const SimCluster*>& scMap,
  const SimClusterCollection& scs) const {

scMap.clear();
for (const auto& sc: scs){
  const auto& sc_hae = sc.hits_and_energies();

  for (const auto& it_sc_hae : sc_hae){
    DetId detid_ = (it_sc_hae.first);
    scMap.emplace(detid_.rawId(),&sc);
    }
  }
}

void NeighborhoodAnalyzer::fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
  const HGCRecHitCollection& rechitsEE,
  const HGCRecHitCollection& rechitsFH,
  const HGCRecHitCollection& rechitsBH) const {
  hitMap.clear();
  int counter = 0;
  for (const auto& hit : rechitsEE) {
  hitMap.emplace(hit.detid(), &hit);
  counter++;
  }

  for (const auto& hit : rechitsFH) {
  hitMap.emplace(hit.detid(), &hit);
  counter++;
  }

  for (const auto& hit : rechitsBH) {
  hitMap.emplace(hit.detid(), &hit);
  counter++;
  }
} // end of EfficiencyStudies::fillHitMap


std::vector<DetId> NeighborhoodAnalyzer::getNeighbors(DetId detid_,HGCalTopology& topoEE,HGCalTopology& topoHF,HGCalTopology& topoHB) const {
  // Get neighbors of detid
  std::vector<DetId> neighbors;
  if (detid_.det() == 8){
    neighbors = topoEE.neighbors(detid_);
  } else if (detid_.det() == 9){
    neighbors =  topoHF.neighbors(detid_);
  } else if (detid_.det() == 10){
    neighbors =  topoHB.neighbors(detid_);
  }
  return neighbors;
}

std::pair<int,double> NeighborhoodAnalyzer::getNeighborHitsAndEnergies(std::vector<DetId> neighbors, std::map<DetId, std::pair<const HGCRecHit*, float>>& hitMap) const {
  std::pair<int,float> result;

  // Check if neighboring detids are associated with a RecHit
  int hitNeighbors = 0;
  float totalEnergy = 0;

  for (auto& neighbor: neighbors){
    auto candidate = hitMap.find(neighbor);
    if (candidate != hitMap.end()){
      hitNeighbors++;
      totalEnergy += candidate->second.first->energy();
    }
  }
  result.first = hitNeighbors;
  result.second = totalEnergy;
  return result;
}

std::string NeighborhoodAnalyzer::getDtype(int detid_) const{
  std::string detector;
  std::string thickness;
  std::string dtype;

  // TODO: if conditions necessary to deal with missing rechits. Not yet possible to determine area of detector properly
  //       Implement it in such a way that we can find the closest detid and from this determine silicon thickness.
  
  if (detid_==0){
    //auto closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(0))->getClosestCell(hits[i].center);
    //layer_ = recHitTools_.getLayerWithOffset(closest_detid); // Can get rid of this in future iterations as I read it out in KFHits anyway
    detector = "Sc";
    thickness = "None";
    dtype = "Sc";
  }
  else if((detid_==8)||(detid_==9)){
    //auto closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(hits[i].detid))->getClosestCellHex(hits[i].center,true);
    //layer_ = recHitTools_.getLayerWithOffset(closest_detid); // Can get rid of this in future iterations as I read it out in KFHits anyway
    detector = "Si";
    //thickness = std::to_string(int(recHitTools_.getSiThickness(closest_detid))); 
    dtype = detector;
  }
  else if (recHitTools_.isSilicon(detid_)){
    detector = "Si";
    thickness = std::to_string(int(recHitTools_.getSiThickness(detid_))); 
    dtype = detector+" "+thickness;
  }
  else{
    detector = "Sc";
    thickness = "None";
    dtype = "Sc";
  } 
  return dtype;
}


// ------------ method called for each event  ------------
void NeighborhoodAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  auto const& eventnr = iEvent.id().event();

  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  recHitTools_.setGeometry(geom);

  edm::Handle<std::vector<CaloParticle>> CaloParticles;
  iEvent.getByToken(caloParticlesToken_, CaloParticles);
  const CaloParticleCollection& cps = *CaloParticles;
  
  edm::Handle<std::vector<SimCluster>> SimClusters;
  iEvent.getByToken(simclustersToken_, SimClusters);
  const SimClusterCollection& scs = *SimClusters;

  edm::Handle<HGCRecHitCollection> recHitHandleEE;  
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);

  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);

  edm::Handle<HGCRecHitCollection> recHitHandleBH;
  iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);

  // Match tracks to simtrack

  reco::RecoToSimCollection const* recSimCollP = nullptr;
  edm::Handle<reco::RecoToSimCollection> recotosimCollectionH;
  iEvent.getByToken(associatormapRtSs[0], recotosimCollectionH);
  recSimCollP = recotosimCollectionH.product();
  reco::RecoToSimCollection const& recSimColl = *recSimCollP;

  edm::Handle<std::vector<SimVertex>> simVerticesHandle;
  iEvent.getByToken(simVerticesToken_, simVerticesHandle);
  std::vector<SimVertex> const& simVertices = *simVerticesHandle;

  CaloParticleSelector cpSelector_ = CaloParticleSelector(0, 100, 1.5, 3,120,280,0,1000000,true,false,false,false,false,{13},-3.2,3.2);

  // Init vars
  const CaloSubdetectorGeometry *geomHB = geom.getSubdetectorGeometry(DetId::Detector(10), ForwardSubdetector::ForwardEmpty);
  auto subgeomHB = static_cast<const HGCalGeometry*>(geomHB);
  auto topoHB = subgeomHB->topology();
  const HGCalDDDConstants* ddd = &(subgeomHB->topology().dddConstants());

  
  const CaloSubdetectorGeometry *geomHF = geom.getSubdetectorGeometry(DetId::Detector(9), ForwardSubdetector::ForwardEmpty);
  auto subgeomHF = static_cast<const HGCalGeometry*>(geomHF);
  auto topoHF = subgeomHF->topology();

  const CaloSubdetectorGeometry *geomEE = geom.getSubdetectorGeometry(DetId::Detector(8), ForwardSubdetector::ForwardEmpty);
  auto subgeomEE = static_cast<const HGCalGeometry*>(geomEE);
  auto topoEE = subgeomEE->topology();

  // Fill Maps

  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);
  //fillSCMap(scMap,scs);
  fillCPMap(cpMap,cps);

  // Pcalohit

  // for (const auto& pCaloHit : iEvent.get(pCaloHitEEToken_)){
  //   std::cout << pCaloHit.id() << "," << pCaloHit.geantTrackId() << ","<< pCaloHit.getPosition().x()<< ","<< pCaloHit.getPosition().y()<< ","<< pCaloHit.getPosition().z() << std::endl;
  // }

  // Loop over Caloparticles and get signal detids

  int obj_id = 0;
  std::cout << cps.size() << std::endl;
  std::vector<DetId> signal_detids;
  std::vector<DetId> all_detids;

  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    const CaloParticle& cp = ((it_cp)); 
    auto gps = it_cp.genParticles();
    auto pdgid = cp.pdgId();
    const SimClusterRefVector& simclusters = cp.simClusters();
    for (const auto& it_simc : simclusters){
      const SimCluster& simc = (*(it_simc));
      const auto& sc_hae = simc.hits_and_energies();

      for (const auto& it_sc_hae : sc_hae){
        DetId detid_ = (it_sc_hae.first);
        all_detids.push_back(detid_);
      }
    }
    //Select only Signal
    if (!cpSelector_(cp,simVertices)) continue;
    int simClusterId = 0;
    for (const auto& it_simc : simclusters){
      const SimCluster& simc = (*(it_simc));
      const auto& sc_hae = simc.hits_and_energies();

      for (const auto& it_sc_hae : sc_hae){
        DetId detid_ = (it_sc_hae.first);
        signal_detids.push_back(detid_);
      }
    }
  }

  auto last = std::unique(all_detids.begin(), all_detids.end());
  all_detids.erase(last, all_detids.end());

  h_tree_evt.push_back(eventnr);
  h_tree_cp.push_back(cps.size());
  h_tree_cp_hits.push_back(all_detids.size());
  h_tree_hits.push_back(((*recHitHandleBH).size() + (*recHitHandleFH).size() + (*recHitHandleEE).size()));

  std::cout << "Nr of CP hits: " << all_detids.size() << std::endl;
  std::cout << "Nr of RecHits: " << ((*recHitHandleBH).size() + (*recHitHandleFH).size() + (*recHitHandleEE).size()) << std::endl;
 

  // Loop over unique signal detids and get neighborhood
  std::sort(signal_detids.begin(), signal_detids.end());
  auto it = std::unique(signal_detids.begin(),signal_detids.end());
  signal_detids.resize( std::distance(signal_detids.begin(),it) );
  
  for (const auto& hit:signal_detids){
    auto layer = recHitTools_.getLayerWithOffset(hit);
    auto eta = recHitTools_.getEta(hit,0);
    auto phi = recHitTools_.getPhi(hit);
    auto energy = hitMap.find(hit)->second->energy();
    auto dtype = getDtype(hit);


    // Fill s_tree
    auto range = cpMap.equal_range(hit.rawId());
    for (auto it = range.first; it!= range.second; it++){
      auto cp = it->second;
      const SimClusterRefVector& simclusters = cp->simClusters();
      for (const auto& it_simc : simclusters){
        const SimCluster& simc = (*(it_simc));
        const auto& sc_haf = simc.hits_and_fractions();
        for (const auto& it_sc_haf : sc_haf){
          DetId detid_ = (it_sc_haf.first);
          if (detid_.rawId()==hit.rawId()){
            s_tree_detid.push_back(hit);
            s_tree_evt.push_back(eventnr);
            s_tree_layer.push_back(layer);
            s_tree_eta.push_back(eta);
            s_tree_phi.push_back(phi);
            s_tree_energy_total.push_back(energy);
            s_tree_dtype.push_back(dtype);
            s_tree_pdgid.push_back(cp->pdgId());
            s_tree_energy_pdgid.push_back(it_sc_haf.second*energy);
          }
        }
      }
    }

    // Fill n_tree

    auto neighbors = getNeighbors(hit,topoEE,topoHF,topoHB);
    //auto neighbors_pair = getNeighborHitsAndEnergies(neighbors, hitMap);
    for (const auto& neighbor:neighbors){
      auto candidate = hitMap.find(neighbor);
      if (candidate == hitMap.end()) continue;


      auto range = cpMap.equal_range(neighbor.rawId());

      // In case there is no CP, choose pdgid 0 and set nenergy as the energy of the candidate hit
      if (range.first == range.second) {
        // Store information on signal detid/detid for wich neighborhood is defined
        n_tree_detid.push_back(hit);
        n_tree_evt.push_back(eventnr);
        n_tree_layer.push_back(layer);
        n_tree_eta.push_back(eta);
        n_tree_phi.push_back(phi);
        n_tree_energy.push_back(energy);
        n_tree_dtype.push_back(dtype);

        // Store information on the neighborhood hit (candidate)
        n_tree_nenergy.push_back(candidate->second->energy());
        n_tree_ndetid.push_back(candidate->second->detid().rawId());
        n_tree_npdgid.push_back(0);
        n_tree_ndtype.push_back(getDtype(candidate->second->detid().rawId()));
      }

      // In case there are CP, loop over CP and store contribution
      for (auto it = range.first; it!= range.second; it++){
        auto cp = it->second;
        const SimClusterRefVector& simclusters = cp->simClusters();
        for (const auto& it_simc : simclusters){
          const SimCluster& simc = (*(it_simc));
          const auto& sc_haf = simc.hits_and_fractions();
          for (const auto& it_sc_haf : sc_haf){
            DetId detid_ = (it_sc_haf.first);
            if (detid_.rawId()==neighbor.rawId()){
              // Store information on signal detid/detid for wich neighborhood is defined
              n_tree_detid.push_back(hit);
              n_tree_evt.push_back(eventnr);
              n_tree_layer.push_back(layer);
              n_tree_eta.push_back(eta);
              n_tree_phi.push_back(phi);
              n_tree_energy.push_back(energy);
              n_tree_dtype.push_back(dtype);

              // Store information on the neighborhood hit (candidate)
              n_tree_nenergy.push_back(it_sc_haf.second*candidate->second->energy());
              n_tree_ndetid.push_back(detid_.rawId());
              n_tree_npdgid.push_back(cp->pdgId());
              n_tree_ndtype.push_back(getDtype(detid_.rawId()));
            } 
          }
        }
      }
    }
  }
  

  // // Loop over unique signal detids and get neighborhood
  // std::sort(signal_detids.begin(), signal_detids.end());
  // auto it = std::unique(signal_detids.begin(),signal_detids.end());
  // signal_detids.resize( std::distance(signal_detids.begin(),it) );
  // for (const auto& hit:signal_detids){
  //   auto neighbors = getNeighbors(hit,topoEE,topoHF,topoHB);
  //   //auto neighbors_pair = getNeighborHitsAndEnergies(neighbors, hitMap);
  //   for (const auto& neighbor:neighbors){
  //     auto range = scMap.equal_range(neighbor.rawId());
  //     if (range.first == range.second) {
  //       std::cout << "No Simcluster found for DetId: " << neighbor.rawId() << std::endl;
  //       continue;
  //     }
  //     for (auto it = range.first; it!= range.second; it++){
  //       auto simclusters = it->second;
  //       const auto& sc_hae = simclusters->hits_and_energies();
  //       for (const auto& it_sc_hae : sc_hae){
  //         DetId detid_ = (it_sc_hae.first);
  //         if (detid_.rawId()==neighbor.rawId()){
  //           std::cout << detid_.rawId() << std::endl;
  //         }
  //       }
  //     }
  //   }
  // }

  h_tree->Fill();
  h_tree_evt.clear();
  h_tree_cp.clear();
  h_tree_cp_hits.clear();
  h_tree_hits.clear();

  n_tree->Fill();
  n_tree_detid.clear();
  n_tree_evt.clear();
  n_tree_layer.clear();
  n_tree_eta.clear();
  n_tree_phi.clear();
  n_tree_energy.clear();
  n_tree_dtype.clear();
  n_tree_ndetid.clear();
  n_tree_npdgid.clear();
  n_tree_nenergy.clear();
  n_tree_ndtype.clear();

  s_tree->Fill();
  s_tree_detid.clear();
  s_tree_evt.clear();
  s_tree_layer.clear();
  s_tree_eta.clear();
  s_tree_phi.clear();
  s_tree_energy_total.clear();
  s_tree_dtype.clear();
  s_tree_pdgid.clear();
  s_tree_energy_pdgid.clear();



#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void NeighborhoodAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void NeighborhoodAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void NeighborhoodAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(NeighborhoodAnalyzer);
