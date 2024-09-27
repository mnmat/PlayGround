// -*- C++ -*-
//
// Package:    PlayGround/SearchWindowAnalyzer
// Class:      SearchWindowAnalyzer
//
/**\class SearchWindowAnalyzer SearchWindowAnalyzer.cc PlayGround/SearchWindowAnalyzer/plugins/SearchWindowAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Wed, 16 Aug 2023 16:53:55 GMT
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

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HGCalReco/interface/TICLLayerTile.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"

#include "Validation/HGCalValidation/interface/CaloParticleSelector.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "DataFormats/HGCalReco/interface/KFHit.h"

#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"



//ROOT includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TProfile.h"
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

class SearchWindowAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SearchWindowAnalyzer(const edm::ParameterSet&);
  ~SearchWindowAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void dumpTiles(const TICLLayerTiles &tiles);
  void dumpTiles(const TICLLayerTiles &tiles,std::vector<const HGCRecHit*>& hitVector);
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  virtual void fetchDetIdsFromSignal(std::vector<DetId>& output,
    std::vector<CaloParticle> cps,
    std::vector<SimVertex> simvertices) const;
  std::map<DetId,int> measurements(
      const DetId &detId, 
      std::vector<const HGCRecHit*>& hitVector,
      const TICLLayerTiles &tiles) const;
  virtual void mergeRecHitCollections(std::vector<const HGCRecHit*>& recHitCollection,
    const HGCRecHitCollection& recHitsEE,
    const HGCRecHitCollection& recHitsFH,
    const HGCRecHitCollection& recHitsBH) const;
  bool isValid(const DetId& detid,
            std::vector<const HGCRecHit*>& hitVector) const;
  void fetchDetIdsFromSearchWindow(std::map<DetId,int>& output, 
                                      std::vector<const HGCRecHit*>& hitVector,
                                      const std::vector<CaloParticle>& cps,
                                      const std::vector<SimVertex>& simVertices,
                                      const TICLLayerTiles& tiles) const;
  void fetchDetIdsFromKFHits(std::map<DetId,int>& output, 
                                      std::vector<const HGCRecHit*>& hitVector,
                                      const std::vector<KFHit>& kfhits,
                                      const TICLLayerTiles& tiles) const;

  // ----------member data ---------------------------
  edm::EDGetTokenT<TICLLayerTiles> recHitTilesToken_;
  edm::EDGetTokenT<TICLLayerTiles> layerTilesToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> lcToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<SimVertex>> simVerticesToken_;
  edm::EDGetTokenT<std::vector<float>> lcMaskToken_;
  edm::EDGetTokenT<std::vector<KFHit>> KFHitsToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;


  static constexpr float etaBinSize = (TICLLayerTiles::constants_type_t::maxEta - TICLLayerTiles::constants_type_t::minEta)/TICLLayerTiles::constants_type_t::nEtaBins;
  static constexpr float phiBinSize = 2*M_PI/TICLLayerTiles::constants_type_t::nPhiBins;
  static constexpr int nPhiBin = TICLLayerTiles::constants_type_t::nPhiBins;
  static constexpr int nEtaBin = TICLLayerTiles::constants_type_t::nEtaBins;

  int eventnr = 0;
  hgcal::RecHitTools rhtools_;

  // Histograms
  TTree *recHitTree;
  TTree *kfHitTree;
  TTree *swTree;

  std::vector<float> recHit_x;
  std::vector<float> recHit_y;
  std::vector<float> recHit_z;
  std::vector<float> recHit_e;
  std::vector<int> recHit_detid;
  std::vector<int> recHit_layer;
  std::vector<int> recHit_evt;
  std::vector<int> recHit_obj_id;
  std::vector<std::string> recHit_dtype;

  std::vector<float> kfHit_x;
  std::vector<float> kfHit_y;
  std::vector<float> kfHit_z;
  std::vector<float> kfHit_e;
  std::vector<int> kfHit_detid;
  std::vector<int> kfHit_layer;
  std::vector<int> kfHit_evt;
  std::vector<int> kfHit_obj_id;
  std::vector<std::string> kfHit_dtype;

  std::vector<int> sw_layer;
  std::vector<float> sw_eta;
  std::vector<float> sw_phi;
  std::vector<int> sw_globalbin;
  std::vector<int> sw_ieta;
  std::vector<int> sw_iphi;
  std::vector<int> sw_evt;
  std::vector<int> sw_size;

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
SearchWindowAnalyzer::SearchWindowAnalyzer(const edm::ParameterSet& iConfig)
    : recHitTilesToken_(consumes<TICLLayerTiles>(iConfig.getParameter<edm::InputTag>("recHitTiles"))),
    layerTilesToken_(consumes<TICLLayerTiles>(iConfig.getParameter<edm::InputTag>("layerTiles"))),
    lcToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("layerClusters"))),
    hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
    hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
    hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
    caloParticlesToken_(consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("caloParticles"))),
    simVerticesToken_(consumes<std::vector<SimVertex>>(iConfig.getParameter<edm::InputTag>("simVertices"))),
    lcMaskToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("lcMask"))),
    KFHitsToken_(consumes<std::vector<KFHit>>(iConfig.getParameter<edm::InputTag>("KFHits"))),
    caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>())
  {

  //dTypes = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};
  //cpHits = {"Simhits", "Rechits"};

  usesResource("TFileService");
  edm::Service<TFileService> file;
  recHitTree = file->make<TTree>("RecHit","RecHit");
  swTree = file->make<TTree>("SearchWindow","SearchWindow");
  kfHitTree = file->make<TTree>("KFHit","KFHit");


  recHitTree->Branch("x", &recHit_x);
  recHitTree->Branch("y", &recHit_y);
  recHitTree->Branch("z", &recHit_z);
  recHitTree->Branch("e", &recHit_e);
  recHitTree->Branch("layer", &recHit_layer);
  recHitTree->Branch("detid", &recHit_detid);
  recHitTree->Branch("dtype", &recHit_dtype);
  recHitTree->Branch("evt", &recHit_evt);
  recHitTree->Branch("obj_id", &recHit_obj_id);

  kfHitTree->Branch("x", &kfHit_x);
  kfHitTree->Branch("y", &kfHit_y);
  kfHitTree->Branch("z", &kfHit_z);
  kfHitTree->Branch("e", &kfHit_e);
  kfHitTree->Branch("layer", &kfHit_layer);
  kfHitTree->Branch("detid", &kfHit_detid);
  kfHitTree->Branch("dtype", &kfHit_dtype);
  kfHitTree->Branch("evt", &kfHit_evt);
  kfHitTree->Branch("obj_id", &kfHit_obj_id);

  swTree->Branch("layer", &sw_layer);
  swTree->Branch("ieta", &sw_ieta);
  swTree->Branch("iphi", &sw_iphi);
  swTree->Branch("globalbin", &sw_globalbin);
  swTree->Branch("eta", &sw_eta);
  swTree->Branch("phi", &sw_phi);
  swTree->Branch("size", &sw_size);
  swTree->Branch("evt",&sw_evt);


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

SearchWindowAnalyzer::~SearchWindowAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty

  
}

//
// member functions
//

void SearchWindowAnalyzer::mergeRecHitCollections(std::vector<const HGCRecHit*>& recHitCollection,
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

void SearchWindowAnalyzer::dumpTiles(const TICLLayerTiles &tiles,
                            std::vector<const HGCRecHit*>& hitMap) {
  constexpr int nEtaBin = TICLLayerTiles::constants_type_t::nEtaBins;
  constexpr int nPhiBin = TICLLayerTiles::constants_type_t::nPhiBins;
  auto lastLayerPerSide = static_cast<int>(rhtools_.lastLayer(false));
  int maxLayer = 2 * lastLayerPerSide - 1;
  for (int layer = 0; layer <= maxLayer; layer++) {
    for (int ieta = 0; ieta <= nEtaBin; ieta++) {
      auto offset = ieta * nPhiBin;
      for (int phi = 0; phi < nPhiBin; phi++) {
        int iphi = ((phi % nPhiBin + nPhiBin) % nPhiBin);
        if (!tiles[layer][offset + iphi].empty()) {
          auto phi = 2*M_PI/nPhiBin*iphi;
          auto eta = (TICLLayerTiles::constants_type_t::maxEta - TICLLayerTiles::constants_type_t::minEta)/nEtaBin*ieta+TICLLayerTiles::constants_type_t::minEta;
          
          for(auto hit: tiles[layer][offset + iphi]) {          
            auto detid = hitMap[hit]->detid();
            std::cout << eventnr << "," << detid.rawId() << "," << layer << "," << offset + iphi << "," << eta << "," << ieta << "," << phi << "," << iphi << std::endl;
          }
        }
      }
    }
  }
}

void SearchWindowAnalyzer::dumpTiles(const TICLLayerTiles &tiles) {
  constexpr int nEtaBin = TICLLayerTiles::constants_type_t::nEtaBins;
  constexpr int nPhiBin = TICLLayerTiles::constants_type_t::nPhiBins;
  auto lastLayerPerSide = static_cast<int>(rhtools_.lastLayer(false));
  int maxLayer = 2 * lastLayerPerSide - 1;
  for (int layer = 0; layer <= maxLayer; layer++) {
    for (int ieta = 0; ieta <= nEtaBin; ieta++) {
      auto offset = ieta * nPhiBin;
      for (int phi = 0; phi < nPhiBin; phi++) {
        int iphi = ((phi % nPhiBin + nPhiBin) % nPhiBin);
        if (!tiles[layer][offset + iphi].empty()) {

          auto phi = 2*M_PI/nPhiBin*iphi;
          auto eta = (TICLLayerTiles::constants_type_t::maxEta - TICLLayerTiles::constants_type_t::minEta)/nEtaBin*ieta+TICLLayerTiles::constants_type_t::minEta;
          sw_layer.push_back(layer);
          sw_ieta.push_back(ieta);
          sw_iphi.push_back(iphi);
          sw_globalbin.push_back(offset+iphi);
          sw_eta.push_back(eta);
          sw_phi.push_back(phi);
          sw_evt.push_back(eventnr);
          sw_size.push_back(tiles[layer][offset+iphi].size());
        }
      }
    }
  }
}

void SearchWindowAnalyzer::fetchDetIdsFromSignal(std::vector<DetId>& output, 
                                      std::vector<CaloParticle> cps,
                                      std::vector<SimVertex> simVertices) const{

  CaloParticleSelector cpSelector_ = CaloParticleSelector(0, 100, 1.5, 3,120,280,0,1000000,true,false,false,false,false,{13},-3.2,3.2);

  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    const CaloParticle& cp = ((it_cp)); 
    if (cpSelector_(cp,simVertices)){
      const SimClusterRefVector& simclusters = cp.simClusters();
      for (const auto& it_simc : simclusters){
        const SimCluster& simc = (*(it_simc));
        const auto& sc_hae = simc.hits_and_energies();
        for (const auto& it_sc_hae : sc_hae){
          output.push_back(it_sc_hae.first);
        }
      }
    }
  }
}

bool SearchWindowAnalyzer::isValid(const DetId& detid,
                                std::vector<const HGCRecHit*>& hitVector) const { 
  auto it = std::find_if(hitVector.begin(), hitVector.end(), [&detid](const HGCRecHit* hit) {
    return hit->detid() == detid;
    });
  return !(it == hitVector.end());
}

std::map<DetId,int> SearchWindowAnalyzer::measurements(
      const DetId &detId, 
      std::vector<const HGCRecHit*>& hitMap,
      const TICLLayerTiles &tiles) const{
  
  // Define output

  std::map<DetId,int> ret;

  // define search window and get bins
  if (detId.rawId() == 0) return ret;
  auto lastLayerPerSide = static_cast<int>(rhtools_.lastLayer(false));
  int layer = rhtools_.getLayerWithOffset(detId)+lastLayerPerSide -1; // Do the correct layre
  float eta = rhtools_.getPosition(detId).eta();
  float phi = rhtools_.getPosition(detId).phi();

  float etaMin = eta - etaBinSize;
  float etaMax = eta + etaBinSize;
  float phiMin = phi - phiBinSize;
  float phiMax = phi + phiBinSize;

  auto bins = tiles[layer].searchBoxEtaPhi(etaMin, etaMax, phiMin, phiMax);

  // loop over candidates
  for (int ieta = bins[0]; ieta <= bins[1]; ieta++) {
    auto offset = ieta * nPhiBin;
    for (int phi = bins[2]; phi <= bins[3]; phi++) {
      int iphi = ((phi % nPhiBin + nPhiBin) % nPhiBin);
      if (!tiles[layer][offset + iphi].empty()) {
        for(auto hit: tiles[layer][offset + iphi]) {          
          auto detid = hitMap[hit]->detid();
          ret.emplace(detid,offset+iphi);
        }
      }
    }
  }

  return ret;
}

void SearchWindowAnalyzer::fetchDetIdsFromSearchWindow(std::map<DetId,int>& output, 
                                      std::vector<const HGCRecHit*>& hitVector,
                                      const std::vector<CaloParticle>& cps,
                                      const std::vector<SimVertex>& simVertices,
                                      const TICLLayerTiles& tiles) const{

  CaloParticleSelector cpSelector_ = CaloParticleSelector(0, 100, 1.5, 3,120,280,0,1000000,true,false,false,false,false,{13},-3.2,3.2);

  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    const CaloParticle& cp = ((it_cp)); 
    if (cpSelector_(cp,simVertices)){
      const SimClusterRefVector& simclusters = cp.simClusters();
      for (const auto& it_simc : simclusters){
        const SimCluster& simc = (*(it_simc));
        const auto& sc_hae = simc.hits_and_energies();
        for (const auto& it_sc_hae : sc_hae){
          std::map<DetId,int> ms = measurements(it_sc_hae.first, hitVector, tiles);
          for (const auto& m : ms){
            output.emplace(m.first, m.second);
          }
        }
      }
    }
  }
}

void SearchWindowAnalyzer::fetchDetIdsFromKFHits(std::map<DetId,int>& output,
                                      std::vector<const HGCRecHit*>& hitVector,
                                      const std::vector<KFHit>& kfhits,
                                      const TICLLayerTiles& tiles) const{

  for (auto& kfhit: kfhits){
    std::map<DetId,int> ms = measurements(kfhit.detid,hitVector,tiles);
    for (const auto& m : ms){
      output.emplace(m.first,m.second);
    }
  }
}


// ------------ method called for each event  ------------
void SearchWindowAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<TICLLayerTiles> recHitTilesHandle;
  iEvent.getByToken(recHitTilesToken_, recHitTilesHandle);

  edm::Handle<HGCRecHitCollection> recHitHandleEE;  
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);

  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);

  edm::Handle<HGCRecHitCollection> recHitHandleBH;
  iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);

  edm::Handle<std::vector<CaloParticle>> CaloParticles;
  iEvent.getByToken(caloParticlesToken_, CaloParticles);

  edm::Handle<std::vector<SimVertex>> simVerticesHandle;
  iEvent.getByToken(simVerticesToken_, simVerticesHandle);

  edm::Handle<std::vector<KFHit>> KFHitsHandle;
  iEvent.getByToken(KFHitsToken_, KFHitsHandle);
  const std::vector<KFHit> &kfhits = *KFHitsHandle;


  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  rhtools_.setGeometry(geom);



  // dump all DetIds

  auto detids = geom.getValidDetIds();
  int counter_si = 0;
  int counter_sc = 0;

  for (auto detid:detids){
    if (!rhtools_.isSilicon(detid)){
      counter_sc++;
    }
    else counter_si ++;
  }

  std::vector<const HGCRecHit*> recHitCollection;
  mergeRecHitCollections(recHitCollection, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);
 
  // Get DetIds from Signal
  std::vector<DetId> signalDetIds;
  fetchDetIdsFromSignal(signalDetIds, *CaloParticles, *simVerticesHandle);

  // Get DetIds from SearchWindow
  std::map<DetId, int> sWMap;
  fetchDetIdsFromSearchWindow(sWMap, recHitCollection, *CaloParticles, *simVerticesHandle, *recHitTilesHandle);

  // Get DetIds from KFHits
  std::map<DetId, int> kfMap;
  fetchDetIdsFromKFHits(kfMap, recHitCollection, kfhits,*recHitTilesHandle);

  // Check how many DetIds are in the search window
  int signalInWindow = 0;
  for (auto &ele : sWMap){
    //auto energy = 
    auto id = ele.first;
    auto globalId = ele.second;
    auto l = rhtools_.getLayerWithOffset(id);

    recHit_x.push_back(rhtools_.getPosition(id).x());
    recHit_y.push_back(rhtools_.getPosition(id).y());
    recHit_z.push_back(rhtools_.getPosition(id).z());
    //e.push_back(energy);
    recHit_detid.push_back(id);
    recHit_layer.push_back(l);
    //dtype.push_back(tmp);
    recHit_evt.push_back(eventnr);
    recHit_obj_id.push_back(globalId);
  }

  for (auto &ele : kfMap){
     //auto energy = 
    auto id = ele.first;
    auto globalId = ele.second;
    auto l = rhtools_.getLayerWithOffset(id);

    kfHit_x.push_back(rhtools_.getPosition(id).x());
    kfHit_y.push_back(rhtools_.getPosition(id).y());
    kfHit_z.push_back(rhtools_.getPosition(id).z());
    //e.push_back(energy);
    kfHit_detid.push_back(id);
    kfHit_layer.push_back(l);
    //dtype.push_back(tmp);
    kfHit_evt.push_back(eventnr);
    kfHit_obj_id.push_back(globalId);   
  }

  dumpTiles(*recHitTilesHandle);
  //dumpTiles(*recHitTilesHandle,recHitCollection);

  eventnr++;

  recHitTree->Fill();
  recHit_x.clear();
  recHit_y.clear();
  recHit_z.clear();
  recHit_detid.clear();
  recHit_layer.clear();
  recHit_evt.clear();
  recHit_obj_id.clear();

  kfHitTree->Fill();
  kfHit_x.clear();
  kfHit_y.clear();
  kfHit_z.clear();
  kfHit_detid.clear();
  kfHit_layer.clear();
  kfHit_evt.clear();
  kfHit_obj_id.clear();

  swTree->Fill();
  sw_eta.clear();
  sw_phi.clear();
  sw_globalbin.clear();
  sw_size.clear();
  sw_ieta.clear();
  sw_iphi.clear();
  sw_layer.clear();
  sw_evt.clear();

  /*

  // Check how many DetIds from the Signal are reconstructed as LCs
  int signalInLcCount = 0;
  for (auto &signal : signalDetIds){
    int layer = rhtools_.getLayerWithOffset(signal)-1;
    h_signal_layer[h][dtype].front()->Fill(layer,1);
    if(std::find(lcDetIds.begin(), lcDetIds.end(), signal) != lcDetIds.end()){
      h_signal_lc_layer[h][dtype].front()->Fill(layer,1);
      signalInLcCount++; 
      if(detIdMask[signal]==0){
        h_signal_filtered_layer[h][dtype].front()->Fill(layer,1);
      }
    }
  }

  // Check how many DetIds are found in the event
  for (auto &detid : hitVector){
    int layer = rhtools_.getLayerWithOffset(detid)-1;
    h_hits_layer[h][dtype].front()->Fill(layer,1);
    if(detIdMask[detid]==0){
      h_hits_filtered_layer[h][dtype].front()->Fill(layer,1);
    }
  }

  // Apply Filter to generate Mask
  std::cout << "Nr of DetIds in Evt: " << hitVector.size() << std::endl;
  std::cout << "Nr of DetIds in LCs: " << lcDetIds.size() << std::endl;
  std::cout << "Nr of DetIds in Signal: " << signalDetIds.size() << std::endl;
  std::cout << "Nr of DetIds from Signal in LCs: " << signalInLcCount << std::endl;
  std::cout << "Nr of DetIds in Search Window: " << signalInWindow << std::endl;

  std::cout << "Nr of LCs: " << (*lcHandle).size() << std::endl;
  std::cout << "Size of LCMask: " << lcMask.size() << std::endl;

  // Count number of unique elements in mask
  std::sort(lcMask.begin(),lcMask.end());
  int count = 0;
  for (std::vector<double>::size_type i = 0; i < lcMask.size(); i++) {
      while (i < lcMask.size() - 1 && lcMask[i] == lcMask[i + 1])
      {
        i++;
      }
    count++;
  } 

  std::cout << "Unique elements in mask: " << lcMask.front() << ", " << lcMask.back() << std::endl; 
  //dumpTiles(*recHitTilesHandle);
  */


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void SearchWindowAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void SearchWindowAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SearchWindowAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(SearchWindowAnalyzer);
