// -*- C++ -*-
//
// Package:    PlayGround/PCaloHitAnalyzer
// Class:      PCaloHitAnalyzer
//
/**\class PCaloHitAnalyzer PCaloHitAnalyzer.cc PlayGround/PCaloHitAnalyzer/plugins/PCaloHitAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Tue, 31 Oct 2023 13:05:54 GMT
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

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/MuonNumbering/interface/MuonSubDetector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

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

class PCaloHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PCaloHitAnalyzer(const edm::ParameterSet&);
  ~PCaloHitAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap, 
    const HGCRecHitCollection& rechitsEE, 
    const HGCRecHitCollection& rechitsFH,
    const HGCRecHitCollection& rechitsBH) const;

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;

  edm::EDGetTokenT<std::vector<PSimHit>> pSimHitToken_;
  edm::EDGetTokenT<std::vector<PCaloHit>> pCaloHitEEToken_;
  edm::EDGetTokenT<std::vector<PCaloHit>> pCaloHitFHToken_;
  edm::EDGetTokenT<std::vector<PCaloHit>> pCaloHitBHToken_;

  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;

  
  edm::ESGetToken<CSCGeometry,MuonGeometryRecord> cscGeomToken_;
  //edm::ESGetToken<TrackerGeometry,TrackerDigiGeometryRecord> tecGeomToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  hgcal::RecHitTools recHitTools_;

  int eventnr =0;

  TTree* tree;

  std::vector<float> pcalo_x;
  std::vector<float> pcalo_y;
  std::vector<float> pcalo_z;
  std::vector<float> pcalo_e;
  std::vector<int> pcalo_detid;
  std::vector<int> pcalo_layer;
  std::vector<int> pcalo_evt;

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
PCaloHitAnalyzer::PCaloHitAnalyzer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
      caloParticlesToken_(consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticles"))),
      pSimHitToken_(consumes<std::vector<PSimHit>>(iConfig.getUntrackedParameter<edm::InputTag>("pSimHits"))),
      pCaloHitEEToken_(consumes<std::vector<PCaloHit>>(iConfig.getUntrackedParameter<edm::InputTag>("pCaloHitsEE"))),
      pCaloHitFHToken_(consumes<std::vector<PCaloHit>>(iConfig.getUntrackedParameter<edm::InputTag>("pCaloHitsFH"))),
      pCaloHitBHToken_(consumes<std::vector<PCaloHit>>(iConfig.getUntrackedParameter<edm::InputTag>("pCaloHitsBH"))),
      hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
      hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
      hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
      //tecGeomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>()),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>())
      {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> file;
  tree = file->make<TTree>("PCaloHit","PCaloHit");

  tree->Branch("x", &pcalo_x);
  tree->Branch("y", &pcalo_y);
  tree->Branch("z", &pcalo_z);
  tree->Branch("e", &pcalo_e);
  tree->Branch("layer", &pcalo_layer);
  tree->Branch("detid", &pcalo_detid);
  tree->Branch("evt", &pcalo_evt);

}

PCaloHitAnalyzer::~PCaloHitAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

void PCaloHitAnalyzer::fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH) const {
  hitMap.clear();
  for (const auto& hit : rechitsEE) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsFH) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsBH) {
    hitMap.emplace(hit.detid(), &hit);
  }
} // end of EfficiencyStudies::fillHitMap


// ------------ method called for each event  ------------
void PCaloHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //const TrackerGeometry &geom = iSetup.getData(tecGeomToken_);
  const CaloGeometry &calogeom = iSetup.getData(caloGeomToken_);
  recHitTools_.setGeometry(calogeom);

  edm::Handle<std::vector<CaloParticle>> CaloParticles;
  iEvent.getByToken(caloParticlesToken_, CaloParticles);
  const CaloParticleCollection& cps = *CaloParticles;

  edm::Handle<HGCRecHitCollection> recHitHandleEE;  
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);

  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);

  edm::Handle<HGCRecHitCollection> recHitHandleBH;
  iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);

  std::map<DetId, const HGCRecHit*> hitMap;
  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);



  for (const auto& track : iEvent.get(tracksToken_)) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();
  }


  
  //for (const auto& pSimHit : iEvent.get(pSimHitToken_)){
    //std::cout << "DetUnitId: " << pSimHit.detUnitId() << std::endl;
    //std::cout << "DetIdPosition: " << geom.idToDet(pSimHit.detUnitId())->position() << std::endl;
    //std::cout << "LocalPosition: "<< pSimHit.localPosition() << std::endl;
  //}
  
  int pcalocounter = 0;
  std::vector<DetId> detids;
  
  //std::cout << "Detid,TrackId, x, y, z" << std::endl;

  for (const auto& pCaloHit : iEvent.get(pCaloHitEEToken_)){
    //std::cout << pCaloHit.id() << "," << pCaloHit.geantTrackId() << ","<< pCaloHit.getPosition().x()<< ","<< pCaloHit.getPosition().y()<< ","<< pCaloHit.getPosition().z() << std::endl;
    
    pcalo_x.push_back(pCaloHit.getPosition().x());
    pcalo_y.push_back(pCaloHit.getPosition().y());
    pcalo_z.push_back(pCaloHit.getPosition().z());
    pcalo_e.push_back(pCaloHit.energy());
    pcalo_layer.push_back(recHitTools_.getLayerWithOffset(pCaloHit.id()));
    pcalo_detid.push_back(pCaloHit.id());
    pcalo_evt.push_back(eventnr);
    //std::cout << "DetIdPosition: " << geom.idToDet(pSimHit.detUnitId())->position() << std::endl;
    auto detid = pCaloHit.id();
    pcalocounter++;
    //if (std::find(detids.begin(), detids.end(), detid) == detids.end()){pcalocounter++;}
    detids.push_back(detid);

  }

  for (const auto& pCaloHit : iEvent.get(pCaloHitFHToken_)){
    //std::cout << pCaloHit.id() << "," << pCaloHit.geantTrackId() << ","<< pCaloHit.getPosition().x()<< ","<< pCaloHit.getPosition().y()<< ","<< pCaloHit.getPosition().z() << std::endl;
    pcalo_x.push_back(pCaloHit.getPosition().x());
    pcalo_y.push_back(pCaloHit.getPosition().y());
    pcalo_z.push_back(pCaloHit.getPosition().z());
    pcalo_e.push_back(pCaloHit.energy());
    pcalo_layer.push_back(recHitTools_.getLayerWithOffset(pCaloHit.id()));
    pcalo_detid.push_back(pCaloHit.id());
    pcalo_evt.push_back(eventnr);

    //std::cout << "DetIdPosition: " << geom.idToDet(pSimHit.detUnitId())->position() << std::endl;
    auto detid = pCaloHit.id();
    pcalocounter++;
    //if (std::find(detids.begin(), detids.end(), detid) == detids.end()){pcalocounter++;}
    detids.push_back(detid);
  }

  for (const auto& pCaloHit : iEvent.get(pCaloHitBHToken_)){
    //std::cout << pCaloHit.id() << "," << pCaloHit.geantTrackId() << ","<< pCaloHit.getPosition().x()<< ","<< pCaloHit.getPosition().y()<< ","<< pCaloHit.getPosition().z() << std::endl;
    pcalo_x.push_back(pCaloHit.getPosition().x());
    pcalo_y.push_back(pCaloHit.getPosition().y());
    pcalo_z.push_back(pCaloHit.getPosition().z());
    pcalo_e.push_back(pCaloHit.energy());
    pcalo_layer.push_back(recHitTools_.getLayerWithOffset(pCaloHit.id()));
    pcalo_detid.push_back(pCaloHit.id());
    pcalo_evt.push_back(eventnr);
    //std::cout << "DetIdPosition: " << geom.idToDet(pSimHit.detUnitId())->position() << std::endl;
    auto detid = pCaloHit.id();
    pcalocounter++;
    //if (std::find(detids.begin(), detids.end(), detid) == detids.end()){pcalocounter++;}
    detids.push_back(detid);
  }

  /*
  for (const auto& id:detids){
    const GlobalPoint &pos = recHitTools_.getPosition(id);
    std::cout << id.rawId() << "," << 1 << "," << pos.x() << "," << pos.y() << "," << pos.z()  << std::endl;
  }
  */
  
  int counter = 0;
  for (const auto& it_cp : cps) {
    const CaloParticle& cp = ((it_cp)); 
    const SimClusterRefVector& simclusters = cp.simClusters();

    for (const auto& it_simc : simclusters){
      const SimCluster& simc = (*(it_simc));
      const auto& sc_hae = simc.hits_and_energies();
      for (const auto& it_sc_hae : sc_hae){
        // Characterize RecHit
        //DetId detid_ = (it_sc_hae.first);
        counter++;
      }
    }
  }
  //std::cout << "Size of PCaloHits: " << pcalocounter << std::endl;
  //std::cout << "Size of CaloParticle: " << counter << std::endl;


  tree->Fill();

  pcalo_x.clear();
  pcalo_y.clear();
  pcalo_z.clear();
  pcalo_e.clear();
  pcalo_detid.clear();
  pcalo_layer.clear();
  pcalo_evt.clear();
  eventnr=eventnr+1;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void PCaloHitAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void PCaloHitAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PCaloHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(PCaloHitAnalyzer);

