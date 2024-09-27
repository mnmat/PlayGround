// -*- C++ -*-
//
// Package:    PlayGround/OutputAnalyzer
// Class:      OutputAnalyzer
//
/**\class OutputAnalyzer OutputAnalyzer.cc PlayGround/OutputAnalyzer/plugins/OutputAnalyzer.cc

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

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"



//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class OutputAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit OutputAnalyzer(const edm::ParameterSet&);
  ~OutputAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  void mergeRecHitCollections(std::vector<const HGCRecHit*>& recHitCollection,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH) const;

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<std::vector<ticl::Trackster>> tracksterToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;


  std::vector<const HGCRecHit*> recHitCollection;
  hgcal::RecHitTools rhtools_;
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
OutputAnalyzer::OutputAnalyzer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
      tracksterToken_(consumes<std::vector<ticl::Trackster>>(iConfig.getUntrackedParameter<edm::InputTag>("tracksters"))),
      hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCEEInput"))),
      hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCFHInput"))),
      hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCBHInput"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()){
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

OutputAnalyzer::~OutputAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

void OutputAnalyzer::mergeRecHitCollections(std::vector<const HGCRecHit*>& recHitCollection,
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


// ------------ method called for each event  ------------
void OutputAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::cout << "Hello" << std::endl;

  using namespace edm;

  edm::Handle<HGCRecHitCollection> ee_hits;
  edm::Handle<HGCRecHitCollection> fh_hits;
  edm::Handle<HGCRecHitCollection> bh_hits;

  iEvent.getByToken(hgcalRecHitsEEToken_, ee_hits);
  iEvent.getByToken(hgcalRecHitsFHToken_, fh_hits);
  iEvent.getByToken(hgcalRecHitsBHToken_, bh_hits);
  mergeRecHitCollections(recHitCollection, *ee_hits, *fh_hits, *bh_hits);

  const CaloGeometry* geom = &iSetup.getData(caloGeomToken_);
  rhtools_.setGeometry(*geom);


  // ----------
  // Tracks
  // ----------

  std::cout << "---------------------" << std::endl;
  std::cout << "Tracks" << std::endl;
  std::cout << "---------------------" << std::endl;

  for (const auto& track : iEvent.get(tracksToken_)) {

    // TrackParameters
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


    // Requires TrackExtra
    
    // Parameteras
    std::cout << "OuterOk: " << track.outerOk() << std::endl;
    std::cout << "InnerOk: " << track.innerOk() << std::endl;
    std::cout << "InnerPos: " << track.innerPosition() << std::endl;
    std::cout << "OuterPos: " << track.outerPosition() << std::endl;
    std::cout << "InnerMom: " << track.innerMomentum() << std::endl;
    std::cout << "OuterMom: " << track.outerMomentum() << std::endl;    
    std::cout << "InnerStateCov: " << track.innerStateCovariance() << std::endl;
    std::cout << "OuterStateCov: " << track.outerStateCovariance() << std::endl;

    // Measurements
    std::cout << "InnerDetId: " << track.innerDetId() << std::endl;
    std::cout << "OuterDetid: " << track.outerDetId() << std::endl;
    std::cout << "Size of RecHits: " << track.recHitsSize() << std::endl;

    // Residuals
    auto residuals = track.residuals();
    std::cout << residuals.residualX(0) << std::endl;
    std::cout << residuals.residualX(1) << std::endl;


    for (auto& h : track.recHits()){
      std::cout << h->rawId() << std::endl;
    }

    std::cout << "Found: " << track.found() << std::endl;
    std::cout << "Lost: " << track.lost() << std::endl;

    // Outer position
    std::cout << "OuterEta: " << track.outerEta() << std::endl;
    std::cout << "OuterTheta: " << track.outerTheta() << std::endl;  
    std::cout << "OuterPhi: " << track.outerPhi() << std::endl; 

  }

  // ----------
  // Tracksters
  // ----------

  std::cout << "---------------------" << std::endl;
  std::cout << "Tracksters" << std::endl;
  std::cout << "---------------------" << std::endl;

  for (const auto& trackster : iEvent.get(tracksterToken_)) {
    // Vertices
    std::cout << "Size of Vertices: " << trackster.vertices().size() << std::endl;
    std::cout << "Indices of Vertices" << std::endl;
    for (auto &v: trackster.vertices()){
      auto detid = recHitCollection[v]->id();
      auto layer = rhtools_.getLayerWithOffset(detid);
      std::cout << layer << ", " <<v << std::endl;
      // Todo: check if RecHits are valid
    }
    // Edges
    std::cout << "Size of Edges: " << trackster.edges().size() << std::endl;
    std::cout << "Edges between Vertices" << std::endl;
    for (auto &e: trackster.edges()){
      std::cout << e[0] << ", " << e[1] << std::endl;
      // TODO: Check if the edges connect next layers/next valid hits
    }

    // Vertex Multiplicity
    std::cout << "Size of Vertex Multiplicity: " << trackster.vertex_multiplicity().size() << std::endl;
    std::cout << "Vertex Multiplicity values" << std::endl;
    for (auto &vm: trackster.vertex_multiplicity()){
      std::cout << vm << std::endl;
    }

    // Seed
    std::cout << "Seed: " << trackster.seedID() << std::endl;
    std::cout << "Seed: " << trackster.seedIndex() << std::endl;

    // Time
    std::cout << "Time: " << trackster.time() << std::endl;
    std::cout << "TimeError: " << trackster.timeError() << std::endl;

    // Energy
    std::cout << "Regressed: " << trackster.regressed_energy() << std::endl;
    std::cout << "Raw: " << trackster.raw_energy() << std::endl;
    std::cout << "Raw EM: " << trackster.raw_em_energy() << std::endl;
    
    // PT
    std::cout << "Raw PT: " << trackster.raw_pt() << std::endl;
    std::cout << "Raw EM PT: " << trackster.raw_em_pt() << std::endl;

    // Shape
    std::cout << "Barycenter: " << trackster.barycenter() << std::endl; 
    std::cout << "Eigenvalues: " << trackster.eigenvalues()[1] << std::endl; 
    std::cout << "Eigenvectors: " << trackster.eigenvectors()[1] << std::endl; 
    std::cout << "Sigmas: " << trackster.sigmas()[1] << std::endl; 
    std::cout << "Sigmas PCA: " << trackster.sigmasPCA()[1] << std::endl; 

    // Probabilities
    //std::cout << "ID Probabilities: " << trackster.id_probabilities() << std::endl; 
  }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void OutputAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void OutputAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void OutputAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(OutputAnalyzer);
