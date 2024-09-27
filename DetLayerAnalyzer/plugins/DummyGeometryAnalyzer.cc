// -*- C++ -*-
//
// Package:    PlayGround/DummyGeometryAnalyzer
// Class:      DummyGeometryAnalyzer
//
/**\class DummyGeometryAnalyzer DummyGeometryAnalyzer.cc PlayGround/DummyGeometryAnalyzer/plugins/DummyGeometryAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Tue, 02 Jul 2024 11:53:40 GMT
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
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"


#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "PlayGround/DetLayerAnalyzer/plugins/DetHGCTracker.h"

#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/GeometrySurface/interface/BoundingBox.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "TrackingTools/DetLayers/interface/GeometricSearchDet.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"

#include "CondFormats/Alignment/interface/Definitions.h"
#include "PlayGround/DetLayerAnalyzer/plugins/HGCGeomSearchDet.h"
#include "PlayGround/DetLayerAnalyzer/plugins/HGCDetLayer.h"

#include "Geometry/CommonTopologies/interface/HGCDiskGeomDet.h"
#include "TrackingTools/DetLayers/interface/DetLayerGeometry.h"

#include <DataFormats/ForwardDetId/interface/HGCSiliconDetId.h>


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class DummyGeometryAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DummyGeometryAnalyzer(const edm::ParameterSet&);
  ~DummyGeometryAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  HGCDetLayer buildLayer(const HGCDiskGeomDet* det);

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::ESGetToken<DetLayerGeometry, RecoGeometryRecord> dummyGeometry_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfieldtoken_;

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
DummyGeometryAnalyzer::DummyGeometryAnalyzer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
      dummyGeometry_(esConsumes<DetLayerGeometry, RecoGeometryRecord>(edm::ESInputTag("","DummyGeometry"))),
      bfieldtoken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();

#endif
  //now do what ever initialization is needed
}

DummyGeometryAnalyzer::~DummyGeometryAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//
 

// ------------ method called for each event  ------------
void DummyGeometryAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;


  auto geom = &iSetup.getData(dummyGeometry_);



  auto id = HGCSiliconDetId(DetId::HGCalEE,1,2,1,0,0,0,0);
  [[maybe_unused]] auto detlayer = geom->idToLayer(id);




  //std::cout << geom->allLayers().size() << std::endl;

  edm::ESHandle<MagneticField> bfield_;
  bfield_ = iSetup.getHandle(bfieldtoken_);
  auto magneticField = bfield_.product();
  
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void DummyGeometryAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void DummyGeometryAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DummyGeometryAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(DummyGeometryAnalyzer);
//*/