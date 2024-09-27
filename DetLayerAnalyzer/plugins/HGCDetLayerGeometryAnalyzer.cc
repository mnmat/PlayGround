// -*- C++ -*-
//
// Package:    PlayGround/HGCDetLayerGeometryAnalyzer
// Class:      HGCDetLayerGeometryAnalyzer
//
/**\class HGCDetLayerGeometryAnalyzer HGCDetLayerGeometryAnalyzer.cc PlayGround/HGCDetLayerGeometryAnalyzer/plugins/HGCDetLayerGeometryAnalyzer.cc

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

#include "PlayGround/DetLayerAnalyzer/plugins/HGCDetLayerGeometry.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class HGCDetLayerGeometryAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit HGCDetLayerGeometryAnalyzer(const edm::ParameterSet&);
  ~HGCDetLayerGeometryAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  HGCDetLayer buildLayer(const HGCDiskGeomDet* det);

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bfieldtoken_;
  edm::ESGetToken<DetHGCTracker,CaloGeometryRecord> DetHGCTrackerToken_;

  hgcal::RecHitTools rhtools_;
  const DetHGCTracker* hgcTracker_;

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
HGCDetLayerGeometryAnalyzer::HGCDetLayerGeometryAnalyzer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      bfieldtoken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
      DetHGCTrackerToken_(esConsumes<DetHGCTracker, CaloGeometryRecord>()) {
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

HGCDetLayerGeometryAnalyzer::~HGCDetLayerGeometryAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

HGCDetLayer HGCDetLayerGeometryAnalyzer::buildLayer(const HGCDiskGeomDet* det){
    // Create HGCGeomSearchDet object
    std::vector<const GeomDet*> gdets;
    gdets.push_back(det);
    const std::vector<const GeomDet*> constVectorRef(gdets);

    std::cout << "z GeomDet: " << constVectorRef.front()->surface().position().z() << std::endl;
    const HGCGeomSearchDet hgcsdet(constVectorRef);
    std::cout << "Done building the HGCGeomSearchDet" << std::endl;

    std::cout << "HGCGeomSearchDet Size of basic components: " << hgcsdet.basicComponents().size() << std::endl;

    // Create HGCDetLayer
    std::vector<const HGCGeomSearchDet*> hgcsdets;
    hgcsdets.push_back(&hgcsdet);
    return HGCDetLayer(hgcsdets);
}
 

// ------------ method called for each event  ------------
void HGCDetLayerGeometryAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;


  const CaloGeometry* geom = &iSetup.getData(caloGeomToken_);
  rhtools_.setGeometry(*geom);

  hgcTracker_ = &iSetup.getData(DetHGCTrackerToken_);

  edm::ESHandle<MagneticField> bfield_;
  bfield_ = iSetup.getHandle(bfieldtoken_);
  auto magneticField = bfield_.product();

  for (const auto& track : iEvent.get(tracksToken_)) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();
  }

  int zside = 1;
  auto geomdets = hgcTracker_->getDiskLayers(zside);

  std::vector<HGCDetLayer> siDetLayers;
  std::vector<HGCDetLayer> scDetLayers;

  int count = 0;

  std::cout << geomdets.size() << std::endl;

  for (const auto& det: geomdets){
      const HGCDiskGeomDet* siDiskGeomDet = det->second; 
      const HGCDiskGeomDet* scDiskGeomDet = det->first; 
      
      std::cout << siDiskGeomDet->subdet() << ", " << siDiskGeomDet->zside() << ", " << siDiskGeomDet->layer() << std::endl; 
      std::cout << siDiskGeomDet->surface().position().z() << std::endl;

      siDetLayers.emplace_back(buildLayer(siDiskGeomDet));    

      if (scDiskGeomDet){
        std::cout << "Scintillator" << std::endl;
        scDetLayers.emplace_back(buildLayer(scDiskGeomDet));
      }
      std::cout << "Layer: " << count << std::endl;
    count++;
    std::cout << "End Loop" << std::endl;
  }
  std::cout << "Done with Loop" << std::endl;

  std::vector<DetLayer*> siLayer;
  std::vector<DetLayer*> scLayer;


  std::cout << "Make DetLayers" << std::endl;
  for (std::vector<HGCDetLayer>::const_iterator it = siDetLayers.begin(); it != siDetLayers.end(); it++){
    siLayer.push_back((DetLayer*)(&(*it)));
  }

  for (std::vector<HGCDetLayer>::const_iterator it = scDetLayers.begin(); it != scDetLayers.end(); it++){
    scLayer.push_back((DetLayer*)(&(*it)));
  }
  std::cout << "Done Making DetLayers" << std::endl;

  
  HGCDetLayerGeometry hgcgeom = HGCDetLayerGeometry();
  hgcgeom.addSiLayers(siLayer);
  hgcgeom.addScintillatorLayers(scLayer);

  std::cout << siDetLayers.size() << std::endl;
  std::cout << scDetLayers.size() << std::endl;

  
//   auto dets = geomdets.front();

//   // Make Structure based on Search Det
//   // MTDSectorForwardLayer
//   // MTDDiskSectorBuilderFromDet.cc

//   std::cout << "HGCDetLayerGeometryAnalyzer: Create HGCGeomSearchDet" << std::endl; 

//   std::vector<const GeomDet*> input;
//   input.push_back(static_cast<const GeomDet*>(dets->second));

//   const std::vector<const GeomDet*>& constVectorRef = input;
//   const HGCGeomSearchDet* test = new HGCGeomSearchDet(constVectorRef);

//   std::cout << "HGCDetLayerGeometryAnalyzer: Done creating HGCGeomSearchDet" << std::endl; 

//   // Build DetLayer from HGCGeomSearchDet
//   std::vector<const HGCGeomSearchDet*> result;
//   result.push_back(test);

//   auto a = test->basicComponents().front();
//   std::cout << typeid(a).name() << "," << typeid(dets).name() << std::endl;
//   a->position();
//   [[maybe_unused]] auto test2 = HGCDetLayer(result);

//   std::cout << test->basicComponents().front()->geographicalId().rawId() << std::endl; 
//   std::cout << test->basicComponents().front()->subDetector() << std::endl; 
//   std::cout << a->position().z() << std::endl;

//   std::cout << dets->second->geographicalId().rawId() << std::endl;
  

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void HGCDetLayerGeometryAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void HGCDetLayerGeometryAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCDetLayerGeometryAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(HGCDetLayerGeometryAnalyzer);
//*/