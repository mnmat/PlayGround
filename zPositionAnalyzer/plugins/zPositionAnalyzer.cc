// -*- C++ -*-
//
// Package:    PlayGround/zPositionAnalyzer
// Class:      zPositionAnalyzer
//
/**\class zPositionAnalyzer zPositionAnalyzer.cc PlayGround/zPositionAnalyzer/plugins/zPositionAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Sat, 02 Dec 2023 17:48:17 GMT
//
//

// system include files
#include <memory>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// topology and geometry
#include "CondFormats/GeometryObjects/interface/PTrackerParameters.h"
#include "Geometry/Records/interface/PTrackerParametersRcd.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeomBuilderFromGeometricDet.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include "Alignment/TrackerAlignment/interface/AlignableTracker.h"
#include "Alignment/TrackerAlignment/interface/AlignableSiStripDet.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class zPositionAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit zPositionAnalyzer(const edm::ParameterSet&);
  ~zPositionAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tTopoToken_;
  const edm::ESGetToken<GeometricDet, IdealGeometryRecord> geomDetToken_;
  const edm::ESGetToken<PTrackerParameters, PTrackerParametersRcd> ptpToken_;
  const edm::ESGetToken<GEMGeometry,MuonGeometryRecord> mgeomToken_;

  // topology and geometry
  const TrackerTopology* trackerTopology;
  const TrackerGeometry* trackerGeometry;
  const GEMGeometry* mGeometry;

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
zPositionAnalyzer::zPositionAnalyzer(const edm::ParameterSet& iConfig)
    : tTopoToken_(esConsumes()),
      geomDetToken_(esConsumes()),
      ptpToken_(esConsumes()),
      mgeomToken_(esConsumes()),
      trackerTopology(0),
      trackerGeometry(0){
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

zPositionAnalyzer::~zPositionAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void zPositionAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Setup Geometry
  trackerTopology = &iSetup.getData(tTopoToken_);
  mGeometry = &iSetup.getData(mgeomToken_);
  edm::ESHandle<GeometricDet> geometricDet = iSetup.getHandle(geomDetToken_);
  edm::ESHandle<PTrackerParameters> trackerParams = iSetup.getHandle(ptpToken_);

  TrackerGeomBuilderFromGeometricDet trackerGeometryBuilder;
  trackerGeometry = trackerGeometryBuilder.build(&(*geometricDet), *trackerParams, trackerTopology);
  

  std::set<float> name; 
  float min = 999;
  float max = 0;
  float r;
  // for (auto& det : trackerGeometry->dets()) {
  //   // std::cout << "There is a Detid" << std::endl;
  //   auto detId = det->geographicalId();
  //   if (detId.det() == DetId::Muon){
  //       std::cout << detId.subdetId() << std::endl;
  //       if (detId.subdetId() == StripSubdetector::TEC){
  //       auto layer = trackerTopology->layer(detId);
  //       const GlobalPoint& globalPosition = trackerGeometry->idToDet(detId)->position();
  //       r = std::sqrt(globalPosition.x()*globalPosition.x() + globalPosition.y()*globalPosition.y());
  //       if (layer == 5){
  //         std::cout << layer << ", " <<  globalPosition.z() << std::endl;
  //         name.emplace(globalPosition.z());
  //         if (r > max){
  //           max = r;
  //         }
  //         else if (r < min){
  //           min = r;
  //         }
  //       }
  //     }
  //   }
  // }

  for (auto& det: mGeometry->dets()){
    auto detid = det->geographicalId();
    auto pos = det->position();
    auto mid = static_cast<GEMDetId>(detid);
    if ((mid.layer() == 0) & (mid.ring()==1) & (mid.station()==1) & (mid.region()==1)){
      //std::cout << mid.chamberId() << ", " << mid.layerId() << std::endl;
      //std::cout << mid.ring() << ", " << mid.station() << ", " << mid.region() << std::endl;
      auto r = std::sqrt(pos.x()*pos.x() + pos.y()*pos.y());
      std::cout << pos.z() << ", " << r << std::endl; 
    }
  }

  // int subdet = 4;
  // int layer = 5;
  // for (auto &pos : name){
  //   std::cout << subdet << "," << layer << "," << pos << "," << max << "," << min << std::endl;
  // }

  // std::cout << "Number of DetIDs" << counter << std::endl;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void zPositionAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void zPositionAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void zPositionAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(zPositionAnalyzer);