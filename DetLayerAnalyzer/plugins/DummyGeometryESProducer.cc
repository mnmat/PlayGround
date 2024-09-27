/** \file
 *
 *  \author N. Amapane - CERN
 *
 *  \modified by R. Radogna & C. Calabria & A. Sharma
 *  \modified by D. Nash
 */
///*




#include "PlayGround/DetLayerAnalyzer/plugins/HGCDetLayerGeometry.h"
#include "TrackingTools/RecoGeometry/interface/RecoGeometryRecord.h"
//#include "Geometry/Records/interface/CaloRecoGeometryRecord.h"
#include "CaloRecoGeometryRecord.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <memory>
#include <iostream>

class CaloRecoGeometryRecord;

class DummyGeometryESProducer : public edm::ESProducer {
public:
  /// Constructor
  DummyGeometryESProducer(const edm::ParameterSet& p);

  /// Produce MuonDeLayerGeometry.
  std::unique_ptr<DetLayerGeometry> produce(const RecoGeometryRecord& record);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  edm::ESGetToken<HGCDetLayerGeometry,CaloRecoGeometryRecord> hgcDetLayerGeometryToken_;
};

using namespace edm;

DummyGeometryESProducer::DummyGeometryESProducer(const edm::ParameterSet& p) {
  auto cc = setWhatProduced(this,"DummyGeometry");
  hgcDetLayerGeometryToken_ = cc.consumes();
}

std::unique_ptr<DetLayerGeometry> DummyGeometryESProducer::produce(const RecoGeometryRecord& iRecord) {
  //const std::string metname = "Muon|RecoMuon|RecoMuonDetLayers|DummyGeometryESProducer";
  const HGCDetLayerGeometry& hgcDetLayerGeometry = iRecord.get(hgcDetLayerGeometryToken_);
  auto dummyGeometry = std::make_unique<DetLayerGeometry>(hgcDetLayerGeometry);
  std::cout << "Done with ESProducer" << std::endl;

  return dummyGeometry;
}

void DummyGeometryESProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  //no parameters are used
  descriptions.addDefault(desc);
}

DEFINE_FWK_EVENTSETUP_MODULE(DummyGeometryESProducer);

//*/