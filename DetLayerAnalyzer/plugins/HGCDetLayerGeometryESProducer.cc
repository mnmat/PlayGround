/** \file
 *
 *  \author N. Amapane - CERN
 *
 *  \modified by R. Radogna & C. Calabria & A. Sharma
 *  \modified by D. Nash
 */

///*

#include "PlayGround/DetLayerAnalyzer/plugins/HGCDetLayerGeometry.h"
#include "PlayGround/DetLayerAnalyzer/plugins/DetHGCTracker.h"
#include "PlayGround/DetLayerAnalyzer/plugins/HGCGeomSearchDet.h"
#include "PlayGround/DetLayerAnalyzer/plugins/HGCDetLayer.h"

#include "Geometry/CommonTopologies/interface/HGCDiskGeomDet.h"
#include "Geometry/Records/interface/CaloRecoGeometryRecord.h"
//#include "CaloRecoGeometryRecord.h"


#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <memory>

class HGCDetLayerGeometryESProducer : public edm::ESProducer {
public:
  /// Constructor
  HGCDetLayerGeometryESProducer(const edm::ParameterSet& p);

  /// Produce MuonDeLayerGeometry.
  std::unique_ptr<HGCDetLayerGeometry> produce(const CaloRecoGeometryRecord& record);
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  std::pair<std::vector<DetLayer*>, std::vector<DetLayer*>> buildLayers(const DetHGCTracker* hgcTracker_) const;
  HGCDetLayer buildLayer(const HGCDiskGeomDet* det) const;

  edm::ESGetToken<DetHGCTracker,CaloGeometryRecord> DetHGCTrackerToken_;
};

using namespace edm;

HGCDetLayerGeometryESProducer::HGCDetLayerGeometryESProducer(const edm::ParameterSet& p) {
  auto cc = setWhatProduced(this);
  DetHGCTrackerToken_ = cc.consumesFrom<DetHGCTracker,CaloGeometryRecord>();
}

std::pair<std::vector<DetLayer*>, std::vector<DetLayer*>> HGCDetLayerGeometryESProducer::buildLayers(const DetHGCTracker* hgcTracker_) const {

  
  int zside = 1;
  auto geomdets = hgcTracker_->getDiskLayers(zside);

  std::vector<HGCDetLayer*> siDetLayers;
  std::vector<HGCDetLayer*> scDetLayers;

  for (const auto& det: geomdets){
    const HGCDiskGeomDet* siDiskGeomDet = det->second; 
    const HGCDiskGeomDet* scDiskGeomDet = det->first; 
    
    std::cout << siDiskGeomDet->subdet() << ", " << siDiskGeomDet->zside() << ", " << siDiskGeomDet->layer() << std::endl; 
    std::cout << siDiskGeomDet->surface().position().z() << std::endl;

    siDetLayers.push_back(new HGCDetLayer(buildLayer(siDiskGeomDet)));    

    if (scDiskGeomDet){
      scDetLayers.push_back(new HGCDetLayer(buildLayer(scDiskGeomDet)));
    }
  }
  std::vector<DetLayer*> siLayer;
  std::vector<DetLayer*> scLayer;

  for (std::vector<HGCDetLayer*>::const_iterator it = siDetLayers.begin(); it != siDetLayers.end(); it++){
    siLayer.push_back((DetLayer*)(*it));
  }
  for (std::vector<HGCDetLayer*>::const_iterator it = scDetLayers.begin(); it != scDetLayers.end(); it++){
    scLayer.push_back((DetLayer*)(*it));
  }

  std::pair<std::vector<DetLayer*>,std::vector<DetLayer*>> ret{siLayer,scLayer};
  return ret;
}

HGCDetLayer HGCDetLayerGeometryESProducer::buildLayer(const HGCDiskGeomDet* det) const{
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

std::unique_ptr<HGCDetLayerGeometry> HGCDetLayerGeometryESProducer::produce(const CaloRecoGeometryRecord& iRecord) {
  //const std::string metname = "Muon|RecoMuon|RecoMuonDetLayers|HGCDetLayerGeometryESProducer";
  auto hgcDetLayerGeometry = std::make_unique<HGCDetLayerGeometry>();
  
  const DetHGCTracker* hgcTracker_ = &(iRecord.get(DetHGCTrackerToken_));

  // Build HGCDetLayers
  std::pair<std::vector<DetLayer*>, std::vector<DetLayer*>> layers = buildLayers(hgcTracker_);

  std::cout << (layers.first).size() << std::endl;
  std::cout << (layers.second).size() << std::endl;

  std::cout << "Fill Layers" << std::endl;
  
  
  int i = 5;
  int* stack_pointer = &i;
  int* heap_pointer = new int(5);

  hgcDetLayerGeometry->addStackPointer(stack_pointer);
  hgcDetLayerGeometry->addHeapPointer(heap_pointer);
  hgcDetLayerGeometry->addSiLayers(layers.first);
  hgcDetLayerGeometry->addScintillatorLayers(layers.second);
  
  std::cout << "Done with ESProducer" << std::endl;

  // 
  // // Build DT layers
  // if (auto dt = record.getHandle(dtToken_)) {
  //   hgcDetLayerGeometry->addDTLayers(MuonDTDetLayerGeometryBuilder::buildLayers(*dt));
  // } else {
  //   LogInfo(metname) << "No DT geometry is available.";
  // }

  // // Build CSC layers
  // if (auto csc = record.getHandle(cscToken_)) {
  //   HGCDetLayerGeometry->addCSCLayers(MuonCSCDetLayerGeometryBuilder::buildLayers(*csc));
  // } else {
  //   LogInfo(metname) << "No CSC geometry is available.";
  // }

  // // Build GEM layers
  // if (auto gem = record.getHandle(gemToken_)) {
  //   HGCDetLayerGeometry->addGEMLayers(MuonGEMDetLayerGeometryBuilder::buildEndcapLayers(*gem));
  // } else {
  //   LogInfo(metname) << "No GEM geometry is available.";
  // }

  // // Build ME0 layers
  // if (auto me0 = record.getHandle(me0Token_)) {
  //   HGCDetLayerGeometry->addME0Layers(MuonME0DetLayerGeometryBuilder::buildEndcapLayers(*me0));
  // } else {
  //   LogDebug(metname) << "No ME0 geometry is available.";
  // }

  // // Build RPC layers
  // if (auto rpc = record.getHandle(rpcToken_)) {
  //   HGCDetLayerGeometry->addRPCLayers(MuonRPCDetLayerGeometryBuilder::buildBarrelLayers(*rpc),
  //                                      MuonRPCDetLayerGeometryBuilder::buildEndcapLayers(*rpc));
  // } else {
  //   LogInfo(metname) << "No RPC geometry is available.";
  // }

  // // Sort layers properly
  // HGCDetLayerGeometry->sortLayers();

  // 

  return hgcDetLayerGeometry;
}

void HGCDetLayerGeometryESProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  //no parameters are used
  descriptions.addDefault(desc);
}

DEFINE_FWK_EVENTSETUP_MODULE(HGCDetLayerGeometryESProducer);

//*/