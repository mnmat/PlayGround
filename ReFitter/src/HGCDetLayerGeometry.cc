/*
#include <PlayGround/ReFitter/interface/HGCDetLayerGeometry.h>
//#include <DataFormats/ForwardDetid/interface/HGCalDetId.h>

#include <FWCore/Utilities/interface/Exception.h>
#include <TrackingTools/DetLayers/interface/DetLayer.h>



HGCDetLayerGeometry::HGCDetLayerGeometry(){};

HGCDetLayerGeometry::~HGCDetLayerGeometry(){};
*/
/*

void HGCDetLayerGeometry::addEELayers(const std::pair<std::vector<DetLayer*>, std::vector<DetLayer*> >& ee_layers){
  for (auto const it : ee_layers.first) {
    eeLayers_fw.push_back(it);
    detLayersMap[makeDetLayerId(it)] = it;
  }

  for (auto const it : ee_layers.second) {
    eeLayers_bw.push_back(it);
    detLayersMap[makeDetLayerId(it)] = it;
  }
}

void HGCDetLayerGeometry::addFHLayers(const std::pair<std::vector<DetLayer*>, std::vector<DetLayer*> >& fh_layers){
  for (auto const it : fh_layers.first) {
    fhLayers_fw.push_back(it);
    detLayersMap[makeDetLayerId(it)] = it;
  }

  for (auto const it : fh_layers.second) {
    fhLayers_bw.push_back(it);
    detLayersMap[makeDetLayerId(it)] = it;
  }
}


void HGCDetLayerGeometry::addBHLayers(const std::pair<std::vector<DetLayer*>, std::vector<DetLayer*> >& bh_layers){
  for (auto const it : bh_layers.first) {
    bhLayers_fw.push_back(it);
    detLayersMap[makeDetLayerId(it)] = it;
  }

  for (auto const it : bh_layers.second) {
    bhLayers_bw.push_back(it);
    detLayersMap[makeDetLayerId(it)] = it;
  }
}



DetId HGCDetLayerGeometry::makeDetLayerId(const DetLayer* detLayer) const {
  if (detLayer->subDetector() == 1) {
    return HGCalDetId(detLayer->basicComponents().front()->geographicalId().rawId());
  } else if (detLayer->subDetector() == 2) {
    return HGCalDetId(detLayer->basicComponents().front()->geographicalId().rawId());
  } else if (detLayer->subDetector() == 3) {
    return HGCalDetId(detLayer->basicComponents().front()->geographicalId().rawId());
  } else
    throw cms::Exception("InvalidModuleIdentification");  // << detLayer->module();
}
*/
