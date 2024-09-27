
///*

/*

*/
#include <DataFormats/ForwardDetId/interface/HGCSiliconDetId.h>
#include <DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h>
#include "PlayGround/DetLayerAnalyzer/plugins/HGCDetLayerGeometry.h"
#include <TrackingTools/DetLayers/interface/DetLayer.h>
#include <DataFormats/DetId/interface/DetId.h>


HGCDetLayerGeometry::HGCDetLayerGeometry() {}

HGCDetLayerGeometry::~HGCDetLayerGeometry() {
    //for (std::vector<const DetLayer*>::const_iterator it = allDetLayers.begin(); it != allDetLayers.end(); ++it){
    //    delete *it;
    //}
}

void HGCDetLayerGeometry::addSiLayers(const std::vector<DetLayer*>& detlayers){
    int layer = 0;
    for (auto const it: detlayers){
        siDetLayers.push_back(it);
        allDetLayers.push_back(it);

        detLayersMap[makeSiDetLayerId(it,layer)] = it;
        layer++;
    }
}

void HGCDetLayerGeometry::addScintillatorLayers(const std::vector<DetLayer*>& detlayers){
    int layer = 0;
    for (auto const it: detlayers){
        scDetLayers.push_back(it);
        allDetLayers.push_back(it);

        detLayersMap[makeScintillatorDetLayerId(it,layer)] = it;
        layer++;
    }
}

DetId HGCDetLayerGeometry::makeSiDetLayerId(const DetLayer* detLayer, const int layer) const{
    auto dtype = (layer < 27)? DetId::HGCalEE:DetId::HGCalHSi;// HGCalEE or HGCalHSi 
    auto test = (detLayer->basicComponents().front()->position().z() > 0)? 1:0;// z-side, +1 or -1
    std::cout << "Basic Components:" << detLayer->basicComponents().front()->position().z() << std::endl;
    int zp = 1;
    return HGCSiliconDetId(dtype,zp,2,layer,0,0,0,0);
}

DetId HGCDetLayerGeometry::makeScintillatorDetLayerId(const DetLayer* detLayer, const int layer) const{
    auto test = (detLayer->basicComponents().front()->position().z() > 0)? 1:0; // (ring < 0) ? 1 : 0; |ring| index (starting from a minimum radius depending on type)
    int ring = 1;
    return HGCScintillatorDetId(0,layer,ring,0);
}

const DetLayer* HGCDetLayerGeometry::idToLayer(const DetId& detId) const{
    HGCSiliconDetId id(detId);
    return siDetLayers[id.layer()]; 
}

#include "FWCore/Utilities/interface/typelookup.h"
TYPELOOKUP_DATA_REG(HGCDetLayerGeometry);


//*/