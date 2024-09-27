
#ifndef DetLayers_HGCDetLayerGeometry_h
#define DetLayers_HGCDetLayerGeometry_h

#include <TrackingTools/DetLayers/interface/DetLayerGeometry.h>
#include <DataFormats/DetId/interface/DetId.h>

#include <map>
class DetLayer;

class HGCDetLayerGeometry : public DetLayerGeometry{
    public:
        HGCDetLayerGeometry();

        ~HGCDetLayerGeometry();

        const std::vector<const DetLayer*>& allLayers() const {return allDetLayers;};

        const DetLayer* idToLayer(const DetId& detId) const override;

        void addSiLayers(const std::vector<DetLayer*>& detlayers);
        void addScintillatorLayers(const std::vector<DetLayer*>& detlayers);

        void addStackPointer(int* pointer){
            stack_integer_pointer = pointer;
        };

        void addHeapPointer(int* pointer){
            heap_integer_pointer = pointer;
        };

    private:
        DetId makeSiDetLayerId(const DetLayer* detLayer, const int layer) const;
        DetId makeScintillatorDetLayerId(const DetLayer* detLayer, const int layer) const;

        std::map<DetId, const DetLayer*> detLayersMap;

        int* stack_integer_pointer;
        int* heap_integer_pointer;

        std::vector<const DetLayer*> allDetLayers;
        std::vector<const DetLayer*> siDetLayers;
        std::vector<const DetLayer*> scDetLayers;
};

#endif
