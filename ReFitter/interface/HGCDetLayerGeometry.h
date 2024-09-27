/*

#ifndef DetLayers_HGCDetLayerGeometry_h
#define DetLayers_HGCDetLayerGeometry_h

#include "DataFormats/DetId/interface/DetId.h"
#include "TrackingTools/DetLayers/interface/DetLayerGeometry.h"
#include <vector>
#include <map>

class DetLayer;

class HGCDetLayerGeometry : public DetLayerGeometry {
	public:
        HGCDetLayerGeometry();
        ~HGCDetLayerGeometry() override;

		const std::vector<const DetLayer*>& allEELayers() const {return eeLayers_all;};
		const std::vector<const DetLayer*>& allFHLayers() const {return fhLayers_all;};
		const std::vector<const DetLayer*>& allBHLayers() const {return bhLayers_all;};
        
        const std::vector<const DetLayer*>& forwardEELayers() const {return eeLayers_fw;};
		const std::vector<const DetLayer*>& backwardEELayers() const {return eeLayers_bw;};
		const std::vector<const DetLayer*>& forwardFHLayers() const {return fhLayers_fw;};
		const std::vector<const DetLayer*>& backwardFHLayers() const {return fhLayers_bw;};
		const std::vector<const DetLayer*>& forwardBHLayers() const {return bhLayers_fw;};
		const std::vector<const DetLayer*>& backwardBHLayers() const {return bhLayers_bw;};

		//const std::vector<const DetLayer*>& idToLayer() const override;
		
	private:
		void addEELayers(const std::pair<std::vector<DetLayer*>, std::vector<DetLayer*> >& ee_layers);
		void addFHLayers(const std::pair<std::vector<DetLayer*>, std::vector<DetLayer*> >& fh_layers);
		void addBHLayers(const std::pair<std::vector<DetLayer*>, std::vector<DetLayer*> >& bh_layers);
		//DetId makeDetLayerId(const DetLayer* detLayer) const;
		void sortLayers();
	
		
		std::vector<const DetLayer*> eeLayers_fw, eeLayers_bw, eeLayers_all;
		std::vector<const DetLayer*> fhLayers_fw, fhLayers_bw, fhLayers_all;
		std::vector<const DetLayer*> bhLayers_fw, bhLayers_bw, bhLayers_all;
	
		std::map<DetId, const DetLayer*> detLayersMap;
};
#endif

*/
