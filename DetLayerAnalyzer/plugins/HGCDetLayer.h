#ifndef HGC_DetLayers_MTDSectorForwardDoubleLayer_H
#define HGC_DetLayers_MTDSectorForwardDoubleLayer_H

#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"
#include "Utilities/BinningTools/interface/BaseBinFinder.h"

class GeomDet;
class HGCGeomSearchDet;

class HGCDetLayer : public ForwardDetLayer {
    public:
        HGCDetLayer(const std::vector<const HGCGeomSearchDet*>& sectors);

        ~HGCDetLayer() override;

        // GeometricSearchDet interface

        const std::vector<const GeomDet*>& basicComponents() const override { return theBasicComponents; }

        const std::vector<const GeometricSearchDet*>& components() const override { return theComponents; }

        // tries closest layer first
        //std::pair<bool, TrajectoryStateOnSurface> compatible(const TrajectoryStateOnSurface&,
        //                                                    const Propagator&,
        //                                                    const MeasurementEstimator&) const override;

        std::vector<DetWithState> compatibleDets(const TrajectoryStateOnSurface& startingState,
                                                const Propagator& prop,
                                                const MeasurementEstimator& est) const override;

        std::vector<DetGroup> groupedCompatibleDets(const TrajectoryStateOnSurface& startingState,
                                                    const Propagator& prop,
                                                    const MeasurementEstimator& est) const override;


        // DetLayer interface
        SubDetector subDetector() const override;


    private:
        std::vector<const HGCGeomSearchDet*> theSectors;
        std::vector<const GeometricSearchDet*> theComponents;  // duplication of the above
        std::vector<const GeomDet*> theBasicComponents;        // All chambers

};
#endif 