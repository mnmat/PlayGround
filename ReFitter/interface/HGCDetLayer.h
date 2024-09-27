/*

#ifndef HGC_DetLayers_HGCForwardLayer_H
#define HGC_DetLayers_HGCForwardLayer_H

#include "TrackingTools/DetLayers/interface/ForwardDetLayer.h"

class HGCDisk;
class GeomDet;

class HGCForwardLayer : public ForwardDetLayer {
public:
  /// Constructor, takes ownership of pointers
  HGCForwardLayer(const std::vector<const HGCDisk*>& sectors);

  ~HGCForwardLayer() override;

  // GeometricSearchDet interface

  const std::vector<const GeomDet*>& basicComponents() const override { return theBasicComps; }

  const std::vector<const GeometricSearchDet*>& components() const override;

  std::vector<DetWithState> compatibleDets(const TrajectoryStateOnSurface& startingState,
                                           const Propagator& prop,
                                           const MeasurementEstimator& est) const override;

  std::vector<DetGroup> groupedCompatibleDets(const TrajectoryStateOnSurface& startingState,
                                              const Propagator& prop,
                                              const MeasurementEstimator& est) const override;

  // DetLayer interface
  SubDetector subDetector() const override;

  // Extension of the interface

  /// Return the vector of sectors
  virtual const std::vector<const HGCDisk*>& sectors() const { return theSectors; }

private:
  std::vector<const HGCDisk*> theSectors;
  std::vector<const GeometricSearchDet*> theComponents;  // duplication of the above
  std::vector<const GeomDet*> theBasicComps;             // All chambers
};
#endif

*/
