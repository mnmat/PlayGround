


#ifndef RecoHGC_DetLayers_HGCGeomSearchDet_H
#define RecoHGC_DetLayers_HGCGeomSearchDet_H

#include "TrackingTools/DetLayers/interface/GeometricSearchDet.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include <ostream>

class GeomDet;

class HGCGeomSearchDet : public GeometricSearchDet {
public:
  using GeometricSearchDet::GeometricSearchDet;

  /// Construct from iterators on GeomDet*
  HGCGeomSearchDet(std::vector<const GeomDet*>::const_iterator first,
               std::vector<const GeomDet*>::const_iterator last);

  /// Construct from a vector of GeomDet*
  HGCGeomSearchDet(const std::vector<const GeomDet*>& dets);

  ~HGCGeomSearchDet() override{};

  // GeometricSearchDet structure

  const std::vector<const GeomDet*>& basicComponents() const override { return theDets; }

  const BoundSurface& surface() const final { return *theDiskS; }

  const std::vector<const GeometricSearchDet*>& components() const override;

  std::pair<bool, TrajectoryStateOnSurface> compatible(const TrajectoryStateOnSurface& ts,
                                                       const Propagator& prop,
                                                       const MeasurementEstimator& est) const override;

  std::vector<DetWithState> compatibleDets(const TrajectoryStateOnSurface& startingState,
                                           const Propagator& prop,
                                           const MeasurementEstimator& est) const override;

  void compatibleDetsV(const TrajectoryStateOnSurface& startingState,
                       const Propagator& prop,
                       const MeasurementEstimator& est,
                       std::vector<DetWithState>& result) const override;

  std::vector<DetGroup> groupedCompatibleDets(const TrajectoryStateOnSurface& startingState,
                                              const Propagator& prop,
                                              const MeasurementEstimator& est) const override;

  // GeometricSearchDet extension

  const BoundDisk& specificSurface() const { return *theDiskS; }

  void compatibleDetsLine(const size_t idetMin,
                          std::vector<DetWithState>& result,
                          const TrajectoryStateOnSurface& tsos,
                          const Propagator& prop,
                          const MeasurementEstimator& est) const;

  size_t hshift(const uint32_t detid, const int horizontalShift) const;
  size_t vshift(const uint32_t detid, const int verticalShift, size_t& closest) const;

protected:
  void setDisk(BoundDisk* diskS) { theDiskS = diskS; }

  bool add(size_t idet,
           std::vector<DetWithState>& result,
           const TrajectoryStateOnSurface& tsos,
           const Propagator& prop,
           const MeasurementEstimator& est) const;

private:
  ReferenceCountingPointer<BoundDisk> theDiskS;
  std::vector<const GeomDet*> theDets;

  void init();
};

std::ostream& operator<<(std::ostream&, const HGCGeomSearchDet&);

#endif