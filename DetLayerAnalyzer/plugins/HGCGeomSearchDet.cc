#include "PlayGround/DetLayerAnalyzer/plugins/HGCGeomSearchDet.h"
#include "PlayGround/DetLayerAnalyzer/plugins/HGCGeomSearchDetBuilderFromDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/DetLayers/interface/MeasurementEstimator.h"

#include <iostream>
#include <iomanip>
#include <vector>

HGCGeomSearchDet::HGCGeomSearchDet(std::vector<const GeomDet*>::const_iterator first,
                           std::vector<const GeomDet*>::const_iterator last)
    : GeometricSearchDet(false), theDets(first, last) {
  std::cout << "HGCGeomSearchDet: Initialize w. Iterator" << std::endl;
  init();
}

HGCGeomSearchDet::HGCGeomSearchDet(const std::vector<const GeomDet*>& vdets)
    : GeometricSearchDet(false), theDets(vdets) {
  std::cout << "HGCGeomSearchDet: Initialize" << std::endl;
  init();
}

void HGCGeomSearchDet::init() {
  // Add here the sector build based on a collection of GeomDets, mimic what done in ForwardDetRingOneZ
  // using the code from tracker BladeShapeBuilderFromDet
  // simple initial version, no sorting for the time being
  setDisk(HGCGeomSearchDetBuilderFromDet()(theDets));
}

const std::vector<const GeometricSearchDet*>& HGCGeomSearchDet::components() const {
  // FIXME dummy impl.
  edm::LogError("HGCDetLayers") << "temporary dummy implementation of HGCGeomSearchDet::components()!!";
  static const std::vector<const GeometricSearchDet*> result;
  return result;
}


std::pair<bool, TrajectoryStateOnSurface> HGCGeomSearchDet::compatible(const TrajectoryStateOnSurface& ts,
                                                              const Propagator& prop,
                                                              const MeasurementEstimator& est) const {
  TrajectoryStateOnSurface ms = prop.propagate(ts, specificSurface());

#ifdef EDM_ML_DEBUG
  LogTrace("HGCDetLayers") << "HGCGeomSearchDet::compatible, sector: \n"
                           << (*this) << "\n  TS at Z,R,phi: " << std::fixed << std::setw(14) << ts.globalPosition().z()
                           << " , " << std::setw(14) << ts.globalPosition().perp() << " , " << std::setw(14)
                           << ts.globalPosition().phi();
  if (ms.isValid()) {
    LogTrace("HGCDetLayers") << " DEST at Z,R,phi: " << std::fixed << std::setw(14) << ms.globalPosition().z() << " , "
                             << std::setw(14) << ms.globalPosition().perp() << " , " << std::setw(14)
                             << ms.globalPosition().phi() << " local Z: " << std::setw(14) << ms.localPosition().z();
  } else {
    LogTrace("HGCDetLayers") << " DEST: not valid";
  }
#endif

  return std::make_pair(ms.isValid() and est.estimate(ms, specificSurface()) != 0, ms);
}

std::vector<GeometricSearchDet::DetWithState> HGCGeomSearchDet::compatibleDets(const TrajectoryStateOnSurface&,
                                   const Propagator&,
                                   const MeasurementEstimator&) const {
  edm::LogError("HGCDetLayers") << "At the moment not a real implementation";
  std::vector<GeometricSearchDet::DetWithState> result;  
  return result;
}

void HGCGeomSearchDet::compatibleDetsV(const TrajectoryStateOnSurface&,
                                   const Propagator&,
                                   const MeasurementEstimator&,
                                   std::vector<GeometricSearchDet::DetWithState>&) const {
  edm::LogError("HGCDetLayers") << "At the moment not a real implementation";
}

std::vector<DetGroup> HGCGeomSearchDet::groupedCompatibleDets(const TrajectoryStateOnSurface& startingState,
                                                     const Propagator& prop,
                                                     const MeasurementEstimator& est) const {
  // FIXME should be implemented to allow returning  overlapping chambers
  // as separate groups!
  edm::LogInfo("HGCDetLayers") << "dummy implementation of MTDDetSector::groupedCompatibleDets()";
  std::vector<DetGroup> result;
  return result;
}