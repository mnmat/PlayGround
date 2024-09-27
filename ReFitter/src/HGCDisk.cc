/*

#include "PlayGround/ReFitter/interface/HGCDisk.h"
#include "PlayGround/ReFitter/interface/HGCDiskBuilderFromDet.h"



HGCDisk::HGCDisk(vector<const GeomDet*>::const_iterator first,
                           vector<const GeomDet*>::const_iterator last,
                           const MTDTopology& topo)
    : GeometricSearchDet(false), theDets(first, last), topo_(&topo) {
  init();
}

HGCDisk::HGCDisk(const vector<const GeomDet*>& vdets, const MTDTopology& topo)
    : GeometricSearchDet(false), theDets(vdets), topo_(&topo) {
  init();
}

void HGCDisk::init() {
  // Add here the sector build based on a collection of GeomDets, mimic what done in ForwardDetRingOneZ
  // using the code from tracker BladeShapeBuilderFromDet
  // simple initial version, no sorting for the time being
  setDisk(HGCDiskBuilderFromDet()(theDets));
}

const vector<const GeometricSearchDet*>& HGCDisk::components() const {
  // FIXME dummy impl.
  edm::LogError("HGCDetLayers") << "temporary dummy implementation of HGCDisk::components()!!";
  static const vector<const GeometricSearchDet*> result;
  return result;
}


pair<bool, TrajectoryStateOnSurface> HGCDisk::compatible(const TrajectoryStateOnSurface& ts,
                                                              const Propagator& prop,
                                                              const MeasurementEstimator& est) const {
  TrajectoryStateOnSurface ms = prop.propagate(ts, specificSurface());

#ifdef EDM_ML_DEBUG
  LogTrace("HGCDetLayers") << "HGCDisk::compatible, sector: \n"
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

  return make_pair(ms.isValid() and est.estimate(ms, specificSurface()) != 0, ms);
}

void HGCDisk::compatibleDets(const TrajectoryStateOnSurface&,
                                   const Propagator&,
                                   const MeasurementEstimator&) const {
  edm::LogError("HGCDetLayers") << "At the moment not a real implementation";
}

void HGCDisk::compatibleDetsV(const TrajectoryStateOnSurface&,
                                   const Propagator&,
                                   const MeasurementEstimator&,
                                   std::vector<DetWithState>&) const {
  edm::LogError("HGCDetLayers") << "At the moment not a real implementation";
}

*/
