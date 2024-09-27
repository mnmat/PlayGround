//#define EDM_ML_DEBUG

#include "PlayGround/DetLayerAnalyzer/plugins/HGCDetLayer.h"
#include "PlayGround/DetLayerAnalyzer/plugins/HGCGeomSearchDet.h"
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <DataFormats/GeometrySurface/interface/DiskSectorBounds.h>
#include <TrackingTools/GeomPropagators/interface/Propagator.h>
#include <TrackingTools/DetLayers/interface/MeasurementEstimator.h>

#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

HGCDetLayer::HGCDetLayer(const std::vector<const HGCGeomSearchDet*>& sectors)
    : ForwardDetLayer(false), theSectors(sectors), theComponents(theSectors.begin(), theSectors.end()) {
  // Initial values for R, Z and Phi bounds

  std::cout << "Entered constructor" << std::endl;
  std::cout << sectors.front()->basicComponents().front()->position().z() << std::endl;

  //float theRmin = sectors.front()->rmin();
  float theRmin = sectors.front()->basicComponents().front()->position().perp();
  float theRmax = theRmin;
  float theZmin = sectors.front()->position().z();
  float theZmax = theZmin;

  std::cout << "Done setting variables" << std::endl;


  // Cache chamber pointers (the basic components_)
  // and find extension in R and Z
  for (const auto& isect : sectors) {
    std::vector<const GeomDet*> tmp2 = isect->basicComponents();
    theBasicComponents.insert(theBasicComponents.end(), tmp2.begin(), tmp2.end());

    theRmin = min(theRmin, isect->specificSurface().innerRadius());
    theRmax = max(theRmax, isect->specificSurface().outerRadius());
    float halfThick = isect->surface().bounds().thickness() / 2.;
    float zCenter = isect->surface().position().z();
    theZmin = min(theZmin, zCenter - halfThick);
    theZmax = max(theZmax, zCenter + halfThick);
  }

  std::cout << "Done with loop" << std::endl;

  // Build surface

  float zPos = (theZmax + theZmin) / 2.;
  PositionType pos(0., 0., zPos);
  RotationType rot;

  setSurface(new BoundDisk(pos, rot, new SimpleDiskBounds(theRmin, theRmax, theZmin - zPos, theZmax - zPos)));

  std::cout << "Done Building surface" << std::endl;

  LogTrace("HGCDetLayers") << "Constructing HGCDetLayer: " << std::fixed << std::setw(14)
                           << basicComponents().size() << " Dets, " << std::setw(14) << theSectors.size()
                           << " Sectors, "
                           << " Z: " << std::setw(14) << specificSurface().position().z() << " R1: " << std::setw(14)
                           << specificSurface().innerRadius() << " R2: " << std::setw(14)
                           << specificSurface().outerRadius();
}

HGCDetLayer::~HGCDetLayer() {
  //for (auto& i : theSectors) {
  //   delete i;
  //}
}


std::vector<GeometricSearchDet::DetWithState> HGCDetLayer::compatibleDets(
    const TrajectoryStateOnSurface& startingState, const Propagator& prop, const MeasurementEstimator& est) const {
  vector<DetWithState> result;
  edm::LogInfo("HGCDetLayers") << "dummy implementation of HGCDetLayer::groupedCompatibleDets()";
  return result;
}

std::vector<DetGroup> HGCDetLayer::groupedCompatibleDets(const TrajectoryStateOnSurface& startingState,
                                                              const Propagator& prop,
                                                              const MeasurementEstimator& est) const {
  // FIXME should return only 1 group
  edm::LogInfo("HGCDetLayers") << "dummy implementation of HGCDetLayer::groupedCompatibleDets()";
  return vector<DetGroup>();
}

GeomDetEnumerators::SubDetector HGCDetLayer::subDetector() const {
  return theBasicComponents.front()->subDetector();
}
