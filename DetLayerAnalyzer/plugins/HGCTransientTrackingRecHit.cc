/** \file
 *
 */

#include "PlayGround/DetLayerAnalyzer/plugins/HGCTransientTrackingRecHit.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/AlignmentPositionError.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Utilities/interface/Exception.h"

#include <map>

typedef HGCTransientTrackingRecHit::HGCRecHitPointer HGCRecHitPointer;
typedef HGCTransientTrackingRecHit::RecHitContainer HGCRecHitContainer;

HGCTransientTrackingRecHit::HGCTransientTrackingRecHit(const GeomDet* geom, const TrackingRecHit* rh)
    : GenericTransientTrackingRecHit(*geom, *rh) {}

HGCTransientTrackingRecHit::HGCTransientTrackingRecHit(const HGCTransientTrackingRecHit& other)
    : GenericTransientTrackingRecHit(*other.det(), *(other.hit())) {}
