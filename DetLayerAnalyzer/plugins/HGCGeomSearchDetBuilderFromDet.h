
#ifndef RecoMTD_DetLayers_HGCGeomSearchDetBuilderFromDet_H
#define RecoMTD_DetLayers_HGCGeomSearchDetBuilderFromDet_H

#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/GeometrySurface/interface/DiskSectorBounds.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include <utility>
#include <vector>
#include <iostream>


class HGCGeomSearchDetBuilderFromDet {
public:
  BoundDisk* operator()(const std::vector<const GeomDet*>& dets) const;
};

#endif

