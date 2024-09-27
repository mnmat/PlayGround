/*

#ifndef HGCDetLayerGeometryBuilder_h
#define HGCDetLayerGeometryBuilder_h

#include <Geometry/HGCalGeometry/interface/HGCalGeometry.h>
#include <vector>

class HGCDetLayerGeometryBuilder {
public:
  /// return.first=forward (+Z), return.second=backward (-Z)
  /// both vectors are sorted inside-out
  static std::pair<std::vector<DetLayer*>, std::vector<DetLayer*> > buildLayers(const HGCalGeometry& geo);

private:
  // Disable constructor - only static access is allowed.
  HGCDetLayerGeometryBuilder() {}
};
#endif

*/
