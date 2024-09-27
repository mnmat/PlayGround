#ifndef HGCTransientTrackingRecHit_h
#define HGCTransientTrackingRecHit_h

/** \class HGCTransientTrackingRecHit
 *
 *  A TransientTrackingRecHit for HGCs.
 *
 *
 *   \author   C. Liu            Purdue University
 *
 *   \modified by C. Calabria    INFN & Universita  Bari
 */

#include "TrackingTools/TransientTrackingRecHit/interface/GenericTransientTrackingRecHit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

class HGCTransientTrackingRecHit final : public GenericTransientTrackingRecHit {
public:
  using HGCRecHitPointer = std::shared_ptr<HGCTransientTrackingRecHit>;
  using ConstHGCRecHitPointer = std::shared_ptr<HGCTransientTrackingRecHit const>;

  //  typedef ReferenceCountingPointer<HGCTransientTrackingRecHit>      HGCRecHitPointer;
  //  typedef ConstReferenceCountingPointer<HGCTransientTrackingRecHit> ConstHGCRecHitPointer;
  typedef std::vector<HGCRecHitPointer> HGCRecHitContainer;
  typedef std::vector<ConstHGCRecHitPointer> ConstHGCRecHitContainer;

  ~HGCTransientTrackingRecHit() override {}

  /// FIXME virtual ConstHGCRecHitContainer specificTransientHits() const;

  static RecHitPointer build(const GeomDet* geom, const TrackingRecHit* rh) {
    return RecHitPointer(new HGCTransientTrackingRecHit(geom, rh));
  }

  static HGCRecHitPointer specificBuild(const GeomDet* geom, const TrackingRecHit* rh) {
    LogDebug("HGC|RecoHGC|HGCDetLayerMeasurements") << "Getting specificBuild" << std::endl;
    return HGCRecHitPointer(new HGCTransientTrackingRecHit(geom, rh));
  }

  void invalidateHit();


  /// Construct from a TrackingRecHit and its GeomDet
  HGCTransientTrackingRecHit(const GeomDet* geom, const TrackingRecHit* rh);

  /// Copy ctor
  HGCTransientTrackingRecHit(const HGCTransientTrackingRecHit& other);

  HGCTransientTrackingRecHit* clone() const override { return new HGCTransientTrackingRecHit(*this); }
};
#endif