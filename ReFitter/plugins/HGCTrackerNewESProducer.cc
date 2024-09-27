// -*- C++ -*-
//
// Package:    PlayGround/HGCTrackerNewESProducer
// Class:      HGCTrackerNewESProducer
//
/**\class HGCTrackerNewESProducer

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Fri, 19 May 2023 09:58:45 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "Geometry/CommonTopologies/interface/HGCDiskGeomDet.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "PlayGround/ReFitter/plugins/DummyHGCTracker.h"

using namespace edm;

class HGCTrackerNewESProducer : public edm::ESProducer {
public:
  HGCTrackerNewESProducer(const edm::ParameterSet&);
  ~HGCTrackerNewESProducer() override;

  using ReturnType = std::unique_ptr<DummyHGCTracker>;

  ReturnType produce(const CaloGeometryRecord&);

private:
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  std::vector<HGCDiskGeomDet *> disksPos_, disksNeg_;
  hgcal::RecHitTools rhtools_;

  std::vector<double> radlen_, xi_;
};

HGCTrackerNewESProducer::HGCTrackerNewESProducer(const edm::ParameterSet& iConfig) {
  //the following line is needed to tell the framework what
  // data is being produced
  auto cc = setWhatProduced(this);

  radlen_ = iConfig.getParameter<std::vector<double>>("radlen");
  xi_ = iConfig.getParameter<std::vector<double>>("xi");
  caloGeomToken_ = cc.consumes<CaloGeometry>();
  //now do what ever other initialization is needed
}

HGCTrackerNewESProducer::~HGCTrackerNewESProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

HGCTrackerNewESProducer::ReturnType HGCTrackerNewESProducer::produce(const CaloGeometryRecord& iRecord) {
  const CaloGeometry* geom = &(iRecord.get(caloGeomToken_));
  rhtools_.setGeometry(*geom);

  auto product = std::make_unique<DummyHGCTracker>(DummyHGCTracker(geom, radlen_, xi_));
  return product;
}

//define this as a plug-in
DEFINE_FWK_EVENTSETUP_MODULE(HGCTrackerNewESProducer);