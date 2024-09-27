// -*- C++ -*-
//
// Package:    PlayGround/TrueTrackAnalyzer
// Class:      TrueTrackAnalyzer
//
/**\class TrueTrackAnalyzer TrueTrackAnalyzer.cc PlayGround/TrueTrackAnalyzer/plugins/TrueTrackAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Mon, 07 Aug 2023 21:43:34 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelectorBase.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//ROOT includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include <TVector.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include <algorithm>
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveStats.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class TrueTrackAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TrueTrackAnalyzer(const edm::ParameterSet&);
  ~TrueTrackAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
  edm::EDGetTokenT<edm::View<reco::Track>> globalMuonTracksToken_;
  edm::EDGetTokenT<edm::View<reco::Track>> standaloneMuonTracksToken_;
  std::vector<edm::InputTag> associators;
  edm::EDGetTokenT<edm::View<reco::Muon>> recoMuonToken_;


  std::vector<edm::EDGetTokenT<reco::SimToRecoCollection>> associatormapStRs;
  std::vector<edm::EDGetTokenT<reco::RecoToSimCollection>> associatormapRtSs;

  TTree* globalMuonTree;
  TTree* standaloneMuonTree;
  TTree* recoMuonTree;


  int eventidx = 0;
  std::string eta_;
  std::string energy_;

  // Global Muons
  std::vector<int> glo_idx;
  std::vector<int> glo_isSignal;
  std::vector<int> glo_isSimTrack;
  std::vector<double> glo_innerMomentum;
  std::vector<double> glo_outerMomentum;
  std::vector<double> glo_innerZ;
  std::vector<double> glo_outerZ;
  std::vector<int> glo_nrRechits;
  std::vector<int> glo_innerOk;
  std::vector<int> glo_outerOk;
  std::vector<int> glo_found;
  std::vector<int> glo_lost;



  // Inner Momentum
  // Outer Momentum
  // Number of Rechits
  // outer ok
  // inner ok
  // found
  // lost

  // Standalone Muons
  std::vector<int> sta_idx;
  std::vector<int> sta_isSignal;
  std::vector<int> sta_isSimTrack;
  std::vector<double> sta_innerMomentum;
  std::vector<double> sta_outerMomentum;
  std::vector<double> sta_innerZ;
  std::vector<double> sta_outerZ;
  std::vector<int> sta_nrRechits;
  std::vector<int> sta_innerOk;
  std::vector<int> sta_outerOk;
  std::vector<int> sta_found;
  std::vector<int> sta_lost;

  // Reco Muons
  std::vector<int> rmu_idx;
  std::vector<int> rmu_isSignal;
  std::vector<int> rmu_isSimTrack;
  std::vector<double> rmu_innerMomentum;
  std::vector<double> rmu_outerMomentum;
  std::vector<double> rmu_innerZ;
  std::vector<double> rmu_outerZ;
  std::vector<int> rmu_nrRechits;
  std::vector<int> rmu_innerOk;
  std::vector<int> rmu_outerOk;
  std::vector<int> rmu_found;
  std::vector<int> rmu_lost;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrueTrackAnalyzer::TrueTrackAnalyzer(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<edm::View<reco::Track>>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
    globalMuonTracksToken_(consumes<edm::View<reco::Track>>(iConfig.getUntrackedParameter<edm::InputTag>("muons"))),
    standaloneMuonTracksToken_(consumes<edm::View<reco::Track>>(iConfig.getUntrackedParameter<edm::InputTag>("standalone"))),
    associators(iConfig.getUntrackedParameter<std::vector<edm::InputTag>>("associators")),
    recoMuonToken_(consumes<edm::View<reco::Muon>>(iConfig.getUntrackedParameter<edm::InputTag>("recoMuons"))),
    eta_(iConfig.getParameter<std::string>("eta")),
    energy_(iConfig.getParameter<std::string>("energy")) {

  //associatormapStRs.push_back(consumes<reco::SimToRecoCollection>(associators[0]));
  associatormapRtSs.push_back(consumes<reco::RecoToSimCollection>(associators[0]));
  associatormapRtSs.push_back(consumes<reco::RecoToSimCollection>(associators[1]));
  associatormapRtSs.push_back(consumes<reco::RecoToSimCollection>(associators[2]));

  usesResource("TFileService");
  edm::Service<TFileService> file;
  globalMuonTree = file->make<TTree>("globalMuons","globalMuons"); // Replace this with an input tree value?
  standaloneMuonTree = file->make<TTree>("standaloneMuon","standaloneMuon"); // Replace this with an input tree value?
  recoMuonTree = file->make<TTree>("recoMuon","recoMuon"); // Replace this with an input tree value?

  // globalMuons
  globalMuonTree->Branch("glo_idx",&glo_idx);
  globalMuonTree->Branch("glo_isSignal", &glo_isSignal);
  globalMuonTree->Branch("glo_isSimTrack", &glo_isSimTrack);
  globalMuonTree->Branch("glo_innerMomentum", &glo_innerMomentum);
  globalMuonTree->Branch("glo_outerMomentum", &glo_outerMomentum);
  globalMuonTree->Branch("glo_innerZ", &glo_innerZ);
  globalMuonTree->Branch("glo_outerZ", &glo_outerZ);
  globalMuonTree->Branch("glo_nrRechits", &glo_nrRechits);
  globalMuonTree->Branch("glo_innerOk", &glo_innerOk);
  globalMuonTree->Branch("glo_outerOk", &glo_outerOk);
  globalMuonTree->Branch("glo_found", &glo_found);
  globalMuonTree->Branch("glo_lost", &glo_lost);

  // standaloneMuons
  standaloneMuonTree->Branch("sta_idx",&sta_idx);
  standaloneMuonTree->Branch("sta_isSignal", &sta_isSignal);
  standaloneMuonTree->Branch("sta_isSimTrack", &sta_isSimTrack);
  standaloneMuonTree->Branch("sta_innerMomentum", &sta_innerMomentum);
  standaloneMuonTree->Branch("sta_outerMomentum", &sta_outerMomentum);
  standaloneMuonTree->Branch("sta_innerZ", &sta_innerZ);
  standaloneMuonTree->Branch("sta_outerZ", &sta_outerZ);
  standaloneMuonTree->Branch("sta_nrRechits", &sta_nrRechits);
  standaloneMuonTree->Branch("sta_innerOk", &sta_innerOk);
  standaloneMuonTree->Branch("sta_outerOk", &sta_outerOk);
  standaloneMuonTree->Branch("sta_found", &sta_found);
  standaloneMuonTree->Branch("sta_lost", &sta_lost);

  // recoMuons
  recoMuonTree->Branch("rmu_idx",&rmu_idx);
  recoMuonTree->Branch("rmu_isSignal", &rmu_isSignal);
  recoMuonTree->Branch("rmu_isSimTrack", &rmu_isSimTrack);
  recoMuonTree->Branch("rmu_innerMomentum", &rmu_innerMomentum);
  recoMuonTree->Branch("rmu_outerMomentum", &rmu_outerMomentum);
  recoMuonTree->Branch("rmu_innerZ", &rmu_innerZ);
  recoMuonTree->Branch("rmu_outerZ", &rmu_outerZ);
  recoMuonTree->Branch("rmu_nrRechits", &rmu_nrRechits);
  recoMuonTree->Branch("rmu_innerOk", &rmu_innerOk);
  recoMuonTree->Branch("rmu_outerOk", &rmu_outerOk);
  recoMuonTree->Branch("rmu_found", &rmu_found);
  recoMuonTree->Branch("rmu_lost", &rmu_lost);

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

TrueTrackAnalyzer::~TrueTrackAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void TrueTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  eventidx +=1;

 // Build tsos from TrackCollection
  edm::Handle<edm::View<reco::Track>> tracks_h;
  iEvent.getByToken(tracksToken_,tracks_h);

  edm::Handle<edm::View<reco::Track>> globalMuonTracks_h;
  iEvent.getByToken(globalMuonTracksToken_,globalMuonTracks_h);

  edm::Handle<edm::View<reco::Track>> standaloneMuonTracks_h;
  iEvent.getByToken(standaloneMuonTracksToken_,standaloneMuonTracks_h);

  edm::Handle<edm::View<reco::Muon>> recoMuons_h;
  iEvent.getByToken(recoMuonToken_, recoMuons_h);

  //edm::Handle<TrackingParticleCollection> trackingParticles;
  //iEvent.getByToken(trackingParticlesToken_, trackingParticles);

  const edm::View<reco::Track> & tkx = *(tracks_h.product()); 
  const edm::View<reco::Track> & gmtkx = *(globalMuonTracks_h.product()); 
  const edm::View<reco::Track> & smtkx = *(standaloneMuonTracks_h.product()); 
  const edm::View<reco::Muon>& rmuons = *(recoMuons_h.product());


  std::cout << "# Tracks: " << tkx.size() << std::endl;
  std::cout << "# Global Muon Tracks: " << gmtkx.size() << std::endl;
  std::cout << "# Standalone Muon Tracks: " << smtkx.size() << std::endl;
  std::cout << "Associators: " << associators.size() << std::endl;
  std::cout << "# RecoMuons: " << rmuons.size() << std::endl;

  // More info on Muon Tracks

  for (auto &t: gmtkx){
    std::cout << "Inner Position: " << t.innerPosition().z() << ", Outer Position: " << t.outerPosition().z() << std::endl; 
  }

  // Load Association Map

  //reco::SimToRecoCollection const* simRecCollPFull = nullptr;
  reco::RecoToSimCollection const* recSimCollP = nullptr;
  reco::RecoToSimCollection const* recSimCollP_Muon = nullptr;
  reco::RecoToSimCollection const* recSimCollP_StandaloneMuon = nullptr;


  //Handle<reco::SimToRecoCollection> simtorecoCollectionH;
  //iEvent.getByToken(associatormapStRs[0], simtorecoCollectionH);
  //simRecCollPFull = simtorecoCollectionH.product();

  edm::Handle<reco::RecoToSimCollection> recotosimCollectionH;
  iEvent.getByToken(associatormapRtSs[0], recotosimCollectionH);
  recSimCollP = recotosimCollectionH.product();

  edm::Handle<reco::RecoToSimCollection> recotosimCollectionH_Muon;
  iEvent.getByToken(associatormapRtSs[1], recotosimCollectionH_Muon);
  recSimCollP_Muon = recotosimCollectionH_Muon.product();
  
  edm::Handle<reco::RecoToSimCollection> recotosimCollectionH_StandaloneMuon;
  iEvent.getByToken(associatormapRtSs[2], recotosimCollectionH_StandaloneMuon);
  recSimCollP_StandaloneMuon = recotosimCollectionH_StandaloneMuon.product();

  // Test ASsociation Map on Global Muon Tracks

  reco::RecoToSimCollection const& recSimColl = *recSimCollP;
  //reco::SimToRecoCollection const& simRecColl = *simRecCollP;

  reco::RecoToSimCollection const& recSimColl_Muon = *recSimCollP_Muon;
  reco::RecoToSimCollection const& recSimColl_StandaloneMuon = *recSimCollP_StandaloneMuon;


  /*
  edm::Handle<View<reco::Track>> trackCollectionH;
  iEvent.getByLabel("generaltracks", trackCollectionH);
  const View<reco::Track> &tC = *(trackCollectionH.product());
  */
 
  /* 
  tpSelector_ = TrackingParticleSelector(tpset.getParameter<double>("ptMin"),
                                    
                                         tpset.getParameter<double>("ptMax"),
                                         tpset.getParameter<double>("minRapidity"),
                                         tpset.getParameter<double>("maxRapidity"),
                                         tpset.getParameter<double>("tip"),
                                         tpset.getParameter<double>("lip"),
                                         tpset.getParameter<int>("minHit"),
                                         tpset.getParameter<bool>("signalOnly"),
                                         tpset.getParameter<bool>("intimeOnly"),
                                         tpset.getParameter<bool>("chargedOnly"),
                                         tpset.getParameter<bool>("stableOnly"),
                                         tpset.getParameter<std::vector<int> >("pdgId"));

  */
  TrackingParticleSelector tpSelector_ = TrackingParticleSelector(0, 100, 1.5, 3,120,280,0,true,false,false,false,{13});

  std::cout << "------------------ General Track -------------------" << std::endl;
  for (edm::View<reco::Track>::size_type i = 0; i < tkx.size(); ++i) {
    RefToBase<reco::Track> track(tracks_h, i);  
    [[maybe_unused]] int isSignal = 0;
    [[maybe_unused]] int isSimTrack = 0;
    if(recSimColl.find(track) != recSimColl.end()){
      isSimTrack = 1;
      std::cout << "Global Muon Track is associated with a SimTrack" << std::endl;

      std::vector<std::pair<TrackingParticleRef, double> > tps;
      tps = recSimColl[track];
      for (auto tp:tps){
        if(tpSelector_(*(tp.first))){
          isSignal=1;
          std::cout << "Global Muon Track is signal track" << std::endl;
          std::cout << "Inner Position: " << tkx[i].innerPosition().z() << ", Outer Position: " << tkx[i].outerPosition().z() << std::endl; 
          std::cout << "Signal track has index: " << i << std::endl;
          std::cout << eventidx << "," << eta_ << "," << energy_ << "," << i << std::endl;
          std::cout << "Signal has pdgId: " << tp.first->pdgId() << std::endl;
          std::cout << "GenParticleIndex: " << tp.first->g4Tracks()[0].genpartIndex() << std::endl;
  
        }
      }
    }

    double inner = std::sqrt(tkx[i].innerMomentum().Mag2());
    double outer = std::sqrt(tkx[i].outerMomentum().Mag2());
    double innerZ = tkx[i].innerPosition().z();
    double outerZ = tkx[i].outerPosition().z();

    // gen_idx.push_back(i);
    // gen_isSignal.push_back(isSignal);
    // gen_isSimTrack.push_back(isSimTrack);
    // gen_innerZ.push_back(innerZ);
    // gen_outerZ.push_back(outerZ);
    // gen_innerMomentum.push_back(inner);
    // gen_outerMomentum.push_back(outer);
    // gen_nrRechits.push_back(gmtkx[i].recHitsSize());
    // gen_innerOk.push_back(gmtkx[i].innerOk());
    // gen_outerOk.push_back(gmtkx[i].outerOk());
    // gen_found.push_back(gmtkx[i].found());
    // gen_lost.push_back(gmtkx[i].lost());
  }




  std::cout << "------------------ Global Muon Track -------------------" << std::endl;
  for (edm::View<reco::Track>::size_type i = 0; i < gmtkx.size(); ++i) {
    RefToBase<reco::Track> track(globalMuonTracks_h, i);  
    int isSignal = 0;
    int isSimTrack = 0;
    if(recSimColl_Muon.find(track) != recSimColl_Muon.end()){
      isSimTrack = 1;
      std::cout << "Global Muon Track is associated with a SimTrack" << std::endl;

      std::vector<std::pair<TrackingParticleRef, double> > tps;
      tps = recSimColl_Muon[track];
      for (auto tp:tps){
        if(tpSelector_(*(tp.first))){
          isSignal=1;
          std::cout << "Global Muon Track is signal track" << std::endl;
          std::cout << "Inner Position: " << gmtkx[i].innerPosition().z() << ", Outer Position: " << gmtkx[i].outerPosition().z() << std::endl; 
          std::cout << "Signal track has index: " << i << std::endl;
          std::cout << eventidx << "," << eta_ << "," << energy_ << "," << i << std::endl;
          std::cout << "Signal has pdgId: " << tp.first->pdgId() << std::endl;
          std::cout << "GenParticleIndex: " << tp.first->g4Tracks()[0].genpartIndex() << std::endl;
        }
      }
    }

    double inner = std::sqrt(gmtkx[i].innerMomentum().Mag2());
    double outer = std::sqrt(gmtkx[i].outerMomentum().Mag2());
    double innerZ = gmtkx[i].innerPosition().z();
    double outerZ = gmtkx[i].outerPosition().z();

    glo_idx.push_back(i);
    glo_isSignal.push_back(isSignal);
    glo_isSimTrack.push_back(isSimTrack);
    glo_innerZ.push_back(innerZ);
    glo_outerZ.push_back(outerZ);
    glo_innerMomentum.push_back(inner);
    glo_outerMomentum.push_back(outer);
    glo_nrRechits.push_back(gmtkx[i].recHitsSize());
    glo_innerOk.push_back(gmtkx[i].innerOk());
    glo_outerOk.push_back(gmtkx[i].outerOk());
    glo_found.push_back(gmtkx[i].found());
    glo_lost.push_back(gmtkx[i].lost());
  }


  std::cout << "------------------ Standalone Muon Track -------------------" << std::endl;
  for (edm::View<reco::Track>::size_type i = 0; i < smtkx.size(); ++i) {
    RefToBase<reco::Track> track(standaloneMuonTracks_h, i);  
    int isSimTrack = 0;
    int isSignal = 0;

    if(recSimColl_StandaloneMuon.find(track) != recSimColl_StandaloneMuon.end()){
      isSimTrack = 0;
      std::cout << "Standalone Muon Track is associated with a SimTrack" << std::endl;

      std::vector<std::pair<TrackingParticleRef, double> > tps;
      tps = recSimColl_StandaloneMuon[track];
      for (auto tp:tps){
        if(tpSelector_(*(tp.first))){
          isSignal=1;
          std::cout << "Standalone Muon Track is signal track" << std::endl;
          std::cout << "Inner Position: " << smtkx[i].innerPosition().z() << ", Outer Position: " << smtkx[i].outerPosition().z() << std::endl; 
          std::cout << "Signal track has index: " << i << std::endl;
          std::cout << eventidx << "," << eta_ << "," << energy_ <<"," << i << std::endl;
          std::cout << "Signal has pdgId: " << tp.first->pdgId() << std::endl;
          std::cout << "GenParticleIndex: " << tp.first->g4Tracks()[0].genpartIndex() << std::endl;
        }
      }
    }

    double inner = std::sqrt(smtkx[i].innerMomentum().Mag2());
    double outer = std::sqrt(smtkx[i].outerMomentum().Mag2());
    double innerZ = smtkx[i].innerPosition().z();
    double outerZ = smtkx[i].outerPosition().z();

    sta_idx.push_back(i);
    sta_isSignal.push_back(isSignal);
    sta_isSimTrack.push_back(isSimTrack);
    sta_innerZ.push_back(innerZ);
    sta_outerZ.push_back(outerZ);
    sta_innerMomentum.push_back(inner);
    sta_outerMomentum.push_back(outer);
    sta_nrRechits.push_back(smtkx[i].recHitsSize());
    sta_innerOk.push_back(smtkx[i].innerOk());
    sta_outerOk.push_back(smtkx[i].outerOk());
    sta_found.push_back(smtkx[i].found());
    sta_lost.push_back(smtkx[i].lost());
  }

  std::cout << "------------------ Reco Muons -------------------" << std::endl;

  const reco::Track& recoMuonTrack = *(rmuons[0].track());
  const reco::Track& recoMuonOuterTrack = *(rmuons[0].outerTrack());

  std::cout << "Inner Position: " << recoMuonTrack.innerPosition().z() << ", Outer Position: " << recoMuonTrack.outerPosition().z() << std::endl; 
  std::cout << "Inner Position: " << recoMuonOuterTrack.innerPosition().z() << ", Outer Position: " << recoMuonOuterTrack.outerPosition().z() << std::endl; 

  globalMuonTree->Fill();
  standaloneMuonTree->Fill();
  recoMuonTree->Fill();


  glo_idx.clear();
  glo_isSignal.clear();
  glo_isSimTrack.clear();
  glo_innerMomentum.clear();
  glo_outerMomentum.clear();
  glo_innerZ.clear();
  glo_outerZ.clear();
  glo_nrRechits.clear();
  glo_innerOk.clear();
  glo_outerOk.clear();
  glo_found.clear();
  glo_lost.clear();

  sta_idx.clear();
  sta_isSignal.clear();
  sta_isSimTrack.clear();
  sta_innerMomentum.clear();
  sta_outerMomentum.clear();
  sta_innerZ.clear();
  sta_outerZ.clear();
  sta_nrRechits.clear();
  sta_innerOk.clear();
  sta_outerOk.clear();
  sta_found.clear();
  sta_lost.clear();

  rmu_idx.clear();
  rmu_isSignal.clear();
  rmu_isSimTrack.clear();
  rmu_innerMomentum.clear();
  rmu_outerMomentum.clear();
  rmu_innerZ.clear();
  rmu_outerZ.clear();
  rmu_nrRechits.clear();
  rmu_innerOk.clear();
  rmu_outerOk.clear();
  rmu_found.clear();
  rmu_lost.clear();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void TrueTrackAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void TrueTrackAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrueTrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrueTrackAnalyzer);