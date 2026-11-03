// Author: Mark Matthewman - mark.matthewman@cern.ch
//
// Truth study of the candidate seed collections for the HGCAL muon Kalman filter:
// generalTracks, reco muons (muons1stStep), standalone muons and GEM/GE0 segments.
// Truth comes from the truth::Graph and its LogicalGraphHitIndex
// (PhysicsTools/TruthInfo) instead of TrackingParticles, but the match definitions
// follow the TrackingParticle associators used by the validation:
//
//   tracks (quickTrackAssociatorByHits, MTV defaults)
//     matched with the truth-graph BranchHitAssociator (SharedHits metric on the
//     tracker channel), keeping the tightest of the best-scoring branches so a muon
//     wins over GEN ancestors with the same hits (as BranchTrackingValidator does).
//     purity = shared hits / valid hits; matched if purity > 0.75 and, for 3-hit
//     tracks, all 3 hits shared. With SimToRecoDenominator = "reco" the sim->reco
//     direction uses the same condition.
//   standalone muons (MuonAssociatorByHits, Validation/RecoMuon settings)
//     matched with the same BranchHitAssociator and tie-break on the muon channel;
//     segments are expanded into their component hits (DT 1D hits on their wire,
//     CSC 2D hits on their layer), the DetIds the muon sim hits carry.
//     purity = shared / valid muon hits, efficiency = shared / sim muon hits.
//     reco->sim (fake rate): purity > 0.75.
//     sim->reco (efficiency): purity > 0.75 and efficiency > 0.5.
//   reco muons
//     global: inner and outer stub both matched to the same particle
//             (tpToGlbMuonAssociation, acceptOneStubMatchings = False);
//     tracker: inner track (tpToTkmuTrackAssociation);
//     STA-only: outer track.
//   GEM segments: no TrackingParticle reference exists; purity on the segment's
//     rechits with its own cut. A segment is local (one chamber), so it is matched on
//     a particle's own hits by HitMatcher rather than through the subtree: with the
//     subtree, a delta ray crossing the same partitions as its muon ties with it.
//
// Differences that remain with respect to the TrackingParticle associators:
//   * a reco hit is shared with a particle if the particle has a sim hit on the
//     same DetUnit (DT: same wire). The TP associators use digi-sim links, i.e.
//     the cluster / wire itself. Identical for isolated muons, looser in dense
//     environments where two particles cross the same module.
//   * sim hit counts (efficiency denominators) are distinct DetUnits of the subtree
//     rather than TrackingParticle::numberOfHits().
//   * tracks and standalone muons are matched to a particle's whole subtree (its own
//     hits plus those of its descendants: delta rays, decay products), where a
//     TrackingParticle has only its own SimTrack's hits. GEM segments use own hits
//     by default; includeDescendantHits = True switches them to the subtree.
//
// One row per object in each tree. Efficiency, fake rate and duplicate rate are
// formed from the flags at plotting time (test/plotMuonSeedTruth.py):
//   efficiency : simMuons with <class>_nMatched >= 1 / simMuons
//   duplicates : simMuons with <class>_nMatched >= 2 / simMuons
//   fake rate  : reco rows with isMatched == 0 / reco rows
//
// The truth is only as complete as the hit index: with a signal-only index (built
// from g4SimHits in RECO) every pileup object is a "fake". Use an index built at
// DIGI with enableTruth (signal + pileup) for fake rates at PU200.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <span>
#include <unordered_map>
#include <vector>

#include "TTree.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/DTRecHit/interface/DTRecHit1D.h"
#include "DataFormats/DTRecHit/interface/DTRecHit1DPair.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "PhysicsTools/TruthInfo/interface/Branch.h"
#include "PhysicsTools/TruthInfo/interface/BranchHitAssociator.h"
#include "PhysicsTools/TruthInfo/interface/RecoHitAdapters.h"
#include "PhysicsTools/TruthInfo/interface/SubgraphHitView.h"
#include "SimDataFormats/TruthInfo/interface/Graph.h"
#include "SimDataFormats/TruthInfo/interface/LogicalGraphHitIndex.h"

namespace {

  constexpr float kNoValue = -999.f;

  void sortUnique(std::vector<uint32_t>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  }

  bool isMuonDetId(uint32_t rawId) { return DetId(rawId).det() == DetId::Muon; }

  bool isGE0(uint32_t rawId) {
    return isMuonDetId(rawId) && DetId(rawId).subdetId() == MuonSubdetId::GEM && GEMDetId(rawId).station() == 0;
  }

  // Granularity at which reco and sim muon hits are compared: the DetUnit, except
  // for DT where sim hits sit on a wire and reco hits on a layer (and left/right
  // ambiguity makes the wire unreliable), so both go to the layer. 0 = not a muon hit.
  uint32_t muonHitKey(uint32_t rawId) {
    if (!isMuonDetId(rawId))
      return 0;
    if (DetId(rawId).subdetId() == MuonSubdetId::DT)
      return DTLayerId(rawId).rawId();
    return rawId;
  }

  uint32_t hgcHitKey(uint32_t rawId) { return rawId; }

  // Chamber DetId of any muon DetId, for the sim-muon summary only.
  uint32_t muonChamberKey(uint32_t rawId) {
    if (!isMuonDetId(rawId))
      return 0;
    switch (DetId(rawId).subdetId()) {
      case MuonSubdetId::DT:
        return DTChamberId(rawId).rawId();
      case MuonSubdetId::CSC:
        return CSCDetId(rawId).chamberId().rawId();
      case MuonSubdetId::RPC:
        return RPCDetId(rawId).chamberId().rawId();
      case MuonSubdetId::GEM:
        return GEMDetId(rawId).chamberId().rawId();
      case MuonSubdetId::ME0:
        return ME0DetId(rawId).chamberId().rawId();
      default:
        return 0;
    }
  }

  // A muon hit and, for segments, its component hits, as MuonAssociatorByHits counts
  // them: a DT 4D segment -> 2D segments -> 1D hits, a CSC/GEM segment -> its 2D hits.
  // Each is given the DetId its sim hits carry in the hit index: the wire for DT (a
  // 1D pair is one hit, not its two left/right halves), the layer for CSC, the roll /
  // eta partition for RPC and GEM.
  void appendMuonRecoHits(TrackingRecHit const& hit, std::vector<truth::RecoHit>& out) {
    if (!hit.isValid())
      return;
    if (auto const* pair = dynamic_cast<DTRecHit1DPair const*>(&hit)) {
      out.push_back(truth::RecoHit{pair->wireId().rawId(), 1.f, 1.f});
      return;
    }
    if (auto const* dt = dynamic_cast<DTRecHit1D const*>(&hit)) {
      out.push_back(truth::RecoHit{dt->wireId().rawId(), 1.f, 1.f});
      return;
    }
    const auto components = hit.recHits();
    if (components.empty()) {
      if (isMuonDetId(hit.geographicalId().rawId()))
        out.push_back(truth::RecoHit{hit.geographicalId().rawId(), 1.f, 1.f});
      return;
    }
    for (auto const* component : components)
      appendMuonRecoHits(*component, out);
  }

  std::vector<truth::RecoHit> muonRecoHits(reco::Track const& track) {
    std::vector<truth::RecoHit> hits;
    for (auto it = track.recHitsBegin(); it != track.recHitsEnd(); ++it)
      appendMuonRecoHits(**it, hits);
    return hits;
  }

  // Result of matching one reco object to the truth graph.
  struct HitMatch {
    int particle = -1;      // best particle, -1 if none shares a hit
    uint32_t shared = 0;    // reco hits shared with it
    uint32_t nReco = 0;     // reco hits considered
    uint32_t nSim = 0;      // its distinct sim hit keys
    float purity = 0.f;     // shared / nReco
    float efficiency = 0.f; // shared / nSim
    int nSharing = 0;       // particles sharing >= mergeThreshold of the reco hits
  };

  // Truth-graph BranchHitAssociator (SharedHits metric) for one reco object: candidates
  // are every particle's subtree. A particle and its ancestors share the same hits and
  // score, so among the best-scoring matches keep the smallest subtree in the channel.
  // With the SharedHits metric, BranchMatch::sharedEnergy holds the number of shared
  // hits (not an energy), score = 1 - shared / reco hits and reverseScore = 1 - shared
  // / branch hits.
  HitMatch bestBranchMatch(truth::BranchHitAssociator const& assoc,
                           truth::SubgraphHitView& view,
                           truth::HitChannel channel,
                           std::vector<truth::RecoHit> const& hits,
                           double mergeThreshold) {
    HitMatch m;
    m.nReco = static_cast<uint32_t>(hits.size());
    if (hits.empty())
      return m;
    const auto matches = assoc.bestBranches(std::span<const truth::RecoHit>(hits));
    if (matches.empty())
      return m;

    const float bestScore = matches.front().score;
    auto best = matches.front();
    std::size_t bestSize = view.subgraphHits(channel, best.rootParticleId).size();
    for (auto const& match : matches) {
      if (match.sharedEnergy >= mergeThreshold * m.nReco)
        ++m.nSharing;
      if (match.score > bestScore)
        continue;
      const std::size_t size = view.subgraphHits(channel, match.rootParticleId).size();
      if (size < bestSize) {
        bestSize = size;
        best = match;
      }
    }
    m.particle = static_cast<int>(best.rootParticleId);
    m.shared = static_cast<uint32_t>(best.sharedEnergy);
    m.nSim = static_cast<uint32_t>(bestSize);
    m.purity = 1.f - best.score;
    m.efficiency = 1.f - best.reverseScore;
    return m;
  }

  // Matches reco hits to particles by shared DetUnits in one hit channel. Candidates
  // are the particles with their own hits in the channel; a reco hit counts as shared
  // with a candidate that has a sim hit with the same key. Each reco hit is counted
  // once, as in the TrackingParticle associators. The best candidate shares the most
  // hits. Ties: with own hits only, the candidate with the most sim keys wins, so a
  // muon beats the delta rays that cross the same DetUnits with it; with descendants
  // included, the fewest keys wins, so a muon beats the ancestors whose subtree
  // contains it. Then the lowest particle id.
  class HitMatcher {
  public:
    using KeyFn = uint32_t (*)(uint32_t);

    HitMatcher(truth::SubgraphHitView& view, truth::HitChannel channel, KeyFn keyFn, bool includeDescendants)
        : preferFewerKeys_(includeDescendants) {
      if (!view.hasChannel(channel))
        return;
      std::vector<uint32_t> keys;
      for (uint32_t p = 0; p < view.nParticles(); ++p) {
        const auto direct = view.directHits(channel, p);
        if (direct.empty())
          continue;
        keys.clear();
        for (auto const& h : includeDescendants ? view.subgraphHits(channel, p) : direct)
          if (const uint32_t k = keyFn(h.detId); k != 0)
            keys.push_back(k);
        sortUnique(keys);
        if (keys.empty())
          continue;
        const auto slot = static_cast<uint32_t>(candidates_.size());
        for (const uint32_t k : keys)
          keyToSlots_[k].push_back(slot);
        candidates_.push_back(Candidate{p, static_cast<uint32_t>(keys.size())});
        nSimOf_[p] = static_cast<uint32_t>(keys.size());
      }
    }

    // Distinct sim hit keys of a particle in this channel (0 if it has none).
    [[nodiscard]] uint32_t nSimKeys(uint32_t particle) const {
      auto it = nSimOf_.find(particle);
      return it != nSimOf_.end() ? it->second : 0;
    }

    [[nodiscard]] HitMatch match(std::vector<uint32_t> const& recoKeys, double mergeThreshold) const {
      HitMatch out;
      out.nReco = static_cast<uint32_t>(recoKeys.size());
      if (recoKeys.empty())
        return out;

      std::unordered_map<uint32_t, uint32_t> sharedBySlot;
      for (const uint32_t k : recoKeys) {
        auto it = keyToSlots_.find(k);
        if (it == keyToSlots_.end())
          continue;
        for (const uint32_t slot : it->second)
          ++sharedBySlot[slot];
      }

      int bestSlot = -1;
      for (auto const& [slot, shared] : sharedBySlot) {
        if (shared >= mergeThreshold * out.nReco)
          ++out.nSharing;
        if (bestSlot >= 0) {
          auto const& best = candidates_[bestSlot];
          auto const& cand = candidates_[slot];
          const bool tighter = preferFewerKeys_ ? cand.nSim < best.nSim : cand.nSim > best.nSim;
          const bool better =
              shared > out.shared ||
              (shared == out.shared && (tighter || (cand.nSim == best.nSim && cand.particle < best.particle)));
          if (!better)
            continue;
        }
        bestSlot = static_cast<int>(slot);
        out.shared = shared;
      }
      if (bestSlot >= 0) {
        auto const& best = candidates_[bestSlot];
        out.particle = static_cast<int>(best.particle);
        out.nSim = best.nSim;
        out.purity = static_cast<float>(out.shared) / static_cast<float>(out.nReco);
        out.efficiency = static_cast<float>(out.shared) / static_cast<float>(best.nSim);
      }
      return out;
    }

  private:
    struct Candidate {
      uint32_t particle;
      uint32_t nSim;
    };

    bool preferFewerKeys_;
    std::vector<Candidate> candidates_;
    std::unordered_map<uint32_t, std::vector<uint32_t>> keyToSlots_;
    std::unordered_map<uint32_t, uint32_t> nSimOf_;
  };

  // --- Tree rows ---------------------------------------------------------------

  struct EventKey {
    uint32_t run = 0;
    uint32_t lumi = 0;
    uint64_t event = 0;
  };

  struct SimMuonRow {
    int particleId = -1;
    int pdgId = 0;
    int charge = 0;
    float pt = 0.f, eta = 0.f, phi = 0.f;
    float vtxRho = kNoValue, vtxZ = kNoValue;
    bool isSignal = false;
    bool hasGen = false;
    int nTrackerSimHits = 0;
    int nMuonSimHits = 0;
    int nMuonSimChambers = 0;
    int nGE0SimHits = 0;
    // generalTracks
    int gt_nMatched = 0;
    int gt_bestIdx = -1;
    float gt_bestPurity = kNoValue, gt_bestPt = kNoValue, gt_bestEta = kNoValue;
    // reco muons (muons1stStep)
    int mu_nMatched = 0;
    int mu_nGlobal = 0;
    int mu_nTracker = 0;
    int mu_nStaOnly = 0;
    // standalone muons
    int sta_nMatched = 0;
    int sta_bestIdx = -1;
    float sta_bestPurity = kNoValue, sta_bestEfficiency = kNoValue, sta_bestPt = kNoValue, sta_bestEta = kNoValue;
    // GEM segments (all stations) and GE0 only
    int seg_nMatched = 0;
    int ge0_nMatched = 0;
    int ge0_bestIdx = -1;
    float ge0_bestPurity = kNoValue;
    // HGCAL Kalman-filter tracks
    int hgcTrack_nMatched = 0;
    int hgcTrack_bestIdx = -1;
    float hgcTrack_bestPurity = kNoValue, hgcTrack_bestEfficiency = kNoValue;
  };

  // Truth fields common to every reco tree.
  struct RecoTruth {
    float purity = 0.f;
    float efficiency = 0.f;
    int nShared = 0;
    int nReco = 0;
    int nSharing = 0;
    int bestParticle = -1;
    int bestPdgId = 0;
    bool bestIsSignal = false;
    bool isMatched = false;     // reco->sim match (fake rate: !isMatched)
    bool isSimMatched = false;  // sim->reco match (counts towards the sim muon's efficiency)
    int matchedSimMuon = -1;    // particle id if isSimMatched to a selected sim muon
    bool isDuplicate = false;   // its sim muon has >= 2 sim-matched objects of this class
  };

  struct TrackRow {
    int idx = -1;
    float pt = 0.f, eta = 0.f, phi = 0.f;
    float outerEta = kNoValue, outerPhi = kNoValue, outerZ = kNoValue;  // outermost hit
    int charge = 0;
    int nValidHits = 0;
    bool highPurity = false;
    bool usedByGlobal = false;
    bool usedByTracker = false;
    RecoTruth truth;
  };

  struct RecoMuonRow {
    int idx = -1;
    float pt = 0.f, eta = 0.f, phi = 0.f;
    bool isGlobal = false, isTracker = false, isStandAlone = false, isStaOnly = false;
    int innerIdx = -1;    // index in generalTracks, -1 if none or another collection
    int matchSource = 0;  // 0 none, 1 inner track, 2 outer track, 3 both stubs (global)
    RecoTruth truth;
  };

  struct StaRow {
    int idx = -1;
    float pt = 0.f, eta = 0.f, phi = 0.f;
    int charge = 0;
    int nValidHits = 0;
    RecoTruth truth;
  };

  struct SegmentRow {
    int idx = -1;
    int region = 0, station = -1, chamber = 0;
    float eta = kNoValue, phi = kNoValue;
    float dirEta = kNoValue, dirPhi = kNoValue;
    int nRecHits = 0;
    float chi2 = 0.f;
    int ndof = 0;
    RecoTruth truth;
  };

  // A track of the HGCAL Kalman filter. Its pt/eta/phi are those of its seed (the KF
  // does not update them yet).
  struct HGCTrackRow {
    int idx = -1;
    float pt = 0.f, eta = 0.f, phi = 0.f;
    int seedTrackIdx = -1;      // seed's index in generalTracks, -1 if the seed was not a tracker track
    int seedBestParticle = -1;  // truth of the seed track (-1 if none or not a tracker track)
    bool seedIsMatched = false;
    bool sameParticleAsSeed = false;  // HGCAL track and seed matched to the same particle
    int nValidHits = 0;
    RecoTruth truth;
  };

  void branchKey(TTree* t, EventKey& k) {
    t->Branch("run", &k.run);
    t->Branch("lumi", &k.lumi);
    t->Branch("event", &k.event);
  }

  void branchTruth(TTree* t, RecoTruth& r) {
    t->Branch("purity", &r.purity);
    t->Branch("efficiency", &r.efficiency);
    t->Branch("nShared", &r.nShared);
    t->Branch("nReco", &r.nReco);
    t->Branch("nSharing", &r.nSharing);
    t->Branch("bestParticle", &r.bestParticle);
    t->Branch("bestPdgId", &r.bestPdgId);
    t->Branch("bestIsSignal", &r.bestIsSignal);
    t->Branch("isMatched", &r.isMatched);
    t->Branch("isSimMatched", &r.isSimMatched);
    t->Branch("matchedSimMuon", &r.matchedSimMuon);
    t->Branch("isDuplicate", &r.isDuplicate);
  }

}  // namespace

class MuonSeedTruthAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MuonSeedTruthAnalyzer(edm::ParameterSet const&);
  void analyze(edm::Event const&, edm::EventSetup const&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  [[nodiscard]] bool selectSimMuon(truth::Graph const& graph, uint32_t p, float& vtxRho, float& vtxZ) const;
  [[nodiscard]] RecoTruth describe(truth::Graph const& graph, HitMatch const& m) const;
  [[nodiscard]] RecoTruth trackTruth(truth::Graph const& graph,
                                     truth::BranchHitAssociator const& assoc,
                                     truth::SubgraphHitView& view,
                                     reco::Track const&) const;
  [[nodiscard]] RecoTruth staTruth(truth::Graph const& graph,
                                   truth::BranchHitAssociator const& assoc,
                                   truth::SubgraphHitView& view,
                                   reco::Track const&) const;
  void setSimMuon(RecoTruth& t, std::unordered_map<uint32_t, std::size_t> const& simRowOf) const;

  const edm::EDGetTokenT<truth::Graph> graphToken_;
  const edm::EDGetTokenT<truth::LogicalGraphHitIndex> hitIndexToken_;
  const edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  const edm::EDGetTokenT<reco::MuonCollection> muonsToken_;
  const edm::EDGetTokenT<reco::TrackCollection> staToken_;
  const edm::EDGetTokenT<GEMSegmentCollection> gemSegmentsToken_;
  const edm::EDGetTokenT<reco::TrackCollection> hgcTracksToken_;
  const edm::EDGetTokenT<std::vector<int>> hgcSeedTrackIdsToken_;
  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;

  const double simMuonPtMin_;
  const double simMuonMinAbsEta_;
  const double simMuonMaxAbsEta_;
  const double simMuonTip_;
  const double simMuonLip_;
  const bool simMuonSignalOnly_;
  const double trackMinAbsEta_;
  const double trackMaxAbsEta_;
  const double trackPtMin_;
  const double trackPurityCut_;
  const bool threeHitTracksAreSpecial_;
  const double staPurityCut_;
  const double staEfficiencyCut_;
  const double segmentPurityCut_;
  const double hgcTrackPurityCut_;
  const double hgcTrackEfficiencyCut_;
  const double mergeThreshold_;
  const bool includeDescendantHits_;

  TTree* eventTree_ = nullptr;
  TTree* simMuonTree_ = nullptr;
  TTree* trackTree_ = nullptr;
  TTree* recoMuonTree_ = nullptr;
  TTree* staTree_ = nullptr;
  TTree* segmentTree_ = nullptr;
  TTree* hgcTrackTree_ = nullptr;

  EventKey key_;
  int nSimMuons_ = 0, nTracks_ = 0, nTracksInWindow_ = 0, nRecoMuons_ = 0, nSta_ = 0, nGEMSegments_ = 0,
      nGE0Segments_ = 0, nHGCTracks_ = 0;
  SimMuonRow simRow_;
  TrackRow trackRow_;
  RecoMuonRow muonRow_;
  StaRow staRow_;
  SegmentRow segRow_;
  HGCTrackRow hgcTrackRow_;

  bool warnedNoMuonChannel_ = false;
  bool warnedSeedTrackIds_ = false;
};

MuonSeedTruthAnalyzer::MuonSeedTruthAnalyzer(edm::ParameterSet const& cfg)
    : graphToken_(consumes<truth::Graph>(cfg.getParameter<edm::InputTag>("graph"))),
      hitIndexToken_(consumes<truth::LogicalGraphHitIndex>(cfg.getParameter<edm::InputTag>("hitIndex"))),
      tracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("tracks"))),
      muonsToken_(consumes<reco::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"))),
      staToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("standAloneMuons"))),
      gemSegmentsToken_(consumes<GEMSegmentCollection>(cfg.getParameter<edm::InputTag>("gemSegments"))),
      hgcTracksToken_(consumes<reco::TrackCollection>(cfg.getParameter<edm::InputTag>("hgcTracks"))),
      hgcSeedTrackIdsToken_(consumes<std::vector<int>>(cfg.getParameter<edm::InputTag>("hgcSeedTrackIds"))),
      gemGeomToken_(esConsumes<GEMGeometry, MuonGeometryRecord>()),
      simMuonPtMin_(cfg.getParameter<double>("simMuonPtMin")),
      simMuonMinAbsEta_(cfg.getParameter<double>("simMuonMinAbsEta")),
      simMuonMaxAbsEta_(cfg.getParameter<double>("simMuonMaxAbsEta")),
      simMuonTip_(cfg.getParameter<double>("simMuonTip")),
      simMuonLip_(cfg.getParameter<double>("simMuonLip")),
      simMuonSignalOnly_(cfg.getParameter<bool>("simMuonSignalOnly")),
      trackMinAbsEta_(cfg.getParameter<double>("trackMinAbsEta")),
      trackMaxAbsEta_(cfg.getParameter<double>("trackMaxAbsEta")),
      trackPtMin_(cfg.getParameter<double>("trackPtMin")),
      trackPurityCut_(cfg.getParameter<double>("trackPurityCut")),
      threeHitTracksAreSpecial_(cfg.getParameter<bool>("threeHitTracksAreSpecial")),
      staPurityCut_(cfg.getParameter<double>("staPurityCut")),
      staEfficiencyCut_(cfg.getParameter<double>("staEfficiencyCut")),
      segmentPurityCut_(cfg.getParameter<double>("segmentPurityCut")),
      hgcTrackPurityCut_(cfg.getParameter<double>("hgcTrackPurityCut")),
      hgcTrackEfficiencyCut_(cfg.getParameter<double>("hgcTrackEfficiencyCut")),
      mergeThreshold_(cfg.getParameter<double>("mergeThreshold")),
      includeDescendantHits_(cfg.getParameter<bool>("includeDescendantHits")) {
  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;

  eventTree_ = fs->make<TTree>("events", "Per-event object counts");
  branchKey(eventTree_, key_);
  eventTree_->Branch("nSimMuons", &nSimMuons_);
  eventTree_->Branch("nTracks", &nTracks_);
  eventTree_->Branch("nTracksInWindow", &nTracksInWindow_);
  eventTree_->Branch("nRecoMuons", &nRecoMuons_);
  eventTree_->Branch("nSta", &nSta_);
  eventTree_->Branch("nGEMSegments", &nGEMSegments_);
  eventTree_->Branch("nGE0Segments", &nGE0Segments_);
  eventTree_->Branch("nHGCTracks", &nHGCTracks_);

  simMuonTree_ = fs->make<TTree>("simMuons", "Selected sim muons (efficiency / duplicate denominator)");
  branchKey(simMuonTree_, key_);
  simMuonTree_->Branch("particleId", &simRow_.particleId);
  simMuonTree_->Branch("pdgId", &simRow_.pdgId);
  simMuonTree_->Branch("charge", &simRow_.charge);
  simMuonTree_->Branch("pt", &simRow_.pt);
  simMuonTree_->Branch("eta", &simRow_.eta);
  simMuonTree_->Branch("phi", &simRow_.phi);
  simMuonTree_->Branch("vtxRho", &simRow_.vtxRho);
  simMuonTree_->Branch("vtxZ", &simRow_.vtxZ);
  simMuonTree_->Branch("isSignal", &simRow_.isSignal);
  simMuonTree_->Branch("hasGen", &simRow_.hasGen);
  simMuonTree_->Branch("nTrackerSimHits", &simRow_.nTrackerSimHits);
  simMuonTree_->Branch("nMuonSimHits", &simRow_.nMuonSimHits);
  simMuonTree_->Branch("nMuonSimChambers", &simRow_.nMuonSimChambers);
  simMuonTree_->Branch("nGE0SimHits", &simRow_.nGE0SimHits);
  simMuonTree_->Branch("gt_nMatched", &simRow_.gt_nMatched);
  simMuonTree_->Branch("gt_bestIdx", &simRow_.gt_bestIdx);
  simMuonTree_->Branch("gt_bestPurity", &simRow_.gt_bestPurity);
  simMuonTree_->Branch("gt_bestPt", &simRow_.gt_bestPt);
  simMuonTree_->Branch("gt_bestEta", &simRow_.gt_bestEta);
  simMuonTree_->Branch("mu_nMatched", &simRow_.mu_nMatched);
  simMuonTree_->Branch("mu_nGlobal", &simRow_.mu_nGlobal);
  simMuonTree_->Branch("mu_nTracker", &simRow_.mu_nTracker);
  simMuonTree_->Branch("mu_nStaOnly", &simRow_.mu_nStaOnly);
  simMuonTree_->Branch("sta_nMatched", &simRow_.sta_nMatched);
  simMuonTree_->Branch("sta_bestIdx", &simRow_.sta_bestIdx);
  simMuonTree_->Branch("sta_bestPurity", &simRow_.sta_bestPurity);
  simMuonTree_->Branch("sta_bestEfficiency", &simRow_.sta_bestEfficiency);
  simMuonTree_->Branch("sta_bestPt", &simRow_.sta_bestPt);
  simMuonTree_->Branch("sta_bestEta", &simRow_.sta_bestEta);
  simMuonTree_->Branch("seg_nMatched", &simRow_.seg_nMatched);
  simMuonTree_->Branch("ge0_nMatched", &simRow_.ge0_nMatched);
  simMuonTree_->Branch("ge0_bestIdx", &simRow_.ge0_bestIdx);
  simMuonTree_->Branch("ge0_bestPurity", &simRow_.ge0_bestPurity);
  simMuonTree_->Branch("hgcTrack_nMatched", &simRow_.hgcTrack_nMatched);
  simMuonTree_->Branch("hgcTrack_bestIdx", &simRow_.hgcTrack_bestIdx);
  simMuonTree_->Branch("hgcTrack_bestPurity", &simRow_.hgcTrack_bestPurity);
  simMuonTree_->Branch("hgcTrack_bestEfficiency", &simRow_.hgcTrack_bestEfficiency);

  trackTree_ = fs->make<TTree>("tracks", "generalTracks inside the eta window");
  branchKey(trackTree_, key_);
  trackTree_->Branch("idx", &trackRow_.idx);
  trackTree_->Branch("pt", &trackRow_.pt);
  trackTree_->Branch("eta", &trackRow_.eta);
  trackTree_->Branch("phi", &trackRow_.phi);
  trackTree_->Branch("outerEta", &trackRow_.outerEta);
  trackTree_->Branch("outerPhi", &trackRow_.outerPhi);
  trackTree_->Branch("outerZ", &trackRow_.outerZ);
  trackTree_->Branch("charge", &trackRow_.charge);
  trackTree_->Branch("nValidHits", &trackRow_.nValidHits);
  trackTree_->Branch("highPurity", &trackRow_.highPurity);
  trackTree_->Branch("usedByGlobal", &trackRow_.usedByGlobal);
  trackTree_->Branch("usedByTracker", &trackRow_.usedByTracker);
  branchTruth(trackTree_, trackRow_.truth);

  recoMuonTree_ = fs->make<TTree>("recoMuons", "Reco muons");
  branchKey(recoMuonTree_, key_);
  recoMuonTree_->Branch("idx", &muonRow_.idx);
  recoMuonTree_->Branch("pt", &muonRow_.pt);
  recoMuonTree_->Branch("eta", &muonRow_.eta);
  recoMuonTree_->Branch("phi", &muonRow_.phi);
  recoMuonTree_->Branch("isGlobal", &muonRow_.isGlobal);
  recoMuonTree_->Branch("isTracker", &muonRow_.isTracker);
  recoMuonTree_->Branch("isStandAlone", &muonRow_.isStandAlone);
  recoMuonTree_->Branch("isStaOnly", &muonRow_.isStaOnly);
  recoMuonTree_->Branch("innerIdx", &muonRow_.innerIdx);
  recoMuonTree_->Branch("matchSource", &muonRow_.matchSource);
  branchTruth(recoMuonTree_, muonRow_.truth);

  staTree_ = fs->make<TTree>("standAloneMuons", "Standalone muon tracks");
  branchKey(staTree_, key_);
  staTree_->Branch("idx", &staRow_.idx);
  staTree_->Branch("pt", &staRow_.pt);
  staTree_->Branch("eta", &staRow_.eta);
  staTree_->Branch("phi", &staRow_.phi);
  staTree_->Branch("charge", &staRow_.charge);
  staTree_->Branch("nValidHits", &staRow_.nValidHits);
  branchTruth(staTree_, staRow_.truth);

  segmentTree_ = fs->make<TTree>("gemSegments", "GEM segments (station 0 = GE0)");
  branchKey(segmentTree_, key_);
  segmentTree_->Branch("idx", &segRow_.idx);
  segmentTree_->Branch("region", &segRow_.region);
  segmentTree_->Branch("station", &segRow_.station);
  segmentTree_->Branch("chamber", &segRow_.chamber);
  segmentTree_->Branch("eta", &segRow_.eta);
  segmentTree_->Branch("phi", &segRow_.phi);
  segmentTree_->Branch("dirEta", &segRow_.dirEta);
  segmentTree_->Branch("dirPhi", &segRow_.dirPhi);
  segmentTree_->Branch("nRecHits", &segRow_.nRecHits);
  segmentTree_->Branch("chi2", &segRow_.chi2);
  segmentTree_->Branch("ndof", &segRow_.ndof);
  branchTruth(segmentTree_, segRow_.truth);

  hgcTrackTree_ = fs->make<TTree>("HGCTrackTree", "HGCAL Kalman-filter tracks");
  branchKey(hgcTrackTree_, key_);
  hgcTrackTree_->Branch("idx", &hgcTrackRow_.idx);
  hgcTrackTree_->Branch("pt", &hgcTrackRow_.pt);
  hgcTrackTree_->Branch("eta", &hgcTrackRow_.eta);
  hgcTrackTree_->Branch("phi", &hgcTrackRow_.phi);
  hgcTrackTree_->Branch("seedTrackIdx", &hgcTrackRow_.seedTrackIdx);
  hgcTrackTree_->Branch("seedBestParticle", &hgcTrackRow_.seedBestParticle);
  hgcTrackTree_->Branch("seedIsMatched", &hgcTrackRow_.seedIsMatched);
  hgcTrackTree_->Branch("sameParticleAsSeed", &hgcTrackRow_.sameParticleAsSeed);
  hgcTrackTree_->Branch("nValidHits", &hgcTrackRow_.nValidHits);
  branchTruth(hgcTrackTree_, hgcTrackRow_.truth);
}

// TrackingParticleSelector-like selection: muon, kinematics, production vertex
// (tip = transverse, lip = longitudinal distance from the origin, in cm).
bool MuonSeedTruthAnalyzer::selectSimMuon(truth::Graph const& graph, uint32_t p, float& vtxRho, float& vtxZ) const {
  auto const& pd = graph.particles()[p];
  if (std::abs(pd.pdgId) != 13 || !pd.hasSim())
    return false;
  const double absEta = std::abs(pd.momentum.eta());
  if (pd.momentum.pt() < simMuonPtMin_ || absEta < simMuonMinAbsEta_ || absEta > simMuonMaxAbsEta_)
    return false;

  // The vertex position prefers SIM (cm) and falls back to GEN (mm), so the cut is
  // applied only to SIM vertices.
  vtxRho = kNoValue;
  vtxZ = kNoValue;
  for (auto const& v : graph.particle(p).productionVertices()) {
    if (!v.hasSim())
      continue;
    vtxRho = std::hypot(v.position().X(), v.position().Y());
    vtxZ = v.position().Z();
    return vtxRho <= simMuonTip_ && std::abs(vtxZ) <= simMuonLip_;
  }
  return true;
}

RecoTruth MuonSeedTruthAnalyzer::describe(truth::Graph const& graph, HitMatch const& m) const {
  RecoTruth r;
  r.purity = m.purity;
  r.efficiency = m.efficiency;
  r.nShared = static_cast<int>(m.shared);
  r.nReco = static_cast<int>(m.nReco);
  r.nSharing = m.nSharing;
  r.bestParticle = m.particle;
  if (m.particle >= 0) {
    r.bestPdgId = graph.particles()[m.particle].pdgId;
    r.bestIsSignal = truth::Branch(&graph, static_cast<uint32_t>(m.particle)).isSignal();
  }
  return r;
}

// BranchHitAssociator on the tracker channel, then the quickTrackAssociatorByHits
// decision: purity > cut both ways, 3-hit tracks need all 3 hits.
RecoTruth MuonSeedTruthAnalyzer::trackTruth(truth::Graph const& graph,
                                            truth::BranchHitAssociator const& assoc,
                                            truth::SubgraphHitView& view,
                                            reco::Track const& track) const {
  const HitMatch m =
      bestBranchMatch(assoc, view, truth::HitChannel::Tracker, truth::recoHits(track), mergeThreshold_);
  RecoTruth r = describe(graph, m);
  r.isMatched = m.particle >= 0 && m.purity > trackPurityCut_ &&
                !(threeHitTracksAreSpecial_ && m.nReco == 3 && m.shared < 3);
  r.isSimMatched = r.isMatched;
  return r;
}

// BranchHitAssociator on the muon channel, then the MuonAssociatorByHits decision:
// purity > cut reco->sim; purity and efficiency > cuts sim->reco.
RecoTruth MuonSeedTruthAnalyzer::staTruth(truth::Graph const& graph,
                                          truth::BranchHitAssociator const& assoc,
                                          truth::SubgraphHitView& view,
                                          reco::Track const& track) const {
  const HitMatch m = bestBranchMatch(assoc, view, truth::HitChannel::Muon, muonRecoHits(track), mergeThreshold_);
  RecoTruth r = describe(graph, m);
  r.isMatched = m.particle >= 0 && m.purity > staPurityCut_;
  r.isSimMatched = r.isMatched && m.efficiency > staEfficiencyCut_;
  return r;
}

void MuonSeedTruthAnalyzer::setSimMuon(RecoTruth& t, std::unordered_map<uint32_t, std::size_t> const& simRowOf) const {
  if (t.isSimMatched && simRowOf.count(static_cast<uint32_t>(t.bestParticle)))
    t.matchedSimMuon = t.bestParticle;
}

void MuonSeedTruthAnalyzer::analyze(edm::Event const& event, edm::EventSetup const& setup) {
  key_.run = event.id().run();
  key_.lumi = event.id().luminosityBlock();
  key_.event = event.id().event();

  auto const& graph = event.get(graphToken_);
  auto const& hitIndex = event.get(hitIndexToken_);
  truth::SubgraphHitView view(hitIndex);

  if (!view.hasChannel(truth::HitChannel::Muon) && !warnedNoMuonChannel_) {
    edm::LogWarning("MuonSeedTruthAnalyzer")
        << "The hit index has no muon channel: standalone, reco-muon (outer track) and GEM segment truth "
           "will all be empty. Build the index with the Muon subdetector enabled.";
    warnedNoMuonChannel_ = true;
  }

  const truth::BranchHitAssociator trackerAssoc(
      hitIndex, {}, truth::BranchHitAssociator::Metric::SharedHits, truth::HitChannel::Tracker);
  const truth::BranchHitAssociator muonAssoc(
      hitIndex, {}, truth::BranchHitAssociator::Metric::SharedHits, truth::HitChannel::Muon);
  // GEM segments only (see the header): a segment is local, so it is matched on own
  // hits rather than through the subtree associator.
  const HitMatcher muonMatcher(view, truth::HitChannel::Muon, muonHitKey, includeDescendantHits_);
  // HGCAL KF tracks: their hits are HGCAL cells, which the hit index keeps (under reco
  // DetIds) in the calorimeter channel. Subtree matching, as for tracks and STAs.
  const HitMatcher hgcTrackMatcher(view, truth::HitChannel::Calo, hgcHitKey, false);

  // --- Sim muons: the efficiency / duplicate denominator ---
  std::vector<SimMuonRow> simRows;
  std::unordered_map<uint32_t, std::size_t> simRowOf;
  for (uint32_t p = 0; p < graph.nParticles(); ++p) {
    SimMuonRow row;
    if (!selectSimMuon(graph, p, row.vtxRho, row.vtxZ))
      continue;
    row.isSignal = truth::Branch(&graph, p).isSignal();
    if (simMuonSignalOnly_ && !row.isSignal)
      continue;

    auto const& pd = graph.particles()[p];
    row.particleId = static_cast<int>(p);
    row.pdgId = pd.pdgId;
    row.charge = pd.pdgId > 0 ? -1 : 1;
    row.pt = pd.momentum.pt();
    row.eta = pd.momentum.eta();
    row.phi = pd.momentum.phi();
    row.hasGen = pd.hasGen();
    row.nTrackerSimHits = static_cast<int>(view.directHits(truth::HitChannel::Tracker, p).size());
    row.nMuonSimHits = static_cast<int>(muonMatcher.nSimKeys(p));
    std::vector<uint32_t> chambers;
    for (auto const& h : view.directHits(truth::HitChannel::Muon, p)) {
      if (const uint32_t c = muonChamberKey(h.detId); c != 0)
        chambers.push_back(c);
      if (isGE0(h.detId))
        ++row.nGE0SimHits;
    }
    sortUnique(chambers);
    row.nMuonSimChambers = static_cast<int>(chambers.size());

    simRowOf[p] = simRows.size();
    simRows.push_back(row);
  }

  // --- generalTracks ---
  auto const tracksHandle = event.getHandle(tracksToken_);
  auto const& tracks = *tracksHandle;
  std::vector<RecoTruth> trackTruths(tracks.size());
  std::vector<TrackRow> trackRows;
  std::vector<int> trackRowOf(tracks.size(), -1);
  nTracksInWindow_ = 0;
  for (std::size_t i = 0; i < tracks.size(); ++i) {
    auto const& track = tracks[i];
    trackTruths[i] = trackTruth(graph, trackerAssoc, view, track);
    setSimMuon(trackTruths[i], simRowOf);

    if (trackTruths[i].matchedSimMuon >= 0) {
      auto& sim = simRows[simRowOf.at(trackTruths[i].matchedSimMuon)];
      ++sim.gt_nMatched;
      if (trackTruths[i].purity > sim.gt_bestPurity) {
        sim.gt_bestIdx = static_cast<int>(i);
        sim.gt_bestPurity = trackTruths[i].purity;
        sim.gt_bestPt = track.pt();
        sim.gt_bestEta = track.eta();
      }
    }

    const double absEta = std::abs(track.eta());
    if (absEta < trackMinAbsEta_ || absEta > trackMaxAbsEta_ || track.pt() < trackPtMin_)
      continue;
    ++nTracksInWindow_;
    TrackRow row;
    row.idx = static_cast<int>(i);
    row.pt = track.pt();
    row.eta = track.eta();
    row.phi = track.phi();
    if (track.extra().isAvailable()) {
      row.outerEta = track.outerPosition().eta();
      row.outerPhi = track.outerPosition().phi();
      row.outerZ = track.outerPosition().z();
    }
    row.charge = track.charge();
    row.nValidHits = track.numberOfValidHits();
    row.highPurity = track.quality(reco::TrackBase::highPurity);
    trackRowOf[i] = static_cast<int>(trackRows.size());
    trackRows.push_back(row);
  }

  // --- Reco muons ---
  auto const& muons = event.get(muonsToken_);
  std::vector<RecoMuonRow> muonRows;
  muonRows.reserve(muons.size());
  for (std::size_t i = 0; i < muons.size(); ++i) {
    auto const& mu = muons[i];
    RecoMuonRow row;
    row.idx = static_cast<int>(i);
    row.pt = mu.pt();
    row.eta = mu.eta();
    row.phi = mu.phi();
    row.isGlobal = mu.isGlobalMuon();
    row.isTracker = mu.isTrackerMuon();
    row.isStandAlone = mu.isStandAloneMuon();
    row.isStaOnly = row.isStandAlone && !row.isGlobal && !row.isTracker;

    RecoTruth inner, outer;
    const bool hasInner = mu.innerTrack().isNonnull();
    const bool hasOuter = mu.outerTrack().isNonnull();
    if (hasInner) {
      if (mu.innerTrack().id() == tracksHandle.id()) {
        row.innerIdx = static_cast<int>(mu.innerTrack().key());
        inner = trackTruths[row.innerIdx];
        if (const int r = trackRowOf[row.innerIdx]; r >= 0) {
          trackRows[r].usedByGlobal |= row.isGlobal;
          trackRows[r].usedByTracker |= row.isTracker;
        }
      } else {
        inner = trackTruth(graph, trackerAssoc, view, *mu.innerTrack());
      }
    }
    if (hasOuter)
      outer = staTruth(graph, muonAssoc, view, *mu.outerTrack());

    if (row.isGlobal && hasInner && hasOuter) {
      // Both stubs must point at the same particle; the combined purity orders
      // candidates like MuonAssociatorByHits' global quality.
      row.matchSource = 3;
      row.truth = inner;
      const int nAll = inner.nReco + outer.nReco;
      const bool same = inner.bestParticle >= 0 && inner.bestParticle == outer.bestParticle;
      row.truth.nShared = same ? inner.nShared + outer.nShared : inner.nShared;
      row.truth.nReco = nAll;
      row.truth.purity = nAll > 0 ? static_cast<float>(row.truth.nShared) / nAll : 0.f;
      row.truth.efficiency = outer.efficiency;
      row.truth.isMatched = same && inner.isMatched && outer.isMatched;
      row.truth.isSimMatched = same && inner.isSimMatched && outer.isSimMatched;
    } else if (hasInner) {
      row.matchSource = 1;
      row.truth = inner;
    } else if (hasOuter) {
      row.matchSource = 2;
      row.truth = outer;
    }
    row.truth.matchedSimMuon = -1;
    setSimMuon(row.truth, simRowOf);

    if (row.truth.matchedSimMuon >= 0) {
      auto& sim = simRows[simRowOf.at(row.truth.matchedSimMuon)];
      ++sim.mu_nMatched;
      sim.mu_nGlobal += row.isGlobal;
      sim.mu_nTracker += row.isTracker;
      sim.mu_nStaOnly += row.isStaOnly;
    }
    muonRows.push_back(row);
  }

  // --- Standalone muons ---
  auto const& stas = event.get(staToken_);
  std::vector<StaRow> staRows;
  staRows.reserve(stas.size());
  for (std::size_t i = 0; i < stas.size(); ++i) {
    auto const& sta = stas[i];
    StaRow row;
    row.idx = static_cast<int>(i);
    row.pt = sta.pt();
    row.eta = sta.eta();
    row.phi = sta.phi();
    row.charge = sta.charge();
    row.nValidHits = sta.numberOfValidHits();
    row.truth = staTruth(graph, muonAssoc, view, sta);
    setSimMuon(row.truth, simRowOf);

    if (row.truth.matchedSimMuon >= 0) {
      auto& sim = simRows[simRowOf.at(row.truth.matchedSimMuon)];
      ++sim.sta_nMatched;
      if (row.truth.purity > sim.sta_bestPurity) {
        sim.sta_bestIdx = row.idx;
        sim.sta_bestPurity = row.truth.purity;
        sim.sta_bestEfficiency = row.truth.efficiency;
        sim.sta_bestPt = row.pt;
        sim.sta_bestEta = row.eta;
      }
    }
    staRows.push_back(row);
  }

  // --- GEM segments ---
  auto const& gemGeom = setup.getData(gemGeomToken_);
  auto const& segments = event.get(gemSegmentsToken_);
  std::vector<SegmentRow> segRows;
  nGE0Segments_ = 0;
  int segIdx = 0;
  for (auto const& seg : segments) {
    const GEMDetId id = seg.gemDetId();
    SegmentRow row;
    row.idx = segIdx++;
    row.region = id.region();
    row.station = id.station();
    row.chamber = id.chamber();
    row.nRecHits = static_cast<int>(seg.specificRecHits().size());
    row.chi2 = seg.chi2();
    row.ndof = seg.degreesOfFreedom();
    if (row.station == 0)
      ++nGE0Segments_;

    // Global position and direction through the segment's own det, falling back to
    // the eta partition of its first rechit.
    if (GeomDet const* det = gemGeom.idToDet(id); det != nullptr) {
      const GlobalPoint pos = det->toGlobal(seg.localPosition());
      const GlobalVector dir = det->toGlobal(seg.localDirection());
      row.eta = pos.eta();
      row.phi = pos.phi();
      row.dirEta = dir.eta();
      row.dirPhi = dir.phi();
    } else if (!seg.specificRecHits().empty()) {
      auto const& rh = seg.specificRecHits().front();
      if (auto const* roll = gemGeom.etaPartition(rh.gemId())) {
        const GlobalPoint pos = roll->toGlobal(rh.localPosition());
        row.eta = pos.eta();
        row.phi = pos.phi();
      }
    }

    std::vector<uint32_t> keys;
    keys.reserve(seg.specificRecHits().size());
    for (auto const& rh : seg.specificRecHits())
      if (rh.isValid())
        keys.push_back(muonHitKey(rh.gemId().rawId()));
    const HitMatch m = muonMatcher.match(keys, mergeThreshold_);
    row.truth = describe(graph, m);
    row.truth.isMatched = m.particle >= 0 && m.purity > segmentPurityCut_;
    row.truth.isSimMatched = row.truth.isMatched;
    setSimMuon(row.truth, simRowOf);

    if (row.truth.matchedSimMuon >= 0) {
      auto& sim = simRows[simRowOf.at(row.truth.matchedSimMuon)];
      ++sim.seg_nMatched;
      if (row.station == 0) {
        ++sim.ge0_nMatched;
        if (row.truth.purity > sim.ge0_bestPurity) {
          sim.ge0_bestIdx = row.idx;
          sim.ge0_bestPurity = row.truth.purity;
        }
      }
    }
    segRows.push_back(row);
  }

  // --- HGCAL Kalman-filter tracks ---
  auto const& hgcTracks = event.get(hgcTracksToken_);
  // One entry per HGCAL track, in the same order: the index of its seed in
  // generalTracks, -1 if the seed was not a tracker track.
  auto const& seedTrackIds = event.get(hgcSeedTrackIdsToken_);
  const bool seedIdsUsable = seedTrackIds.size() == hgcTracks.size();
  if (!seedIdsUsable && !warnedSeedTrackIds_) {
    edm::LogWarning("MuonSeedTruthAnalyzer")
        << "HGCSeedingTrackIds has " << seedTrackIds.size() << " entries for " << hgcTracks.size()
        << " HGCAL tracks; the seed indices are ignored. They must be filled in parallel with HGCALTracks.";
    warnedSeedTrackIds_ = true;
  }

  std::vector<HGCTrackRow> hgcTrackRows;
  hgcTrackRows.reserve(hgcTracks.size());
  for (std::size_t i = 0; i < hgcTracks.size(); ++i) {
    auto const& hgcTrack = hgcTracks[i];
    HGCTrackRow row;
    row.idx = static_cast<int>(i);
    row.pt = hgcTrack.pt();
    row.eta = hgcTrack.eta();
    row.phi = hgcTrack.phi();

    // The seed's own truth, when it was a generalTracks track.
    if (seedIdsUsable) {
      const int seed = seedTrackIds[i];
      if (seed >= 0 && static_cast<std::size_t>(seed) < tracks.size()) {
        row.seedTrackIdx = seed;
        row.seedBestParticle = trackTruths[seed].bestParticle;
        row.seedIsMatched = trackTruths[seed].isMatched;
      }
    }

    // Valid hits only: layers without a compatible hit are stored as invalid hits on
    // the disk.
    std::vector<uint32_t> keys;
    for (auto it = hgcTrack.recHitsBegin(); it != hgcTrack.recHitsEnd(); ++it) {
      TrackingRecHit const& hit = **it;
      if (hit.isValid())
        keys.push_back(hgcHitKey(hit.geographicalId().rawId()));
    }
    row.nValidHits = static_cast<int>(keys.size());

    const HitMatch m = hgcTrackMatcher.match(keys, mergeThreshold_);
    row.truth = describe(graph, m);
    row.truth.isMatched = m.particle >= 0 && m.purity > hgcTrackPurityCut_;
    row.truth.isSimMatched = row.truth.isMatched && m.efficiency > hgcTrackEfficiencyCut_;
    setSimMuon(row.truth, simRowOf);
    row.sameParticleAsSeed = row.seedBestParticle >= 0 && row.seedBestParticle == row.truth.bestParticle;

    if (row.truth.matchedSimMuon >= 0) {
      auto& sim = simRows[simRowOf.at(row.truth.matchedSimMuon)];
      ++sim.hgcTrack_nMatched;
      if (row.truth.purity > sim.hgcTrack_bestPurity) {
        sim.hgcTrack_bestIdx = row.idx;
        sim.hgcTrack_bestPurity = row.truth.purity;
        sim.hgcTrack_bestEfficiency = row.truth.efficiency;
      }
    }
    hgcTrackRows.push_back(row);
  }


  // --- Duplicate flags need the complete per-sim-muon counts ---
  auto isDuplicate = [&](RecoTruth const& t, auto countOf) {
    return t.matchedSimMuon >= 0 && countOf(simRows[simRowOf.at(t.matchedSimMuon)]) >= 2;
  };

  for (std::size_t i = 0; i < tracks.size(); ++i) {
    if (const int r = trackRowOf[i]; r >= 0) {
      trackRow_ = trackRows[r];
      trackRow_.truth = trackTruths[i];
      trackRow_.truth.isDuplicate = isDuplicate(trackTruths[i], [](SimMuonRow const& s) { return s.gt_nMatched; });
      trackTree_->Fill();
    }
  }
  for (auto const& row : muonRows) {
    muonRow_ = row;
    muonRow_.truth.isDuplicate = isDuplicate(row.truth, [](SimMuonRow const& s) { return s.mu_nMatched; });
    recoMuonTree_->Fill();
  }
  for (auto const& row : staRows) {
    staRow_ = row;
    staRow_.truth.isDuplicate = isDuplicate(row.truth, [](SimMuonRow const& s) { return s.sta_nMatched; });
    staTree_->Fill();
  }
  for (auto const& row : segRows) {
    segRow_ = row;
    segRow_.truth.isDuplicate = isDuplicate(row.truth, [&row](SimMuonRow const& s) {
      return row.station == 0 ? s.ge0_nMatched : s.seg_nMatched;
    });
    segmentTree_->Fill();
  }

  for (auto const& row : hgcTrackRows) {
    hgcTrackRow_ = row;
    hgcTrackRow_.truth.isDuplicate = isDuplicate(row.truth, [](SimMuonRow const& s) { return s.hgcTrack_nMatched; });
    hgcTrackTree_->Fill();
  }
  for (auto const& row : simRows) {
    simRow_ = row;
    simMuonTree_->Fill();
  }

  nSimMuons_ = static_cast<int>(simRows.size());
  nHGCTracks_ = static_cast<int>(hgcTracks.size());
  nTracks_ = static_cast<int>(tracks.size());
  nRecoMuons_ = static_cast<int>(muons.size());
  nSta_ = static_cast<int>(stas.size());
  nGEMSegments_ = static_cast<int>(segRows.size());
  eventTree_->Fill();
}

void MuonSeedTruthAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("graph", edm::InputTag("truthLogicalGraphProducer"));
  desc.add<edm::InputTag>("hitIndex", edm::InputTag("truthLogicalGraphHitIndexProducer"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("muons", edm::InputTag("muons1stStep"));
  desc.add<edm::InputTag>("standAloneMuons", edm::InputTag("standAloneMuons"));
  desc.add<edm::InputTag>("gemSegments", edm::InputTag("gemSegments"));
  desc.add<edm::InputTag>("hgcTracks", edm::InputTag("ticlTrackstersKalmanFilter", "HGCALTracks"));
  desc.add<edm::InputTag>("hgcSeedTrackIds", edm::InputTag("ticlTrackstersKalmanFilter", "HGCSeedingTrackIds"))
      ->setComment("Per HGCAL track: index of its seed in `tracks`, -1 if the seed was not a tracker track.");
  // Sim-muon selection, as TPmu in Validation/RecoMuon (|eta| left open for HGCAL).
  desc.add<double>("simMuonPtMin", 0.9);
  desc.add<double>("simMuonMinAbsEta", 0.);
  desc.add<double>("simMuonMaxAbsEta", 5.);
  desc.add<double>("simMuonTip", 3.5)->setComment("Max production-vertex rho [cm] (SIM vertices only).");
  desc.add<double>("simMuonLip", 30.)->setComment("Max production-vertex |z| [cm] (SIM vertices only).");
  desc.add<bool>("simMuonSignalOnly", true)->setComment("Only muons from the signal interaction (bx 0, event 0).");

  desc.add<double>("trackMinAbsEta", 1.3)->setComment("Tracks outside the window are matched but not stored.");
  desc.add<double>("trackMaxAbsEta", 4.5);
  desc.add<double>("trackPtMin", 0.);

  desc.add<double>("trackPurityCut", 0.75)->setComment("quickTrackAssociatorByHits Cut_RecoToSim / Purity_SimToReco.");
  desc.add<bool>("threeHitTracksAreSpecial", true)->setComment("3-hit tracks must share all 3 hits.");
  desc.add<double>("staPurityCut", 0.75)->setComment("MuonAssociatorByHits PurityCut_muon (Validation/RecoMuon).");
  desc.add<double>("staEfficiencyCut", 0.5)
      ->setComment("MuonAssociatorByHits EfficiencyCut_muon (Validation/RecoMuon), sim->reco only.");
  desc.add<double>("segmentPurityCut", 0.5)->setComment("GEM segments; no TrackingParticle reference exists.");
  desc.add<double>("hgcTrackPurityCut", 0.5)
      ->setComment("HGCAL KF tracks: min shared / valid hits. No TrackingParticle reference exists.");
  desc.add<double>("hgcTrackEfficiencyCut", 0.)
      ->setComment(
          "HGCAL KF tracks, sim->reco only: min shared / the particle's HGCAL cells (subtree). The KF takes at "
          "most one cell per layer, so this fraction is small by construction; 0 disables the cut.");
  desc.add<double>("mergeThreshold", 0.3)
      ->setComment("Min shared fraction for a particle to count in nSharing (>= 2 -> merged object).");
  desc.add<bool>("includeDescendantHits", false)
      ->setComment(
          "GEM segments only (tracks and standalone muons always use the BranchHitAssociator subtree). False: "
          "a particle owns only its own SimTrack's hits, as a TrackingParticle. True: also its descendants'.");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(MuonSeedTruthAnalyzer);
