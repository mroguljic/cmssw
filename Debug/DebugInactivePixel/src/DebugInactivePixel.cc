#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/SiPixelRawData/interface/SiPixelRawDataError.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include <iostream>
#include <vector>

class DebugInactivePixel : public edm::one::EDAnalyzer<> {
public:
  explicit DebugInactivePixel(const edm::ParameterSet& cfg) :
    tracksLabel_(cfg.getParameter<edm::InputTag>("tracksLabel")),
    refitTracksLabel_(cfg.exists("refitTracksLabel") ? cfg.getParameter<edm::InputTag>("refitTracksLabel") : edm::InputTag("TrackRefitter")),
    rawErrLabel_(cfg.getParameter<edm::InputTag>("rawErrLabel")),
    tracksToken_(consumes<std::vector<reco::Track>>(tracksLabel_)),
    refitTracksToken_(consumes<std::vector<reco::Track>>(refitTracksLabel_)),
    printedRawErrCheck_(false), haveRawErr_(false)
  {

    if (!rawErrLabel_.encode().empty()) {
      rawErrToken_ = consumes<edm::DetSetVector<SiPixelRawDataError>>(rawErrLabel_);
      haveRawErr_ = true;
    }

    std::cout << "[DebugInactivePixel]: tracksLabel='" << tracksLabel_.encode()
              << "' refitTracksLabel='" << refitTracksLabel_.encode()
              << "' rawErrLabel='" << rawErrLabel_.encode() << "'\n";
    if (haveRawErr_) std::cout << "[DebugInactivePixel]: will check raw-error collection presence on first event.\n";
  }

  void analyze(const edm::Event& e, const edm::EventSetup& es) override {
    using namespace edm;

    // One-time check for raw-error collection validity
    if (!printedRawErrCheck_ && haveRawErr_) {
      Handle<edm::DetSetVector<SiPixelRawDataError>> errH;
      e.getByToken(rawErrToken_, errH);
      std::cout << "[DebugInactivePixel]: rawErr collection ('" << rawErrLabel_.encode() << "') present? " << errH.isValid() << "\n";
      printedRawErrCheck_ = true;
    }

    // --- Print HitPattern for refitted tracks ---
    Handle<std::vector<reco::Track>> refitTracksH;
    e.getByToken(refitTracksToken_, refitTracksH);
    if(refitTracksH.isValid()) {
      for(size_t i=0;i<refitTracksH->size();++i) {
        const reco::Track &t = (*refitTracksH)[i];
        std::cout << "[REFIT] track["<<i<<"] pt="<<t.pt()<<" eta="<<t.eta()
                  <<" nValidPix="<<t.hitPattern().numberOfValidPixelHits()
                  <<" nMissingInner="<<t.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)
                  <<" nMissingOuter="<<t.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)
                  <<std::endl;

        using HP = reco::HitPattern;
        const std::vector<HP::HitCategory> cats = {HP::TRACK_HITS, HP::MISSING_INNER_HITS, HP::MISSING_OUTER_HITS};
        for (auto cat : cats) {
          int n = t.hitPattern().numberOfAllHits(cat);
          const char* catName = (cat == HP::TRACK_HITS) ? "TRACK_HITS" : (cat == HP::MISSING_INNER_HITS) ? "MISSING_INNER_HITS" : "MISSING_OUTER_HITS";
          for (int ip = 0; ip < n; ++ip) {
            const uint16_t word = t.hitPattern().getHitPattern(cat, ip);
            const uint32_t ht = HP::getHitType(word);
            const char* htStr = (ht == HP::VALID) ? "VALID" : (ht == HP::MISSING) ? "MISSING" : (ht == HP::INACTIVE) ? "INACTIVE" : "BAD";
            std::string detStr = "UNKNOWN";
            if (HP::pixelBarrelHitFilter(word)) detStr = "PixelBarrel";
            else if (HP::pixelEndcapHitFilter(word)) detStr = "PixelEndcap";
            else continue; // skip non-pixel hits
            const int layer = HP::getLayer(word);
            std::cout << "[REFIT] category=" << catName << " pos=" << ip << " type=" << htStr
                      << " det=" << detStr << " layer=" << layer << std::endl;
          }
        }
      }
    }

    // --- Print HitPattern for stored tracks ---
    Handle<std::vector<reco::Track>> tracksH;
    e.getByToken(tracksToken_, tracksH);
    if(tracksH.isValid()) {
      for(size_t i=0;i<tracksH->size();++i) {
        const reco::Track &t = (*tracksH)[i];
        std::cout << "[STORED] track["<<i<<"] pt="<<t.pt()<<" eta="<<t.eta()
                  <<" nValidPix="<<t.hitPattern().numberOfValidPixelHits()
                  <<" nMissingInner="<<t.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)
                  <<" nMissingOuter="<<t.hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)
                  <<std::endl;

        using HP = reco::HitPattern;
        const std::vector<HP::HitCategory> cats = {HP::TRACK_HITS, HP::MISSING_INNER_HITS, HP::MISSING_OUTER_HITS};
        for (auto cat : cats) {
          int n = t.hitPattern().numberOfAllHits(cat);
          const char* catName = (cat == HP::TRACK_HITS) ? "TRACK_HITS" : (cat == HP::MISSING_INNER_HITS) ? "MISSING_INNER_HITS" : "MISSING_OUTER_HITS";
          for (int ip = 0; ip < n; ++ip) {
            const uint16_t word = t.hitPattern().getHitPattern(cat, ip);
            const uint32_t ht = HP::getHitType(word);
            const char* htStr = (ht == HP::VALID) ? "VALID" : (ht == HP::MISSING) ? "MISSING" : (ht == HP::INACTIVE) ? "INACTIVE" : "BAD";
            std::string detStr = "UNKNOWN";
            if (HP::pixelBarrelHitFilter(word)) detStr = "PixelBarrel";
            else if (HP::pixelEndcapHitFilter(word)) detStr = "PixelEndcap";
            else continue; // skip non-pixel hits
            const int layer = HP::getLayer(word);
            std::cout << "[STORED] category=" << catName << " pos=" << ip << " type=" << htStr
                      << " det=" << detStr << " layer=" << layer << std::endl;
          }
        }
      }
    }
    std::cout << "-------------------------" << std::endl;
  }

private:
  edm::InputTag tracksLabel_;
  edm::InputTag refitTracksLabel_;
  edm::InputTag rawErrLabel_;

  edm::EDGetTokenT<std::vector<reco::Track>> tracksToken_;
  edm::EDGetTokenT<std::vector<reco::Track>> refitTracksToken_;
  edm::EDGetTokenT<edm::DetSetVector<SiPixelRawDataError>> rawErrToken_;

  bool printedRawErrCheck_;
  bool haveRawErr_;
};

DEFINE_FWK_MODULE(DebugInactivePixel);
