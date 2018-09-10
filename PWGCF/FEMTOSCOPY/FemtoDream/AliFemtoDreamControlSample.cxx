#include "AliFemtoDreamControlSample.h"
#include "AliLog.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMath.h"

ClassImp(AliFemtoDreamControlSample)

AliFemtoDreamControlSample::AliFemtoDreamControlSample()
    : fHists(nullptr),
      fPDGParticleSpecies(),
      fMultBins(),
      fRandom(),
      fPi(TMath::Pi()),
      fSpinningDepth(0),
      fStravinsky(false),
      fDeltaEtaMax(0.f),
      fDeltaPhiMax(0.f) {
  fRandom.SetSeed(0);
}

AliFemtoDreamControlSample::AliFemtoDreamControlSample(
    const AliFemtoDreamControlSample &samp)
    : fHists(samp.fHists),
      fPDGParticleSpecies(samp.fPDGParticleSpecies),
      fMultBins(samp.fMultBins),
      fRandom(),
      fPi(TMath::Pi()),
      fSpinningDepth(samp.fSpinningDepth),
      fStravinsky(false),
      fDeltaEtaMax(samp.fDeltaEtaMax),
      fDeltaPhiMax(samp.fDeltaPhiMax) {
  fRandom.SetSeed(0);
}

AliFemtoDreamControlSample::AliFemtoDreamControlSample(
    AliFemtoDreamCollConfig *conf, bool minimalBooking)
    : fHists(new AliFemtoDreamCorrHists(conf, minimalBooking)),
      fPDGParticleSpecies(conf->GetPDGCodes()),
      fMultBins(conf->GetMultBins()),
      fRandom(),
      fPi(TMath::Pi()),
      fSpinningDepth(0),
      fStravinsky(conf->GetDoStravinsky()),
      fDeltaEtaMax(conf->GetDeltaEtaMax()),
      fDeltaPhiMax(conf->GetDeltaPhiMax()) {
  if (fStravinsky) {
    fSpinningDepth = 1;
  } else {
    fSpinningDepth = conf->GetSpinningDepth();
  }
  fRandom.SetSeed(0);
}

AliFemtoDreamControlSample& AliFemtoDreamControlSample::operator=(
    const AliFemtoDreamControlSample& samp) {
  if (this == &samp) {
    return *this;
  }
  this->fHists = samp.fHists;
  this->fPDGParticleSpecies = samp.fPDGParticleSpecies;
  this->fMultBins = samp.fMultBins;
  this->fRandom.SetSeed(0);
  this->fPi = TMath::Pi();
  this->fSpinningDepth = samp.fSpinningDepth;
  this->fDeltaEtaMax = samp.fDeltaEtaMax;
  this->fDeltaPhiMax = samp.fDeltaPhiMax;
  return *this;
}

AliFemtoDreamControlSample::~AliFemtoDreamControlSample() {
}

void AliFemtoDreamControlSample::SetEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles, float mult) {
  const int iMult = FindBin(mult);
  float RelativeK = 0.f;
  int HistCounter = 0;
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  for (auto itSpec1 = Particles.begin(); itSpec1 != Particles.end();
      ++itSpec1) {
    auto itPDGPar2 = fPDGParticleSpecies.begin();
    itPDGPar2 += itSpec1 - Particles.begin();
    for (auto itSpec2 = itSpec1; itSpec2 != Particles.end(); ++itSpec2) {
      fHists->FillPartnersSE(HistCounter, itSpec1->size(), itSpec2->size());
      // Now loop over the actual Particles and correlate them
      for (auto itPart1 = itSpec1->begin(); itPart1 != itSpec1->end();
          ++itPart1) {
        std::vector<AliFemtoDreamBasePart>::iterator itPart2;
        if (itSpec1 == itSpec2) {
          itPart2 = itPart1 + 1;
        } else {
          itPart2 = itSpec2->begin();
        }
        while (itPart2 != itSpec2->end()) {

          // Delta eta - Delta phi* cut
          if (ComputeDeltaEta(*itPart1, *itPart2) < fDeltaEtaMax) {
            ++itPart2;
            continue;
          }
          if (ComputeDeltaPhi(*itPart1, *itPart2) < fDeltaPhiMax) {
            ++itPart2;
            continue;
          }

          // correlated sample
          RelativeK = RelativePairMomentum(itPart1->GetMomentum(), *itPDGPar1,
                                           itPart2->GetMomentum(), *itPDGPar2);
          fHists->FillSameEventDist(HistCounter, RelativeK);
          if (fHists->GetDoMultBinning()) {
            fHists->FillSameEventMultDist(HistCounter, iMult + 1, RelativeK);
          }
          for (int i = 0; i < fSpinningDepth; ++i) {
            // randomized sample - who is the father???
            RelativeK = RelativePairMomentum(itPart1->GetMomentum(), *itPDGPar1,
                                             itPart2->GetMomentum(), *itPDGPar2,
                                             true);
            fHists->FillMixedEventDist(HistCounter, RelativeK);
            if (fHists->GetDoMultBinning()) {
              fHists->FillMixedEventMultDist(HistCounter, iMult + 1, RelativeK);
            }
          }
          ++itPart2;
        }
      }
      ++HistCounter;
      itPDGPar2++;
    }
    itPDGPar1++;
  }
}

float AliFemtoDreamControlSample::RelativePairMomentum(TVector3 Part1Momentum,
                                                       int PDGPart1,
                                                       TVector3 Part2Momentum,
                                                       int PDGPart2,
                                                       bool random) {
  if (PDGPart1 == 0 || PDGPart2 == 0) {
    AliError("Invalid PDG Code");
  }
  float results = 0.;
  TLorentzVector SPtrack, TPProng, trackSum, SPtrackCMS, TPProngCMS;
  // Even if the Daughter tracks were switched up during PID doesn't play a role
  // here cause we are
  // only looking at the mother mass
  SPtrack.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(), Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  TPProng.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(), Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());

  if (random) {
    // Do the randomization here
    if (fStravinsky) {
      if (fRandom.Uniform() < 0.5) {
        SPtrack.SetPhi(SPtrack.Phi() + fPi);
      } else {
        TPProng.SetPhi(TPProng.Phi() + fPi);
      }
    } else {
      SPtrack.SetPhi(SPtrack.Phi() + fRandom.Uniform(2 * fPi));
      TPProng.SetPhi(TPProng.Phi() + fRandom.Uniform(2 * fPi));
    }
  }
  trackSum = SPtrack + TPProng;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  SPtrackCMS = SPtrack;
  TPProngCMS = TPProng;

  SPtrackCMS.Boost(-betax, -betay, -betaz);
  TPProngCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK;

  trackRelK = SPtrackCMS - TPProngCMS;
  results = 0.5 * trackRelK.P();
  return results;
}

int AliFemtoDreamControlSample::FindBin(float Multiplicity) {
  int binCounter = fMultBins.size();
  for (std::vector<int>::reverse_iterator itBin = fMultBins.rbegin();
      itBin != fMultBins.rend(); ++itBin) {
    binCounter--;
    if (Multiplicity >= *itBin) {
      break;
    }
  }
  return binCounter;
}

float AliFemtoDreamControlSample::ComputeDeltaEta(
    AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2) {
  float eta1 = part1.GetEta().at(0);
  float eta2 = part2.GetEta().at(0);
  return std::abs(eta1 - eta2);
}

float AliFemtoDreamControlSample::ComputeDeltaPhi(
    AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2) {
  std::vector<float> Phirad1 = part1.GetPhiAtRaidius().at(0);
  std::vector<float> Phirad2 = part2.GetPhiAtRaidius().at(0);
  std::vector<float> radVector;
  float dphi = 999.f;
  for (int iRad = 0; iRad < 9; ++iRad) {
    float currentdphi = std::abs(Phirad1.at(iRad) - Phirad2.at(iRad));
    if(currentdphi < dphi) dphi = currentdphi;
  }
  return dphi;
}
