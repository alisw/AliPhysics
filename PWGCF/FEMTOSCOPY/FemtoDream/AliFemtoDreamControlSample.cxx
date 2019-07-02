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
      fmode(kNone),
      fSpinningDepth(0),
      fCorrelationRange(0.),
      fDeltaEtaMax(0.f),
      fDeltaPhiMax(0.f),
      fDoDeltaEtaDeltaPhiCut(false) {
  fRandom.SetSeed(0);
}

AliFemtoDreamControlSample::AliFemtoDreamControlSample(
    const AliFemtoDreamControlSample &samp)
    : fHists(samp.fHists),
      fPDGParticleSpecies(samp.fPDGParticleSpecies),
      fMultBins(samp.fMultBins),
      fRandom(),
      fPi(TMath::Pi()),
      fmode(samp.fmode),
      fSpinningDepth(samp.fSpinningDepth),
      fCorrelationRange(samp.fCorrelationRange),
      fDeltaEtaMax(samp.fDeltaEtaMax),
      fDeltaPhiMax(samp.fDeltaPhiMax),
      fDoDeltaEtaDeltaPhiCut(samp.fDoDeltaEtaDeltaPhiCut) {
  fRandom.SetSeed(0);
}

AliFemtoDreamControlSample::AliFemtoDreamControlSample(
    AliFemtoDreamCollConfig *conf, bool minimalBooking)
    : fHists(new AliFemtoDreamCorrHists(conf, minimalBooking)),
      fPDGParticleSpecies(conf->GetPDGCodes()),
      fMultBins(conf->GetMultBins()),
      fRandom(),
      fPi(TMath::Pi()),
      fmode(kNone),
      fSpinningDepth(0),
      fCorrelationRange(0.),
      fDeltaEtaMax(conf->GetDeltaEtaMax()),
      fDeltaPhiMax(conf->GetDeltaPhiMax()),
      fDoDeltaEtaDeltaPhiCut(conf->GetDoDeltaEtaDeltaPhiCut()) {
//  if (fStravinsky) {
//    fSpinningDepth = 1;
//  } else {
//    fSpinningDepth = conf->GetSpinningDepth();
//  }
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
  this->fDoDeltaEtaDeltaPhiCut = samp.fDoDeltaEtaDeltaPhiCut;
  return *this;
}

AliFemtoDreamControlSample::~AliFemtoDreamControlSample() {
}

void AliFemtoDreamControlSample::SetEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles, float mult) {
  const int iMult = FindBin(mult);
  int HistCounter = 0;
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  auto itSpec1 = Particles.begin();
  bool sameParticle = false;
  while (itSpec1 != Particles.end()) {
    auto itSpec2 = itSpec1;
    auto itPDGPar2 = itPDGPar1;
    while (itSpec2 != Particles.end()) {
      fHists->FillPartnersSE(HistCounter, itSpec1->size(), itSpec2->size());
      if ((itSpec1->size() + itSpec2->size()) > 1) {  // for single particle pairs, this is pointless
        if (itSpec1 == itSpec2) {
          sameParticle = true;
        } else {
          sameParticle = false;
        }
        CorrelatedSample(*itSpec1, *itPDGPar1, *itSpec2, *itPDGPar1,
                         sameParticle, iMult, HistCounter);
        UncorrelatedSample(*itSpec1, *itPDGPar1, *itSpec2, *itPDGPar1,
                           sameParticle, iMult, HistCounter);
      }
      ++itSpec2;
      ++HistCounter;
      ++itPDGPar2;
    }
    itPDGPar1++;
    itSpec1++;
  }
}

float AliFemtoDreamControlSample::RelativePairMomentum(TVector3 Part1Momentum,
                                                       int PDGPart1,
                                                       TVector3 Part2Momentum,
                                                       int PDGPart2,
                                                       UncorrelatedMode mode) {
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
  // Do the randomization here
  if (mode == kStravinsky) {
    if (fRandom.Uniform() < 0.5) {
      SPtrack.SetPhi(SPtrack.Phi() + fPi);
    } else {
      TPProng.SetPhi(TPProng.Phi() + fPi);
    }
  } else if (mode == kStravinsky) {
    SPtrack.SetPhi(SPtrack.Phi() + fRandom.Uniform(2 * fPi));
    TPProng.SetPhi(TPProng.Phi() + fRandom.Uniform(2 * fPi));
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

void AliFemtoDreamControlSample::CorrelatedSample(
    std::vector<AliFemtoDreamBasePart> &part1, int &PDGPart1,
    std::vector<AliFemtoDreamBasePart> &part2, int &PDGPart2, bool SameParticle,
    int Mult, int HistCounter) {
  float RelativeK = 0;
  auto itPart1 = part1.begin();
  while (itPart1 != part1.end()) {
    auto itPart2 = SameParticle ? itPart1 + 1 : part2.begin();
    while (itPart2 != part2.end()) {
      if (fDoDeltaEtaDeltaPhiCut) {
        if (ComputeDeltaEta(*itPart1, *itPart2) < fDeltaEtaMax) {
          ++itPart2;
          continue;
        }
        if (ComputeDeltaPhi(*itPart1, *itPart2) < fDeltaPhiMax) {
          ++itPart2;
          continue;
        }
        RelativeK = RelativePairMomentum(itPart1->GetMomentum(), PDGPart1,
                                         itPart2->GetMomentum(), PDGPart2);
        fHists->FillSameEventDist(HistCounter, RelativeK);
        if (fHists->GetDoMultBinning()) {
          fHists->FillSameEventMultDist(HistCounter, Mult + 1, RelativeK);
        }
      }
      itPart2++;
    }
    itPart1++;
  }
}

void AliFemtoDreamControlSample::UncorrelatedSample(
    std::vector<AliFemtoDreamBasePart> &part1, int &PDGPart1,
    std::vector<AliFemtoDreamBasePart> &part2, int &PDGPart2, bool SameParticle,
    int Mult, int HistCounter) {
  if (fmode == kPhiSpin || fmode == kStravinsky) {
    PhiSpinning(part1, PDGPart1, part2, PDGPart2, SameParticle,Mult, HistCounter);
  } else if (fmode == kCorrelatedPhi) {
    for (int iSpin = 0; iSpin < fSpinningDepth; ++iSpin) {
      LimitedPhiSpinning(part1, PDGPart1, part2, PDGPart2, SameParticle,Mult, HistCounter);
    }
  } else {
    Error("AliFemtoDreamControlSample::UncorrelatedSample",
          "None implemented mode\n");
  }
  return;
}

void AliFemtoDreamControlSample::PhiSpinning(
    std::vector<AliFemtoDreamBasePart> &part1, int PDGPart1,
    std::vector<AliFemtoDreamBasePart> &part2, int PDGPart2, bool SameParticle,
    int Mult, int HistCounter) {
  float RelativeK = 0;
  auto itPart1 = part1.begin();
  while (itPart1 != part1.end()) {
    auto itPart2 = SameParticle ? itPart1 + 1 : part2.begin();
    while (itPart2 != part2.end()) {
      if (fDoDeltaEtaDeltaPhiCut) {
        if (ComputeDeltaEta(*itPart1, *itPart2) < fDeltaEtaMax) {
          ++itPart2;
          continue;
        }
        if (ComputeDeltaPhi(*itPart1, *itPart2) < fDeltaPhiMax) {
          ++itPart2;
          continue;
        }
        for (int i = 0; i < fSpinningDepth; ++i) {
          // randomized sample - who is the father???
          RelativeK = RelativePairMomentum(itPart1->GetMomentum(), PDGPart1,
                                           itPart2->GetMomentum(), PDGPart2,
                                           fmode);
          fHists->FillMixedEventDist(HistCounter, RelativeK);
          if (fHists->GetDoMultBinning()) {
            fHists->FillMixedEventMultDist(HistCounter, Mult + 1, RelativeK);
          }
        }
      }
      itPart2++;
    }
    itPart1++;
  }
}

void AliFemtoDreamControlSample::LimitedPhiSpinning(
    std::vector<AliFemtoDreamBasePart> &part1, int PDGPart1,
    std::vector<AliFemtoDreamBasePart> &part2, int PDGPart2, bool SameParticle,
    int Mult, int HistCounter) {
  //copy input to not touch the initial state & make it to pointers to do some nasty things.
  std::vector<AliFemtoDreamBasePart*> CopyPart1;
  std::vector<AliFemtoDreamBasePart*> CopyPart2;

  for (auto it : part1) {
    CopyPart1.emplace_back(new AliFemtoDreamBasePart(it));
  }
  // randomize the sample: add references to all particles to a vector
  std::vector<AliFemtoDreamBasePart*> RandomizeMe;
  for (auto it : CopyPart1) {
    RandomizeMe.push_back(it);
  }
  // Only add the second particle, in case it is a different particle
  if (!SameParticle) {
    for (auto it : part2) {
      CopyPart2.emplace_back(new AliFemtoDreamBasePart(it));
    }
    for (auto it : CopyPart2) {
      RandomizeMe.push_back(it);
    }
  } else {  // if the particles are the same this should point to the same vector as particle 1
    CopyPart2 = CopyPart1;
  }
  Randomizer(RandomizeMe);
  // calculate the relative momentum
  float RelativeK = 0;
  auto itPart1 = CopyPart1.begin();
  while (itPart1 != CopyPart1.end()) {
    auto itPart2 = SameParticle ? itPart1 + 1 : CopyPart2.begin();
    while (itPart2 != CopyPart2.end()) {
      if (fDoDeltaEtaDeltaPhiCut) {
        if (ComputeDeltaEta(*itPart1, *itPart2) < fDeltaEtaMax) {
          ++itPart2;
          continue;
        }
        if (ComputeDeltaPhi(*itPart1, *itPart2) < fDeltaPhiMax) {
          ++itPart2;
          continue;
        }
        for (int i = 0; i < fSpinningDepth; ++i) {
          // randomized sample - who is the father???
          RelativeK = RelativePairMomentum((*itPart1)->GetMomentum(), PDGPart1,
                                           (*itPart2)->GetMomentum(), PDGPart2,
                                           kNone);
          fHists->FillMixedEventDist(HistCounter, RelativeK);
          if (fHists->GetDoMultBinning()) {
            fHists->FillMixedEventMultDist(HistCounter, Mult + 1, RelativeK);
          }
        }
      }
      itPart2++;
    }
    itPart1++;
  }
}

void AliFemtoDreamControlSample::Randomizer(
    std::vector<AliFemtoDreamBasePart*> &part) {
  if (part.size() < 2) {
    Warning("AliFemtoDreamControlSample::Randomizer",
            "this should not happen, check the vector size\n");
  }
  auto itPart1 = part.begin();
  auto itPart2 = itPart1;  // dummy vale to start the while
  while (itPart1 != part.end()) {  // this doersn't work
    itPart2 = itPart1 + 1;  //pick the next particle
    if (itPart2 == part.end()) {  // if this is out of range, finish the while
      break;
    }
    // Get the momenta
    TVector3 Part1Mom = (*itPart1)->GetMomentum();
    TVector3 Part2Mom = (*itPart2)->GetMomentum();
    Part2Mom += Part1Mom;
    // give the first particle a phi kick
    Part1Mom.SetPhi(
        Part2Mom.Phi()
            + fRandom.Uniform(-fCorrelationRange, fCorrelationRange));
    Part2Mom -= Part1Mom;
    (*itPart1)->SetMomentum(Part1Mom);
    (*itPart2)->SetMomentum(Part2Mom);
    itPart2 += 2;
  }
}

float AliFemtoDreamControlSample::ComputeDeltaEta(
    AliFemtoDreamBasePart *part1, AliFemtoDreamBasePart *part2) {
  float eta1 = part1->GetEta().at(0);
  float eta2 = part2->GetEta().at(0);
  return std::abs(eta1 - eta2);
}

float AliFemtoDreamControlSample::ComputeDeltaPhi(
    AliFemtoDreamBasePart *part1, AliFemtoDreamBasePart *part2) {
  std::vector<float> Phirad1 = part1->GetPhiAtRaidius().at(0);
  std::vector<float> Phirad2 = part2->GetPhiAtRaidius().at(0);
  std::vector<float> radVector;
  float dphi = 999.f;
  for (unsigned int iRad = 0; iRad < Phirad1.size(); ++iRad) {
    float currentdphi = std::abs(Phirad1.at(iRad) - Phirad2.at(iRad));
    if (currentdphi < dphi)
      dphi = currentdphi;
  }
  return dphi;
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
  for (unsigned int iRad = 0; iRad < Phirad1.size(); ++iRad) {
    float currentdphi = std::abs(Phirad1.at(iRad) - Phirad2.at(iRad));
    if (currentdphi < dphi)
      dphi = currentdphi;
  }
  return dphi;
}
