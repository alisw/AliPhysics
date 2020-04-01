#include "AliFemtoDreamControlSample.h"
#include "AliLog.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TMath.h"

ClassImp(AliFemtoDreamControlSample)

AliFemtoDreamControlSample::AliFemtoDreamControlSample()
    : fHigherMath(nullptr),
      fPDGParticleSpecies(),
      fRejPairs(),
      fMultBins(),
      fRandom(),
      fPi(TMath::Pi()),
      fmode(AliFemtoDreamCollConfig::kNone),
      fSpinningDepth(0),
      fCorrelationRange(0.),
      fDeltaEtaMax(0.f),
      fDeltaPhiMax(0.f),
      fDoDeltaEtaDeltaPhiCut(false),
      fMult(0),
      fCent(0) {
  fRandom.SetSeed(0);
}

AliFemtoDreamControlSample::AliFemtoDreamControlSample(
    const AliFemtoDreamControlSample &samp)
    : fHigherMath(samp.fHigherMath),
      fPDGParticleSpecies(samp.fPDGParticleSpecies),
      fRejPairs(samp.fRejPairs),
      fMultBins(samp.fMultBins),
      fRandom(),
      fPi(TMath::Pi()),
      fmode(samp.fmode),
      fSpinningDepth(samp.fSpinningDepth),
      fCorrelationRange(samp.fCorrelationRange),
      fDeltaEtaMax(samp.fDeltaEtaMax),
      fDeltaPhiMax(samp.fDeltaPhiMax),
      fDoDeltaEtaDeltaPhiCut(samp.fDoDeltaEtaDeltaPhiCut),
      fMult(samp.fMult),
      fCent(samp.fCent) {
  fRandom.SetSeed(0);
}

AliFemtoDreamControlSample::AliFemtoDreamControlSample(
    AliFemtoDreamCollConfig *conf)
    : fHigherMath(new AliFemtoDreamHigherPairMath(conf, conf->GetMinimalBookingSample())),
      fPDGParticleSpecies(conf->GetPDGCodes()),
      fRejPairs(conf->GetClosePairRej()),
      fMultBins(conf->GetMultBins()),
      fRandom(),
      fPi(TMath::Pi()),
      fmode(conf->GetControlMode()),
      fSpinningDepth(conf->GetSpinningDepth()),
      fCorrelationRange(conf->GetCorrelationRange()),
      fDeltaEtaMax(conf->GetDeltaEtaMax()),
      fDeltaPhiMax(conf->GetDeltaPhiMax()),
      fDoDeltaEtaDeltaPhiCut(false),
      fMult(0),
      fCent(0) {
  fRandom.SetSeed(0);
}

AliFemtoDreamControlSample& AliFemtoDreamControlSample::operator=(
    const AliFemtoDreamControlSample& samp) {
  if (this == &samp) {
    return *this;
  }
  this->fHigherMath = samp.fHigherMath;
  this->fPDGParticleSpecies = samp.fPDGParticleSpecies;
  this->fRejPairs = samp.fRejPairs;
  this->fMultBins = samp.fMultBins;
  this->fRandom.SetSeed(0);
  this->fPi = TMath::Pi();
  this->fSpinningDepth = samp.fSpinningDepth;
  this->fDeltaEtaMax = samp.fDeltaEtaMax;
  this->fDeltaPhiMax = samp.fDeltaPhiMax;
  this->fDoDeltaEtaDeltaPhiCut = samp.fDoDeltaEtaDeltaPhiCut;
  this->fMult = samp.fMult;
  this->fCent = samp.fCent;
  return *this;
}

AliFemtoDreamControlSample::~AliFemtoDreamControlSample() {
}

void AliFemtoDreamControlSample::SetEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
    AliFemtoDreamEvent *evt) {
  fMult = FindBin(evt->GetMultiplicity());
  fCent = evt->GetV0MCentrality();
  fHigherMath->SetBField(evt->GetMagneticField());
  int HistCounter = 0;
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  auto itSpec1 = Particles.begin();
  bool sameParticle = false;
  size_t minimumSize;
  while (itSpec1 != Particles.end()) {
    auto itSpec2 = itSpec1;
    auto itPDGPar2 = itPDGPar1;
    while (itSpec2 != Particles.end()) {
      fHigherMath->FillPairCounterSE(HistCounter, itSpec1->size(),
                                     itSpec2->size());
      if (itSpec1 == itSpec2) {
        sameParticle = true;
        minimumSize = 1;
      } else {
        minimumSize = 0;
        sameParticle = false;
      }
      //if it is the same particle, then we need at least two to pair them 2*2 is the minimum
      if ((itSpec1->size() > minimumSize) && (itSpec2->size() > minimumSize)) {  // for single particle pairs, this is pointless
        CorrelatedSample(*itSpec1, *itPDGPar1, *itSpec2, *itPDGPar2,
                         sameParticle, HistCounter);
        UncorrelatedSample(*itSpec1, *itPDGPar1, *itSpec2, *itPDGPar2,
                           sameParticle, HistCounter);
      }
      ++itSpec2;
      ++HistCounter;
      ++itPDGPar2;
    }
    itPDGPar1++;
    itSpec1++;
  }
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
    int HistCounter) {
  float RelativeK = 0;
  auto itPart1 = part1.begin();
  bool CPR = fRejPairs.at(HistCounter);
  while (itPart1 != part1.end()) {
    auto itPart2 = SameParticle ? itPart1 + 1 : part2.begin();
    while (itPart2 != part2.end()) {
      TLorentzVector PartOne, PartTwo;
      PartOne.SetXYZM(
          itPart1->GetMomentum().X(), itPart1->GetMomentum().Y(),
          itPart1->GetMomentum().Z(),
          TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
      PartTwo.SetXYZM(
          itPart2->GetMomentum().X(), itPart2->GetMomentum().Y(),
          itPart2->GetMomentum().Z(),
          TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
      float RelativeK = fHigherMath->RelativePairMomentum(PartOne, PartTwo);
      if (!fHigherMath->PassesPairSelection(HistCounter, *itPart1, *itPart2,
                                           RelativeK, true, false)) {
        ++itPart2;
        continue;
      }
      RelativeK = fHigherMath->FillSameEvent(HistCounter, fMult, fCent,
                                             *itPart1, PDGPart1,
                                             *itPart2, PDGPart2);
      fHigherMath->MassQA(HistCounter, RelativeK, *itPart1, *itPart2);
      fHigherMath->SEDetaDPhiPlots(HistCounter, *itPart1, PDGPart1, *itPart2,
                                   PDGPart2, RelativeK, true);

      itPart2++;
    }
    itPart1++;
  }
}

void AliFemtoDreamControlSample::UncorrelatedSample(
    std::vector<AliFemtoDreamBasePart> &part1, int &PDGPart1,
    std::vector<AliFemtoDreamBasePart> &part2, int &PDGPart2, bool SameParticle,
    int HistCounter) {
  if (fmode == AliFemtoDreamCollConfig::kPhiSpin
      || fmode == AliFemtoDreamCollConfig::kStravinsky) {
    PhiSpinning(part1, PDGPart1, part2, PDGPart2, SameParticle, HistCounter);
  } else if (fmode == AliFemtoDreamCollConfig::kCorrelatedPhi) {
    for (int iSpin = 0; iSpin < fSpinningDepth; ++iSpin) {
      LimitedPhiSpinning(part1, PDGPart1, part2, PDGPart2, SameParticle,
                         HistCounter);
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
    int HistCounter) {
  float RelativeK = 0;
  auto itPart1 = part1.begin();
  bool CPR = fRejPairs.at(HistCounter);
  while (itPart1 != part1.end()) {
    auto itPart2 = SameParticle ? itPart1 + 1 : part2.begin();
    while (itPart2 != part2.end()) {
      //in this case it does NOT work as intended since the new phi of the particle has to
      //be considered for this cut
      TLorentzVector PartOne, PartTwo;
      PartOne.SetXYZM(
          itPart1->GetMomentum().X(), itPart1->GetMomentum().Y(),
          itPart1->GetMomentum().Z(),
          TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
      PartTwo.SetXYZM(
          itPart2->GetMomentum().X(), itPart2->GetMomentum().Y(),
          itPart2->GetMomentum().Z(),
          TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
      float RelativeK = fHigherMath->RelativePairMomentum(PartOne, PartTwo);
      if (!fHigherMath->PassesPairSelection(HistCounter, *itPart1, *itPart2,
                                            RelativeK, false, true)) {
        ++itPart2;
        continue;
      }
      for (int i = 0; i < fSpinningDepth; ++i) {
        // randomized sample - who is the father???
        RelativeK = fHigherMath->FillMixedEvent(HistCounter, fMult, fCent,
                                                *itPart1,
                                                PDGPart1,
                                                *itPart2,
                                                PDGPart2, fmode);
        fHigherMath->MEDetaDPhiPlots(HistCounter, *itPart1, PDGPart1, *itPart2,
                                     PDGPart2, RelativeK, true);
      }
      itPart2++;
    }
    itPart1++;
  }
}

void AliFemtoDreamControlSample::LimitedPhiSpinning(
    std::vector<AliFemtoDreamBasePart> &part1, int PDGPart1,
    std::vector<AliFemtoDreamBasePart> &part2, int PDGPart2, bool SameParticle,
    int HistCounter) {
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
  bool CPR = fRejPairs.at(HistCounter);
  auto itPart1 = CopyPart1.begin();
  while (itPart1 != CopyPart1.end()) {
    auto itPart2 = SameParticle ? itPart1 + 1 : CopyPart2.begin();
    while (itPart2 != (SameParticle ? CopyPart1 : CopyPart2).end()) {
      //in this case it does NOT work as intended since the new phi of the particle has to
      //be considered for this cut
      TLorentzVector PartOne, PartTwo;
      PartOne.SetXYZM(
          (*itPart1)->GetMomentum().X(), (*itPart1)->GetMomentum().Y(),
          (*itPart1)->GetMomentum().Z(),
          TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
      PartTwo.SetXYZM(
          (*itPart2)->GetMomentum().X(), (*itPart2)->GetMomentum().Y(),
          (*itPart2)->GetMomentum().Z(),
          TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
      float RelativeK = fHigherMath->RelativePairMomentum(PartOne, PartTwo);
      if (!fHigherMath->PassesPairSelection(HistCounter, *(*itPart2), *(*itPart2),
                                            RelativeK, false, true)) {
        ++itPart2;
        continue;
      }
      // randomized sample - who is the father???
      RelativeK = fHigherMath->FillMixedEvent(HistCounter, fMult, fCent,
                                              *(*itPart1),
                                              PDGPart1,
                                              *(*itPart2),
                                              PDGPart2,
                                              AliFemtoDreamCollConfig::kNone);
      fHigherMath->MEDetaDPhiPlots(HistCounter, *(*itPart1), PDGPart1, *(*itPart2),
                                   PDGPart2, RelativeK, true);
      itPart2++;
    }
    itPart1++;
  }
  for (auto iDel : RandomizeMe) {
    delete iDel;
  }
  RandomizeMe.clear();
//  it was checked with massif that there is no memory leak and even if the pointer and its members
//  still exist.
//  for (auto iTest : CopyPart1) {
//    if (iTest){
//      Warning("LimitedPhiSpinning",
//              "This should be deleted, think about your code \n");
//    }
//  }
  CopyPart1.clear();
  CopyPart2.clear();
  return;
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
        Part1Mom.Phi()
            + fRandom.Uniform(-fCorrelationRange, fCorrelationRange));
    Part2Mom -= Part1Mom;
    (*itPart1)->SetMomentum(0, Part1Mom);
    (*itPart2)->SetMomentum(0, Part2Mom);
    itPart1 += 2;
  }
}
