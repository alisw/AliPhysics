/*
 * AliFemtoDreamHigherPairMath.cxx
 *
 *  Created on: Jul 3, 2019
 *      Author: schmollweger
 */

#include <AliFemtoDreamHigherPairMath.h>
#include "TMath.h"
#include "TDatabasePDG.h"
static const float piHi = TMath::Pi();

AliFemtoDreamHigherPairMath::AliFemtoDreamHigherPairMath(
    AliFemtoDreamCollConfig *conf, bool minBooking)
    : fHists(new AliFemtoDreamCorrHists(conf, minBooking)),
      fWhichPairs(conf->GetWhichPairs()),
      fBField(-99.),
      fRejPairs(conf->GetClosePairRej()),
      fDoDeltaEtaDeltaPhiCut(false),
      fDeltaPhiEtaMax(conf->GetSqDeltaPhiEtaMax()),
      fRandom(),
      fPi(TMath::Pi()) {
  fRandom.SetSeed(0);
  TDatabasePDG::Instance()->AddParticle("deuteron", "deuteron", 1.8756134,
                                        kTRUE, 0.0, 1, "Nucleus", 1000010020);
  TDatabasePDG::Instance()->AddAntiParticle("anti-deuteron", -1000010020);
}

AliFemtoDreamHigherPairMath::~AliFemtoDreamHigherPairMath() {
}

AliFemtoDreamHigherPairMath::AliFemtoDreamHigherPairMath(
    const AliFemtoDreamHigherPairMath& samp)
    : fHists(samp.fHists),
      fWhichPairs(samp.fWhichPairs),
      fBField(-99.),
      fRejPairs(samp.fRejPairs),
      fDoDeltaEtaDeltaPhiCut(samp.fRejPairs),
      fDeltaPhiEtaMax(samp.fDeltaPhiEtaMax),
      fRandom(),
      fPi(TMath::Pi()) {
  fRandom.SetSeed(0);
}
AliFemtoDreamHigherPairMath& AliFemtoDreamHigherPairMath::operator=(
    const AliFemtoDreamHigherPairMath& math) {
  if (this == &math) {
    return *this;
  }
  fHists = math.fHists;
  fWhichPairs = math.fWhichPairs;
  fBField = math.fBField;
  fDeltaPhiEtaMax = math.fDeltaPhiEtaMax;
  fRandom = math.fRandom;
  fPi = TMath::Pi();
  return *this;
}

bool AliFemtoDreamHigherPairMath::PassesPairSelection(
    int iHC, AliFemtoDreamBasePart& part1, AliFemtoDreamBasePart& part2,
    float RelativeK, bool SEorME, bool Recalculate) {
  bool outBool = true;
  //Method calculates the average separation between two tracks
  //at different radii within the TPC and rejects pairs which a
  //too low separation
  if (Recalculate) {
    RecalculatePhiStar(part1);
    RecalculatePhiStar(part2);
  }
  bool pass = true;
  bool CPR = fRejPairs.at(iHC);

  //only need to do this, if we wanna fill the plots or wanna do the CPR ...
  if ((CPR && fDoDeltaEtaDeltaPhiCut) || fHists->GetEtaPhiPlots()) {
    pass = DeltaEtaDeltaPhi(iHC, part1, part2, SEorME, RelativeK);
  }
  return pass;
}

void AliFemtoDreamHigherPairMath::RecalculatePhiStar(
    AliFemtoDreamBasePart &part) {
//Only use this for Tracks, this was implemented for the case where particles
//were added a random phi to obtain an uncorrelated sample. This should be extended
//to the daughter tracks. For some reason the momentum is stored in a TVector3,
//and not a std::vector<TVector3>. Should be changed in time, conceptually the
//difficulty is to make this work in a proper way after the mother,
//e.g. the Lambda was phi shifted.
//
  if (fBField < -95) {
    AliWarning(
        "BField was most probably not set! PhiStar Calculation meaningless. \n");
  }
  part.ResizePhiAtRadii(0);
  static float TPCradii[9] = { 85., 105., 125., 145., 165., 185., 205., 225.,
      245. };
  std::vector<float> tmpVec;
  auto phi0 = part.GetMomentum().Phi();
  float pt = part.GetMomentum().Pt();
  float chg = part.GetCharge().at(0);
  for (int radius = 0; radius < 9; radius++) {
    tmpVec.push_back(
        phi0
            - TMath::ASin(
                0.1 * chg * fBField * 0.3 * TPCradii[radius] * 0.01
                    / (2. * pt)));
  }
  part.SetPhiAtRadius(tmpVec);
}

float AliFemtoDreamHigherPairMath::FillSameEvent(int iHC, int Mult, float cent,
                                                 TVector3 Part1Momentum,
                                                 int PDGPart1,
                                                 TVector3 Part2Momentum,
                                                 int PDGPart2) {
  if (PDGPart1 == 0 || PDGPart2 == 0) {
    AliError("Invalid PDG Code");
  }
  bool fillHists = fWhichPairs.at(iHC);
  TLorentzVector PartOne, PartTwo;
// Even if the Daughter tracks were switched up during PID doesn't play a role
// here cause we are
// only looking at the mother mass
  PartOne.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(), Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  PartTwo.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(), Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());

  float RelativeK = RelativePairMomentum(PartOne, PartTwo);
  fHists->FillSameEventDist(iHC, RelativeK);
  if (fHists->GetDoMultBinning()) {
    fHists->FillSameEventMultDist(iHC, Mult + 1, RelativeK);
  }
  if (fillHists && fHists->GetDoCentBinning()) {
    fHists->FillSameEventCentDist(iHC, cent, RelativeK);
  }
  if (fillHists && fHists->GetDokTBinning()) {
    fHists->FillSameEventkTDist(iHC, RelativePairkT(PartOne, PartTwo),
                                RelativeK, cent);
  }
  if (fillHists && fHists->GetDomTBinning()) {
    fHists->FillSameEventmTDist(iHC, RelativePairmT(PartOne, PartTwo),
                                RelativeK);
  }
  if (fillHists && fHists->GetDokTandMultBinning()) {
    fHists->FillSameEventkTandMultDist(iHC, RelativePairkT(PartOne, PartTwo),
                                       RelativeK, Mult + 1);
  }
  if (fillHists && fHists->GetDoPtQA()) {
    fHists->FillPtQADist(iHC, RelativeK, Part1Momentum.Pt(),
                         Part2Momentum.Pt());
    fHists->FillPtSEOneQADist(iHC, Part1Momentum.Pt(), Mult + 1);
    fHists->FillPtSETwoQADist(iHC, Part2Momentum.Pt(), Mult + 1);
  }
  return RelativeK;
}

void AliFemtoDreamHigherPairMath::MassQA(int iHC, float RelK,
                                         AliFemtoDreamBasePart &part1,
                                         AliFemtoDreamBasePart &part2) {
  if (fWhichPairs.at(iHC) && fHists->GetDoMassQA()) {
    fHists->FillMassQADist(iHC, RelK, part1.GetInvMass(), part2.GetInvMass());
    fHists->FillPairInvMassQAD(iHC, part1, part2);
  }
}

float AliFemtoDreamHigherPairMath::FillMixedEvent(
    int iHC, int Mult, float cent, TVector3 Part1Momentum, int PDGPart1,
    TVector3 Part2Momentum, int PDGPart2,
    AliFemtoDreamCollConfig::UncorrelatedMode mode) {
  if (PDGPart1 == 0 || PDGPart2 == 0) {
    AliError("Invalid PDG Code");
  }
  bool fillHists = fWhichPairs.at(iHC);
  TLorentzVector PartOne, PartTwo;
// Even if the Daughter tracks were switched up during PID doesn't play a role
// here cause we are
// only looking at the mother mass
  PartOne.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(), Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  PartTwo.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(), Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
// Do the randomization here
  if (mode == AliFemtoDreamCollConfig::kStravinsky) {
    if (fRandom.Uniform() < 0.5) {
      PartOne.SetPhi(PartOne.Phi() + fPi);
    } else {
      PartTwo.SetPhi(PartTwo.Phi() + fPi);
    }
  } else if (mode == AliFemtoDreamCollConfig::kPhiSpin) {
    PartOne.SetPhi(PartOne.Phi() + fRandom.Uniform(2 * fPi));
    PartTwo.SetPhi(PartTwo.Phi() + fRandom.Uniform(2 * fPi));
  }
  float RelativeK = RelativePairMomentum(PartOne, PartTwo);
  fHists->FillMixedEventDist(iHC, RelativeK);
  if (fHists->GetDoMultBinning()) {
    fHists->FillMixedEventMultDist(iHC, Mult + 1, RelativeK);
  }
  if (fillHists && fHists->GetDoCentBinning()) {
    fHists->FillMixedEventCentDist(iHC, cent, RelativeK);
  }
  if (fillHists && fHists->GetDokTBinning()) {
    fHists->FillMixedEventkTDist(iHC, RelativePairkT(PartOne, PartTwo),
                                 RelativeK, cent);
  }
  if (fillHists && fHists->GetDomTBinning()) {
    fHists->FillMixedEventmTDist(iHC, RelativePairmT(PartOne, PartTwo),
                                 RelativeK);
  }
  if (fillHists && fHists->GetDokTandMultBinning()) {
    fHists->FillMixedEventkTandMultDist(iHC, RelativePairkT(PartOne, PartTwo),
                                        RelativeK, Mult + 1);
  }
  if (fillHists && fHists->GetDoPtQA()) {
    fHists->FillPtMEOneQADist(iHC, Part1Momentum.Pt(), Mult + 1);
    fHists->FillPtMETwoQADist(iHC, Part2Momentum.Pt(), Mult + 1);
  }
  return RelativeK;
}

void AliFemtoDreamHigherPairMath::SEDetaDPhiPlots(int iHC,
                                                  AliFemtoDreamBasePart &part1,
                                                  int PDGPart1,
                                                  AliFemtoDreamBasePart &part2,
                                                  int PDGPart2, float RelativeK,
                                                  bool recalculate) {
  if (fWhichPairs.at(iHC) && fHists->GetDodPhidEtaPlots()) {
    float deta = part1.GetEta().at(0) - part2.GetEta().at(0);
    float dphi = part1.GetPhi().at(0) - part2.GetPhi().at(0);
    float mT = 0;
    if (fHists->GetDodPhidEtamTPlots()) {
      TLorentzVector PartOne, PartTwo;
      PartOne.SetXYZM(part1.GetMCMomentum().X(), part1.GetMCMomentum().Y(),
                      part1.GetMCMomentum().Z(),
                      TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
      PartTwo.SetXYZM(part2.GetMCMomentum().X(), part2.GetMCMomentum().Y(),
                      part2.GetMCMomentum().Z(),
                      TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
      mT = RelativePairmT(PartOne, PartTwo);
    }
    if (dphi < 0) {
      fHists->FilldPhidEtaSE(iHC, dphi + 2 * TMath::Pi(), deta, mT);
    } else {
      fHists->FilldPhidEtaSE(iHC, dphi, deta, mT);
    }
  }
}

void AliFemtoDreamHigherPairMath::MEDetaDPhiPlots(int iHC,
                                                  AliFemtoDreamBasePart &part1,
                                                  int PDGPart1,
                                                  AliFemtoDreamBasePart &part2,
                                                  int PDGPart2, float RelativeK,
                                                  bool recalculate) {
  if (fWhichPairs.at(iHC) && fHists->GetDodPhidEtaPlots()) {
    float deta = part1.GetEta().at(0) - part2.GetEta().at(0);
    float dphi = part1.GetPhi().at(0) - part2.GetPhi().at(0);
    float mT = 0;
    if (fHists->GetDodPhidEtamTPlots()) {
      TLorentzVector PartOne, PartTwo;
      PartOne.SetXYZM(part1.GetMCMomentum().X(), part1.GetMCMomentum().Y(),
                      part1.GetMCMomentum().Z(),
                      TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
      PartTwo.SetXYZM(part2.GetMCMomentum().X(), part2.GetMCMomentum().Y(),
                      part2.GetMCMomentum().Z(),
                      TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
      mT = RelativePairmT(PartOne, PartTwo);
    }
    if (dphi < 0) {
      fHists->FilldPhidEtaME(iHC, dphi + 2 * TMath::Pi(), deta, mT);
    } else {
      fHists->FilldPhidEtaME(iHC, dphi, deta, mT);
    }
  }
}

void AliFemtoDreamHigherPairMath::SEMomentumResolution(
    int iHC, AliFemtoDreamBasePart* part1, int PDGPart1,
    AliFemtoDreamBasePart* part2, int PDGPart2, float RelativeK) {
  if (fWhichPairs[iHC] && fHists->GetObtainMomentumResolution()) {
    //It is sufficient to do this in Mixed events, which allows
    //to increase the statistics. The Resolution of the tracks and therefore
    //of the pairs does not change event by event.
    //Now we only want to use the momentum of particles we are after, hence
    //we check the PDG Code!
    //This should not be used in the rotated sample, since the MC Truth is not rotated
    TLorentzVector PartOne, PartTwo;
    // Even if the Daughter tracks were switched up during PID doesn't play a role
    // here cause we are
    // only looking at the mother mass
    PartOne.SetXYZM(part1->GetMCMomentum().X(), part1->GetMCMomentum().Y(),
                    part1->GetMCMomentum().Z(),
                    TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
    PartTwo.SetXYZM(part2->GetMCMomentum().X(), part2->GetMCMomentum().Y(),
                    part2->GetMCMomentum().Z(),
                    TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
    float RelKTrue = RelativePairMomentum(PartOne, PartTwo);
    fHists->FillMomentumResolutionSEAll(iHC, RelKTrue, RelativeK);
    if ((PDGPart1 == TMath::Abs(part1->GetMCPDGCode()))
        && ((PDGPart2 == TMath::Abs(part2->GetMCPDGCode())))) {
      fHists->FillMomentumResolutionSE(iHC, RelKTrue, RelativeK);
    }
  }
}

void AliFemtoDreamHigherPairMath::MEMomentumResolution(
    int iHC, AliFemtoDreamBasePart* part1, int PDGPart1,
    AliFemtoDreamBasePart* part2, int PDGPart2, float RelativeK) {
  if (fWhichPairs[iHC] && fHists->GetObtainMomentumResolution()) {
    //It is sufficient to do this in Mixed events, which allows
    //to increase the statistics. The Resolution of the tracks and therefore
    //of the pairs does not change event by event.
    //Now we only want to use the momentum of particles we are after, hence
    //we check the PDG Code!
    //This should not be used in the rotated sample, since the MC Truth is not rotated
    TLorentzVector PartOne, PartTwo;
    // Even if the Daughter tracks were switched up during PID doesn't play a role
    // here cause we are
    // only looking at the mother mass
    PartOne.SetXYZM(part1->GetMCMomentum().X(), part1->GetMCMomentum().Y(),
                    part1->GetMCMomentum().Z(),
                    TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
    PartTwo.SetXYZM(part2->GetMCMomentum().X(), part2->GetMCMomentum().Y(),
                    part2->GetMCMomentum().Z(),
                    TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
    float RelKTrue = RelativePairMomentum(PartOne, PartTwo);
    fHists->FillMomentumResolutionMEAll(iHC, RelKTrue, RelativeK);
    if ((PDGPart1 == TMath::Abs(part1->GetMCPDGCode()))
        && ((PDGPart2 == TMath::Abs(part2->GetMCPDGCode())))) {
      fHists->FillMomentumResolutionME(iHC, RelKTrue, RelativeK);
    }
  }
}

float AliFemtoDreamHigherPairMath::RelativePairMomentum(
    AliFemtoDreamBasePart *part1, const int pdg1, AliFemtoDreamBasePart *part2,
    const int pdg2) {
  TLorentzVector PartOne, PartTwo;
  PartOne.SetXYZM(part1->GetMomentum().X(), part1->GetMomentum().Y(),
                  part1->GetMomentum().Z(),
                  TDatabasePDG::Instance()->GetParticle(pdg1)->Mass());
  PartTwo.SetXYZM(part2->GetMomentum().X(), part2->GetMomentum().Y(),
                  part2->GetMomentum().Z(),
                  TDatabasePDG::Instance()->GetParticle(pdg2)->Mass());
  return RelativePairMomentum(PartOne, PartTwo);
}

float AliFemtoDreamHigherPairMath::RelativePairMomentum(
    TLorentzVector &PartOne, TLorentzVector &PartTwo) {
  TLorentzVector trackSum = PartOne + PartTwo;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  TLorentzVector PartOneCMS = PartOne;
  TLorentzVector PartTwoCMS = PartTwo;

  PartOneCMS.Boost(-betax, -betay, -betaz);
  PartTwoCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;

  return 0.5 * trackRelK.P();
}

float AliFemtoDreamHigherPairMath::RelativePairkT(AliFemtoDreamBasePart *part1,
                                                  const int pdg1,
                                                  AliFemtoDreamBasePart *part2,
                                                  const int pdg2) {
  TLorentzVector PartOne, PartTwo;
  PartOne.SetXYZM(part1->GetMomentum().X(), part1->GetMomentum().Y(),
                  part1->GetMomentum().Z(),
                  TDatabasePDG::Instance()->GetParticle(pdg1)->Mass());
  PartTwo.SetXYZM(part2->GetMomentum().X(), part2->GetMomentum().Y(),
                  part2->GetMomentum().Z(),
                  TDatabasePDG::Instance()->GetParticle(pdg2)->Mass());
  return RelativePairkT(PartOne, PartTwo);
}

float AliFemtoDreamHigherPairMath::RelativePairkT(TLorentzVector &PartOne,
                                                  TLorentzVector &PartTwo) {
  float results = 0.;
  TLorentzVector trackSum;
  trackSum = PartOne + PartTwo;
  results = 0.5 * trackSum.Pt();
  return results;
}

float AliFemtoDreamHigherPairMath::RelativePairmT(AliFemtoDreamBasePart *part1,
                                                  const int pdg1,
                                                  AliFemtoDreamBasePart *part2,
                                                  const int pdg2) {
  TLorentzVector PartOne, PartTwo;
  PartOne.SetXYZM(part1->GetMomentum().X(), part1->GetMomentum().Y(),
                  part1->GetMomentum().Z(),
                  TDatabasePDG::Instance()->GetParticle(pdg1)->Mass());
  PartTwo.SetXYZM(part2->GetMomentum().X(), part2->GetMomentum().Y(),
                  part2->GetMomentum().Z(),
                  TDatabasePDG::Instance()->GetParticle(pdg2)->Mass());
  return RelativePairmT(PartOne, PartTwo);
}

float AliFemtoDreamHigherPairMath::RelativePairmT(TLorentzVector &PartOne,
                                                  TLorentzVector &PartTwo) {
  float results = 0.;
  TLorentzVector trackSum;
  trackSum = PartOne + PartTwo;
  float pairKT = 0.5 * trackSum.Pt();
  float averageMass = 0.5 * (PartOne.M() + PartTwo.M());
  results = TMath::Sqrt(pow(pairKT, 2.) + pow(averageMass, 2.));
  return results;
}

bool AliFemtoDreamHigherPairMath::DeltaEtaDeltaPhi(int Hist,
                                                   AliFemtoDreamBasePart &part1,
                                                   AliFemtoDreamBasePart &part2,
                                                   bool SEorME, float relk) {
  bool pass = true;
  // if nDaug == 1 => Single Track, else decay
  unsigned int DoThisPair = fWhichPairs.at(Hist);
  unsigned int nDaug1 = (unsigned int) DoThisPair / 10;
  if (nDaug1 > 9) {
    AliWarning("you are doing something wrong \n");
  }
  if (nDaug1 != part1.GetPhiAtRaidius().size()) {
    TString outMessage = TString::Format(
        "Your number of Daughters 1 (%u) and Radii 1 (%u) do not correspond \n",
        nDaug1, part1.GetPhiAtRaidius().size());
    AliWarning(outMessage.Data());
  }
  unsigned int nDaug2 = (unsigned int) DoThisPair % 10;

  if (nDaug2 != part2.GetPhiAtRaidius().size()) {
    TString outMessage = TString::Format(
        "Your number of Daughters 2 (%u) and Radii 2 (%u) do not correspond \n",
        nDaug2, part2.GetPhiAtRaidius().size());
    AliWarning(outMessage.Data());
  }
  std::vector<float> eta1 = part1.GetEta();
  std::vector<float> eta2 = part2.GetEta();

  for (unsigned int iDaug1 = 0; iDaug1 < nDaug1; ++iDaug1) {
    std::vector<float> PhiAtRad1 = part1.GetPhiAtRaidius().at(iDaug1);
    float etaPar1;
    if (nDaug1 == 1) {
      etaPar1 = eta1.at(0);
    } else {
      etaPar1 = eta1.at(iDaug1 + 1);
    }
    for (unsigned int iDaug2 = 0; iDaug2 < part2.GetPhiAtRaidius().size();
        ++iDaug2) {
      std::vector<float> phiAtRad2 = part2.GetPhiAtRaidius().at(iDaug2);
      float etaPar2;
      if (nDaug2 == 1) {
        etaPar2 = eta2.at(0);
      } else {
        etaPar2 = eta2.at(iDaug2 + 1);
      }
      float deta = etaPar1 - etaPar2;
      const int size =
          (PhiAtRad1.size() > phiAtRad2.size()) ?
              phiAtRad2.size() : PhiAtRad1.size();
      float dphiAvg = 0;
      for (int iRad = 0; iRad < size; ++iRad) {
        float dphi = PhiAtRad1.at(iRad) - phiAtRad2.at(iRad);
        dphiAvg += dphi;
        if (dphi > piHi) {
          dphi += -piHi * 2;
        } else if (dphi < -piHi) {
          dphi += piHi * 2;
        }
        dphi = TVector2::Phi_mpi_pi(dphi);
        if (pass) {
          if (dphi * dphi + deta * deta < fDeltaPhiEtaMax) {
            pass = false;
          }
        }
        if (fWhichPairs.at(Hist)) {
          if (SEorME) {
            fHists->FillEtaPhiAtRadiiSE(Hist, 3 * iDaug1 + iDaug2, iRad, dphi,
                                        deta, relk);
          } else {
            fHists->FillEtaPhiAtRadiiME(Hist, 3 * iDaug1 + iDaug2, iRad, dphi,
                                        deta, relk);
          }
        }
      }
      //fill dPhi avg
      if (fWhichPairs.at(Hist)) {
        if (SEorME) {
          fHists->FillEtaPhiAverageSE(Hist, 3 * iDaug1 + iDaug2,
                                      dphiAvg / (float) size, deta, true);
          if (pass) {
            fHists->FillEtaPhiAverageSE(Hist, 3 * iDaug1 + iDaug2,
                                        dphiAvg / (float) size, deta, false);
          }
        } else {
          fHists->FillEtaPhiAverageME(Hist, 3 * iDaug1 + iDaug2,
                                      dphiAvg / (float) size, deta, true);
          if (pass) {
            fHists->FillEtaPhiAverageME(Hist, 3 * iDaug1 + iDaug2,
                                        dphiAvg / (float) size, deta, false);
          }
        }
      }
    }
  }
  return pass;
}
