/*
 * AliFemtoPPbpbLamZVtxMultContainer.cxx
 *
 *  Created on: Aug 30, 2017
 *      Author: gu74req
 */
//#include "AliLog.h"
#include <iostream>
#include "AliFemtoDreamZVtxMultContainer.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TVector2.h"

ClassImp(AliFemtoDreamPartContainer)
static const float piHi = TMath::Pi();
AliFemtoDreamZVtxMultContainer::AliFemtoDreamZVtxMultContainer()
    : fPartContainer(0),
      fPDGParticleSpecies(0),
      fWhichPairs(),
      fRejPairs(),
      fDeltaEtaMax(0.f),
      fDeltaPhiMax(0.f),
      fDeltaPhiEtaMax(0.f),
      fDoDeltaEtaDeltaPhiCut(false) {
  TDatabasePDG::Instance()->AddParticle("deuteron", "deuteron", 1.8756134,
                                        kTRUE, 0.0, 1, "Nucleus", 1000010020);
  TDatabasePDG::Instance()->AddAntiParticle("anti-deuteron", -1000010020);
}

AliFemtoDreamZVtxMultContainer::AliFemtoDreamZVtxMultContainer(
    AliFemtoDreamCollConfig *conf)
    : fPartContainer(conf->GetNParticles(),
                     AliFemtoDreamPartContainer(conf->GetMixingDepth())),
      fPDGParticleSpecies(conf->GetPDGCodes()),
      fWhichPairs(conf->GetWhichPairs()),
      fRejPairs(conf->GetClosePairRej()),
      fDeltaEtaMax(conf->GetDeltaEtaMax()),
      fDeltaPhiMax(conf->GetDeltaPhiMax()),
      fDeltaPhiEtaMax(
          fDeltaPhiMax * fDeltaPhiMax + fDeltaEtaMax * fDeltaEtaMax),
      fDoDeltaEtaDeltaPhiCut(conf->GetDoDeltaEtaDeltaPhiCut()) {
  TDatabasePDG::Instance()->AddParticle("deuteron", "deuteron", 1.8756134,
                                        kTRUE, 0.0, 1, "Nucleus", 1000010020);
  TDatabasePDG::Instance()->AddAntiParticle("anti-deuteron", -1000010020);
}

AliFemtoDreamZVtxMultContainer::~AliFemtoDreamZVtxMultContainer() {
  // TODO Auto-generated destructor stub
}

void AliFemtoDreamZVtxMultContainer::SetEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles) {
  //This method sets the particles of an event only in the case, that
  //more than one particle was identified, to avoid empty events.
  //  if (Particles->size()!=fParticleSpecies){
  //    TString errMessage = Form("Number of Input Particlese (%d) doese not"
  //        "correspond to the Number of particles Set (%d)",Particles->size(),
  //        fParticleSpecies);
  //    AliFatal(errMessage.Data());
  //  } else {
  std::vector<std::vector<AliFemtoDreamBasePart>>::iterator itInput = Particles
      .begin();
  std::vector<AliFemtoDreamPartContainer>::iterator itContainer = fPartContainer
      .begin();
  while (itContainer != fPartContainer.end()) {
    if (itInput->size() > 0) {
      itContainer->SetEvent(*itInput);
    }
    ++itInput;
    ++itContainer;
  }
  //  }
}
void AliFemtoDreamZVtxMultContainer::PairParticlesSE(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
    AliFemtoDreamCorrHists *ResultsHist, int iMult, float cent) {
  float RelativeK = 0;
  int HistCounter = 0;
  //First loop over all the different Species
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  for (auto itSpec1 = Particles.begin(); itSpec1 != Particles.end();
      ++itSpec1) {
    auto itPDGPar2 = fPDGParticleSpecies.begin();
    itPDGPar2 += itSpec1 - Particles.begin();
    for (auto itSpec2 = itSpec1; itSpec2 != Particles.end(); ++itSpec2) {
      ResultsHist->FillPartnersSE(HistCounter, itSpec1->size(),
                                  itSpec2->size());
      //Now loop over the actual Particles and correlate them
      unsigned int DoThisPair = fWhichPairs.at(HistCounter);
      bool fillHists = DoThisPair > 0 ? true : false;
      bool CPR = fRejPairs.at(HistCounter);
      for (auto itPart1 = itSpec1->begin(); itPart1 != itSpec1->end();
          ++itPart1) {
        AliFemtoDreamBasePart part1 = *itPart1;
        std::vector<AliFemtoDreamBasePart>::iterator itPart2;
        if (itSpec1 == itSpec2) {
          itPart2 = itPart1 + 1;
        } else {
          itPart2 = itSpec2->begin();
        }
        while (itPart2 != itSpec2->end()) {
          AliFemtoDreamBasePart part2 = *itPart2;
          // Delta eta - Delta phi* cut
          if (fDoDeltaEtaDeltaPhiCut && CPR) {
            if (!RejectClosePairs(part1, part2)) {
              ++itPart2;
              continue;
            }
          }
          RelativeK = RelativePairMomentum(itPart1->GetMomentum(), *itPDGPar1,
                                           itPart2->GetMomentum(), *itPDGPar2);

          if (fillHists && ResultsHist->GetEtaPhiPlots()) {
            DeltaEtaDeltaPhi(HistCounter, part1, part2, true, ResultsHist,
                             RelativeK);
          }
          if (fillHists && ResultsHist->GetDodPhidEtaPlots()) {
            float deta = itPart1->GetEta().at(0) - itPart2->GetEta().at(0);
            float dphi = itPart1->GetPhi().at(0) - itPart2->GetPhi().at(0);
            float mT =
                ResultsHist->GetDodPhidEtamTPlots() ?
                    RelativePairmT(itPart1->GetMomentum(), *itPDGPar1,
                                   itPart2->GetMomentum(), *itPDGPar2) :
                    0;
            if (dphi < 0) {
              ResultsHist->FilldPhidEtaSE(HistCounter, dphi + 2 * TMath::Pi(),
                                          deta, mT);
            } else {
              ResultsHist->FilldPhidEtaSE(HistCounter, dphi, deta, mT);
            }
          }

          ResultsHist->FillSameEventDist(HistCounter, RelativeK);
          if (ResultsHist->GetDoMultBinning()) {
            ResultsHist->FillSameEventMultDist(HistCounter, iMult + 1,
                                               RelativeK);
          }
          if (fillHists && ResultsHist->GetDoCentBinning()) {
            ResultsHist->FillSameEventCentDist(HistCounter, cent, RelativeK);
          }
          if (fillHists && ResultsHist->GetDokTBinning()) {
            ResultsHist->FillSameEventkTDist(
                HistCounter,
                RelativePairkT(itPart1->GetMomentum(), *itPDGPar1,
                               itPart2->GetMomentum(), *itPDGPar2),
                RelativeK, cent);
          }
          if (fillHists && ResultsHist->GetDomTBinning()) {
            ResultsHist->FillSameEventmTDist(
                HistCounter,
                RelativePairmT(itPart1->GetMomentum(), *itPDGPar1,
                               itPart2->GetMomentum(), *itPDGPar2),
                RelativeK);
          }
          if (fillHists && ResultsHist->GetDoPtQA()) {
            ResultsHist->FillPtQADist(HistCounter, RelativeK, itPart1->GetPt(),
                                      itPart2->GetPt());
          }
          if (fillHists && ResultsHist->GetDoMassQA()) {
            ResultsHist->FillMassQADist(HistCounter, RelativeK,
                                        itPart1->GetInvMass(),
                                        itPart2->GetInvMass());
	    ResultsHist->FillPairInvMassQAD(HistCounter, part1, part2);

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

void AliFemtoDreamZVtxMultContainer::PairParticlesME(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
    AliFemtoDreamCorrHists *ResultsHist, int iMult, float cent) {
  float RelativeK = 0;
  int HistCounter = 0;
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  //First loop over all the different Species
  for (auto itSpec1 = Particles.begin(); itSpec1 != Particles.end();
      ++itSpec1) {
    //We dont want to correlate the particles twice. Mixed Event Dist. of
    //Particle1 + Particle2 == Particle2 + Particle 1
    int SkipPart = itSpec1 - Particles.begin();
    auto itPDGPar2 = fPDGParticleSpecies.begin() + SkipPart;
    for (auto itSpec2 = fPartContainer.begin() + SkipPart;
        itSpec2 != fPartContainer.end(); ++itSpec2) {
      if (itSpec1->size() > 0) {
        ResultsHist->FillEffectiveMixingDepth(HistCounter,
                                              (int) itSpec2->GetMixingDepth());
      }
      unsigned int DoThisPair = fWhichPairs.at(HistCounter);
      bool fillHists = DoThisPair > 0 ? true : false;
      bool CPR = fRejPairs.at(HistCounter);
      for (int iDepth = 0; iDepth < (int) itSpec2->GetMixingDepth(); ++iDepth) {
        std::vector<AliFemtoDreamBasePart> ParticlesOfEvent = itSpec2->GetEvent(
            iDepth);
        ResultsHist->FillPartnersME(HistCounter, itSpec1->size(),
                                    ParticlesOfEvent.size());
        for (auto itPart1 = itSpec1->begin(); itPart1 != itSpec1->end();
            ++itPart1) {
          AliFemtoDreamBasePart part1 = *itPart1;
          for (auto itPart2 = ParticlesOfEvent.begin();
              itPart2 != ParticlesOfEvent.end(); ++itPart2) {
            AliFemtoDreamBasePart part2 = *itPart2;
            // Delta eta - Delta phi* cut
            if (fDoDeltaEtaDeltaPhiCut && CPR) {
              if (!RejectClosePairs(part1, part2)) {
                continue;
              }
            }
            RelativeK = RelativePairMomentum(itPart1->GetMomentum(), *itPDGPar1,
                                             itPart2->GetMomentum(),
                                             *itPDGPar2);
            if (fillHists && ResultsHist->GetEtaPhiPlots()) {
              DeltaEtaDeltaPhi(HistCounter, part1, part2, false, ResultsHist,
                               RelativeK);
            }
            if (fillHists && ResultsHist->GetDodPhidEtaPlots()) {
              float deta = itPart1->GetEta().at(0) - itPart2->GetEta().at(0);
              float dphi = itPart1->GetPhi().at(0) - itPart2->GetPhi().at(0);
              float mT =
                  ResultsHist->GetDodPhidEtamTPlots() ?
                      RelativePairmT(itPart1->GetMomentum(), *itPDGPar1,
                                     itPart2->GetMomentum(), *itPDGPar2) :
                      0;
              if (dphi < 0) {
                ResultsHist->FilldPhidEtaME(HistCounter, dphi + 2 * TMath::Pi(),
                                            deta, mT);
              } else {
                ResultsHist->FilldPhidEtaME(HistCounter, dphi, deta, mT);
              }
            }

            ResultsHist->FillMixedEventDist(HistCounter, RelativeK);
            if (ResultsHist->GetDoMultBinning()) {
              ResultsHist->FillMixedEventMultDist(HistCounter, iMult + 1,
                                                  RelativeK);
            }
            if (fillHists && ResultsHist->GetDoCentBinning()) {
              ResultsHist->FillMixedEventCentDist(HistCounter, cent, RelativeK);
            }
            if (fillHists && ResultsHist->GetDokTBinning()) {
              ResultsHist->FillMixedEventkTDist(
                  HistCounter,
                  RelativePairkT(itPart1->GetMomentum(), *itPDGPar1,
                                 itPart2->GetMomentum(), *itPDGPar2),
                  RelativeK, cent);
            }
            if (fillHists && ResultsHist->GetDomTBinning()) {
              ResultsHist->FillMixedEventmTDist(
                  HistCounter,
                  RelativePairmT(itPart1->GetMomentum(), *itPDGPar1,
                                 itPart2->GetMomentum(), *itPDGPar2),
                  RelativeK);
            }
            if (fillHists && ResultsHist->GetObtainMomentumResolution()) {
              //It is sufficient to do this in Mixed events, which allows
              //to increase the statistics. The Resolution of the tracks and therefore
              //of the pairs does not change event by event.
              //Now we only want to use the momentum of particles we are after, hence
              //we check the PDG Code!
              if ((*itPDGPar1 == TMath::Abs(itPart1->GetMCPDGCode()))
                  && ((*itPDGPar2 == TMath::Abs(itPart2->GetMCPDGCode())))) {
                float RelKTrue = RelativePairMomentum(itPart1->GetMCMomentum(),
                                                      *itPDGPar1,
                                                      itPart2->GetMCMomentum(),
                                                      *itPDGPar2);
                ResultsHist->FillMomentumResolution(HistCounter, RelKTrue,
                                                    RelativeK);
              }
            }
          }
        }
      }
      ++HistCounter;
      ++itPDGPar2;
    }
    ++itPDGPar1;
  }
}
float AliFemtoDreamZVtxMultContainer::RelativePairMomentum(
    TVector3 Part1Momentum, int PDGPart1, TVector3 Part2Momentum,
    int PDGPart2) {
  if (PDGPart1 == 0 || PDGPart2 == 0) {
    AliError("Invalid PDG Code");
  }
  float results = 0.;
  TLorentzVector SPtrack, TPProng, trackSum, SPtrackCMS, TPProngCMS;
  //Even if the Daughter tracks were switched up during PID doesn't play a role here cause we are
  //only looking at the mother mass
  SPtrack.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(), Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  TPProng.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(), Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
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
float AliFemtoDreamZVtxMultContainer::RelativePairkT(TVector3 Part1Momentum,
                                                     int PDGPart1,
                                                     TVector3 Part2Momentum,
                                                     int PDGPart2) {
  if (PDGPart1 == 0 || PDGPart2 == 0) {
    AliError("Invalid PDG Code");
  }
  float results = 0.;
  TLorentzVector SPtrack, TPProng, trackSum;
  //Even if the Daughter tracks were switched up during PID doesn't play a role here cause we are
  //only looking at the mother mass
  SPtrack.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(), Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  TPProng.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(), Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
  trackSum = SPtrack + TPProng;

  results = 0.5 * trackSum.Pt();
  return results;
}
float AliFemtoDreamZVtxMultContainer::RelativePairmT(TVector3 Part1Momentum,
                                                     int PDGPart1,
                                                     TVector3 Part2Momentum,
                                                     int PDGPart2) {
  if (PDGPart1 == 0 || PDGPart2 == 0) {
    AliError("Invalid PDG Code");
  }
  float results = 0.;
  TLorentzVector SPtrack, TPProng, trackSum;
  //Even if the Daughter tracks were switched up during PID doesn't play a role here cause we are
  //only looking at the mother mass
  SPtrack.SetXYZM(Part1Momentum.X(), Part1Momentum.Y(), Part1Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass());
  TPProng.SetXYZM(Part2Momentum.X(), Part2Momentum.Y(), Part2Momentum.Z(),
                  TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
  trackSum = SPtrack + TPProng;
  float pairKT = 0.5 * trackSum.Pt();
  float averageMass = 0.5
      * (TDatabasePDG::Instance()->GetParticle(PDGPart1)->Mass()
          + TDatabasePDG::Instance()->GetParticle(PDGPart2)->Mass());
  results = TMath::Sqrt(pow(pairKT, 2.) + pow(averageMass, 2.));
  return results;
}

void AliFemtoDreamZVtxMultContainer::DeltaEtaDeltaPhi(
    int Hist, AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2,
    bool SEorME, AliFemtoDreamCorrHists *ResultsHist, float relk) {
  //used to check for track splitting/merging
  //this function only produces meaningful results for track with x Daughter
  //looking at this quantity makes only sense anyways for Track - Track not
  //for v0 - v0 ...

  unsigned int DoThisPair = fWhichPairs.at(Hist);
  unsigned int nDaug1 = (unsigned int) DoThisPair / 10;
  if (nDaug1 > 9) {
    AliWarning("you are doing something wrong \n");
  }
  unsigned int nDaug2 = (unsigned int) DoThisPair % 10;
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
        if (SEorME) {
          ResultsHist->FillEtaPhiAtRadiiSE(Hist, 3 * iDaug1 + iDaug2, iRad,
                                           dphi, deta, relk);
        } else {
          ResultsHist->FillEtaPhiAtRadiiME(Hist, 3 * iDaug1 + iDaug2, iRad,
                                           dphi, deta, relk);
        }
      }
      //fill dPhi avg
      if (SEorME) {
        ResultsHist->FillEtaPhiAverageSE(Hist, 3 * iDaug1 + iDaug2,
                                         dphiAvg / (float) size, deta, false);
      } else {
        ResultsHist->FillEtaPhiAverageME(Hist, 3 * iDaug1 + iDaug2,
                                         dphiAvg / (float) size, deta, false);
      }
    }
  }
  return;
}

float AliFemtoDreamZVtxMultContainer::ComputeDeltaEta(
    AliFemtoDreamBasePart &part1, AliFemtoDreamBasePart &part2) {
  float eta1 = part1.GetEta().at(0);
  float eta2 = part2.GetEta().at(0);
  return std::abs(eta1 - eta2);
}

float AliFemtoDreamZVtxMultContainer::ComputeDeltaPhi(
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

bool AliFemtoDreamZVtxMultContainer::RejectClosePairs(
    AliFemtoDreamBasePart& part1, AliFemtoDreamBasePart& part2) {
  bool outBool = true;
  //Method calculates the average separation between two tracks
  //at different radii within the TPC and rejects pairs which a
  //too low separation
  unsigned int nDaug1 = part1.GetPhiAtRaidius().size();
  unsigned int nDaug2 = part2.GetPhiAtRaidius().size();
  // if nDaug == 1 => Single Track, else decay
  std::vector<float> eta1 = part1.GetEta();
  std::vector<float> eta2 = part2.GetEta();
  for (unsigned int iDaug1 = 0; iDaug1 < nDaug1 && outBool; ++iDaug1) {
    std::vector<float> PhiAtRad1 = part1.GetPhiAtRaidius().at(iDaug1);
    float etaPar1;
    if (nDaug1 == 1) {
      etaPar1 = eta1.at(0);
    } else {
      etaPar1 = eta1.at(iDaug1 + 1);
    }
    for (unsigned int iDaug2 = 0;
        iDaug2 < part2.GetPhiAtRaidius().size() && outBool; ++iDaug2) {
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
      for (int iRad = 0; iRad < size; ++iRad) {
        float dphi = PhiAtRad1.at(iRad) - phiAtRad2.at(iRad);
        if (dphi > piHi) {
          dphi += -piHi * 2;
        } else if (dphi < -piHi) {
          dphi += piHi * 2;
        }
        dphi = TVector2::Phi_mpi_pi(dphi);
        if (dphi * dphi + deta * deta < fDeltaPhiEtaMax) {
          outBool = false;
          break;
        }
      }
    }
  }
  return outBool;
}
