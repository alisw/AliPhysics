/*
 * AliFemtoDreamPairCleaner.cxx
 *
 *  Created on: 8 Nov 2017
 *      Author: bernhardhohlweger
 */

#include <iostream>
#include "AliFemtoDreamPairCleaner.h"
ClassImp(AliFemtoDreamPairCleaner)
AliFemtoDreamPairCleaner::AliFemtoDreamPairCleaner()
    : fMinimalBooking(false),
      fCounter(0),
      fParticles(),
      fHists(0) {
}

AliFemtoDreamPairCleaner::AliFemtoDreamPairCleaner(
    const AliFemtoDreamPairCleaner& cleaner)
    : fMinimalBooking(cleaner.fMinimalBooking),
      fCounter(0),
      fParticles(),
      fHists(cleaner.fHists) {
}

AliFemtoDreamPairCleaner::AliFemtoDreamPairCleaner(int nTrackDecayChecks,
                                                   int nDecayDecayChecks,
                                                   bool MinimalBooking)
    : fMinimalBooking(MinimalBooking),
      fCounter(0),
      fParticles(),
      fHists(nullptr) {
  if (!fMinimalBooking) {
    fHists = new AliFemtoDreamPairCleanerHists(nTrackDecayChecks,
                                               nDecayDecayChecks);
  }
}

AliFemtoDreamPairCleaner::~AliFemtoDreamPairCleaner() {
  if (fHists) {
    delete fHists;
  }
}

AliFemtoDreamPairCleaner& AliFemtoDreamPairCleaner::operator=(
    const AliFemtoDreamPairCleaner& cleaner) {
  if (this == &cleaner) {
    return *this;
  }
  this->fMinimalBooking = cleaner.fMinimalBooking;
  this->fCounter = cleaner.fCounter;
  this->fParticles = cleaner.fParticles;
  this->fHists = cleaner.fHists;
  return *this;
}

void AliFemtoDreamPairCleaner::CleanTrackAndDecay(
    std::vector<AliFemtoDreamBasePart> *Tracks,
    std::vector<AliFemtoDreamBasePart> *Decay, int histnumber) {
  int counter = 0;
  for (auto itTrack = Tracks->begin(); itTrack != Tracks->end(); ++itTrack) {
    //std::cout  << "New Track" << std::endl;
    for (auto itDecay = Decay->begin(); itDecay != Decay->end(); ++itDecay) {
      if (itDecay->UseParticle()) {
        //std::cout  << "New v0" << std::endl;
        std::vector<int> IDTrack = itTrack->GetIDTracks();
        std::vector<int> IDDaug = itDecay->GetIDTracks();
        for (auto itIDs = IDDaug.begin(); itIDs != IDDaug.end(); ++itIDs) {
          //std::cout <<"ID of Track: "<<IDTrack.at(0)<<" IDs of Daughter: "
          //              <<*itIDs<<'\n';
          if (*itIDs == IDTrack.at(0)) {
            itDecay->SetUse(false);
            counter++;
          }
        }
      } else {
        continue;
      }
    }
  }
  if (!fMinimalBooking)
    fHists->FillDaughtersSharedTrack(histnumber, counter);
}
void AliFemtoDreamPairCleaner::CleanDecayAndDecay(
    std::vector<AliFemtoDreamBasePart> *Decay1,
    std::vector<AliFemtoDreamBasePart> *Decay2, int histnumber) {
  int counter = 0;
  for (auto itDecay1 = Decay1->begin(); itDecay1 != Decay1->end(); ++itDecay1) {
    if (itDecay1->UseParticle()) {
      for (auto itDecay2 = Decay2->begin(); itDecay2 != Decay2->end();
          ++itDecay2) {
        if (itDecay2->UseParticle()) {
          std::vector<int> IDDaug1 = itDecay1->GetIDTracks();
          std::vector<int> IDDaug2 = itDecay2->GetIDTracks();
          for (auto itID1s = IDDaug1.begin(); itID1s != IDDaug1.end();
              ++itID1s) {
            for (auto itID2s = IDDaug2.begin(); itID2s != IDDaug2.end();
                ++itID2s) {
              if (*itID1s == *itID2s) {
                if (itDecay1->GetCPA() < itDecay2->GetCPA()) {
                  itDecay1->SetUse(false);
                  counter++;
                } else {
                  itDecay2->SetUse(false);
                  counter++;
                }
              }
            }
          }
        } else {
          continue;
        }
      }
    } else {
      continue;
    }
  }
  if (!fMinimalBooking)
    fHists->FillDaughtersSharedDaughter(histnumber, counter);
}

void AliFemtoDreamPairCleaner::CleanDecay(
    std::vector<AliFemtoDreamBasePart> *Decay, int histnumber) {
  int counter = 0;
  for (std::vector<AliFemtoDreamBasePart>::iterator itDecay1 = Decay->begin();
      itDecay1 != Decay->end(); ++itDecay1) {
    if (itDecay1->UseParticle()) {
      //std::cout  << "New Particle 1" << std::endl;
      for (auto itDecay2 = itDecay1 + 1; itDecay2 != Decay->end(); ++itDecay2) {
        if (itDecay2->UseParticle()) {
          //std::cout  << "New Particle 2" << std::endl;
          std::vector<int> IDDaug1 = itDecay1->GetIDTracks();
          std::vector<int> IDDaug2 = itDecay2->GetIDTracks();
          for (auto itID1s = IDDaug1.begin(); itID1s != IDDaug1.end();
              ++itID1s) {
            for (auto itID2s = IDDaug2.begin(); itID2s != IDDaug2.end();
                ++itID2s) {
              //std::cout <<"ID of Daug v01: "<<*itID1s<<" ID of Daug v02: "
              //                  <<*itID2s<<'\n';
              if (*itID1s == *itID2s) {
                if (itDecay1->GetCPA() < itDecay2->GetCPA()) {
                  itDecay1->SetUse(false);
                  counter++;
                } else {
                  itDecay2->SetUse(false);
                  counter++;
                }
              }
            }
          }
        } else {
          continue;
        }
      }
    } else {
      continue;
    }
  }
  if (!fMinimalBooking)
    fHists->FillDaughtersSharedDaughter(histnumber, counter);
}

void AliFemtoDreamPairCleaner::StoreParticle(
    std::vector<AliFemtoDreamBasePart> Particles) {
  std::vector<AliFemtoDreamBasePart> tmpParticles;
  for (auto itPart : Particles) {
    if (itPart.UseParticle()) {
      tmpParticles.push_back(itPart);
      fCounter++;
    }
  }
  fParticles.push_back(tmpParticles);
}
void AliFemtoDreamPairCleaner::ResetArray() {
  fCounter = 0;
  fParticles.clear();
}

float AliFemtoDreamPairCleaner::RelativePairMomentum(TVector3 Part1Momentum,
                                                     int PDGPart1,
                                                     TVector3 Part2Momentum,
                                                     int PDGPart2) {
  if (PDGPart1 == 0 || PDGPart2 == 0) {
    printf("Invalid PDG Code");
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
