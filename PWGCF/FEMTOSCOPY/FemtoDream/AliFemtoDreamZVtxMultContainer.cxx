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
AliFemtoDreamZVtxMultContainer::AliFemtoDreamZVtxMultContainer()
    : fPartContainer(0),
      fPDGParticleSpecies(0),
      fWhichPairs(){
}

AliFemtoDreamZVtxMultContainer::AliFemtoDreamZVtxMultContainer(
    AliFemtoDreamCollConfig *conf)
    : fPartContainer(conf->GetNParticles(),
                     AliFemtoDreamPartContainer(conf->GetMixingDepth())),
      fPDGParticleSpecies(conf->GetPDGCodes()),
      fWhichPairs(conf->GetWhichPairs()){
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
    AliFemtoDreamHigherPairMath *HigherMath, int iMult, float cent) {
  int HistCounter = 0;
  //First loop over all the different Species
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  for (auto itSpec1 = Particles.begin(); itSpec1 != Particles.end();
      ++itSpec1) {
    auto itPDGPar2 = fPDGParticleSpecies.begin();
    itPDGPar2 += itSpec1 - Particles.begin();
    for (auto itSpec2 = itSpec1; itSpec2 != Particles.end(); ++itSpec2) {
      HigherMath->FillPairCounterSE(HistCounter, itSpec1->size(),
                                    itSpec2->size());
      //Now loop over the actual Particles and correlate them
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
          TLorentzVector PartOne, PartTwo;
          PartOne.SetXYZM(
              itPart1->GetMomentum().X(), itPart1->GetMomentum().Y(),
              itPart1->GetMomentum().Z(),
              TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass());
          PartTwo.SetXYZM(
              itPart2->GetMomentum().X(), itPart2->GetMomentum().Y(),
              itPart2->GetMomentum().Z(),
              TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass());
          float RelativeK = HigherMath->RelativePairMomentum(PartOne, PartTwo);
          if (!HigherMath->PassesPairSelection(HistCounter, *itPart1, *itPart2,
                                               RelativeK, true, false)) {
            ++itPart2;
            continue;
          }
          RelativeK = HigherMath->FillSameEvent(HistCounter, iMult, cent,
                                                part1,
                                                *itPDGPar1,
                                                part2,
                                                *itPDGPar2);
          HigherMath->MassQA(HistCounter, RelativeK, *itPart1, *itPDGPar1,
                                                     *itPart2, *itPDGPar2);
          HigherMath->SEDetaDPhiPlots(HistCounter, *itPart1, *itPDGPar1,
                                      *itPart2, *itPDGPar2, RelativeK, false);
          HigherMath->SEMomentumResolution(HistCounter, &(*itPart1), *itPDGPar1,
                                           &(*itPart2), *itPDGPar2, RelativeK);
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
    AliFemtoDreamHigherPairMath *HigherMath, int iMult, float cent) {
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
        HigherMath->FillEffectiveMixingDepth(HistCounter,
                                             (int) itSpec2->GetMixingDepth());
      }
      for (int iDepth = 0; iDepth < (int) itSpec2->GetMixingDepth(); ++iDepth) {
        std::vector<AliFemtoDreamBasePart> ParticlesOfEvent = itSpec2->GetEvent(
            iDepth);
        HigherMath->FillPairCounterME(HistCounter, itSpec1->size(),
                                      ParticlesOfEvent.size());
        for (auto itPart1 = itSpec1->begin(); itPart1 != itSpec1->end();
            ++itPart1) {
          for (auto itPart2 = ParticlesOfEvent.begin();
              itPart2 != ParticlesOfEvent.end(); ++itPart2) {

            TLorentzVector PartOne, PartTwo;
            PartOne.SetXYZM(
                itPart1->GetMomentum().X(), itPart1->GetMomentum().Y(),
                itPart1->GetMomentum().Z(),
                TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass());
            PartTwo.SetXYZM(
                itPart2->GetMomentum().X(), itPart2->GetMomentum().Y(),
                itPart2->GetMomentum().Z(),
                TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass());
            float RelativeK = HigherMath->RelativePairMomentum(PartOne, PartTwo);
            if (!HigherMath->PassesPairSelection(HistCounter, *itPart1, *itPart2,
                                                 RelativeK, false, false)) {
              continue;
            }
            RelativeK = HigherMath->FillMixedEvent(
                HistCounter, iMult, cent, *itPart1, *itPDGPar1,
                *itPart2, *itPDGPar2,
                AliFemtoDreamCollConfig::kNone);

            HigherMath->MEMassQA(HistCounter, RelativeK, *itPart1, *itPDGPar1,
                                                         *itPart2, *itPDGPar2);
            HigherMath->MEDetaDPhiPlots(HistCounter, *itPart1, *itPDGPar1,
                                        *itPart2, *itPDGPar2, RelativeK, false);
            HigherMath->MEMomentumResolution(HistCounter, &(*itPart1),
                                             *itPDGPar1, &(*itPart2),
                                             *itPDGPar2, RelativeK);
          }
        }
      }
      ++HistCounter;
      ++itPDGPar2;
    }
    ++itPDGPar1;
  }
}
