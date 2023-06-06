/*
 * AliFemtoPPbpbLamZVtxMultContainer.cxx
 *
 *  Created on: Aug 30, 2017
 *      Author: gu74req
 */
//#include "AliLog.h"
#include <map>
#include <utility>

#include <iostream>
#include "AliFemtoDreamZVtxMultContainer.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TVector2.h"
#include "TTree.h"

ClassImp(AliFemtoDreamPartContainer)
AliFemtoDreamZVtxMultContainer::AliFemtoDreamZVtxMultContainer()
    : fPartContainer(0),
      fPDGParticleSpecies(0),
      fWhichPairs(),
      fSummedPtLimit1(0.0),
      fSummedPtLimit2(999.0){
}

AliFemtoDreamZVtxMultContainer::AliFemtoDreamZVtxMultContainer(
    AliFemtoDreamCollConfig *conf)
    : fPartContainer(conf->GetNParticles(),
                     AliFemtoDreamPartContainer(conf->GetMixingDepth())),
      fPDGParticleSpecies(conf->GetPDGCodes()),
      fWhichPairs(conf->GetWhichPairs()),
      fSummedPtLimit1(conf->GetSummedPtLimit1()),
      fSummedPtLimit2(conf->GetSummedPtLimit2()){
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
    AliFemtoDreamHigherPairMath *HigherMath, int iMult, float cent,
    std::map<std::pair<int, int>, TTree *> *kStarsSE) {
  int HistCounter = 0;
  //First loop over all the different Species
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  for (auto itSpec1 = Particles.begin(); itSpec1 != Particles.end();
      ++itSpec1) {
    int iParticles1 = std::distance(Particles.begin(), itSpec1);

    auto itPDGPar2 = fPDGParticleSpecies.begin();
    itPDGPar2 += itSpec1 - Particles.begin();
    for (auto itSpec2 = itSpec1; itSpec2 != Particles.end(); ++itSpec2) {
      int iParticles2 = std::distance(Particles.begin(), itSpec2);

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
        auto itStartPart2 = itPart2;
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
          if (!HigherMath->PassesMDPairSelection(*itPart1, *itPart2)) {
            ++itPart2;
            continue;
          }
          RelativeK = HigherMath->FillSameEvent(HistCounter, iMult, cent,
                                                part1,
                                                *itPDGPar1,
                                                part2,
                                                *itPDGPar2,
						fSummedPtLimit1,
						fSummedPtLimit2);

          // save the list of k* for each pair (required for Dmeson-LF analyses)
          if (kStarsSE) {
            auto combo = std::pair<int, int>({iParticles1, iParticles2});
            if (kStarsSE->find(combo) != kStarsSE->end()) { // ignore uninteresting pairs
              auto tree = kStarsSE->at(std::pair<int, int>({iParticles1, iParticles2}));

              // load event info
              int mult = part2.GetMult();
              float zvtx = part1.GetZVtx();

              // pair
              bool is_oldpcrm = part1.IsRemovedByOldPC() || part2.IsRemovedByOldPC();
              bool is_newpcrm = part1.IsRemovedByNewPC() || part2.IsRemovedByNewPC();
              bool is_crosspcrm = part1.IsRemovedByCrossPC() || part2.IsRemovedByCrossPC();

              // auto p1 = TLorentzVector(part1.GetMomentum(), TMath::Sqrt(part1.GetInvMass() * part1.GetInvMass() + part1.GetMomentum().Mag2()));
              // auto p2 = TLorentzVector(part2.GetMomentum(), TMath::Sqrt(part2.GetInvMass() * part2.GetInvMass() + part2.GetMomentum().Mag2()));
              // float inv_mass = (p1 + p2).Mag();

              // auto p1pdg = TLorentzVector(part1.GetMomentum(), TMath::Sqrt(0.139 * 0.139 + part1.GetMomentum().Mag2()));
              // auto p2pdg = TLorentzVector(part2.GetMomentum(), TMath::Sqrt(2.010 * 2.010 + part2.GetMomentum().Mag2()));
              // float inv_masspdg = (p1pdg + p2pdg).Mag();
              
              // load dmeson info
              int heavy_mult = part2.GetParticleMult();
              float heavy_invmass = part2.GetInvMass();
              float heavy_pt = part2.GetPt();
              float heavy_eta = part2.GetEta()[0];
              int heavy_origin = part2.GetParticleOrigin();
              std::vector<int>  *heavy_daus = new std::vector<int>(part2.GetIDTracks());
              float heavy_softpion_px = part2.GetSoftPionPx();
              float heavy_softpion_py = part2.GetSoftPionPy();
              float heavy_softpion_pz = part2.GetSoftPionPz();
              float heavy_bkgscore = part2.GetBkgScore();
              float heavy_promptscore = part2.GetPromptScore();
              int heavy_d0label = part2.GetDzeroLabel();

              // load light info
              int light_mult = part1.GetParticleMult();
              float light_px = part1.GetPx();
              float light_py = part1.GetPy();
              float light_pz = part1.GetPz();
              float light_eta = part1.GetEta()[0];
              float light_nsigtpc = part1.GetNSigTPC();
              float light_nsigtof = part1.GetNSigTOF();
              int light_ncls = part1.GetNCls();
              int light_ncrossed = part1.GetNCrossedRows();
              float light_dcaz = part1.GetDCAZ();
              float light_dcaxy = part1.GetDCAXY();
              int light_label = part1.GetID();
              int light_pdg = part1.GetPDGCode();
              int light_origin = part1.GetParticleOrigin();
              bool light_isprim = part1.GetIsPrim();
              int light_motherPdg = part1.GetMotherPDG();

              // event
              if (tree->FindBranch("mult")) tree->SetBranchAddress("mult", &mult);
              if (tree->FindBranch("vz")) tree->SetBranchAddress("vz", &zvtx);

              // pair
              if (tree->FindBranch("kStar")) tree->SetBranchAddress("kStar", &RelativeK);
              if (tree->FindBranch("is_oldpcrm")) tree->SetBranchAddress("is_oldpcrm", &is_oldpcrm);
              if (tree->FindBranch("is_newpcrm")) tree->SetBranchAddress("is_newpcrm", &is_newpcrm);
              if (tree->FindBranch("is_crosspcrm")) tree->SetBranchAddress("is_crosspcrm", &is_crosspcrm);
              // if (tree->FindBranch("inv_mass")) tree->SetBranchAddress("inv_mass", &inv_mass);
              // if (tree->FindBranch("inv_masspdg")) tree->SetBranchAddress("inv_masspdg", &inv_masspdg);


              // heavy particle
              if (tree->FindBranch("heavy_mult")) tree->SetBranchAddress("heavy_mult", &heavy_mult);
              if (tree->FindBranch("heavy_invmass")) tree->SetBranchAddress("heavy_invmass", &heavy_invmass);
              if (tree->FindBranch("heavy_pt")) tree->SetBranchAddress("heavy_pt", &heavy_pt);
              if (tree->FindBranch("heavy_eta")) tree->SetBranchAddress("heavy_eta", &heavy_eta);
              if (tree->FindBranch("heavy_origin")) tree->SetBranchAddress("heavy_origin", &heavy_origin);
              if (tree->FindBranch("heavy_daus")) tree->SetBranchAddress("heavy_daus", &heavy_daus);
              if (tree->FindBranch("heavy_softpion_px")) tree->SetBranchAddress("heavy_softpion_px", &heavy_softpion_px);
              if (tree->FindBranch("heavy_softpion_py")) tree->SetBranchAddress("heavy_softpion_py", &heavy_softpion_py);
              if (tree->FindBranch("heavy_softpion_pz")) tree->SetBranchAddress("heavy_softpion_pz", &heavy_softpion_pz);
              if (tree->FindBranch("heavy_bkg_score")) tree->SetBranchAddress("heavy_bkg_score", &heavy_bkgscore);
              if (tree->FindBranch("heavy_prompt_score")) tree->SetBranchAddress("heavy_prompt_score", &heavy_promptscore);
              if (tree->FindBranch("heavy_d0label")) tree->SetBranchAddress("heavy_d0label", &heavy_d0label);

              // light particle
              if (tree->FindBranch("light_mult")) tree->SetBranchAddress("light_mult", &light_mult);
              if (tree->FindBranch("light_px")) tree->SetBranchAddress("light_px", &light_px);
              if (tree->FindBranch("light_py")) tree->SetBranchAddress("light_py", &light_py);
              if (tree->FindBranch("light_pz")) tree->SetBranchAddress("light_pz", &light_pz);
              if (tree->FindBranch("light_eta")) tree->SetBranchAddress("light_eta", &light_eta);
              if (tree->FindBranch("light_nsigtpc")) tree->SetBranchAddress("light_nsigtpc", &light_nsigtpc);
              if (tree->FindBranch("light_nsigtof")) tree->SetBranchAddress("light_nsigtof", &light_nsigtof);
              if (tree->FindBranch("light_ncls")) tree->SetBranchAddress("light_ncls", &light_ncls);
              if (tree->FindBranch("light_ncrossed")) tree->SetBranchAddress("light_ncrossed", &light_ncrossed);
              if (tree->FindBranch("light_dcaz")) tree->SetBranchAddress("light_dcaz", &light_dcaz);
              if (tree->FindBranch("light_dcaxy")) tree->SetBranchAddress("light_dcaxy", &light_dcaxy);
              if (tree->FindBranch("light_label")) tree->SetBranchAddress("light_label", &light_label);
              if (tree->FindBranch("light_pdg")) tree->SetBranchAddress("light_pdg", &light_pdg);
              if (tree->FindBranch("light_origin")) tree->SetBranchAddress("light_origin", &light_origin);
              if (tree->FindBranch("light_isprim")) tree->SetBranchAddress("light_isprim", &light_isprim);
              if (tree->FindBranch("light_motherpdg")) tree->SetBranchAddress("light_motherpdg", &light_motherPdg);
              tree->Fill();
            }
          }

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
  AliFemtoDreamHigherPairMath *HigherMath,
  int iMult,
  float cent,
  std::map<std::pair<int, int>, TTree *> *kStarsME,
  bool usePart2Buffer){
    if (usePart2Buffer) PairParticlesMEPart2Buffer( Particles, HigherMath, iMult, cent, kStarsME);
    else PairParticlesMEPart1Buffer( Particles, HigherMath, iMult, cent, kStarsME);
  }

void AliFemtoDreamZVtxMultContainer::PairParticlesMEPart2Buffer(
    std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
    AliFemtoDreamHigherPairMath *HigherMath, int iMult, float cent,
    std::map<std::pair<int, int>, TTree *> *kStarsME) {
  int HistCounter = 0;
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  //First loop over all the different Species
  for (auto itSpec1 = Particles.begin(); itSpec1 != Particles.end();
      ++itSpec1) {
    int iParticles1 = std::distance(Particles.begin(), itSpec1);
    //We dont want to correlate the particles twice. Mixed Event Dist. of
    //Particle1 + Particle2 == Particle2 + Particle 1
    int SkipPart = itSpec1 - Particles.begin();
    auto itPDGPar2 = fPDGParticleSpecies.begin() + SkipPart;
    for (auto itSpec2 = fPartContainer.begin() + SkipPart;
        itSpec2 != fPartContainer.end(); ++itSpec2) {
      int iParticles2 = std::distance(fPartContainer.begin(), itSpec2);
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

            if (kStarsME) {
              auto combo = std::pair<int, int>({iParticles1, iParticles2});
              if (kStarsME->find(combo) != kStarsME->end()) { // ignore uninteresting pairs
                auto tree = kStarsME->at(std::pair<int, int>({iParticles1, iParticles2}));

                // load event info
                int mult = itPart2->GetMult();
                float zvtx = itPart2->GetZVtx();
                
                // pair
                bool is_oldpcrm = itPart1->IsRemovedByOldPC() || itPart2->IsRemovedByOldPC();
                bool is_newpcrm = itPart1->IsRemovedByNewPC() || itPart2->IsRemovedByNewPC();
                bool is_crosspcrm = itPart1->IsRemovedByCrossPC() || itPart2->IsRemovedByCrossPC();
              
                // load dmeson info
                int heavy_mult = itPart2->GetParticleMult();
                float heavy_invmass = itPart2->GetInvMass();
                float heavy_pt = itPart2->GetPt();
                float heavy_eta = itPart2->GetEta()[0];
                int heavy_origin = itPart2->GetParticleOrigin();
                std::vector<Int_t>  * heavy_daus = new std::vector<int>(itPart2->GetIDTracks());
                float heavy_softpion_px = itPart2->GetSoftPionPx();
                float heavy_softpion_py = itPart2->GetSoftPionPy();
                float heavy_softpion_pz = itPart2->GetSoftPionPz();
                float heavy_bkgscore = itPart2->GetBkgScore();
                float heavy_promptscore = itPart2->GetPromptScore();
                int heavy_d0label = itPart2->GetDzeroLabel();

                // load light flavor info
                int light_mult = itPart1->GetParticleMult();
                float light_px = itPart1->GetPx();
                float light_py = itPart1->GetPy();
                float light_pz = itPart1->GetPz();
                float light_eta = itPart1->GetEta()[0];
                float light_nsigtpc = itPart1->GetNSigTPC();
                float light_nsigtof = itPart1->GetNSigTOF();
                int light_ncls = itPart1->GetNCls();
                int light_ncrossed = itPart1->GetNCrossedRows();
                float light_dcaz = itPart1->GetDCAZ();
                float light_dcaxy = itPart1->GetDCAXY();
                int light_label = itPart1->GetID();
                int light_pdg = itPart1->GetPDGCode();
                int light_origin = itPart1->GetParticleOrigin();
                bool light_isprim = itPart1->GetIsPrim();
                int light_motherPdg = itPart1->GetMotherPDG();

                // event
                if (tree->FindBranch("mult")) tree->SetBranchAddress("mult", &mult);
                if (tree->FindBranch("vz")) tree->SetBranchAddress("vz", &zvtx);

                // pair
                if (tree->FindBranch("kStar")) tree->SetBranchAddress("kStar", &RelativeK);
                if (tree->FindBranch("is_oldpcrm")) tree->SetBranchAddress("is_oldpcrm", &is_oldpcrm);
                if (tree->FindBranch("is_newpcrm")) tree->SetBranchAddress("is_newpcrm", &is_newpcrm);
                if (tree->FindBranch("is_crosspcrm")) tree->SetBranchAddress("is_crosspcrm", &is_crosspcrm);
                // if (tree->FindBranch("inv_mass")) tree->SetBranchAddress("inv_mass", &inv_mass);
                // if (tree->FindBranch("inv_masspdg")) tree->SetBranchAddress("inv_masspdg", &inv_masspdg);

                // heavy particle
                if (tree->FindBranch("heavy_mult")) tree->SetBranchAddress("heavy_mult", &heavy_mult);
                if (tree->FindBranch("heavy_invmass")) tree->SetBranchAddress("heavy_invmass", &heavy_invmass);
                if (tree->FindBranch("heavy_pt")) tree->SetBranchAddress("heavy_pt", &heavy_pt);
                if (tree->FindBranch("heavy_eta")) tree->SetBranchAddress("heavy_eta", &heavy_eta);
                if (tree->FindBranch("heavy_origin")) tree->SetBranchAddress("heavy_origin", &heavy_origin);
                if (tree->FindBranch("heavy_daus")) tree->SetBranchAddress("heavy_daus", &heavy_daus);
                if (tree->FindBranch("heavy_softpion_px")) tree->SetBranchAddress("heavy_softpion_px", &heavy_softpion_px);
                if (tree->FindBranch("heavy_softpion_py")) tree->SetBranchAddress("heavy_softpion_py", &heavy_softpion_py);
                if (tree->FindBranch("heavy_softpion_pz")) tree->SetBranchAddress("heavy_softpion_pz", &heavy_softpion_pz);
                if (tree->FindBranch("heavy_bkg_score")) tree->SetBranchAddress("heavy_bkg_score", &heavy_bkgscore);
                if (tree->FindBranch("heavy_prompt_score")) tree->SetBranchAddress("heavy_prompt_score", &heavy_promptscore);
                if (tree->FindBranch("heavy_d0label")) tree->SetBranchAddress("heavy_d0label", &heavy_d0label);

                // light particle
                if (tree->FindBranch("light_mult")) tree->SetBranchAddress("light_mult", &light_mult);
                if (tree->FindBranch("light_px")) tree->SetBranchAddress("light_px", &light_px);
                if (tree->FindBranch("light_py")) tree->SetBranchAddress("light_py", &light_py);
                if (tree->FindBranch("light_pz")) tree->SetBranchAddress("light_pz", &light_pz);
                if (tree->FindBranch("light_eta")) tree->SetBranchAddress("light_eta", &light_eta);
                if (tree->FindBranch("light_nsigtpc")) tree->SetBranchAddress("light_nsigtpc", &light_nsigtpc);
                if (tree->FindBranch("light_nsigtof")) tree->SetBranchAddress("light_nsigtof", &light_nsigtof);
                if (tree->FindBranch("light_ncls")) tree->SetBranchAddress("light_ncls", &light_ncls);
                if (tree->FindBranch("light_ncrossed")) tree->SetBranchAddress("light_ncrossed", &light_ncrossed);
                if (tree->FindBranch("light_dcaz")) tree->SetBranchAddress("light_dcaz", &light_dcaz);
                if (tree->FindBranch("light_dcaxy")) tree->SetBranchAddress("light_dcaxy", &light_dcaxy);
                if (tree->FindBranch("light_label")) tree->SetBranchAddress("light_label", &light_label);
                if (tree->FindBranch("light_pdg")) tree->SetBranchAddress("light_pdg", &light_pdg);
                if (tree->FindBranch("light_origin")) tree->SetBranchAddress("light_origin", &light_origin);
                if (tree->FindBranch("light_isprim")) tree->SetBranchAddress("light_isprim", &light_isprim);
                if (tree->FindBranch("light_motherpdg")) tree->SetBranchAddress("light_motherpdg", &light_motherPdg);
                tree->Fill();
              }
            }
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

void AliFemtoDreamZVtxMultContainer::PairParticlesMEPart1Buffer(
  std::vector<std::vector<AliFemtoDreamBasePart>> &Particles,
  AliFemtoDreamHigherPairMath *HigherMath,
  int iMult,
  float cent,
  std::map<std::pair<int, int>, TTree *> *kStarsME){
  
  int HistCounter = 0;
  auto itPDGPar1 = fPDGParticleSpecies.begin();
  //First loop over all the different Species
  for (auto itSpec_to1 = fPartContainer.begin(); itSpec_to1 != fPartContainer.end(); ++itSpec_to1) {
    //We dont want to correlate the particles twice. Mixed Event Dist. of
    //Particle1 + Particle2 == Particle2 + Particle 1
    
    int SkipPart = itSpec_to1 - fPartContainer.begin();
    int iParticles1 = std::distance(fPartContainer.begin(), itSpec_to1);
    auto itPDGPar2 = fPDGParticleSpecies.begin() + SkipPart;

    for (auto itSpec_to2 = Particles.begin() + SkipPart; itSpec_to2 != Particles.end(); ++itSpec_to2) {
      int iParticles2 = std::distance(Particles.begin(), itSpec_to2);

      if (itSpec_to2->size() > 0) {
        HigherMath->FillEffectiveMixingDepth(HistCounter, (int) itSpec_to1->GetMixingDepth());
      }

      for (int iDepth = 0; iDepth < (int) itSpec_to1->GetMixingDepth(); ++iDepth) {
        std::vector<AliFemtoDreamBasePart> ParticlesOfEvent = itSpec_to1->GetEvent(iDepth);
        HigherMath->FillPairCounterME(HistCounter, itSpec_to2->size(), ParticlesOfEvent.size());
        
        for (auto itPart1 = ParticlesOfEvent.begin(); itPart1 != ParticlesOfEvent.end(); ++itPart1) {
          for (auto itPart2 = itSpec_to2->begin(); itPart2 != itSpec_to2->end(); ++itPart2) {
            
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
            
            // continue;
            RelativeK = HigherMath->FillMixedEvent(
                HistCounter, iMult, cent, *itPart1, *itPDGPar1,
                *itPart2, *itPDGPar2,
                AliFemtoDreamCollConfig::kNone);

            if (kStarsME) {
              auto combo = std::pair<int, int>({iParticles1, iParticles2});
              if (kStarsME->find(combo) != kStarsME->end()) { // ignore uninteresting pairs
                auto tree = kStarsME->at(std::pair<int, int>({iParticles1, iParticles2}));

                // load event info
                int mult = itPart2->GetMult();
                float zvtx = itPart2->GetZVtx();
                
                // pair
                bool is_oldpcrm = itPart1->IsRemovedByOldPC() || itPart2->IsRemovedByOldPC();
                bool is_newpcrm = itPart1->IsRemovedByNewPC() || itPart2->IsRemovedByNewPC();
                bool is_crosspcrm = itPart1->IsRemovedByCrossPC() || itPart2->IsRemovedByCrossPC();
              
                // load dmeson info
                int heavy_mult = itPart2->GetParticleMult();
                float heavy_invmass = itPart2->GetInvMass();
                float heavy_pt = itPart2->GetPt();
                float heavy_eta = itPart2->GetEta()[0];
                int heavy_origin = itPart2->GetParticleOrigin();
                std::vector<Int_t>  * heavy_daus = new std::vector<int>(itPart2->GetIDTracks());
                float heavy_softpion_px = itPart2->GetSoftPionPx();
                float heavy_softpion_py = itPart2->GetSoftPionPy();
                float heavy_softpion_pz = itPart2->GetSoftPionPz();
                float heavy_bkgscore = itPart2->GetBkgScore();
                float heavy_promptscore = itPart2->GetPromptScore();
                int heavy_d0label = itPart2->GetDzeroLabel();

                // load light flavor info
                int light_mult = itPart1->GetParticleMult();
                float light_px = itPart1->GetPx();
                float light_py = itPart1->GetPy();
                float light_pz = itPart1->GetPz();
                float light_eta = itPart1->GetEta()[0];
                float light_nsigtpc = itPart1->GetNSigTPC();
                float light_nsigtof = itPart1->GetNSigTOF();
                int light_ncls = itPart1->GetNCls();
                int light_ncrossed = itPart1->GetNCrossedRows();
                float light_dcaz = itPart1->GetDCAZ();
                float light_dcaxy = itPart1->GetDCAXY();
                int light_label = itPart1->GetID();
                int light_pdg = itPart1->GetPDGCode();
                int light_origin = itPart1->GetParticleOrigin();
                bool light_isprim = itPart1->GetIsPrim();
                int light_motherPdg = itPart1->GetMotherPDG();

                // event
                if (tree->FindBranch("mult")) tree->SetBranchAddress("mult", &mult);
                if (tree->FindBranch("vz")) tree->SetBranchAddress("vz", &zvtx);

                // pair
                if (tree->FindBranch("kStar")) tree->SetBranchAddress("kStar", &RelativeK);
                if (tree->FindBranch("is_oldpcrm")) tree->SetBranchAddress("is_oldpcrm", &is_oldpcrm);
                if (tree->FindBranch("is_newpcrm")) tree->SetBranchAddress("is_newpcrm", &is_newpcrm);
                if (tree->FindBranch("is_crosspcrm")) tree->SetBranchAddress("is_crosspcrm", &is_crosspcrm);
                // if (tree->FindBranch("inv_mass")) tree->SetBranchAddress("inv_mass", &inv_mass);
                // if (tree->FindBranch("inv_masspdg")) tree->SetBranchAddress("inv_masspdg", &inv_masspdg);

                // heavy particle
                if (tree->FindBranch("heavy_mult")) tree->SetBranchAddress("heavy_mult", &heavy_mult);
                if (tree->FindBranch("heavy_invmass")) tree->SetBranchAddress("heavy_invmass", &heavy_invmass);
                if (tree->FindBranch("heavy_pt")) tree->SetBranchAddress("heavy_pt", &heavy_pt);
                if (tree->FindBranch("heavy_eta")) tree->SetBranchAddress("heavy_eta", &heavy_eta);
                if (tree->FindBranch("heavy_origin")) tree->SetBranchAddress("heavy_origin", &heavy_origin);
                if (tree->FindBranch("heavy_daus")) tree->SetBranchAddress("heavy_daus", &heavy_daus);
                if (tree->FindBranch("heavy_softpion_px")) tree->SetBranchAddress("heavy_softpion_px", &heavy_softpion_px);
                if (tree->FindBranch("heavy_softpion_py")) tree->SetBranchAddress("heavy_softpion_py", &heavy_softpion_py);
                if (tree->FindBranch("heavy_softpion_pz")) tree->SetBranchAddress("heavy_softpion_pz", &heavy_softpion_pz);
                if (tree->FindBranch("heavy_bkg_score")) tree->SetBranchAddress("heavy_bkg_score", &heavy_bkgscore);
                if (tree->FindBranch("heavy_prompt_score")) tree->SetBranchAddress("heavy_prompt_score", &heavy_promptscore);
                if (tree->FindBranch("heavy_d0label")) tree->SetBranchAddress("heavy_d0label", &heavy_d0label);

                // light particle
                if (tree->FindBranch("light_mult")) tree->SetBranchAddress("light_mult", &light_mult);
                if (tree->FindBranch("light_px")) tree->SetBranchAddress("light_px", &light_px);
                if (tree->FindBranch("light_py")) tree->SetBranchAddress("light_py", &light_py);
                if (tree->FindBranch("light_pz")) tree->SetBranchAddress("light_pz", &light_pz);
                if (tree->FindBranch("light_eta")) tree->SetBranchAddress("light_eta", &light_eta);
                if (tree->FindBranch("light_nsigtpc")) tree->SetBranchAddress("light_nsigtpc", &light_nsigtpc);
                if (tree->FindBranch("light_nsigtof")) tree->SetBranchAddress("light_nsigtof", &light_nsigtof);
                if (tree->FindBranch("light_ncls")) tree->SetBranchAddress("light_ncls", &light_ncls);
                if (tree->FindBranch("light_ncrossed")) tree->SetBranchAddress("light_ncrossed", &light_ncrossed);
                if (tree->FindBranch("light_dcaz")) tree->SetBranchAddress("light_dcaz", &light_dcaz);
                if (tree->FindBranch("light_dcaxy")) tree->SetBranchAddress("light_dcaxy", &light_dcaxy);
                if (tree->FindBranch("light_label")) tree->SetBranchAddress("light_label", &light_label);
                if (tree->FindBranch("light_pdg")) tree->SetBranchAddress("light_pdg", &light_pdg);
                if (tree->FindBranch("light_origin")) tree->SetBranchAddress("light_origin", &light_origin);
                if (tree->FindBranch("light_isprim")) tree->SetBranchAddress("light_isprim", &light_isprim);
                if (tree->FindBranch("light_motherpdg")) tree->SetBranchAddress("light_motherpdg", &light_motherPdg);
                tree->Fill();
              }
            }
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