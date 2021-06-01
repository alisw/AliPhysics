/************************************************************************************
 * Copyright (C) 2018, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <iostream>
#include <sstream>
#include <string>

#include <TClonesArray.h>
#include <THistManager.h>
#include <TList.h>
#include <TString.h>

#include "AliAnalysisManager.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskEmcalPatchMultCorr.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalPatchMultCorr)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalPatchMultCorr::AliAnalysisTaskEmcalPatchMultCorr():
  AliAnalysisTaskEmcal(),
  fHistos(nullptr)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetUseAliAnaUtils(true);
  SetMakeGeneralHistograms(true);
}

AliAnalysisTaskEmcalPatchMultCorr::AliAnalysisTaskEmcalPatchMultCorr(const char *name):
  AliAnalysisTaskEmcal(name, true),
  fHistos(nullptr)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetUseAliAnaUtils(true);
  SetMakeGeneralHistograms(true);
}

AliAnalysisTaskEmcalPatchMultCorr::~AliAnalysisTaskEmcalPatchMultCorr(){

}

void AliAnalysisTaskEmcalPatchMultCorr::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fHistos = new THistManager(Form("histos_%s", GetName()));
  fHistos->CreateTH2("hMultCorrGAJE", "Multiplicity correlation firing gamma/jet patches; Mult JE; Mult GA", 100., 0., 100., 100., 0., 100.);
  fHistos->CreateTH1("hCountIsolatedGA", "Number of GA patches not contained by a JE patch; Isolated GA patches; Number of events", 100, 0., 100.);
  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);
}

bool AliAnalysisTaskEmcalPatchMultCorr::Run(){
  // Only study EJ1 events
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE)) return false; 
  if(!TString(fInputEvent->GetFiredTriggerClasses()).Contains("EJ1")) return false;
  auto ejepatches = SelectPatches(fTriggerPatchInfo, 255, true),
       egapatches = SelectPatches(fTriggerPatchInfo, 115, false),
       isolatedGA = GetIsolatedGammaPatches(egapatches, ejepatches);
  fHistos->FillTH2("hMultCorrGAJE", ejepatches.size(), egapatches.size());
  fHistos->FillTH1("hCountIsolatedGA", isolatedGA.size());
  return true;
}

const std::vector<const AliEMCALTriggerPatchInfo *> AliAnalysisTaskEmcalPatchMultCorr::SelectPatches(const TClonesArray *patchcontainer, double threshold, bool isEJE) const{
  std::vector<const AliEMCALTriggerPatchInfo *> selpatches;
  for(auto p : *patchcontainer){
    auto testpatch = static_cast<AliEMCALTriggerPatchInfo *>(p);
    if((isEJE && !testpatch->IsJetHighRecalc()) || (!isEJE && !testpatch->IsGammaHighRecalc())) continue;
    if(testpatch->GetADCAmp() < threshold) continue;
    selpatches.emplace_back(testpatch);
  } 
  return selpatches;  
}

const std::vector<const AliEMCALTriggerPatchInfo *> AliAnalysisTaskEmcalPatchMultCorr::GetIsolatedGammaPatches(const std::vector<const AliEMCALTriggerPatchInfo *> &gapaches, const std::vector<const AliEMCALTriggerPatchInfo *> &jepatches) const{
  std::vector<const AliEMCALTriggerPatchInfo *> isolated;
  for(auto ga : gapaches) {
    bool hasOverlap = false;
    for(auto je : jepatches){
      if(CheckGAJEOverlap(ga, je)) {
        hasOverlap = true;
        break;
      }
    }
    if(!hasOverlap) isolated.emplace_back(ga);
  }
  return isolated;
}

bool AliAnalysisTaskEmcalPatchMultCorr::CheckGAJEOverlap(const AliEMCALTriggerPatchInfo *gapatch, const AliEMCALTriggerPatchInfo *jepatch) const{
  int garowmin = gapatch->GetRowStart(), garowmax = gapatch->GetRowStart() + gapatch->GetPatchSize() - 1,
      gacolmin = gapatch->GetColStart(), gacolmax = gapatch->GetColStart() + gapatch->GetPatchSize() - 1,
      jerowmin = jepatch->GetRowStart(), jerowmax = jepatch->GetRowStart() + jepatch->GetPatchSize() - 1,
      jecolmin = jepatch->GetColStart(), jecolmax = jepatch->GetColStart() + jepatch->GetPatchSize() - 1;
  // Check whether any of the 4 edges of the GA patch is contianed within the square spanned by the jet patch
  return InSquare(gacolmin, garowmin, jecolmin, jerowmin, jecolmax, jerowmax) || 
         InSquare(gacolmax, garowmin, jecolmin, jerowmin, jecolmax, jerowmax) || 
         InSquare(gacolmin, garowmax, jecolmin, jerowmin, jecolmax, jerowmax) || 
         InSquare(gacolmax, garowmax, jecolmin, jerowmin, jecolmax, jerowmax);
}

bool AliAnalysisTaskEmcalPatchMultCorr::InSquare(int px, int py, int sxlow, int sylow, int sxhigh, int syhigh) const {
  // all edge coordinates are part of the square
  return px >= sxlow && px <= sxhigh && py >= sylow && py <= syhigh;
}

AliAnalysisTaskEmcalPatchMultCorr *AliAnalysisTaskEmcalPatchMultCorr::AddTaskEmcalPatchMultCorr(const char *name){
  auto *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    std::cerr << "No analysis manager found. Returning ..." << std::endl;
    return nullptr;
  }

  auto *patchtask = new AliAnalysisTaskEmcalPatchMultCorr(Form("PatchMultCorrTask_%s", name));
  mgr->AddTask(patchtask);

  std::stringstream outfilename;
  outfilename << mgr->GetCommonFileName() << ":PatchMultCorrResults_" << name;
  mgr->ConnectInput(patchtask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(patchtask, 1, mgr->CreateContainer(Form("PatchMultCorrHists_%s", name), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename.str().data()));
  
  return patchtask;
}