/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <TClonesArray.h>
#include <THashList.h>
#include <THistManager.h>

#include "AliAnalysisUtils.h"
#include "AliEMCALTriggerPatchInfo.h"
#include <AliEmcalTriggerRejectionMaker.h>
#include "AliInputEventHandler.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::AliEmcalTriggerRejectionMaker)
/// \endcond

using namespace PWG::EMCAL;

AliEmcalTriggerRejectionMaker::AliEmcalTriggerRejectionMaker() :
    AliAnalysisTaskEmcal(),
    fHistos(nullptr),
    fSelectTrigger(AliVEvent::kINT7)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliEmcalTriggerRejectionMaker::AliEmcalTriggerRejectionMaker(const char *name) :
    AliAnalysisTaskEmcal(name, kTRUE),
    fHistos(nullptr),
    fSelectTrigger(AliVEvent::kINT7)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetMakeGeneralHistograms(true);
}

AliEmcalTriggerRejectionMaker::~AliEmcalTriggerRejectionMaker() {
}

void AliEmcalTriggerRejectionMaker::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if(!fAliAnalysisUtils) fAliAnalysisUtils = new AliAnalysisUtils;

  fHistos = new THistManager("histMaxPatch");
  fHistos->CreateTH1("hTrueEventCount", "Maximum energy patch in the event", 1, 0.5, 1.5);

  const std::map<TString, TString> triggers {
    {"EGAOffline", "offline EGA"}, {"EJEOffline", "offline EJE"}, {"DGAOffline", "offline DGA"}, {"DJEOffline", "offline DJE"},
    {"EGARecalc", "recalc EGA"}, {"EJERecalc", "recalc EJE"}, {"DGARecalc", "recalc DGA"}, {"DJERecalc", "recalc DJE"},
    {"EG1Online", "online EG1"}, {"EG2Online", "online EG2"}, {"DG1Online", "online DG1"}, {"DG2Online", "online DG2"},
    {"EJ1Online", "online EJ1"}, {"EJ2Online", "online EJ2"}, {"DJ1Online", "online DJ1"}, {"DJ2Online", "online DJ2"},
  };
  // Calibrated FEE energy
  for(const auto &t : triggers)
    fHistos->CreateTH1("hPatchEnergyMax" + t.first, "Energy spectrum of the maximum " + t.second + " patch", 2000, 0., 200.);

  // Online ADC counts
  for(const auto &t : triggers)
    fHistos->CreateTH1("hPatchADCMax" + t.first, "ADC spectrum of the maximum " + t.second + " patch", 2049, -0.5, 2048.5);

  // ADC vs energy
  for(const auto &t : triggers)
    fHistos->CreateTH2("hPatchADCvsEnergyMax" + t.first, "ADC vs. Energy of the maximum " + t.second + " patch", 300, 0., 3000, 200, 0., 200.);

  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);

  PostData(1, fOutput);
}

Bool_t AliEmcalTriggerRejectionMaker::IsEventSelected(){
  AliDebugStream(2) << GetName() << ": Using custom event selection method" << std::endl;
  if(!fTriggerPatchInfo){
    AliErrorStream() << GetName() << ": Trigger patch container not found but required" << std::endl;
    return false;
  }
  if(!(fInputHandler->IsEventSelected() & fSelectTrigger)) return false;
  if(fTriggerPattern.Length()){
    TString triggerstring = InputEvent()->GetFiredTriggerClasses();
    if(!triggerstring.Contains(fTriggerPattern)) return false;
  }
  AliDebugStream(3) << GetName() << "Event is selected for the given trigger" << std::endl;

  // Generall event quality cuts
  // The vertex cut also contains cuts on the number
  // of contributors and the position in z
  AliDebugStream(3) << GetName() << ": Applying vertex selection" << std::endl;
  if(fAliAnalysisUtils){
    if(!fAliAnalysisUtils->IsVertexSelected2013pA(InputEvent())) return false;
    if(fAliAnalysisUtils->IsPileUpEvent(InputEvent())) return false;
    AliDebugStream(3) << GetName() << ": Vertex selection passed" << std::endl;
  }

  AliDebugStream(2) << GetName() << "Event selected" << std::endl;
  return true;
}

Bool_t AliEmcalTriggerRejectionMaker::Run(){
  fHistos->FillTH1("hTrueEventCount", 1);

  const AliEMCALTriggerPatchInfo *currentpatch(nullptr),
       *maxOfflineEGA(nullptr),
       *maxOfflineEJE(nullptr),
       *maxOfflineDGA(nullptr),
       *maxOfflineDJE(nullptr),
       *maxRecalcEGA(nullptr),
       *maxRecalcEJE(nullptr),
       *maxRecalcDGA(nullptr),
       *maxRecalcDJE(nullptr),
       *maxOnlineEG1(nullptr),
       *maxOnlineEG2(nullptr),
       *maxOnlineDG1(nullptr),
       *maxOnlineDG2(nullptr),
       *maxOnlineEJ1(nullptr),
       *maxOnlineEJ2(nullptr),
       *maxOnlineDJ1(nullptr),
       *maxOnlineDJ2(nullptr);

  // Find the maximum patch for each cathegory
  for(auto patchiter : *fTriggerPatchInfo){
    currentpatch = static_cast<AliEMCALTriggerPatchInfo *>(patchiter);

    // Offline patches - make cut on energy
    if(currentpatch->IsOfflineSimple()){
      if(currentpatch->IsGammaHighSimple()){
        if(currentpatch->IsEMCal()){
          if(!maxOfflineEGA || currentpatch->GetPatchE() > maxOfflineEGA->GetPatchE())
            maxOfflineEGA = currentpatch;
        } else {
          if(!maxOfflineDGA || currentpatch->GetPatchE() > maxOfflineDGA->GetPatchE())
            maxOfflineDGA = currentpatch;
        }
      } else if(currentpatch->IsJetHighSimple()) {
        if(currentpatch->IsEMCal()){
          if(!maxOfflineEJE || currentpatch->GetPatchE() > maxOfflineEJE->GetPatchE())
            maxOfflineEJE = currentpatch;
        } else {
          if(!maxOfflineDJE || currentpatch->GetPatchE() > maxOfflineDJE->GetPatchE())
            maxOfflineDJE = currentpatch;
        }
      }
    }

    // Recalc patches - make cut on FastOR ADC
    if(currentpatch->IsRecalc()){
      if(currentpatch->IsGammaHighRecalc()){
        if(currentpatch->IsEMCal()){
          if(!maxRecalcEGA || currentpatch->GetADCAmp() > maxRecalcEGA->GetADCAmp())
            maxRecalcEGA = currentpatch;
        } else {
          if(!maxRecalcDGA || currentpatch->GetADCAmp() > maxRecalcDGA->GetADCAmp())
            maxRecalcDGA = currentpatch;
        }
      } else if(currentpatch->IsJetHighSimple()) {
        if(currentpatch->IsEMCal()){
          if(!maxRecalcEJE || currentpatch->GetADCAmp() > maxRecalcEJE->GetADCAmp())
            maxRecalcEJE = currentpatch;
        } else {
          if(!maxRecalcDJE || currentpatch->GetADCAmp() > maxRecalcDJE->GetADCAmp())
            maxRecalcDJE = currentpatch;
        }
      }
    }

    // Online patches
    if(currentpatch->IsGammaHigh() || currentpatch->IsGammaLow()){
      if(currentpatch->IsEMCal()){
        if(currentpatch->IsGammaHigh()){
          if(!maxOnlineEG1 || currentpatch->GetADCAmp() > maxOnlineEG1->GetADCAmp())
            maxOnlineEG1 = currentpatch;
        }
        if(currentpatch->IsGammaLow()){
          if(!maxOnlineEG2 || currentpatch->GetADCAmp() > maxOnlineEG2->GetADCAmp())
            maxOnlineEG2 = currentpatch;
        }
      } else {
        if(currentpatch->IsGammaHigh()){
          if(!maxOnlineDG1 || currentpatch->GetADCAmp() > maxOnlineDG1->GetADCAmp())
            maxOnlineDG1 = currentpatch;
        }
        if(currentpatch->IsGammaLow()){
          if(!maxOnlineDG2 || currentpatch->GetADCAmp() > maxOnlineDG2->GetADCAmp())
            maxOnlineDG2 = currentpatch;
        }
      }
    }

    if(currentpatch->IsJetHigh() || currentpatch->IsJetLow()){
      if(currentpatch->IsEMCal()){
        if(currentpatch->IsJetHigh()){
          if(!maxOnlineEJ1 || currentpatch->GetADCAmp() > maxOnlineEJ1->GetADCAmp())
            maxOnlineEJ1 = currentpatch;
        }
        if(currentpatch->IsJetLow()){
          if(!maxOnlineEJ2 || currentpatch->GetADCAmp() > maxOnlineEJ2->GetADCAmp())
            maxOnlineEJ2 = currentpatch;
        }
      } else {
        if(currentpatch->IsJetHigh()){
          if(!maxOnlineDG1 || currentpatch->GetADCAmp() > maxOnlineDJ1->GetADCAmp())
            maxOnlineDG1 = currentpatch;
        }
        if(currentpatch->IsJetLow()){
          if(!maxOnlineDJ2 || currentpatch->GetADCAmp() > maxOnlineDJ2->GetADCAmp())
            maxOnlineDJ2 = currentpatch;
        }
      }
    }
  }

  std::function<void (const AliEMCALTriggerPatchInfo *, const TString &)> FillHistos = [this](const AliEMCALTriggerPatchInfo * testpatch, const TString & triggername){
    fHistos->FillTH1("hPatchEnergyMax" + triggername, testpatch ? testpatch->GetPatchE() : 0.);
    fHistos->FillTH1("hPatchADCMax" + triggername, testpatch ? testpatch->GetADCAmp() : 0.);
    fHistos->FillTH2("hPatchADCvsEnergyMax" + triggername, testpatch ? testpatch->GetADCAmp() : 0, testpatch ? testpatch->GetPatchE() : 0.);
  };

  FillHistos(maxOfflineEGA, "EGAOffline");
  FillHistos(maxOfflineEJE, "EJEOffline");
  FillHistos(maxOfflineDGA, "DGAOffline");
  FillHistos(maxOfflineDJE, "DJEOffline");
  FillHistos(maxRecalcEGA, "EGARecalc");
  FillHistos(maxRecalcEJE, "EJERecalc");
  FillHistos(maxRecalcDGA, "DGARecalc");
  FillHistos(maxRecalcDJE, "DJERecalc");
  FillHistos(maxOnlineEG1, "EG1Online");
  FillHistos(maxOnlineEG2, "EG2Online");
  FillHistos(maxOnlineDG1, "DG1Online");
  FillHistos(maxOnlineDG2, "DG2Online");
  FillHistos(maxOnlineEJ1, "EJ1Online");
  FillHistos(maxOnlineEJ2, "EJ2Online");
  FillHistos(maxOnlineDJ1, "DJ1Online");
  FillHistos(maxOnlineDJ2, "DJ2Online");

  return true;

} /* namespace PWGJE::EMCALJetTasks */
