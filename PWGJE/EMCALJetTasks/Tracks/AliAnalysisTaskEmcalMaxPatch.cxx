/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <TClonesArray.h>
#include <THashList.h>
#include <THistManager.h>

#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskEmcalMaxPatch.h"
#include "AliEMCalTriggerWeightHandler.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEmcalMaxPatch)
/// \endcond

namespace EMCalTriggerPtAnalysis {

AliAnalysisTaskEmcalMaxPatch::AliAnalysisTaskEmcalMaxPatch() :
    AliAnalysisTaskEmcal(),
    fWeightHandler(nullptr),
    fHistos(nullptr),
    fSelectTrigger(AliVEvent::kINT7)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
}

AliAnalysisTaskEmcalMaxPatch::AliAnalysisTaskEmcalMaxPatch(const char *name) :
    AliAnalysisTaskEmcal(name, kTRUE),
    fWeightHandler(nullptr),
    fHistos(nullptr),
    fSelectTrigger(AliVEvent::kINT7)
{
  SetCaloTriggerPatchInfoName("EmcalTriggers");
  SetMakeGeneralHistograms(true);
}

AliAnalysisTaskEmcalMaxPatch::~AliAnalysisTaskEmcalMaxPatch() {
}

void AliAnalysisTaskEmcalMaxPatch::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if(!fAliAnalysisUtils) fAliAnalysisUtils = new AliAnalysisUtils;

  fHistos = new THistManager("histMaxPatch");
  fHistos->CreateTH1("hTrueEventCount", "Maximum energy patch in the event", 1, 0.5, 1.5);

  const std::map<std::string, std::string> triggers {
    {"EGAOffline", "offline EGA"}, {"EJEOffline", "offline EJE"}, {"EGARecalc", "recalc EGA"},
    {"EJERecalc", "recalc EJE"}, {"EG1Online", "online EG1"}, {"EG2Online", "online EG2"},
    {"DG1Online", "online DG1"}, {"DG2Online", "online DG2"}, {"EJ1Online", "online EJ1"},
    {"EJ2Online", "online EJ2"}, {"DJ1Online", "online DJ1"}, {"DJ2Online", "online DJ2"},
  };
  // Calibrated FEE energy
  for(const auto &t : triggers)
    fHistos->CreateTH1(Form("hPatchEnergyMax%s", t.first.c_str()), Form("Energy spectrum of the maximum %s patch", t.second.c_str()), 2000, 0., 200.);

  // Online ADC counts
  for(const auto &t : triggers)
    fHistos->CreateTH1(Form("hPatchADCMax%s", t.first.c_str()), Form("ADC spectrum of the maximum %s patch", t.second.c_str()), 2049, -0.5, 2048.5);

  // ADC vs energy
  for(const auto &t : triggers)
    fHistos->CreateTH2(Form("hPatchADCvsEnergyMax%s", t.first.c_str()), Form("ADC vs. Energy of the maximum %s patch", t.second.c_str()), 300, 0., 3000, 200, 0., 200.);

  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalMaxPatch::IsEventSelected(){
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

Bool_t AliAnalysisTaskEmcalMaxPatch::Run(){
  fHistos->FillTH1("hTrueEventCount", 1);

  const AliEMCALTriggerPatchInfo *currentpatch(nullptr),
       *maxOfflineEGA(nullptr),
       *maxOfflineEJE(nullptr),
       *maxRecalcEGA(nullptr),
       *maxRecalcEJE(nullptr),
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
        if(!maxOfflineEGA || currentpatch->GetPatchE() > maxOfflineEGA->GetPatchE())
          maxOfflineEGA = currentpatch;
      } else if(currentpatch->IsJetHighSimple()) {
        if(!maxOfflineEJE || currentpatch->GetPatchE() > maxOfflineEJE->GetPatchE())
          maxOfflineEJE = currentpatch;
      }
    }

    // Recalc patches - make cut on FastOR ADC
    if(currentpatch->IsRecalc()){
      if(currentpatch->IsGammaHighRecalc()){
        if(!maxRecalcEGA || currentpatch->GetADCAmp() > maxRecalcEGA->GetADCAmp())
          maxRecalcEGA = currentpatch;
      } else if(currentpatch->IsJetHighSimple()) {
        if(!maxRecalcEJE || currentpatch->GetADCAmp() > maxRecalcEJE->GetADCAmp())
          maxRecalcEJE = currentpatch;
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

  std::function<void (const AliEMCALTriggerPatchInfo *, const std::string &)> FillHistos = [this](const AliEMCALTriggerPatchInfo * testpatch, const std::string & triggername){
    fHistos->FillTH1(Form("hPatchEnergyMax%s", triggername.c_str()), testpatch ? testpatch->GetPatchE() : 0.);
    fHistos->FillTH1(Form("hPatchADCMax%s", triggername.c_str()), testpatch ? testpatch->GetPatchE() : 0.);
    fHistos->FillTH2(Form("hPatchADCvsEnergyMax%s", triggername.c_str()), testpatch ? testpatch->GetPatchE() : 0., testpatch ? testpatch->GetADCAmp() : 0);
  };

  FillHistos(maxOfflineEGA, "EGAOffline");
  FillHistos(maxOfflineEJE, "EJEOffline");
  FillHistos(maxRecalcEGA, "EGARecalc");
  FillHistos(maxRecalcEJE, "EJERecalc");
  FillHistos(maxOnlineEG1, "EG1Online");
  FillHistos(maxOnlineEG2, "EG2Online");
  FillHistos(maxOnlineDG1, "DG1Online");
  FillHistos(maxOnlineDG2, "DG2Online");
  FillHistos(maxOnlineEJ1, "EJ1Online");
  FillHistos(maxOnlineEJ2, "EJ2Online");
  FillHistos(maxOnlineDJ1, "DJ1Online");
  FillHistos(maxOnlineDJ2, "DJ2Online");

  return true;
}

} /* namespace EMCalTriggerPtAnalysis */
