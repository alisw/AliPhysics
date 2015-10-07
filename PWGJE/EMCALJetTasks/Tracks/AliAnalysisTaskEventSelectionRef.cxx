/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <TClonesArray.h>
#include <THashList.h>
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliEmcalTriggerPatchInfo.h"
#include "AliEMCalHistoContainer.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"

#include "AliAnalysisTaskEventSelectionRef.h"

#if __cplusplus < 201103L
#ifndef nullptr
#define nullptr NULL
#endif
#endif

ClassImp(EMCalTriggerPtAnalysis::AliAnalysisTaskEventSelectionRef)

namespace EMCalTriggerPtAnalysis {

AliAnalysisTaskEventSelectionRef::AliAnalysisTaskEventSelectionRef():
  AliAnalysisTaskSE(),
  fAnalysisUtils(nullptr),
  fHistos(nullptr)
{
}

AliAnalysisTaskEventSelectionRef::AliAnalysisTaskEventSelectionRef(const char *name):
  AliAnalysisTaskSE(name),
  fAnalysisUtils(nullptr),
  fHistos(nullptr)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEventSelectionRef::~AliAnalysisTaskEventSelectionRef() {
  if(fAnalysisUtils) delete fAnalysisUtils;
}

void AliAnalysisTaskEventSelectionRef::UserCreateOutputObjects(){
  fAnalysisUtils = new AliAnalysisUtils;

  fHistos = new AliEMCalHistoContainer("Ref");

  TString triggers[6] = {"MB", "EMC7", "EJ1", "EJ2", "EG1", "EG2"};
  for(TString *trgit = triggers; trgit < triggers + sizeof(triggers)/sizeof(TString); ++trgit){
    fHistos->CreateTH1(Form("hEventCount%sBeforeEventSelection", trgit->Data()), Form("Event count for trigger %s before event selection", trgit->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hEventCount%sBeforeOfflineTrigger", trgit->Data()), Form("Event count for trigger %s before offline selection", trgit->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hEventCount%sAfterOfflineTrigger", trgit->Data()), Form("Event count for trigger %s before offline selection", trgit->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hVertexTrigger%sBeforeEventSelection", trgit->Data()), Form("Vertex Distribution for trigger %s before event selection", trgit->Data()), 400, -40, 40);
    fHistos->CreateTH1(Form("hVertexTrigger%sBeforeOfflineTrigger", trgit->Data()), Form("Vertex Distribution for trigger %s before offline trigger", trgit->Data()), 400, -40, 40);
    fHistos->CreateTH1(Form("hVertexTrigger%sAfterOfflineTrigger", trgit->Data()), Form("Vertex Distribution for trigger %s after offline trigger", trgit->Data()), 400, -40, 40);
  }

  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEventSelectionRef::UserExec(Option_t *){
  // Select event
  TClonesArray *triggerpatches = static_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
  TString triggerstring = fInputEvent->GetFiredTriggerClasses();
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1"),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2"),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1"),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2"),
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("CEMC7");
  if(!(isMinBias || isEMC7 || isEG1 || isEG2 || isEJ1 || isEJ2)) return;
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return;
  if(fInputEvent->IsA() == AliESDEvent::Class() && fAnalysisUtils->IsFirstEventInChunk(fInputEvent)) return;
  bool isSelected  = kTRUE;
  if(!fAnalysisUtils->IsVertexSelected2013pA(fInputEvent)) isSelected = kFALSE;       // Apply new vertex cut
  if(fAnalysisUtils->IsPileUpEvent(fInputEvent)) isSelected = kFALSE;       // Apply new vertex cut
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) isSelected = kFALSE;

  // Fill Event counter and reference vertex distributions for the different trigger classes
  if(isMinBias){
    FillEventCounterHists("MB", vtx->GetZ(), isSelected, true);
  }
  if(isEMC7){
    FillEventCounterHists("EMC7", vtx->GetZ(), isSelected, IsOfflineSelected(kCPREJ1, triggerpatches));
  }
  if(isEJ2){
    FillEventCounterHists("EJ2", vtx->GetZ(), isSelected, IsOfflineSelected(kCPREJ2, triggerpatches));
  }
  if(isEJ1){
    FillEventCounterHists("EJ1", vtx->GetZ(), isSelected, IsOfflineSelected(kCPREG1, triggerpatches));
  }
  if(isEG2){
    FillEventCounterHists("EG2", vtx->GetZ(), isSelected, IsOfflineSelected(kCPREG2, triggerpatches));
  }
  if(isEG1){
    FillEventCounterHists("EG1", vtx->GetZ(), isSelected, IsOfflineSelected(kCPREL0, triggerpatches));
  }

  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEventSelectionRef::FillEventCounterHists(
    const char *triggerclass,
    double vtxz,
    bool isSelected,
    bool isOfflineSelected
)
{
  // Fill reference distribution for the primary vertex before any z-cut
  fHistos->FillTH1(Form("hVertexTrigger%sBeforeEventSelection", triggerclass), vtxz);
  fHistos->FillTH1(Form("hEventCount%sBeforeEventSelection", triggerclass), 1.);
  if(isSelected){
    // Fill Event counter and reference vertex distributions after event selection
    fHistos->FillTH1(Form("hEventCount%sBeforeOfflineTrigger", triggerclass), 1);
    fHistos->FillTH1(Form("hVertexTrigger%sBeforeOfflineTrigger", triggerclass), vtxz);
    if(isOfflineSelected){
      fHistos->FillTH1(Form("hEventCount%sAfterOfflineTrigger", triggerclass), 1);
      fHistos->FillTH1(Form("hVertexTrigger%sAfterOfflineTrigger", triggerclass), vtxz);
    }
  }
}

Bool_t AliAnalysisTaskEventSelectionRef::IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const {
  if(fOfflineEnergyThreshold[trgcls] < 0) return true;
  bool isSingleShower = ((trgcls == kCPREL0) || (trgcls == kCPREG1) || (trgcls == kCPREG2));
  int nfound = 0;
  AliEmcalTriggerPatchInfo *patch = NULL;
  for(TIter patchIter = TIter(triggerpatches).Begin(); patchIter != TIter::End(); ++patchIter){
    patch = static_cast<AliEmcalTriggerPatchInfo *>(*patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if(isSingleShower){
     if(!patch->IsGammaLowSimple()) continue;
    } else {
      if(!patch->IsJetLowSimple()) continue;
    }
    if(patch->GetPatchE() > fOfflineEnergyThreshold[trgcls]) nfound++;
  }
  return nfound > 0;
}



} /* namespace EMCalTriggerPtAnalysis */
