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
  fHistos->CreateTH1("hPatchEnergyMaxEGA", "Energy spectrum of the maximum EGA patch", 5000, 0., 500.);
  fHistos->CreateTH1("hPatchEnergyMaxEJE", "Energy spectrum of the maximum EJE patch", 5000, 0., 500.);
  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalMaxPatch::IsEventSelected(){
  AliDebugStream(2) << GetName() << ": Using custom event selection method" << std::endl;
  if(!fTriggerPatchInfo){
    AliErrorStream() << GetName() << ": Trigger patch container not found but required" << std::endl;
    return false;
  }
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;
  if(fTriggerPattern.Length()){
    TString triggerstring = InputEvent()->GetFiredTriggerClasses();
    if(!triggerstring.Contains(fTriggerPattern)) return false;
  }
  AliDebugStream(3) << GetName() << "Event is an INT7 event" << std::endl;

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

  double maxEGA(0), maxEJE(0);
  int nEGA(0), nEJE(0);
  AliEMCALTriggerPatchInfo *currentpatch(nullptr);
  for(TIter patchiter = TIter(this->fTriggerPatchInfo).Begin(); patchiter != TIter::End(); ++patchiter){
    currentpatch = static_cast<AliEMCALTriggerPatchInfo *>(*patchiter);
    if(!currentpatch->IsOfflineSimple()) continue;
    if(currentpatch->IsGammaHighSimple()){
      nEGA++;
      maxEGA = currentpatch->GetPatchE() > maxEGA ? currentpatch->GetPatchE() : maxEGA;
    } else if(currentpatch->IsJetHighSimple()) {
      nEJE++;
      maxEJE = currentpatch->GetPatchE() > maxEJE ? currentpatch->GetPatchE() : maxEJE;
    }
  }

  if(nEGA) fHistos->FillTH1("hPatchEnergyMaxEGA", maxEGA);
  if(nEJE) fHistos->FillTH1("hPatchEnergyMaxEJE", maxEJE);

  return true;
}

} /* namespace EMCalTriggerPtAnalysis */
