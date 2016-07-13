/*
 * AliAnalysisTaskEmcalTriggerPosition.cxx
 *
 *  Created on: 11.07.2016
 *      Author: markusfasel
 */
#include <TClonesArray.h>
#include <THistManager.h>
#include <TString.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliEMCALTriggerPatchInfo.h"

#include "AliAnalysisTaskEmcalTriggerPosition.h"

ClassImp(AliAnalysisTaskEmcalTriggerPosition)

AliAnalysisTaskEmcalTriggerPosition::AliAnalysisTaskEmcalTriggerPosition():
  AliAnalysisTaskSE(),
  fHistos(nullptr),
  fThresholdEG1(0.)
{

}

AliAnalysisTaskEmcalTriggerPosition::AliAnalysisTaskEmcalTriggerPosition(const char *name):
  AliAnalysisTaskSE(name),
  fHistos(nullptr),
  fThresholdEG1(0.)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEmcalTriggerPosition::~AliAnalysisTaskEmcalTriggerPosition() {
}

void AliAnalysisTaskEmcalTriggerPosition::UserCreateOutputObjects(){
  fHistos = new THistManager("EG1position");

  fHistos->CreateTH2("posRecalcEG1", "Position of recalc EG1 patches", 48, -0.5, 47.5, 104, -0.5, 103.5);
  fHistos->CreateTH2("posOnlineEG1", "Position of recalc EG1 patches", 48, -0.5, 47.5, 104, -0.5, 103.5);

  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEmcalTriggerPosition::UserExec(Option_t *opt){
  // Select only EG1 events
  TString firedtriggers(InputEvent()->GetFiredTriggerClasses());
  if(!firedtriggers.Contains("EG1")) return;
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kEMCEGA)) return;

  TClonesArray *triggerpatches = dynamic_cast<TClonesArray *>(InputEvent()->FindListObject("EmcalTriggers"));
  if(!triggerpatches)
    AliErrorStream() << "Trigger patch container EmcalTriggers not found in task " << GetName() << std::endl;

  AliEMCALTriggerPatchInfo *currentpatch(nullptr);
  for(TIter patchiter = TIter(triggerpatches).Begin(); patchiter != TIter::End(); ++patchiter){
    currentpatch = static_cast<AliEMCALTriggerPatchInfo *>(*patchiter);
    if(currentpatch->GetPatchSize() != 2) continue;
    // Reject L0 patches
    if(currentpatch->IsLevel0() || currentpatch->IsLevel0Recalc() || currentpatch->IsLevel0Simple()) continue;
    if(currentpatch->GetADCAmp() > this->fThresholdEG1)
      fHistos->FillTH1("posRecalcEG1", currentpatch->GetColStart(), currentpatch->GetRowStart());
    if(currentpatch->IsGammaHigh())
      fHistos->FillTH1("posOnlineEG1", currentpatch->GetColStart(), currentpatch->GetRowStart());
  }

  PostData(1, fHistos->GetListOfHistograms());
}

