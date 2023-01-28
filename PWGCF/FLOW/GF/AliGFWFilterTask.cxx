#include "AliGFWFilterTask.h"

AliGFWFilterTask::AliGFWFilterTask():
 AliAnalysisTaskSE(),
 fAOD(0),
 fOutList(0),
 fFilter(0),
 fEmbedES(kTRUE),
 fAddQA(kTRUE),
 fTrigger(0),
 fExtendV0MAcceptance(kFALSE),
 fCustomCuts({}),
 fDisableDefaultCuts(kFALSE),
 fStandardChi2Cut(GFWFlags::klTPCchi2PC25),
 fPtMin(0.2),
 fPtMax(3.0),
 fEtaMin(-0.8),
 fEtaMax(0.8)
{}
//_____________________________________________________________________________
AliGFWFilterTask::AliGFWFilterTask(const char* name):
 AliAnalysisTaskSE(name),
 fAOD(0),
 fOutList(0),
 fFilter(0),
 fEmbedES(kTRUE),
 fAddQA(kTRUE),
 fTrigger(0),
 fExtendV0MAcceptance(kFALSE),
 fCustomCuts({}),
 fDisableDefaultCuts(kFALSE),
 fStandardChi2Cut(GFWFlags::klTPCchi2PC25),
 fPtMin(0.2),
 fPtMax(3.0),
 fEtaMin(-0.8),
 fEtaMax(0.8)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliGFWFilterTask::~AliGFWFilterTask()
{}
//_____________________________________________________________________________
void AliGFWFilterTask::UserCreateOutputObjects()
{
    OpenFile(1);
    // creating output lists
    fOutList = new TList();
    fOutList->SetOwner(kTRUE);
    fFilter = new AliGFWFilter();
    if(!fDisableDefaultCuts || !fCustomCuts.size()) fFilter->CreateStandardCutMasks(fStandardChi2Cut); //If no custom cuts are set, then create standard cuts nevertheless
    for(auto lcut: fCustomCuts) fFilter->AddCustomCuts(kFALSE,lcut.first, lcut.second);
    fFilter->SetPt(fPtMin,fPtMax);
    fFilter->SetEta(fEtaMin,fEtaMax);
    if(fEmbedES) {
      fEventCuts = new AliEventCuts();
      fEventCuts->AddQAplotsToList(fOutList,kTRUE);
      fEventCuts->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
      if(fTrigger) fEventCuts->OverrideAutomaticTriggerSelection(fTrigger,kTRUE);
      if(fExtendV0MAcceptance) {
        fEventCuts->OverrideCentralityFramework(1);
        fEventCuts->SetCentralityEstimators("V0M","CL0");
        fEventCuts->SetCentralityRange(0.f,101.f);
      }
      fFilter->SetEventCuts(fEventCuts);
    };
    PostData(1, fOutList);
}
//_____________________________________________________________________________
void AliGFWFilterTask::SetupLWCuts(Bool_t fb96, Bool_t fb768, Bool_t use25Chi)
{
  if(!fb96 && !fb768) return;
  fDisableDefaultCuts = kTRUE;
  GFWFlags::kLocalTrackFlags l_standardTPCChiCut = use25Chi?klTPCchi2PC25:klTPCchi2PC40;
  if(fb96)  AddCustomCuts(klVtxZ10 + klEventCuts, klFB96 + klDCAz20 + klDCAxy2011 + l_standardTPCChiCut + klNTPCcls70);
  if(fb768) AddCustomCuts(klVtxZ10 + klEventCuts, klFB768 + l_standardTPCChiCut + klNTPCcls70);
}
//_____________________________________________________________________________
void AliGFWFilterTask::NotifyRun() {
  if(fEmbedES) {
    AliAODEvent *lEv = dynamic_cast<AliAODEvent*>(InputEvent());
    fEventCuts->AcceptEvent(lEv);
    fEventCuts->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
    if(fTrigger) fEventCuts->OverrideAutomaticTriggerSelection(fTrigger,kTRUE);
    if(fExtendV0MAcceptance) {
      fEventCuts->OverrideCentralityFramework(1);
      fEventCuts->SetCentralityEstimators("V0M","CL0");
      fEventCuts->SetCentralityRange(0.f,101.f);
    }
  }
}
//_____________________________________________________________________________
void AliGFWFilterTask::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;
    fFilter->CheckEvent(fAOD);
    TObject *obA = fAOD->FindListObject("GFWFlags");
    if(!obA)
      fAOD->AddObject(fFilter->GetFlags());
    PostData(1, fOutList);
}
//_____________________________________________________________________________
void AliGFWFilterTask::Terminate(Option_t *)
{}
//_____________________________________________________________________________
