#include "AliGFWFilterTask.h"

AliGFWFilterTask::AliGFWFilterTask():
 AliAnalysisTaskSE(),
 fAOD(0),
 fOutList(0),
 fFilter(0),
 fEmbedES(kTRUE),
 fAddQA(kTRUE),
 fTrigger(0),
 fExtendV0MAcceptance(kFALSE)
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
 fExtendV0MAcceptance(kFALSE)
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
    fFilter->CreateCutMasks();
    if(fEmbedES) {
      fEventCuts = new AliEventCuts();
      fEventCuts->AddQAplotsToList(fOutList,kTRUE);
      fEventCuts->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
      if(fTrigger) fEventCuts->OverrideCentralityFramework(fTrigger);
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
void AliGFWFilterTask::NotifyRun() {
  if(fEmbedES) {
    AliAODEvent *lEv = dynamic_cast<AliAODEvent*>(InputEvent());
    fEventCuts->AcceptEvent(lEv);
    fEventCuts->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
    if(fTrigger) {
        fEventCuts->OverrideCentralityFramework(fTrigger);
    };
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
