#include "AliMultEventClassifierTask.h"
#include "AliAODMultEventClass.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliESDEvent.h"

//____________________________________________________________________
AliMultEventClassifierTask::AliMultEventClassifierTask()
  : AliAnalysisTaskSE(),
    fClassifier(),
    fData(),
    fList(0)
{}

//____________________________________________________________________
AliMultEventClassifierTask::AliMultEventClassifierTask(const char* name)
  : AliAnalysisTaskSE(name),
    fClassifier("dummy"),
    fData(),
    fList(0)
{
  DefineOutput(1,TList::Class());
}
//____________________________________________________________________
AliMultEventClassifierTask::AliMultEventClassifierTask(const
						       AliMultEventClassifierTask& o)
  : AliAnalysisTaskSE(o),
    fClassifier(o.fClassifier),
    fData(o.fData),
    fList(0)
{
  if (o.fList) fList = static_cast<TList*>(o.fList->Clone());
}
//____________________________________________________________________
AliMultEventClassifierTask&
AliMultEventClassifierTask::operator=(const AliMultEventClassifierTask& o)
{
  if (&o == this) return *this;

  AliAnalysisTaskSE::operator=(o);

  fClassifier = o.fClassifier;
  fData       = o.fData;

  if (fList) { delete fList; fList = 0; }
  if (o.fList) fList = static_cast<TList*>(o.fList->Clone());

  if (fList) fClassifier.CreateOutputObjects(fList);

  return *this;
}

//____________________________________________________________________
Bool_t
AliMultEventClassifierTask::Connect(const char* sumFile)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskForwardMult", "No analysis manager to connect to.");
    return false;
  }   

  // Add to the manager 
  mgr->AddTask(this);
  
  // Create and connect output containers 
  TString sumOut;
  if (sumFile && sumFile[0] != '\0') sumOut = sumFile;
  // If the string is null or 'default' connect to standard output file 
  if (sumOut.IsNull() || sumOut.EqualTo("default", TString::kIgnoreCase)) 
    sumOut = AliAnalysisManager::GetCommonFileName();

  // Always connect input 
  mgr->ConnectInput(this, 0, mgr->GetCommonInputContainer());

  // Connect sum list unless the output 'none' is specified
  if (!sumOut.EqualTo("none", TString::kIgnoreCase)) {
    TString sumName(Form("%s%s", GetName(), "Sums"));
    AliAnalysisDataContainer* sumCon = 
      mgr->CreateContainer(sumName, TList::Class(), 
			   AliAnalysisManager::kOutputContainer, sumOut);
    mgr->ConnectOutput(this, 1, sumCon);
  }
  
  return true;
}

//____________________________________________________________________
void
AliMultEventClassifierTask::UserCreateOutputObjects()
{
  fList = new TList;
  fList->SetOwner();
  fList->SetName("MultClass");

  fClassifier.CreateOutputObjects(fList);

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());

  if (ah) {
    TObject* data = &fData;
    ah->AddBranch(fData.ClassName(), &data);
  }
  
  PostData(1, fList);
}
//____________________________________________________________________
void
AliMultEventClassifierTask::UserExec(Option_t*)
{
  fData.Clear();

  LoadBranches();

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliWarning("No ESD event found for input event");
    return;
  }

  fClassifier.Process(esd, &fData);

  AliAnalysisManager* am = AliAnalysisManager::GetAnalysisManager();
  AliAODHandler*      ah = 
    dynamic_cast<AliAODHandler*>(am->GetOutputEventHandler());
  if (ah) ah->SetFillAOD(kTRUE);
  
  PostData(1,fList);
}
//____________________________________________________________________
void
AliMultEventClassifierTask::Terminate(Option_t*)
{}
//____________________________________________________________________
//
// EOF
// 
