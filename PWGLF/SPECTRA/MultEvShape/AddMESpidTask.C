#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "TError.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWGLF/SPECTRA/MultEvShape/AliMESbaseTask.h"
#include "PWGLF/SPECTRA/MultEvShape/AliMESpidTask.h"
#endif

void AddMESpidTask(AliAnalysisDataContainer **ci, Bool_t mc)
{
  Info("AddMESpidTask", Form("Inputs : [0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\"", 
      ci[0]?ci[0]->GetName():"none", 
      ci[1]?ci[1]->GetName():"none", 
      ci[2]?ci[2]->GetName():"none", 
      ci[3]?ci[3]->GetName():"none"));

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager(); 
  AliMESpidTask *pid = new AliMESpidTask("PID");
  mgr->AddTask(pid);
  pid->SetPostProcess(kFALSE);
  pid->SetMCdata(mc);
  pid->SetDebugLevel(0);
  mgr->ConnectInput(pid, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container
  mgr->ConnectInput(pid, AliMESbaseTask::kEventInfo, ci[0]); // connect event info
  mgr->ConnectInput(pid, AliMESbaseTask::kTracks, ci[1]);    // connect track info container
  if(mc){
    mgr->ConnectInput(pid, AliMESbaseTask::kMCeventInfo, ci[2]); // connect MC event info
    mgr->ConnectInput(pid, AliMESbaseTask::kMCtracks, ci[3]);    // connect MC tracks container
  }
  mgr->ConnectOutput(pid, AliMESbaseTask::kQA, mgr->CreateContainer("pidQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:MES", mgr->GetCommonFileName())));
}

