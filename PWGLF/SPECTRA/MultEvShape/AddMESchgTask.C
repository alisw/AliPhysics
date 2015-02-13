#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "TError.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWGLF/SPECTRA/MultEvShape/AliMESbaseTask.h"
#include "PWGLF/SPECTRA/MultEvShape/AliMESchgTask.h"
#endif

void AddMESchgTask(AliAnalysisDataContainer **ci, Bool_t mc)
{
  Info("AddMESchgTask", Form("Inputs : [0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\"", 
      ci[0]?ci[0]->GetName():"none", 
      ci[1]?ci[1]->GetName():"none", 
      ci[2]?ci[2]->GetName():"none", 
      ci[3]?ci[3]->GetName():"none"));

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager(); 
  AliMESchgTask *chg = new AliMESchgTask("PID");
  mgr->AddTask(chg);
  chg->SetPostProcess(kFALSE);
  chg->SetMCdata(mc);
  chg->SetDebugLevel(1);
  mgr->ConnectInput(chg, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container
  mgr->ConnectInput(chg, AliMESbaseTask::kEventInfo, ci[0]); // connect event info
  mgr->ConnectInput(chg, AliMESbaseTask::kTracks, ci[1]);    // connect track info container
  if(mc){
    mgr->ConnectInput(chg, AliMESbaseTask::kMCeventInfo, ci[2]); // connect MC event info
    mgr->ConnectInput(chg, AliMESbaseTask::kMCtracks, ci[3]);    // connect MC tracks container
  }
  mgr->ConnectOutput(chg, AliMESbaseTask::kQA, mgr->CreateContainer("chgQA", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:MES", mgr->GetCommonFileName())));
}

