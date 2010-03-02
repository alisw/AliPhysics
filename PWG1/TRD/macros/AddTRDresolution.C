#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/macros/AliTRDperformanceTrain.h"
#include "PWG1/TRD/AliTRDresolution.h"
#include "PWG1/TRD/AliTRDclusterResolution.h"
#include "PWG1/TRD/AliTRDalignmentTask.h"
#endif

#include "PWG1/TRD/macros/helper.C"
void AddTRDresolution(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  Int_t map = ParseOptions(trd);

  //AliLog::SetClassDebugLevel("AliTRDresolution", 5);  
  AliTRDresolution *task(NULL);
  if(!TSTBIT(map, kResolution)) return;
  mgr->AddTask(task = new AliTRDresolution((char*)"TRDresolution"));
  task->SetMCdata(mgr->GetMCtruthEventHandler());
  task->SetPostProcess(kFALSE);
  task->SetDebugLevel(0);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());  
  mgr->ConnectInput(task, 1, ci[0]);
  mgr->ConnectOutput(task,1, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.Performance.root"));

  // Create output containers for calibration tasks
  const Int_t nc = 4;
  const Char_t *cn[nc] = {"Cl", "Trklt", "MC_Cl", "MC_Trklt"}; 
  AliAnalysisDataContainer *co[] = {0x0, 0x0, 0x0, 0x0};
  for(Int_t ic = 0; ic<nc; ic++){
    co[ic] = mgr->CreateContainer(Form("%s%s", task->GetName(), cn[ic]), TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(task, 1+ic, co[ic]);
  }

  // Cluster Error Parameterization
  if(TSTBIT(map, kClErrParam)){
    //AliLog::SetClassDebugLevel("AliTRDclusterResolution", 5);  

    AliTRDclusterResolution *taskCl(NULL);
    mgr->AddTask(taskCl = new AliTRDclusterResolution((char*)"ESD", "ESD Cluster error parameterization"));
    taskCl->SetExB();
    mgr->ConnectInput(taskCl,  0, mgr->GetCommonInputContainer());  
    mgr->ConnectInput(taskCl,  1, co[0]);
    mgr->ConnectOutput(taskCl, 1, mgr->CreateContainer(taskCl->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.CalibClErrParam.root"));

    mgr->AddTask(taskCl = new AliTRDclusterResolution((char*)"MC", "MC Cluster error parameterization"));
    taskCl->SetExB();
    mgr->ConnectInput(taskCl,  0, mgr->GetCommonInputContainer());  
    mgr->ConnectInput(taskCl,  1, co[2]);
    mgr->ConnectOutput(taskCl, 1, mgr->CreateContainer(taskCl->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.CalibClErrParam.root"));
  }

  // TRD alignment
  if(TSTBIT(map, kAlignment)){
    AliTRDalignmentTask *taskAl = 0x0;
    mgr->AddTask(taskAl = new AliTRDalignmentTask((char*)"TRDalignment"));
    taskAl->SetDebugLevel(0);
    mgr->ConnectInput(taskAl,  0, mgr->GetCommonInputContainer());  
    mgr->ConnectInput(taskAl,  1, ci[0]);  
    mgr->ConnectOutput(taskAl, 1, mgr->CreateContainer(Form("h%s", taskAl->GetName()), TObjArray::Class(), AliAnalysisManager::kExchangeContainer));
    mgr->ConnectOutput(taskAl, 2, mgr->CreateContainer(task->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", task->GetName())));
  }
}

