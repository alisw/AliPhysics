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
void AddTRDresolution(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **ci)
{
  Int_t map = ParseOptions(trd);
  if(!TSTBIT(map, kResolution)) return;
  printf("AddTRDresolution <- [0]=\"%s\" [1]=\"%s\" [2]=\"%s\"\n", ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName());

  //AliLog::SetClassDebugLevel("AliTRDrecoTask", 2);
  //AliLog::SetClassDebugLevel("AliTRDresolution", 2);
  AliTRDresolution *res(NULL);
  const Char_t *suffix[]={"", "SA", "K"};
  for(Int_t itq=0; itq<1/*3*/; itq++){
    mgr->AddTask(res = new AliTRDresolution(Form("TRDresolution%s", suffix[itq])));
    res->SetMCdata(mgr->GetMCtruthEventHandler());
    res->SetPostProcess(kFALSE);
    res->SetDebugLevel(0);
    //if(itq==0) res->SetSegmentationLevel(1);
    res->SetNameId(suffix[itq]);
    mgr->ConnectInput(res, 0, mgr->GetCommonInputContainer());  
    mgr->ConnectInput(res, 1, ci[itq]);
    mgr->ConnectOutput(res,1, mgr->CreateContainer(res->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.Performance.root"));
  
    // Create output containers for calibration tasks
    AliAnalysisDataContainer *co(NULL);
    co = mgr->CreateContainer(Form("%sCl2Trk%s", res->GetName(), suffix), TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(res, AliTRDresolution::kClToTrk, co);
    co = mgr->CreateContainer(Form("%sCl2MC%s", res->GetName(), suffix), TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(res, AliTRDresolution::kClToMC, co);
//     co = mgr->CreateContainer(Form("%sTrklt2Trk%s", res->GetName(), suffix), TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
//     mgr->ConnectOutput(res, AliTRDresolution::kTrkltToTrk, co);
//     co = mgr->CreateContainer(Form("%sTrklt2MC%s", res->GetName(), suffix), TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
//     mgr->ConnectOutput(res, AliTRDresolution::kTrkltToMC, co);
    
    TObjArray *coa = mgr->GetContainers();
    // Cluster Error Parameterization
    if(TSTBIT(map, kClErrParam)){
      AliTRDclusterResolution *taskCl(NULL);
/*      mgr->AddTask(taskCl = new AliTRDclusterResolution((char*)"ClErrCalibESD"));
      taskCl->SetExB();
      taskCl->SetDebugLevel(0);
      //AliLog::SetClassDebugLevel("AliTRDclusterResolution", 5);  
  
      mgr->ConnectInput(taskCl,  0, mgr->GetCommonInputContainer()); 
      mgr->ConnectInput(taskCl,  1, (AliAnalysisDataContainer*)coa->FindObject(Form("%sCl2Trk%s", res->GetName(), suffix)));
      mgr->ConnectOutput(taskCl, 1, mgr->CreateContainer(taskCl->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.CalibClErrParam.root"));*/
  
      AliLog::SetClassDebugLevel("AliTRDclusterResolution", 5);  
      mgr->AddTask(taskCl = new AliTRDclusterResolution((char*)"ClErrCalibMC"));
      taskCl->SetExB();
      taskCl->SetDebugLevel(0);
      mgr->ConnectInput(taskCl,  0, mgr->GetCommonInputContainer());  
      mgr->ConnectInput(taskCl,  1, (AliAnalysisDataContainer*)coa->FindObject(Form("%sCl2MC%s", res->GetName(), suffix)));
      mgr->ConnectOutput(taskCl, 1, mgr->CreateContainer(taskCl->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.CalibClErrParam.root"));
    }
  }

  // TRD alignment
  if(TSTBIT(map, kAlignment)){
    TObjArray *coa = mgr->GetContainers();
    AliTRDalignmentTask *taskAlign(NULL);
    mgr->AddTask(taskAlign = new AliTRDalignmentTask((char*)"TRDalignment"));
    taskAlign->SetDebugLevel(0);
    //AliLog::SetClassDebugLevel("AliTRDalignmentTask", 5);  
    mgr->ConnectInput(taskAlign,  0, mgr->GetCommonInputContainer());  
    mgr->ConnectInput(taskAlign,  1, (AliAnalysisDataContainer*)coa->FindObject("TRDresolutionCl2Trk"));  
    mgr->ConnectOutput(taskAlign, 1, mgr->CreateContainer(Form("h%s", taskAlign->GetName()), TObjArray::Class(), AliAnalysisManager::kExchangeContainer));
    mgr->ConnectOutput(taskAlign, 2, mgr->CreateContainer(taskAlign->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", taskAlign->GetName())));
  }
}

