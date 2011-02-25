#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "TError.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "TRD/AliTRDgeometry.h"
#include "PWG1/TRD/AliTRDpwg1Helper.h"
#include "PWG1/TRD/AliTRDresolution.h"
#include "PWG1/TRD/AliTRDclusterResolution.h"
#include "PWG1/TRD/AliTRDalignmentTask.h"
#endif

void AddTRDresolution(AliAnalysisManager *mgr, Int_t map, AliAnalysisDataContainer **ci)
{
  Info("AddTRDresolution", Form("[0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\"", ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName(), ci[3]->GetName()));
  AliAnalysisDataContainer *evInfoContainer = ci[3];

  //AliLog::SetClassDebugLevel("AliTRDrecoTask", 2);
  //AliLog::SetClassDebugLevel("AliTRDresolution", 2);
  AliTRDresolution *res(NULL);
  const Char_t *suffix[]={"", "SA", "K"};
  for(Int_t itq=0; itq<1/*3*/; itq++){
    mgr->AddTask(res = new AliTRDresolution(Form("TRDresolution%s", suffix[itq])));
    res->SetMCdata(mgr->GetMCtruthEventHandler());
    res->SetPostProcess(kFALSE);
    res->SetDebugLevel(0);
    if(itq==0) res->SetSegmentationLevel(AliTRDresolution::kDetector);
    // use these settings if you know what you are doing !
    //res->SetTrackRefit(); 
    //res->SetPtThreshold(0.);
    res->SetNameId(suffix[itq]);
    mgr->ConnectInput(res, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container 
    mgr->ConnectInput(res, 1, ci[itq]);                        // conect track info container
    mgr->ConnectInput(res, 2, evInfoContainer);                // conect event info container
    mgr->ConnectOutput(res,1, mgr->CreateContainer(res->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance", mgr->GetCommonFileName())));
  
    // Create output containers for calibration tasks
    AliAnalysisDataContainer *co(NULL);
    co = mgr->CreateContainer(Form("%sCl2Trk%s", res->GetName(), suffix[itq]), TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(res, AliTRDresolution::kClToTrk, co);
    co = mgr->CreateContainer(Form("%sCl2MC%s", res->GetName(), suffix[itq]), TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(res, AliTRDresolution::kClToMC, co);
    
    TObjArray *coa = mgr->GetContainers();
    // Cluster Error Parameterization
    if(TESTBIT(map, AliTRDpwg1Helper::kClErrParam)){
      AliTRDclusterResolution *taskCl(NULL);
      AliLog::SetClassDebugLevel("AliTRDclusterResolution", 2);
      for(Int_t idet(10); idet<11/*AliTRDgeometry::kNdet*/; idet++){
        mgr->AddTask(taskCl = new AliTRDclusterResolution(Form("ClErrCalib%03d", idet)));
        taskCl->SetCalibrationRegion(idet);
        taskCl->SetDebugLevel(0);
        mgr->ConnectInput(taskCl,  0, mgr->GetCommonInputContainer());  // connect main (ESD) container
        mgr->ConnectInput(taskCl,  1, (AliAnalysisDataContainer*)coa->FindObject(Form("%sCl2Trk%s", res->GetName(), suffix[itq])));
        mgr->ConnectInput(taskCl,  2, evInfoContainer);
        mgr->ConnectOutput(taskCl, 1, mgr->CreateContainer(taskCl->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Calibration", mgr->GetCommonFileName())));
        if(mgr->GetMCtruthEventHandler()){
          mgr->AddTask(taskCl = new AliTRDclusterResolution(Form("ClErrCalibMC%03d", idet)));
          taskCl->SetCalibrationRegion(idet);
          taskCl->SetDebugLevel(0);
          mgr->ConnectInput(taskCl,  0, mgr->GetCommonInputContainer());  
          mgr->ConnectInput(taskCl,  1, (AliAnalysisDataContainer*)coa->FindObject(Form("%sCl2MC%s", res->GetName(), suffix[itq])));
          mgr->ConnectInput(taskCl,  2, evInfoContainer);
          mgr->ConnectOutput(taskCl, 1, mgr->CreateContainer(taskCl->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Calibration", mgr->GetCommonFileName())));
        }
      }
    }
  }

  // TRD alignment
  if(TESTBIT(map, AliTRDpwg1Helper::kAlignment)){
    AliTRDalignmentTask *taskAlign(NULL);
    mgr->AddTask(taskAlign = new AliTRDalignmentTask((char*)"TRDalignment"));
    taskAlign->SetDebugLevel(0);
    //AliLog::SetClassDebugLevel("AliTRDalignmentTask", 5);  
    mgr->ConnectInput(taskAlign,  0, mgr->GetCommonInputContainer());  
    mgr->ConnectInput(taskAlign,  1, ci[0]);
    mgr->ConnectInput(taskAlign,  2, evInfoContainer);
    mgr->ConnectOutput(taskAlign, 1, mgr->CreateContainer(Form("h%s", taskAlign->GetName()), TObjArray::Class(), AliAnalysisManager::kExchangeContainer));
    mgr->ConnectOutput(taskAlign, 2, mgr->CreateContainer(taskAlign->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Alignment",mgr->GetCommonFileName())));
  }
}

