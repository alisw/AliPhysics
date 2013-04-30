#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "TError.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWGPP/TRD/AliTRDpwgppHelper.h"
#include "PWGPP/TRD/AliTRDcheckPID.h"
#include "PWGPP/TRD/AliTRDpidRefMaker.h"
#include "PWGPP/TRD/AliTRDpidRefMakerNN.h"
#include "PWGPP/TRD/AliTRDpidRefMakerLQ.h"
#endif

void AddTRDcheckPID(AliAnalysisManager *mgr, Int_t map, AliAnalysisDataContainer **ci, AliAnalysisDataContainer **co)
{
  Info("AddTRDcheckPID", "[0]=\"%s\" [1]=\"%s\" [2]=\"%s\"", 
       ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName());

  AliTRDcheckPID *pid(NULL);
  mgr->AddTask(pid = new AliTRDcheckPID((char*)"TRDcheckPID"));
  //AliLog::SetClassDebugLevel("AliTRDcheckPID", 5);  
  pid->SetDebugLevel(0);
  pid->SetMCdata(mgr->GetMCtruthEventHandler());

  // define PID exchange container
  co[0] = mgr->CreateContainer("InfoPID", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput (pid, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container
  mgr->ConnectInput (pid, 1, ci[0]);                          // connect barrel tracks container
  mgr->ConnectInput (pid, 2, ci[1]);                          // connect event info container
  mgr->ConnectInput (pid, 3, ci[2]);                          // connect V0s container
  mgr->ConnectOutput(pid, 1, mgr->CreateContainer(pid->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
  mgr->ConnectOutput(pid, 2, co[0]);

  if(TESTBIT(map, AliTRDpwgppHelper::kPIDRefMaker)){

    //AliLog::SetClassDebugLevel("AliTRDpidRefMaker", 3);
    //AliLog::SetClassDebugLevel("AliTRDpidRefMakerNN", 3);
    //AliLog::SetClassDebugLevel("AliTRDpidRefMakerLQ", 3);
  
    // TRD pid reference maker NN
    AliTRDpidRefMaker *ref(NULL);
    mgr->AddTask(ref = new AliTRDpidRefMakerNN((char*)"TRDrefMakerNN"));
    ref->SetDebugLevel(3);
    ref->SetMCdata(mgr->GetMCtruthEventHandler());
    ref->SetFriends(kTRUE);
    mgr->ConnectInput( ref, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container
    mgr->ConnectInput( ref, 1, ci[0]);                          // connect barrel tracks container
    mgr->ConnectInput( ref, 2, ci[1]);                          // connect event info container
    mgr->ConnectInput( ref, 3, ci[2]);                          // connect V0s container
    mgr->ConnectInput( ref, 4, co[0]);                          // connect pid Info container
    mgr->ConnectOutput(ref, 1, mgr->CreateContainer("MonitorNN", TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Calibration",mgr->GetCommonFileName())));
    mgr->ConnectOutput(ref, 2, mgr->CreateContainer(ref->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Calibration", mgr->GetCommonFileName())));
  
    // TRD pid reference maker LQ 
    mgr->AddTask(ref = new AliTRDpidRefMakerLQ((char*)"TRDrefMakerLQ"));
    ref->SetDebugLevel(3);
    ref->SetMCdata(mgr->GetMCtruthEventHandler());
    ref->SetFriends(kTRUE);
    mgr->ConnectInput(ref, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container
    mgr->ConnectInput(ref, 1, ci[0]);                          // connect barrel tracks container
    mgr->ConnectInput(ref, 2, ci[1]);                          // connect event info container
    mgr->ConnectInput(ref, 3, ci[2]);                          // connect V0s container
    mgr->ConnectInput(ref, 4, co[0]);                          // connect pid Info container
    mgr->ConnectOutput(ref, 1, mgr->CreateContainer("MonitorLQ", TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Calibration", mgr->GetCommonFileName())));
    mgr->ConnectOutput(ref, 2, mgr->CreateContainer(ref->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Calibration", mgr->GetCommonFileName())));
    mgr->ConnectOutput(ref, 3, mgr->CreateContainer("PDF", TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Calibration", mgr->GetCommonFileName())));
  }
}

