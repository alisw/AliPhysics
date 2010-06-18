#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/macros/AliTRDperformanceTrain.h"
#include "PWG1/TRD/AliTRDcheckPID.h"
#include "PWG1/TRD/AliTRDpidRefMaker.h"
#include "PWG1/TRD/AliTRDpidRefMakerNN.h"
#include "PWG1/TRD/AliTRDpidRefMakerLQ.h"
#endif

#include "PWG1/TRD/macros/helper.C"
void AddTRDcheckPID(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  Int_t map = ParseOptions(trd);
  if(!TSTBIT(map, kCheckPID)) return;
  printf("AddTRDcheckPID <- [0]=\"%s\" [1]=\"%s\"\n", ci[0]->GetName(), ci[1]->GetName());

  AliTRDcheckPID *pid(NULL);
  mgr->AddTask(pid = new AliTRDcheckPID((char*)"checkPID"));
  //AliLog::SetClassDebugLevel("AliTRDcheckPID", 5);  
  pid->SetDebugLevel(0);
  pid->SetMCdata(mgr->GetMCtruthEventHandler());

  // define PID exchange container
  AliAnalysisDataContainer *ce = mgr->CreateContainer("InfoPID", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
  mgr->ConnectInput (pid, 0, mgr->GetCommonInputContainer());
  mgr->ConnectInput (pid, 1, ci[0]);
  mgr->ConnectInput (pid, 2, ci[1]);
  mgr->ConnectOutput(pid, 1, mgr->CreateContainer(pid->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
  mgr->ConnectOutput(pid, 2, ce);

  if(!TSTBIT(map, kPIDRefMaker)) return;

  //AliLog::SetClassDebugLevel("AliTRDpidRefMaker", 3);
  //AliLog::SetClassDebugLevel("AliTRDpidRefMakerNN", 3);
  //AliLog::SetClassDebugLevel("AliTRDpidRefMakerLQ", 3);

  // TRD pid reference maker NN
  AliTRDpidRefMaker *ref(NULL);
  mgr->AddTask(ref = new AliTRDpidRefMakerNN((char*)"refMakerNN"));
  ref->SetDebugLevel(3);
  ref->SetMCdata(mgr->GetMCtruthEventHandler());
  ref->SetFriends(kTRUE);
  mgr->ConnectInput( ref, 0, mgr->GetCommonInputContainer());
  mgr->ConnectInput( ref, 1, ci[0]);
  mgr->ConnectInput( ref, 2, ci[1]);
  mgr->ConnectInput( ref, 3, ce);
  mgr->ConnectOutput(ref, 1, mgr->CreateContainer("MonitorNN", TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD.CalibPIDrefMaker",mgr->GetCommonFileName())));
  mgr->ConnectOutput(ref, 2, mgr->CreateContainer(ref->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD.CalibPIDrefMaker", mgr->GetCommonFileName())));

  // TRD pid reference maker LQ 
  mgr->AddTask(ref = new AliTRDpidRefMakerLQ((char*)"refMakerLQ"));
  ref->SetDebugLevel(3);
  ref->SetMCdata(mgr->GetMCtruthEventHandler());
  ref->SetFriends(kTRUE);
  mgr->ConnectInput(ref, 0, mgr->GetCommonInputContainer());
  mgr->ConnectInput(ref, 1, ci[0]);
  mgr->ConnectInput(ref, 2, ci[1]);
  mgr->ConnectInput(ref, 3, ce);
  mgr->ConnectOutput(ref, 1, mgr->CreateContainer("MonitorLQ", TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD.CalibPIDrefMaker", mgr->GetCommonFileName())));
  mgr->ConnectOutput(ref, 2, mgr->CreateContainer(ref->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD.CalibPIDrefMaker", mgr->GetCommonFileName())));
  mgr->ConnectOutput(ref, 3, mgr->CreateContainer("PDF", TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD.CalibPIDrefMakerLQ", mgr->GetCommonFileName())));
}

