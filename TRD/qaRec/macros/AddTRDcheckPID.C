#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "TRD/qaRec/macros/AliTRDperformanceTrain.h"
#include "TRD/qaRec/AliTRDcheckPID.h"
#include "TRD/qaRec/AliTRDpidRefMaker.h"
#endif


void AddTRDcheckPID(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  Int_t map = ParseOptions(trd);
  if(!(TSTBIT(map, kCheckPID))) return;

  AliTRDcheckPID *pid = 0x0;
  mgr->AddTask(pid = new AliTRDcheckPID());
  pid->SetDebugLevel(0);
  pid->SetMCdata(mgr->GetMCtruthEventHandler());
  mgr->ConnectInput(pid, 0, ci[0]);
  mgr->ConnectOutput(pid, 0, mgr->CreateContainer(pid->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.Performance.root"));

  // TRD pid reference 
  if(!(TSTBIT(map, kPIDRefMaker))) return;
  AliTRDpidRefMaker *ref = 0x0; 
  mgr->AddTask(ref = new AliTRDpidRefMaker());
  ref->SetDebugLevel(0);
  ref->SetMCdata(mgr->GetMCtruthEventHandler());
  
  // Create containers for input/output
  mgr->ConnectInput( ref, 0, ci[0]);
  mgr->ConnectInput( ref, 1, ci[2]);
  mgr->ConnectOutput(ref, 0, mgr->CreateContainer(ref->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", ref->GetName())));
  
  // network container
  AliAnalysisDataContainer *co[] = {0x0, 0x0};
  co[0] = mgr->CreateContainer(Form("%sNN", ref->GetName()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%sNN.root", ref->GetName()));
  // likelihood container
  co[1] = mgr->CreateContainer(Form("%sLQ", ref->GetName()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%sLQ.root", ref->GetName()));
  mgr->ConnectOutput(ref, 1, co[0]);
  mgr->ConnectOutput(ref, 2, co[1]);
}

