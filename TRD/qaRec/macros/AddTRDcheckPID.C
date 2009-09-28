#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "TRD/qaRec/macros/AliTRDperformanceTrain.h"
#include "TRD/qaRec/AliTRDcheckPID.h"
#include "TRD/qaRec/AliTRDpidRefMakerNN.h"
#include "TRD/qaRec/AliTRDpidRefMakerLQ.h"
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

  if(!(TSTBIT(map, kPIDRefMaker))) return;
  // TRD pid reference 
  AliTRDpidRefMaker *ref = 0x0; 
  // Neural network PID
  mgr->AddTask(ref = new AliTRDpidRefMakerNN());
  ref->SetDebugLevel(3);
  AliLog::SetClassDebugLevel("AliTRDpidRefMakerNN", 3);
  ref->SetMCdata(mgr->GetMCtruthEventHandler());
  mgr->ConnectInput( ref, 0, ci[0]);
  mgr->ConnectInput( ref, 1, ci[2]);
  mgr->ConnectOutput(ref, 0, mgr->CreateContainer(ref->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", ref->GetName())));
  mgr->ConnectOutput(ref, 1, mgr->CreateContainer(Form("%sNN", ref->GetName()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.%s.root", ref->GetName())));

  // Multidimensional Likelihood PID
  mgr->AddTask(ref = new AliTRDpidRefMakerLQ());
  ref->SetDebugLevel(3);
  AliLog::SetClassDebugLevel("AliTRDpidRefMakerLQ", 3);
  ref->SetMCdata(mgr->GetMCtruthEventHandler());
  mgr->ConnectInput( ref, 0, ci[0]);
  mgr->ConnectInput( ref, 1, ci[2]);
  mgr->ConnectOutput(ref, 0, mgr->CreateContainer(ref->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", ref->GetName())));
  mgr->ConnectOutput(ref, 1, mgr->CreateContainer(Form("%sLQ", ref->GetName()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.%s.root", ref->GetName())));
}

