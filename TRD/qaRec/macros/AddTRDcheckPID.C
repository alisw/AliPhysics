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
  if(TSTBIT(map, kCheckPID)){
    AliTRDcheckPID *pid = 0x0;
    mgr->AddTask(pid = new AliTRDcheckPID());
    pid->SetDebugLevel(0);
    pid->SetMCdata(mgr->GetMCtruthEventHandler());
    mgr->ConnectInput(pid, 0, ci[0]);
    mgr->ConnectOutput(pid, 0, mgr->CreateContainer(pid->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.Performance.root"));
  }

  if(TSTBIT(map, kPIDRefMakerNN)){
    // TRD NN pid reference maker 
    AliTRDpidRefMakerNN *ref = new AliTRDpidRefMakerNN(); 
    //AliLog::SetClassDebugLevel("AliTRDpidRefMakerNN", 4);
    // Neural network PID
    mgr->AddTask(ref);
    ref->SetDebugLevel(3);
    AliLog::SetClassDebugLevel("AliTRDpidRefMakerNN", 3);
    ref->SetMCdata(mgr->GetMCtruthEventHandler());
    mgr->ConnectInput( ref, 0, ci[0]);
    mgr->ConnectInput( ref, 1, ci[2]);
    mgr->ConnectOutput(ref, 0, mgr->CreateContainer(Form("Moni%s", ref->GetName()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", ref->GetName())));
    mgr->ConnectOutput(ref, 1, mgr->CreateContainer(ref->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", ref->GetName())));
  }

  // Multidimensional Likelihood PID
  if(TSTBIT(map, kPIDRefMakerLQ)){
    AliTRDpidRefMakerLQ *reflq = new AliTRDpidRefMakerLQ(); 
    mgr->AddTask(reflq);
    reflq->SetDebugLevel(3);
    //AliLog::SetClassDebugLevel("AliTRDpidRefMakerLQ", 3);
    reflq->SetMCdata(mgr->GetMCtruthEventHandler());
    mgr->ConnectInput( reflq, 0, ci[0]);
    mgr->ConnectInput( reflq, 1, ci[2]);
    mgr->ConnectOutput(reflq, 0, mgr->CreateContainer(Form("Moni%s", reflq->GetName()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", reflq->GetName())));
    mgr->ConnectOutput(reflq, 1, mgr->CreateContainer(reflq->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", reflq->GetName())));
  }
}

