#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/macros/AliTRDperformanceTrain.h"
#include "PWG1/TRD/AliTRDcheckPID.h"
#include "PWG1/TRD/AliTRDpidRefMaker.h"
#include "PWG1/TRD/AliTRDpidRefMakerLQ.h"
#endif

#include "PWG1/TRD/macros/helper.C"
void AddTRDcheckPID(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  AliLog::SetClassDebugLevel("AliTRDcheckPID", 5);  
  Int_t map = ParseOptions(trd);
  AliAnalysisDataContainer *ce(NULL);
  if(TSTBIT(map, kCheckPID)){
    AliTRDcheckPID *pid = 0x0;
    mgr->AddTask(pid = new AliTRDcheckPID((char*)"checkPID"));
    pid->SetDebugLevel(0);
    pid->SetMCdata(mgr->GetMCtruthEventHandler());
    mgr->ConnectInput (pid, 0, mgr->GetCommonInputContainer());
    mgr->ConnectInput (pid, 1, ci[0]);
    mgr->ConnectOutput(pid, 1, mgr->CreateContainer(pid->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.Performance.root"));

    // define PID exchange container
    ce = mgr->CreateContainer("InfoPID", TObjArray::Class(), AliAnalysisManager::kExchangeContainer);
    mgr->ConnectOutput(pid, 2, ce);
  }  

  if(TSTBIT(map, kPIDRefMaker)){
    // TRD pid reference maker 
    AliTRDpidRefMaker *ref = new AliTRDpidRefMaker(); 
    mgr->AddTask(ref);
    ref->SetDebugLevel(3);
    AliLog::SetClassDebugLevel("AliTRDpidRefMaker", 3);
    ref->SetMCdata(mgr->GetMCtruthEventHandler());
    ref->SetFriends(kTRUE);

    // link basic ref maker
    mgr->ConnectInput( ref, 1, ci[0]);
    mgr->ConnectInput( ref, 2, ci[2]);
    if(ce) mgr->ConnectInput( ref, 2, ce);
    mgr->ConnectOutput(ref, 1, mgr->CreateContainer(Form("Moni%s", ref->GetName()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", ref->GetName())));
    mgr->ConnectOutput(ref, 2, mgr->CreateContainer(ref->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", ref->GetName())));

    // TRD pid reference maker LQ 
    AliTRDpidRefMakerLQ *lq = new AliTRDpidRefMakerLQ(); 
    mgr->AddTask(lq);
    lq->SetDebugLevel(3);
    AliLog::SetClassDebugLevel("AliTRDpidRefMakerLQ", 3);
    lq->SetMCdata(mgr->GetMCtruthEventHandler());
    lq->SetFriends(kTRUE);
    mgr->ConnectInput(lq, 1, ci[0]);
    mgr->ConnectInput(lq, 0, ci[2]);
    if(ce) mgr->ConnectInput(lq, 2, ce);
    mgr->ConnectOutput(lq, 1, mgr->CreateContainer(Form("Moni%s", lq->GetName()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", ref->GetName())));
    mgr->ConnectOutput(lq, 2, mgr->CreateContainer(lq->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", ref->GetName())));
    mgr->ConnectOutput(lq, 3, mgr->CreateContainer("PDF", TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", lq->GetName())));
  }
}

