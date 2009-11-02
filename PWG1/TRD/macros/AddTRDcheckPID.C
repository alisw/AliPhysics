#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/macros/AliTRDperformanceTrain.h"
#include "PWG1/TRD/AliTRDcheckPID.h"
#include "PWG1/TRD/AliTRDpidRefMaker.h"
#endif

#include "PWG1/TRD/macros/helper.C"
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

  if(TSTBIT(map, kPIDRefMaker)){
    // TRD pid reference maker 
    AliTRDpidRefMaker *ref = new AliTRDpidRefMaker(); 
    mgr->AddTask(ref);
    ref->SetDebugLevel(3);
    AliLog::SetClassDebugLevel("AliTRDpidRefMaker", 3);
    ref->SetMCdata(mgr->GetMCtruthEventHandler());
    ref->SetFriends(kTRUE);
    mgr->ConnectInput( ref, 0, ci[0]);
    mgr->ConnectInput( ref, 1, ci[2]);
    mgr->ConnectOutput(ref, 0, mgr->CreateContainer(Form("Moni%s", ref->GetName()), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", ref->GetName())));
    mgr->ConnectOutput(ref, 1, mgr->CreateContainer(ref->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Calib%s.root", ref->GetName())));
  }
}

