#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "TError.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWGPP/TRD/AliTRDv0Monitor.h"
#endif

void AddTRDv0Monitor(AliAnalysisManager *mgr, Int_t /*map*/, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  Info("AddTRDv0Monitor", "[0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\"", 
       ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName(), ci[3]->GetName());

  AliTRDv0Monitor *v0Mon(NULL);;
  mgr->AddTask(v0Mon = new AliTRDv0Monitor((char*)"TRDv0Monitor"));
  v0Mon->SetDebugLevel(0);
  //AliLog::SetClassDebugLevel("AliTRDpidRefMaker", 3);
  v0Mon->SetMCdata(mgr->GetMCtruthEventHandler());
  v0Mon->SetFriends(kTRUE);
  //v0Mon->SetSource(AliTRDpidRefMaker::kV0,AliTRDpidRefMaker::kRec);
  mgr->ConnectInput( v0Mon, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container
  mgr->ConnectInput( v0Mon, 1, ci[0]);                          // connect barrel tracks container
  mgr->ConnectInput( v0Mon, 2, ci[1]);                          // connect event info container
  mgr->ConnectInput( v0Mon, 3, ci[2]);                          // connect V0s container
  mgr->ConnectInput( v0Mon, 4, ci[3]);                          // connect pid Info container

  mgr->ConnectOutput(v0Mon, 1, mgr->CreateContainer(v0Mon->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
  //mgr->ConnectOutput(v0Mon, 2, mgr->CreateContainer(v0Mon->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD.CalibPIDrefMaker", mgr->GetCommonFileName())));
}

