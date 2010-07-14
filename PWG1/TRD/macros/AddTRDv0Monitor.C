#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "TTree.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/macros/AliTRDperformanceTrain.h"
#include "PWG1/TRD/AliTRDv0Monitor.h"
#endif

#include "PWG1/TRD/macros/helper.C"
void AddTRDv0Monitor(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  Int_t map = ParseOptions(trd);
  if(!TSTBIT(map, kV0Monitor)) return;
  printf("AddTRDv0Monitor <- [0]=\"%s\" [1]=\"%s\" [2]=\"%s\"\n", ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName());

  AliTRDv0Monitor *v0Mon = new AliTRDv0Monitor("v0Monitor", "v0Monitor");
  mgr->AddTask(v0Mon);
  v0Mon->SetDebugLevel(0);
  //AliLog::SetClassDebugLevel("AliTRDpidRefMaker", 3);
  v0Mon->SetMCdata(mgr->GetMCtruthEventHandler());
  v0Mon->SetFriends(kTRUE);
  //v0Mon->SetSource(AliTRDpidRefMaker::kV0,AliTRDpidRefMaker::kRec);
  mgr->ConnectInput( v0Mon, 0, mgr->GetCommonInputContainer());
  mgr->ConnectInput( v0Mon, 1, ci[0]);
  mgr->ConnectInput( v0Mon, 2, ci[1]);
  mgr->ConnectInput( v0Mon, 3, ci[2]);

  mgr->ConnectOutput(v0Mon, 1, mgr->CreateContainer(v0Mon->GetName(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
  //mgr->ConnectOutput(v0Mon, 2, mgr->CreateContainer(v0Mon->GetName(), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD.CalibPIDrefMaker", mgr->GetCommonFileName())));
}

