#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "TRD/qaRec/macros/AliTRDperformanceTrain.h"
#include "TRD/qaRec/AliTRDcheckDET.h"
#include "TRD/qaRec/AliTRDcalibration.h"
#endif

void AddTRDcheckDET(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  Int_t map = ParseOptions(trd);
  if(!(TSTBIT(map, kCheckDET))) return;

  AliTRDcheckDET *task = 0x0;
  mgr->AddTask(task = new AliTRDcheckDET());
  task->SetDebugLevel(0);
  task->SetMCdata(mgr->GetMCtruthEventHandler());
  
  // Create containers for input/output
  mgr->ConnectInput( task, 0, ci[0]);
  mgr->ConnectInput( task, 1, ci[1]);
  mgr->ConnectOutput(task, 0, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, "TRD.Performance.root"));
  

  // CALIBRATION
  if(!(TSTBIT(map, kCalibration))) return;
  AliTRDcalibration *ctask = 0x0;
  mgr->AddTask(ctask = new AliTRDcalibration());
  ctask->SetLow(0);
  ctask->SetHigh(30);
  ctask->SetFillZero(kFALSE);
  ctask->SetDebugLevel(0);

  // Create containers for input/output
  mgr->ConnectInput(task, 0, ci[0]);
  mgr->ConnectOutput(task, 0, mgr->CreateContainer(ctask->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", ctask->GetName())));
}

