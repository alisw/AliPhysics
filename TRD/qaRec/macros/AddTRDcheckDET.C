#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "TRD/qaRec/macros/AliTRDperformanceTrain.h"
#include "TRD/qaRec/AliTRDcheckDET.h"
#include "TRD/qaRec/AliTRDcalibration.h"
#endif

#include "TRD/qaRec/macros/helper.C"
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
  ctask->SetHisto2d(kTRUE);
  ctask->SetVector2d(kTRUE);
  ctask->SetVdriftLinear(kTRUE);
  ctask->SetNz(0,0);
  ctask->SetNrphi(0,0);
  ctask->SetNz(0,1);
  ctask->SetNrphi(0,1);
  ctask->SetNz(0,2);
  ctask->SetNrphi(0,2);
  ctask->SetLow(0);
  ctask->SetHigh(30);
  ctask->SetFillZero(kFALSE);
  ctask->SetDebugLevel(1);

  // Create containers for input/output
  mgr->ConnectInput(ctask, 0, ci[0]);
  mgr->ConnectOutput(ctask, 0, mgr->CreateContainer(ctask->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("TRD.Task%s.root", ctask->GetName())));
}
