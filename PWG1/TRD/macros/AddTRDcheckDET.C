#if ! defined (__CINT__) || defined (__MAKECINT__)
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "PWG1/TRD/macros/AliTRDperformanceTrain.h"
#include "PWG1/TRD/AliTRDcheckDET.h"
#include "PWG1/TRD/AliTRDcalibration.h"
#endif

#include "PWG1/TRD/macros/helper.C"
void AddTRDcheckDET(AliAnalysisManager *mgr, Char_t *trd, AliAnalysisDataContainer **ci/*, AliAnalysisDataContainer **co*/)
{
  Int_t map = ParseOptions(trd);
  if(!(TSTBIT(map, kCheckDET))) return;
  printf("AddTRDcheckDET <- [0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\"\n", ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName(), ci[3]->GetName());

  //AliLog::SetClassDebugLevel("AliTRDcheckDET", 5);
  AliTRDcheckDET *task(NULL);
  mgr->AddTask(task = new AliTRDcheckDET((char*)"checkDET"));
  task->UseClustersOutsideChamber();
  task->SetDebugLevel(0);
  task->SetMCdata(mgr->GetMCtruthEventHandler());
  
  // Create containers for input/output
  mgr->ConnectInput ( task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectInput ( task, 1, ci[1]);
  mgr->ConnectInput ( task, 2, ci[0]);
  mgr->ConnectOutput( task, 1, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
  

  // CALIBRATION
  if(!(TSTBIT(map, kCalibration))) return;
  AliTRDcalibration *ctask = 0x0;
  mgr->AddTask(ctask = new AliTRDcalibration((char*)"calibration"));
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
  mgr->ConnectInput(ctask,  0, mgr->GetCommonInputContainer());
  mgr->ConnectInput(ctask,  1, ci[1]);
  mgr->ConnectOutput(ctask, 1, mgr->CreateContainer(ctask->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD.Calib%s", mgr->GetCommonFileName(),ctask->GetName())));
}
