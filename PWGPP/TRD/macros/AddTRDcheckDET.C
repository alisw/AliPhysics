void AddTRDcheckDET(AliAnalysisDataContainer **ci, Int_t cal)
{
  Info("AddTRDcheckDET", "[0]=\"%s\" [1]=\"%s\" [2]=\"%s\" [3]=\"%s\" [4]=\"%s\"", ci[0]->GetName(), ci[1]->GetName(), ci[2]->GetName(), ci[3]->GetName(), ci[4]->GetName());
  AliAnalysisDataContainer *evInfoContainer = ci[3];

  //AliLog::SetClassDebugLevel("AliTRDcheckDET", 5);
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) return;
  AliTRDcheckDET *task(NULL);
  mgr->AddTask(task = new AliTRDcheckDET((char*)"TRDcheckDET"));
  task->UseClustersOutsideChamber();
  task->SetDebugLevel(0);
  task->SetMCdata((Bool_t)mgr->GetMCtruthEventHandler());
  
  // Create containers for input/output
  Int_t trackStatus = 0; // barrel tracks
//                    = 1; // kink tracks
//                    = 2; // SA tracks
  mgr->ConnectInput ( task, 0, mgr->GetCommonInputContainer()); // connect main (ESD) container
  mgr->ConnectInput ( task, 1, ci[trackStatus]);                // conect track info container
  mgr->ConnectInput ( task, 2, evInfoContainer);                // conect event info container
  mgr->ConnectInput ( task, 3, ci[4]);                          // conect onl.tracklets container
  mgr->ConnectInput ( task, 4, ci[5]);                          // conect clusters container
  mgr->ConnectOutput( task, 1, mgr->CreateContainer(task->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Performance",mgr->GetCommonFileName())));
  

  // CALIBRATION
  if(!cal) return;
  AliTRDcalibration *ctask(NULL);
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
  mgr->ConnectInput(ctask,  1, ci[0]);
  mgr->ConnectOutput(ctask, 1, mgr->CreateContainer(ctask->GetName(), TObjArray::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TRD_Calibration", mgr->GetCommonFileName())));
}
