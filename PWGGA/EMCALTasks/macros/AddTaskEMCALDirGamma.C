AliAnalysisTask *AddTaskEMCALDirGamma(const UInt_t triggermask = AliVEvent::kMB,
                                      Bool_t mcmode = 0,
                                      const char geoname[] = "EMCAL_COMPLETEV1",
                                      Bool_t qf = 0,
                                      Int_t ncells = 2)
{
  //gSystem->Load("libPWGGAGammaConv.so");

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALDirGamma", "No analysis manager to connect to.");
    return NULL;
  }  

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEMCALDirGamma", "This task requires an input event handler");
    return NULL;
  }

//  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
//  Bool_t isMC=kFALSE; // kTRUE in case of MC
//  AddTaskPIDResponse(isMC);
  
  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  //connect input V0Reader
//  mgr->ConnectInput(fV0ReaderV1,0,cinput);
  
  
  // Create the task and configure it.
  //===========================================================================
  AliAnalysisTaskEMCALDirGamma* task = new  AliAnalysisTaskEMCALDirGamma("EMCALDirGammaTask");
  task->SelectCollisionCandidates(triggermask);
  task->SetUseQualFlag(qf);
  task->SetAsymMax1(0.3);
  task->SetAsymMax2(0.7);
  task->SetMinEcc(-10);
  task->SetNminCells(ncells);
  task->SetGeoName(geoname);
  task->SetMcMode(mcmode);
  task->SetAddedSignal(0);
  task->SetFillNtuple(0);
  mgr->AddTask(task);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosEMCALDirGamma",
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (task, 1, coutput1 );
  
  return task;
  
}
