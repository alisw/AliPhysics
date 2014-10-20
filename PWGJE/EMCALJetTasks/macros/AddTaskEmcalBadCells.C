AliAnalysisTaskEmcalBadCells *AddTaskEmcalBadCells(
						   const char *CentEst             = "V0A",
						   Int_t       pSel                = AliVEvent::kINT7,
						   TString     kEmcalCellsName     = "emcalCells",
						   TString     tag                 = ""
) {

  // #### Define manager and data container names
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEmcalBadCells", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalBadCells", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName = Form("EmcalBadCells%s",tag.Data());

  //Configure DiJet task
  AliAnalysisTaskEmcalBadCells *task = new AliAnalysisTaskEmcalBadCells(wagonName.Data(),kTRUE);

  task->SetCaloCellsName(kEmcalCellsName.Data());
  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  AliAnalysisDataContainer *coutput1 = 0x0;

  TString containerName1 = Form("%s",wagonName.Data());

  TString outputfile = Form("%s:%s",AliAnalysisManager::GetCommonFileName(),wagonName.Data());

  coutput1 = mgr->CreateContainer(containerName1, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  mgr->ConnectOutput(task,1,coutput1);
  
  return task;

}
