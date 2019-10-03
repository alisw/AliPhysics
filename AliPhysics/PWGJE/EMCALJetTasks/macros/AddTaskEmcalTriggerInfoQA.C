// $Id$

AliAnalysisTaskEmcalTriggerInfoQA* AddTaskEmcalTriggerInfoQA(
  const char *triggersName     = "EmcalTriggers",
  const char *triggerSetupName = "EmcalTriggerSetup",
  const char *cellsName           = 0,
  const char *taskName            = "AliEmcalTriggerInfoQA"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalTriggerInfoQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AddTaskEmcalTriggerInfoQA", "This task requires an input event handler");
    return NULL;
  }

  if (!cellsName) {
    cellsName = new char[100];

    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strcpy(cellsName,"EMCALCells");
      ::Info("AddTaskEmcalTriggerInfoQA", Form( "ESD analysis, cellsName = \"%s\"", cellsName ));
    }
    else {
      strcpy(cellsName,"emcalCells");
      ::Info("AddTaskEmcalTriggerInfoQA", Form( "AOD analysis, cellsName = \"%s\"", cellsName ));
    }
  }
 
   //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliAnalysisTaskEmcalTriggerInfoQA *eTask = new AliAnalysisTaskEmcalTriggerInfoQA(taskName);
  eTask->SetCaloTriggerPatchInfoName(triggersName);
  eTask->SetCaloTriggerSetupInfoName(triggerSetupName);
  eTask->SetCaloCellsName(cellsName);
  eTask->SetAnaType(AliAnalysisTaskEmcal::kEMCAL);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
	TString listName = Form("ListEmcalTriggerInfoQA" );
	TString fileName = Form("%s:EmcalTriggerInfoQA", AliAnalysisManager::GetCommonFileName());

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );
  
 	AliAnalysisDataContainer *coutput = mgr->CreateContainer( listName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
	mgr->ConnectOutput( eTask, 1, coutput);

 return eTask;
}
