// $Id$

AliEmcalTriggerMaker* AddTaskEmcalTriggerMaker(
  const char *triggersOutName     = "EmcalTriggers",
  const char *triggerSetupOutName = "EmcalTriggerSetup",
  const char *cellsName           = 0,
  const char *triggersName        = 0,
  const char *taskName            = "AliEmcalTriggerMaker",
  int jetLowA                     = 0,
  int jetLowB                     = 0,
  int jetLowC                     = 0,
  int jetHighA                    = 0,
  int jetHighB                    = 0,
  int jetHighC                    = 0
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalTriggerMaker", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AddTaskEmcalTriggerMaker", "This task requires an input event handler");
    return NULL;
  }

  if (!triggersName) {
    triggersName = new char[100];

    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strcpy(triggersName,"EMCALTrigger");
      ::Info("AddTaskEmcalTriggerMaker", Form( "ESD analysis, triggersName = \"%s\"", triggersName ));
    }
    else {
      strcpy(triggersName,"emcalTrigger");
      ::Info("AddTaskEmcalTriggerMaker", Form( "AOD analysis, triggersName = \"%s\"", triggersName ));
    }
  }
  if (!cellsName) {
    cellsName = new char[100];

    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strcpy(cellsName,"EMCALCells");
      ::Info("AddTaskEmcalTriggerMaker", Form( "ESD analysis, cellsName = \"%s\"", cellsName ));
    }
    else {
      strcpy(cellsName,"emcalCells");
      ::Info("AddTaskEmcalTriggerMaker", Form( "AOD analysis, cellsName = \"%s\"", cellsName ));
    }
  }
  char *v0Name;
  v0Name = new char[100];
  if (evhand->InheritsFrom("AliESDInputHandler")) {
    strcpy(v0Name,"AliESDVZERO");
    ::Info("AddTaskEmcalTriggerMaker", Form( "ESD analysis, v0Name = \"%s\"", v0Name ));
  }
  else {
    strcpy(v0Name,"AliAODVZERO");
    ::Info("AddTaskEmcalTriggerMaker", Form( "AOD analysis, v0Name = \"%s\"", cellsName ));
  }
 
   //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalTriggerMaker *eTask = new AliEmcalTriggerMaker(taskName);
  eTask->SetCaloTriggersName(triggersName);
  eTask->SetCaloTriggersOutName(triggersOutName);
  eTask->SetCaloTriggerSetupOutName(triggerSetupOutName);
  eTask->SetCaloCellsName(cellsName);
  eTask->SetV0InName(v0Name);
  eTask->SetTriggerThresholdJetLow( jetLowA, jetLowB, jetLowC );
  eTask->SetTriggerThresholdJetHigh( jetHighA, jetHighB, jetHighC );

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );
  
  return eTask;
}
