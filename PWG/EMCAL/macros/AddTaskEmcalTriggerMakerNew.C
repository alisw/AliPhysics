AliEmcalTriggerMakerTask* AddTaskEmcalTriggerMakerNew(
  const char *triggersOutName     = "EmcalTriggers",
  const char *cellsName           = "",
  const char *triggersName        = "",
  bool doQA                       = kFALSE,
  const char* suffix              = ""
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalTriggerMaker", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AddTaskEmcalTriggerMaker", "This task requires an input event handler");
    return 0;
  }

  TString taskName("AliEmcalTriggerMakerTask"), strtriggersOutName(triggersOutName);
  if (!TString(suffix).IsNull()) {
    taskName += TString::Format("_%s", suffix);
    strtriggersOutName += TString::Format("_%s", suffix);
  }

  // Check if the task already exists, if yes only return the pointer
  AliEmcalTriggerMakerTask* eTask = 0;
  if ((eTask = dynamic_cast<AliEmcalTriggerMakerTask *>(mgr->GetTask(taskName)))) return eTask;

  TString strTriggersName(triggersName);
  TString strCellsName(cellsName);

  if(strTriggersName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strTriggersName = "EMCALTrigger";
      ::Info("AddTaskEmcalTriggerMaker", "ESD analysis, triggersName = \"%s\"", strTriggersName.Data());
    }
    else {
      strTriggersName = "emcalTrigger";
      ::Info("AddTaskEmcalTriggerMaker", "AOD analysis, triggersName = \"%s\"", strTriggersName.Data());
    }
  }

  if(strCellsName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strCellsName = "EMCALCells";
      ::Info("AddTaskEmcalTriggerMaker", "ESD analysis, cellsName = \"%s\"", strCellsName.Data());
    }
    else {
      strCellsName = "emcalCells";
      ::Info("AddTaskEmcalTriggerMaker", "AOD analysis, cellsName = \"%s\"", strCellsName.Data());
    }
  }

  TString v0Name;
  if (evhand->InheritsFrom("AliESDInputHandler")) {
    v0Name = "AliESDVZERO";
    ::Info("AddTaskEmcalTriggerMaker", "ESD analysis, v0Name = \"%s\"", v0Name.Data());
  }
  else {
    v0Name = "AliAODVZERO";
    ::Info("AddTaskEmcalTriggerMaker", "AOD analysis, v0Name = \"%s\"", v0Name.Data());
  }
 
   //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  eTask = new AliEmcalTriggerMakerTask(taskName, doQA);
  eTask->SetCaloTriggersName(strTriggersName.Data());
  eTask->SetCaloTriggersOutName(strtriggersOutName.Data());
  eTask->SetCaloCellsName(strCellsName.Data());
  eTask->SetV0InName(v0Name);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );

  if(doQA){
    TString commonoutput = mgr->GetCommonFileName();
    commonoutput += ":TriggerQA";
    mgr->ConnectOutput(eTask, 1, mgr->CreateContainer("TriggerQA", TList::Class(), AliAnalysisManager::kOutputContainer, commonoutput.Data()));
  }

  return eTask;
}
