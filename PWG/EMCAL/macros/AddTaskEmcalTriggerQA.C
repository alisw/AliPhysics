AliEmcalTriggerQATask* AddTaskEmcalTriggerQA(
  const char* triggerPatchesName  = "EmcalTriggers",
  const char* cellsName           = 0,
  const char* triggersName        = 0,
  const char* suffix              = "")
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalTriggerQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AddTaskEmcalTriggerQA", "This task requires an input event handler");
    return NULL;
  }

   //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString taskName("AliEmcalTriggerQATask");
  if (strcmp(suffix,"") != 0) {
    taskName += "_";
    taskName += suffix;
  }
  AliEmcalTriggerQATask* eTask = new AliEmcalTriggerQATask(taskName);
  eTask->SetTriggerPatchesName(triggerPatchesName);
  eTask->Set2015CaloTriggerNames();
  eTask->GetTriggerQA()->EnablePatchType(AliEMCALTriggerQA::kOnlinePatch);
  eTask->GetTriggerQA()->EnablePatchType(AliEMCALTriggerQA::kOfflinePatch);
  eTask->GetTriggerQA()->EnablePatchType(AliEMCALTriggerQA::kRecalcPatch);
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

  eTask->SetCaloTriggersName(strTriggersName);
  eTask->SetCaloCellsName(strCellsName);
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );

  TString commonoutput = mgr->GetCommonFileName();
  TString contOutName(Form("%s_histos", taskName.Data()));
  mgr->ConnectOutput(eTask, 1, mgr->CreateContainer(contOutName, TList::Class(), AliAnalysisManager::kOutputContainer, commonoutput.Data()));

  return eTask;
}
