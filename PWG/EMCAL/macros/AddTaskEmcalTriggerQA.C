/// \file AddTaskEmcalTriggerQA.C
/// \brief This is the AddTask macro of AliEmcalTriggerQATask
///
/// This macro is used in a LEGO train to add an instance
/// of AliEmcalTriggerQATask in the analysis manager.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Apr 4, 2016

AliEmcalTriggerQATask* AddTaskEmcalTriggerQA(
    const char* triggerPatchesName  = "EmcalTriggers",
    const char* cellsName           = 0,
    const char* triggersName        = 0,
    Int_t       nCentBins           = 0,
    Bool_t      online              = kFALSE,
    const char* suffix              = "")
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AliEmcalTriggerQATask", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AliEmcalTriggerQATask", "This task requires an input event handler");
    return 0;
  }

   //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString taskName("AliEmcalTriggerQATask");
  if (strcmp(suffix,"") != 0) {
    taskName += "_";
    taskName += suffix;
  }
  AliEmcalTriggerQATask* eTask = new AliEmcalTriggerQATask(taskName, nCentBins, online);
  eTask->SetTriggerPatchesName(triggerPatchesName);
  eTask->GetTriggerQA()->EnablePatchType(AliEMCALTriggerQA::kOnlinePatch);
  eTask->GetTriggerQA()->EnablePatchType(AliEMCALTriggerQA::kOfflinePatch);
  eTask->GetTriggerQA()->EnablePatchType(AliEMCALTriggerQA::kRecalcPatch);
  TString strTriggersName(triggersName);
  TString strCellsName(cellsName);
  if(strTriggersName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strTriggersName = "EMCALTrigger";
      ::Info("AddTaskEmcalTriggerQA", "ESD analysis, triggersName = \"%s\"", strTriggersName.Data());
    }
    else {
      strTriggersName = "emcalTrigger";
      ::Info("AddTaskEmcalTriggerQA", "AOD analysis, triggersName = \"%s\"", strTriggersName.Data());
    }
  }
  if(strCellsName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strCellsName = "EMCALCells";
      ::Info("AddTaskEmcalTriggerQA", "ESD analysis, cellsName = \"%s\"", strCellsName.Data());
    }
    else {
      strCellsName = "emcalCells";
      ::Info("AddTaskEmcalTriggerQA", "AOD analysis, cellsName = \"%s\"", strCellsName.Data());
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
