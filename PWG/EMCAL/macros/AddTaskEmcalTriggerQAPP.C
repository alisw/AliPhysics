/// \file AddTaskEmcalTriggerQAPP.C
/// \brief This is the AddTask macro of AliEmcalTriggerQATaskPP
///
/// This macro is used in a LEGO train to add an instance
/// of AliEmcalTriggerQATaskPP in the analysis manager.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Feb 23, 2016

AliEmcalTriggerQATaskPP* AddTaskEmcalTriggerQAPP(
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
    ::Error("AliEmcalTriggerQATaskPP", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    ::Error("AliEmcalTriggerQATaskPP", "This task requires an input event handler");
    return 0;
  }

   //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString taskName("AliEmcalTriggerQATaskPP");
  if (strcmp(suffix,"") != 0) {
    taskName += "_";
    taskName += suffix;
  }
  AliEmcalTriggerQATaskPP* eTask = new AliEmcalTriggerQATaskPP(taskName);
  eTask->SetTriggerPatchesName(triggerPatchesName);
  eTask->GetTriggerQA()->EnablePatchType(AliEmcalTriggerQAPP::kOnlinePatch);
  eTask->GetTriggerQA()->EnablePatchType(AliEmcalTriggerQAPP::kOfflinePatch);
  eTask->GetTriggerQA()->EnablePatchType(AliEmcalTriggerQAPP::kRecalcPatch);
  TString strTriggersName(triggersName);
  TString strCellsName(cellsName);

  if(strTriggersName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strTriggersName = "EMCALTrigger";
      ::Info("AliEmcalTriggerQATaskPP", "ESD analysis, triggersName = \"%s\"", strTriggersName.Data());
    }
    else {
      strTriggersName = "emcalTrigger";
      ::Info("AliEmcalTriggerQATaskPP", "AOD analysis, triggersName = \"%s\"", strTriggersName.Data());
    }
  }

  if(strCellsName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strCellsName = "EMCALCells";
      ::Info("AliEmcalTriggerQATaskPP", "ESD analysis, cellsName = \"%s\"", strCellsName.Data());
    }
    else {
      strCellsName = "emcalCells";
      ::Info("AliEmcalTriggerQATaskPP", "AOD analysis, cellsName = \"%s\"", strCellsName.Data());
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
