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
  int jetHighC                    = 0,
  int gammaLowA                  = 0,
  int gammaLowB                  = 0,
  int gammaLowC                  = 0,
  int gammaHighA                  = 0,
  int gammaHighB                  = 0,
  int gammaHighC                  = 0,
  bool useOldBitConfig            = kFALSE,
  bool doQA                       = kFALSE
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

  // Check if the task already exists, if yes only return the pointer
  AliEmcalTriggerMaker *eTask(NULL);
  if((eTask = dynamic_cast<AliEmcalTriggerMaker *>(mgr->GetTask(taskName)))) return eTask;

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

  char *v0Name;
  v0Name = new char[100];
  if (evhand->InheritsFrom("AliESDInputHandler")) {
    strcpy(v0Name,"AliESDVZERO");
    ::Info("AddTaskEmcalTriggerMaker", "ESD analysis, v0Name = \"%s\"", v0Name);
  }
  else {
    strcpy(v0Name,"AliAODVZERO");
    ::Info("AddTaskEmcalTriggerMaker", "AOD analysis, v0Name = \"%s\"", v0Name);
  }
 
   //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  eTask = new AliEmcalTriggerMaker(taskName, doQA);
  eTask->SetCaloTriggersName(strTriggersName.Data());
  eTask->SetCaloTriggersOutName(triggersOutName);
  eTask->SetCaloTriggerSetupOutName(triggerSetupOutName);
  eTask->SetCaloCellsName(strCellsName.Data());
  eTask->SetV0InName(v0Name);
  eTask->SetTriggerThresholdJetLow( jetLowA, jetLowB, jetLowC );
  eTask->SetTriggerThresholdJetHigh( jetHighA, jetHighB, jetHighC );
  eTask->SetTriggerThresholdGammaLow( gammaLowA, gammaLowB, gammaLowC );
  eTask->SetTriggerThresholdGammaHigh( gammaHighA, gammaHighB, gammaHighC );
  if (useOldBitConfig)
    eTask->SetUseTriggerBitConfig(AliEmcalTriggerMaker::kOldConfig);

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
