// $Id$

AliEmcalTriggerMaker* AddTaskEmcalTriggerMaker(const char *configurationString)
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

  // Handle wagon configuration
  AliEMCALConfiguration defaultConfiguration("triggerMakerDefault");
  defaultConfiguration.AddParam("triggersOutName", new AliEMCALConfigurationValueString("EmcalTriggers"));
  defaultConfiguration.AddParam("cellsName", new AliEMCALConfigurationValueString(""));
  defaultConfiguration.AddParam("triggersName", new AliEMCALConfigurationValueString(""));
  defaultConfiguration.AddParam("taskName", new AliEMCALConfigurationValueString("AliEmcalTriggerMaker"));
  defaultConfiguration.AddParam("jetLowA", new AliEMCALConfigurationValueInt(0));
  defaultConfiguration.AddParam("jetLowB", new AliEMCALConfigurationValueInt(0));
  defaultConfiguration.AddParam("jetLowC", new AliEMCALConfigurationValueInt(0));
  defaultConfiguration.AddParam("jetHighA", new AliEMCALConfigurationValueInt(0));
  defaultConfiguration.AddParam("jetHighB", new AliEMCALConfigurationValueInt(0));
  defaultConfiguration.AddParam("jetHighC", new AliEMCALConfigurationValueInt(0));
  defaultConfiguration.AddParam("doQA", new AliEMCALConfigurationValueBool(kFALSE));
  AliEMCALConfiguration userConfiguration("userConfig");
  userCofiguration.Build(configurationString)
  AliEMCALConfigurationMatcher combinedConfiguration(userConfig, defaultConfig);

  AliEMCALConfigurationValueString *strTriggersName = static_cast<AliEMCALConfigurationValueString *>(combinedConfiguration->GetValue("triggersName")),
  				   *strCellsName = static_cast<AliEMCALConfigurationValueString *>(combinedConfiguration->GetValue("cellsName"));

  if(strTriggersName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strTriggersName = "EMCALTrigger";
      ::Info("AddTaskEmcalTriggerMaker", Form( "ESD analysis, triggersName = \"%s\"", strTriggersName->GetValue() ));
    }
    else {
      strTriggersName = "emcalTrigger";
      ::Info("AddTaskEmcalTriggerMaker", Form( "AOD analysis, triggersName = \"%s\"", strTriggersName->GetValue() ));
    }
  }

  if(strCellsName.IsNull()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strCellsName = "EMCALCells";
      ::Info("AddTaskEmcalTriggerMaker", Form( "ESD analysis, cellsName = \"%s\"", strCellsName->GetValue() ));
    }
    else {
      strCellsName = "emcalCells";
      ::Info("AddTaskEmcalTriggerMaker", Form( "AOD analysis, cellsName = \"%s\"", strCellsName->GetValue() ));
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
    ::Info("AddTaskEmcalTriggerMaker", Form( "AOD analysis, v0Name = \"%s\"", v0Name ));
  }
 
   //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalTriggerMaker *eTask = new AliEmcalTriggerMaker(
			(static_cast<AliEMCALConfigurationValueString *>(combinedConfiguration->GetValue("taskName")))->GetValue(), 
			(static_cast<AliEMCALConfigurationValueBool *>(combinedConfiguration->GetValue("doQA")))->GetValue());
  eTask->SetCaloTriggersName(strTriggersName->GetValue());
  eTask->SetCaloTriggersOutName((static_cast<AliEMCALConfigurationValueString *>(combinedConfiguration->GetValue("triggersOutName"))->GetValue());
  eTask->SetCaloTriggerSetupOutName((static_cast<AliEMCALConfigurationValueString *>(combinedConfiguration->GetValue("triggerSetupOutName")))->GetValue());
  eTask->SetCaloCellsName(strCellsName->GetValue());
  eTask->SetV0InName(v0Name);
  AliEMCALConfigurationValueInt *jetLowA = static_cast<AliEMCALConfigurationValueInt *>(combinedConfiguration->GetValue("jetLowA")),
  				*jetLowB = static_cast<AliEMCALConfigurationValueInt *>(combinedConfiguration->GetValue("jetLowB")),
  				*jetLowC = static_cast<AliEMCALConfigurationValueInt *>(combinedConfiguration->GetValue("jetLowC")),
  				*jetHighA = static_cast<AliEMCALConfigurationValueInt *>(combinedConfiguration->GetValue("jetHighA")),
  				*jetHighB = static_cast<AliEMCALConfigurationValueInt *>(combinedConfiguration->GetValue("jetHighB")),
  				*jetHighC = static_cast<AliEMCALConfigurationValueInt *>(combinedConfiguration->GetValue("jetHighC")),
  eTask->SetTriggerThresholdJetLow( jetLowA->GetValue(), jetLowB->GetValue(), jetLowC->GetValue() );
  eTask->SetTriggerThresholdJetHigh( jetHighA->GetValue(), jetHighB->GetValue(), jetHighC->GetValue() );

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
