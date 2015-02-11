// $Id$

AliEmcalTriggerMaker* AddTaskEmcalTriggerMakerJSON(const char *configurationString)
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

  // Handle wagon configuration
  // Definition of possible parameters
  AliEMCALConfiguration defaultConfiguration("triggerMakerDefault");
  defaultConfiguration.AddParam("triggersOutName", new AliJSONString("EmcalTriggers"));
  defaultConfiguration.AddParam("triggerSetupOutName", new AliJSONString("EmcalTriggerSetup"));
  defaultConfiguration.AddParam("cellsName", new AliJSONString(""));
  defaultConfiguration.AddParam("triggersName", new AliJSONString(""));
  defaultConfiguration.AddParam("taskName", new AliJSONString("AliEmcalTriggerMaker"));
  defaultConfiguration.AddParam("jetLowA", new AliJSONInt(0));
  defaultConfiguration.AddParam("jetLowB", new AliJSONInt(0));
  defaultConfiguration.AddParam("jetLowC", new AliJSONInt(0));
  defaultConfiguration.AddParam("jetHighA", new AliJSONInt(0));
  defaultConfiguration.AddParam("jetHighB", new AliJSONInt(0));
  defaultConfiguration.AddParam("jetHighC", new AliJSONInt(0));
  defaultConfiguration.AddParam("gammaLowA", new AliJSONInt(0));
  defaultConfiguration.AddParam("gammaLowB", new AliJSONInt(0));
  defaultConfiguration.AddParam("gammaLowC", new AliJSONInt(0));
  defaultConfiguration.AddParam("gammaHighA", new AliJSONInt(0));
  defaultConfiguration.AddParam("gammaHighB", new AliJSONInt(0));
  defaultConfiguration.AddParam("gammaHighC", new AliJSONInt(0));
  defaultConfiguration.AddParam("doQA", new AliJSONBool(kFALSE));
  AliEMCALConfiguration userConfiguration("userConfig");
  userConfiguration.Build(configurationString);
  AliEMCALConfigurationMatcher combinedConfiguration(&userConfiguration, &defaultConfiguration);

  defaultConfiguration.Print();
  userConfiguration.Print();

  TString strTriggersName = (static_cast<AliJSONString *>(combinedConfiguration.GetValue("triggersName")))->GetValue(),
  				   strCellsName = (static_cast<AliJSONString *>(combinedConfiguration.GetValue("cellsName")))->GetValue();

  if(!strTriggersName.Length()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strTriggersName = "EMCALTrigger";
      ::Info("AddTaskEmcalTriggerMaker", Form( "ESD analysis, triggersName = \"%s\"", strTriggersName.Data() ));
    }
    else {
      strTriggersName = "emcalTrigger";
      ::Info("AddTaskEmcalTriggerMaker", Form( "AOD analysis, triggersName = \"%s\"", strTriggersName.Data() ));
    }
  }

  if(!strCellsName.Length()) {
    if (evhand->InheritsFrom("AliESDInputHandler")) {
      strCellsName = "EMCALCells";
      ::Info("AddTaskEmcalTriggerMaker", Form( "ESD analysis, cellsName = \"%s\"", strCellsName.Data() ));
    }
    else {
      strCellsName = "emcalCells";
      ::Info("AddTaskEmcalTriggerMaker", Form( "AOD analysis, cellsName = \"%s\"", strCellsName.Data() ));
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
  AliJSONBool *doQA = static_cast<AliJSONBool *>(combinedConfiguration.GetValue("doQA"));
  TString taskname = (static_cast<AliJSONString *>(combinedConfiguration.GetValue("taskName")))->GetValue();

  eTask = new AliEmcalTriggerMaker(taskname.Data(), doQA->GetValue());
  eTask->SetCaloTriggersName(strTriggersName.Data());
  eTask->SetCaloTriggersOutName((static_cast<AliJSONString *>(combinedConfiguration.GetValue("triggersOutName")))->GetValue());
  eTask->SetCaloTriggerSetupOutName((static_cast<AliJSONString *>(combinedConfiguration.GetValue("triggerSetupOutName")))->GetValue());
  eTask->SetCaloCellsName(strCellsName.Data());
  eTask->SetV0InName(v0Name);
  AliJSONInt *jetLowA = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("jetLowA")),
  				*jetLowB = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("jetLowB")),
  				*jetLowC = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("jetLowC")),
  				*jetHighA = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("jetHighA")),
  				*jetHighB = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("jetHighB")),
  				*jetHighC = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("jetHighC")),
  				*gammaLowA = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("gammaLowA")),
  				*gammaLowB = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("gammaLowB")),
  				*gammaLowC = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("gammaLowC")),
  				*gammaHighA = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("gammaHighA")),
  				*gammaHighB = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("gammaHighB")),
  				*gammaHighC = static_cast<AliJSONInt *>(combinedConfiguration.GetValue("gammaHighC"));
  eTask->SetTriggerThresholdJetLow( jetLowA->GetValue(), jetLowB->GetValue(), jetLowC->GetValue() );
  eTask->SetTriggerThresholdJetHigh( jetHighA->GetValue(), jetHighB->GetValue(), jetHighC->GetValue() );
  eTask->SetTriggerThresholdJetLow( gammaLowA->GetValue(), gammaLowB->GetValue(), gammaLowC->GetValue() );
  eTask->SetTriggerThresholdJetHigh( gammaHighA->GetValue(), gammaHighB->GetValue(), gammaHighC->GetValue() );

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
