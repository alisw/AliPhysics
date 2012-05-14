// $Id$

AliAnalysisTaskSE *AddTaskEMCALTender()
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALTender", "No analysis manager to connect to.");
    return NULL;
  }

  AliVEventHandler *evhand = mgr->GetInputEventHandler();

  // Create the task and configure it.
  //===========================================================================

  AliAnalysisTaskSE *ana = 0;
  AliEMCALTenderSupply *EMCALSupply = 0;

  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/ConfigEmcalTenderSupply.C");

  if (evhand->InheritsFrom("AliESDInputHandler")) {
    EMCALSupply = ConfigEmcalTenderSupply(kTRUE);

    AliTender* alitender = new  AliTender("AliTender");
    alitender->AddSupply(EMCALSupply);
    ana = alitender;
  }
  else if (evhand->InheritsFrom("AliAODInputHandler")) {
    EMCALSupply = ConfigEmcalTenderSupply(kFALSE);

    AliEmcalTenderTask* emcaltender = new  AliEmcalTenderTask("AliEmcalTenderTask");
    emcaltender->SetEMCALTenderSupply(EMCALSupply);
    ana = emcaltender;
  }
  else {
    ::Error("AddTaskEMCALTender", "Input event handler not recognized, AOD/ESD expected. Returning...");
    return NULL;
  }

  mgr->AddTask(ana);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("tender_event", 
                         AliESDEvent::Class(), 
                         AliAnalysisManager::kExchangeContainer,
                         "default_tender");
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1 );
   
  return ana;
}
