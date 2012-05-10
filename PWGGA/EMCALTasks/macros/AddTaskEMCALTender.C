// $Id$

AliTender *AddTaskEMCALTender()

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALTender", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Create the task and configure it.
  //===========================================================================
  AliTender* ana = new  AliTender("AliTender");
  
  Bool_t ismc = (mgr->GetMCtruthEventHandler() != NULL);

  mgr->AddTask(ana);

  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/ConfigEmcalTenderSupply.C");
  AliEMCALTenderSupply *EMCALSupply = ConfigEmcalTenderSupply(kTRUE);
 
  ana->AddSupply(EMCALSupply);

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
