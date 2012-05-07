// $Id$

AliEmcalTenderTask *AddTaskEmcalAodTender(const char *geoname="EMCAL_COMPLETEV1", const char* datatype="pp")
{
  // Parameters: geoname = "EMCAL_FIRSTYEARV1" or "EMCAL_COMPLETEV1" or ""

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskTrgContam", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Create the task and configure it.
  //===========================================================================
  AliEmcalTenderTask* ana = new  AliEmcalTenderTask("AliEmcalTenderTask");
  
  Bool_t ismc = (mgr->GetMCtruthEventHandler() != NULL);

  mgr->AddTask(ana);

  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/ConfigEmcalTenderSupply.C");
  AliEMCALTenderSupply *EMCALSupply = ConfigEmcalTenderSupply(geoname, kFALSE);

  ana->SetEMCALTenderSupply(EMCALSupply);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  /*AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("tender_event", 
                         AliESDEvent::Class(), 
                         AliAnalysisManager::kExchangeContainer,
                         "default_tender");
  */
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  //mgr->ConnectOutput (ana, 1, coutput1 );
   
  return ana;
}
