// $Id$

AliEmcalTenderTask *AddTaskEmcalAodTender()
{
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
  AliEMCALTenderSupply *EMCALSupply = ConfigEmcalTenderSupply(kFALSE);

  ana->SetEMCALTenderSupply(EMCALSupply);

  // Get and connect common input/output containers via the manager as below
  //==============================================================================

  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
   
  return ana;
}
