AliEmcalIsolatedPhotonsTask* AddTaskEmcalIsolatedPhotons(
						       const char *ntracks            = "Tracks",
						       const char *nclusters          = "CaloClusters",
						       const char *njets              = "Jets",
						       const Int_t AODtrackFilterBit  = 256  // hybrid LHC11h tracks
                                                       )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalIsolatedPhotons", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalIsolatedPhotons", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliEmcalIsolatedPhotonsTask* phTask = new AliEmcalIsolatedPhotonsTask("aiolaIsoPhotons");
  phTask->SetTracksName(ntracks);
  phTask->SetClusName(nclusters);
  phTask->SetJetsName(njets);
  phTask->SetAODFilterBit(AODtrackFilterBit); // global hybrids for LHC11h

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(phTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histosEMCALIsoPhoton_aiola", 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (phTask, 0,  cinput1 );
  mgr->ConnectOutput (phTask, 1, coutput1 );

  return phTask;
  
}
