AliMCHFParticleSelector* AddTaskMCHFParticleSelector(
  const char *outname    = "mcparticles",
  Bool_t      nk         = kFALSE,
  Bool_t      ch         = kFALSE,
  Double_t    etamax     = 1,
  Bool_t      physPrim   = kTRUE
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskMCHFParticleSelector", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskMCHFParticleSelector", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString name("AliMCHFParticleSelector_");
  name += outname;
  AliMCHFParticleSelector *eTask = new AliMCHFParticleSelector(name);
  eTask->SetParticlesOutName(outname);
  eTask->SetRejectNK(nk);
  eTask->SetChargedMC(ch);
  eTask->SetEtaMax(etamax);
  eTask->SetOnlyPhysPrim(physPrim);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(eTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  mgr->ConnectInput  (eTask, 0,  cinput1 );
  
  return eTask;
}
