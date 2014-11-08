// $Id$

AliAnalysisTaskEmcalHighMultTrigger* AddTaskEmcalHighMultTrigger(
  const char *ntracks            = "Tracks",
  const char *nPatches           = "EmcalPatches32x32",
  Int_t       nCentBins          = 1,
  const char *taskname           = "HighMultTrigger"
)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalHighMultTrigger", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalHighMultTrigger", "This task requires an input event handler");
    return NULL;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name = Form("%s_%s",taskname,nPatches);

  Printf("name: %s",name.Data());

  AliAnalysisTaskEmcalHighMultTrigger* task = new AliAnalysisTaskEmcalHighMultTrigger(name);
  task->SetCentRange(0.,100.);
  task->SetNCentBins(nCentBins);
  task->SetCaloTriggerPatchInfoName(nPatches);
  task->SetVzRange(-10.,10.);
  task->SetUseAliAnaUtils(kTRUE,kTRUE);//kFALSE);

  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  trackCont->SetClassName("AliVTrack");

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (task, 0,  cinput1 );
  mgr->ConnectOutput (task, 1, coutput1 );

  return task;
}

