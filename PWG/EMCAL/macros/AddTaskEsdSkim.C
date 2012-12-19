// $Id$

AliEsdSkimTask* AddTaskEsdSkim(
  Bool_t tof       = kFALSE,
  Bool_t emc       = kFALSE,
  Bool_t emt       = kFALSE,
  Bool_t phc       = kFALSE,
  Bool_t pht       = kFALSE,
  Bool_t clus      = kFALSE,
  Bool_t tracks    = kFALSE,
  Bool_t mtracks   = kFALSE,
  Bool_t ptracks   = kFALSE,
  Bool_t remcov    = kFALSE,
  Bool_t rescov    = kFALSE,
  Bool_t sbytes    = kFALSE,
  Bool_t phosclus  = kFALSE,
  Bool_t emcclus   = kFALSE,
  Bool_t muons     = kFALSE,
  const char *tname = "Tracks",
  const char *filename = "AliSkimmedESD.root"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalSetup", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalSetup", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  AliEsdSkimTask *task = new AliEsdSkimTask("EmcalSkimTask");
  task->SetDoClus(clus);
  task->SetDoEmC(emc);
  task->SetDoEmT(emt);
  task->SetDoMiniTracks(mtracks);
  task->SetDoMuonTracks(muons);
  task->SetDoPhC(phc);
  task->SetDoPhT(pht);
  task->SetDoPicoTracks(ptracks);
  task->SetDoSaveBytes(sbytes);
  task->SetDoTof(tof);
  task->SetDoTracks(tracks);
  task->SetEmcalClusOnly(emcclus);
  task->SetPhosClusOnly(phosclus);
  task->SetRemoveCP(remcov);
  task->SetResetCov(rescov);
  task->SetTracks(tname);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("EsdSkimTree"),
                                                           TTree::Class(),   
                                                           AliAnalysisManager::kOutputContainer, 
                                                           filename);
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}
