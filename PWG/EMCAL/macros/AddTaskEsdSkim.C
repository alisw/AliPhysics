// $Id$

AliEsdSkimTask* AddTaskEsdSkim(
  Bool_t tof       = kTRUE,
  Bool_t emc       = kTRUE,
  Bool_t emt       = kFALSE,
  Bool_t phc       = kTRUE,
  Bool_t pht       = kFALSE,
  Bool_t clus      = kTRUE,
  Bool_t tracks    = kTRUE,
  Bool_t atracks   = kTRUE,
  Bool_t mtracks   = kFALSE,
  Bool_t ptracks   = kFALSE,
  Bool_t v0s       = kFALSE,
  Bool_t mult      = kTRUE,
  Bool_t fmd       = kFALSE,
  Bool_t muons     = kFALSE,
  Bool_t remcov    = kFALSE,
  Bool_t rescov    = kFALSE,
  Bool_t sbytes    = kFALSE,
  Bool_t phosclus  = kFALSE,
  Bool_t emcclus   = kFALSE,
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
  task->SetDoAllTracks(atracks);
  task->SetDoV0s(v0s);
  task->SetDoMult(mult);
  task->SetDoFmd(fmd);
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
