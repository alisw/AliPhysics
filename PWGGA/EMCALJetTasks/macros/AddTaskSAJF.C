// $Id$

AliAnalysisTaskSAJF* AddTaskSAJF(
  const char *taskname           = "AliAnalysisTaskSAJF",
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nktjets            = "KtJets",
  const char *nembjets           = "EmbJets",
  const char *nrho               = "Rho",
  Double_t    jetradius          = 0.4,
  Double_t    ptcut              = 0.15,
  Double_t    jetpartcut         = 10,
  UInt_t      type               = AliAnalysisTaskSAJF::kTPC
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskSAJF", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskSAJF", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliAnalysisTaskSAJF* phTask = new AliAnalysisTaskSAJF(taskname);
  phTask->SetAnaType(type);
  phTask->SetTracksName(ntracks);
  phTask->SetClusName(nclusters);
  phTask->SetJetsName(njets);
  phTask->SetKtJetsName(nktjets);
  phTask->SetEmbJetsName(nembjets);
  phTask->SetRhoName(nrho);
  phTask->SetPtCut(ptcut);
  phTask->SetJetRadius(jetradius);
  phTask->SetPtCutJetPart(jetpartcut);
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(phTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(taskname);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (phTask, 0,  cinput1 );
  mgr->ConnectOutput (phTask, 1, coutput1 );

  return phTask;
}
