// $Id: AddTaskSAQA.C 56113 2012-05-01 21:39:27Z loizides $

AliAnalysisTaskSAQA* AddTaskSAQA(
  const char *taskname           = "AliAnalysisTaskSAQA",
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *ntrgclusters       = "ClustersL1GAMMAFEE",
  Double_t    jetradius          = 0.4,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.1,
  UInt_t      type               = AliAnalysisTaskEmcal::kTPC
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskSAQA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskSAQA", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  AliAnalysisTaskSAQA* qaTask = new AliAnalysisTaskSAQA(taskname);
  qaTask->SetTracksName(ntracks);
  qaTask->SetClusName(nclusters);
  qaTask->SetJetsName(njets);
  qaTask->SetTrgClusName(ntrgclusters);
  qaTask->SetJetRadius(jetradius);
  qaTask->SetJetPtCut(jetptcut);
  qaTask->SetJetAreaCut(jetareacut);
  qaTask->SetAnaType(type);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(qaTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;

  TString contName(taskname);
  contName += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (qaTask, 0,  cinput1 );
  mgr->ConnectOutput (qaTask, 1, coutput1 );

  return qaTask;
}
