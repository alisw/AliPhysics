// $Id$

AliAnalysisTaskSAQA* AddTaskSAQA(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.557,
  Double_t    trackptcut         = 0.15,
  Double_t    clusptcut          = 0.30,
  UInt_t      type               = AliAnalysisTaskEmcal::kTPC,
  const char *taskname           = "AliAnalysisTaskSAQA"
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
  
  TString name(taskname);
  if (strcmp(ntracks,"")) {
    name += "_";
    name += ntracks;
  }
  if (strcmp(nclusters,"")) {
    name += "_";
    name += nclusters;
  }
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (type == AliAnalysisTaskEmcal::kTPC) 
    name += "_TPC";
  else if (type == AliAnalysisTaskEmcal::kEMCAL) 
    name += "_EMCAL";
  else if (type == AliAnalysisTaskEmcal::kUser) 
    name += "_USER";

  AliAnalysisTaskSAQA* qaTask = new AliAnalysisTaskSAQA(name);
  qaTask->SetTracksName(ntracks);
  qaTask->SetClusName(nclusters);
  qaTask->SetJetsName(njets);
  qaTask->SetJetRadius(jetradius);
  qaTask->SetJetPtCut(jetptcut);
  qaTask->SetPercAreaCut(jetareacut);
  qaTask->SetTrackPtCut(trackptcut);
  qaTask->SetClusPtCut(clusptcut);
  qaTask->SetAnaType(type);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(qaTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;

  TString contName(name);
  contName += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (qaTask, 0,  cinput1 );
  mgr->ConnectOutput (qaTask, 1, coutput1 );

  return qaTask;
}
