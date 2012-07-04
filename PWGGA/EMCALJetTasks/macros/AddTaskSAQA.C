// $Id: AddTaskSAQA.C 56113 2012-05-01 21:39:27Z loizides $

AliAnalysisTaskSAQA* AddTaskSAQA(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *ntrgclusters       = "",
  Double_t    jetradius          = 0.4,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.8,
  Double_t    ptcut              = 0.15,
  Double_t    jetBiasTrack       = 5,
  Double_t    jetBiasClus        = 5,
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
  name += "_";
  name += njets;
  name += "_Track";
  name += jetBiasTrack;
  name += "_Clus";
  name += jetBiasClus;
  name += "_R0";
  name += floor(jetradius*10+0.5);
  name += "_PtCut";
  name += floor(ptcut*1000+0.5);
  name += "_";
  if (type == AliAnalysisTaskEmcal::kTPC) 
    name += "TPC";
  else if (type == AliAnalysisTaskEmcal::kEMCAL) 
    name += "EMCAL";
  else if (type == AliAnalysisTaskEmcal::kTPCSmall) 
    name += "TPCSmall";
  else if (type == AliAnalysisTaskEmcal::kEMCALOnly) 
    name += "EMCALOnly";
  AliAnalysisTaskSAQA* qaTask = new AliAnalysisTaskSAQA(name);
  qaTask->SetTracksName(ntracks);
  qaTask->SetClusName(nclusters);
  qaTask->SetJetsName(njets);
  qaTask->SetTrgClusName(ntrgclusters);
  qaTask->SetJetRadius(jetradius);
  qaTask->SetJetPtCut(jetptcut);
  qaTask->SetPercAreaCut(jetareacut);
  qaTask->SetPtCut(ptcut);
  qaTask->SetPtBiasJetTrack(jetBiasTrack);
  qaTask->SetPtBiasJetClus(jetBiasClus);
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
