// $Id$

AliAnalysisTaskDeltaPt* AddTaskDeltaPt(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nembtracks         = "TracksEmbedded",
  const char *nembclusters       = "CaloClustersEmbedded",
  const char *nembjets           = "EmbJets",
  const char *nrandtracks        = "TracksRandomized",
  const char *nrandclusters      = "CaloClustersRandomized",
  const char *nrho               = "Rho",
  Double_t    jetradius          = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.557,
  Double_t    trackptcut         = 0.15,
  Double_t    clusptcut          = 0.30,
  UInt_t      type               = AliAnalysisTaskEmcal::kTPC,
  const char *taskname           = "AliAnalysisTaskDeltaPt"
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
  TString name;
  if (strcmp(ntracks, "") == 0 && strcmp(nclusters, "") == 0) 
    name = Form("%s_%s_R0%d",taskname,nrho,(Int_t)floor(jetradius*100+0.5));
  else if (strcmp(ntracks, "") == 0) 
    name = Form("%s_%s_%s_R0%d",taskname,nclusters,nrho,(Int_t)floor(jetradius*100+0.5));
  else if (strcmp(nclusters, "") == 0) 
    name = Form("%s_%s_%s_R0%d",taskname,ntracks,nrho,(Int_t)floor(jetradius*100+0.5));
  else
    name = Form("%s_%s_%s_%s_R0%d",taskname,ntracks,nclusters,nrho,(Int_t)floor(jetradius*100+0.5));

  if (type == AliAnalysisTaskEmcal::kTPC) 
    name += "_TPC";
  else if (type == AliAnalysisTaskEmcal::kEMCAL) 
    name += "_EMCAL";
  else if (type == AliAnalysisTaskEmcal::kUser) 
    name += "_USER";

  AliAnalysisTaskDeltaPt* jetTask = new AliAnalysisTaskDeltaPt(name);
  jetTask->SetAnaType(type);
  jetTask->SetTracksName(ntracks);
  jetTask->SetClusName(nclusters);
  jetTask->SetJetsName(njets);
  jetTask->SetEmbTracksName(nembtracks);
  jetTask->SetEmbClusName(nembclusters);
  jetTask->SetEmbJetsName(nembjets);
  jetTask->SetRandTracksName(nrandtracks);
  jetTask->SetRandClusName(nrandclusters);
  jetTask->SetRhoName(nrho);
  jetTask->SetClusPtCut(clusptcut);
  jetTask->SetTrackPtCut(trackptcut);
  jetTask->SetJetRadius(jetradius);
  jetTask->SetJetPtCut(jetptcut);
  jetTask->SetPercAreaCut(jetareacut);
  
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer();
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput(jetTask, 0, cinput1);
  mgr->ConnectOutput(jetTask, 1, coutput1);
  
  return jetTask;
}
