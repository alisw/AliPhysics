// $Id$

AliAnalysisTaskSAJF* AddTaskSAJF(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nembtracks         = "TracksEmbedded",
  const char *nembclusters       = "CaloClustersEmbedded",
  const char *nembjets           = "EmbJets",
  const char *nrandtracks        = "TracksRandomized",
  const char *nrandclusters      = "CaloClustersRandomized",
  const char *nrho               = "Rho",
  Double_t    jetradius          = 0.4,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.4,
  Double_t    ptcut              = 0.15,
  Double_t    jetBiasTrack       = 5,
  Double_t    jetBiasClus        = 5,
  UInt_t      type               = AliAnalysisTaskEmcal::kTPC,
  const char *taskname           = "AliAnalysisTaskSAJF"
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

  TString name(taskname);
  name += "_";
  name += njets;
  name += "_";
  name += nrho;
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
  AliAnalysisTaskSAJF* jetTask = new AliAnalysisTaskSAJF(name);
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
  jetTask->SetPtCut(ptcut);
  jetTask->SetJetRadius(jetradius);
  jetTask->SetJetPtCut(jetptcut);
  jetTask->SetJetAreaCut(jetareacut);
  jetTask->SetPtBiasJetTrack(jetBiasTrack);
  jetTask->SetPtBiasJetClus(jetBiasClus);
  
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );
  
  return jetTask;
}
