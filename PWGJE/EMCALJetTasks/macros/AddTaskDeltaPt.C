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
  Double_t    jetareacut         = 0.557,
  Double_t    trackptcut         = 0.15,
  Double_t    clusptcut          = 0.30,
  const char *type               = "TPC",
  const char *taskname           = "AliAnalysisTaskDeltaPt"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskDeltaPt", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskDeltaPt", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString name;
  if (strcmp(ntracks, "") == 0 && strcmp(nclusters, "") == 0) 
    name = Form("%s_%s_R0%d_%s",taskname,nrho,(Int_t)floor(jetradius*100+0.5),type);
  else if (strcmp(ntracks, "") == 0) 
    name = Form("%s_%s_%s_R0%d_%s",taskname,nclusters,nrho,(Int_t)floor(jetradius*100+0.5),type);
  else if (strcmp(nclusters, "") == 0) 
    name = Form("%s_%s_%s_R0%d_%s",taskname,ntracks,nrho,(Int_t)floor(jetradius*100+0.5),type);
  else
    name = Form("%s_%s_%s_%s_R0%d_%s",taskname,ntracks,nclusters,nrho,(Int_t)floor(jetradius*100+0.5),type);

  AliAnalysisTaskDeltaPt* jetTask = new AliAnalysisTaskDeltaPt(name);
  jetTask->SetConeRadius(jetradius);
  jetTask->SetRhoName(nrho,-1);
  if (strcmp(type,"TPC")==0) 
    jetTask->SetConeEtaPhiTPC();
  else if (strcmp(type,"EMCAL")==0) 
    jetTask->SetConeEtaPhiEMCAL();

  AliParticleContainer *partCont = jetTask->AddParticleContainer(ntracks);
  if (partCont) {
    partCont->SetName("Tracks");
    partCont->SetParticlePtCut(trackptcut);
  }

  AliClusterContainer *clusCont = jetTask->AddClusterContainer(nclusters);
  if (clusCont) {
    clusCont->SetName("CaloClusters");
    clusCont->SetClusPtCut(clusptcut);
  }

  AliJetContainer *jetCont = jetTask->AddJetContainer(njets,type,jetradius);
  if (jetCont) {
    jetCont->SetName("Jets");
    jetCont->SetPercAreaCut(jetareacut);
    jetCont->SetRhoName(nrho);
    jetCont->ConnectParticleContainer(partCont);
    jetCont->ConnectClusterContainer(clusCont);
  }

  AliParticleContainer *embPartCont = jetTask->AddParticleContainer(nembtracks);
  if (embPartCont) {
    embPartCont->SetName("EmbTracks");
    embPartCont->SetParticlePtCut(trackptcut);
  }

  AliClusterContainer *embClusCont = jetTask->AddClusterContainer(nembclusters);
  if (embClusCont) {
    embClusCont->SetName("EmbClusters");
    embClusCont->SetClusPtCut(clusptcut);
  }

  AliJetContainer *embJetCont = jetTask->AddJetContainer(nembjets,type,jetradius);
  if (embJetCont) {
    embJetCont->SetName("EmbJets");
    embJetCont->SetPercAreaCut(jetareacut);
    embJetCont->SetRhoName(nrho);
    embJetCont->ConnectParticleContainer(embPartCont);
    embJetCont->ConnectClusterContainer(embClusCont);
  }

  AliParticleContainer *randPartCont = jetTask->AddParticleContainer(nrandtracks);
  if (randPartCont) {
    randPartCont->SetName("RandTracks");
    randPartCont->SetParticlePtCut(trackptcut);
  }

  AliClusterContainer *randClusCont = jetTask->AddClusterContainer(nrandclusters);    
  if (randClusCont) {
    randClusCont->SetName("RandClusters");
    randClusCont->SetClusPtCut(clusptcut);
  }
  
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
