// $Id$

AliAnalysisTaskDijetHadron* AddTaskDijetHadron(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nMCtracks          = "TracksMC",
  const char *nMCclusters        = "CaloClustersMC",
  const char *nMCjets            = "JetsMC",
  const char *nembtracks         = "TracksEmbedded",
  const char *nembclusters       = "CaloClustersEmbedded",
  const char *nembjets           = "EmbJets",
  const char *nrandtracks        = "TracksRandomized",
  const char *nrandclusters      = "CaloClustersRandomized",
  const char *nPbPbrho           = "Rho",
  const char *nMCrho             = "RhoMC",
  const char *nEMBrho            = "RhoEMB",
  Double_t    jetradius          = 0.2,
  Double_t    leadinghadron1     = 0.0,
  Double_t    leadinghadron2     = 3.0,
  Double_t    leadinghadron3     = 5.0,
  Double_t    jet1pt1            = 10.0,
  Double_t    jet1pt2            = 20.0,
  Double_t    jet1pt3            = 30.0,
  Double_t    jet2pt1            = 10.0,
  Double_t    jet2pt2            = 20.0,
  Double_t    jet2pt3            = 30.0,
  Double_t    jetareacut         = 0.557,
  Double_t    trackptcut         = 0.15,
  Double_t    clusptcut          = 0.30,
  const char *type               = "TPC",
  const char *taskname           = "AliAnalysisTaskDijetHadron"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskDijetHadron", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskDijetHadron", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString name;
  if (strcmp(ntracks, "") == 0 && strcmp(nclusters, "") == 0) 
    name = Form("%s_%s_R0%d_%s",taskname,nPbPbrho,(Int_t)floor(jetradius*100+0.5),type);
  else if (strcmp(ntracks, "") == 0) 
    name = Form("%s_%s_%s_R0%d_%s",taskname,nclusters,nPbPbrho,(Int_t)floor(jetradius*100+0.5),type);
  else if (strcmp(nclusters, "") == 0) 
    name = Form("%s_%s_%s_R0%d_%s",taskname,ntracks,nPbPbrho,(Int_t)floor(jetradius*100+0.5),type);
  else
    name = Form("%s_%s_%s_%s_R0%d_%s",taskname,ntracks,nclusters,nPbPbrho,(Int_t)floor(jetradius*100+0.5),type);

  AliAnalysisTaskDijetHadron* jetTask = new AliAnalysisTaskDijetHadron(name);
  jetTask->SetConeRadius(jetradius);
  jetTask->SetLeadingHadronPtThreshold1(leadinghadron1);
  jetTask->SetLeadingHadronPtThreshold2(leadinghadron2);
  jetTask->SetLeadingHadronPtThreshold3(leadinghadron3);
  jetTask->SetJet1PtThreshold1(jet1pt1);
  jetTask->SetJet1PtThreshold2(jet1pt2);
  jetTask->SetJet1PtThreshold3(jet1pt3);
  jetTask->SetJet2PtThreshold1(jet2pt1);
  jetTask->SetJet2PtThreshold2(jet2pt2);
  jetTask->SetJet2PtThreshold3(jet2pt3);
  jetTask->SetRhoName(nPbPbrho,-1);
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
    jetCont->SetRhoName(nPbPbrho);
    jetCont->ConnectParticleContainer(partCont);
    jetCont->ConnectClusterContainer(clusCont);
  }

  AliParticleContainer *MCpartCont = jetTask->AddParticleContainer(nMCtracks);
  if (partCont) {
    MCpartCont->SetName("MCTracks");
    MCpartCont->SetParticlePtCut(trackptcut);
  }

  AliClusterContainer *MCclusCont = jetTask->AddClusterContainer(nMCclusters);
  if (clusCont) {
    MCclusCont->SetName("MCCaloClusters");
    MCclusCont->SetClusPtCut(clusptcut);
  }

  AliJetContainer *MCjetCont = jetTask->AddJetContainer(nMCjets,type,jetradius);
  if (jetCont) {
    MCjetCont->SetName("MCJets");
    MCjetCont->SetPercAreaCut(jetareacut);
    MCjetCont->SetRhoName(nMCrho);
    MCjetCont->ConnectParticleContainer(MCpartCont);
    MCjetCont->ConnectClusterContainer(MCclusCont);
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
    embJetCont->SetRhoName(nEMBrho);
    embJetCont->ConnectParticleContainer(embPartCont);
    embJetCont->ConnectClusterContainer(embClusCont);
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
