// $Id$

AliJetResponseMaker* AddTaskJetResponseMaker(
  const char *ntracks1           = "Tracks",
  const char *nclusters1         = "CaloClusters",
  const char *njets1             = "Jets",
  const char *nrho1              = "Rho",
  Double_t    jetradius1         = 0.2,
  const char *ntracks2           = "MCParticles",
  const char *nclusters2         = "",
  const char *njets2             = "MCJets",
  const char *nrho2              = "",
  Double_t    jetradius2         = 0.2,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.557,
  Double_t    jetBias            = 5,
  Int_t       biasType           = 0,   //  0 = charged, 1 = neutral, 2 = both
  UInt_t      matching           = AliJetResponseMaker::kGeometrical,
  Double_t    maxDistance1       = 0.25,
  Double_t    maxDistance2       = 0.25,
  const char *cutType            = "TPC",
  Int_t       ptHardBin          = -999,
  Double_t    minCent            = -999,
  Double_t    maxCent            = -999,
  const char *taskname           = "AliJetResponseMaker",
  Bool_t      biggerMatrix       = kFALSE,
  AliJetResponseMaker* address   = 0,
  Double_t    nefmincut          = -10,
  Double_t    nefmaxcut          = 10,
  Int_t       jetTagging         = 0,
  Double_t    maxTrackPt         = 100
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJetResponseMaker", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetResponseMaker", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("%s_%s_%s_Bias%d_BiasType%d_%s",taskname,njets1,njets2,(Int_t)floor(jetBias),biasType,cutType));

  if (minCent != -999 && maxCent != -999) 
    name += Form("_Cent%d_%d", (Int_t)floor(minCent), (Int_t)floor(maxCent));

  if (ptHardBin != -999) 
    name += Form("_PtHard%d", ptHardBin);
  if (nefmaxcut<1.0)
    name += Form("_NEF%d", (Int_t)(100*nefmaxcut));
    
  AliJetResponseMaker* jetTask = address;
  if (jetTask)
    new (jetTask) AliJetResponseMaker(name);
  else
    jetTask = new AliJetResponseMaker(name);

  AliParticleContainer *trackCont1 = jetTask->AddParticleContainer(ntracks1);
  AliClusterContainer *clusCont1 = jetTask->AddClusterContainer(nclusters1);
  AliJetContainer *jetCont1 = jetTask->AddJetContainer(njets1, cutType, jetradius1);
  jetCont1->SetRhoName(nrho1);
  jetCont1->SetLeadingHadronType(biasType);
  jetCont1->SetPtBiasJetTrack(jetBias);
  jetCont1->SetPtBiasJetClus(jetBias);
  jetCont1->SetJetPtCut(jetptcut);
  jetCont1->SetPercAreaCut(jetareacut);
  jetCont1->SetIsParticleLevel(kFALSE);
  jetCont1->ConnectParticleContainer(trackCont1);
  jetCont1->ConnectClusterContainer(clusCont1);
  jetCont1->SetNEFCut(nefmincut,nefmaxcut);
  jetCont1->SetFlavourCut(jetTagging);
  jetCont1->SetMaxTrackPt(maxTrackPt);

    
  AliParticleContainer *trackCont2 = jetTask->AddParticleContainer(ntracks2);
  trackCont2->SetParticlePtCut(0);
  AliClusterContainer *clusCont2 = jetTask->AddClusterContainer(nclusters2);
  AliJetContainer *jetCont2 = jetTask->AddJetContainer(njets2, cutType, jetradius2);
  jetCont2->SetRhoName(nrho2);
  jetCont2->SetLeadingHadronType(biasType);
  jetCont2->SetPtBiasJetTrack(jetBias);
  jetCont2->SetPtBiasJetClus(jetBias);
  jetCont2->SetJetPtCut(jetptcut);
  jetCont2->SetPercAreaCut(jetareacut);
  jetCont2->SetIsParticleLevel(kTRUE);
  jetCont2->ConnectParticleContainer(trackCont2);
  jetCont2->ConnectClusterContainer(clusCont2);
  jetCont2->SetFlavourCut(jetTagging);
  jetCont2->SetMaxTrackPt(1000); // disable default 100 GeV/c track cut for particle level jets

    
  jetTask->SetMatching(matching, maxDistance1, maxDistance2);
  jetTask->SetVzRange(-10,10);
  jetTask->SetIsPythia(kTRUE);
  jetTask->SetPtHardBin(ptHardBin);
  jetTask->SetCentRange(minCent,maxCent);

  if (biggerMatrix) 
    jetTask->SetHistoBins(1000,0,500);
  
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
