// AddTaskEmcalJetSpectraQA.C

AliAnalysisHFjetTagHFE* AddTaskHFjetTagHFE(
  const char *ntracks            = "usedefault",
  const char *nclusters          = "usedefault",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Double_t    jetradius          = 0.3,
  Double_t    jetptcut           = 1,
  Double_t    jetareacut         = 0.2,
  const char *cutType            = "TPCfid",
  Int_t       leadhadtype        = 0,
  const char *suffix             = "",
  Bool_t     iMC                 = kFALSE,
  Bool_t     iNarrowEta          = kFALSE
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskHFjetTagHFE", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskHFjetTagHFE", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(ntracks);
  TString clusName(nclusters);

  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  TString name("AliAnalysisHFjetTagHFE");
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  name += "_";
  name += cutType;

  if (strcmp(suffix,"")) {
    name += "_";
    name += suffix;
  }

  AliAnalysisHFjetTagHFE* jetTask = new AliAnalysisHFjetTagHFE(name);
  jetTask->SetVzRange(-10,10);
  jetTask->SetNeedEmcalGeom(kFALSE);

  Double_t JetEta = 0.9-jetradius;
  if(iNarrowEta)JetEta = 0.6-jetradius;  // reference eta is EMC acc

  cout << "<----------- JetEta =  " << JetEta << endl;
  jetTask->SetJetEtaCut(JetEta);

  AliTrackContainer* trackCont = 0;

  //if (trackName == "mcparticles") {
    //AliMCParticleContainer* mcpartCont = jetTask->AddMCParticleContainer(trackName);
    //mcpartCont->SelectPhysicalPrimaries(kTRUE);
  //}
  //else if (trackName == "tracks" || trackName == "Tracks") {
  if (trackName == "tracks" || trackName == "Tracks") {
    //AliTrackContainer* trackCont = jetTask->AddTrackContainer(trackName);
    trackCont = jetTask->AddTrackContainer(trackName);
    trackCont->SetFilterHybridTracks(kTRUE);
  }
  else if (!trackName.IsNull()) {
    jetTask->AddParticleContainer(trackName);
  }

 /*
  AliParticleContainer *partCont = jetTask->GetParticleContainer(0);
  if (partCont) {
    partCont->SetParticlePtCut(trackPtCut);
  }
*/ 

  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(clusName);
  if (clusterCont) {
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  //AliJetContainer *jetCont = jetTask->AddJetContainer(njets, cutType, jetradius);
  AliJetContainer* jetCont = jetTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, jetradius, AliJetContainer::kTPCfid, "Jet");
  if (jetCont) {
    //jetCont->SetRhoName(nrho);
    jetCont->SetRhoName("Rho");
    cout << "Name of Rho " << jetCont->GetRhoName() << endl;
    //if(jetradius==0.3)jetareacut=0.2;
    jetareacut = jetradius*jetradius*TMath::Pi()*0.6;
    cout << "jetradius = " << jetradius << " ; jetareacut = " << jetareacut << endl;  
    //+++ jetCont->SetPercAreaCut(jetareacut);
    jetCont->SetJetAreaCut(jetareacut);
    jetCont->SetJetPtCut(jetptcut);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
    jetCont->SetLeadingHadronType(leadhadtype);
    jetCont->SetMaxTrackPt(1000);
    jetCont->SetZLeadingCut(0.98,0.98);
    }

   if(iMC)
     {
     //AliTrackContainer* trackContMC = jetTask->AddTrackContainer("mcparticles");
      AliMCParticleContainer* trackContMC = jetTask->AddMCParticleContainer("mcparticles");
      //AliJetContainer* jetContMC = jetTask->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, jetradius, AliJetContainer::kTPCfid, "JetMC");
      //AliJetContainer* jetContMC = jetTask->AddJetContainer("JetMC_AKTChargedR030_mcparticles_pT0150_pt_scheme");
      AliJetContainer* jetContMC;
      if(jetradius==0.3)jetContMC = jetTask->AddJetContainer("JetMC_AKTChargedR030_mcparticles_pT0150_pt_scheme");
      if(jetradius==0.2)jetContMC = jetTask->AddJetContainer("JetMC_AKTChargedR020_mcparticles_pT0150_pt_scheme");
      if(jetradius==0.4)jetContMC = jetTask->AddJetContainer("JetMC_AKTChargedR040_mcparticles_pT0150_pt_scheme");
      if(jetradius==0.6)jetContMC = jetTask->AddJetContainer("JetMC_AKTChargedR060_mcparticles_pT0150_pt_scheme");
     
      if (jetContMC) {
      //jetCont->SetRhoName(nrho);
      //if(jetradius==0.3)jetareacut=0.2;
      jetContMC->SetJetAreaCut(jetareacut);
      jetContMC->SetJetPtCut(jetptcut);
      jetContMC->ConnectParticleContainer(trackContMC);
      jetContMC->ConnectClusterContainer(clusterCont);
      jetContMC->SetLeadingHadronType(leadhadtype);
      jetContMC->SetMaxTrackPt(1000);
      jetContMC->SetZLeadingCut(0.98,0.98);
     }
   }
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
