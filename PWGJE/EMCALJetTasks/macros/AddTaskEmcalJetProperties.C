AliAnalysisTaskEmcalJetProperties* AddTaskEmcalJetProperties(
							     //const char *ntracks            = "mcparticles",
							     const char *ntracks            = "usedefault",
							     const char *nclusters          = "usedefault",
							     const char* ncells             = "usedefault",
							     const char *suffix             = ""
							     )
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      ::Error("AddTaskEmcalJetProperties", "No analysis manager to connect to.");
      return 0;
    }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
    {
      ::Error("AddTaskEmcalJetProperties", "This task requires an input event handler");
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
  TString cellName(ncells);
  
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

  if (cellName == "usedefault") {
    if (dataType == kESD) {
      cellName = "EMCALCells";
    }
    else if (dataType == kAOD) {
      cellName = "emcalCells";
    }
    else {
      cellName = "";
    }
  }

  TString name("AliAnalysisTaskEmcalJetProperties");
  if (!trackName.IsNull()) {
    name += "_";
    name += trackName;
  }
  if (!clusName.IsNull()) {
    name += "_";
    name += clusName;
  }
  if (!cellName.IsNull()) {
    name += "_";
    name += cellName;
  }
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

  AliAnalysisTaskEmcalJetProperties* sampleTask = new AliAnalysisTaskEmcalJetProperties(name);
  sampleTask->SetCaloCellsName(cellName);
  sampleTask->SetVzRange(-10,10);

  if (trackName == "mcparticles") {
    AliMCParticleContainer* mcpartCont = sampleTask->AddMCParticleContainer(trackName);
  }
  else if (trackName == "tracks" || trackName == "Tracks") {
    AliTrackContainer* trackCont = sampleTask->AddTrackContainer(trackName);
  }
  else if (!trackName.IsNull()) {
    sampleTask->AddParticleContainer(trackName);
  }
  sampleTask->AddClusterContainer(clusName);
  
  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(sampleTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(),
      TList::Class(),AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (sampleTask, 0,  cinput1 );
  mgr->ConnectOutput (sampleTask, 1, coutput1 );

  return sampleTask;
}

//*****************

AliAnalysisTaskEmcalJetProperties* AddTaskEmcalJetProperties(const char * njetsKtRec, 
							     const char * njetsKtGen, 
							     const char * njetsAktRec, 
							     const char * njetsAktGen, 
							     const Double_t R,
							     const char * type,
							     //const char *ntracks,
							     //const char *nclusters,
							     Int_t pSel
							     )
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalJetTagger","No analysis manager found.");
      return 0;
    }
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalJetTagger", "This task requires an input event handler");
      return NULL;
    }
  
  TString wagonName = Form("JetTagger_%s_%s_TC",njetsAktRec,njetsAktGen);

  //Configure jet tagger task
  AliAnalysisTaskEmcalJetProperties *task = new AliAnalysisTaskEmcalJetProperties(wagonName);
  
  AliTrackContainer *trackCont  = task->AddTrackContainer("tracks");
  AliClusterContainer *clusterCont = task->AddClusterContainer("caloClusters");
  AliMCParticleContainer* mcpartCont = task->AddMCParticleContainer("mcparticles");
  
  task->SetJetContainerKtRec(0);
  task->SetJetContainerKtGen(1);
  task->SetJetContainerAktRec(2);
  task->SetJetContainerAktGen(3);

  TString strType(type);
  AliJetContainer *jetContKtRec = task->AddJetContainer(njetsKtRec,strType,R);
  if(jetContKtRec) {
    //jetContKtRec->SetRhoName(nrhoBase);
    jetContKtRec->ConnectParticleContainer(trackCont);
    jetContKtRec->ConnectClusterContainer(clusterCont);
    jetContKtRec->SetMaxTrackPt(10000.);
  }

  AliJetContainer *jetContKtGen = task->AddJetContainer(njetsKtGen,strType,R);
  if(jetContKtGen) {
    //jetContKtGen->SetRhoName(nrhoBase);
    jetContKtGen->ConnectParticleContainer(mcpartCont);
    jetContKtGen->ConnectClusterContainer(clusterCont);
    jetContKtGen->SetMaxTrackPt(10000.);
  }
  AliJetContainer *jetContAktRec = task->AddJetContainer(njetsAktRec,strType,R);
  if(jetContAktRec) {
    //jetContAktRec->SetRhoName(nrhoBase);
    jetContAktRec->ConnectParticleContainer(trackCont);
    jetContAktRec->ConnectClusterContainer(clusterCont);
    jetContAktRec->SetMaxTrackPt(10000.);
  }
  AliJetContainer *jetContAktGen = task->AddJetContainer(njetsAktGen,strType,R);
  if(jetContAktGen) {
    //jetContAktGen->SetRhoName(nrhoBase);
    jetContAktGen->ConnectParticleContainer(mcpartCont);
    jetContAktGen->ConnectClusterContainer(clusterCont);
    jetContAktGen->SetMaxTrackPt(10000.);
  }

  for(Int_t i=0; i<4; i++) {
    task->SetPercAreaCut(0.6, i); //keep?
  }
  
  task->SelectCollisionCandidates(pSel);
  task->SetCaloCellsName("emcalCells");
  task->SetVzRange(-10,10);

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName(wagonName);
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);

  return task;  
}




