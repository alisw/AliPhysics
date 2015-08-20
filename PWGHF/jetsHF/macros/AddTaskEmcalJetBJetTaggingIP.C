AliAnalysisTaskEmcalJetBJetTaggingIP* AddTaskEmcalJetBJetTaggingIP(const char *ntracks            = "PicoTracks",
								   const char *nclusters           = "",
								   const char *njets              = "Jets",

								   const char *nrho               = "Rho",
								   Double_t jetradius =0.4,
								   const char * type = "TPC",	   
								   const char *taskname           = "AliAnalysisTaskEmcalJetBJetTaggingIP",
								   Bool_t isMC = kFALSE,
								   const char *njetsMC              = "Jets",
								   const char *nrhoMC               = "RhoMC",
								   Bool_t doTrackQAEvent =kTRUE,
								   Bool_t doTrackQAJet =kTRUE,
								   Bool_t doBackgroundFluctuations =kTRUE,
								   Bool_t doCorrectPt =kTRUE)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetSample", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalJetSample", "This task requires an input event handler");
    return NULL;
  }

  
  TString name(taskname);
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  if (strcmp(type,"")) {
    name += "_";
    name += type;
  }

  Printf("name: %s",name.Data());
  AliAnalysisTaskEmcalJetBJetTaggingIP* jetTask = new AliAnalysisTaskEmcalJetBJetTaggingIP(taskname);
  AliParticleContainer *trackCont  = jetTask->AddParticleContainer(ntracks);
  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);

  TString strType(type);
  AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);
    if(jetCont) {
    jetCont->SetRhoName(nrho);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
  }

    if(isMC)
      {
  AliJetContainer *jetContMC = jetTask->AddJetContainer(njetsMC,strType,jetradius);
    if(jetContMC) {
    jetContMC->SetRhoName(nrhoMC);
    jetContMC->SetIsParticleLevel(kTRUE);
    jetContMC->SetMaxTrackPt(1000);
    }


      
      }
    //-------------------------------------------------------
    //  Configure analysis task
    //-------------------------------------------------------
    jetTask->SetMC(isMC);
    jetTask->SetIsPythia(isMC);
    jetTask->SetDoTrackQA(doTrackQAEvent);
    jetTask->SetDoTrackQAConstituent(doTrackQAJet);
    jetTask->SetDoBackgroundFluctuations(doBackgroundFluctuations);
    jetTask->SetUseCorrectedPt(doCorrectPt);

    
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
