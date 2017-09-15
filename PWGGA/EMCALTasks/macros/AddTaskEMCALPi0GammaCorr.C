// $Id$
AliAnalysisTaskEMCALPi0GammaCorr* AddTaskEMCALPi0GammaCorr(Bool_t isMC=kFALSE)
{  
  std::cout << "Beggining AddTaskEMCALPi0GammaCorr" << std::endl;
  if(isMC) std::cout << "Will analyze MC" << std::endl;
  UInt_t      evtTriggerType         = AliVEvent::kEMCEGA; //AliVEvent::kAnyINT,// AliVEvent::kEMCEGA,//..use this type of events to combine gammas(trigger) with hadrons
  UInt_t      evtMixingType          = AliVEvent::kAnyINT;//..use only this type of events to fill your mixed event pool with tracks

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)  ::Error("AddTask", "No analysis manager to connect to.");
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) ::Error("AddTask", "This task requires an input event handler");
 
  //-------------------------------------------------------
  // Built the name of the Task together
  //-------------------------------------------------------
  const char *taskname               = "AliAnalysisTask";
  TString combinedName;
  combinedName.Form("%s_histos",taskname);
  TString contName(combinedName);
  //-------------------------------------------------------
  // Init the task and do settings
  //------------------------------------------------------
  std::cout << "----#About to start AliAnalysisTaskEMCALPi0GammaCorr task " << std::endl;
  AliAnalysisTaskEMCALPi0GammaCorr* AnalysisTask = new AliAnalysisTaskEMCALPi0GammaCorr(kTRUE);
  //AnalysisTask->SetPeriod(period.Data());
  AnalysisTask->SetMC(isMC);

  AliEmcalMCTrackSelector *mcPartTask = NULL;
  if(isMC){ 
      std::cout << "AddTask MC Track Selector " << std::endl;
      AddTaskMCTrackSelector("mcparticles", kFALSE, kFALSE, 1, kFALSE); //this is needed only in mC
  }
  std::cout << "----Adding cluster container, track container" << std::endl;
  AnalysisTask->AddClusterContainer("usedefault");
  
  AliTrackContainer* trackCont = AnalysisTask->AddTrackContainer("usedefault");
  if(!trackCont) ::Error("AddTask", "Track container not found" ) ;
  trackCont->SetName("ForCorrelation");
  trackCont->SetFilterHybridTracks(kTRUE); //gives me Hyprid tracks
  const char* periodstr = "LHC13d";
 
  AliTrackContainer* trackContMatching = AnalysisTask->AddTrackContainer("usedefault");
  if(!trackContMatching) ::Error("AddTask", "Track container not found" ) ;
  trackContMatching->SetName("ForMatching");
  trackContMatching->SetTrackFilterType(AliEmcalTrackSelection::kTPCOnlyTracks);  


  AliTrackContainer* trackITSOnly  = AnalysisTask->AddTrackContainer("usedefault");
  trackITSOnly->SetName("ITSOnly");
  trackITSOnly->SetTrackFilterType(AliEmcalTrackSelection::kCustomTrackFilter);
  AliESDtrackCuts* myCuts = AliESDtrackCuts::GetStandardITSPureSATrackCuts2010();
  trackITSOnly->AddTrackCuts(myCuts);
  
  if(isMC){
    std::cout << "Adding MC particle container \n" << std::endl;
    AliMCParticleContainer *mcpcont = AnalysisTask->AddMCParticleContainer("mcparticles");
    if(!mcpcont) ::Error("AddTask", "MC Particle container not found");
    mcpcont->SetName("mcparticles");
  }
  else 
  {
    std::cout << " Not adding MC particle container " << std::endl;
  }

  //-------------------------------------------------------
  // Add some selection criteria
  //-------------------------------------------------------
  if(!isMC) AnalysisTask->SetOffTrigger(evtTriggerType|evtMixingType); //..select only evets of type evtTriggerType and evtMixingType
    //..for Run1 pPb
  //AnalysisTask->SetUseManualEvtCuts(kTRUE);
  AnalysisTask->SetUseAliAnaUtils(kTRUE); //this does automatically some vertex selection and pileup suppression
  AnalysisTask->SetVzRange(-10,10);
  if(!isMC) AnalysisTask->SetCentRange(0,100.0);
  //..new task for run2
  //AnalysisTask->SetUseNewCentralityEstimation(kFALSE); //maybe this is what is required
  
  AnalysisTask->GetTrackContainer("ForMatching")->SetParticlePtCut(0.20);
  AnalysisTask->GetTrackContainer("ForCorrelation")->SetParticlePtCut(1.0);

  if(AnalysisTask->GetClusterContainer(0))
  {
    std::cout << "Setting cuts for clusters" << std::endl;
    AnalysisTask->GetClusterContainer(0)->SetClusECut(0);                
    AnalysisTask->GetClusterContainer(0)->SetClusPtCut(0.30);       
    AnalysisTask->GetClusterContainer(0)->SetClusUserDefEnergyCut(AliVCluster::kHadCorr,0);
    AnalysisTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }
  AnalysisTask->SetNeedEmcalGeom(kTRUE);
  AnalysisTask->SetSavePool(kFALSE);
  AnalysisTask->SetEvtTriggerType(evtTriggerType);   
  AnalysisTask->SetEvtMixType(evtMixingType);       

  mgr->AddTask(AnalysisTask);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(),TList::Class(),
	  	  	  	  	  	  	  	  	  	  	  	  	   AliAnalysisManager::kOutputContainer,
		  	  	  	  	  	  	  	  	  	  	  	  	   Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (AnalysisTask, 0,  cinput1 );
  mgr->ConnectOutput (AnalysisTask, 1, coutput1 );
 
  return AnalysisTask;
}
