// $Id$
AliAnalysisTaskEMCALPi0GammaCorr* AddTaskEMCALPi0GammaCorr(
  UInt_t      evtTriggerType         = AliVEvent::kEMCEGA, //AliVEvent::kAnyINT,// AliVEvent::kEMCEGA,//..use this type of events to combine gammas(trigger) with hadrons
  UInt_t      evtMixingType          = AliVEvent::kAnyINT,//..use only this type of events to fill your mixed event pool with tracks
  Double_t    trackptcut             = 0.15,              //..
  Double_t    clusptcut              = 0.15,              //..
  Bool_t      SavePool               = kFALSE,                 //..saves a mixed event pool to the output event
  Bool_t      isMC                   = kFALSE, 
  const char *period                 = "lhc13d", 
  const char *trackName              = "usedefault",
  const char *clusName               = "usedefault",
  const char *taskname               = "AliAnalysisTask",
)
{  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)  ::Error("AddTask", "No analysis manager to connect to.");
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) ::Error("AddTask", "This task requires an input event handler");
 
  if(trackName=="usedefault")trackName = "usedefault"; //it was Tracks "tracks"
  if(clusName =="usedefault")clusName  = "usedefault";//"caloClusters"; //CaloClusters"; //caloClusters
  //-------------------------------------------------------
  // Built the name of the Task together
  //-------------------------------------------------------
  TString combinedName;
  combinedName.Form("%s_histos",taskname);
  TString contName(combinedName);
  //-------------------------------------------------------
  // Init the task and do settings
  //------------------------------------------------------
  std::cout << "----#About to start AliAnalysisTaskEMCALPi0GammaCorr task " << std::endl;
  AliAnalysisTaskEMCALPi0GammaCorr* AnalysisTask = new AliAnalysisTaskEMCALPi0GammaCorr(kTRUE);
  //AnalysisTask->SetPeriod(period.Data());
  std::cout << "----Adding cluster container, track container" << std::endl;
  AnalysisTask->AddClusterContainer("usedefault");
  AliTrackContainer* trackCont = AnalysisTask->AddTrackContainer("Tracks");
  trackCont->SetName("ForCorrelation");
  trackCont->SetFilterHybridTracks(kTRUE); //gives me Hyprid tracks
  const char* periodstr = "LHC13d";
 
  AliTrackContainer* trackContMatching = AnalysisTask->AddTrackContainer("Tracks");
  trackContMatching->SetName("ForMatching");
  trackContMatching->SetTrackFilterType(AliEmcalTrackSelection::kTPCOnlyTracks);  

  //-------------------------------------------------------
  // Add some selection criteria
  //-------------------------------------------------------
  if(!isMC) AnalysisTask->SetOffTrigger(evtTriggerType|evtMixingType); //..select only evets of type evtTriggerType and evtMixingType
  
    //..for Run1 pPb
  AnalysisTask->SetUseManualEvtCuts(kTRUE);
  AnalysisTask->SetUseAliAnaUtils(kTRUE); //this does automatically some vertex selection and pileup suppression
  AnalysisTask->SetVzRange(-10,10);
  if(!isMC) AnalysisTask->SetCentRange(0,100.0);
  //..new task for run2
  AnalysisTask->SetUseNewCentralityEstimation(kFALSE); //maybe this is what is required
  
  AnalysisTask->GetTrackContainer("ForMatching")->SetParticlePtCut(trackptcut);
  AnalysisTask->GetTrackContainer("ForCorrelation")->SetParticlePtCut(1.0);

  if(AnalysisTask->GetClusterContainer(clusName))
  {
      std::cout << "Setting cuts for clusters" << std::endl;
	  AnalysisTask->GetClusterContainer(clusName)->SetClusECut(0);                
	  AnalysisTask->GetClusterContainer(clusName)->SetClusPtCut(clusptcut);       
	  AnalysisTask->GetClusterContainer(clusName)->SetClusUserDefEnergyCut(AliVCluster::kHadCorr,0);
	  AnalysisTask->GetClusterContainer(clusName)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  AnalysisTask->SetNeedEmcalGeom(kTRUE);
  AnalysisTask->SetSavePool(SavePool);
  AnalysisTask->SetEvtTriggerType(evtTriggerType);   
  AnalysisTask->SetEvtMixType(evtMixingType);       

  //Whether or not running on MC: 
  AnalysisTask->SetMC(isMC);
  
  mgr->AddTask(AnalysisTask);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(),TList::Class(),
	  	  	  	  	  	  	  	  	  	  	  	  	   AliAnalysisManager::kOutputContainer,
		  	  	  	  	  	  	  	  	  	  	  	  	   Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput  (AnalysisTask, 0,  cinput1 );
  mgr->ConnectOutput (AnalysisTask, 1, coutput1 );
  //mgr->ConnectOutput (AnalysisTask, 2, coutput2 );
  return AnalysisTask;
}
