AliAnalysisTaskLK0Spectra * AddTaskLK0Spectra(const char * outfilename, Int_t ibin, Int_t iMCAnalysis = 0,  Bool_t usePID = kTRUE) {


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPhysicsSelection", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPhysicsSelection", "This task requires an input event handler");
    return NULL;
  }
	
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  if (inputDataType != "ESD") {
    Printf("ERROR! This task can only run on ESDs!");
  }

  // Configure analysis
  //===========================================================================
  // Int_t nbMinTPCclusters = 80;
  // Int_t lCollidingSystems = 1; 
  // TString fAnalysisType = "ESD";
    //    TString lAnalysisPidMode  = "withPID";
    // TString lAnalysisCut      = "no";    
    //Int_t iMCAnalysis = 0;
     
    AliESDtrackCuts * myTracksCuts = new AliESDtrackCuts();
     myTracksCuts->SetRequireTPCRefit(kTRUE);
     //myTracksCuts->SetRequireITSRefit(kTRUE);
	 myTracksCuts->SetRequireITSRefit(kFALSE);
     // myTracksCuts->SetMinNClustersTPC(nbMinTPCclusters);
     myTracksCuts->SetMinNCrossedRowsTPC(70);
     myTracksCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
	
	char taskName[15];
	sprintf(taskName,"TaskLambdaK0_%d",ibin);
	
	
  AliAnalysisTaskLK0Spectra * task = new AliAnalysisTaskLK0Spectra(taskName);
 
	task->SetAnalysisType("ESD");
    task->SetAnalysisMC(iMCAnalysis); // 0 or 1
    task->SetCollidingSystems(1); // 0 =pp, 1=AA
    task->SetAnalysisCut("no");
    task->SetQASelector(kFALSE);
    if(usePID) 
		task->SetUsePID("withPID"); // withPID or withoutPID
    else
		task->SetUsePID("withoutPID"); // withPID or withoutPID
	task->SetArmenterosCut(0.2);
    task->SetTrackCuts(myTracksCuts);
   
    mgr->AddTask(task);
      
    // set centrality
	
	AliAnalysisCentralitySelector * centrBin = new AliAnalysisCentralitySelector();
	centrBin->SetIsMC(iMCAnalysis);
	centrBin->SetCentralityEstimator("V0M"); // Todo: add parameter to macro?
	
	if(ibin == 0){
    centrBin->SetCentralityBin(0,5);
    task->SetCentralitySelector(centrBin);
    }
	
    if(ibin == 1){
    centrBin->SetCentralityBin(5,10);
    task->SetCentralitySelector(centrBin);
    }
	
    if(ibin == 2){
    centrBin->SetCentralityBin(10,20);
    task->SetCentralitySelector(centrBin);
    }
	
    if(ibin == 3){
    centrBin->SetCentralityBin(20,30);
    task->SetCentralitySelector(centrBin);
    }
	
    if(ibin == 4){
    centrBin->SetCentralityBin(30,40);
    task->SetCentralitySelector(centrBin);
    }
	
    if(ibin == 5){
    centrBin->SetCentralityBin(40,50);
    task->SetCentralitySelector(centrBin);
    }
	
    if(ibin == 6){
    centrBin->SetCentralityBin(50,60);
    task->SetCentralitySelector(centrBin);
    }
	
    if(ibin == 7){
    centrBin->SetCentralityBin(60,70);
    task->SetCentralitySelector(centrBin);
    }
	
    if(ibin == 8){
    centrBin->SetCentralityBin(70,80);
    task->SetCentralitySelector(centrBin);
    }
	
    if(ibin == 9){
    centrBin->SetCentralityBin(80,90);
    task->SetCentralitySelector(centrBin);
    }
	
    if(ibin == 10){
    centrBin->SetCentralityBin(0,90);
    task->SetCentralitySelector(centrBin);
    }
	
    TString outfilenameCentr = outfilename;
    outfilenameCentr.ReplaceAll(".root",Form("_%2.2d.root",ibin));

    AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
   
	sprintf(taskName,"clambdak0Histo_%2.2d",ibin);
	
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(taskName, TList::Class(),AliAnalysisManager::kOutputContainer, Form("%s:lambdak0", AliAnalysisManager::GetCommonFileName()));
    //AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("clambdak0Centr_%2.2d",ibin), AliAnalysisCentralitySelector::Class(),AliAnalysisManager::kOutputContainer, outfilenameCentr.Data());
    //AliAnalysisDataContainer *output_cuts = mgr->CreateContainer(Form("cuts_%2.2d",ibin), AliESDtrackCuts::Class(), AliAnalysisManager::kOutputContainer, outfilenameCentr.Data()); 
	
	mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,coutput1);
    //mgr->ConnectOutput(task,2,coutput2);
    //mgr->ConnectOutput(task,3,output_cuts);
    

  
  // TODO:
  // IO into folders in a file?

  // Set I/O

  return task;
}   


