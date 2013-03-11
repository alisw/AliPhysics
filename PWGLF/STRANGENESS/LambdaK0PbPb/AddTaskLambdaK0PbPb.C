AliAnalysisTaskPerformanceStrange ** AddTaskLambdaK0PbPb(const char * outfilename = "lambdak0.root", Int_t &nbin, Int_t binMin, Int_t binMax, Int_t isMCAnalysis = 0, Bool_t usePID = 1) {


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


  AliAnalysisCentralitySelector * centrSelector = new AliAnalysisCentralitySelector();
  centrSelector->SetIsMC(isMCAnalysis);
  centrSelector->SetCentralityEstimator("V0M"); 
  // Configure analysis
  //===========================================================================
     
  AliESDtrackCuts * myTracksCuts = new AliESDtrackCuts();
  myTracksCuts->SetRequireTPCRefit(kTRUE);
  //          myTracksCuts->SetRequireITSRefit(kTRUE);
  myTracksCuts->SetRequireITSRefit(kFALSE);
  // myTracksCuts->SetMinNClustersTPC(nbMinTPCclusters);
  myTracksCuts->SetMinNCrossedRowsTPC(70);
  myTracksCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
 
  nbin = binMax - binMin + 1;

  AliAnalysisTaskPerformanceStrange ** task = new AliAnalysisTaskPerformanceStrange*[nbin];
  Int_t itask = -1;
  for(Int_t ibin = binMin; ibin <= binMax; ibin++){  
    itask++;

    task[itask] = new AliAnalysisTaskPerformanceStrange("TaskLambdaK0");
    cout << "Booking " << ibin << "  "<< itask << " " << task[itask] <<endl;

    task[ibin]->SetAnalysisType("ESD");
    cout << "1" << endl;
    task[ibin]->SetAnalysisMC(isMCAnalysis); // 0 or 1
    cout << "2" << endl;
    task[ibin]->SetCollidingSystems(2); // 0 =pp, 1=AA  2=pA
    cout << "3" << endl;
    task[ibin]->SetAnalysisCut("no");
    cout << "4" << endl;
    task[ibin]->SetQASelector(kFALSE);  // Todo -> put trees for QA
    cout<< "5" << endl;
    if(usePID) 
      task[ibin]->SetUsePID("withPID"); // withPID or withoutPID
    else
      task[ibin]->SetUsePID("withoutPID"); // withPID or withoutPID
    cout << "5" << endl;
    task[itask]->SetTrackCuts(myTracksCuts);
   
    mgr->AddTask(task[itask]);
  
    // // Set Cuts
    // if (!esdTrackCuts)
    //   {
    //     printf("ERROR: esdTrackCuts could not be created\n");
    //     return;
    //   }  
    // task->SetTrackCuts(esdTrackCuts);
    
    // set centrality
    AliAnalysisCentralitySelector * centrBin = (AliAnalysisCentralitySelector*) centrSelector->Clone();
    //centrBin->SetCentralityBin(ibin);
    if(ibin == 0){
      centrBin->SetCentralityBin(0,5);
      task[itask]->SetCentralitySelector(centrBin);
    }
    if(ibin == 1){
      centrBin->SetCentralityBin(5,10);
      task[itask]->SetCentralitySelector(centrBin);
    }
    if(ibin == 2){
      centrBin->SetCentralityBin(10,20);
      task[itask]->SetCentralitySelector(centrBin);
    }
    if(ibin == 3){
      centrBin->SetCentralityBin(20,30);
      task[itask]->SetCentralitySelector(centrBin);
    }
    if(ibin == 4){
      centrBin->SetCentralityBin(30,40);
      task[itask]->SetCentralitySelector(centrBin);
    }
    if(ibin == 5){
      centrBin->SetCentralityBin(40,50);
      task[itask]->SetCentralitySelector(centrBin);
    }
    if(ibin == 6){
      centrBin->SetCentralityBin(50,60);
      task[itask]->SetCentralitySelector(centrBin);
    }
    if(ibin == 7){
      centrBin->SetCentralityBin(60,70);
      task[itask]->SetCentralitySelector(centrBin);
    }
    if(ibin == 8){
      centrBin->SetCentralityBin(70,80);
      task[itask]->SetCentralitySelector(centrBin);
    }
    if(ibin == 9){
      centrBin->SetCentralityBin(80,90);
      task[itask]->SetCentralitySelector(centrBin);
    }
    if(ibin == 10){
      centrBin->SetCentralityBin(0,90);
      task[itask]->SetCentralitySelector(centrBin);
    }
    TString outfilenameCentr = outfilename;
    outfilenameCentr.ReplaceAll(".root",Form("_%2.2d.root",ibin));

    AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("clambdak0Histo_%2.2d",ibin), TList::Class(),AliAnalysisManager::kOutputContainer, outfilenameCentr.Data());
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("clambdak0Centr_%2.2d",ibin), AliAnalysisCentralitySelector::Class(),AliAnalysisManager::kOutputContainer, outfilenameCentr.Data());
    AliAnalysisDataContainer *output_cuts = mgr->CreateContainer(Form("cuts_%2.2d",ibin), AliESDtrackCuts::Class(), AliAnalysisManager::kOutputContainer, outfilenameCentr.Data()); 
    mgr->ConnectInput (task[itask], 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[itask],1,coutput1);
    mgr->ConnectOutput(task[itask],2,coutput2);
    mgr->ConnectOutput(task[itask],3,output_cuts);
    

  }
  // TODO:
  // IO into folders in a file?

  // Set I/O

  return task;
}   


