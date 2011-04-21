AliAnalysisTaskPerformanceStrange ** AddTaskLambdaK0PbPb(const char * outfilename, AliAnalysisCentralitySelector * centr, Int_t &nbin, Int_t binMin, Int_t binMax, Int_t iMCAnalysis = 0) {


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
    Int_t nbMinTPCclusters = 80;
    Int_t lCollidingSystems = 1; 
    TString fAnalysisType = "ESD";
    TString lAnalysisPidMode  = "withoutPID";
    TString lAnalysisCut      = "no";    
    //Int_t iMCAnalysis = 0;
     
    AliESDtrackCuts * myTracksCuts = new AliESDtrackCuts();
     myTracksCuts->SetRequireTPCRefit(kTRUE);
     myTracksCuts->SetMinNClustersTPC(nbMinTPCclusters);
 
     nbin = binMax - binMin + 1;

  AliAnalysisTaskPerformanceStrange ** task = new AliAnalysisTaskPerformanceStrange*[nbin];
  Int_t itask = -1;
  for(Int_t ibin = binMin; ibin <= binMax; ibin++){  
    itask++;

    task[itask] = new AliAnalysisTaskPerformanceStrange("TaskLambdaK0");
    cout << "Booking " << ibin << "  "<< itask << " " << task[itask] <<endl;
 
    task[itask]->SetCollidingSystems(lCollidingSystems);
    task[itask]->SetAnalysisType(fAnalysisType);
    task[itask]->SetAnalysisMC(iMCAnalysis);
    task[itask]->SetAnalysisCut(lAnalysisCut);
    task[itask]->SetUsePID(lAnalysisPidMode);
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
    AliAnalysisCentralitySelector * centrBin = (AliAnalysisCentralitySelector*) centr->Clone();
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


