AliAnalysisTaskPerformanceStrange ** AddTaskLambdaK0PbPb(const char * outfilename, AliAnalysisCentralitySelector * centr, Int_t &nbin) {


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
    
  Int_t binMin = 0; // FIXME: settable? different percentiles?
  Int_t binMax = 10;
  nbin = binMax - binMin + 1;

  AliAnalysisTaskPerformanceStrange ** task = new AliAnalysisTaskPerformanceStrange*[nbin];
  Int_t itask = -1;
  for(Int_t ibin = binMin; ibin <= binMax; ibin++){  
    itask++;

    task[itask] = new AliAnalysisTaskPerformanceStrange("TaskLambdaK0");
    cout << "Booking " << ibin << "  "<< itask << " " << task[itask] <<endl;
    
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
    centrBin->SetCentralityBin(ibin);
    task[itask]->SetCentralitySelector(centrBin);

    TString outfilenameCentr = outfilename;
    outfilenameCentr.ReplaceAll(".root",Form("_%2.2d.root",ibin));

    AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("clambdak0Histo_%2.2d",ibin), TList::Class(),AliAnalysisManager::kOutputContainer, outfilenameCentr.Data());
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("clambdak0Centr_%2.2d",ibin), AliAnalysisCentralitySelector::Class(),AliAnalysisManager::kOutputContainer, outfilenameCentr.Data());

    mgr->ConnectInput (task[itask], 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[itask],1,coutput1);
    mgr->ConnectOutput(task[itask],2,coutput2);

  }
  // TODO:
  // IO into folders in a file?

  // Set I/O

  return task;
}   


