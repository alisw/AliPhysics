AliAnalysisTaskPerformanceStrangeAOD ** AddTaskLambdaK0PbPbAOD(const char * outfilename, Int_t &nbin, Int_t binMin, Int_t binMax, Int_t iMCAnalysis = 0) {


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
  
  if (inputDataType != "AOD") {
    Printf("ERROR! This task can only run on AODs!");
  }

  // Configure analysis
  //===========================================================================
  // Int_t nbMinTPCclusters = 80;
  // Int_t lCollidingSystems = 1; 
  // TString fAnalysisType = "ESD";
    //    TString lAnalysisPidMode  = "withPID";
    // TString lAnalysisCut      = "no";    
    //Int_t iMCAnalysis = 0;
     
     nbin = binMax - binMin + 1;

  AliAnalysisTaskPerformanceStrangeAOD ** task = new AliAnalysisTaskPerformanceStrangeAOD*[nbin];
  Int_t itask = -1;
  for(Int_t ibin = binMin; ibin <= binMax; ibin++){  
    itask++;

    task[itask] = new AliAnalysisTaskPerformanceStrangeAOD("TaskLambdaK0");
    cout << "Booking " << ibin << "  "<< itask << " " << task[itask] <<endl;
   
    mgr->AddTask(task[itask]);
 
   if(ibin == 0){
      task[itask]->SetCentMin(0);
      task[itask]->SetCentMax(5);
    }
    if(ibin == 1){
      task[itask]->SetCentMin(5);
      task[itask]->SetCentMax(10);
    }
    if(ibin == 2){
      task[itask]->SetCentMin(10);
      task[itask]->SetCentMax(20);
    }
    if(ibin == 3){
      task[itask]->SetCentMin(20);
      task[itask]->SetCentMax(30);
    }
    if(ibin == 4){
      task[itask]->SetCentMin(30);
      task[itask]->SetCentMax(40);
    }
    if(ibin == 5){
      task[itask]->SetCentMin(40);
      task[itask]->SetCentMax(50);
    }
    if(ibin == 6){
      task[itask]->SetCentMin(50);
      task[itask]->SetCentMax(60);
    }
    if(ibin == 7){
      task[itask]->SetCentMin(60);
      task[itask]->SetCentMax(70);
    }
    if(ibin == 8){
      task[itask]->SetCentMin(70);
      task[itask]->SetCentMax(80);
    }
    if(ibin == 9){
      task[itask]->SetCentMin(80);
      task[itask]->SetCentMax(90);
    }
    if(ibin == 10){
      task[itask]->SetCentMin(0);
      task[itask]->SetCentMax(90);
    }
    TString outfilenameCentr = outfilename;
    outfilenameCentr.ReplaceAll(".root",Form("_%2.2d.root",ibin));

    AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("clambdak0Histo_%2.2d",ibin), TList::Class(),AliAnalysisManager::kOutputContainer, outfilenameCentr.Data());
    //    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("clambdak0Centr_%2.2d",ibin), AliAnalysisCentralitySelector::Class(),AliAnalysisManager::kOutputContainer, outfilenameCentr.Data());
    //    AliAnalysisDataContainer *output_cuts = mgr->CreateContainer(Form("cuts_%2.2d",ibin), AliESDtrackCuts::Class(), AliAnalysisManager::kOutputContainer, outfilenameCentr.Data()); 
    mgr->ConnectInput (task[itask], 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[itask],1,coutput1);
    //  mgr->ConnectOutput(task[itask],2,coutput2);
    // mgr->ConnectOutput(task[itask],3,output_cuts);
    

  }
  // TODO:
  // IO into folders in a file?

  // Set I/O

  return task;
}   


