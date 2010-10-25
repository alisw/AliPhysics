AliAnalysisTaskMinijet* AddTaskMinijet(TString format="esd",Bool_t useMC = kFALSE, TString kGridDataSet="LHC10e")
{

  //starting with periode LHC10e, there are also events triggered with High Mult trigger
  //add stand alone task in case these events are availale
  Bool_t IsHighMult= (!kGridDataSet.CompareTo("LHC10e") || !kGridDataSet.CompareTo("LHC10f"));
  //extent with new periods (LHC10g,h,...)

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) mgr = new AliAnalysisManager("Analysis train");

  // Set cuts (used for ESDs) 
  //===========================================================================
  AliESDtrackCuts* esdTrackCutsITSTPC=0x0;
  if(!format.CompareTo("esd")){
    esdTrackCutsITSTPC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
  }

  // Configure tasks
  //===========================================================================
   
  //first task for min bias events
  AliAnalysisTaskMinijet *task = new AliAnalysisTaskMinijet("AliAnalysisTaskMinijet Min bias");
  task->UseMC(useMC);
  task->SetRadiusCut(0.7);
  task->SetTriggerPtCut(0.7);
  task->SetAssociatePtCut(0.4);
  task->SetEventAxis(1); //1=random track
  task->SetDebugLevel(0);
  task->SetMaxVertexZ(10.);
  if(!format.CompareTo("esd")){
    task->SetCuts(esdTrackCutsITSTPC);
    task->SetMode(0);//0=reading ESDs
    task->SelectCollisionCandidates(AliVEvent::kMB);//MB
  }

  //second task for high multipliciy events
  AliAnalysisTaskMinijet *taskHM =0x0;
  if(IsHighMult && !format.CompareTo("esd")){
    taskHM  = new AliAnalysisTaskMinijet("AliAnalysisTaskMinijet HighMult");
    taskHM->UseMC(useMC);
    taskHM->SetRadiusCut(0.7);
    taskHM->SetTriggerPtCut(0.7);
    taskHM->SetAssociatePtCut(0.4);
    taskHM->SetEventAxis(1); //1=random track
    taskHM->SetDebugLevel(0);
    taskHM->SetMaxVertexZ(10.);
    taskHM->SetCuts(esdTrackCutsITSTPC);
    taskHM->SetMode(0);//0=reading ESDs
    taskHM->SelectCollisionCandidates(AliVEvent::kHighMult);//high mult triggered event
  }

  if (!format.CompareTo("aod")){
    task->SetMode(1);// 1 = reading AODs
  }

  //create output container and add task to mgr
  //===========================================================================
  AliAnalysisDataContainer *output1 = 
    mgr->CreateContainer("chist", TList::Class(), 
			 AliAnalysisManager::kOutputContainer, 
			 "EventAxis.root");
  // add task to the manager
  mgr->AddTask(task);
  
  //connect input and output
  mgr->ConnectInput(task, 0,  mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, output1);


  //===========================================================================
  //do same for high mult task if necassary
  if(!format.CompareTo("esd") && IsHighMult){
    //create output container
    AliAnalysisDataContainer *outputHM = 0x0;
    outputHM =  mgr->CreateContainer("chist_HM", TList::Class(), 
				     AliAnalysisManager::kOutputContainer, 
				     "EventAxis.root");
    // add task to the manager
    mgr->AddTask(taskHM);
    
    //connect input and output
    mgr->ConnectInput (taskHM, 0,  mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskHM, 1, outputHM);
  }
  
  return task;

}

