
//const Bool_t anaType   = 1;//0 HD; 1 UU;
//const Bool_t OnlineSearchV0 = kFALSE; // True=online, false=offline
//const TString whichV0 = "Lambda"; //Lambda or Kaon
//const TString whichV0region = "signal";//signal or sideband region
//const TString whichMixing = "Tracklets"; //V0M or #Tracklets
//const int whichfilterbit = 96;
//const int whichfilterbit = 128;

AliAnalysisTaskPLFemto *AliAnalysisTaskPLFemto(Int_t system=0/*0=pp,1=PbPb*/,
					       Bool_t theMCon=kFALSE,Bool_t OnlineSearchV0=kFALSE,
					       TString whichV0="Lambda",TString whichV0region="signal")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) 
    {
      ::Error("AddTaskPLFemtoSpectra", "No analysis manager to connect to.");
      return NULL;
    }  


  //Specify cut values for systematic checks:
  //Default
  AliFemtoCutValues::systematics DefaultValues = AliFemtoCutValues::kDefault;
  //Pt threshold for protons varied by 20% down:
  AliFemtoCutValues::systematics ProtonPtValueDown = AliFemtoCutValues::kProtonVariationLowerPtThresholdDown;
  //Pt threshold for protons varied by 20% up:
  AliFemtoCutValues::systematics ProtonPtValueUp = AliFemtoCutValues::kProtonVariationLowerPtThresholdUp;



  printf("CREATE TASK\n");
  // create the task

  //AliAnalysisTaskPLFemto *task[2];

  AliAnalysisTaskPLFemto *task = new AliAnalysisTaskPLFemto("AliAnalysisPLFemto",OnlineSearchV0,whichV0,whichV0region);//calls the constructor of the class
  task->SelectCollisionCandidates(AliVEvent::kMB);
  //task->SelectCollisionCandidates(AliVEvent::kINT7+AliVEvent::kMB);
  //task->SelectCollisionCandidates(AliVEvent::kINT7+AliVEvent::kMB+AliVEvent::kINT8);

  task->SetSystematics(DefaultValues);
  //task->SetSystematics(PtProtonLowValue);
  //task->SetSystematics(kDefault);
  //task->SetTrackCuts(DefaultValues);
  //task->SetCuts(cutVariationPtProtonLow);
  task->SetDebugLevel(0);
  task->SetMC(theMCon);

  mgr->AddTask(task);



  // Create and connect containers for input/output

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWGCF_PLFemto";
  //outputfile += usercomment;


  // ------ input data ------
  TString input = "infemto";
  //input += usercomment;
  TString output1 = "Evtinfo";
  //output1 += usercomment;
  TString output2 = "SPdir";//directory for single particle quantities
  //output2 += usercomment;
  TString output3 = "PIDdir";//directory for PID quantities
  //output3 += usercomment;
  TString output4 = "TPdir";//directory for two particle quantities
  //output4 += usercomment;

  //AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *cinput0  =  mgr->CreateContainer(input,TChain::Class(),AliAnalysisManager::kInputContainer);

 // ----- output data -----

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(output1,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(output2,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(output3,TList::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(output4,TList::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);



  return task;
}



