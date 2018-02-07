
//const Bool_t anaType   = 1;//0 HD; 1 UU;
const Bool_t OnlineSearchV0 = kFALSE; // True=online, false=offline
const TString whichV0 = "Lambda"; //Lambda or Kaon
const TString whichV0region = "signal";//signal or sideband region
const TString whichMixing = "Tracklets"; //V0M or #Tracklets
//const int whichfilterbit = 96;
const int whichfilterbit = 128;

AliAnalysisTaskPLFemto *AliAnalysisTaskPLFemto(Int_t system=0/*0=pp,1=PbPb*/,
					       
					       Bool_t theMCon=kFALSE)
  
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) 
    {
      ::Error("AddTaskPLFemtoSpectra", "No analysis manager to connect to.");
      return NULL;
    }  
  
  //AliVEventHandler *inputEvents = 0;
  
  // Event handlers:
  /*
  if(!theMCon)
    {
      inputEvents = new AliAODInputHandler();
      mgr->SetInputEventHandler(inputEvents);
    }
  else //Monte Carlo
    {
      inputEvents = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(inputEvents);
    }
  */
  //CREATE THE TASK

  printf("CREATE TASK\n");
  // create the task

  AliAnalysisTaskPLFemto *task = new AliAnalysisTaskPLFemto("AliAnalysisPLFemto",OnlineSearchV0,whichV0,whichV0region,whichfilterbit);//calls the constructor of the class
  task->SelectCollisionCandidates(AliVEvent::kMB);
  //task->SelectCollisionCandidates(AliVEvent::kINT7+AliVEvent::kMB);
  //task->SelectCollisionCandidates(AliVEvent::kINT7+AliVEvent::kMB+AliVEvent::kINT8);
  //AliAnalysisTaskPLFemto *task = new AliAnalysisTaskPLFemto("AliAnalysisPLFemto");//calls the constructor of the class

  task->SetMC(theMCon);
  task->SetDebugLevel(0);
  mgr->AddTask(task);



  // Create and connect containers for input/output
  
  //usercomment = "_" + usercomment;  

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
  //TString output5 = "MC_dir";//directory for purely Monte Carlo output (doesn't work at the moment)
  //output5 += usercomment;

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



