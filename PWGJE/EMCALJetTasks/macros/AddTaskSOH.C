// $Id$

AliAnalysisTaskSOH* AddTaskSOH(const char *name, Bool_t ScaleFactorP=kFALSE, Bool_t ClusterP=kFALSE, Double_t Zvtx=10)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskSOH", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskSOH", "This task requires an input event handler");
    return NULL;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliAnalysisTaskSOH *taskSOH = new AliAnalysisTaskSOH(Form("QA_soh_%s",name));

  AliESDtrackCuts *esdTrackCuts = 0x0;
  AliESDtrackCuts *hybridTrackCuts1 = 0x0;
  AliESDtrackCuts *hybridTrackCuts2 = 0x0;
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");

  esdTrackCuts = CreateTrackCutsPWGJE(10001006);
  hybridTrackCuts1 = CreateTrackCutsPWGJE(1006);
  hybridTrackCuts2 = CreateTrackCutsPWGJE(10041006);
  
  taskSOH->SetEsdTrackCuts(esdTrackCuts);
  taskSOH->SetHybridTrackCuts1(hybridTrackCuts1);
  taskSOH->SetHybridTrackCuts2(hybridTrackCuts2);

  taskSOH->SetMcProcess(kTRUE);
  taskSOH->SetTrackProcess(kTRUE);
  taskSOH->SetSFProcess(ScaleFactorP);
  taskSOH->SetClusterProcess(ClusterP);
  taskSOH->SetZvtx(Zvtx);

  // Add task(s)
  mgr->AddTask(taskSOH); 

  // ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);
  
  // MC truth handler
  AliMCEventHandler* mcEvtHdl = new AliMCEventHandler();
  mcEvtHdl->SetReadTR(kTRUE);
  mgr->SetMCtruthEventHandler(mcEvtHdl); 

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutputpt = mgr->CreateContainer(name, TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults_soh.root");

  // Connect input/output
  mgr->ConnectInput(taskSOH, 0, cinput);
  mgr->ConnectOutput(taskSOH, 1, coutputpt);

  return taskSOH;
}
