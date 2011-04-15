/*
  Template of calibration/filtering  macro using ESD

  Example:
  .L $ALICE_ROOT/ANALYSIS/macros/runCalibTrain.C
  runCalibTrain(105160);

*/

void runCalibTrain(TString runNumberString, const char *inFileName = "AliESDs.root")
{
  gROOT->Macro("LoadLibraries.C");
  gROOT->LoadMacro("ConfigCalibTrain.C");
  gROOT->LoadMacro("AddTaskFilterFriend.C");
  gROOT->LoadMacro("AddTaskFilterFriendSecond.C");
  gROOT->LoadMacro("AddTaskAddObject.C");

  // detector tasks
  
  gROOT->LoadMacro("AddTaskTPCCalib.C");
  
  AliLog::SetClassDebugLevel("AliESDEvent",19);
  TChain *chain = new TChain("esdTree");
  
  // Steering input chain
  
  chain->Add(inFileName);
  Int_t runNumber = runNumberString.Atoi();
  printf("runNumber from runCalibTrain = %d\n",runNumber);
  ConfigCalibTrain(runNumber, "raw://");
  
  AliAnalysisManager *mgr  = new AliAnalysisManager("ESD to ESD", "Analysis Manager");
  // mgr->SetDebugLevel(3);
  
  // Input
  
  AliESDInputHandler* inpHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler  (inpHandler);
  
  // Output
  
  AliESDHandler* esdHandler   = new AliESDHandler();
  mgr->SetOutputEventHandler(esdHandler);
  esdHandler->SetOutputFileName("AliESDfriends_v1.root");
  // Steering Tasks
  
  AliAnalysisTaskFilterFriend* filter = AddTaskFilterFriend();
  AliAnalysisTaskFilterFriendSecond* filter2 = AddTaskFilterFriendSecond();
  AliAnalysisTaskAddObject* add = AddTaskAddObject();
  
  // Detector Tasks
  
  AliAnalysisTask* tTPC = AddTaskTPCCalib(runNumber);
  
  // Run the analysis
  
  if (!mgr->InitAnalysis()) {
    printf("Analysis cannot be started, returning\n");
    return;
  }
  
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain);
  
  return;
}
