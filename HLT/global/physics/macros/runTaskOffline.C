void runTaskOffline()
{
  // load analysis framework
  gSystem->Load("libANALYSISalice");

  gROOT->LoadMacro("$ALICE_ROOT/PWGUD/macros/CreateESDChain.C");

  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$(ALICE_ROOT)/include
  // or in each macro
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gROOT->Macro("$ALICE_ROOT/PWGPP/CalibMacros/CPass0/LoadLibraries.C");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  // Add ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kTRUE);
  esdH->SetActiveBranches("ESDfriend");

  // Register input handler to manager
  mgr->SetInputEventHandler(esdH);

  // Create task

  //gROOT->LoadMacro("AliAnalysisTaskPt.cxx+g");
  AliAnalysisTaskPt *task = new AliAnalysisTaskPt("TaskPt");
  task->SetUseFriends(kTRUE);
  
  // Add task
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("chist", TList::Class(),    AliAnalysisManager::kOutputContainer, "Pt.ESD.1.root");

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 0, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  TChain *chain1 = new TChain("esdTree");
  chain1->Add("AliESDs.root");

  mgr->StartAnalysis("local", chain1);
}
