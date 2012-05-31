void runAcordeAna()
{
  // load analysis framework

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("files.txt",1);
  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$(ALICE_ROOT)/include
  // or in each macro

  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);

  // Create task

  gROOT->LoadMacro("AliAnalysisTaskAcorde.cxx+g");
  AliAnalysisTask *task = new AliAnalysisTaskAcorde("TaskAcordeTest");

  // Add task
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("treeCosmic",TTree::Class(),AliAnalysisManager::kOutputContainer,"acordeOutput1.root");
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("treeList",TList::Class(),AliAnalysisManager::kOutputContainer,"acordeOutput2.root");

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 0, coutput);
  mgr->ConnectOutput(task,1,coutput1);
  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  mgr->StartAnalysis("local", chain);
}
