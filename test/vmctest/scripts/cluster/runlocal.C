void runlocal() {
  TStopwatch timer;
  timer.Start();


  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include"); 
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  
  //Enable the needed package
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

 // Create chain of input files
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("file.txt", 3);
 
  //ANALYSIS PART
  gROOT->LoadMacro("AliAnalysisTaskCluster.cxx++g");

    
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  
  //AliVEventHandler* esdH = new AliESDInputHandler;
  AliESDInputHandler* esdH = new AliESDInputHandler;
  esdH->SetReadFriends(0);
  mgr->SetInputEventHandler(esdH);  

  AliMCEventHandler* mcH = new AliMCEventHandler();
  mcH->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(mcH);


 //____________________________________________//
  // event selection 

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
  physSelTask->GetPhysicsSelection()->SetAnalyzeMC(); 

  AliPhysicsSelection* physSel = physSelTask->GetPhysicsSelection();
  physSel->AddBackgroundIdentification(new AliBackgroundSelection());

  
  //____________________________________________//
  // analysis task and esd track cuts
  AliAnalysisTaskCluster *task1 = new AliAnalysisTaskCluster("AliAnalysisTaskCluster");

  AliESDtrackCuts* esdTrackCutsL1 = new AliESDtrackCuts("AliESDtrackCuts","test");
  esdTrackCutsL1->SetMaxDCAToVertexXY(3.);
  esdTrackCutsL1->SetMaxDCAToVertexZ(3.);
  esdTrackCutsL1->SetAcceptKinkDaughters(kFALSE);
  

  task1->SetCuts(esdTrackCutsL1);
  task1->SelectCollisionCandidates();


  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = 
    mgr->CreateContainer("cchain",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("chist1",TList::Class(),AliAnalysisManager::kOutputContainer,
			 "Cluster.local.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task1,1,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
  //mgr->StartAnalysis("local",chain, 100,200);//startevent, nevents

  timer.Stop();
  timer.Print();
}

