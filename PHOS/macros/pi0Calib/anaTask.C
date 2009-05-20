void anaTask()
{

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("PHOS/*","local://./");

  // load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice"); //AliAnalysisTaskSE

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("list.txt", 300);

  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$ALICE_ROOT/include
  // or in each macro
  // gSystem->AddIncludePath("$ALICE_ROOT/include");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0Calib");



  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  // Create task

  gROOT->LoadMacro("AliAnalysisTaskPi0CalibSelection.cxx+g");
  AliAnalysisTaskPi0CalibSelection *task = new AliAnalysisTaskPi0CalibSelection("Pi0CalibSelection");
  //task->SetClusterMinEnergy(0.4); 

  // Add task
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("histos", TList::Class(),  AliAnalysisManager::kOutputContainer, "histos.root");

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(10);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  mgr->StartAnalysis("local", chain);
}
