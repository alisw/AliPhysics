void anaTask()
{
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB"); 
  AliCDBManager::Instance()->SetSpecificStorage("PHOS/*","local://./");
  
  // load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice"); //AliAnalysisTaskSE
  gSystem->Load("libPWG4CaloCalib");
  
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("list.txt", 15);
  
  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$ALICE_ROOT/include
  // or in each macro
  // gSystem->AddIncludePath("$ALICE_ROOT/include");
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("PHOSPi0Calib");
  
  //Input event handler
  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);
  
  //Output event handler
  AliAODHandler* aodoutHandler   = new AliAODHandler();
  aodoutHandler->SetOutputFileName("aod.root");
  mgr->SetOutputEventHandler(aodoutHandler);
  
  // ESD filter task
  gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/AddTaskESDfilter.C");
  AliAnalysisTaskESDfilter *esdfilter = AddTaskESDfilter(kFALSE);
  
  // Calibration task 
  //
  AliAnalysisTaskPHOSPi0CalibSelection *task = new AliAnalysisTaskPHOSPi0CalibSelection("PHOSPi0CalibSelection");
  task->SetClusterMinEnergy(0.4); 
  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histos", TList::Class(),  AliAnalysisManager::kOutputContainer, "PHOShistos.root");
  
  // Connect input/output
  mgr->ConnectInput   (task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput2);
  
  // Enable debug printouts
  //mgr->SetDebugLevel(10);
  
  if (!mgr->InitAnalysis())
    return;
  
  mgr->PrintStatus();
  
  mgr->StartAnalysis("local", chain);
}
