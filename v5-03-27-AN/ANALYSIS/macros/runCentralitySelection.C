void runCentralitySelection(const char *mode="local")
{
  // Load common libraries
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
  gSystem->Load("libCORRFW");   
  gSystem->Load("libPWGHFbase");   
  gSystem->Load("libPWGmuon");   
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");

  // filename 
  TChain* chain = new TChain("esdTree");
  chain->AddFile("/home/alberica/analysis/centrality/data/alice/sim/LHC10a12/104157/998/root_archive.zip#AliESDs.root");
 
 // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisManager");
  mgr->SetDebugLevel(10);
  
  // My task
  gROOT->LoadMacro("AliCentralitySelectionTask.cxx++g");   
  AliCentralitySelectionTask *task = new AliCentralitySelectionTask("CentralitySelection"); 
  task->SetPercentileFile("test_AliCentralityBy1D.root");
  task->SetCentralityMethod("V0");
  mgr->AddTask(task);

  // My dummy task
  gROOT->LoadMacro("AliDummy.cxx++g");   
  AliDummy *dummytask = new AliDummy("Dummy"); 
  mgr->AddTask(dummytask);



  AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);

  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  // Create containers for input/output
  mgr->ConnectInput (task,0, mgr->GetCommonInputContainer());
  mgr->ConnectInput (dummytask,0, mgr->GetCommonInputContainer());
  //  mgr->ConnectOutput(task,0, mgr->GetCommonOutputContainer());  

  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if(!mgr->InitAnalysis()) return;

  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis(mode, chain);

};
