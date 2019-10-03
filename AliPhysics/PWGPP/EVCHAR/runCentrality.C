void runCentrality(const char * type = "a", const char *mode="grid")
{
  // Load common libraries
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  // form filename 
  TString filenameStr = Form("%s",type);
  filenameStr = TString("LHC10g2")+filenameStr;
  const char * filename = filenameStr.Data();
  
  // Create and configure the alien handler plugin
  gROOT->LoadMacro("CreateAlienHandlerAOD.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandlerAOD(filename);  
  if(!alienHandler) return;
 
 // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisManager");

  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
  
  // My task
  gROOT->LoadMacro("AliAnalysisTaskCentrality.cxx++g");   
  AliAnalysisTaskCentrality *task = new AliAnalysisTaskCentrality("CentralityTask");  
  // Writing (or not) output tree
  //task->SetTreeFilling(writeTree);
  task->SetMCInput();
  mgr->AddTask(task);

  AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);

  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);
  
  // Physics selection
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  //AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(kFALSE);
  // Selecting collision candidates
  //task->SelectCollisionCandidates();

  // Create containers for input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutput1",TList::Class(),
                                            AliAnalysisManager::kOutputContainer, 
  					    "cenHistos.root");
  mgr->ConnectOutput(task, 1, coutput1);
  
//   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("coutput2",TTree::Class(),
//   					    AliAnalysisManager::kOutputContainer,
//         				    "cenTree.root");
// //   coutput2->SetSpecialOutput();  
//   mgr->ConnectOutput(task, 2, coutput2);

  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if(!mgr->InitAnalysis()) return;

  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis(mode);

};
