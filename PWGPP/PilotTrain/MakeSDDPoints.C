

void MakeSDDPoints(){
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");

// Common packages
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  TChain *esdChain = new TChain("esdTree");
  esdChain->AddFile("/home/prino/alice/FirstPP/104044/09000104044018.10/AliESDs.root");
 
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliESDInputHandlerRP* esdH = new AliESDInputHandlerRP();
  mgr->SetInputEventHandler(esdH);

  gROOT->LoadMacro("AliAnalysisTaskSDDRP.cxx+g");
  AliAnalysisTaskSDDRP *task= new AliAnalysisTaskSDDRP();
  task->SetGeometryFile("/home/prino/alice/FirstPP/chunk1/geometry.root");
  task->SetRunNumber(104044);
  mgr->AddTask(task);
  mgr->SetDebugLevel(2);


  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutputRP",TList::Class(),AliAnalysisManager::kOutputContainer,"/home/prino/alice/test/SDDPoints.root");



  mgr->ConnectInput(task,0,cinput1);
  mgr->ConnectOutput(task,0,coutput1);


  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", esdChain);
  }   
  delete mgr;
}   



