void RunZDCTowerCalib(char *path="/alice/cern.ch/user/d/defalco/data", 
		    Bool_t isLocal=kFALSE, 
		    Double_t adcmin=7){
  // load analysis framework
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice") ;

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  AliAnalysisGrid *alienHandler = 0;    // only used in GRID analysis   
  TChain* chain = 0;
  if (!isLocal) {
    gROOT->LoadMacro("CreateAlienHandler.C"); // Create and configure the alien handler plugin
    alienHandler = CreateAlienHandler(path);  
    if (!alienHandler) return;
  }
  else { 
    char filename[500]; 
    sprintf (filename,"%s/AliESDs.root",path); 
    chain = new TChain("esdTree"); 
    chain->Add(filename); 
    chain->Add(filename); 
  }

  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/ZDC");

  // Create the analysis manager

  AliAnalysisManager *mgr = new AliAnalysisManager("ZNCalibManager");
  if (!isLocal) mgr->SetGridHandler(alienHandler);

  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  // Create task

  gROOT->LoadMacro("AliZDCTowerCalibTask.cxx+g");
  AliAnalysisTask *task = new AliZDCTowerCalibTask("TaskZDCCalib");
  ((AliZDCTowerCalibTask*) task)->SetADCMin(adcmin); 

  // Add task
  mgr->AddTask(task);

  // Create containers for input/output

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
//   AliAnalysisDataContainer *coutput = mgr->CreateContainer("chist", TH1::Class(),    AliAnalysisManager::kOutputContainer, "Pt.ESD.1.root");

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
//   mgr->ConnectOutput(task, 0, coutput);

  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  if (isLocal) mgr->StartAnalysis("local", chain);
  else mgr->StartAnalysis("grid");

}
