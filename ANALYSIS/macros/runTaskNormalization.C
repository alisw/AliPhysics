
void runTaskNormalization(const char * incollection,const char * filename = "LHC09b12_7TeV_0.5T_norm.root", Bool_t isMC = 1,Int_t nev =123456789) {

  
  // Load libraries
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;
  gSystem->Load("libCORRFW") ;
  gSystem->Load("libITSbase") ;
  gSystem->Load("libPWG0base") ;


  // chain
  TChain* analysisChain = 0;
  analysisChain = new TChain("esdTree");
  if (TString(incollection).Contains(".root")){
    analysisChain->Add(incollection);
  }
  else if (TString(incollection).Contains("xml")){
    TGrid::Connect("alien://");
    TAlienCollection * coll = TAlienCollection::Open (incollection);
    while(coll->Next()){
      analysisChain->Add(TString("alien://")+coll->GetLFN());
    }
  } else {
    ifstream file_collect(incollection);
    TString line;
    while (line.ReadLine(file_collect) ) {
      analysisChain->Add(line.Data());
    }
  }
  analysisChain->GetListOfFiles()->Print();


  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  //  mgr->SetDebugLevel(3);
  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler; 

  mgr->SetInputEventHandler(esdH);
	
  if(isMC) {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mc->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mc);
  }
  // assign simple task
//   gROOT->LoadMacro("AliCollisionNormalization.cxx++g");   
//   gROOT->LoadMacro("AliCollisionNormalizationTask.cxx++g");   
  //____________________________________________//
  // Physics selection
  gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,1,!isMC); // Use Physics Selection. Enable computation of BG if is not MC
  //  task->SelectCollisionCandidates(); /// This should be disabled, at least for MC: we need all the events
  physSelTask->GetPhysicsSelection()->SetBin0Callback("TaskNormalization");

  // assign simple task
  AliCollisionNormalizationTask * task = new AliCollisionNormalizationTask("TaskNormalization");
  //  task->SetMC();
  task->SetMC(isMC);
  mgr->AddTask(task);



  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();	
  mgr->ConnectInput(task,0,cinput1);


  
  // Attach output
  cOutput = mgr->CreateContainer("Norm", TList::Class(), AliAnalysisManager::kOutputContainer,filename);
  mgr->ConnectOutput(task, 1, cOutput);      
	
  if (!mgr->InitAnalysis()) return;
	
  mgr->PrintStatus();
  mgr->StartAnalysis("local",analysisChain,nev);

}
