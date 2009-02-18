
Bool_t AliPriorsTask()
{
  
  TBenchmark benchmark;
  benchmark.Start("AliPriorsTask");

  AliLog::SetGlobalDebugLevel(0);

  Load() ; //load the required libraries

  TChain * analysisChain = new TChain("esdTree");
  analysisChain->Add("AliESDs.root");

  Int_t nIter = 5;

  // create the task
  AliPriorsTask *task = new AliPriorsTask("Task for Priors");
  Double_t priors[5]={0.2,0.2,0.2,0.2,0.2};
  task->SetPriors(priors);
  task->SetNiterations(nIter);

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestIteration");

  AliMCEventHandler*  mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);


   AliESDInputHandler* esdHandler = new AliESDInputHandler();
   mgr->SetInputEventHandler(esdHandler);


  // Create and connect containers for input/output
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  //slot 0 : default output tree (by default handled by AliAnalysisTaskSE)
  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("ctree0", TTree::Class(),AliAnalysisManager::kOutputContainer,"output.root");
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TH1I::Class(),AliAnalysisManager::kOutputContainer,"output.root");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cgraph0", TH2D::Class(),AliAnalysisManager::kOutputContainer,"output.root");


  cinput0->SetData(analysisChain);

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput0);
  mgr->ConnectOutput(task,0,coutput0);
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);



  AliEventPoolLoop* loop = new AliEventPoolLoop(nIter);
  loop->SetChain(analysisChain);
  mgr->SetEventPool(loop);

  //mgr->SetDebugLevel(1); 
  printf("READY TO RUN\n");
  //RUN !!!
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("mixing",analysisChain);
    //mgr->StartAnalysis("local",analysisChain);
    
  }

  benchmark.Stop("AliPriorsTask");
  benchmark.Show("AliPriorsTask");

  return kTRUE ;
}

void Load() {

  //load the required aliroot libraries
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;

  //compile online the task class
  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ROOTSYS/include");
  gROOT->LoadMacro("./AliPriorsTask.cxx++g");
}
