Bool_t AliHMPIDAnalysisTask()
{
  
  TBenchmark benchmark;
  benchmark.Start("AliHMPIDAnalysisTask");

  AliLog::SetGlobalDebugLevel(0);

  Load() ; //load the required libraries

  TChain * analysisChain ;
   analysisChain = new TChain("esdTree");
    //here put your input data path 
   analysisChain->Add("AliESDs.root");


  Info("AliHMPIDAnalysisTask",Form("N events %d",(Int_t)analysisChain->GetEntries()));

  // create the task
  AliHMPIDAnalysisTask *task = new AliHMPIDAnalysisTask("CosmicAnalysisTask");
  task->SetTrigger(2);

 
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

 
  AliESDInputHandler* esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);

  // Create and connect containers for input/output

  //------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();

  // ----- output data -----
  
  //slot 0 : default output tree (by default handled by AliAnalysisTaskSE)
//  AliAnalysisDataContainer *coutput0 = mgr->CreateContainer("ctree0", TTree::Class(),AliAnalysisManager::kOutputContainer,"Houtput.root");

  //now comes user's output objects :
  
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist0", TH1I::Class(),AliAnalysisManager::kOutputContainer,"Houtput.root");
  // output list of histos
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("clist0", TList::Class(),AliAnalysisManager::kOutputContainer,"Houtput.root");
 
  cinput0->SetData(analysisChain);

  mgr->AddTask(task);
  mgr->ConnectInput(task,0,cinput0);
//  mgr->ConnectOutput(task,0,coutput0);
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
 

  //RUN !!!
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",analysisChain);
  }

  benchmark.Stop("AliHMPIDAnalysisTask");
  benchmark.Show("AliHMPIDAnalysisTask");

  return kTRUE ;


}

void Load() {

  //load the required aliroot libraries
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;

  //compile online the task class
  gSystem->SetIncludePath("-I. -I$HOME/HMPID -I$ALICE_ROOT/include -I$ROOTSYS/include");
  gROOT->LoadMacro("./AliHMPIDAnalysisTask.cxx+");
}
