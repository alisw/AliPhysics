//--------------------------------------------------------------------------
// Base macro for submitting single muon analysis.
//
// In case it is not run with full aliroot, it needs to have in the working directory:
//  - STEERBase.par
//  - AOD.par
//  - ESD.par
//  - ANALYSIS.par
//  - ANALYSISalice.par
//  - PWG3muon.par
//
// The macro reads ESDs and outputs file:
// - SingleMuESD.root
//--------------------------------------------------------------------------

void RunSingleMuonAnalysisFromESD(Bool_t local = kFALSE) {

  TStopwatch timer;
  timer.Start();

  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libAOD");
  gSystem->Load("libESD");
  gSystem->Load("libPWG3muon.so");

  TChain* chain = new TChain("esdTree");

  if (!local) {
    printf("*** Connect to AliEn ***\n");
    TGrid::Connect("alien://");
    TAlienCollection* coll = TAlienCollection::Open("wn.xml");
    TGridResult* result = coll->GetGridResult("",0,0);
    for(Int_t i = 0; i < result->GetEntries(); i++) {
      printf("TURL = %s \n",result->GetKey(i,"turl"));
      chain->Add(result->GetKey(i,"turl"));
    }
  } else {
    chain->Add("/path_1_to/AliESDs.root");
    chain->Add("/path_2_to/AliESDs.root");
    chain->Add("/path_3_to/AliESDs.root");
  }

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  

  //____________________________________________//
  // ntuple task
  AliAnalysisTaskSingleMuESD *task = new AliAnalysisTaskSingleMuESD("TaskSingleMuESD");
  // force values for the trigger mask, for version where 
  // GetFiredTriggerClasses() does not work
  //task->SetTriggerType("MUON");
  mgr->AddTask(task);

  // Create containers for input/output

  // input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  // output
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("ctree", TNtuple::Class(),AliAnalysisManager::kOutputContainer,"SingleMuESD.root");
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("chist", TH1F::Class(),AliAnalysisManager::kOutputContainer,"SingleMuESD.root");

  //____________________________________________//
  mgr->ConnectInput(task,0,cinput);
  mgr->ConnectOutput(task,0,coutput1);
  mgr->ConnectOutput(task,1,coutput2);

  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
  }
  
  timer.Stop();
  timer.Print();
}

