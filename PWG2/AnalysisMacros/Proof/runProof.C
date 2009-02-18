void runProof() {
  TStopwatch timer;
  timer.Start();
  
  runProofESD("AliAnalysisTaskPt.cxx++");

  timer.Stop();
  timer.Print();
}

//==========================================//
void runProofESD(const char *file) {
  //the next line should point to the local $ALICE_ROOT
  //that contains the latest ANALYSIS developments
  printf("****** Connect to PROOF *******\n");
  TProof::Open("proof://lxb6046.cern.ch"); 
  gProof->SetParallel(1);

  // Enable the Analysis Package
  gProof->UploadPackage("STEERBase.par");
  gProof->EnablePackage("STEERBase");
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");
  
  // You should get this macro and the txt file from:
  // http://aliceinfo.cern.ch/Offline/Analysis/CAF/
  gROOT->LoadMacro("CreateESDChain.C");
  TChain* chain = 0x0;
  chain = CreateESDChain("ESD82XX_30K.txt",10);
  chain->SetBranchStatus("*Calo*",0);

  gProof->Load(file);
 
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);  
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTaskPt *task1 = new AliAnalysisTaskPt("TaskPt");
  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TH1::Class(),AliAnalysisManager::kOutputContainer,"Pt.ESD.root");
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->SetDebugLevel(2);
  mgr->StartAnalysis("proof",chain);
}
