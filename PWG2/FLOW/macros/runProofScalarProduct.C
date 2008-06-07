void runProofScalarProduct() {
  TStopwatch timer;
  timer.Start();

  printf("*** Connect to PROOF ***\n");
  TProof::Open("snelling@lxb6046.cern.ch");

  gProof->UploadPackage("STEERBase.par");
  gProof->EnablePackage("STEERBase");
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");
  gProof->UploadPackage("ANALYSISalice.par");
  gProof->EnablePackage("ANALYSISalice");
  gProof->UploadPackage("PWG2AOD.par");
  gProof->EnablePackage("PWG2AOD");
  gProof->UploadPackage("PWG2flow.par");
  gProof->EnablePackage("PWG2flow");
  gSystem->SetIncludePath("-I$ROOTSYS/include -I./PWG2flow/FLOW -I./ESD -I./AOD -I./ANALYSIS -I./PWG2AOD/AOD");
  gProof->AddIncludePath("./PWG2AOD/AOD");
  gProof->AddIncludePath("./PWG2flow/FLOW");


  //  TChain *chain = 0x0;

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliESDInputHandler* esdH = new AliESDInputHandler;
  esdH->SetInactiveBranches("FMD CaloCluster");
  mgr->SetInputEventHandler(esdH);  
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTaskScalarProduct *task1 = new AliAnalysisTaskScalarProduct("TaskScalarProduct");

  mgr->AddTask(task1);



  // Create chain of input files
  gROOT->LoadMacro("CreateESDChain.C");
  chain = CreateESDChain("ESD82XX_30K.txt",200);


  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  //  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist1", TList::Class(),AliAnalysisManager::kOutputContainer,"OutputFromCumulantAnlysisESD.root");
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist1", TList::Class(),AliAnalysisManager::kOutputContainer,"outputFromScalarProductAnalysisESD.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("proof",chain);

  timer.Stop();
  timer.Print();
}

