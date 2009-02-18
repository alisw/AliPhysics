void runProof(int dataType=0, const char *dataSource="ESD82XX_30K.txt") {
  // Run on PROOF with data from:
  // dataType =
  // 0 - a local file list
  // 1 - a PROOF dataset
  // 
  // for dataTpye = 0
  // dataSource is the list file name
  // for dataType = 1
  // dataSource is the PROOF dataset name

  TStopwatch timer;
  timer.Start();

  printf("*** Connect to PROOF ***\n");
  // ****
  // You have to change this to Your own username !!!!
  TProof::Open("akisiel@lxb6046.cern.ch");
  //
  // ****

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
  gProof->UploadPackage("PWG2femtoscopy.par");
  gProof->EnablePackage("PWG2femtoscopy");
  gProof->UploadPackage("PWG2femtoscopyUser.par");
  gProof->EnablePackage("PWG2femtoscopyUser");
  gSystem->SetIncludePath("-I$ROOTSYS/include -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser -I./ESD -I./AOD -I./ANALYSIS -I./PWG2AOD/AOD");
  gProof->AddIncludePath("./PWG2AOD/AOD");
  gProof->AddIncludePath("./PWG2femtoscopy/FEMTOSCOPY/AliFemto");
  gProof->AddIncludePath("./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser");

  gProof->Load("ConfigFemtoAnalysis.C++g");
  gProof->Load("AliAnalysisTaskFemto.cxx++g");

  TChain *chain = 0x0;

  if (dataType == 0) {
    gROOT->LoadMacro("CreateESDChain.C");
    chain = CreateESDChain("ESD82XX_30K.txt",200);
  }

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliESDInputHandler* esdH = new AliESDInputHandler;
  esdH->SetInactiveBranches("FMD CaloCluster");
  mgr->SetInputEventHandler(esdH);  
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTaskFemto *task1 = new AliAnalysisTaskFemto("TaskFemto");

  mgr->AddTask(task1);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clist1", TList::Class(),AliAnalysisManager::kOutputContainer,"Femto.ESD.root");
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if (dataType == 0)
    mgr->StartAnalysis("proof",chain);
  else if (dataType == 1)
    mgr->StartAnalysis("proof",dataSource, -1, 0);

  timer.Stop();
  timer.Print();
}

