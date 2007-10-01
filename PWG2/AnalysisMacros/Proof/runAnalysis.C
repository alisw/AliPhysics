void runAnalysis() {
  TProof::Open("lxb6046.cern.ch");

  // Enable the Analysis Package
  gProof->UploadPackage("ESD.par");
  gProof->EnablePackage("ESD");
  gProof->UploadPackage("AOD.par");
  gProof->EnablePackage("AOD");
  gProof->UploadPackage("ANALYSIS.par");
  gProof->EnablePackage("ANALYSIS");

  gROOT->LoadMacro("CreateESDChain.C");
  TChain* chain = CreateESDChain("ESD100_110_v4.txt", 100);

  gProof->Load("AliAnalysisTaskPt.cxx+");

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  // Add Pt task
  AliAnalysisTask *task1 = new AliAnalysisTaskPt("TaskPt");
  mgr->AddTask(task1);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist2", TH1::Class(),    AliAnalysisManager::kOutputContainer, "Pt.ESD.1.root");

  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);

  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  mgr->StartAnalysis("proof",chain);
}
