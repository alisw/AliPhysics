void RunAliTPCCalibKrTask(const char* list="KrClusters_250508.txt",Bool_t bProof = kFALSE)
{

  if(bProof) {
    TProof *proof = TProof::Open("jacek@gsiaf.gsi.de");
    gProof->GetManager()->SetROOTVersion("5.18/00a");

    // Proof Enable Libraries
    gROOT->LoadMacro("ProofEnableAliRoot.C");
    ProofEnableAliRoot("/d/alice11/jacek/alice/x86_64/AliRoot/HEAD");
  }

  //
  // Create chain of input files
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx++");
  AliXRDPROOFtoolkit tool;

  // -- Make chain of files
  TChain * chain = tool.MakeChain(list,"Kr","",2,0);
  chain->SetBranchStatus("Cl.fCluster",kFALSE);
  //
  // Create the analysis manager
  //
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  // Calibration component 
  AliTPCCalibKr *calibObj = new AliTPCCalibKr;
  //calibObj->SetASide(kFALSE);

  // Add task 
  AliTPCCalibKrTask *task = new AliTPCCalibKrTask;
  task->SetInputChain(chain);
  task->SetTPCCalibKr(calibObj);
  mgr->AddTask(task); 

  // Attach input
  cInput  = mgr->CreateContainer("cInput", TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer,"outHistFile.root");
  mgr->ConnectOutput(task, 0, cOutput);

  // Run analysis
  mgr->InitAnalysis();
  mgr->PrintStatus();
  if(bProof) mgr->StartAnalysis("proof", chain);
  else mgr->StartAnalysis("local", chain);
}
