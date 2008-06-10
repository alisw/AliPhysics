void runProof(const char * dataset = "/PWG0/COMMON/run30000X_10TeV_0.5T", Long64_t nentries=1234567890, Long64_t firstentry=0)
{
  // Connect to Proof
  TProof::Open("lxb6046");

  // Upload and enable packages: please use the correct version!
  gProof->UploadPackage("AF-v4-12");
  gProof->EnablePackage("AF-v4-12");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AliTPCComparison");

  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  // Enable MC event handler
  AliVEventHandler* handler = new AliMCEventHandler;
  mgr->SetMCtruthEventHandler(handler);

  // Create task
  gProof->Load("AliMCComparisonTrack.cxx++g");
  gProof->Load("AliTPCComparisonTask.cxx++g");
  AliAnalysisTask *task = new AliTPCComparisonTask("AliTPCComparisonTask");

  // Add task
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("coutput", TList::Class(), AliAnalysisManager::kOutputContainer, "Hist.root");

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);


  // Enable debug printouts
  mgr->SetDebugLevel(3);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  mgr->StartAnalysis("proof",dataset,nentries,firstentry);
}
