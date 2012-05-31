void runProofT0analysis(const char * dataset = "/COMMON/COMMON/LHC09a4_10TeV_200k#esdTree",Long64_t nentries=20000, Long64_t firstentry=0)
{
// Connect to Proof
  TProof::Open("proof://alla@alicecaf.cern.ch"); 
  //TProof::Open("lxb6046");

  // Upload and enable packages: please use the correct version!
  gProof->UploadPackage("AF-v4-16");
  gProof->EnablePackage("AF-v4-16");
  gProof->ShowDataSets();      
 
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AliT0MultiplicityTask");

  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  // Enable MC event handler
  AliVEventHandler* handler = new AliMCEventHandler;
  mgr->SetMCtruthEventHandler(handler);

  // Create task
  //  gProof->Load("AliMCComparisonTrack.cxx++g");
  gProof->Load("AliT0MultiplicityTask.cxx++g");
  AliAnalysisTask *task = new AliT0MultiplicityTask("AliT0MultiplicityTask");

  // Add task
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = 
    mgr->CreateContainer("cchain", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput = 
    mgr->CreateContainer("coutput", TList::Class(), 
    AliAnalysisManager::kOutputContainer, "MultHist.root");

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

