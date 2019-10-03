void runProofFwdDetsQA(const char * dataset = "/COMMON/COMMON/LHC09a4_run8101X",Long64_t nentries=100000, Long64_t firstentry=0)
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");

  // Connect to Proof
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  TProof::Open("cheshkov:PWG0@alicecaf");

  // Upload and enable packages: please use the correct version!
  gProof->UploadPackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-16-Release/AF-v4-16");
  gProof->EnablePackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-16-Release/AF-v4-16");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AliAnaFwdDetsQA");

  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  // Enable MC event handler
  AliVEventHandler* handler = new AliMCEventHandler;
  mgr->SetMCtruthEventHandler(handler);

  // Create task

  gProof->Load(Form("%s/PWGPP/AliAnaFwdDetsQA.cxx++g",
		    gSystem->Getenv("ALICE_ROOT")));
  AliAnalysisTask *task = new AliAnaFwdDetsQA("AliAnaFwdDetsQA");

  // Add task
  mgr->AddTask(task);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = 
    mgr->CreateContainer("coutput", TList::Class(), 
    AliAnalysisManager::kOutputContainer, "FwdDetsQA.root");

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);


  // Enable debug printouts
  mgr->SetDebugLevel(3);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  TFileCollection *proofColl = gProof->GetDataSet(dataset);
  TChain *chain = new TChain("esdTree");
  chain->AddFileInfoList((TCollection*)(proofColl->GetList()));
  mgr->StartAnalysis("proof", chain, nentries, firstentry);

    //  mgr->StartAnalysis("proof",dataset,nentries,firstentry);
}

