void run(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE, const char* option = "")
{
  if (aProof)
  {
    TProof::Open("lxb6046");

    // Enable the needed package
    gProof->UploadPackage("STEERBase");
    gProof->EnablePackage("STEERBase");
    gProof->UploadPackage("ESD");
    gProof->EnablePackage("ESD");
    gProof->UploadPackage("AOD");
    gProof->EnablePackage("AOD");
    gProof->UploadPackage("ANALYSIS");
    gProof->EnablePackage("ANALYSIS");
    gProof->UploadPackage("PWG0base");
    gProof->EnablePackage("PWG0base");

    //gProof->AddInput(new TNamed("PROOF_Packetizer", "TPacketizer"));
  }
  else
  {
    gSystem->Load("libVMC");
    gSystem->Load("libTree");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libPWG0base");
  }

  // Create chain of input files
  gROOT->LoadMacro("../CreateESDChain.C");
  chain = CreateESDChain(data, nRuns, offset);

  // Create the analysis manager
  mgr = new AliAnalysisManager("testAnalysis");

  // selection of esd tracks
  gROOT->ProcessLine(".L CreateCuts.C");
  AliESDtrackCuts* esdTrackCuts = CreateTrackCuts();
  if (!esdTrackCuts)
  {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  }

  TString taskName("AliMultiplicityTask.cxx+");
  if (aDebug)
    taskName += "+g";

  // Create, add task
  if (aProof) {
    gProof->Load(taskName);
  } else
    gROOT->Macro(taskName);

  task = new AliMultiplicityTask;
  task->SetTrackCuts(esdTrackCuts);
  mgr->AddTask(task);

  // Enable MC event handler
  AliMCEventHandler* handler = new AliMCEventHandler;
  handler->SetReadTR(kFALSE);
  //mgr->SetMCtruthEventHandler(handler);

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);

  // Attach input
  cInput  = mgr->CreateContainer("cInput", TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer);
  //cOutput->SetDataOwned(kTRUE);
  mgr->ConnectOutput(task, 0, cOutput);

  // Enable debug printouts
  if (aDebug)
    mgr->SetDebugLevel(2);

  // Run analysis
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis((aProof) ? "proof" : "local", chain);
}

void runAll()
{
  run("part1.txt", 1000000, 0, kFALSE, kTRUE);
  gSystem->Exec("mv multiplicityMC.root multiplicityMC_1.root");
  run("part2.txt", 1000000, 0, kFALSE, kTRUE);
  gSystem->Exec("mv multiplicityMC.root multiplicityMC_2.root");
}
