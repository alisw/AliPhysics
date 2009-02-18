void run(Char_t* data, Int_t nRuns=20, Int_t offset=0, Bool_t aDebug = kFALSE, Bool_t aProof = kFALSE, Bool_t mc = kFALSE, const char* option = "")
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

  TString taskName("AliCutTask.cxx+");
  if (aDebug)
    taskName += "+g";

  // Create, add task
  if (aProof) {
    gProof->Load(taskName);
  } else
    gROOT->Macro(taskName);

  task = new AliCutTask;

  AliPWG0Helper::AnalysisMode analysisMode = AliPWG0Helper::kTPC;
  task->SetAnalysisMode(analysisMode);

  if (analysisMode != AliPWG0Helper::kSPD)
  {
    // selection of esd tracks
    gROOT->ProcessLine(".L ../CreateStandardCuts.C");
    AliESDtrackCuts* esdTrackCuts = CreateTrackCuts(analysisMode);
    if (!esdTrackCuts)
    {
      printf("ERROR: esdTrackCuts could not be created\n");
      return;
    }

    task->SetTrackCuts(esdTrackCuts);
  }

  mgr->AddTask(task);

  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  
  if (mc) {
    task->EnableSecondaryStudy();
    // Enable MC event handler
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // Attach input
  cInput  = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cInput);

  // Attach output
  cOutput = mgr->CreateContainer("cOutput", TList::Class(), AliAnalysisManager::kOutputContainer);
  mgr->ConnectOutput(task, 0, cOutput);

  // Enable debug printouts
  if (aDebug)
    mgr->SetDebugLevel(2);

  // Run analysis
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis((aProof) ? "proof" : "local", chain);
}
