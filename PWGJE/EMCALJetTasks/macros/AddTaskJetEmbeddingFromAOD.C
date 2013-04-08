// $Id: AddTaskJetEmbeddingFromAOD.C  $

TObjArray* GenerateFileList(const char* list, Int_t nFiles);

AliJetEmbeddingFromAODTask* AddTaskJetEmbeddingFromAOD(
  const char     *tracksName    = "Tracks",
  const char     *clusName      = "",
  const char     *cellsName     = "EMCALCells",
  const char     *MCPartName    = "",
  const char     *fileList      = "files.txt",
  const char     *aodTreeName   = "aodTree",
  const char     *aodTracksName = "tracks",
  const char     *aodClusName   = "",
  const char     *aodCellsName  = "emcalCells",
  const char     *aodMCPartName = "",
  const char     *runperiod     = "lhc11h",
  Bool_t          includeNoITS  = kTRUE,
  Double_t        minCent       = 0,
  Double_t        maxCent       = 10,
  UInt_t          mask          = AliVEvent::kAny,
  const Bool_t    copyArray     = kTRUE,  
  const Bool_t    makeQA        = kFALSE,
  Int_t           nFiles        = 1234567890,
  const char     *taskName      = "JetEmbeddingFromAODTask"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJetEmbeddingFromAOD", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetEmbeddingFromAOD", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliJetEmbeddingFromAODTask *jetEmb = new AliJetEmbeddingFromAODTask(taskName,makeQA);
  jetEmb->SetTracksName(tracksName);
  jetEmb->SetClusName(clusName);
  jetEmb->SetCellsName(cellsName);
  jetEmb->SetMCParticlesName(MCPartName);
  jetEmb->SetFileList(GenerateFileList(fileList, nFiles));
  jetEmb->SetAODTreeName(aodTreeName);
  jetEmb->SetAODTracksName(aodTracksName);
  jetEmb->SetAODClusName(aodClusName);
  jetEmb->SetAODCellsName(aodCellsName);
  jetEmb->SetAODMCParticlesName(aodMCPartName);
  jetEmb->SetCentralityRange(minCent, maxCent);
  jetEmb->SetTriggerMask(mask);
  jetEmb->SetCopyArray(copyArray);
  jetEmb->SetNClusters(1);
  jetEmb->SetMarkMC(0);

  jetEmb->SetIncludeNoITS(includeNoITS);
  TString runPeriod(runperiod);
  runPeriod.ToLower();
  if (runPeriod == "lhc11h" || runPeriod == "lhc12a15e") {
    jetEmb->SetAODfilterBits(256,512); // hybrid tracks for LHC11h and LHC12a15e
  }
  else if (runPeriod == "lhc11a" || runPeriod == "lhc12a15a") {
    jetEmb->SetAODfilterBits(256,16); // hybrid tracks for LHC11a and LHC12a15a
  }
  else {
    if (!runPeriod.IsNull())
      ::Warning("Run period %s not known. It will use IsHybridGlobalConstrainedGlobal.", runPeriod.Data());
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetEmb);
    
  // Create containers for input/output
  mgr->ConnectInput(jetEmb, 0, mgr->GetCommonInputContainer());

  if (makeQA) {
    TString contName = taskName;
    contName += "_histos";
    AliAnalysisDataContainer *outc = mgr->CreateContainer(contName,
							  TList::Class(),
							  AliAnalysisManager::kOutputContainer,
							  "AnalysisResults.root");
    mgr->ConnectOutput(jetEmb, 1, outc);
  }

  return jetEmb;
}

TObjArray* GenerateFileList(const char* list, Int_t nFiles)
{
  TObjArray *array = new TObjArray(9999);

  TString myList = list;
  if (myList.Contains("alien:///")) {
    TFile::Cp(myList,"file:./list.txt");
    myList = "./list.txt";
  }

  // Open the input stream
  ifstream in;
  in.open(myList.Data());

  Int_t count = 0;

  // Read the input list of files and add them to the chain
  TString line;
  while (in.good()) {
    if (nFiles != 1234567890) {
      if (count >= nFiles)
	break;
    }

    in >> line;

    if (line.Length() == 0)
      continue;

    TObjString *aodFile = new TObjString(line);
    array->Add(aodFile);
    
    count++;
  }

  return array;
}
