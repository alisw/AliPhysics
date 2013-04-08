// $Id: AddTaskJetEmbeddingFromPYTHIA.C  $

THashTable* GenerateFileTable(const char* list);

AliJetEmbeddingFromPYTHIATask* AddTaskJetEmbeddingFromPYTHIA(
  const char     *tracksName    = "Tracks",
  const char     *clusName      = "",
  const char     *cellsName     = "EMCALCells",
  const char     *MCPartName    = "",
  const char     *simpath       = "/alice/sim/2012/LHC12a15e",
  Int_t           nPtHard       = 11,
  Double_t       *ptHardScaling = 0,
  const char     *aodTreeName   = "aodTree",
  const char     *aodTracksName = "tracks",
  const char     *aodClusName   = "",
  const char     *aodCellsName  = "emcalCells",
  const char     *aodMCPartName = "mcparticles",
  const char     *runperiod     = "lhc12a15e",
  Bool_t          includeNoITS  = kTRUE,
  Double_t        minCent       = -1,
  Double_t        maxCent       = -1,
  UInt_t          mask          = AliVEvent::kAny,
  Double_t        minJetPt      = 0,
  const Bool_t    copyArray     = kTRUE,  
  const Bool_t    makeQA        = kFALSE,
  const char     *fileTable     = "",
  const char     *taskName      = "JetEmbeddingFromPYTHIATask"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJetEmbeddingFromPYTHIA", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetEmbeddingFromPYTHIA", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliJetEmbeddingFromPYTHIATask *jetEmb = new AliJetEmbeddingFromPYTHIATask(taskName,makeQA);
  jetEmb->SetTracksName(tracksName);
  jetEmb->SetClusName(clusName);
  jetEmb->SetCellsName(cellsName);
  jetEmb->SetMCParticlesName(MCPartName);
  jetEmb->SetAODTreeName(aodTreeName);
  jetEmb->SetAODTracksName(aodTracksName);
  jetEmb->SetAODClusName(aodClusName);
  jetEmb->SetAODCellsName(aodCellsName);
  jetEmb->SetAODMCParticlesName(aodMCPartName);
  jetEmb->SetCentralityRange(minCent, maxCent);
  jetEmb->SetTriggerMask(mask);
  jetEmb->SetCopyArray(copyArray);
  jetEmb->SetJetMinPt(minJetPt);
  jetEmb->SetNClusters(1);
  jetEmb->SetMarkMC(0);

  if (strcmp(fileTable, "") != 0)
    jetEmb->SetFileTable(GenerateFileTable(fileTable));

  jetEmb->SetIncludeNoITS(includeNoITS);
  TString runPeriod(runperiod);
  runPeriod.ToLower();
  if (runPeriod == "lhc12a15e") {
    jetEmb->SetAODfilterBits(256,512);
  }
  else if (runPeriod == "lhc12a15a") {
    jetEmb->SetAODfilterBits(256,16);
  }
  else {
    if (!runPeriod.IsNull())
      ::Warning("Run period %s not known. It will use IsHybridGlobalConstrainedGlobal.", runPeriod.Data());
  }

  jetEmb->SetPYTHIAPath(simpath);

  if (nPtHard > 0) {
    if (ptHardScaling==0) {
      ptHardScaling = new Double_t[nPtHard];
      for (Int_t i = 0; i < nPtHard; i++)
	ptHardScaling[i] = 1;
    }
    jetEmb->SetPtHardBinScaling(nPtHard, ptHardScaling);
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

THashTable* GenerateFileTable(const char* list)
{
  THashTable *table = new THashTable();

  TString myList = list;
  if (myList.Contains("alien:///")) {
    TFile::Cp(myList,"file:./list.txt");
    myList = "./list.txt";
  }

  // Open the input stream
  ifstream in;
  in.open(myList.Data());

  // Read the input list of files and add them to the chain
  TString line;
  while (in.good()) {
    in >> line;

    if (line.Length() == 0)
      continue;

    TObjString *aodFile = new TObjString(line);
    table->Add(aodFile);
  }

  return table;
}
