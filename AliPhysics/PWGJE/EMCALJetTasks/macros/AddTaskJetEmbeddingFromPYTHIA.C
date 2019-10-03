// $Id$

THashTable* GenerateFileTable(const char* list);

AliJetEmbeddingFromPYTHIATask* AddTaskJetEmbeddingFromPYTHIA(
  const char     *tracksName    = "Tracks",
  const char     *clusName      = "",
  const char     *cellsName     = "EMCALCells",
  const char     *MCPartName    = "",
  const char     *simpath       = "alien:///alice/sim/2012/LHC12a15e_fix/%d/%d/AOD149/%04d/AliAOD.root",
  Int_t           nPtHard       = 0,
  Double_t       *ptHardScaling = 0,
  const char     *aodTreeName   = "aodTree",
  const char     *aodTracksName = "tracks",
  const char     *aodClusName   = "",
  const char     *aodCellsName  = "emcalCells",
  const char     *aodMCPartName = "mcparticles",
  const char     *runperiod     = "lhc12a15e",
  Bool_t          includeNoITS  = kFALSE,
  Double_t        minCent       = -1,
  Double_t        maxCent       = -1,
  UInt_t          mask          = 0,
  Double_t        minJetPt      = 0,
  const Bool_t    copyArray     = kFALSE,  
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
      ::Warning("AddTaskJetEmbeddingFromPYTHIA","Run period %s not known. It will use IsHybridGlobalConstrainedGlobal.", runPeriod.Data());
  }

  jetEmb->SetPYTHIAPath(simpath);

  if (nPtHard == 0 || ptHardScaling == 0) { // the pt hard bin scaling was not provided, use the default for LHC12a15e_fix AOD149
    ::Warning("AddTaskJetEmbeddingFromPYTHIA","The pT hard bin scaling has not been provided, will use the default for LHC12a15e_fix AOD149!", runPeriod.Data());
    nPtHard = 11;
    ptHardScaling = new Double_t[nPtHard];
    ptHardScaling [0]  = 0;
    ptHardScaling [1]  = 5.135193e-05;
    ptHardScaling [2]  = 5.859497e-06;
    ptHardScaling [3]  = 4.444755e-07;
    ptHardScaling [4]  = 4.293118e-08;
    ptHardScaling [5]  = 5.154750e-09;
    ptHardScaling [6]  = 6.958612e-10;
    ptHardScaling [7]  = 1.149828e-10;
    ptHardScaling [8]  = 2.520137e-11;
    ptHardScaling [9]  = 6.222240e-12;
    ptHardScaling [10] = 2.255832e-12;
  }
  jetEmb->SetPtHardBinScaling(nPtHard, ptHardScaling);

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
