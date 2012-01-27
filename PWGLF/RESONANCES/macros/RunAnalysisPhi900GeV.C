//
// This is an example of steering macro for running RSN analysis task
// locally with a collection of files written in a text file
//
// Allowed inputs:
// - ESD + MC (MC is included automatically)
// - AOD
//
// All settings are specified as arguments in the main macro,
// which is at the end of this script.
//
void RunAnalysisPhi900GeV
(
  Int_t       nReadFiles   = 2,
  Int_t       nSkipFiles   = 0,
  const char *dataType     = "900GeV_pass4_sim",
  const char *inputSource  = "list.txt",
  const char *outName1     = "Phi900GeV_all.root",
  const char *outName2     = "Phi900GeV_true.root",
  const char *outName3     = "Phi900GeV_info.root"
)
{
  // convert the last argument into a BOOL variable
  Bool_t isMC = kTRUE;
  if (!strcmp(dataType, "900GeV_pass4_data")) isMC = kFALSE;
  if (!strcmp(dataType, "7TeV_pass1_data")) isMC = kFALSE;

  // message on aliroot version
  cout << "*** ALIROOT PATH = " << gSystem->Getenv("ALICE_ROOT") << " ***" << endl;
  cout << "*** MC " << (isMC ? "" : "NOT") << " INCLUDED ***" << endl;

  // check extension of input to distinguish between XML and TXT
  TString sInput(inputSource);
  sInput.ToLower();
  Bool_t isTXT = (!strcmp(sInput(sInput.Length() - 3, 3).Data(), "txt"));
  cout << "Input = " << (isTXT ? "TXT" : "XML") << endl;

  // load compiled libraries (for aliroot session)
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG2resonances.so");

  // if input is XML, connect to AliEn
  if (!isTXT) TGrid::Connect("alien://");

  // create analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MyTaskManager");

  // create handlers for input
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdHandler);
  // if required, create also MC handler
  if (isMC)
  {
    AliMCEventHandler *mcHandler  = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
  }
  
  // add event selection for data
  gROOT->LoadMacro("AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);

  // add task macro
  AliRsnAnalysisPhi900GeV *task = new AliRsnAnalysisPhi900GeV("taskphi900gev");
  task->SelectCollisionCandidates();
  task->SetTPCparams(isMC);
  task->SetMaxChi2(4.0);
  task->SetMaxDCAr(0.5);
  task->SetMaxDCAz(3.0);
  task->SetMinNTPC(80);
  task->SetUseMC(kFALSE);
  if (!strcmp(dataType, "900GeV_pass4_data"))
  {
    task->SetTOFESD(kFALSE);
    task->SetTOFSigma(130.0);
    task->SetTOFSettings(AliTOFT0makerANA::kPass4);
  }
  if (!strcmp(dataType, "7TeV_pass1_data"))
  {
    task->SetTOFESD(kFALSE);
    task->SetTOFSigma(130.0);
    task->SetTOFSettings(AliTOFT0makerANA::kPass4);
  }
  else if (!strcmp(dataType, "900GeV_pass4_sim"))
  {
    task->SetTOFESD(kTRUE);
    task->SetTOFSigma(130.0);
    task->SetTOFSettings(AliTOFT0makerANA::kNone);
  }
  mgr->AddTask(task);
  
  // create containers for input/output
  AliAnalysisDataContainer *out1 = mgr->CreateContainer("tracks", TTree::Class(), AliAnalysisManager::kOutputContainer, outName1);
  AliAnalysisDataContainer *out2 = mgr->CreateContainer("rsn"   , TTree::Class(), AliAnalysisManager::kOutputContainer, outName2);
  AliAnalysisDataContainer *out3 = mgr->CreateContainer("info"  , TList::Class(), AliAnalysisManager::kOutputContainer, outName3);
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, out1);
  mgr->ConnectOutput(task, 2, out2);
  mgr->ConnectOutput(task, 3, out3);

  // create TChain of input events
  TChain *analysisChain = 0x0;
  if (isTXT) analysisChain = CreateChainFromText(inputSource, "esdTree", nReadFiles, nSkipFiles);
  else       analysisChain = CreateChainFromXML (inputSource, "esdTree", nReadFiles, nSkipFiles);

  // start analysis
  if (!analysisChain)
  {
    Error("runLocal", "Analysis chain not properly initialized");
    return;
  }
  mgr->InitAnalysis();
  mgr->PrintStatus();
  if (isTXT) mgr->StartAnalysis("local", analysisChain);
  else       mgr->StartAnalysis("alien", analysisChain);
}

//_________________________________________________________________________________________________
TChain* CreateChainFromXML
(const char *xmlFileName, const char *treeName, Int_t nread, Int_t nskip)
{
//
// Create a TChain with all required files listed into an XML collection.
// Necessary to run analysis in AliEn jobs.
// ---
// Arguments:
//  - xmlFileName = input list
//  - treeName    = "esdTree" or "aodTree"
//  - nread       = how many files to read (0 = all)
//  - nskip       = how many files to skip from beginning
//

  // if nread argument is 0, it is disabled
  if (nread == 0) nread = 1000000000;

  // initialize output object
  TChain *chain = new TChain(treeName);

  // initialize the AliEn collection
  TAlienCollection *myCollection = TAlienCollection::Open(xmlFileName);
  if (!myCollection)
  {
    Error("CreateChainFromXML", "Cannot create an AliEn collection from %s", xmlFileName);
    return 0x0;
  }

  // loop on collection
  myCollection->Reset();
  while (myCollection->Next())
  {
    // skip until reached required number of offset
    if (nskip > 0) {--nskip; continue;}

    // stop if required number of read files is reached
    // otherwise update the counter
    if (nread <= 0) break;
    nread--;

    // recovery file and add it
    Info("CreateChainFromXML", Form("Adding: %s", myCollection->GetTURL("")));
    chain->Add(myCollection->GetTURL(""));
  }

  return chain;
}

//_________________________________________________________________________________________________
TChain* CreateChainFromText(const char *fileName, const char *treeName, Int_t nread, Int_t nskip)
{
//
// Create a TChain with all required files listed into a text file.
// Necessary to run analysis in local jobs.
// ---
// Arguments:
//  - xmlFileName = input file list
//  - treeName    = "esdTree" or "aodTree"
//  - nread       = how many files to read (0 = all)
//  - nskip       = how many files to skip from beginning
//

  // if third argument is 0, it is interpreted
  // as "read all lines"
  Bool_t readAll = (nread <= 0);
  
  // initialize output object
  TChain* target = new TChain(treeName);
  
  // open text file
  ifstream fileIn(fileName);
  
  // loop on collection
  TString line;
  while (fileIn.good())
  {
    fileIn >> line;
    if (line.IsNull()) continue;
    
    // skip until reached required number of offset
    if (nskip > 0) {--nskip; continue;}
    
    // stop if required number of read files is reached
    // otherwise update the counter
    if (!readAll && nread <= 0) break;
    nread--;
    
    // add file
    Info("CreateChainFromText", "Adding '%s'", line.Data());
    target->Add(line.Data());
  }
  
  return target;
}
