//
// This is an example steering macro for running RSN analysis task
// locally with a collection of files specified in a text file:
//
// Inputs:
//   - nReadFiles  = number of files to process from the list
//   - nSkipFiles  = how many lines to be skipped when reading the list
//   - addTaskName = name of the macro to add the RSN analysis task
//                   (assumed to have inside it a function named like the file)
//   - inputSource = name of the file containing all the inputs
//                   ---> to run on a local collection, the collection file 
//                        must contain on each line the full path 
//                        of one input file and it must have the ".txt" extension
//                   ---> to run on an AliEn collection, the collection file must be an XML
//                        file collection like those built from the "find -x" method in aliensh.
//   - dataLabel   = a label which is used to know what kind of data are being read
//                   (it is propagated to the 'addTask' macro for eventual setting up of something
//   - outName     = name for the file with RSN package outputs (without ROOT extension)
//
// Notes:
//   - in case the source is an ESD, and if inputs are a MC production
//     the MC input handler is created by default
// 
//
// In principle, the user should never modify this macro. 
//
void runLocal
(
  Int_t       nReadFiles  = 0,
  Int_t       nSkipFiles  = 0,
  const char *addTaskName = "AddAnalysisTaskRsnTest.C",
  const char *inputSource = "/home/pulvir/analysis/resonances/LHC2010-7TeV-phi/alien+plugin/sim.txt",
  const char *dataLabel   = "7TeV_pass2_sim_ESD",
  const char *outName     = "rsn-test"
)
{
  
  // convert the last argument into a BOOL variable
  TString strDataLabel(dataLabel);
  Bool_t isESD = strDataLabel.Contains("ESD");
  Bool_t isAOD = strDataLabel.Contains("AOD");
  Bool_t isSim = strDataLabel.Contains("sim");   
  
  //AliLog::SetGlobalDebugLevel(AliLog::kDebug+1);
  //AliLog::SetClassDebugLevel("AliRsnCutPID", AliLog::kDebug+2);
  //AliLog::SetClassDebugLevel("AliRsnPair", AliLog::kDebug+2);
  //AliLog::SetClassDebugLevel("AliRsnPairFunctions", AliLog::kDebug+2);
  //AliLog::SetClassDebugLevel("AliRsnValue", AliLog::kDebug+4);

  // check extension of input to distinguish between XML and TXT
  TString sInput(inputSource);
  sInput.ToLower();
  Bool_t isTXT = (!strcmp(sInput(sInput.Length() - 3, 3).Data(), "txt"));
  cout << "Input = " << (isTXT ? "TXT" : "XML") << endl;

  // load compiled libraries (for aliroot session)
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG2resonances.so");
  
  // if input is XML, connect to AliEn
  if (!isTXT) TGrid::Connect("alien://");

  // create analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("taskRsnTest");
  mgr->SetCommonFileName(Form("%s.root", outName));
  
  // create input handler
  if (isESD)
  {
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdHandler);
    // if possible, create also MC handler
    if (isSim)
    {
      AliMCEventHandler *mcHandler  = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcHandler);
    }
  }
  else if (isAOD)
  {
    AliAODInputHandler *aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);
  }
  else
  {
    ::Error("Required an ESD or AOD input data set");
    return;
  }
  
  // add event selection for data
  gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isSim);
  
  // split the argument for macros in order to add many tasks
  TString    sList(addTaskName);
  TObjArray *list = sList.Tokenize(":");
  for (Int_t i = 0; i < list->GetEntries(); i++)
  {
    TObjString *os  = (TObjString*)list->At(i);
    TString     str = os->GetString();
    gROOT->ProcessLine(Form(".x %s(\"%s\")", str.Data(), dataLabel));
  }

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
Bool_t LoadPars(const char *parList, const char *path)
{
//
// Load PAR libraries locally
// ---
// Arguments:
//  - parList = list of PARs without extension, separated by ':'
//  - path    = path where PARs are stored
//

  // store position of working directory
  TString ocwd = gSystem->WorkingDirectory();

  // tokenize list
  TString     pars(parList);
  TObjArray  *array = pars.Tokenize(":");

  // loop on list
  TObjString *ostr;
  TString     str;
  Char_t      parName[200], parFile[200];
  for (Int_t i = 0; i < array->GetEntriesFast(); i++)
  {
    ostr = (TObjString*) array->At(i);
    str = ostr->GetString();
    sprintf(parName, "%s", str.Data());
    sprintf(parFile, "%s/%s.par", path, str.Data());

    // check that file exists
    if (!gSystem->AccessPathName(parFile))
    {
      // explode tar-ball and enter it
      gROOT->ProcessLine(Form(".! tar xzf %s", parFile));
      gSystem->ChangeDirectory(Form("%s", parName));
      // checks for BUILD.sh and execute it
      if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh"))
      {
        ::Info("", ">> Building PARs: %s", parName);
        if (gSystem->Exec("PROOF-INF/BUILD.sh"))
        {
          ::Error("LoadPars", Form("BUILD.sh error for '%s'", parFile));
          gSystem->ChangeDirectory(ocwd);
          return kFALSE;
        }
      }
      // check and execute SETUP.C
      if (!gSystem->AccessPathName("PROOF-INF/SETUP.C"))
      {
        ::Info("", ">> Setting up PARs: %s", parName);
        if (gROOT->Macro("PROOF-INF/SETUP.C"))
        {
          Error("LoadPars", Form("SETUP.C error for '%s'", parFile));
          gSystem->ChangeDirectory(ocwd);
          return kFALSE;
        }
      }
    }
    else
    {
      Error("LoadParsLocal", Form("File '%s' not found", parFile));
      return kFALSE;
    }
  }

  gSystem->ChangeDirectory(ocwd);
  return kTRUE;
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
