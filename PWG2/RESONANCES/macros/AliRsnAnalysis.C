//=============================================================================
//
// *** AliRsnAnalysis.C ***
//
// This macro contains all stuff which is required to run
// a resonance analysis in whatever environment.
// It contains a pool of macroes which do all the steps to set it up
// and run it with the appropriate configurations.
// All the possibilities are made available through enumerations,
// with the same style of the Config.C file used for simulation.
//
//=============================================================================

#include <TProof.h>
#include <TGrid.h>

enum Rsn_DataSource
{
  kTextFile,        // text file with paths of all used files (local analysis)
  kXMLCollection,   // XML collection of files (local/AliEn analysis)
  kDataset          // CAF dataset
};

enum Rsn_Environment
{
  kLocal,   // local analysis with a locally defined TXT list of files
  kAlien,   // AliEn analysis with an XML collection
  kProof    // CAF analysis with a DataSet
};

static Bool_t            cleanPars = kFALSE;
static Bool_t            cleanPWG2resonances = kFALSE;
static const char*       pathPar = "/home/pulvir/ALICE/ALIROOT/par-files/afs_proof";
static const char*       listPar = "STEERBase:ESD:AOD:ANALYSIS:ANALYSISalice:PWG2resonances";
static const char*       proofConnection = "localhost";

// Uncomment these lines for an example run with a local list of local files
//static const char*       inputSource = "local.txt";
//static Rsn_DataSource    listType = kTextFile;
//static Rsn_Environment   environment = kLocal;

// Uncomment these lines for an example run with PROOF
static const char*       inputSource = "/COMMON/COMMON/LHC08c11_10TeV_0.5T";
static Rsn_DataSource    listType    = kDataset;
static Rsn_Environment   environment = kProof;

// Uncomment these lines for an example run with AliEn XML collection
//static const char*       inputSource = "wn.xml";
//static Rsn_DataSource    listType    = kXMLCollection;
//static Rsn_Environment   environment = kAlien;

static Int_t             nReadFiles = 5;
static Int_t             nSkippedFiles = 0;
static TString           treeName = "esdTree";

//_________________________________________________________________________________________________
Bool_t IsProof()
{
//
// Check if the environment is PROOF
//
  return (environment == kProof);
}

//_________________________________________________________________________________________________
Bool_t IsAlien()
{
//
// Check if the environment is Alien
//
  return (environment == kAlien);
}

//_________________________________________________________________________________________________
Bool_t CleanPars(const char *pars)
{
//
// Cleans existing PAR library directories.
// The string argument must contain their names
// separated by colons (':').
// Returns kTRUE if everything went well, otherwise returns kFALSE.
//

  TString list(pars);
  TObjArray* array = list.Tokenize(":");
  TObjString *ostr;
  TString str;

  for (Int_t i = 0; i < array->GetEntriesFast(); i++)
  {
    ostr = (TObjString *) array->At(i);
    str = ostr->GetString();
    Info("", ">> Cleaning PARs: %s", str.Data());
    if (IsProof()) {
      if (!gProof) {
        Error("CleanPars", "gProof object not initialized");
        return kFALSE;
      }
      gProof->ClearPackage(Form("%s.par", str.Data()));
    }
    else {
      gSystem->Exec(Form("rm -Rf %s/", str.Data()));
    }
  }

  return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t LoadPars(TString pars)
{
//
// Loads required PAR libraries.
// The string argument must contain their names
// separated by colons (':').
// Returns kTRUE if everything went well, otherwise returns kFALSE.
//

  TObjArray *array = pars.Tokenize(":");
  TObjString *ostr;
  TString str;

  for (Int_t i = 0; i < array->GetEntriesFast(); i++) {
    ostr = (TObjString *) array->At(i);
    str = ostr->GetString();
    Info("", ">> Creating PARs: %s", str.Data());
    if (IsProof()) {
       if (!gProof) {
        Error("CleanPars", "gProof object not initialized");
        return kFALSE;
      }
      if (!str.CompareTo("PWG2resonances")) {
        gProof->UploadPackage(Form("%s.par", str.Data()));
      }
      else {
          gProof->UploadPackage(Form("%s/%s.par", pathPar, str.Data()));
      }
      gProof->EnablePackage(str.Data());
    }
    else {
      // check that file exists...
      if (gSystem->AccessPathName(Form("%s.par", str.Data()))) {
        Error("", " !!! File not found");
        return kFALSE;
      }
      // ...explode tar-ball...
      gROOT->ProcessLine(Form(".! tar xzf %s.par", str.Data()));
      // ...go to unpacked directory of PAR library...
      TString ocwd = gSystem->WorkingDirectory();
      gSystem->ChangeDirectory(Form("%s", str.Data()));
      // ...checks for BUILD.sh and execute it...
      if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
        if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
          Error("", " !!! Build error");
          return kFALSE;
        }
      }
      // ...check and execute SETUP.C...
      if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
        if (gROOT->Macro("PROOF-INF/SETUP.C")) {
          Error("", " !!! Set-up error");
          return kFALSE;
        }
      }
      // ...return to main dir and finish
      gSystem->ChangeDirectory(ocwd);
    }
  }

  return kTRUE;
}

//_________________________________________________________________________________________________
TChain* CreateChainFromTXT()
{
//
// Create a TChain of files to be analyzed from a text file.
// This text file must contain in each line the FULL PATH of an analyzed file.
// Last argument specifies the TTree name.
//

  TChain *chain = new TChain(treeName.Data());

  ifstream fileIn(inputSource);
  Int_t count = 0, skip = nSkippedFiles;

  TString line;
  while (fileIn.good()) {
    fileIn >> line;
    if (line.IsNull()) continue;
    if (skip > 0) {
      --skip;
      continue;
    }
    if (nReadFiles > 0 && count++ >= nReadFiles) break;
    chain->Add(Form("%s", line.Data()));
  }
  fileIn.close();

  return chain;
}

//_________________________________________________________________________________________________
TChain* CreateChainFromXML()
{
//
// Create a TChain of files to be analyzed from a XML file.
// Last argument specifies the TTree name.
//

  if (!gGrid) {
    Error("CreateChainFromXML", "gGrid object not instantiated");
    return 0;
  }

  TChain *chain = new TChain(treeName.Data());
  TAlienCollection *myCollection = TAlienCollection::Open(inputSource);
  if (!myCollection)
  {
    Error("CreateChainFromXML", "Cannot create an AliEn collection from %s", collectionFile);
    return 0x0;
  }

  // initialize a counter to check the number of read files
  Int_t nfiles = 0, skip = nSkippedFiles;
  TString filename;
  myCollection->Reset();
  while (myCollection->Next())
  {
    if (skip > 0) { --skip; continue; }
    if (nReadFiles > 0 && nfiles >= nReadFiles) break;
    // char fileName[255];
    // sprintf(fileName, "%s", myCollection->GetTURL(""));
    filename = myCollection->GetTURL("");
    chain->Add(filename.Data());
    nfiles++;
  }

  return chain;
}

//_________________________________________________________________________________________________
TChain* CreateChainFromDataset()
{
//
// Create a TChain of files to be analyzed from a CAF dataset.
// Last argument specifies the TTree name.
//

  if (!gProof) {
    Error("CreateChainFromDataset", "gProof object not initialized");
    return 0;
  }

  TFileCollection *fc = gProof->GetDataSet(inputSource)->GetStagedSubset();
  TIter iter(fc->GetList());

  TChain* target = new TChain(treeName.Data());

  TFileInfo* fileInfo = 0;
  Int_t nfiles = 0, skip = nSkippedFiles;
  while ((fileInfo = dynamic_cast<TFileInfo*> (iter()))) {
    if (fileInfo->GetFirstUrl()) {
      if (skip > 0) { --skip; continue; }
      if (nReadFiles > 0 && nfiles >= nReadFiles) break;
      target->Add(fileInfo->GetFirstUrl()->GetUrl());
      nfiles++;
    }
  }

  return target;
}

//_________________________________________________________________________________________________
Int_t AliRsnAnalysis(const char *addMacro = "AddRsnAnalysisTask.C")
{
//
// Main macro (named after the macro filename).
// It initializes the job, creates the AnalysisManager and calls the config macro which
// fills this AnalysisManager with all AnalysisTask objects for resonance analysis.
// ---
// The argument is the name of the macro used for configuration
// ---
// Return values:
// - 0 = successful execution
// - 1 = problem while cleanind PAR libraries
// - 2 = problem while loading PAR libraries
// - 3 = invalid TTree name (allowed: esdTree, aodTree)
// - 4 = data source not found
// - 5 = wrong definition of source type (allowed any member of enum Rsn_DataSource)
// - 6 = TChain initialization returned a NULL object
// - 7 = error while initializing AliAnalysisManager object
//

  // connect to PROOF if required
  if (IsProof()) TProof::Open(proofConnection);
  else if (IsAlien()) TGrid::Connect("alien://");

  //
  // *** SETUP PAR LIBRARIES **********************************************************************
  //  - libraries are cleaned if requested
  //

  if (cleanPars) {
    if (!CleanPars(listPar)) {
      Error("AliRsnAnalysis", "Error while cleaning PAR libraries");
      return 1;
    }
  }
  if (cleanPWG2resonances) CleanPars("PWG2resonances");
  if (!LoadPars(listPar)) {
    Error("AliRsnAnalysis", "Error while initializing PAR libraries");
    return 2;
  }

  // set up log level
  AliLog::SetGlobalLogLevel(AliLog::kError);

  // *** CREATE TChain OF INPUT EVENTS ************************************************************

  // preliminary check #1: is the tree name OK?
  if (treeName != "esdTree" && treeName != "aodTree") {
    Error("AliRsnAnalysis", "Invalid tree name specified");
    return 3;
  }

  // preliminary check #2: does the input source exist?
  // (valid only for local files and XML collections)
  if (listType == kTextFile || listType == kXMLCollection) {
    Long_t id, size, flags, modtime;
    if (gSystem->GetPathInfo(inputSource, &id, &size, &flags, &modtime)) {
      Error("AliRsnAnalysis", "Input source not found");
      return 4;
    }
  }

  // create TChain
  TChain *analysisChain = 0;
  switch (listType) {
    case kTextFile:
      analysisChain = CreateChainFromTXT();
      break;
    case kXMLCollection:
      analysisChain = CreateChainFromXML();
      break;
    case kDataset:
      analysisChain = CreateChainFromDataset();
      break;
    default:
      Error("AliRsnAnalysis", "Input type not supported");
      return 5;
  }
  if (!analysisChain) {
    Error("AliRsnAnalysis", "Analysis TChain not properly initialized");
    return 6;
  }

  //
  // *** INITIALIZE ANALYSIS MANAGER **************************************************************
  //

  AliAnalysisManager *mgr = new AliAnalysisManager("ResonanceAnalysis");
  if (!mgr) {
    Error("AliRsnAnalysis", "Error occurred while initializing AnalysisManager");
    return 7;
  }

  // set ESD and MC input handlers
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  AliAODInputHandler *aodHandler = new AliAODInputHandler();
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  esdHandler->SetInactiveBranches("FMD CaloCluster");
  if (!treeName.CompareTo("esdTree")) {
    mgr->SetInputEventHandler(esdHandler);
  }
  else if (!treeName.CompareTo("aodTree")) {
    mgr->SetInputEventHandler(aodHandler);
  }
  else {
    Error("AliRsnAnalysis", "Input tree type not supported");
    return 8;
  }
  mgr->SetMCtruthEventHandler(mcHandler); 

  // load configuration macro and uses it to initialize this object
  gROOT->LoadMacro(addMacro);
  AddRsnAnalysisTask();

  // initialize analysis and run it
  TString strEnv = (IsProof() ? "proof" : "local");
  TStopwatch timer;
  timer.Start();
  mgr->InitAnalysis();
  mgr->StartAnalysis(strEnv.Data(), analysisChain);
  timer.Stop();
  timer.Print();

  return 0;
}
