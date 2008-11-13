#ifndef AliRsnUtils_cxx
#define AliRsnUtils_cxx

class AliRsnUtils
{
  public:
    enum EType
    {
      kLocal = 0,
      kProof,
      kAlien,
      kLastIndex
    };

    enum EDataType
    {
      kTxt = 0,
      kDataSet,
      kXmlCollection,
      kDataLastIndex
    };

    AliRsnUtils(EType type, TString where, TString reset = "", TString rootVer = "");
    AliRsnUtils(EType type, Bool_t shouldRun = kTRUE, Bool_t shouldCopy = kTRUE);
    ~AliRsnUtils() {}

    Bool_t Connect();
    Bool_t CleanPackages(TString pars = "");
    Bool_t LoadPars(TString pars = "STEERBase:ESD:AOD:ANALYSIS:ANALYSISalice:PWG2resonances");
    Bool_t SetInputData(EDataType type = kTxt, TString source = "ESD.txt", TString treeName = "esdTree",
                        TString port = "", Int_t numfiles = 1000000, Int_t filesSkip = 0);
    Bool_t Run(TString macro, Long64_t numEvents = 1, Long64_t numEventsSkip = 0);

    void Print();

    void PrintInfoString(TString s);
    void PrintWarningString(TString s);
    void PrintErrorString(TString s);

    Bool_t IsConnected() const { return fIsConnected; }

    void ShouldCopy(const Bool_t& theValue) { fShouldCopy = theValue; }

    void SetAlienRunFile(const TString& theValue) { fAlienRunFile = theValue;}

    void AddFilesToAlien(const TString& theValue);

    void SetAlienProjectDir(const TString& theValue) { fAlienProjectDir = theValue; }
    TString GetAlienProjectDir() const { return fAlienProjectDir; }

    TAlienJDL* GetJDL() const;

    void ShouldRunAlienJob(const Bool_t& theValue) { fShouldRunAlienJob = theValue; }


  private:

    EType      fCurrentType;              // current analysis type
    EDataType  fCurrentDataType;          // current data type
    Bool_t     fIsConnected;              // flag if we are connected
    Bool_t     fShouldEnd;                // flag if macro should be terminated(in case error)
    TString    fConnectInfo[3];           // connection info (all 3 are use in proof mode)
    TString    fDataInfo[kLastIndex];     // data info
    TString    fTreeName;                 // tree name
    TString    fParNames;                 // pars name
    TChain*    fChain;                    // chain which will be processed

    // AliEn Part
    TString    fAlienProjectDir;          // alien Project directory
    TString    fAlienOutputDir;           // alien output directory
    TString    fAlienRunFile;             // run file name (default is AliRsnAlien.C)
    TString    fFilesToAlien;             // list of files to copy to alien additional to
    // pars. for example (PWG2resonances.C)
    TAlienJDL* fJDL;                      // jdl class which submits job
    Bool_t     fShouldCopy;               // flag if files should be copied to alien
    Bool_t     fShouldRunAlienJob;        // flag if job will be submited

    // connect functions
    Bool_t     ConnectProof();
    Bool_t     ConnectAliEn();

    // clean functions
    Bool_t     CleanParsLocal(TString pars);
    Bool_t     CleanParsProof(TString pars);

    // load functions
    Bool_t     LoadParsLocal(TString pars);
    Bool_t     LoadParsProof(TString pars);

    // process par
    Bool_t     ProcessPAR(TString pars);

    // gets chain from TXT file
    TChain* CreateChainFromTXT(TString chainName = "esdTree", TString aDataDir = "ESDfiles.txt",Int_t aRuns = 200, Int_t offset = 0, TString portNum = "1094");

    // gets chain from xml collection file
    TChain* CreateChainFromCollection(TString chainName = "esdTree", TString collname = "my.xml", Int_t maxFiles = 1000, Int_t skipfiles = 0);


    // copy to AliEn function
    Bool_t CopyFilesToAliEn();

    // generates macro which will be run on alien
    Bool_t CreateRunAliEnFile(TString macro);

    // runs alien job
    Bool_t RunAliEnJob();

    // copies one file to alien
    Bool_t CopyFileToAlien(TString file,Bool_t overwrite = kTRUE);

    ClassDef(AliRsnUtils, 1)
};

#endif

ClassImp(AliRsnUtils)

//________________________________________________________________________
AliRsnUtils::AliRsnUtils(EType type, TString where, TString reset, TString rootVer) :
    fCurrentType(type),
    fCurrentDataType(kTxt),
    fIsConnected(kFALSE),
    fShouldEnd(kFALSE),
    fChain(0),
    fJDL(0),
    fShouldCopy(kTRUE),
    fShouldRunAlienJob(kTRUE),
    fAlienRunFile("AliRsnAlien.C")
{
//
// Constructor
//
  // setting up parameters
  fConnectInfo[0] = where;
  fConnectInfo[1] = reset;
  fConnectInfo[2] = rootVer;

  // connecting
  Connect();
}

//________________________________________________________________________
AliRsnUtils::AliRsnUtils(EType type, Bool_t shouldRun, Bool_t shouldCopy) :
    fCurrentType(type),
    fCurrentDataType(kTxt),
    fIsConnected(kFALSE),
    fShouldEnd(kFALSE),
    fChain(0),
    fJDL(0),
    fShouldRunAlienJob(shouldRun),
    fShouldCopy(shouldCopy),
    fAlienRunFile("AliRsnAlien.C")
{

  // connecting
  Connect();
}

//________________________________________________________________________
void AliRsnUtils::PrintInfoString(const char *msg)
{
  // prints message as Info
  Info("", msg);
}

//________________________________________________________________________
void AliRsnUtils::PrintWarningString(const char *msg)
{
  // prints message as Warning
  Warning("", msg);
}

//________________________________________________________________________
void AliRsnUtils::PrintErrorString(const char *msg)
{
  // prints message as Error
  Error("", msg);
}

//________________________________________________________________________
Bool_t AliRsnUtils::Connect()
{
// connect to local, proof, alien if needed
  switch (fCurrentType)
  {
    case kLocal:
      PrintInfoString("Running Local...");
      fIsConnected = kTRUE;
      break;
    case kProof:
      PrintInfoString(Form("Running Proof at %s...", fConnectInfo[0].Data()));
      fIsConnected = ConnectProof();
      if (!IsConnected())
      {
        PrintErrorString(Form("Not connected to %s !!!", fConnectInfo[0].Data()));
      }
      break;
    case kAlien:
      PrintInfoString(Form("Running AliEn..."));
      fIsConnected = ConnectAliEn();
      if (!IsConnected())
      {
        PrintErrorString("Not connected to AliEn !!!");
      }
      break;
    default:
      PrintErrorString("Wrong type");
      fIsConnected = kFALSE;
      return kFALSE;
  }

  if (!IsConnected())
  {
    PrintErrorString("Connection failed. Aborting...");
    return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRsnUtils::CleanPackages(TString pars)
{
// cleans packages
  if (fShouldEnd) return kFALSE;
  switch (fCurrentType)
  {
    case kLocal:
      if (!IsConnected()) return kFALSE;
      return CleanParsLocal(pars);
      break;
    case kProof:
      if (!IsConnected()) return kFALSE;
      return CleanParsProof(pars);
      break;
    case kAlien:
      return kTRUE;
      break;
    default:
      PrintErrorString("Wrong type");
      return kFALSE;
  }
}

//________________________________________________________________________
Bool_t AliRsnUtils::LoadPars(TString pars)
{
// loads par files
  if (fShouldEnd) return kFALSE;
  switch (fCurrentType)
  {
    case kLocal:
      if (!IsConnected()) return kFALSE;
      return LoadParsLocal(pars);
      break;
    case kProof:
      if (!IsConnected()) return kFALSE;
      return LoadParsProof(pars);
      break;
    case kAlien:
      if (!pars.IsNull())
      {
        fParNames = pars;
        pars.ReplaceAll(":",".par:");
        pars += ".par";
        AddFilesToAlien(pars);
      }
      return kTRUE;
      break;
    default:
      PrintErrorString("Wrong type");
      return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRsnUtils::SetInputData
(EDataType type, TString source, TString treeName, TString port, Int_t numfiles, Int_t filesSkip)
{
// sets input which will be used in Run() function
//

  fDataInfo[type] = source;

  switch (type)
  {
    case kTxt:
      fCurrentDataType = type;
      fChain = CreateChainFromTXT(treeName, fDataInfo[type], numfiles, filesSkip, port);
      if (!fChain) return kFALSE;
      return kTRUE;
      break;
    case kDataSet:
      fCurrentDataType = type;
      break;
    case kXmlCollection:
      fCurrentDataType = type;
      if (fCurrentType == AliRsnUtils::kLocal) ConnectAliEn();
      fChain = CreateChainFromCollection(treeName, fDataInfo[type], numfiles, filesSkip);
      fTreeName = treeName;
      break;
    default:
      PrintErrorString("Wrong type");
      return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRsnUtils::Run(TString macro, Long64_t numEvents, Long64_t numEventsSkip)
{
// runs job
  if (fCurrentType == AliRsnUtils::kAlien)
  {
    Bool_t runFileOK = CreateRunAliEnFile(macro);
    Bool_t copyOK = CopyFilesToAliEn();
    Bool_t runJobOK = RunAliEnJob();
    return runFileOK&&copyOK&&runJobOK;
  }

  gROOT->LoadMacro(macro.Data());
  macro.ReplaceAll(".C","");
  AliAnalysisManager *mgr = (AliAnalysisManager *) gROOT->ProcessLine(Form("%s();", macro.Data()));


  if (mgr->InitAnalysis())
  {
    mgr->PrintStatus();

    switch (fCurrentDataType)
    {
      case kTxt:
        if (!fChain) PrintErrorString("fChain is null");
        TString mode = "local";
        if (fCurrentType == AliRsnUtils::kProof) mode="proof";
        return mgr->StartAnalysis(mode, fChain,numEvents,numEventsSkip);
        break;
      case kDataSet:
        return mgr->StartAnalysis("proof", fDataInfo[kDataSet].Data(),numEvents,numEventsSkip);
        break;
      case kXmlCollection:
        if (!fChain) PrintErrorString("fChain is null");
        return mgr->StartAnalysis("local", fChain,numEvents,numEventsSkip);
        break;
      default:
        PrintErrorString("Wrong type");
        return kFALSE;

    }
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRsnUtils::ConnectProof()
{
// connects to proof
//
// if one uses AliRsnUtils::kProof constructor
// AliRsnUtils ( AliRsnUtils::kProof,TString where = "mvala@lxb6046.cern.ch
// " , TString reset = "RESET", TString rootVer = "HEADXYZ" );
//
// where : means where to connect
// reset : "RESET" will do gProof->Reset(where.Data());
// rootVer : "HEADXYZ" will use different version of root, in this case "HEADXYZ"
//
  if (fConnectInfo[0].IsNull())
  {
    PrintErrorString("Proof is empty!!!");
    return kFALSE;
  }

  if (!fConnectInfo[1].IsNull())
    TProof::Reset(fConnectInfo[0].Data());

  if (!fConnectInfo[2].IsNull())
    TProof::Mgr(fConnectInfo[0].Data())->SetROOTVersion(fConnectInfo[2].Data());

  return TProof::Open(fConnectInfo[0].Data());
}

//________________________________________________________________________
Bool_t AliRsnUtils::ConnectAliEn()
{
// connects to alien
  TGrid::Connect("alien://");
  if (!gGrid)
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRsnUtils::CleanParsLocal(TString pars)
{
// just delete directories which are in pars variable
//
// pars = "STEERBase:ESD"
// will delete STEERBase and ESD directories
//
  TObjArray* array = pars.Tokenize(":");
  TObjString *str;
  TString strr;
  for (Int_t i=0;i< array->GetEntriesFast();i++)
  {
    TObjString *str = (TObjString *) array->At(i);
    strr = str->GetString();
    gSystem->Exec(Form("rm -Rf %s/",strr.Data()));
  }

}

//________________________________________________________________________
Bool_t AliRsnUtils::CleanParsProof(TString pars)
{
// clear par files in proof (divided by ":"). example "STEERBase:ESD"
// will clear STEERBase.par and ESD.par
// if pars is null clears all packages

  TObjArray* array = (TObjArray*) pars.Tokenize(":");

  if (pars.IsNull())
  {
    gProof->ClearPackages();
    return kTRUE;
  }
  TObjString *str;
  TString strr;
  for (Int_t i=0;i< array->GetEntriesFast();i++)
  {
    str = (TObjString *) array->At(i);
    strr = str->GetString();
    gProof->ClearPackage(strr.Data());
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRsnUtils::LoadParsLocal(TString pars)
{
// loads pars localy
  TObjArray* array = pars.Tokenize(":");
  TObjString *str;
  TString strr;
  for (Int_t i=0;i< array->GetEntriesFast();i++)
  {
    str = (TObjString *) array->At(i);
    strr = str->GetString();
    if (!ProcessPAR(strr.Data()))
    {
      return kFALSE;
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRsnUtils::LoadParsProof(TString pars)
{
// load pars on proof
  if (!IsConnected()) return kFALSE;

  TObjArray* array = pars.Tokenize(":");
  TObjString *str;
  TString strr;
  for (Int_t i=0;i< array->GetEntriesFast();i++)
  {
    if (fShouldEnd) return kFALSE;
    str = (TObjString *) array->At(i);
    strr = str->GetString();
    gSystem->Exec(Form("rm -Rf %s",strr.Data()));
    if (gProof->UploadPackage(Form("%s.par", strr.Data())))
    {
      PrintErrorString(Form("Error uploading %s package!!!",strr));
      fShouldEnd = kTRUE;
      return kFALSE;
    }
    if (gProof->EnablePackage(strr))
    {
      PrintErrorString(Form("Error loading %s package!!!",strr));
      fShouldEnd = kTRUE;
      return kFALSE;
    }
    fShouldEnd = kFALSE;
  }
  gProof->ShowEnabledPackages();
  return kTRUE;
}

//________________________________________________________________________
void AliRsnUtils::Print()
{
// prints number of events
  if (fChain)
    PrintInfoString(Form("Number of events %d",fChain->GetEntries()));
}

//________________________________________________________________________
Bool_t AliRsnUtils::CopyFilesToAliEn()
{
// copy needed files to alien
  if (!IsConnected()) return kFALSE;

  if (!fShouldCopy) return kTRUE;

  if (fAlienProjectDir.IsNull())
  {
    fAlienProjectDir = Form("%sRSNTASK/01",gGrid->GetHomeDirectory());
    PrintWarningString(Form("Alien Project directory is not defined.Setting \"%s\" as project directory ...",fAlienProjectDir.Data()));
  }
  gGrid->Rmdir(Form("%s",fAlienProjectDir.Data()));
  gGrid->Mkdir(Form("%s",fAlienProjectDir.Data()),"-p");


  TObjArray* array = fFilesToAlien.Tokenize(":");
  TObjString *str;
  TString strr;
  for (Int_t i=0;i< array->GetEntriesFast();i++)
  {
    str = (TObjString *) array->At(i);
    strr = str->GetString();
    if (!CopyFileToAlien(Form("%s",strr.Data()) ,kFALSE)) return kFALSE;
  }
  if (CopyFileToAlien(fAlienRunFile.Data(),kFALSE)) return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRsnUtils::CreateRunAliEnFile(TString macro)
{
// generated AliRsnAlien.C macro which will be executed in alien
  if (!IsConnected()) return kFALSE;
  TString file = fAlienRunFile.Data();
  ofstream outFile(file.Data());
  if (!outFile)
  {
    PrintErrorString(Form("Cannot open file %s",file.Data()));
    return kFALSE;
  }
  file.ReplaceAll(".C","");
  TString tmp;
  tmp = Form("Int_t %s() {",file.Data());
  outFile << tmp.Data() << endl;

  tmp = "\tgROOT->LoadMacro (\"PWG2resonancesUtils.C\");";
  outFile << tmp.Data() << endl;

  tmp = "\tAliRsnUtils *utils = new AliRsnUtils (AliRsnUtils::kLocal);";
  outFile << tmp.Data() << endl;

  tmp = Form("\tif (!utils->LoadPars(\"%s\")) return -1;",fParNames.Data());
//   tmp = Form("\tif (!utils->LoadPars(\"PWG2resonances\")) return -1;",fParNames.Data());
  outFile << tmp.Data() << endl;

  tmp = Form("\tif (!(utils->SetInputData(AliRsnUtils::kXmlCollection,\"wn.xml\",\"%s\"))) return -1;",fTreeName.Data());
  outFile << tmp.Data() << endl;

  tmp = Form("\tutils->Run(\"%s\",10000000);",macro.Data());
  outFile << tmp.Data() << endl;

  tmp = Form("return 0;",file.Data());
  outFile << tmp.Data() << endl;

  tmp = Form("}",file.Data());
  outFile << tmp.Data() << endl;

  outFile.close();

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRsnUtils::RunAliEnJob()
{
// runs alien job

  fJDL->SetArguments(Form("%s",fAlienRunFile.Data()));

  TObjArray* array = fFilesToAlien.Tokenize(":");
  TObjString *str;
  TString strr;
  for (Int_t i=0;i< array->GetEntriesFast();i++)
  {
    str = (TObjString *) array->At(i);
    strr = str->GetString().Data();
    fJDL->AddToInputSandbox(Form("LF:%s/%s",fAlienProjectDir.Data(),strr.Data())) ;
  }
  fJDL->AddToInputSandbox(Form("LF:%s/%s",fAlienProjectDir.Data(), fAlienRunFile.Data())) ;

  if (fShouldRunAlienJob == kFALSE)
    PrintInfoString(Form("\n%s",fJDL->Generate().Data()));
  else
  {
    TGridJob* job = gGrid->Submit(fJDL->Generate());
    if (job == 0)
    {
      Error("SubmitTest", "submitting failed");
      return kFALSE;
    }
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRsnUtils::CopyFileToAlien(TString file, Bool_t overwrite)
{
// copy one file to alien
  if (fAlienProjectDir.IsNull())
  {
    PrintErrorString("AliEn project directory does not exist...");
    return kFALSE;
  }

  if (overwrite)
    gGrid->Rm(Form("alien://%s/%s",fAlienProjectDir.Data(),file.Data()));

  Info("",Form("Copy from %s to %s",file.Data(),Form("alien:://%s/%s",fAlienProjectDir.Data(),file.Data())));


  return TFile::Cp(file.Data(),Form("alien:://%s/%s",fAlienProjectDir.Data(),file.Data()));
}

//________________________________________________________________________
Bool_t AliRsnUtils::ProcessPAR(TString pars)
{
// process one par

  TString pararchivenameS(pars);
  PrintInfoString(Form("==== Processing %s package LOCAL =====",pararchivenameS.Data()));

  if (gSystem->AccessPathName(Form("%s.par",pararchivenameS.Data())))
  {
    PrintErrorString(Form("File %s.par was not found !!!",pararchivenameS.Data()));
    return kFALSE;
  }

  gROOT->ProcessLine(Form(".! tar xzf %s.par", pararchivenameS.Data()));
  TString ocwd = gSystem->WorkingDirectory();
//     gSystem->ChangeDirectory(Form("%s/%s",dirname,pararchivenameS.Data()));
  gSystem->ChangeDirectory(Form("%s",pararchivenameS.Data()));
  // check for BUILD.sh and execute
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh"))
  {
    PrintInfoString(Form("==== Building %s package =====",pararchivenameS.Data()));
    if (gSystem->Exec("PROOF-INF/BUILD.sh"))
    {
      EPrintErrorString("Cannot Build the PAR Archive! - Abort!");
      return kFALSE;
    }
  }
  // check for SETUP.C and execute
  if (!gSystem->AccessPathName("PROOF-INF/SETUP.C"))
  {
    PrintInfoString(Form("==== Running SETUP.C of %s package =====",pararchivenameS.Data()));
    gROOT->Macro("PROOF-INF/SETUP.C");
    if (gROOT->Macro("PROOF-INF/SETUP.C"))
    {
      PrintErrorString("Cannot SETUP the PAR Archive! - Abort!");
      return kFALSE;
    }
  }

  gSystem->ChangeDirectory(ocwd);

  PrintInfoString(Form("==== All OK for %s package LOCAL =====",pararchivenameS.Data()));
  return kTRUE;

}

//________________________________________________________________________
TChain* AliRsnUtils::CreateChainFromTXT
(TString chainName, TString aDataDir, Int_t aRuns, Int_t offset, TString portNum)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...

  PrintInfoString(Form("Loading the chain %s", chainName.Data()));
  if (aDataDir.IsNull())
  {
    PrintErrorString(Form("aDataDir not found."));
    return kFALSE;
  }

  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
  {
    PrintErrorString(Form("%s not found.", aDataDir.Data()));
    return kFALSE;
  }

  TChain *chain = new TChain(chainName.Data());

  // Open the input stream
  ifstream fileIn;
  fileIn.open(aDataDir.Data());

  Int_t count = 0;

  // Read the input list of files and add them to the chain
  TString esdfile;
  while (fileIn.good())
  {
    fileIn >> esdfile;
    if (esdfile.IsNull()) continue;
    if (offset > 0)
    {
      --offset;
      continue;
    }
    if (count++ == aRuns) break;

    // add esd file
    TString esdFileWithPort(esdfile);
    esdFileWithPort.ReplaceAll("root://","");
    esdFileWithPort.ReplaceAll("//",Form(":%s//",portNum.Data()));
    // PrintInfoString (Form ("Adding root://%s",esdFileWithPort.Data()));
    chain->Add(Form("root://%s",esdFileWithPort.Data()));
  }
  fileIn.close();

  PrintInfoString(Form("Loading the chain %s DONE.", chainName.Data()));
  return chain;
}

//________________________________________________________________________
TChain* AliRsnUtils::CreateChainFromCollection
(TString chainName, TString collname, Int_t maxFiles, Int_t skipfiles)
{
// gets chain from xml collection file

  TChain *chain = new TChain(chainName.Data());
  TAlienCollection *myCollection = TAlienCollection::Open(collname.Data());
  if (!myCollection)
  {
    PrintErrorString(Form("Cannot create an AliEn collection from %s", collectionFile));
    return 0x0;
  }

  // initialize a counter to check the number of read files
  Int_t nfiles = 0;
  TString filename;
  myCollection->Reset();
  while (myCollection->Next())
  {
    if (skipfiles > 0) { --skipfiles; continue; }
    if (maxFiles > 0 && nfiles >= maxFiles) break;
    // char fileName[255];
    // sprintf(fileName, "%s", myCollection->GetTURL(""));
    filename = myCollection->GetTURL("");
//     PrintInfoString(Form("Adding file '%s'", filename.Data()));
    chain->Add(filename.Data());
    nfiles++;
  }

  return chain;
}

//________________________________________________________________________
void AliRsnUtils::AddFilesToAlien(const TString &fileName)
{
// adds file to alien's list which will be copied to alien
  if (fFilesToAlien.IsNull()) fFilesToAlien += fileName;
  else
  {
    fFilesToAlien += ":";
    fFilesToAlien += fileName;
  }
}

//________________________________________________________________________
TAlienJDL* AliRsnUtils::GetJDL() const
{
// gets JDL
  if (!fJDL) fJDL = (TAlienJDL*) gGrid->GetJDLGenerator();
  return fJDL;
}
