/* $Id:  $ */

#ifndef __aodtoolsC__
#define __aodtoolsC__

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TChain.h>
#include <TEnv.h>
#include <TError.h>
#include <TFile.h>
#include <TProof.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTimeStamp.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TAlienCollection.h>

class AODProdParameters : public TObject
{
public:
  AODProdParameters() :
    TObject(),
    runmode("esd"),
    runtype("local"),
    branches("ESD:AliESDRun.,AliESDHeader.,AliESDZDC.,AliESDFMD.,AliESDVZERO."
             ",SPDVertex.,TPCVertex.,PrimaryVertex.,AliMultiplicity.,Tracks"
             ",EMCALCells.,PHOSCells."),
    libs("STEERBase,ESD,AOD,ANALYSIS,ANALYSISalice,CORRFW,EMCALUtils,PWG4CaloCalib,PWG4PartCorrBase,PWG4PartCorrDep"),
    pars(),
    esdfiles("esdfiles.txt"),
    nfiles(3),
    noffset(0),
    xmlfile("collection.xml"),
    dlevel(0),
    doAOD(1),
    doPS(1),
    doCS(1),
    doMC(0)
  {}
  TString runmode;
  TString runtype;
  TString branches;
  TString libs;
  TString pars;
  TString esdfiles;
  Int_t   nfiles;
  Int_t   noffset;
  TString xmlfile;
  Int_t   dlevel;
  Bool_t  doAOD;
  Bool_t  doPS;
  Bool_t  doCS;
  Bool_t  doMC;
  ClassDef(AODProdParameters,1) // AOD production parameters
};
#endif

void LoadLibraries(const AODProdParameters &params);
void SetupPar(const char *pararchivename, Bool_t mflag=kFALSE);
TChain *CreateChain(const AODProdParameters &params);
TChain *CreateESDChain(const char* aDataDir = "ESDfiles.txt", 
                       Int_t aRuns          = 20, 
                       Int_t offset         = 0, 
                       Bool_t addFileName   = kFALSE, 
                       Bool_t addFriend     = kFALSE, 
                       const char* check    = 0);


//_________________________________________________________________________________________________
TChain *CreateChain(const AODProdParameters &params)
{
  if (params.runtype == "grid") {
    ::Error("LoadLibraries", "GRID mode not supported!");

    if (!TFile::Open(params.xmlfile)) {
      ::Error("CreateChain","Collection file %s not found",params.xmlfile.Data());
      return 0;
    }

    TGrid *grid = TGrid::Connect("alien://");
    if (!grid)
      return 0;

    TGridCollection *collection = (TGridCollection*)TAlienCollection::Open(params.xmlfile);
    if (!collection) {
      ::Error("CreateChain","Could not find %s", params.xmlfile.Data());
      return 0;
    }

    TChain *chain = new TChain("esdTree");
    TGridResult *result = collection->GetGridResult("",0 ,0);
    for (Int_t index = 0; index < result->GetEntries(); ++index) {
      TString alienURL(result->GetKey(index, "turl")); 
      Info("CreateChain","Added %s", alienURL.Data());
      chain->Add(alienURL) ; 
    }
    return chain;
  } 

  if (params.runmode == "aod") {
    ::Error("CreateChain", "AOD files not yet supported in local mode");
    return 0;
  }

  return CreateESDChain(params.esdfiles,params.nfiles,params.noffset);
}

//_________________________________________________________________________________________________
TChain* CreateESDChain(const char *aDataDir,
                       Int_t aRuns,
                       Int_t offset,
                       Bool_t addFileName,
                       Bool_t addFriend,
                       const char *check)
{
  // Creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...
  //
  // if addFileName is true the list only needs to contain the directories that contain the AliESDs.root files
  // if addFriend is true a file AliESDfriends.root is expected in the same directory and added to the chain as friend
  // if check is != 0 the files that work are written back into the textfile with the name check

  if (!aDataDir)
    return 0;

  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime)) {
    ::Error("CreateESDChain", "%s not found.", aDataDir);
    return 0;
  }

  TChain *chain = new TChain("esdTree");
  TChain *chainFriend = 0;
  
  if (addFriend)
    chainFriend = new TChain("esdFriendTree");

  if (flags & 2) {
    TString execDir(gSystem->pwd());
    TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
    TList* dirList            = baseDir->GetListOfFiles();
    Int_t nDirs               = dirList->GetEntries();
    gSystem->cd(execDir);

    Int_t count = 0;
    for (Int_t iDir=0; iDir<nDirs; ++iDir) {
      TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
      if (!presentDir || 
          !presentDir->IsDirectory() || 
          strcmp(presentDir->GetName(), ".") == 0 || 
          strcmp(presentDir->GetName(), "..") == 0)
        continue;
      if (offset > 0) {
        --offset;
        continue;
      }
      if (count++ == aRuns)
        break;
      TString presentDirName(aDataDir);
      presentDirName += "/";
      presentDirName += presentDir->GetName();
      chain->Add(presentDirName + "/AliESDs.root/esdTree");
    }
  } else {
    // Open the input stream
    ifstream in;
    in.open(aDataDir);

    ofstream outfile;
    if (check)
      outfile.open(check);

    Int_t count = 0;
    // Read the input list of files and add them to the chain
    TString line;
    while (in.good()) {
      in >> line;
      if (line.Length() == 0)
        continue;
      if (offset > 0) {
        offset--;
        continue;
      }
      TString esdFile(line);
      if (esdFile.BeginsWith("#"))
        continue;
      if (count++ == aRuns)
        break;
      if (esdFile.EndsWith(".zip"))
          esdFile.Append("#AliESDs.root");
      if (addFileName)
        esdFile += "/AliESDs.root";
        
      TString esdFileFriend(esdFile);
      esdFileFriend.ReplaceAll("AliESDs.root", "AliESDfriends.root");
        
      if (check) {
        TFile *file = TFile::Open(esdFile);
        if (!file)
          continue;
        file->Close();

        if (chainFriend) {
          TFile *file2 = TFile::Open(esdFileFriend);
          if (!file2)
            continue;
          file2->Close();
        }
        
        outfile << line.Data() << endl;
        Info("CreateESDChain", "Checked file %s", line.Data());
      }        
      // add esd file
      chain->Add(esdFile);
      if (chainFriend)
        chainFriend->Add(esdFileFriend);
    }
    in.close();
    if (check)
      outfile.close();
  }
  if (chainFriend)
    chain->AddFriend(chainFriend);
  return chain;
}

//_________________________________________________________________________________________________
void LoadLibraries(const AODProdParameters &params) 
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libGeom");
  gSystem->Load("libGui");
  gSystem->Load("libXMLParser");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libMinuit");
  gSystem->Load("libMinuit2");
  gSystem->Load("libProof");
  gSystem->Load("libPhysics");
  gSystem->Load("libMatrix.so");
  if (params.runtype=="grid") {
    gSystem->Load("libNetx.so") ; 
    gSystem->Load("libRAliEn.so"); 
  }

  if (gSystem->Getenv("ALICE_ROOT")) {
    gSystem->AddIncludePath("-I$ALICE_ROOT");
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->AddIncludePath("-I$ALICE_ROOT/PWG2/FLOW/AliFlowCommon");
  }

  if (gSystem->Getenv("TMPDIR")) 
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));

  if (1) {
    delete gRandom;
    gRandom = new TRandom3(0);
    gSystem->Setenv("CONFIG_SEED","0");
  }

  if (params.runtype == "proof") {
    ::Error("LoadLibraries", "PROOF mode not supported!");
    return;
  } 

  TObjArray *libs=params.libs.Tokenize(",");
  if (libs) {
    Int_t n = libs->GetEntries();
    for (Int_t i=0;i<n;++i) {
      TString lib(libs->At(i)->GetName());
      if (lib.Length()<=0)
        continue;
      Bool_t buildpar=0;
      if (params.pars.Contains(lib))
        buildpar=1;
      if (buildpar) {
        SetupPar(lib);
      } else {
        gSystem->Load(Form("lib%s.so",lib.Data()));
      }
      Info("LoadLibraries","Sucessfully loaded %s",lib.Data());
    }
  }
  delete libs;
}

//_________________________________________________________________________________________________
void SetupPar(const char *pararchivename, Bool_t mflag)
{
  // Load par files, create analysis libraries.
  // If mflag then create par file.

  TString cdir(Form("%s", gSystem->WorkingDirectory())) ; 
  TString parpar(Form("%s.par", pararchivename)) ; 

  if (mflag && gSystem->AccessPathName(parpar.Data())) {
    gSystem->ChangeDirectory(gSystem->Getenv("ALICE_ROOT")) ;
    TString processline(Form(".! make %s", parpar.Data())) ; 
    gROOT->ProcessLine(processline.Data()) ;
    gSystem->ChangeDirectory(cdir) ; 
    processline = Form(".! mv $ALICE_ROOT/%s .", parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()) ;
  } 

  if (gSystem->AccessPathName(pararchivename)) {  
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }
  
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);
  
  if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
    ::Info("SetupPar", "Building %s\n", pararchivename);
    
    if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
      ::Error("SetupPar","Can not build %s! - Aborting", pararchivename);
      return;
    }
  }

  // check for SETUP.C and execute
  if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
    ::Info("SetupPar", "Loading %s\n", pararchivename);
    gROOT->Macro("PROOF-INF/SETUP.C");
  }
  
  gSystem->ChangeDirectory(ocwd.Data());
}
#endif
