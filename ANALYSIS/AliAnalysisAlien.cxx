/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Author: Mihaela Gheata, 01/09/2008

//==============================================================================
//   AliAnalysisAlien - AliEn utility class. Provides interface for creating
// a personalized JDL, finding and creating a dataset.
//==============================================================================

#include "Riostream.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TGridCollection.h"
#include "TGridJDL.h"
#include "TGridJobStatusList.h"
#include "TGridJobStatus.h"
#include "TFileMerger.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisAlien.h"

ClassImp(AliAnalysisAlien)

//______________________________________________________________________________
AliAnalysisAlien::AliAnalysisAlien()
                 :AliAnalysisGrid(),
                  fGridJDL(NULL),
                  fPrice(0),
                  fTTL(0),
                  fSplitMaxInputFileNumber(0),
                  fMaxInitFailed(0),
                  fMasterResubmitThreshold(0),
                  fNtestFiles(0),
                  fNrunsPerMaster(0),
                  fMaxMergeFiles(0),
                  fNsubmitted(0),
                  fProductionMode(0),
                  fRunNumbers(),
                  fExecutable(),
                  fExecutableCommand(),
                  fArguments(),
                  fAnalysisMacro(),
                  fAnalysisSource(),
                  fAdditionalLibs(),
                  fSplitMode(),
                  fAPIVersion(),
                  fROOTVersion(),
                  fAliROOTVersion(),
                  fUser(),
                  fGridWorkingDir(),
                  fGridDataDir(),
                  fDataPattern(),
                  fGridOutputDir(),
                  fOutputArchive(),
                  fOutputFiles(),
                  fInputFormat(),
                  fDatasetName(),
                  fJDLName(),
		            fMergeExcludes(),
                  fIncludePath(),
                  fCloseSE(),
                  fFriendChainName(),
                  fJobTag(),
                  fInputFiles(0),
                  fPackages(0)
{
// Dummy ctor.
   SetDefaults();
}

//______________________________________________________________________________
AliAnalysisAlien::AliAnalysisAlien(const char *name)
                 :AliAnalysisGrid(name),
                  fGridJDL(NULL),
                  fPrice(0),
                  fTTL(0),
                  fSplitMaxInputFileNumber(0),
                  fMaxInitFailed(0),
                  fMasterResubmitThreshold(0),
                  fNtestFiles(0),
                  fNrunsPerMaster(0),
                  fMaxMergeFiles(0),
                  fNsubmitted(0),
                  fProductionMode(0),
                  fRunNumbers(),
                  fExecutable(),
                  fExecutableCommand(),
                  fArguments(),
                  fAnalysisMacro(),
                  fAnalysisSource(),
                  fAdditionalLibs(),
                  fSplitMode(),
                  fAPIVersion(),
                  fROOTVersion(),
                  fAliROOTVersion(),
                  fUser(),
                  fGridWorkingDir(),
                  fGridDataDir(),
                  fDataPattern(),
                  fGridOutputDir(),
                  fOutputArchive(),
                  fOutputFiles(),
                  fInputFormat(),
                  fDatasetName(),
                  fJDLName(),
                  fMergeExcludes(),
                  fIncludePath(),
                  fCloseSE(),
                  fFriendChainName(),
                  fJobTag(),
                  fInputFiles(0),
                  fPackages(0)
{
// Default ctor.
   SetDefaults();
}

//______________________________________________________________________________
AliAnalysisAlien::AliAnalysisAlien(const AliAnalysisAlien& other)
                 :AliAnalysisGrid(other),
                  fGridJDL(NULL),
                  fPrice(other.fPrice),
                  fTTL(other.fTTL),
                  fSplitMaxInputFileNumber(other.fSplitMaxInputFileNumber),
                  fMaxInitFailed(other.fMaxInitFailed),
                  fMasterResubmitThreshold(other.fMasterResubmitThreshold),
                  fNtestFiles(other.fNtestFiles),
                  fNrunsPerMaster(other.fNrunsPerMaster),
                  fMaxMergeFiles(other.fMaxMergeFiles),
                  fNsubmitted(other.fNsubmitted),
                  fProductionMode(other.fProductionMode),
                  fRunNumbers(other.fRunNumbers),
                  fExecutable(other.fExecutable),
                  fExecutableCommand(other.fExecutableCommand),
                  fArguments(other.fArguments),
                  fAnalysisMacro(other.fAnalysisMacro),
                  fAnalysisSource(other.fAnalysisSource),
                  fAdditionalLibs(other.fAdditionalLibs),
                  fSplitMode(other.fSplitMode),
                  fAPIVersion(other.fAPIVersion),
                  fROOTVersion(other.fROOTVersion),
                  fAliROOTVersion(other.fAliROOTVersion),
                  fUser(other.fUser),
                  fGridWorkingDir(other.fGridWorkingDir),
                  fGridDataDir(other.fGridDataDir),
                  fDataPattern(other.fDataPattern),
                  fGridOutputDir(other.fGridOutputDir),
                  fOutputArchive(other.fOutputArchive),
                  fOutputFiles(other.fOutputFiles),
                  fInputFormat(other.fInputFormat),
                  fDatasetName(other.fDatasetName),
                  fJDLName(other.fJDLName),
                  fMergeExcludes(other.fMergeExcludes),
                  fIncludePath(other.fIncludePath),
                  fCloseSE(other.fCloseSE),
                  fFriendChainName(other.fFriendChainName),
                  fJobTag(other.fJobTag),
                  fInputFiles(0),
                  fPackages(0)
{
// Copy ctor.
   fGridJDL = (TGridJDL*)gROOT->ProcessLine("new TAlienJDL()");
   fRunRange[0] = other.fRunRange[0];
   fRunRange[1] = other.fRunRange[1];
   if (other.fInputFiles) {
      fInputFiles = new TObjArray();
      TIter next(other.fInputFiles);
      TObject *obj;
      while ((obj=next())) fInputFiles->Add(new TObjString(obj->GetName()));
      fInputFiles->SetOwner();
   }   
   if (other.fPackages) {
      fPackages = new TObjArray();
      TIter next(other.fPackages);
      TObject *obj;
      while ((obj=next())) fPackages->Add(new TObjString(obj->GetName()));
      fPackages->SetOwner();
   }   
}

//______________________________________________________________________________
AliAnalysisAlien::~AliAnalysisAlien()
{
// Destructor.
   if (fGridJDL) delete fGridJDL;
   if (fInputFiles) delete fInputFiles;
   if (fPackages) delete fPackages;
}   

//______________________________________________________________________________
AliAnalysisAlien &AliAnalysisAlien::operator=(const AliAnalysisAlien& other)
{
// Assignment.
   if (this != &other) {
      AliAnalysisGrid::operator=(other);
      fGridJDL = (TGridJDL*)gROOT->ProcessLine("new TAlienJDL()");
      fPrice                   = other.fPrice;
      fTTL                     = other.fTTL;
      fSplitMaxInputFileNumber = other.fSplitMaxInputFileNumber;
      fMaxInitFailed           = other.fMaxInitFailed;
      fMasterResubmitThreshold = other.fMasterResubmitThreshold;
      fNtestFiles              = other.fNtestFiles;
      fNrunsPerMaster          = other.fNrunsPerMaster;
      fMaxMergeFiles           = other.fMaxMergeFiles;
      fNsubmitted              = other.fNsubmitted;
      fProductionMode          = other.fProductionMode;
      fRunNumbers              = other.fRunNumbers;
      fExecutable              = other.fExecutable;
      fExecutableCommand       = other.fExecutableCommand;
      fArguments               = other.fArguments;
      fAnalysisMacro           = other.fAnalysisMacro;
      fAnalysisSource          = other.fAnalysisSource;
      fAdditionalLibs          = other.fAdditionalLibs;
      fSplitMode               = other.fSplitMode;
      fAPIVersion              = other.fAPIVersion;
      fROOTVersion             = other.fROOTVersion;
      fAliROOTVersion          = other.fAliROOTVersion;
      fUser                    = other.fUser;
      fGridWorkingDir          = other.fGridWorkingDir;
      fGridDataDir             = other.fGridDataDir;
      fDataPattern             = other.fDataPattern;
      fGridOutputDir           = other.fGridOutputDir;
      fOutputArchive           = other.fOutputArchive;
      fOutputFiles             = other.fOutputFiles;
      fInputFormat             = other.fInputFormat;
      fDatasetName             = other.fDatasetName;
      fJDLName                 = other.fJDLName;
      fMergeExcludes           = other.fMergeExcludes;
      fIncludePath             = other.fIncludePath;
      fCloseSE                 = other.fCloseSE;
      fFriendChainName         = other.fFriendChainName;
      fJobTag                  = other.fJobTag;
      if (other.fInputFiles) {
         fInputFiles = new TObjArray();
         TIter next(other.fInputFiles);
         TObject *obj;
         while ((obj=next())) fInputFiles->Add(new TObjString(obj->GetName()));
         fInputFiles->SetOwner();
      }   
      if (other.fPackages) {
         fPackages = new TObjArray();
         TIter next(other.fPackages);
         TObject *obj;
         while ((obj=next())) fPackages->Add(new TObjString(obj->GetName()));
         fPackages->SetOwner();
      }   
   }
   return *this;
}

//______________________________________________________________________________
void AliAnalysisAlien::AddIncludePath(const char *path)
{
// Add include path in the remote analysis macro.
   TString p(path);
   if (p.Contains("-I")) fIncludePath += Form("%s ", path);
   else                  fIncludePath += Form("-I%s ", path);
}

//______________________________________________________________________________
void AliAnalysisAlien::AddRunNumber(Int_t run)
{
// Add a run number to the list of runs to be processed.
   if (fRunNumbers.Length()) fRunNumbers += " ";
   fRunNumbers += Form("%d", run);
}   

//______________________________________________________________________________
void AliAnalysisAlien::AddRunNumber(const char* run)
{
// Add a run number to the list of runs to be processed.
   if (fRunNumbers.Length()) fRunNumbers += " ";
   fRunNumbers += run;
}   

//______________________________________________________________________________
void AliAnalysisAlien::AddDataFile(const char *lfn)
{
// Adds a data file to the input to be analysed. The file should be a valid LFN
// or point to an existing file in the alien workdir.
   if (!fInputFiles) fInputFiles = new TObjArray();
   fInputFiles->Add(new TObjString(lfn));
}
   
//______________________________________________________________________________
Bool_t AliAnalysisAlien::Connect()
{
// Try to connect to AliEn. User needs a valid token and /tmp/gclient_env_$UID sourced.
   if (gGrid && gGrid->IsConnected()) return kTRUE;
   if (!gSystem->Getenv("alien_API_USER")) {
      Error("Connect", "Make sure you:\n 1. Have called: alien-token-init <username> today\n 2. Have sourced /tmp/gclient_env_%s",
            gSystem->Getenv("UID"));
      return kFALSE;
   }         
   if (!gGrid) {
      Info("Connect", "Trying to connect to AliEn ...");
      TGrid::Connect("alien://");
   }
   if (!gGrid || !gGrid->IsConnected()) {
      Error("Connect", "Did not managed to connect to AliEn. Make sure you have a valid token.");
      return kFALSE;
   }  
   fUser = gGrid->GetUser();
   Info("Connect", "\n#####   Connected to AliEn as user %s. Setting analysis user to <%s>", fUser.Data(), fUser.Data());
   return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisAlien::CdWork()
{
// Check validity of alien workspace. Create directory if possible.
   if (!Connect()) {
      Error("CdWork", "Alien connection required");
      return;
   } 
   TString homedir = gGrid->GetHomeDirectory();
   TString workdir = homedir + fGridWorkingDir;
   if (DirectoryExists(workdir)) {
      gGrid->Cd(workdir);
      return;
   }   
   // Work directory not existing - create it
   gGrid->Cd(homedir);
   if (gGrid->Mkdir(workdir)) {
      gGrid->Cd(fGridWorkingDir);
      Info("CreateJDL", "\n#####   Created alien working directory %s", fGridWorkingDir.Data());
   } else {
      Warning("CreateJDL", "Working directory %s cannot be created.\n Using %s instead.",
              workdir.Data(), homedir.Data());
      fGridWorkingDir = "";
   }          
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::CheckInputData()
{
// Check validity of input data. If necessary, create xml files.
   if (!fInputFiles && !fRunNumbers.Length() && !fRunRange[0]) {
      if (!fGridDataDir.Length()) {
         Error("CkeckInputData", "AliEn path to base data directory must be set.\n = Use: SetGridDataDir()");
         return kFALSE;
      }
      Info("CheckInputData", "Analysis will make a single xml for base data directory %s",fGridDataDir.Data());
      return kTRUE;
   }
   // Process declared files
   Bool_t is_collection = kFALSE;
   Bool_t is_xml = kFALSE;
   Bool_t use_tags = kFALSE;
   Bool_t checked = kFALSE;
   CdWork();
   TString file;
   TString workdir = gGrid->GetHomeDirectory();
   workdir += fGridWorkingDir;
   if (fInputFiles) {
      TObjString *objstr;
      TIter next(fInputFiles);
      while ((objstr=(TObjString*)next())) {
         file = workdir;
         file += "/";
         file += objstr->GetString();
         // Store full lfn path
         if (FileExists(file)) objstr->SetString(file);
         else {
            file = objstr->GetName();
            if (!FileExists(objstr->GetName())) {
               Error("CheckInputData", "Data file %s not found or not in your working dir: %s",
                     objstr->GetName(), workdir.Data());
               return kFALSE;
            }         
         }
         Bool_t iscoll, isxml, usetags;
         CheckDataType(file, iscoll, isxml, usetags);
         if (!checked) {
            checked = kTRUE;
            is_collection = iscoll;
            is_xml = isxml;
            use_tags = usetags;
            TObject::SetBit(AliAnalysisGrid::kUseTags, use_tags);
         } else {
            if ((iscoll != is_collection) || (isxml != is_xml) || (usetags != use_tags)) {
               Error("CheckInputData", "Some conflict was found in the types of inputs");
               return kFALSE;
            } 
         }
      }
   }
   // Process requested run numbers
   if (!fRunNumbers.Length() && !fRunRange[0]) return kTRUE;
   // Check validity of alien data directory
   if (!fGridDataDir.Length()) {
      Error("CkeckInputData", "AliEn path to base data directory must be set.\n = Use: SetGridDataDir()");
      return kFALSE;
   }
   if (!DirectoryExists(fGridDataDir)) {
      Error("CheckInputData", "Data directory %s not existing.", fGridDataDir.Data());
      return kFALSE;
   }
   if (is_collection) {
      Error("CheckInputData", "You are using raw AliEn collections as input. Cannot process run numbers.");
      return kFALSE;   
   }
   
   if (checked && !is_xml) {
      Error("CheckInputData", "Cannot mix processing of full runs with non-xml files");
      return kFALSE;   
   }
   // Check validity of run number(s)
   TObjArray *arr;
   TObjString *os;
   Int_t nruns = 0;
   TString schunk;
   TString path;
   if (!checked) {
      checked = kTRUE;
      use_tags = fDataPattern.Contains("tag");
      TObject::SetBit(AliAnalysisGrid::kUseTags, use_tags);
   }   
   if (use_tags != fDataPattern.Contains("tag")) {
      Error("CheckInputData", "Cannot mix input files using/not using tags");
      return kFALSE;
   }
   if (fRunNumbers.Length()) {
      Info("CheckDataType", "Using supplied run numbers (run ranges are ignored)");
      arr = fRunNumbers.Tokenize(" ");
      TIter next(arr);
      while ((os=(TObjString*)next())) {
         path = Form("%s/%s ", fGridDataDir.Data(), os->GetString().Data());
         if (!DirectoryExists(path)) {
            Warning("CheckInputData", "Run number %s not found in path: <%s>", os->GetString().Data(), path.Data());
            continue;
         }
         path = Form("%s/%s.xml", workdir.Data(),os->GetString().Data());
         TString msg = "\n#####   file: ";
         msg += path;
         msg += " type: xml_collection;";
         if (use_tags) msg += " using_tags: Yes";
         else          msg += " using_tags: No";
         Info("CheckDataType", msg.Data());
         if (fNrunsPerMaster<2) {
            AddDataFile(Form("%s.xml", os->GetString().Data()));
         } else {
            nruns++;
            if (((nruns-1)%fNrunsPerMaster) == 0) {
               schunk = os->GetString();
            }   
            if ((nruns%fNrunsPerMaster)!=0 && os!=arr->Last()) continue;
            schunk += Form("_%s.xml", os->GetString().Data());
            AddDataFile(schunk);
         }   
      }
      delete arr;   
   } else {
      Info("CheckDataType", "Using run range [%d, %d]", fRunRange[0], fRunRange[1]);
      for (Int_t irun=fRunRange[0]; irun<=fRunRange[1]; irun++) {
         path = Form("%s/%d ", fGridDataDir.Data(), irun);
         if (!DirectoryExists(path)) {
//            Warning("CheckInputData", "Run number %d not found in path: <%s>", irun, path.Data());
            continue;
         }
         path = Form("%s/%d.xml", workdir.Data(),irun);
         TString msg = "\n#####   file: ";
         msg += path;
         msg += " type: xml_collection;";
         if (use_tags) msg += " using_tags: Yes";
         else          msg += " using_tags: No";
         Info("CheckDataType", msg.Data());
         if (fNrunsPerMaster<2) {
            AddDataFile(Form("%d.xml",irun));
         } else {
            nruns++;
            if (((nruns-1)%fNrunsPerMaster) == 0) {
               schunk = Form("%d", irun);
            }
            if ((nruns%fNrunsPerMaster)!=0 && irun != fRunRange[1]) continue;
            schunk += Form("_%d.xml",  irun);
            AddDataFile(schunk);
         }   
      }
   }
   return kTRUE;      
}   

//______________________________________________________________________________
Bool_t AliAnalysisAlien::CreateDataset(const char *pattern)
{
// Create dataset for the grid data directory + run number.
   if (TestBit(AliAnalysisGrid::kOffline)) return kFALSE;
   if (!Connect()) {
      Error("CreateDataset", "Cannot create dataset with no grid connection");
      return kFALSE;
   }   

   // Cd workspace
   CdWork();
   TString workdir = gGrid->GetHomeDirectory();
   workdir += fGridWorkingDir;

   // Compose the 'find' command arguments
   TString command;
   TString options = "-x collection ";
   if (TestBit(AliAnalysisGrid::kTest)) options += Form("-l %d ", fNtestFiles);
   TString conditions = "";
   
   TString file;
   TString path;
   Int_t nruns = 0;
   TString schunk;
   TGridCollection *cbase=0, *cadd=0;
   if (!fRunNumbers.Length() && !fRunRange[0]) {
      if (fInputFiles && fInputFiles->GetEntries()) return kTRUE;
      // Make a single data collection from data directory.
      path = fGridDataDir;
      if (!DirectoryExists(path)) {
         Error("CreateDataset", "Path to data directory %s not valid",fGridDataDir.Data());
         return kFALSE;
      }   
//      CdWork();
      if (TestBit(AliAnalysisGrid::kTest)) file = "wn.xml";
      else file = Form("%s.xml", gSystem->BaseName(path));
      if (gSystem->AccessPathName(file)) {
         command = "find ";
         command += options;
         command += path;
         command += " ";
         command += pattern;
         command += conditions;
         Printf("command: %s", command.Data());
         TGridResult *res = gGrid->Command(command);
         if (res) delete res;
         // Write standard output to file
         gROOT->ProcessLine(Form("gGrid->Stdout(); > %s", file.Data()));
      }   
      if (!TestBit(AliAnalysisGrid::kTest) && !FileExists(file)) {
         // Copy xml file to alien space
         TFile::Cp(Form("file:%s",file.Data()), Form("alien://%s/%s",workdir.Data(), file.Data()));
         if (!FileExists(file)) {
            Error("CreateDataset", "Command %s did NOT succeed", command.Data());
            return kFALSE;
         }
         // Update list of files to be processed.
      }
      AddDataFile(Form("%s/%s", workdir.Data(), file.Data()));
      return kTRUE;
   }   
   // Several runs
   if (fRunNumbers.Length()) {
      TObjArray *arr = fRunNumbers.Tokenize(" ");
      TObjString *os;
      TIter next(arr);
      while ((os=(TObjString*)next())) {
         path = Form("%s/%s ", fGridDataDir.Data(), os->GetString().Data());
         if (!DirectoryExists(path)) continue;
//         CdWork();
         if (TestBit(AliAnalysisGrid::kTest)) file = "wn.xml";
         else file = Form("%s.xml", os->GetString().Data());
         // If local collection file does not exist, create it via 'find' command.
         if (gSystem->AccessPathName(file)) {
            command = "find ";
            command += options;
            command += path;
            command += pattern;
            command += conditions;
            TGridResult *res = gGrid->Command(command);
            if (res) delete res;
            // Write standard output to file
            gROOT->ProcessLine(Form("gGrid->Stdout(); > %s", file.Data()));
         }   
         if (TestBit(AliAnalysisGrid::kTest)) break;
         // Check if there is one run per master job.
         if (fNrunsPerMaster<2) {
            if (FileExists(file)) {
               Info("CreateDataset", "\n#####   Dataset %s exist. Skipping creation...", file.Data());
               continue;
            }        
            // Copy xml file to alien space
            TFile::Cp(Form("file:%s",file.Data()), Form("alien://%s/%s",workdir.Data(), file.Data()));
            if (!FileExists(file)) {
               Error("CreateDataset", "Command %s did NOT succeed", command.Data());
               delete arr;
               return kFALSE;
            }
         } else {
            nruns++;
            if (((nruns-1)%fNrunsPerMaster) == 0) {
               schunk = os->GetString();
               cbase = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"%s\", 1000000);",file.Data()));
            } else {
               cadd = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"%s\", 1000000);",file.Data()));
               Printf("   Merging collection <%s> into masterjob input...", file.Data());
               cbase->Add(cadd);
               delete cadd;
            }
            if ((nruns%fNrunsPerMaster)!=0 && os!=arr->Last()) {
               continue;
            }   
            schunk += Form("_%s.xml", os->GetString().Data());
            if (FileExists(schunk)) {
               Info("CreateDataset", "\n#####   Dataset %s exist. Skipping creation...", schunk.Data());
               continue;
            }        
            Printf("Exporting merged collection <%s> and copying to AliEn.", schunk.Data());
            cbase->ExportXML(Form("file://%s", schunk.Data()),kFALSE,kFALSE, schunk, "Merged runs");
            TFile::Cp(Form("file:%s",file.Data()), Form("alien://%s/%s",workdir.Data(), file.Data()));
            if (!FileExists(schunk)) {
               Error("CreateDataset", "Copy command did NOT succeed for %s", schunk.Data());
               delete arr;
               return kFALSE;
            }
         }   
      }   
      delete arr;
   } else {
      // Process a full run range.
      for (Int_t irun=fRunRange[0]; irun<=fRunRange[1]; irun++) {
         path = Form("%s/%d ", fGridDataDir.Data(), irun);
         if (!DirectoryExists(path)) continue;
//         CdWork();
         if (TestBit(AliAnalysisGrid::kTest)) file = "wn.xml";
         else file = Form("%d.xml", irun);
         if (FileExists(file) && fNrunsPerMaster<2 && !TestBit(AliAnalysisGrid::kTest)) {
            Info("CreateDataset", "\n#####   Dataset %s exist. Skipping creation...", file.Data());
//            gGrid->Rm(file); 
            continue;
         }
         // If local collection file does not exist, create it via 'find' command.
         if (gSystem->AccessPathName(file)) {
            command = "find ";
            command += options;
            command += path;
            command += pattern;
            command += conditions;
            TGridResult *res = gGrid->Command(command);
            if (res) delete res;
            // Write standard output to file
            gROOT->ProcessLine(Form("gGrid->Stdout(); > %s", file.Data()));
         }   
         if (TestBit(AliAnalysisGrid::kTest)) break;
         // Check if there is one run per master job.
         if (fNrunsPerMaster<2) {
            if (FileExists(file)) {
               Info("CreateDataset", "\n#####   Dataset %s exist. Skipping creation...", file.Data());
               continue;
            }        
            // Copy xml file to alien space
            TFile::Cp(Form("file:%s",file.Data()), Form("alien://%s/%s",workdir.Data(), file.Data()));
            if (!FileExists(file)) {
               Error("CreateDataset", "Command %s did NOT succeed", command.Data());
               return kFALSE;
            }
         } else {
            nruns++;
            // Check if the collection for the chunk exist locally.
            Int_t nchunk = (nruns-1)/fNrunsPerMaster;
            if (FileExists(fInputFiles->At(nchunk)->GetName())) continue;
            Printf("   Merging collection <%s> into %d runs chunk...",file.Data(),fNrunsPerMaster);
            if (((nruns-1)%fNrunsPerMaster) == 0) {
               schunk = Form("%d", irun);
               cbase = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"%s\", 1000000);",file.Data()));
            } else {
               cadd = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"%s\", 1000000);",file.Data()));
               cbase->Add(cadd);
               delete cadd;
            }
            if ((nruns%fNrunsPerMaster)!=0 && irun!=fRunRange[1]) {
               continue;
            }   
            schunk += Form("_%d.xml", irun);
            if (FileExists(schunk)) {
               Info("CreateDataset", "\n#####   Dataset %s exist. Skipping creation...", schunk.Data());
               continue;
            }        
            Printf("Exporting merged collection <%s> and copying to AliEn.", schunk.Data());
            cbase->ExportXML(Form("file://%s", schunk.Data()),kFALSE,kFALSE, schunk, "Merged runs");
            if (FileExists(schunk)) {
               Info("CreateDataset", "\n#####   Dataset %s exist. Skipping copy...", schunk.Data());
               continue;
            }   
            TFile::Cp(Form("file:%s",schunk.Data()), Form("alien://%s/%s",workdir.Data(), schunk.Data()));
            if (!FileExists(schunk)) {
               Error("CreateDataset", "Copy command did NOT succeed for %s", schunk.Data());
               return kFALSE;
            }
         }   
      }
   }      
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::CreateJDL()
{
// Generate a JDL file according to current settings. The name of the file is 
// specified by fJDLName.
   Bool_t error = kFALSE;
   TObjArray *arr = 0;
   Bool_t copy = kTRUE;
   if (TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   Bool_t generate = kTRUE;
   if (TestBit(AliAnalysisGrid::kTest) || TestBit(AliAnalysisGrid::kSubmit)) generate = kFALSE;
   if (!Connect()) {
      Error("CreateJDL", "Alien connection required");
      return kFALSE;
   }   
   // Check validity of alien workspace
   CdWork();
   TString workdir = gGrid->GetHomeDirectory();
   workdir += fGridWorkingDir;
   if (generate) {
      TObjString *os;
      if (!fInputFiles) {
         Error("CreateJDL()", "Define some input files for your analysis.");
         error = kTRUE;
      }
      // Compose list of input files   
      // Check if output files were defined
      if (!fOutputFiles.Length()) {
         Error("CreateJDL", "You must define at least one output file");
         error = kTRUE;
      }   
      // Check if an output directory was defined and valid
      if (!fGridOutputDir.Length()) {
         Error("CreateJDL", "You must define AliEn output directory");
         error = kTRUE;
      } else {
         if (!fGridOutputDir.Contains("/")) fGridOutputDir = Form("%s/%s", workdir.Data(), fGridOutputDir.Data());
         if (!DirectoryExists(fGridOutputDir)) {
            if (gGrid->Mkdir(fGridOutputDir)) {
               Info("CreateJDL", "\n#####   Created alien output directory %s", fGridOutputDir.Data());
            } else {
               Error("CreateJDL", "Could not create alien output directory %s", fGridOutputDir.Data());
               // error = kTRUE;
            }
         }
         gGrid->Cd(workdir);
      }   
      // Exit if any error up to now
      if (error) return kFALSE;   
      // Set JDL fields
      fGridJDL->SetValue("User", Form("\"%s\"", fUser.Data()));
      fGridJDL->SetExecutable(fExecutable);
//      fGridJDL->SetTTL((UInt_t)fTTL);
      fGridJDL->SetValue("TTL", Form("\"%d\"", fTTL));
      if (fMaxInitFailed > 0) 
         fGridJDL->SetValue("MaxInitFailed", Form("\"%d\"",fMaxInitFailed));
      if (fSplitMaxInputFileNumber > 0) 
         fGridJDL->SetValue("SplitMaxInputFileNumber", Form("\"%d\"", fSplitMaxInputFileNumber));
      if (fSplitMode.Length()) 
         fGridJDL->SetValue("Split", Form("\"%s\"", fSplitMode.Data()));
//      fGridJDL->SetSplitMode(fSplitMode, (UInt_t)fSplitMaxInputFileNumber);
      if (fAliROOTVersion.Length())  
         fGridJDL->AddToPackages("AliRoot", fAliROOTVersion);
      if (fROOTVersion.Length()) 
         fGridJDL->AddToPackages("ROOT", fROOTVersion);
      if (fAPIVersion.Length()) 
         fGridJDL->AddToPackages("APISCONFIG", fAPIVersion);
      fGridJDL->SetInputDataListFormat(fInputFormat);
      fGridJDL->SetInputDataList("wn.xml");
      fGridJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(), fAnalysisMacro.Data()));
      TString analysisFile = fExecutable;
      analysisFile.ReplaceAll(".sh", ".root");
      fGridJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(),analysisFile.Data()));
      if (IsUsingTags() && !gSystem->AccessPathName("ConfigureCuts.C"))
         fGridJDL->AddToInputSandbox(Form("LF:%s/ConfigureCuts.C", workdir.Data()));
      if (fAdditionalLibs.Length()) {
         arr = fAdditionalLibs.Tokenize(" ");
         TIter next(arr);
         while ((os=(TObjString*)next())) {
            if (os->GetString().Contains(".so")) continue;
            fGridJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(), os->GetString().Data()));
         }   
         delete arr;   
      }
      if (fPackages) {
         TIter next(fPackages);
         TObject *obj;
         while ((obj=next()))
            fGridJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(), obj->GetName()));
      }
      if (fOutputArchive.Length()) {
         arr = fOutputArchive.Tokenize(" ");
         TIter next(arr);
         while ((os=(TObjString*)next()))
         if (!os->GetString().Contains("@") && fCloseSE.Length())
            fGridJDL->AddToOutputArchive(Form("%s@%s",os->GetString().Data(), fCloseSE.Data())); 
         else
            fGridJDL->AddToOutputArchive(os->GetString());
         delete arr;
      }      
      arr = fOutputFiles.Tokenize(" ");
      TIter next(arr);
      while ((os=(TObjString*)next())) {
         // Ignore ouputs in jdl that are also in outputarchive
         TString sout = os->GetString();
         if (sout.Index("@")>0) sout.Remove(sout.Index("@"));
         if (fOutputArchive.Contains(sout)) continue;
         if (!os->GetString().Contains("@") && fCloseSE.Length())
            fGridJDL->AddToOutputSandbox(Form("%s@%s",os->GetString().Data(), fCloseSE.Data())); 
         else
            fGridJDL->AddToOutputSandbox(os->GetString());
      }   
      delete arr;
//      fGridJDL->SetPrice((UInt_t)fPrice);
      fGridJDL->SetValue("Price", Form("\"%d\"", fPrice));
      TString validationScript = fExecutable;
      validationScript.ReplaceAll(".sh", "_validation.sh");
      fGridJDL->SetValidationCommand(Form("%s/%s", workdir.Data(),validationScript.Data()));
      if (fMasterResubmitThreshold) fGridJDL->SetValue("MasterResubmitThreshold", Form("\"%d%%\"", fMasterResubmitThreshold));
      // Write a jdl with 2 input parameters: collection name and output dir name.
      WriteJDL(copy);
   }
   // Copy jdl to grid workspace   
   if (copy) {
      if (fAdditionalLibs.Length()) {
         arr = fAdditionalLibs.Tokenize(" ");
         TObjString *os;
         TIter next(arr);
         while ((os=(TObjString*)next())) {
            if (os->GetString().Contains(".so")) continue;
            Info("CreateJDL", "\n#####   Copying dependency: <%s> to your alien workspace", os->GetString().Data());
            if (FileExists(os->GetString())) gGrid->Rm(os->GetString());
            TFile::Cp(Form("file:%s",os->GetString().Data()), Form("alien://%s/%s", workdir.Data(), os->GetString().Data()));
         }   
         delete arr;   
      }
      if (fPackages) {
         TIter next(fPackages);
         TObject *obj;
         while ((obj=next())) {
            Info("CreateJDL", "\n#####   Copying dependency: <%s> to your alien workspace", obj->GetName());
            TFile::Cp(Form("file:%s",obj->GetName()), Form("alien://%s/%s", workdir.Data(), obj->GetName()));
         }   
      }      
   } 
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::WriteJDL(Bool_t copy)
{
// Writes one or more JDL's corresponding to findex. If findex is negative,
// all run numbers are considered in one go (jdl). For non-negative indices
// they correspond to the indices in the array fInputFiles.
   if (!fInputFiles) return kFALSE;
   TObjString *os;
   TString workdir = gGrid->GetHomeDirectory();
   workdir += fGridWorkingDir;
   
   if (!fRunNumbers.Length() && !fRunRange[0]) {
      // One jdl with no parameters in case input data is specified by name.
      TIter next(fInputFiles);
      while ((os=(TObjString*)next()))
         fGridJDL->AddToInputDataCollection(Form("LF:%s,nodownload", os->GetString().Data()));
      fGridJDL->SetOutputDirectory(Form("%s/#alien_counter_03i#", fGridOutputDir.Data()));
   } else {
      // One jdl to be submitted with 2 input parameters: data collection name and output dir prefix
      fGridJDL->AddToInputDataCollection(Form("LF:%s/$1,nodownload", workdir.Data()));
      fGridJDL->SetOutputDirectory(Form("%s/$2/#alien_counter_03i#", fGridOutputDir.Data()));
   }
      

   // Generate the JDL as a string
   TString sjdl = fGridJDL->Generate();
   Int_t index;
   index = sjdl.Index("Executable");
   if (index >= 0) sjdl.Insert(index, "\n# This is the startup script\n");
   index = sjdl.Index("Split ");
   if (index >= 0) sjdl.Insert(index, "\n# We split per storage element\n");
   index = sjdl.Index("SplitMaxInputFileNumber");
   if (index >= 0) sjdl.Insert(index, "\n# We want each subjob to get maximum this number of input files\n");
   index = sjdl.Index("InputDataCollection");
   if (index >= 0) sjdl.Insert(index, "# Input xml collections\n");
   index = sjdl.Index("InputFile");
   if (index >= 0) sjdl.Insert(index, "\n# List of input files to be uploaded to wn's\n");
   index = sjdl.Index("InputDataList ");
   if (index >= 0) sjdl.Insert(index, "\n# Collection to be processed on wn\n");
   index = sjdl.Index("InputDataListFormat");
   if (index >= 0) sjdl.Insert(index, "\n# Format of input data\n");
   index = sjdl.Index("Price");
   if (index >= 0) sjdl.Insert(index, "\n# AliEn price for this job\n");
   index = sjdl.Index("Requirements");
   if (index >= 0) sjdl.Insert(index, "\n# Additional requirements for the computing element\n");
   index = sjdl.Index("Packages");
   if (index >= 0) sjdl.Insert(index, "\n# Packages to be used\n");
   index = sjdl.Index("User =");
   if (index >= 0) sjdl.Insert(index, "\n# AliEn user\n");
   index = sjdl.Index("TTL");
   if (index >= 0) sjdl.Insert(index, "\n# Time to live for the job\n");
   index = sjdl.Index("OutputFile");
   if (index >= 0) sjdl.Insert(index, "\n# List of output files to be registered\n");
   index = sjdl.Index("OutputDir");
   if (index >= 0) sjdl.Insert(index, "\n# Output directory\n");
   index = sjdl.Index("OutputArchive");
   if (index >= 0) sjdl.Insert(index, "\n# Files to be archived\n");
   index = sjdl.Index("MaxInitFailed");
   if (index >= 0) sjdl.Insert(index, "\n# Maximum number of first failing jobs to abort the master job\n");
   index = sjdl.Index("MasterResubmitThreshold");
   if (index >= 0) sjdl.Insert(index, "\n# Resubmit failed jobs until DONE rate reaches this percentage\n");
   sjdl.ReplaceAll("ValidationCommand", "Validationcommand");
   index = sjdl.Index("Validationcommand");
   if (index >= 0) sjdl.Insert(index, "\n# Validation script to be run for each subjob\n");
   sjdl.ReplaceAll("\"LF:", "\n   \"LF:");
   sjdl.ReplaceAll("(member", "\n   (member");
   sjdl.ReplaceAll("\",\"VO_", "\",\n   \"VO_");
   sjdl.ReplaceAll("{", "{\n   ");
   sjdl.ReplaceAll("};", "\n};");
   sjdl.ReplaceAll("{\n   \n", "{\n");
   sjdl.ReplaceAll("\n\n", "\n");
   sjdl.ReplaceAll("OutputDirectory", "OutputDir");
   sjdl += "JDLVariables = \n{\n   \"Packages\",\n   \"OutputDir\"\n};\n";
   sjdl.Prepend(Form("Jobtag = {\n   \"comment:%s\"\n};\n", fJobTag.Data()));
   index = sjdl.Index("JDLVariables");
   if (index >= 0) sjdl.Insert(index, "\n# JDL variables\n");
   // Write jdl to file
   ofstream out;
   out.open(fJDLName.Data(), ios::out);
   if (out.bad()) {
      Error("CreateJDL", "Bad file name: %s", fJDLName.Data());
      return kFALSE;
   }
   out << sjdl << endl;

   // Copy jdl to grid workspace   
   if (!copy) {
      Info("CreateJDL", "\n#####   You may want to review jdl:%s and analysis macro:%s before running in <submit> mode", fJDLName.Data(), fAnalysisMacro.Data());
   } else {
      Info("CreateJDL", "\n#####   Copying JDL file <%s> to your AliEn output directory", fJDLName.Data());
      TString locjdl = Form("%s/%s", fGridOutputDir.Data(),fJDLName.Data());
      if (fProductionMode)
         locjdl = Form("%s/%s", workdir.Data(),fJDLName.Data());
      if (FileExists(locjdl)) gGrid->Rm(locjdl);
      TFile::Cp(Form("file:%s",fJDLName.Data()), Form("alien://%s", locjdl.Data()));
   } 
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::FileExists(const char *lfn)
{
// Returns true if file exists.
   if (!gGrid) return kFALSE;
   TGridResult *res = gGrid->Ls(lfn);
   if (!res) return kFALSE;
   TMap *map = dynamic_cast<TMap*>(res->At(0));
   if (!map) {
      delete res;
      return kFALSE;
   }   
   TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("name"));
   if (!objs || !objs->GetString().Length()) {
      delete res;
      return kFALSE;
   }
   delete res;   
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::DirectoryExists(const char *dirname)
{
// Returns true if directory exists. Can be also a path.
   if (!gGrid) return kFALSE;
   // Check if dirname is a path
   TString dirstripped = dirname;
   dirstripped = dirstripped.Strip();
   dirstripped = dirstripped.Strip(TString::kTrailing, '/');
   TString dir = gSystem->BaseName(dirstripped);
   dir += "/";
   TString path = gSystem->DirName(dirstripped);
   TGridResult *res = gGrid->Ls(path, "-F");
   if (!res) return kFALSE;
   TIter next(res);
   TMap *map;
   TObject *obj;
   while ((map=dynamic_cast<TMap*>(next()))) {
      obj = map->GetValue("name");
      if (!obj) break;
      if (dir == obj->GetName()) {
         delete res;
         return kTRUE;
      }
   }
   delete res;
   return kFALSE;
}      

//______________________________________________________________________________
void AliAnalysisAlien::CheckDataType(const char *lfn, Bool_t &is_collection, Bool_t &is_xml, Bool_t &use_tags)
{
// Check input data type.
   is_collection = kFALSE;
   is_xml = kFALSE;
   use_tags = kFALSE;
   if (!gGrid) {
      Error("CheckDataType", "No connection to grid");
      return;
   }
   is_collection = IsCollection(lfn);
   TString msg = "\n#####   file: ";
   msg += lfn;
   if (is_collection) {
      msg += " type: raw_collection;";
   // special treatment for collections
      is_xml = kFALSE;
      // check for tag files in the collection
      TGridResult *res = gGrid->Command(Form("listFilesFromCollection -z -v %s",lfn), kFALSE);
      if (!res) {
         msg += " using_tags: No (unknown)";
         Info("CheckDataType", msg.Data());
         return;
      }   
      const char* typeStr = res->GetKey(0, "origLFN");
      if (!typeStr || !strlen(typeStr)) {
         msg += " using_tags: No (unknown)";
         Info("CheckDataType", msg.Data());
         return;
      }   
      TString file = typeStr;
      use_tags = file.Contains(".tag");
      if (use_tags) msg += " using_tags: Yes";
      else          msg += " using_tags: No";
      Info("CheckDataType", msg.Data());
      return;
   }
   TString slfn(lfn);
   slfn.ToLower();
   is_xml = slfn.Contains(".xml");
   if (is_xml) {
   // Open xml collection and check if there are tag files inside
      msg += " type: xml_collection;";
      TGridCollection *coll = (TGridCollection*)gROOT->ProcessLine(Form("TAlienCollection::Open(\"alien://%s\",1);",lfn));
      if (!coll) {
         msg += " using_tags: No (unknown)";
         Info("CheckDataType", msg.Data());
         return;
      }   
      TMap *map = coll->Next();
      if (!map) {
         msg += " using_tags: No (unknown)";
         Info("CheckDataType", msg.Data());
         return;
      }   
      map = (TMap*)map->GetValue("");
      TString file;
      if (map && map->GetValue("name")) file = map->GetValue("name")->GetName();
      use_tags = file.Contains(".tag");
      delete coll;
      if (use_tags) msg += " using_tags: Yes";
      else          msg += " using_tags: No";
      Info("CheckDataType", msg.Data());
      return;
   }
   use_tags = slfn.Contains(".tag");
   if (slfn.Contains(".root")) msg += " type: root file;";
   else                        msg += " type: unhnown file;";
   if (use_tags) msg += " using_tags: Yes";
   else          msg += " using_tags: No";
   Info("CheckDataType", msg.Data());
}

//______________________________________________________________________________
void AliAnalysisAlien::EnablePackage(const char *package)
{
// Enables a par file supposed to exist in the current directory.
   TString pkg(package);
   pkg.ReplaceAll(".par", "");
   pkg += ".par";
   if (gSystem->AccessPathName(pkg)) {
      Error("EnablePackage", "Package %s not found", pkg.Data());
      return;
   }
   if (!TObject::TestBit(AliAnalysisGrid::kUsePars))
      Info("EnablePackage", "AliEn plugin will use .par packages");
   TObject::SetBit(AliAnalysisGrid::kUsePars, kTRUE);
   if (!fPackages) {
      fPackages = new TObjArray();
      fPackages->SetOwner();
   }
   fPackages->Add(new TObjString(pkg));
}      

//______________________________________________________________________________
const char *AliAnalysisAlien::GetJobStatus(Int_t jobidstart, Int_t lastid, Int_t &nrunning, Int_t &nwaiting, Int_t &nerror, Int_t &ndone)
{
// Get job status for all jobs with jobid>jobidstart.
   static char mstatus[20];
   mstatus[0] = '\0';
   nrunning = 0;
   nwaiting = 0;
   nerror   = 0;
   ndone    = 0;
   TGridJobStatusList *list = gGrid->Ps("");
   if (!list) return mstatus;
   Int_t nentries = list->GetSize();
   TGridJobStatus *status;
   Int_t pid;
   for (Int_t ijob=0; ijob<nentries; ijob++) {
      status = (TGridJobStatus *)list->At(ijob);
      pid = gROOT->ProcessLine(Form("atoi(((TAlienJobStatus*)0x%lx)->GetKey(\"queueId\"));", (ULong_t)status));
      if (pid<jobidstart) continue;
      if (pid == lastid) {
         gROOT->ProcessLine(Form("sprintf((char*)0x%lx,((TAlienJobStatus*)0x%lx)->GetKey(\"status\"));",(ULong_t)mstatus, (ULong_t)status));
      }   
      switch (status->GetStatus()) {
         case TGridJobStatus::kWAITING:
            nwaiting++; break;
         case TGridJobStatus::kRUNNING:
            nrunning++; break;
         case TGridJobStatus::kABORTED:
         case TGridJobStatus::kFAIL:
         case TGridJobStatus::kUNKNOWN:
            nerror++; break;
         case TGridJobStatus::kDONE:
            ndone++;
      }
   }
   list->Delete();
   delete list;
   return mstatus;
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::IsCollection(const char *lfn) const
{
// Returns true if file is a collection. Functionality duplicated from
// TAlien::Type() because we don't want to directly depend on TAlien.
   if (!gGrid) {
      Error("IsCollection", "No connection to grid");
      return kFALSE;
   }
   TGridResult *res = gGrid->Command(Form("type -z %s",lfn),kFALSE);
   if (!res) return kFALSE;
   const char* typeStr = res->GetKey(0, "type");
   if (!typeStr || !strlen(typeStr)) return kFALSE;
   if (!strcmp(typeStr, "collection")) return kTRUE;
   delete res;
   return kFALSE;
}   

//______________________________________________________________________________
void AliAnalysisAlien::Print(Option_t *) const
{
// Print current plugin settings.
   Printf("### AliEn analysis plugin current settings ###");
   Printf("=   Production mode:______________________________ %d", fProductionMode);
   Printf("=   Version of API requested: ____________________ %s", fAPIVersion.Data());
   Printf("=   Version of ROOT requested: ___________________ %s", fROOTVersion.Data());
   Printf("=   Version of AliRoot requested: ________________ %s", fAliROOTVersion.Data());
   if (fUser.Length()) 
   Printf("=   User running the plugin: _____________________ %s", fUser.Data());
   Printf("=   Grid workdir relative to user $HOME: _________ %s", fGridWorkingDir.Data());
   Printf("=   Grid output directory relative to workdir: ___ %s", fGridOutputDir.Data());
   Printf("=   Data base directory path requested: __________ %s", fGridDataDir.Data());
   Printf("=   Data search pattern: _________________________ %s", fDataPattern.Data());
   Printf("=   Input data format: ___________________________ %s", fInputFormat.Data());
   if (fRunNumbers.Length()) 
   Printf("=   Run numbers to be processed: _________________ %s", fRunNumbers.Data());
   if (fRunRange[0])
   Printf("=   Run range to be processed: ___________________ %d-%d", fRunRange[0], fRunRange[1]);
   if (!fRunRange[0] && !fRunNumbers.Length()) {
      TIter next(fInputFiles);
      TObject *obj;
      TString list;
      while ((obj=next())) list += obj->GetName();
      Printf("=   Input files to be processed: _________________ %s", list.Data());
   }
   if (TestBit(AliAnalysisGrid::kTest))
   Printf("=   Number of input files used in test mode: _____ %d", fNtestFiles);
   Printf("=   List of output files to be registered: _______ %s", fOutputFiles.Data());
   Printf("=   List of outputs going to be archived: ________ %s", fOutputArchive.Data());
   Printf("=   List of outputs that should not be merged: ___ %s", fMergeExcludes.Data());
   Printf("=====================================================================");
   Printf("=   Job price: ___________________________________ %d", fPrice);
   Printf("=   Time to live (TTL): __________________________ %d", fTTL);
   Printf("=   Max files per subjob: ________________________ %d", fSplitMaxInputFileNumber);
   if (fMaxInitFailed>0) 
   Printf("=   Max number of subjob fails to kill: __________ %d", fMaxInitFailed);
   if (fMasterResubmitThreshold>0) 
   Printf("=   Resubmit master job if failed subjobs >_______ %d", fMasterResubmitThreshold);
   if (fNrunsPerMaster>0)
   Printf("=   Number of runs per master job: _______________ %d", fNrunsPerMaster);
   Printf("=   Number of files in one chunk to be merged: ___ %d", fMaxMergeFiles);
   Printf("=   Name of the generated execution script: ______ %s",fExecutable.Data());
   if (fArguments.Length()) 
   Printf("=   Arguments for the execution script: __________ %s",fArguments.Data());
   Printf("=   Name of the generated analysis macro: ________ %s",fAnalysisMacro.Data());
   Printf("=   User analysis files to be deployed: __________ %s",fAnalysisSource.Data());
   Printf("=   Additional libs to be loaded or souces to be compiled runtime: <%s>",fAdditionalLibs.Data());
   Printf("=   Master jobs split mode: ______________________ %s",fSplitMode.Data());
   if (fDatasetName)
   Printf("=   Custom name for the dataset to be created: ___ %s", fDatasetName.Data());
   Printf("=   Name of the generated JDL: ___________________ %s", fJDLName.Data());
   if (fIncludePath.Data())
   Printf("=   Include path for runtime task compilation: ___ %s", fIncludePath.Data());
   if (fCloseSE.Length())
   Printf("=   Force job outputs to storage element: ________ %s", fCloseSE.Data());
   if (fFriendChainName.Length())
   Printf("=   Open friend chain file on worker: ____________ %s", fFriendChainName.Data());
   if (fPackages) {
      TIter next(fPackages);
      TObject *obj;
      TString list;
      while ((obj=next())) list += obj->GetName();
      Printf("=   Par files to be used: ________________________ %s", list.Data());
   }   
}

//______________________________________________________________________________
void AliAnalysisAlien::SetDefaults()
{
// Set default values for everything. What cannot be filled will be left empty.
   if (fGridJDL) delete fGridJDL;
   fGridJDL = (TGridJDL*)gROOT->ProcessLine("new TAlienJDL()");
   fPrice                      = 1;
   fTTL                        = 30000;
   fSplitMaxInputFileNumber    = 100;
   fMaxInitFailed              = 0;
   fMasterResubmitThreshold    = 0;
   fNtestFiles                 = 10;
   fRunRange[0]                = 0;
   fRunRange[1]                = 0;
   fNrunsPerMaster             = 1;
   fMaxMergeFiles              = 100;
   fRunNumbers                 = "";
   fExecutable                 = "analysis.sh";
   fExecutableCommand          = "root -b -q";
   fArguments                  = "";
   fAnalysisMacro              = "myAnalysis.C";
   fAnalysisSource             = "";
   fAdditionalLibs             = "";
   fSplitMode                  = "se";
   fAPIVersion                 = "";
   fROOTVersion                = "";
   fAliROOTVersion             = "";
   fUser                       = "";  // Your alien user name
   fGridWorkingDir             = "";
   fGridDataDir                = "";  // Can be like: /alice/sim/PDC_08a/LHC08c9/
   fDataPattern                = "*AliESDs.root";  // Can be like: *AliESDs.root, */pass1/*AliESDs.root, ...
   fFriendChainName            = "";
   fGridOutputDir              = "output";
   fOutputArchive              = "log_archive.zip:stdout,stderr root_archive.zip:*.root";
   fOutputFiles                = "";  // Like "AliAODs.root histos.root"
   fInputFormat                = "xml-single";
   fJDLName                    = "analysis.jdl";
   fJobTag                     = "Automatically generated analysis JDL";
   fMergeExcludes              = "";
}   

//______________________________________________________________________________
Bool_t AliAnalysisAlien::MergeOutputs()
{
// Merge analysis outputs existing in the AliEn space.
   if (TestBit(AliAnalysisGrid::kTest)) return kTRUE;
   if (TestBit(AliAnalysisGrid::kOffline)) return kFALSE;
   if (!Connect()) {
      Error("MergeOutputs", "Cannot merge outputs without grid connection. Terminate will NOT be executed");
      return kFALSE;
   }   
   // Get the output path
   if (!fGridOutputDir.Contains("/")) fGridOutputDir = Form("/%s/%s/%s", gGrid->GetHomeDirectory(), fGridWorkingDir.Data(), fGridOutputDir.Data());
   if (!DirectoryExists(fGridOutputDir)) {
      Error("MergeOutputs", "Grid output directory %s not found. Terminate() will NOT be executed", fGridOutputDir.Data());
      return kFALSE;
   }
   if (!fOutputFiles.Length()) {
      Error("MergeOutputs", "No output file names defined. Are you running the right AliAnalysisAlien configuration ?");
      return kFALSE;
   }   
   TObjArray *list = fOutputFiles.Tokenize(" ");
   TIter next(list);
   TObjString *str;
   TString command;
   TString output_file;
   TString output_chunk;
   TString previous_chunk;
   Int_t count_chunk = 0;
   Int_t count_zero = fMaxMergeFiles;
   Bool_t merged = kTRUE;
   while((str=(TObjString*)next())) {
      output_file = str->GetString();
      Int_t index = output_file.Index("@");
      if (index > 0) output_file.Remove(index);
      // Skip already merged outputs
      if (!gSystem->AccessPathName(output_file)) {
         Info("MergeOutputs", "Output file <%s> found. Not merging again.", output_file.Data());
         continue;
      }   
      if (fMergeExcludes.Length() &&
          fMergeExcludes.Contains(output_file.Data())) continue;
      // Perform a 'find' command in the output directory, looking for registered outputs    
      command = Form("find %s/ *%s", fGridOutputDir.Data(), output_file.Data());
      Printf("command: %s", command.Data());
      TGridResult *res = gGrid->Command(command);
      if (!res) continue;
      TFileMerger *fm = 0;
      TIter nextmap(res);
      TMap *map = 0;
      previous_chunk = "";
      count_chunk = 0;
      // Check if there is a merge operation to resume
      output_chunk = output_file;
      output_chunk.ReplaceAll(".root", "_*.root");
      if (!gSystem->Exec(Form("ls %s", output_chunk.Data()))) {
         while (1) {
            for (Int_t counter=0; counter<fMaxMergeFiles; counter++) map = (TMap*)nextmap();
            if (!map) {
               Error("MergeOutputs", "Cannot resume merging for <%s>, nentries=%d", output_file.Data(), res->GetSize());
               delete res;
               return kFALSE;
            }
            output_chunk = output_file;
            output_chunk.ReplaceAll(".root", Form("_%04d.root", count_chunk));
            Printf("%s", output_chunk.Data());
            count_chunk++;
            if (gSystem->AccessPathName(output_chunk)) continue;
            // Merged file with chunks up to <count_chunk> found
            Printf("Resume merging of <%s> from <%s>", output_file.Data(), output_chunk.Data());
            previous_chunk = output_chunk;
            break;
         }
      }
      count_zero = fMaxMergeFiles;
      while ((map=(TMap*)nextmap())) {
      // Loop 'find' results and get next LFN
         if (count_zero == fMaxMergeFiles) {
            // First file in chunk - create file merger and add previous chunk if any.
            fm = new TFileMerger(kFALSE);
            fm->SetFastMethod(kTRUE);
            if (previous_chunk.Length()) fm->AddFile(previous_chunk.Data());
            output_chunk = output_file;
            output_chunk.ReplaceAll(".root", Form("_%04d.root", count_chunk));
         }
         // If last file found, put merged results in the output file
         if (map == res->Last()) output_chunk = output_file;
         TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
         if (!objs || !objs->GetString().Length()) {
            // Nothing found - skip this output
            delete res;
            delete fm;
            break;
         } 
         // Add file to be merged and decrement chunk counter.
         fm->AddFile(objs->GetString());
         count_zero--;
         if (count_zero==0 || map == res->Last()) {            
            fm->OutputFile(output_chunk);
            if (!fm->GetMergeList() || !fm->GetMergeList()->GetSize()) {
            // Nothing found - skip this output
               Warning("MergeOutputs", "No <%s> files found.", output_file.Data());
               delete res;
               delete fm;
               break;
            }
            // Merge the outputs, then go to next chunk      
            if (!fm->Merge()) {
               Error("MergeOutputs", "Could not merge all <%s> files", output_file.Data());
               delete res;
               delete fm;
               merged = kFALSE;
               break;
            } else {
               Info("MergeOutputs", "\n#####   Merged %d output files to <%s>", fm->GetMergeList()->GetSize(), output_chunk.Data());
               gSystem->Unlink(previous_chunk);
            }
            if (map == res->Last()) {
               delete res;
               delete fm;
               break;
            }      
            count_chunk++;
            count_zero = fMaxMergeFiles;
            previous_chunk = output_chunk;
         }
      }
   } 
   if (!merged) {
      Error("MergeOutputs", "Terminate() will  NOT be executed");
   }  
   return merged;
}   

//______________________________________________________________________________
void AliAnalysisAlien::SetDefaultOutputs(Bool_t flag)
{
// Use the output files connected to output containers from the analysis manager
// rather than the files defined by SetOutputFiles
   if (flag && !TObject::TestBit(AliAnalysisGrid::kDefaultOutputs))
      Info("SetDefaultOutputs", "Plugin will use the output files taken from \
      analysis manager");
   TObject::SetBit(AliAnalysisGrid::kDefaultOutputs, flag);
}
      
//______________________________________________________________________________
Bool_t AliAnalysisAlien::StartAnalysis(Long64_t /*nentries*/, Long64_t /*firstEntry*/)
{
// Start remote grid analysis.
   
   // Check if output files have to be taken from the analysis manager
   if (TestBit(AliAnalysisGrid::kDefaultOutputs)) {
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      if (!mgr || !mgr->IsInitialized()) {
         Error("StartAnalysis", "You need an initialized analysis manager for this");
         return kFALSE;
      }
      fOutputFiles = "";
      TIter next(mgr->GetOutputs());
      AliAnalysisDataContainer *output;
      while ((output=(AliAnalysisDataContainer*)next())) {
         const char *filename = output->GetFileName();
         if (!(strcmp(filename, "default"))) {
            if (!mgr->GetOutputEventHandler()) continue;
            filename = mgr->GetOutputEventHandler()->GetOutputFileName();
         }
         if (fOutputFiles.Length()) fOutputFiles += " ";
         fOutputFiles += filename;
      }
      // Add extra files registered to the analysis manager
      if (mgr->GetExtraFiles().Length()) {
         if (fOutputFiles.Length()) fOutputFiles += " ";
         fOutputFiles += mgr->GetExtraFiles();
      }   
   }
//   if (!fCloseSE.Length()) fCloseSE = gSystem->Getenv("alien_CLOSE_SE");
   if (TestBit(AliAnalysisGrid::kOffline)) {
      Info("StartAnalysis","\n##### OFFLINE MODE ##### Files to be used in GRID are produced but not copied \
      \n                         there nor any job run. You can revise the JDL and analysis \
      \n                         macro then run the same in \"submit\" mode.");
   } else if (TestBit(AliAnalysisGrid::kTest)) {
      Info("StartAnalysis","\n##### LOCAL MODE #####   Your analysis will be run locally on a subset of the requested \
      \n                         dataset.");
   } else if (TestBit(AliAnalysisGrid::kSubmit)) {
      Info("StartAnalysis","\n##### SUBMIT MODE #####  Files required by your analysis are copied to your grid working \
      \n                         space and job submitted.");
   } else if (TestBit(AliAnalysisGrid::kMerge)) {
      Info("StartAnalysis","\n##### MERGE MODE #####   The registered outputs of the analysis will be merged");
      return kTRUE;
   } else {
      Info("StartAnalysis","\n##### FULL ANALYSIS MODE ##### Producing needed files and submitting your analysis job...");   
   }   
      
   if (!Connect()) {
      Error("StartAnalysis", "Cannot start grid analysis without grid connection");
      return kFALSE;
   }
   Print();   
   if (!CheckInputData()) {
      Error("StartAnalysis", "There was an error in preprocessing your requested input data");
      return kFALSE;
   }   
   CreateDataset(fDataPattern);
   WriteAnalysisFile();   
   WriteAnalysisMacro();
   WriteExecutable();
   WriteValidationScript();
   if (!CreateJDL()) return kFALSE;
   if (TestBit(AliAnalysisGrid::kOffline)) return kFALSE;
   if (TestBit(AliAnalysisGrid::kTest)) {
      // Locally testing the analysis
      Info("StartAnalysis", "\n_______________________________________________________________________ \
      \n   Running analysis script in a daughter shell as on a worker node \
      \n_______________________________________________________________________");
      TObjArray *list = fOutputFiles.Tokenize(" ");
      TIter next(list);
      TObjString *str;
      TString output_file;
      while((str=(TObjString*)next())) {
         output_file = str->GetString();
         Int_t index = output_file.Index("@");
         if (index > 0) output_file.Remove(index);         
         if (!gSystem->AccessPathName(output_file)) gSystem->Exec(Form("rm %s", output_file.Data()));
      }
      delete list;
      gSystem->Exec(Form("bash %s 2>stderr", fExecutable.Data()));
      TString validationScript = fExecutable;
      validationScript.ReplaceAll(".sh", "_validation.sh");
      gSystem->Exec(Form("bash %s",validationScript.Data()));
//      gSystem->Exec("cat stdout");
      return kFALSE;
   }
   // Check if submitting is managed by LPM manager
   if (fProductionMode) {
      TString prodfile = fJDLName;
      prodfile.ReplaceAll(".jdl", ".prod");
      WriteProductionFile(prodfile);
      Info("StartAnalysis", "Job submitting is managed by LPM. Rerun in terminate mode after jobs finished.");
      return kFALSE;
   }   
   // Submit AliEn job(s)
   gGrid->Cd(fGridOutputDir);
   TGridResult *res;
   TString jobID = "";
   if (!fRunNumbers.Length() && !fRunRange[0]) {
      // Submit a given xml or a set of runs
      res = gGrid->Command(Form("submit %s", fJDLName.Data()));
      Printf("*************************** %s",Form("submit %s", fJDLName.Data()));
      if (res) {
         const char *cjobId = res->GetKey(0,"jobId");
         if (!cjobId) {
            Error("StartAnalysis", "Your JDL %s could not be submitted", fJDLName.Data());
            return kFALSE;
         } else {
            Info("StartAnalysis", "\n_______________________________________________________________________ \
            \n#####   Your JDL %s was successfully submitted. \nTHE JOB ID IS: %s \
            \n_______________________________________________________________________",
                   fJDLName.Data(), cjobId);
            jobID = cjobId;      
         }          
         delete res;
      }   
   } else {
      // Submit for a range of enumeration of runs.
      Submit();
   }   
         
   Info("StartAnalysis", "\n#### STARTING AN ALIEN SHELL FOR YOU. EXIT WHEN YOUR JOB %s HAS FINISHED. #### \
   \n You may exit at any time and terminate the job later using the option <terminate> \
   \n ##################################################################################", jobID.Data());
   gSystem->Exec("aliensh");
   return kTRUE;
}

//______________________________________________________________________________
void AliAnalysisAlien::Submit()
{
// Submit all master jobs.
   Int_t nmasterjobs = fInputFiles->GetEntries();
   Long_t tshoot = gSystem->Now();
   if (!fNsubmitted) SubmitNext();
   while (fNsubmitted < nmasterjobs) {
      Long_t now = gSystem->Now();
      if ((now-tshoot)>30000) {
         tshoot = now;
         SubmitNext();
      }   
   }
}

//______________________________________________________________________________
void AliAnalysisAlien::SubmitNext()
{
// Submit next bunch of master jobs if the queue is free.
   static Bool_t iscalled = kFALSE;
   static Int_t firstmaster = 0;
   static Int_t lastmaster = 0;
   static Int_t npermaster  = 0;
   if (iscalled) return;
   iscalled = kTRUE;
   Int_t nrunning=0, nwaiting=0, nerror=0, ndone=0;
   Int_t ntosubmit = 0;
   TGridResult *res;
   TString jobID = "";
   if (!fNsubmitted) ntosubmit = 1;
   else {
      TString status = GetJobStatus(firstmaster, lastmaster, nrunning, nwaiting, nerror, ndone);
      Printf("=== master %d: %s", lastmaster, status.Data());
      // If last master not split, just return
      if (status != "SPLIT") {iscalled = kFALSE; return;}
      // No more than 100 waiting jobs
      if (nwaiting>100) {iscalled = kFALSE; return;}
      npermaster = (nrunning+nwaiting+nerror+ndone)/fNsubmitted;      
      if (npermaster) ntosubmit = (100-nwaiting)/npermaster;
      Printf("=== WAITING(%d) RUNNING(%d) DONE(%d) OTHER(%d) NperMaster=%d => to submit %d jobs", 
             nwaiting, nrunning, ndone, nerror, npermaster, ntosubmit);
   }
   Int_t nmasterjobs = fInputFiles->GetEntries();
   for (Int_t i=0; i<ntosubmit; i++) {
      // Submit for a range of enumeration of runs.
      if (fNsubmitted>=nmasterjobs) {iscalled = kFALSE; return;}
      TString query;
      query = Form("submit %s %s %03d", fJDLName.Data(), fInputFiles->At(fNsubmitted)->GetName(), fNsubmitted);
      Printf("********* %s",query.Data());
      res = gGrid->Command(query);
      if (res) {
         TString cjobId1 = res->GetKey(0,"jobId");
         if (!cjobId1.Length()) {
            Error("StartAnalysis", "Your JDL %s could not be submitted", fJDLName.Data());
            iscalled = kFALSE;
            return;
         } else {
            Info("StartAnalysis", "\n_______________________________________________________________________ \
            \n#####   Your JDL %s submitted (%d to go). \nTHE JOB ID IS: %s \
            \n_______________________________________________________________________",
                fJDLName.Data(), nmasterjobs-fNsubmitted-1, cjobId1.Data());
            jobID += cjobId1;
            jobID += " ";
            lastmaster = cjobId1.Atoi();
            if (!firstmaster) firstmaster = lastmaster;
            fNsubmitted++;
         }          
         delete res;
      }   
   }
   iscalled = kFALSE;
}

//______________________________________________________________________________
void AliAnalysisAlien::WriteAnalysisFile()
{
// Write current analysis manager into the file <analysisFile>
   TString analysisFile = fExecutable;
   analysisFile.ReplaceAll(".sh", ".root");
   if (!TestBit(AliAnalysisGrid::kSubmit)) {  
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      if (!mgr || !mgr->IsInitialized()) {
         Error("WriteAnalysisFile", "You need an initialized analysis manager for this");
         return;
      }
      // Check analysis type
      TObject *handler;
      if (mgr->GetMCtruthEventHandler()) TObject::SetBit(AliAnalysisGrid::kUseMC);
      handler = (TObject*)mgr->GetInputEventHandler();
      if (handler) {
         if (handler->InheritsFrom("AliESDInputHandler")) TObject::SetBit(AliAnalysisGrid::kUseESD);
         if (handler->InheritsFrom("AliAODInputHandler")) TObject::SetBit(AliAnalysisGrid::kUseAOD);
      }
      TDirectory *cdir = gDirectory;
      TFile *file = TFile::Open(analysisFile, "RECREATE");
      if (file) {
         mgr->Write();
         delete file;
      }
      if (cdir) cdir->cd();
      Info("WriteAnalysisFile", "\n#####   Analysis manager: %s wrote to file <%s>\n", mgr->GetName(),analysisFile.Data());
   }   
   Bool_t copy = kTRUE;
   if (TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      workdir += fGridWorkingDir;
      Info("CreateJDL", "\n#####   Copying file <%s> containing your initialized analysis manager to your alien workspace", analysisFile.Data());
      if (FileExists(analysisFile)) gGrid->Rm(analysisFile);
      TFile::Cp(Form("file:%s",analysisFile.Data()), Form("alien://%s/%s", workdir.Data(),analysisFile.Data()));
   }   
}

//______________________________________________________________________________
void AliAnalysisAlien::WriteAnalysisMacro()
{
// Write the analysis macro that will steer the analysis in grid mode.
   if (!TestBit(AliAnalysisGrid::kSubmit)) {  
      ofstream out;
      out.open(fAnalysisMacro.Data(), ios::out);
      if (!out.good()) {
         Error("WriteAnalysisMacro", "could not open file %s for writing", fAnalysisMacro.Data());
         return;
      }
      TString func = fAnalysisMacro;
      TString type = "ESD";
      TString comment = "// Analysis using ";
      if (TObject::TestBit(AliAnalysisGrid::kUseESD)) comment += "ESD";
      if (TObject::TestBit(AliAnalysisGrid::kUseAOD)) {
         type = "AOD";
         comment += "AOD";
      }   
      if (type!="AOD" && fFriendChainName!="") {
         Error("WriteAnalysisMacro", "Friend chain can be attached only to AOD");
         return;
      }
      if (TObject::TestBit(AliAnalysisGrid::kUseMC)) comment += "/MC";
      else comment += " data";
      out << "const char *anatype = \"" << type.Data() << "\";" << endl << endl;
      func.ReplaceAll(".C", "");
      out << "void " << func.Data() << "()" << endl; 
      out << "{" << endl;
      out << comment.Data() << endl;
      out << "// Automatically generated analysis steering macro executed in grid subjobs" << endl << endl;
      out << "   TStopwatch timer;" << endl;
      out << "   timer.Start();" << endl << endl;
      out << "// load base root libraries" << endl;
      out << "   gSystem->Load(\"libTree\");" << endl;
      out << "   gSystem->Load(\"libGeom\");" << endl;
      out << "   gSystem->Load(\"libVMC\");" << endl;
      out << "   gSystem->Load(\"libPhysics\");" << endl << endl;
      out << "// Load analysis framework libraries" << endl;
      if (!fPackages) {
         out << "   gSystem->Load(\"libSTEERBase\");" << endl;
         out << "   gSystem->Load(\"libESD\");" << endl;
         out << "   gSystem->Load(\"libAOD\");" << endl;
         out << "   gSystem->Load(\"libANALYSIS\");" << endl;
         out << "   gSystem->Load(\"libANALYSISalice\");" << endl;
         out << "   gSystem->Load(\"libCORRFW\");" << endl << endl;
      } else {
         TIter next(fPackages);
         TObject *obj;
         TString pkgname;
         Bool_t hasSTEERBase = kFALSE;
         Bool_t hasESD = kFALSE;
         Bool_t hasAOD = kFALSE;
         Bool_t hasANALYSIS = kFALSE;
         Bool_t hasANALYSISalice = kFALSE;
         Bool_t hasCORRFW = kFALSE;
         while ((obj=next())) {
            pkgname = obj->GetName();
            if (pkgname == "STEERBase" ||
                pkgname == "STEERBase.par") hasSTEERBase = kTRUE;
            if (pkgname == "ESD" ||
                pkgname == "ESD.par")       hasESD = kTRUE;
            if (pkgname == "AOD" ||
                pkgname == "AOD.par")       hasAOD = kTRUE;
            if (pkgname == "ANALYSIS" ||
                pkgname == "ANALYSIS.par")  hasANALYSIS = kTRUE;
            if (pkgname == "ANALYSISalice" ||
                pkgname == "ANALYSISalice.par") hasANALYSISalice = kTRUE;
            if (pkgname == "CORRFW" ||
                pkgname == "CORRFW.par")    hasCORRFW = kTRUE;
         }   
         if (!hasSTEERBase) out << "   gSystem->Load(\"libSTEERBase\");" << endl;
         else out << "   if (!SetupPar(\"STEERBase\")) return;" << endl;
         if (!hasESD)       out << "   gSystem->Load(\"libESD\");" << endl;
         else out << "   if (!SetupPar(\"ESD\")) return;" << endl;
         if (!hasAOD)       out << "   gSystem->Load(\"libAOD\");" << endl;
         else out << "   if (!SetupPar(\"AOD\")) return;" << endl;
         if (!hasANALYSIS)  out << "   gSystem->Load(\"libANALYSIS\");" << endl;
         else out << "   if (!SetupPar(\"ANALYSIS\")) return;" << endl;
         if (!hasANALYSISalice)   out << "   gSystem->Load(\"libANALYSISalice\");" << endl;
         else out << "   if (!SetupPar(\"ANALYSISalice\")) return;" << endl;
         if (!hasCORRFW)    out << "   gSystem->Load(\"libCORRFW\");" << endl << endl;
         else out << "   if (!SetupPar(\"CORRFW\")) return;" << endl << endl;
         out << "// Compile other par packages" << endl;
         next.Reset();
         while ((obj=next())) {
            pkgname = obj->GetName();
            if (pkgname == "STEERBase" ||
                pkgname == "STEERBase.par" ||
                pkgname == "ESD" ||
                pkgname == "ESD.par" ||
                pkgname == "AOD" ||
                pkgname == "AOD.par" ||
                pkgname == "ANALYSIS" ||
                pkgname == "ANALYSIS.par" ||
                pkgname == "ANALYSISalice" ||
                pkgname == "ANALYSISalice.par" ||
                pkgname == "CORRFW" ||
                pkgname == "CORRFW.par") continue;
            out << "   if (!SetupPar(\"" << obj->GetName() << "\")) return;" << endl;
         }   
      }   
      out << "// include path" << endl;
      if (fIncludePath.Length()) out << "   gSystem->AddIncludePath(\"" << fIncludePath.Data() << "\");" << endl;
      out << "   gSystem->AddIncludePath(\"-I$ALICE_ROOT/include\");" << endl << endl;
      if (fAdditionalLibs.Length()) {
         out << "// Add aditional AliRoot libraries" << endl;
         TObjArray *list = fAdditionalLibs.Tokenize(" ");
         TIter next(list);
         TObjString *str;
         while((str=(TObjString*)next())) {
            if (str->GetString().Contains(".so"))
               out << "   gSystem->Load(\"" << str->GetString().Data() << "\");" << endl;
         }
         if (list) delete list;
      }
      out << endl;
      out << "// analysis source to be compiled at runtime (if any)" << endl;
      if (fAnalysisSource.Length()) {
         TObjArray *list = fAnalysisSource.Tokenize(" ");
         TIter next(list);
         TObjString *str;
         while((str=(TObjString*)next())) {
            out << "   gROOT->ProcessLine(\".L " << str->GetString().Data() << "+g\");" << endl;
         }   
         if (list) delete list;
      }
      out << endl;
      out << "// connect to AliEn and make the chain" << endl;
      out << "   if (!TGrid::Connect(\"alien://\")) return;" << endl;
      if (IsUsingTags()) {
         out << "   TChain *chain = CreateChainFromTags(\"wn.xml\", anatype);" << endl << endl;
      } else {
         if(fFriendChainName!="AliAOD.VertexingHF.root") {
            out << "   TChain *chain = CreateChain(\"wn.xml\", anatype);" << endl << endl;    
         } else {
            out << "   // Check if the macro to create the chain was provided" << endl;
            out << "   if (gSystem->AccessPathName(\"MakeAODInputChain.C\")) {" << endl;
            out << "      ::Error(\"" << func.Data() << "\", \"File MakeAODInputChain.C not provided. Aborting.\");" << endl;
            out << "      return;" << endl;
            out << "   }" << endl;
            out << "   gROOT->LoadMacro(\"MakeAODInputChain.C\");" << endl;
            out << "   TChain *chain = MakeAODInputChain(\"wn.xml\",\"none\");" << endl << endl;
         }  
      }   
      out << "// read the analysis manager from file" << endl;
      TString analysisFile = fExecutable;
      analysisFile.ReplaceAll(".sh", ".root");
      out << "   TFile *file = TFile::Open(\"" << analysisFile << "\");" << endl;
      out << "   if (!file) return;" << endl; 
      out << "   TIter nextkey(file->GetListOfKeys());" << endl;
      out << "   AliAnalysisManager *mgr = 0;" << endl;
      out << "   TKey *key;" << endl;
      out << "   while ((key=(TKey*)nextkey())) {" << endl;
      out << "      if (!strcmp(key->GetClassName(), \"AliAnalysisManager\"))" << endl;
      out << "         mgr = (AliAnalysisManager*)file->Get(key->GetName());" << endl;
      out << "   };" << endl;
      out << "   if (!mgr) {" << endl;
      out << "      ::Error(\"" << func.Data() << "\", \"No analysis manager found in file" << analysisFile <<"\");" << endl;
      out << "      return;" << endl;
      out << "   }" << endl << endl;
      out << "   mgr->PrintStatus();" << endl;
      out << "   mgr->StartAnalysis(\"localfile\", chain);" << endl;
      out << "   timer.Stop();" << endl;
      out << "   timer.Print();" << endl;
      out << "}" << endl << endl;
      if (IsUsingTags()) {
         out << "TChain* CreateChainFromTags(const char *xmlfile, const char *type=\"ESD\")" << endl;
         out << "{" << endl;
         out << "// Create a chain using tags from the xml file." << endl;
         out << "   TAlienCollection* coll = TAlienCollection::Open(xmlfile);" << endl;
         out << "   if (!coll) {" << endl;
         out << "      ::Error(\"CreateChainFromTags\", \"Cannot create an AliEn collection from %s\", xmlfile);" << endl;
         out << "      return NULL;" << endl;
         out << "   }" << endl;
         out << "   TGridResult* tagResult = coll->GetGridResult(\"\",kFALSE,kFALSE);" << endl;
         out << "   AliTagAnalysis *tagAna = new AliTagAnalysis(type);" << endl;
         out << "   tagAna->ChainGridTags(tagResult);" << endl << endl;
         out << "   AliRunTagCuts      *runCuts = new AliRunTagCuts();" << endl;
         out << "   AliLHCTagCuts      *lhcCuts = new AliLHCTagCuts();" << endl;
         out << "   AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();" << endl;
         out << "   AliEventTagCuts    *evCuts  = new AliEventTagCuts();" << endl;
         out << "   // Check if the cuts configuration file was provided" << endl;
         out << "   if (!gSystem->AccessPathName(\"ConfigureCuts.C\")) {" << endl;
         out << "      gROOT->LoadMacro(\"ConfigureCuts.C\");" << endl;
         out << "      ConfigureCuts(runCuts, lhcCuts, detCuts, evCuts);" << endl;
         out << "   }" << endl;
         if (fFriendChainName=="") {
            out << "   TChain *chain = tagAna->QueryTags(runCuts, lhcCuts, detCuts, evCuts);" << endl;
         } else {
            out << "   TString tmpColl=\"tmpCollection.xml\";" << endl;
            out << "   tagAna->CreateXMLCollection(tmpColl.Data(),runCuts, lhcCuts, detCuts, evCuts);" << endl;
            out << "   TChain *chain = CreateChain(tmpColl.Data(),type);" << endl;
         }
         out << "   if (!chain || !chain->GetNtrees()) return NULL;" << endl;
         out << "   chain->ls();" << endl;
         out << "   return chain;" << endl;
         out << "}" << endl << endl;
         if (gSystem->AccessPathName("ConfigureCuts.C")) {
            TString msg = "\n#####   You may want to provide a macro ConfigureCuts.C with a method:\n";
            msg += "   void ConfigureCuts(AliRunTagCuts *runCuts,\n";
            msg += "                      AliLHCTagCuts *lhcCuts,\n";
            msg += "                      AliDetectorTagCuts *detCuts,\n";
            msg += "                      AliEventTagCuts *evCuts)";
            Info("WriteAnalysisMacro", msg.Data());
         }
      } 
      if (!IsUsingTags() || fFriendChainName!="") {
         out <<"//________________________________________________________________________________" << endl;
         out << "TChain* CreateChain(const char *xmlfile, const char *type=\"ESD\")" << endl;
         out << "{" << endl;
         out << "// Create a chain using url's from xml file" << endl;
         out << "   TString treename = type;" << endl;
         out << "   treename.ToLower();" << endl;
         out << "   treename += \"Tree\";" << endl;
         out << "   printf(\"***************************************\\n\");" << endl;
         out << "   printf(\"    Getting chain of trees %s\\n\", treename.Data());" << endl;
         out << "   printf(\"***************************************\\n\");" << endl;
         out << "   TAlienCollection *coll = TAlienCollection::Open(xmlfile);" << endl;
         out << "   if (!coll) {" << endl;
         out << "      ::Error(\"CreateChain\", \"Cannot create an AliEn collection from %s\", xmlfile);" << endl;
         out << "      return NULL;" << endl;
         out << "   }" << endl;
         out << "   TChain *chain = new TChain(treename);" << endl;
         if(fFriendChainName!="") {
            out << "   TChain *chainFriend = new TChain(treename);" << endl;
         }
         out << "   coll->Reset();" << endl;
         out << "   while (coll->Next()) {" << endl;
         out << "      chain->Add(coll->GetTURL(\"\"));" << endl;
         if(fFriendChainName!="") {
            out << "      TString fileFriend=coll->GetTURL(\"\");" << endl;
            out << "      fileFriend.ReplaceAll(\"AliAOD.root\",\""<<fFriendChainName.Data()<<"\");" << endl;
            out << "      fileFriend.ReplaceAll(\"AliAODs.root\",\""<<fFriendChainName.Data()<<"\");" << endl;
            out << "      chainFriend->Add(fileFriend.Data());" << endl;
         }
         out << "   }" << endl;
         out << "   if (!chain->GetNtrees()) {" << endl;
         out << "      ::Error(\"CreateChain\", \"No tree found from collection %s\", xmlfile);" << endl;
         out << "      return NULL;" << endl;
         out << "   }" << endl;
         if(fFriendChainName!="") {
            out << "   chain->AddFriend(chainFriend);" << endl;
         }
         out << "   return chain;" << endl;
         out << "}" << endl << endl;
      }   
      if (fPackages) {
         out <<"//________________________________________________________________________________" << endl;
         out << "Bool_t SetupPar(const char *package) {" << endl;
         out << "// Compile the package and set it up." << endl;
         out << "   TString pkgdir = package;" << endl;
         out << "   pkgdir.ReplaceAll(\".par\",\"\");" << endl;
         out << "   gSystem->Exec(Form(\"tar xvzf %s.par\", pkgdir.Data()));" << endl;
         out << "   TString cdir = gSystem->WorkingDirectory();" << endl;
         out << "   gSystem->ChangeDirectory(pkgdir);" << endl;
         out << "   // Check for BUILD.sh and execute" << endl;
         out << "   if (!gSystem->AccessPathName(\"PROOF-INF/BUILD.sh\")) {" << endl;
         out << "      printf(\"*******************************\\n\");" << endl;
         out << "      printf(\"*** Building PAR archive    ***\\n\");" << endl;
         out << "      printf(\"*******************************\\n\");" << endl;
         out << "      if (gSystem->Exec(\"PROOF-INF/BUILD.sh\")) {" << endl;
         out << "         ::Error(\"SetupPar\", \"Cannot build par archive %s\", pkgdir.Data());" << endl;
         out << "         gSystem->ChangeDirectory(cdir);" << endl;
         out << "         return kFALSE;" << endl;
         out << "      }" << endl;
         out << "   } else {" << endl;
         out << "      ::Error(\"SetupPar\",\"Cannot access PROOF-INF/BUILD.sh for package %s\", pkgdir.Data());" << endl;
         out << "      gSystem->ChangeDirectory(cdir);" << endl;
         out << "      return kFALSE;" << endl;
         out << "   }" << endl;
         out << "   // Check for SETUP.C and execute" << endl;
         out << "   if (!gSystem->AccessPathName(\"PROOF-INF/SETUP.C\")) {" << endl;
         out << "      printf(\"*******************************\\n\");" << endl;
         out << "      printf(\"***    Setup PAR archive    ***\\n\");" << endl;
         out << "      printf(\"*******************************\\n\");" << endl;
         out << "      gROOT->Macro(\"PROOF-INF/SETUP.C\");" << endl;
         out << "   } else {" << endl;
         out << "      ::Error(\"SetupPar\",\"Cannot access PROOF-INF/SETUP.C for package %s\", pkgdir.Data());" << endl;
         out << "      gSystem->ChangeDirectory(cdir);" << endl;
         out << "      return kFALSE;" << endl;
         out << "   }" << endl;
         out << "   // Restore original workdir" << endl;
         out << "   gSystem->ChangeDirectory(cdir);" << endl;
         out << "   return kTRUE;" << endl;
         out << "}" << endl;
      }
      Info("WriteAnalysisMacro", "\n#####   Analysis macro to run on worker nodes <%s> written",fAnalysisMacro.Data());
   }   
   Bool_t copy = kTRUE;
   if (TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      workdir += fGridWorkingDir;
      if (FileExists(fAnalysisMacro)) gGrid->Rm(fAnalysisMacro);
      if (IsUsingTags() && !gSystem->AccessPathName("ConfigureCuts.C")) {
         if (FileExists("ConfigureCuts.C")) gGrid->Rm("ConfigureCuts.C");
         Info("WriteAnalysisMacro", "\n#####   Copying cuts configuration macro: <ConfigureCuts.C> to your alien workspace");
         TFile::Cp("file:ConfigureCuts.C", Form("alien://%s/ConfigureCuts.C", workdir.Data()));
      }   
      Info("WriteAnalysisMacro", "\n#####   Copying analysis macro: <%s> to your alien workspace", fAnalysisMacro.Data());
      TFile::Cp(Form("file:%s",fAnalysisMacro.Data()), Form("alien://%s/%s", workdir.Data(), fAnalysisMacro.Data()));
   }
}

//______________________________________________________________________________
void AliAnalysisAlien::WriteExecutable()
{
// Generate the alien executable script.
   if (!TestBit(AliAnalysisGrid::kSubmit)) {  
      ofstream out;
      out.open(fExecutable.Data(), ios::out);
      if (out.bad()) {
         Error("WriteExecutable", "Bad file name for executable: %s", fExecutable.Data());
         return;
      }
      out << "#!/bin/bash" << endl;
      out << "echo \"=========================================\"" << endl; 
      out << "echo \"############## PATH : ##############\"" << endl;
      out << "echo $PATH" << endl;
      out << "echo \"############## LD_LIBRARY_PATH : ##############\"" << endl;
      out << "echo $LD_LIBRARY_PATH" << endl;
      out << "echo \"############## ROOTSYS : ##############\"" << endl;
      out << "echo $ROOTSYS" << endl;
      out << "echo \"############## which root : ##############\"" << endl;
      out << "which root" << endl;
      out << "echo \"############## ALICE_ROOT : ##############\"" << endl;
      out << "echo $ALICE_ROOT" << endl;
      out << "echo \"############## which aliroot : ##############\"" << endl;
      out << "which aliroot" << endl;
      out << "echo \"=========================================\"" << endl << endl;
//      if (TestBit(AliAnalysisGrid::kTest)) out << "root ";
      out << fExecutableCommand << " "; 
      out << fAnalysisMacro.Data() << endl << endl;
      out << "echo \"======== " << fAnalysisMacro.Data() << " finished ========\"" << endl;
   }   
   Bool_t copy = kTRUE;
   if (TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      TString bindir = Form("%s/bin", workdir.Data());
      if (!DirectoryExists(bindir)) gGrid->Mkdir(bindir);
      workdir += fGridWorkingDir;
      TString executable = Form("%s/bin/%s", gGrid->GetHomeDirectory(), fExecutable.Data());
      if (FileExists(executable)) gGrid->Rm(executable);
      Info("CreateJDL", "\n#####   Copying executable file <%s> to your AliEn bin directory", fExecutable.Data());
      TFile::Cp(Form("file:%s",fExecutable.Data()), Form("alien://%s", executable.Data()));
   } 
}

//______________________________________________________________________________
void AliAnalysisAlien::WriteProductionFile(const char *filename) const
{
// Write the production file to be submitted by LPM manager. The format is:
// First line: full_path_to_jdl estimated_no_subjobs_per_master
// Next lines: full_path_to_dataset XXX (XXX is a string)
// To submit, one has to: submit jdl XXX for all lines
   ofstream out;
   out.open(filename, ios::out);
   if (out.bad()) {
      Error("WriteProductionFile", "Bad file name: %s", filename);
      return;
   }
   TString workdir = gGrid->GetHomeDirectory();
   workdir += fGridWorkingDir;
   Int_t njobspermaster = 1000*fNrunsPerMaster/fSplitMaxInputFileNumber;
   TString locjdl = Form("%s/%s", workdir.Data(),fJDLName.Data());
   out << locjdl << " " << njobspermaster << endl;
   Int_t nmasterjobs = fInputFiles->GetEntries();
   for (Int_t i=0; i<nmasterjobs; i++) {
      out << Form("%s", fInputFiles->At(i)->GetName()) << " " << Form("%03d", i) << endl;
   }
   Info("WriteProductionFile", "\n#####   Copying production file <%s> to your work directory", filename);
   TFile::Cp(Form("file:%s",filename), Form("alien://%s/%s", workdir.Data(),filename));   
}

//______________________________________________________________________________
void AliAnalysisAlien::WriteValidationScript()
{
// Generate the alien validation script.
   // Generate the validation script
   TObjString *os;
   TString validationScript = fExecutable;
   validationScript.ReplaceAll(".sh", "_validation.sh");
   if (!Connect()) {
      Error("WriteValidationScript", "Alien connection required");
      return;
   }
   TString out_stream = "";
   if (!TestBit(AliAnalysisGrid::kTest)) out_stream = " >> stdout";
   if (!TestBit(AliAnalysisGrid::kSubmit)) {  
      ofstream out;
      out.open(validationScript, ios::out);
      out << "#!/bin/bash" << endl;
      out << "##################################################" << endl;
      out << "validateout=`dirname $0`" << endl;
      out << "validatetime=`date`" << endl;
      out << "validated=\"0\";" << endl;
      out << "error=0" << endl;
      out << "if [ -z $validateout ]" << endl;
      out << "then" << endl;
      out << "    validateout=\".\"" << endl;
      out << "fi" << endl << endl;
      out << "cd $validateout;" << endl;
      out << "validateworkdir=`pwd`;" << endl << endl;
      out << "echo \"*******************************************************\"" << out_stream << endl;
      out << "echo \"* Automatically generated validation script           *\""  << out_stream << endl;
      out << "" << endl;
      out << "echo \"* Time:    $validatetime \""  << out_stream << endl;
      out << "echo \"* Dir:     $validateout\""  << out_stream << endl;
      out << "echo \"* Workdir: $validateworkdir\""  << out_stream << endl;
      out << "echo \"* ----------------------------------------------------*\""  << out_stream << endl;
      out << "ls -la ./"  << out_stream << endl;
      out << "echo \"* ----------------------------------------------------*\""  << out_stream << endl << endl;
      out << "##################################################" << endl;

      out << "" << endl;
      out << "parArch=`grep -Ei \"Cannot Build the PAR Archive\" stderr`" << endl;
      out << "segViol=`grep -Ei \"Segmentation violation\" stderr`" << endl;
      out << "segFault=`grep -Ei \"Segmentation fault\" stderr`" << endl;
      out << "" << endl;

      out << "if [ ! -f stderr ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"* ########## Job not validated - no stderr  ###\" " << out_stream << endl;
      out << "   echo \"Error = $error\" " << out_stream << endl;
      out << "fi" << endl;

      out << "if [ \"$parArch\" != \"\" ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"* ########## Job not validated - PAR archive not built  ###\" " << out_stream << endl;
      out << "   echo \"$parArch\" " << out_stream << endl;
      out << "   echo \"Error = $error\" " << out_stream << endl;
      out << "fi" << endl;

      out << "if [ \"$segViol\" != \"\" ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"* ########## Job not validated - Segment. violation  ###\" " << out_stream << endl;
      out << "   echo \"$segViol\" " << out_stream << endl;
      out << "   echo \"Error = $error\" " << out_stream << endl;
      out << "fi" << endl;

      out << "if [ \"$segFault\" != \"\" ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"* ########## Job not validated - Segment. fault  ###\" " << out_stream << endl;
      out << "   echo \"$segFault\" " << out_stream << endl;
      out << "   echo \"Error = $error\" " << out_stream << endl;
      out << "fi" << endl;

      // Part dedicated to the specific analyses running into the train

      TObjArray *arr = fOutputFiles.Tokenize(" ");
      TIter next1(arr);
      TString output_file;
      while ((os=(TObjString*)next1())) { 
         output_file = os->GetString();
         Int_t index = output_file.Index("@");
         if (index > 0) output_file.Remove(index);
         out << "if ! [ -f " << output_file.Data() << " ] ; then" << endl;
         out << "   error=1" << endl;
         out << "   echo \"Output file(s) not found. Job FAILED !\""  << out_stream << endl;
         out << "   echo \"Output file(s) not found. Job FAILED !\" >> stderr" << endl;
         out << "fi" << endl;
      }   
      delete arr;
      out << "if ! [ -f outputs_valid ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"Output files were not validated by the analysis manager\" >> stdout" << endl;
      out << "   echo \"Output files were not validated by the analysis manager\" >> stderr" << endl;
      out << "fi" << endl;
      
      out << "if [ $error = 0 ] ; then" << endl;
      out << "   echo \"* ----------------   Job Validated  ------------------*\""  << out_stream << endl;
      out << "fi" << endl;

      out << "echo \"* ----------------------------------------------------*\""  << out_stream << endl;
      out << "echo \"*******************************************************\""  << out_stream << endl;
      out << "cd -" << endl;
      out << "exit $error" << endl;
   }    
   Bool_t copy = kTRUE;
   if (TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      workdir += fGridWorkingDir;
      Info("CreateJDL", "\n#####   Copying validation script <%s> to your AliEn working space", validationScript.Data());
      if (FileExists(validationScript)) gGrid->Rm(validationScript);
      TFile::Cp(Form("file:%s",validationScript.Data()), Form("alien://%s/%s", workdir.Data(),validationScript.Data()));
   } 
}
