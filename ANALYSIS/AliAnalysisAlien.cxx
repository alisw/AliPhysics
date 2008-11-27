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
#include "TFileMerger.h"
#include "AliAnalysisManager.h"
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
                  fRunNumbers(),
                  fExecutable(),
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
                  fRunNumbers(),
                  fExecutable(),
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
                  fRunNumbers(other.fRunNumbers),
                  fExecutable(other.fExecutable),
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
                  fInputFiles(0),
                  fPackages(0)
{
// Copy ctor.
   fGridJDL = (TGridJDL*)gROOT->ProcessLine("new TAlienJDL()");
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
      fRunNumbers              = other.fRunNumbers;
      fExecutable              = other.fExecutable;
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
void AliAnalysisAlien::AddRunNumber(Int_t run)
{
// Add a run number to the list of runs to be processed.
   if (fRunNumbers.Length()) fRunNumbers += " ";
   fRunNumbers += Form("%d", run);
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
   if (!gGrid->Cd(workdir)) {
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
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::CheckInputData()
{
// Check validity of input data. If necessary, create xml files.
   if (!fInputFiles && !fRunNumbers.Length()) {
       Error("CheckInputData", "You have to specify either a set of run numbers or some existing grid files. Use AddRunNumber()/AddDataFile().");
      return kFALSE;
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
   if (!fRunNumbers.Length()) return kTRUE;
   // Check validity of alien data directory
   if (!fGridDataDir.Length()) {
      Error("CkeckInputData", "AliEn path to base data directory must be set.\n = Use: SetGridDataDir()");
      return kFALSE;
   }
   if (!gGrid->Cd(fGridDataDir)) {
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
      arr = fRunNumbers.Tokenize(" ");
      TIter next(arr);
      while ((os=(TObjString*)next())) {
         path = Form("%s/%s ", fGridDataDir.Data(), os->GetString().Data());
         if (!gGrid->Cd(path)) {
            Error("CheckInputData", "Run number %s not found in path: %s", os->GetString().Data(), path.Data());
            return kFALSE;
         }
         path = Form("%s/%s.xml", workdir.Data(),os->GetString().Data());
         TString msg = "\n#####   file: ";
         msg += path;
         msg += " type: xml_collection;";
         if (use_tags) msg += " using_tags: Yes";
         else          msg += " using_tags: No";
         Info("CheckDataType", msg.Data());
         AddDataFile(path);
      }
      delete arr;   
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
   if (TestBit(AliAnalysisGrid::kTest)) options += "-l 10 ";
   TString conditions = "";
   
   TString file;
   TString path;
   if (!fRunNumbers.Length()) return kTRUE;   
   // Several runs
   TObjArray *arr = fRunNumbers.Tokenize(" ");
   TObjString *os;
   TIter next(arr);
   while ((os=(TObjString*)next())) {
      path = Form("%s/%s ", fGridDataDir.Data(), os->GetString().Data());
      if (TestBit(AliAnalysisGrid::kTest)) file = "wn.xml";
      else file = Form("%s.xml", os->GetString().Data());
      if (FileExists(file) && !TestBit(AliAnalysisGrid::kTest)) {
         Info("CreateDataset", "\n#####   Removing previous dataset %s", file.Data());
         gGrid->Rm(file); 
      }
      command = "find ";
      command += options;
      command += path;
      command += pattern;
//      conditions = Form(" > %s", file.Data());
      command += conditions;
      TGridResult *res = gGrid->Command(command);
      if (res) delete res;
      // Write standard output to file
      gROOT->ProcessLine(Form("gGrid->Stdout(); > %s", file.Data()));
      if (TestBit(AliAnalysisGrid::kTest)) break;
      // Copy xml file to alien space
      TFile::Cp(Form("file:%s",file.Data()), Form("alien://%s/%s",workdir.Data(), file.Data()));
      if (!FileExists(file)) {
         Error("CreateDataset", "Command %s did NOT succeed", command.Data());
         delete arr;
         return kFALSE;
      }
   }   
   delete arr;
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
         if (!gGrid->Cd(fGridOutputDir)) {
            if (gGrid->Mkdir(fGridOutputDir)) {
               Info("CreateJDL", "\n#####   Created alien output directory %s", fGridOutputDir.Data());
            } else {
               Error("CreateJDL", "Could not create alien output directory %s", fGridOutputDir.Data());
               error = kTRUE;
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
      if (fInputFiles) {
         TIter next(fInputFiles);
         while ((os=(TObjString*)next()))
            fGridJDL->AddToInputDataCollection(Form("LF:%s,nodownload", os->GetString().Data()));
      }      
      fGridJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(), fAnalysisMacro.Data()));
      fGridJDL->AddToInputSandbox(Form("LF:%s/analysis.root", workdir.Data()));
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
            fGridJDL->AddToOutputArchive(os->GetString().Data());
         delete arr;
      }      
      fGridJDL->SetOutputDirectory(Form("%s/%s/#alien_counter_03i#", workdir.Data(), fGridOutputDir.Data()));
      arr = fOutputFiles.Tokenize(" ");
      TIter next(arr);
      while ((os=(TObjString*)next())) fGridJDL->AddToOutputSandbox(os->GetString());
      delete arr;
//      fGridJDL->SetPrice((UInt_t)fPrice);
      fGridJDL->SetValue("Price", Form("\"%d\"", fPrice));
      fGridJDL->SetValidationCommand(Form("%s/validate.sh", workdir.Data()));
      if (fMasterResubmitThreshold) fGridJDL->SetValue("MasterResubmitThreshold", Form("\"%d%%\"", fMasterResubmitThreshold));
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
      index = sjdl.Index("User");
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
      sjdl.Prepend("JobTag = \"Automatically generated analysis JDL\";\n");
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
   }
   // Copy jdl to grid workspace   
   if (!copy) {
      Info("CreateJDL", "\n#####   You may want to review jdl:%s and analysis macro:%s before running in <submit> mode", fJDLName.Data(), fAnalysisMacro.Data());
   } else {
      Info("CreateJDL", "\n#####   Copying JDL file <%s> to your AliEn working space", fJDLName.Data());
      if (FileExists(fJDLName)) gGrid->Rm(fJDLName);
      TFile::Cp(Form("file:%s",fJDLName.Data()), Form("alien://%s/%s", workdir.Data(), fJDLName.Data()));
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
Bool_t AliAnalysisAlien::FileExists(const char *lfn) const
{
// Returns true if file exists.
   if (!gGrid) {
      Error("FileExists", "No connection to grid");
      return kFALSE;
   }
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
   fRunNumbers                 = "";
   fExecutable                 = "analysis.sh";
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
   fGridOutputDir              = "output";
   fOutputArchive              = "log_archive.zip:stdout,stderr root_archive.zip:*.root";
   fOutputFiles                = "";  // Like "AliAODs.root histos.root"
   fInputFormat                = "xml-single";
   fJDLName                    = "analysis.jdl";
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
   TString output = Form("/%s/%s/%s", gGrid->GetHomeDirectory(), fGridWorkingDir.Data(), fGridOutputDir.Data());
   if (!gGrid->Cd(output)) output = Form("/%s/%s", gGrid->GetHomeDirectory(), fGridOutputDir.Data());
   if (!gGrid->Cd(output)) {
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
   Bool_t merged = kTRUE;
   while((str=(TObjString*)next())) {
      output_file = str->GetString();
      Int_t index = output_file.Index("@");
      if (index > 0) output_file.Remove(index);
      command = Form("find %s/ *%s", output.Data(), output_file.Data());
      printf("command: %s\n", command.Data());
      TGridResult *res = gGrid->Command(command);
      if (!res) continue;
      TFileMerger *fm = 0;
      TIter nextmap(res);
      TMap *map;
      while ((map=(TMap*)nextmap())) {
         TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
         if (!objs || !objs->GetString().Length()) {
            delete res;
            continue;
         }   
         if (!fm) {
            fm = new TFileMerger(kFALSE);
            fm->SetFastMethod(kTRUE);
            fm->OutputFile(output_file);
         }
         fm->AddFile(objs->GetString());   
      }
      if (!fm || !fm->GetMergeList() || !fm->GetMergeList()->GetSize()) {
         Warning("MergeOutputs", "No <%s> files found.", output_file.Data());
         merged = kFALSE;
         delete res;
         continue;
      }
      if (!fm->Merge()) {
         Error("MergeOutputs", "Could not merge all <%s> files", output_file.Data());
         merged = kFALSE;
      } else {
         Info("MergeOutputs", "\n#####   Merged %d output files <%s>", fm->GetMergeList()->GetSize(), output_file.Data());
      }   
      delete fm;
      delete res;
   } 
   if (!merged) {
      Error("MergeOutputs", "Terminate() will  NOT be executed");
   }  
   return merged;
}   

//______________________________________________________________________________
void AliAnalysisAlien::StartAnalysis(Long64_t /*nentries*/, Long64_t /*firstEntry*/)
{
// Start remote grid analysis.
   
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
      return;
   } else {
      Info("StartAnalysis","\n##### FULL ANALYSIS MODE ##### Producing needed files and submitting your analysis job...");   
   }   
      
   if (!Connect()) {
      Error("StartAnalysis", "Cannot start grid analysis without grid connection");
      return;
   }   
   if (!CheckInputData()) {
      Error("StartAnalysis", "There was an error in preprocessing your requested input data");
      return;
   }   
   CreateDataset(fDataPattern);
   WriteAnalysisFile();   
   WriteAnalysisMacro();
   WriteExecutable();
   WriteValidationScript();
   if (!CreateJDL()) return;
   if (TestBit(AliAnalysisGrid::kOffline)) return;
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
         gSystem->Exec(Form("rm %s", output_file.Data()));
      }
      delete list;
      gSystem->Exec(Form("bash %s 2>stderr", fExecutable.Data()));
      gSystem->Exec("bash validate.sh");
//      gSystem->Exec("cat stdout");
      return;
   }
   // Submit AliEn job
   CdWork();
   TGridResult *res = gGrid->Command(Form("submit %s", fJDLName.Data()));
   TString jobID = "";
   if (res) {
      const char *cjobId = res->GetKey(0,"jobId");
      if (!cjobId) {
         Error("StartAnalysis", "Your JDL %s could not be submitted", fJDLName.Data());
         return;
      } else {
         Info("StartAnalysis", "\n_______________________________________________________________________ \
         \n#####   Your JDL %s was successfully submitted. \nTHE JOB ID IS: %s \
      \n_______________________________________________________________________",
                fJDLName.Data(), cjobId);
         jobID = cjobId;      
      }          
      delete res;
   }   
   Info("StartAnalysis", "\n#### STARTING AN ALIEN SHELL FOR YOU. EXIT WHEN YOUR JOB %s HAS FINISHED. #### \
   \n You may exit at any time and terminate the job later using the option <terminate> \
   \n ##################################################################################", jobID.Data());
   gGrid->Shell();
}

//______________________________________________________________________________
void AliAnalysisAlien::WriteAnalysisFile()
{
// Write current analysis manager into the file analysis.root
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
      TFile *file = TFile::Open("analysis.root", "RECREATE");
      if (file) {
         mgr->Write();
         delete file;
      }
      if (cdir) cdir->cd();
      Info("WriteAnalysisFile", "\n#####   Analysis manager: %s wrote to file <analysis.root>\n", mgr->GetName());
   }   
   Bool_t copy = kTRUE;
   if (TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      workdir += fGridWorkingDir;
      Info("CreateJDL", "\n#####   Copying file <analysis.root> containing your initialized analysis manager to your alien workspace");
      if (FileExists("analysis.root")) gGrid->Rm("analysis.root");
      TFile::Cp("file:analysis.root", Form("alien://%s/analysis.root", workdir.Data()));
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
      if (TObject::TestBit(AliAnalysisGrid::kUseMC)) comment += "/MC";
      else comment += " data";
      out << "const char *anatype = \"" << type.Data() << "\";" << endl << endl;
      func.ReplaceAll(".C", "");
      out << "void " << func.Data() << "()" << endl; 
      out << "{" << endl;
      out << comment.Data() << endl;
      out << "// Automatically generated analysis steering macro executed in grid subjobs" << endl << endl;
      out << "// load base root libraries" << endl;
      out << "   gSystem->Load(\"libTree\");" << endl;
      out << "   gSystem->Load(\"libGeom\");" << endl;
      out << "   gSystem->Load(\"libVMC\");" << endl;
      out << "   gSystem->Load(\"libPhysics\");" << endl << endl;
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
      if (!fPackages) {
         out << "// Load analysis framework libraries" << endl;
         out << "   gSystem->Load(\"libSTEERBase\");" << endl;
         out << "   gSystem->Load(\"libESD\");" << endl;
         out << "   gSystem->Load(\"libAOD\");" << endl;
         out << "   gSystem->Load(\"libANALYSIS\");" << endl;
         out << "   gSystem->Load(\"libANALYSISalice\");" << endl << endl;
         out << "// include path (remove if using par files)" << endl;
         out << "   gROOT->ProcessLine(\".include $ALICE_ROOT/include\");" << endl << endl;
      } else {
         out << "// Compile all par packages" << endl;
         TIter next(fPackages);
         TObject *obj;
         while ((obj=next())) 
            out << "   if (!SetupPar(\"" << obj->GetName() << "\")) return;" << endl;
      }   
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
         out << "   TChain *chain = CreateChain(\"wn.xml\", anatype);" << endl << endl;      
      }   
      out << "// read the analysis manager from file" << endl;
      out << "   TFile *file = TFile::Open(\"analysis.root\");" << endl;
      out << "   if (!file) return;" << endl; 
      out << "   TIter nextkey(file->GetListOfKeys());" << endl;
      out << "   AliAnalysisManager *mgr = 0;" << endl;
      out << "   TKey *key;" << endl;
      out << "   while ((key=(TKey*)nextkey())) {" << endl;
      out << "      if (!strcmp(key->GetClassName(), \"AliAnalysisManager\"))" << endl;
      out << "         mgr = (AliAnalysisManager*)file->Get(key->GetName());" << endl;
      out << "   };" << endl;
      out << "   if (!mgr) {" << endl;
      out << "      ::Error(\"" << func.Data() << "\", \"No analysis manager found in file analysis.root\");" << endl;
      out << "      return;" << endl;
      out << "   }" << endl << endl;
      out << "   mgr->PrintStatus();" << endl;
      out << "   mgr->StartAnalysis(\"localfile\", chain);" << endl;
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
         out << "   TChain *chain = tagAna->QueryTags(runCuts, lhcCuts, detCuts, evCuts);" << endl;
         out << "   if (!chain || !chain->GetNtrees()) return NULL;" << endl;
         out << "   chain->ls();" << endl;
         out << "   return chain;" << endl;
         out << "}" << endl;
         if (gSystem->AccessPathName("ConfigureCuts.C")) {
            TString msg = "\n#####   You may want to provide a macro ConfigureCuts.C with a method:\n";
            msg += "   void ConfigureCuts(AliRunTagCuts *runCuts,\n";
            msg += "                      AliLHCTagCuts *lhcCuts,\n";
            msg += "                      AliDetectorTagCuts *detCuts,\n";
            msg += "                      AliEventTagCuts *evCuts)";
            Info("WriteAnalysisMacro", msg.Data());
         }
      } else {
         out << "TChain* CreateChain(const char *xmlfile, const char *type=\"ESD\")" << endl;
         out << "{" << endl;
         out << "// Create a chain using url's from xml file" << endl;
         out << "   TString treename = type;" << endl;
         out << "   treename.ToLower();" << endl;
         out << "   treename += \"Tree\";" << endl;
         out << "   printf(\"***************************************\");" << endl;
         out << "   printf(\"    Getting chain of trees %s\\n\", treename);" << endl;
         out << "   printf(\"***************************************\");" << endl;
         out << "   TAlienCollection *coll = TAlienCollection::Open(xmlfile);" << endl;
         out << "   if (!coll) {" << endl;
         out << "      ::Error(\"CreateChain\", \"Cannot create an AliEn collection from %s\", xmlfile);" << endl;
         out << "      return NULL;" << endl;
         out << "   }" << endl;
         out << "   TChain *chain = new TChain(treename);" << endl;
         out << "   coll->Reset();" << endl;
         out << "   while (coll->Next()) chain->Add(coll->GetTURL(\"\"));" << endl;
         out << "   if (!chain->GetNtrees()) {" << endl;
         out << "      ::Error(\"CreateChain\", \"No tree found from collection %s\", xmlfile);" << endl;
         out << "      return NULL;" << endl;
         out << "   }" << endl;
         out << "   return chain;" << endl;
         out << "}" << endl;
      }   
      if (fPackages) {
         out << "Bool_t SetupPar(const char *package) {" << endl;
         out << "// Compile the package and set it up." << endl;
         out << "   TString pkgdir = package;" << endl;
         out << "   pkgdir.ReplaceAll(\".par\",\"\");" << endl;
         out << "   gSystem->Exec(Form(\"tar xvzf %s\", package));" << endl;
         out << "   TString cdir = gSystem->WorkingDirectory();" << endl;
         out << "   gSystem->ChangeDirectory(pkgdir);" << endl;
         out << "   // Check for BUILD.sh and execute" << endl;
         out << "   if (!gSystem->AccessPathName(\"PROOF-INF/BUILD.sh\")) {" << endl;
         out << "      printf(\"*******************************\\n\");" << endl;
         out << "      printf(\"*** Building PAR archive    ***\\n\");" << endl;
         out << "      printf(\"*******************************\\n\");" << endl;
         out << "      if (gSystem->Exec(\"PROOF-INF/BUILD.sh\")) {" << endl;
         out << "         ::Error(\"SetupPar\", \"Cannot build par archive %s\", package);" << endl;
         out << "         gSystem->ChangeDirectory(cdir);" << endl;
         out << "         return kFALSE;" << endl;
         out << "      }" << endl;
         out << "   } else {" << endl;
         out << "      ::Error(\"SetupPar\",\"Cannot access PROOF-INF/BUILD.sh for package %s\", package);" << endl;
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
         out << "      ::Error(\"SetupPar\",\"Cannot access PROOF-INF/SETUP.C for package %s\", package);" << endl;
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
         Error("CreateJDL", "Bad file name for executable: %s", fExecutable.Data());
         return;
      }
      out << "#!/bin/bash" << endl;
      out << "export GCLIENT_SERVER_LIST=\"pcapiserv04.cern.ch:10000|pcapiserv05.cern.ch:10000|pcapiserv06.cern.ch:10000|pcapiserv07.cern.ch:10000\"" << endl;
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
      out << "root -b -q "; 
      out << fAnalysisMacro.Data() << endl << endl;
      out << "echo \"======== " << fAnalysisMacro.Data() << " finished ========\"" << endl;
   }   
   Bool_t copy = kTRUE;
   if (TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      workdir += fGridWorkingDir;
      TString executable = Form("%s/bin/%s", gGrid->GetHomeDirectory(), fExecutable.Data());
      if (FileExists(executable)) gGrid->Rm(executable);
      Info("CreateJDL", "\n#####   Copying executable file <%s> to your AliEn bin directory", fExecutable.Data());
      TFile::Cp(Form("file:%s",fExecutable.Data()), Form("alien://%s", executable.Data()));
   } 
}

//______________________________________________________________________________
void AliAnalysisAlien::WriteValidationScript()
{
// Generate the alien validation script.
   // Generate the validation script
   TObjString *os;
   if (!Connect()) {
      Error("WriteValidationScript", "Alien connection required");
      return;
   }
   TString out_stream = "";
   if (!TestBit(AliAnalysisGrid::kTest)) out_stream = " >> stdout";
   if (!TestBit(AliAnalysisGrid::kSubmit)) {  
      ofstream out;
      out.open("validate.sh", ios::out);
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
      Info("CreateJDL", "\n#####   Copying validation script <validate.sh> to your AliEn working space");
      if (FileExists("validate.sh")) gGrid->Rm("validate.sh");
      TFile::Cp("file:validate.sh", Form("alien://%s/validate.sh", workdir.Data()));
   } 
}
