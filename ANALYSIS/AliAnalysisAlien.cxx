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

#include "AliAnalysisAlien.h"

#include "Riostream.h"
#include "TEnv.h"
#include "TKey.h"
#include "TBits.h"
#include "TError.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMacro.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TGridCollection.h"
#include "TGridJDL.h"
#include "TGridJobStatusList.h"
#include "TGridJobStatus.h"
#include "TFileMerger.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCfg.h"
#include "AliVEventHandler.h"
#include "AliAnalysisDataContainer.h"
#include "AliMultiInputEventHandler.h"

ClassImp(AliAnalysisAlien)
#if 0
;
#endif  

namespace {
  Bool_t copyLocal2Alien(const char* where, const char* loc, const char* rem)
  {
    TString sl(Form("file:%s", loc));
    TString sr(Form("alien://%s", rem));
    Bool_t ret = TFile::Cp(sl, sr);
    if (!ret) { 
      Warning(where, "Failed to copy %s to %s", sl.Data(), sr.Data());
    }
    return ret;
  }
}
    
//______________________________________________________________________________
AliAnalysisAlien::AliAnalysisAlien()
                 :AliAnalysisGrid(),
                  fGridJDL(NULL),
                  fMergingJDL(NULL),
                  fPrice(0),
                  fTTL(0),
                  fSplitMaxInputFileNumber(0),
                  fMaxInitFailed(0),
                  fMasterResubmitThreshold(0),
                  fNtestFiles(0),
                  fNrunsPerMaster(0),
                  fMaxMergeFiles(0),
                  fMaxMergeStages(0),
                  fNsubmitted(0),
                  fProductionMode(0),
                  fOutputToRunNo(0),
                  fMergeViaJDL(0),
                  fFastReadOption(0),
                  fOverwriteMode(1),
                  fNreplicas(2),
                  fNproofWorkers(0),
                  fNproofWorkersPerSlave(0),
                  fProofReset(0),
                  fRunNumbers(),
                  fExecutable(),
                  fExecutableCommand(),
                  fArguments(),
                  fExecutableArgs(),
                  fAnalysisMacro(),
                  fAnalysisSource(),
                  fValidationScript(),
                  fAdditionalRootLibs(),
                  fAdditionalLibs(),
                  fSplitMode(),
                  fAPIVersion(),
                  fROOTVersion(),
                  fAliROOTVersion(),
                  fExternalPackages(),
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
                  fTerminateFiles(),
		            fMergeExcludes(),
                  fIncludePath(),
                  fCloseSE(),
                  fFriendChainName(),
                  fJobTag(),
                  fOutputSingle(),
                  fRunPrefix(),
                  fProofCluster(),
                  fProofDataSet(),
                  fFileForTestMode(),
                  fRootVersionForProof(),
                  fAliRootMode(),
                  fMergeDirName(),
                  fInputFiles(0),
                  fPackages(0),
                  fModules(0),
                  fProofParam()
{
// Dummy ctor.
   SetDefaults();
}

//______________________________________________________________________________
AliAnalysisAlien::AliAnalysisAlien(const char *name)
                 :AliAnalysisGrid(name),
                  fGridJDL(NULL),
                  fMergingJDL(NULL),
                  fPrice(0),
                  fTTL(0),
                  fSplitMaxInputFileNumber(0),
                  fMaxInitFailed(0),
                  fMasterResubmitThreshold(0),
                  fNtestFiles(0),
                  fNrunsPerMaster(0),
                  fMaxMergeFiles(0),
                  fMaxMergeStages(0),
                  fNsubmitted(0),
                  fProductionMode(0),
                  fOutputToRunNo(0),
                  fMergeViaJDL(0),
                  fFastReadOption(0),
                  fOverwriteMode(1),
                  fNreplicas(2),
                  fNproofWorkers(0),
                  fNproofWorkersPerSlave(0),
                  fProofReset(0),
                  fRunNumbers(),
                  fExecutable(),
                  fExecutableCommand(),
                  fArguments(),
                  fExecutableArgs(),
                  fAnalysisMacro(),
                  fAnalysisSource(),
                  fValidationScript(),
                  fAdditionalRootLibs(),
                  fAdditionalLibs(),
                  fSplitMode(),
                  fAPIVersion(),
                  fROOTVersion(),
                  fAliROOTVersion(),
                  fExternalPackages(),
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
                  fTerminateFiles(),
                  fMergeExcludes(),
                  fIncludePath(),
                  fCloseSE(),
                  fFriendChainName(),
                  fJobTag(),
                  fOutputSingle(),
                  fRunPrefix(),
                  fProofCluster(),
                  fProofDataSet(),
                  fFileForTestMode(),
                  fRootVersionForProof(),
                  fAliRootMode(),
                  fMergeDirName(),
                  fInputFiles(0),
                  fPackages(0),
                  fModules(0),
                  fProofParam()
{
// Default ctor.
   SetDefaults();
}

//______________________________________________________________________________
AliAnalysisAlien::AliAnalysisAlien(const AliAnalysisAlien& other)
                 :AliAnalysisGrid(other),
                  fGridJDL(NULL),
                  fMergingJDL(NULL),
                  fPrice(other.fPrice),
                  fTTL(other.fTTL),
                  fSplitMaxInputFileNumber(other.fSplitMaxInputFileNumber),
                  fMaxInitFailed(other.fMaxInitFailed),
                  fMasterResubmitThreshold(other.fMasterResubmitThreshold),
                  fNtestFiles(other.fNtestFiles),
                  fNrunsPerMaster(other.fNrunsPerMaster),
                  fMaxMergeFiles(other.fMaxMergeFiles),
                  fMaxMergeStages(other.fMaxMergeStages),
                  fNsubmitted(other.fNsubmitted),
                  fProductionMode(other.fProductionMode),
                  fOutputToRunNo(other.fOutputToRunNo),
                  fMergeViaJDL(other.fMergeViaJDL),
                  fFastReadOption(other.fFastReadOption),
                  fOverwriteMode(other.fOverwriteMode),
                  fNreplicas(other.fNreplicas),
                  fNproofWorkers(other.fNproofWorkers),
                  fNproofWorkersPerSlave(other.fNproofWorkersPerSlave),
                  fProofReset(other.fProofReset),
                  fRunNumbers(other.fRunNumbers),
                  fExecutable(other.fExecutable),
                  fExecutableCommand(other.fExecutableCommand),
                  fArguments(other.fArguments),
                  fExecutableArgs(other.fExecutableArgs),
                  fAnalysisMacro(other.fAnalysisMacro),
                  fAnalysisSource(other.fAnalysisSource),
                  fValidationScript(other.fValidationScript),
                  fAdditionalRootLibs(other.fAdditionalRootLibs),
                  fAdditionalLibs(other.fAdditionalLibs),
                  fSplitMode(other.fSplitMode),
                  fAPIVersion(other.fAPIVersion),
                  fROOTVersion(other.fROOTVersion),
                  fAliROOTVersion(other.fAliROOTVersion),
                  fExternalPackages(other.fExternalPackages),
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
                  fTerminateFiles(other.fTerminateFiles),
                  fMergeExcludes(other.fMergeExcludes),
                  fIncludePath(other.fIncludePath),
                  fCloseSE(other.fCloseSE),
                  fFriendChainName(other.fFriendChainName),
                  fJobTag(other.fJobTag),
                  fOutputSingle(other.fOutputSingle),
                  fRunPrefix(other.fRunPrefix),
                  fProofCluster(other.fProofCluster),
                  fProofDataSet(other.fProofDataSet),
                  fFileForTestMode(other.fFileForTestMode),
                  fRootVersionForProof(other.fRootVersionForProof),
                  fAliRootMode(other.fAliRootMode),
                  fMergeDirName(other.fMergeDirName),
                  fInputFiles(0),
                  fPackages(0),
                  fModules(0),
                  fProofParam()
{
// Copy ctor.
   fGridJDL = (TGridJDL*)gROOT->ProcessLine("new TAlienJDL()");
   fMergingJDL = (TGridJDL*)gROOT->ProcessLine("new TAlienJDL()");
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
   if (other.fModules) {
      fModules = new TObjArray();
      fModules->SetOwner();
      TIter next(other.fModules);
      AliAnalysisTaskCfg *mod, *crt;
      while ((crt=(AliAnalysisTaskCfg*)next())) {
         mod = new AliAnalysisTaskCfg(*crt);
         fModules->Add(mod);
      }
   }   
}

//______________________________________________________________________________
AliAnalysisAlien::~AliAnalysisAlien()
{
// Destructor.
   delete fGridJDL;
   delete fMergingJDL;
   delete fInputFiles;
   delete fPackages;
   delete fModules;
   fProofParam.DeleteAll();
}   

//______________________________________________________________________________
AliAnalysisAlien &AliAnalysisAlien::operator=(const AliAnalysisAlien& other)
{
// Assignment.
   if (this != &other) {
      AliAnalysisGrid::operator=(other);
      fGridJDL = (TGridJDL*)gROOT->ProcessLine("new TAlienJDL()");
      fMergingJDL = (TGridJDL*)gROOT->ProcessLine("new TAlienJDL()");
      fPrice                   = other.fPrice;
      fTTL                     = other.fTTL;
      fSplitMaxInputFileNumber = other.fSplitMaxInputFileNumber;
      fMaxInitFailed           = other.fMaxInitFailed;
      fMasterResubmitThreshold = other.fMasterResubmitThreshold;
      fNtestFiles              = other.fNtestFiles;
      fNrunsPerMaster          = other.fNrunsPerMaster;
      fMaxMergeFiles           = other.fMaxMergeFiles;
      fMaxMergeStages          = other.fMaxMergeStages;
      fNsubmitted              = other.fNsubmitted;
      fProductionMode          = other.fProductionMode;
      fOutputToRunNo           = other.fOutputToRunNo;
      fMergeViaJDL             = other.fMergeViaJDL;
      fFastReadOption          = other.fFastReadOption;
      fOverwriteMode           = other.fOverwriteMode;
      fNreplicas               = other.fNreplicas;
      fNproofWorkers           = other.fNproofWorkers;
      fNproofWorkersPerSlave   = other.fNproofWorkersPerSlave;
      fProofReset              = other.fProofReset;
      fRunNumbers              = other.fRunNumbers;
      fExecutable              = other.fExecutable;
      fExecutableCommand       = other.fExecutableCommand;
      fArguments               = other.fArguments;
      fExecutableArgs          = other.fExecutableArgs;
      fAnalysisMacro           = other.fAnalysisMacro;
      fAnalysisSource          = other.fAnalysisSource;
      fValidationScript        = other.fValidationScript;
      fAdditionalRootLibs      = other.fAdditionalRootLibs;
      fAdditionalLibs          = other.fAdditionalLibs;
      fSplitMode               = other.fSplitMode;
      fAPIVersion              = other.fAPIVersion;
      fROOTVersion             = other.fROOTVersion;
      fAliROOTVersion          = other.fAliROOTVersion;
      fExternalPackages        = other.fExternalPackages;
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
      fTerminateFiles          = other.fTerminateFiles;
      fMergeExcludes           = other.fMergeExcludes;
      fIncludePath             = other.fIncludePath;
      fCloseSE                 = other.fCloseSE;
      fFriendChainName         = other.fFriendChainName;
      fJobTag                  = other.fJobTag;
      fOutputSingle            = other.fOutputSingle;
      fRunPrefix               = other.fRunPrefix;
      fProofCluster            = other.fProofCluster;
      fProofDataSet            = other.fProofDataSet;
      fFileForTestMode         = other.fFileForTestMode;
      fRootVersionForProof     = other.fRootVersionForProof;
      fAliRootMode             = other.fAliRootMode;
      fMergeDirName            = other.fMergeDirName;
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
      if (other.fModules) {
         fModules = new TObjArray();
         fModules->SetOwner();
         TIter next(other.fModules);
         AliAnalysisTaskCfg *mod, *crt;
         while ((crt=(AliAnalysisTaskCfg*)next())) {
            mod = new AliAnalysisTaskCfg(*crt);
            fModules->Add(mod);
         }
      }   
   }
   return *this;
}

//______________________________________________________________________________
void AliAnalysisAlien::AddModule(AliAnalysisTaskCfg *module)
{
// Adding a module. Checks if already existing. Becomes owned by this.
   if (!module) return;
   if (GetModule(module->GetName())) {
      Error("AddModule", "A module having the same name %s already added", module->GetName());
      return;
   }
   if (!fModules) {
      fModules = new TObjArray();
      fModules->SetOwner();
   }
   fModules->Add(module);
}

//______________________________________________________________________________
void AliAnalysisAlien::AddModules(TObjArray *list)
{
// Adding a list of modules. Checks if already existing. Becomes owned by this.
   TIter next(list);
   AliAnalysisTaskCfg *module;
   while ((module = (AliAnalysisTaskCfg*)next())) AddModule(module);
}   

//______________________________________________________________________________
Bool_t AliAnalysisAlien::CheckDependencies()
{
// Check if all dependencies are satisfied. Reorder modules if needed.
   Int_t nmodules = GetNmodules();
   if (!nmodules) {
      Warning("CheckDependencies", "No modules added yet to check their dependencies");
      return kTRUE;
   }   
   AliAnalysisTaskCfg *mod = 0;
   AliAnalysisTaskCfg *dep = 0;
   TString depname;
   Int_t i, j, k;
   for (i=0; i<nmodules; i++) {
      mod = (AliAnalysisTaskCfg*) fModules->At(i);
      Int_t ndeps = mod->GetNdeps();
      Int_t istart = i;
      for (j=0; j<ndeps; j++) {
         depname = mod->GetDependency(j);
         dep = GetModule(depname);
         if (!dep) {
            Error("CheckDependencies","Dependency %s not added for module %s",
                   depname.Data(), mod->GetName());
            return kFALSE;
         }
         if (dep->NeedsDependency(mod->GetName())) {
            Error("CheckDependencies","Modules %s and %s circularly depend on each other",
                   mod->GetName(), dep->GetName());
            return kFALSE;
         }                  
         Int_t idep = fModules->IndexOf(dep);
         // The dependency task must come first
         if (idep>i) {
            // Remove at idep and move all objects below up one slot
            // down to index i included.
            fModules->RemoveAt(idep);
            for (k=idep-1; k>=i; k--) fModules->AddAt(fModules->RemoveAt(k),k+1);
            fModules->AddAt(dep, i++);
         }
         //Redo from istart if dependencies were inserted
         if (i>istart) i=istart-1;
      }
   }
   return kTRUE;
}      

//______________________________________________________________________________
AliAnalysisManager *AliAnalysisAlien::CreateAnalysisManager(const char *name, const char *filename)
{
// Create the analysis manager and optionally execute the macro in filename.
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (mgr) return mgr;
   mgr = new AliAnalysisManager(name);
   mgr->SetGridHandler((AliAnalysisGrid*)this);
   if (strlen(filename)) {
      TString line = gSystem->ExpandPathName(filename);
      line.Prepend(".x ");
      gROOT->ProcessLine(line.Data());
   }
   return mgr;
}      
      
//______________________________________________________________________________
Int_t AliAnalysisAlien::GetNmodules() const
{
// Get number of modules.
   if (!fModules) return 0;
   return fModules->GetEntries();
}

//______________________________________________________________________________
AliAnalysisTaskCfg *AliAnalysisAlien::GetModule(const char *name)
{
// Get a module by name.
   if (!fModules) return 0;
   return (AliAnalysisTaskCfg*)fModules->FindObject(name);
}
   
//______________________________________________________________________________
Bool_t AliAnalysisAlien::LoadModule(AliAnalysisTaskCfg *mod)
{
// Load a given module.
   if (mod->IsLoaded()) return kTRUE;
   Int_t ndeps = mod->GetNdeps();
   TString depname;
   for (Int_t j=0; j<ndeps; j++) {
      depname = mod->GetDependency(j);
      AliAnalysisTaskCfg *dep = GetModule(depname);
      if (!dep) {
         Error("LoadModule","Dependency %s not existing for module %s",
                depname.Data(), mod->GetName());
         return kFALSE;
      }
      if (!LoadModule(dep)) {
         Error("LoadModule","Dependency %s for module %s could not be loaded",
                depname.Data(), mod->GetName());
         return kFALSE;
      }
   }
   // Load libraries for the module
   if (!mod->CheckLoadLibraries()) {
      Error("LoadModule", "Cannot load all libraries for module %s", mod->GetName());
      return kFALSE;
   }
   // Execute the macro
   if (mod->ExecuteMacro()<0) {
      Error("LoadModule", "Executing the macro %s with arguments: %s for module %s returned a negative value",
             mod->GetMacroName(), mod->GetMacroArgs(), mod->GetName());
      return kFALSE;
   }
   // Configure dependencies
   if (mod->GetConfigMacro() && mod->ExecuteConfigMacro()<0) {
      Error("LoadModule", "There was an error executing the deps config macro %s for module %s",
            mod->GetConfigMacro()->GetTitle(), mod->GetName());
      return kFALSE;
   }
   // Adjust extra libraries
   Int_t nlibs = mod->GetNlibs();
   TString lib;
   for (Int_t i=0; i<nlibs; i++) {
      lib = mod->GetLibrary(i);
      if (fAdditionalLibs.Contains(lib)) continue;
      lib = Form("lib%s.so", lib.Data());
      if (!fAdditionalLibs.IsNull()) fAdditionalLibs += " ";
      fAdditionalLibs += lib;
   }
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::GenerateTest(const char *name, const char *modname)
{
// Generate test macros for a single module or for the full train.
   fAdditionalLibs = "";
   if (strlen(modname)) {
      if (!CheckDependencies()) return kFALSE;
      AliAnalysisTaskCfg *mod = GetModule(modname);
      if (!mod) {
         Error("GenerateTest", "cannot generate test for inexistent module %s", modname);
         return kFALSE;
      }
      if (!LoadModule(mod)) return kFALSE;
   } else if (!LoadModules()) return kFALSE;
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr->InitAnalysis()) return kFALSE;
   mgr->PrintStatus();
   SetLocalTest(kTRUE);
   Int_t productionMode = fProductionMode;
   SetProductionMode();
   TString macro = fAnalysisMacro;
   TString executable = fExecutable;
   TString validation = fValidationScript;
   TString execCommand = fExecutableCommand;
   SetAnalysisMacro(Form("%s.C", name));
   SetExecutable(Form("%s.sh", name));
   SetExecutableCommand("aliroot -b -q ");
   SetValidationScript(Form("%s_validation.sh", name));
   WriteAnalysisFile();   
   WriteAnalysisMacro();
   WriteExecutable();
   WriteValidationScript();   
   SetLocalTest(kFALSE);
   SetProductionMode(productionMode);
   fAnalysisMacro = macro;
   fExecutable = executable;
   fExecutableCommand = execCommand;
   fValidationScript = validation;
   return kTRUE;   
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::LoadModules()
{
// Load all modules by executing the AddTask macros. Checks first the dependencies.
   fAdditionalLibs = "";
   Int_t nmodules = GetNmodules();
   if (!nmodules) {
      Warning("LoadModules", "No module to be loaded");
      return kTRUE;
   }   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("LoadModules", "No analysis manager created yet. Use CreateAnalysisManager first.");
      return kFALSE;
   }   
   if (!CheckDependencies()) return kFALSE;
   nmodules = GetNmodules();
   AliAnalysisTaskCfg *mod;
   for (Int_t imod=0; imod<nmodules; imod++) {
      mod = (AliAnalysisTaskCfg*)fModules->At(imod);
      if (!LoadModule(mod)) return kFALSE;
   }
   return kTRUE;
}      

//______________________________________________________________________________
void AliAnalysisAlien::SetRunPrefix(const char *prefix)
{
// Set the run number format. Can be a prefix or a format like "%09d"
   fRunPrefix = prefix;
   if (!fRunPrefix.Contains("%")) fRunPrefix += "%d";
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
   fRunNumbers += Form(fRunPrefix.Data(), run);
}   

//______________________________________________________________________________
void AliAnalysisAlien::AddRunList(const char* runList)
{
// Add several runs into the list of runs; they are expected to be separated by a blank character.  
  TString    sList = runList;
  TObjArray *list  = sList.Tokenize(" ");
  Int_t n = list->GetEntries();
  for (Int_t i = 0; i < n; i++) {
    TObjString *os = (TObjString*)list->At(i);
    AddRunNumber(os->GetString().Atoi());
  }
  delete list;
}

//______________________________________________________________________________
void AliAnalysisAlien::AddRunNumber(const char* run)
{
// Add a run number to the list of runs to be processed.
   TString runs = run;
   TObjString *os;
   TObjArray *arr = runs.Tokenize(" ");
   TIter next(arr);
   TString prefix; 
   prefix.Append(fRunPrefix, fRunPrefix.Index("%d"));
   while ((os=(TObjString*)next())){
       if (fRunNumbers.Length()) fRunNumbers += " ";
       fRunNumbers += Form("%s%s", prefix.Data(), os->GetString().Data());
   }
   delete arr;
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
void AliAnalysisAlien::AddExternalPackage(const char *package)
{
// Adds external packages w.r.t to the default ones (root,aliroot and gapi)
   if (fExternalPackages) fExternalPackages += " ";
   fExternalPackages += package;
}   
      
//______________________________________________________________________________
Bool_t AliAnalysisAlien::Connect()
{
// Try to connect to AliEn. User needs a valid token and /tmp/gclient_env_$UID sourced.
   if (gGrid && gGrid->IsConnected()) return kTRUE;
   if (fProductionMode) return kTRUE;
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
   if (gGrid->Mkdir(workdir, "-p")) {
      gGrid->Cd(fGridWorkingDir);
      Info("CdWork", "\n#####   Created alien working directory %s", fGridWorkingDir.Data());
   } else {
      Warning("CdWork", "Working directory %s cannot be created.\n Using %s instead.",
              workdir.Data(), homedir.Data());
      fGridWorkingDir = "";
   }          
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::CheckFileCopy(const char *alienpath)
{
// Check if file copying is possible.
   if (fProductionMode) return kTRUE;
   if (!Connect()) {
      Error("CheckFileCopy", "Not connected to AliEn. File copying cannot be tested.");
      return kFALSE;
   }
   Info("CheckFileCopy", "Checking possibility to copy files to your AliEn home directory... \
        \n +++ NOTE: You can disable this via: plugin->SetCheckCopy(kFALSE);");
   // Check if alien_CLOSE_SE is defined
   TString closeSE = gSystem->Getenv("alien_CLOSE_SE");
   if (!closeSE.IsNull()) {
      Info("CheckFileCopy", "Your current close storage is pointing to: \
           \n      alien_CLOSE_SE = \"%s\"", closeSE.Data());
   } else {
      Warning("CheckFileCopy", "Your current close storage is empty ! Depending on your location, file copying may fail.");
   }        
   // Check if grid directory exists.
   if (!DirectoryExists(alienpath)) {
      Error("CheckFileCopy", "Alien path %s does not seem to exist", alienpath);
      return kFALSE;
   }
   TFile f("plugin_test_copy", "RECREATE");
   // User may not have write permissions to current directory 
   if (f.IsZombie()) {
      Error("CheckFileCopy", "Cannot create local test file. Do you have write access to current directory: <%s> ?",
            gSystem->WorkingDirectory());
      return kFALSE;
   }
   f.Close();
   if (FileExists(Form("alien://%s/%s",alienpath, f.GetName()))) gGrid->Rm(Form("alien://%s/%s",alienpath, f.GetName()));
   if (!TFile::Cp(f.GetName(), Form("alien://%s/%s",alienpath, f.GetName()))) {
      Error("CheckFileCopy", "Cannot copy files to Alien destination: <%s> This may be temporary, or: \
           \n# 1. Make sure you have write permissions there. If this is the case: \
           \n# 2. Check the storage availability at: http://alimonitor.cern.ch/stats?page=SE/table \
           \n#    Do:           export alien_CLOSE_SE=\"working_disk_SE\" \
           \n#    To make this permanent put in in your .bashrc (in .alienshrc is not enough) \
           \n#    Redo token:   rm /tmp/x509up_u$UID then: alien-token-init <username>", alienpath);
      gSystem->Unlink(f.GetName());
      return kFALSE;
   }   
   gSystem->Unlink(f.GetName());
   gGrid->Rm(Form("%s%s",alienpath,f.GetName()));
   Info("CheckFileCopy", "### ...SUCCESS ###");
   return kTRUE;
}   

//______________________________________________________________________________
Bool_t AliAnalysisAlien::CheckInputData()
{
// Check validity of input data. If necessary, create xml files.
   if (fProductionMode) return kTRUE;
   if (!fInputFiles && !fRunNumbers.Length() && !fRunRange[0]) {
      if (!fGridDataDir.Length()) {
         Error("CkeckInputData", "AliEn path to base data directory must be set.\n = Use: SetGridDataDir()");
         return kFALSE;
      }
      if (fMergeViaJDL) {
         Error("CheckInputData", "Merging via jdl works only with run numbers, run range or provided xml");
         return kFALSE;
      }   
      Info("CheckInputData", "Analysis will make a single xml for base data directory %s",fGridDataDir.Data());
      if (fDataPattern.Contains("tag") && TestBit(AliAnalysisGrid::kTest))
         TObject::SetBit(AliAnalysisGrid::kUseTags, kTRUE); // ADDED (fix problem in determining the tag usage in test mode) 
      return kTRUE;
   }
   // Process declared files
   Bool_t isCollection = kFALSE;
   Bool_t isXml = kFALSE;
   Bool_t useTags = kFALSE;
   Bool_t checked = kFALSE;
   if (!TestBit(AliAnalysisGrid::kTest)) CdWork();
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
            isCollection = iscoll;
            isXml = isxml;
            useTags = usetags;
            TObject::SetBit(AliAnalysisGrid::kUseTags, useTags);
         } else {
            if ((iscoll != isCollection) || (isxml != isXml) || (usetags != useTags)) {
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
   if (isCollection) {
      Error("CheckInputData", "You are using raw AliEn collections as input. Cannot process run numbers.");
      return kFALSE;   
   }
   
   if (checked && !isXml) {
      Error("CheckInputData", "Cannot mix processing of full runs with non-xml files");
      return kFALSE;   
   }
   // Check validity of run number(s)
   TObjArray *arr;
   TObjString *os;
   TString format;
   Int_t nruns = 0;
   TString schunk, schunk2;
   TString path;
   if (!checked) {
      checked = kTRUE;
      useTags = fDataPattern.Contains("tag");
      TObject::SetBit(AliAnalysisGrid::kUseTags, useTags);
   }   
   if (useTags != fDataPattern.Contains("tag")) {
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
         if (useTags) msg += " using_tags: Yes";
         else          msg += " using_tags: No";
         Info("CheckDataType", "%s", msg.Data());
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
         format = Form("%%s/%s ", fRunPrefix.Data());
         path = Form(format.Data(), fGridDataDir.Data(), irun);
         if (!DirectoryExists(path)) {
            continue;
         }
         format = Form("%%s/%s.xml", fRunPrefix.Data());
         path = Form(format.Data(), workdir.Data(),irun);
         TString msg = "\n#####   file: ";
         msg += path;
         msg += " type: xml_collection;";
         if (useTags) msg += " using_tags: Yes";
         else          msg += " using_tags: No";
         Info("CheckDataType", "%s", msg.Data());
         if (fNrunsPerMaster<2) {
            format = Form("%s.xml", fRunPrefix.Data());
            AddDataFile(Form(format.Data(),irun));
         } else {
            nruns++;
            if (((nruns-1)%fNrunsPerMaster) == 0) {
               schunk = Form(fRunPrefix.Data(),irun);
            }
            format = Form("_%s.xml", fRunPrefix.Data());
            schunk2 = Form(format.Data(), irun);
            if ((nruns%fNrunsPerMaster)!=0 && irun != fRunRange[1]) continue;
            schunk += schunk2;
            AddDataFile(schunk);
         }   
      }
      if (!fInputFiles) {
         schunk += schunk2;
         AddDataFile(schunk);
      }   
   }
   return kTRUE;      
}   

//______________________________________________________________________________
Bool_t AliAnalysisAlien::CopyLocalDataset(const char *griddir, const char *pattern, Int_t nfiles, const char *output, const char *anchorfile)
{
// Copy data from the given grid directory according a pattern and make a local
// dataset.
   if (!Connect()) {
      Error("CopyLocalDataset", "Cannot copy local dataset with no grid connection");
      return kFALSE;
   }
   if (!DirectoryExists(griddir)) {
      Error("CopyLocalDataset", "Data directory %s not existing.", griddir);
      return kFALSE;
   }
   TString command = Form("find -z -l %d %s %s", nfiles, griddir, pattern);
   printf("Running command: %s\n", command.Data());
   TGridResult *res = gGrid->Command(command);
   Int_t nfound = res->GetEntries();
   if (!nfound) {
      Error("CopyLocalDataset", "No file found in <%s> having pattern <%s>", griddir, pattern);
      return kFALSE;
   }
   printf("... found %d files. Copying locally ...\n", nfound);
   // Copy files locally
   ofstream out;
   out.open(output, ios::out);
   TMap *map;
   TString turl, dirname, filename, temp;
   TString cdir = gSystem->WorkingDirectory();
   gSystem->MakeDirectory("data");
   gSystem->ChangeDirectory("data");
   for (Int_t i=0; i<nfound; i++) {
      map = (TMap*)res->At(i);
      turl = map->GetValue("turl")->GetName();
      filename = gSystem->BaseName(turl.Data());
      dirname = gSystem->DirName(turl.Data());
      dirname = gSystem->BaseName(dirname.Data());
      gSystem->MakeDirectory(dirname);
      if (TFile::Cp(turl, Form("file:./%s/%s", dirname.Data(), filename.Data()))) {
         if (strlen(anchorfile)) filename = Form("%s#%s", filename.Data(), anchorfile);
         out << cdir << "/data/" << Form("%s/%s", dirname.Data(), filename.Data()) << endl;
      }
   }
   gSystem->ChangeDirectory(cdir);
   delete res;
   return kTRUE;
}   

//______________________________________________________________________________
Bool_t AliAnalysisAlien::CreateDataset(const char *pattern)
{
// Create dataset for the grid data directory + run number.
   const Int_t gMaxEntries = 15000;
   if (fProductionMode || TestBit(AliAnalysisGrid::kOffline)) return kTRUE;
   if (!Connect()) {
      Error("CreateDataset", "Cannot create dataset with no grid connection");
      return kFALSE;
   }   

   // Cd workspace
   if (!TestBit(AliAnalysisGrid::kTest)) CdWork();
   TString workdir = gGrid->GetHomeDirectory();
   workdir += fGridWorkingDir;

   // Compose the 'find' command arguments
   TString format;
   TString command;
   TString options = "-x collection ";
   if (TestBit(AliAnalysisGrid::kTest)) options += Form("-l %d ", fNtestFiles);
   else options += Form("-l %d ", gMaxEntries);  // Protection for the find command
   TString conditions = "";
   Int_t nstart = 0;
   Int_t ncount = 0;
   Int_t stage = 0;
   TString file;
   TString path;
   Int_t nruns = 0;
   TString schunk, schunk2;
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
      while (1) {
         ncount = 0;
         stage++;
         if (gSystem->AccessPathName(file) || TestBit(AliAnalysisGrid::kTest) || fOverwriteMode) {
            command = "find ";
            command += Form("%s -o %d ",options.Data(), nstart);
            command += path;
            command += " ";
            command += pattern;
            command += conditions;
            printf("command: %s\n", command.Data());
            TGridResult *res = gGrid->Command(command);
            if (res) delete res;
            // Write standard output to file
            gROOT->ProcessLine(Form("gGrid->Stdout(); > __tmp%d__%s", stage, file.Data()));
            Bool_t hasGrep = (gSystem->Exec("grep --version 2>/dev/null > /dev/null")==0)?kTRUE:kFALSE;
            Bool_t nullFile = kFALSE;
            if (!hasGrep) {
                Warning("CreateDataset", "'grep' command not available on this system - cannot validate the result of the grid 'find' command");
            } else {
               nullFile = (gSystem->Exec(Form("grep -c /event __tmp%d__%s 2>/dev/null > __tmp__",stage,file.Data()))==0)?kFALSE:kTRUE;
               if (nullFile) {
                  Error("CreateDataset","Dataset %s produced by the previous find command is empty !", file.Data());
                  gSystem->Exec("rm -f __tmp*");
                  return kFALSE;
               }
               TString line;
               ifstream in;
               in.open("__tmp__");
               in >> line;
               in.close();
               gSystem->Exec("rm -f __tmp__");
               ncount = line.Atoi();
            }         
         }
         if (ncount == gMaxEntries) {
            Info("CreateDataset", "Dataset %s has more than 15K entries. Trying to merge...", file.Data());
            cadd = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"__tmp%d__%s\", 1000000);",stage,file.Data()));
            if (!cbase) cbase = cadd;
            else {
               cbase->Add(cadd);
               delete cadd;
            }   
            nstart += ncount;
         } else {
            if (cbase) {
               cadd = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"__tmp%d__%s\", 1000000);",stage,file.Data()));
               printf("... please wait - TAlienCollection::Add() scales badly...\n");
               cbase->Add(cadd);
               delete cadd;
               cbase->ExportXML(Form("file://%s", file.Data()),kFALSE,kFALSE, file, "Merged entries for a run");
               delete cbase; cbase = 0;               
            } else {
               TFile::Cp(Form("__tmp%d__%s",stage, file.Data()), file.Data());
            }
            gSystem->Exec("rm -f __tmp*");   
            Info("CreateDataset", "Created dataset %s with %d files", file.Data(), nstart+ncount);
            break;
         }
      }
      Bool_t fileExists = FileExists(file);
      if (!TestBit(AliAnalysisGrid::kTest) && (!fileExists || fOverwriteMode)) {
         // Copy xml file to alien space
         if (fileExists) gGrid->Rm(file);
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
   Bool_t nullResult = kTRUE;
   if (fRunNumbers.Length()) {
      TObjArray *arr = fRunNumbers.Tokenize(" ");
      TObjString *os;
      TIter next(arr);
      while ((os=(TObjString*)next())) {
         nstart = 0;
         stage = 0;
         path = Form("%s/%s/ ", fGridDataDir.Data(), os->GetString().Data());
         if (!DirectoryExists(path)) continue;
//         CdWork();
         if (TestBit(AliAnalysisGrid::kTest)) file = "wn.xml";
         else file = Form("%s.xml", os->GetString().Data());
         // If local collection file does not exist, create it via 'find' command.
         while (1) {
            ncount = 0;
            stage++;
            if (gSystem->AccessPathName(file) || TestBit(AliAnalysisGrid::kTest) || fOverwriteMode) {
               command = "find ";
               command +=  Form("%s -o %d ",options.Data(), nstart);
               command += path;
               command += pattern;
               command += conditions;
               TGridResult *res = gGrid->Command(command);
               if (res) delete res;
               // Write standard output to file
               gROOT->ProcessLine(Form("gGrid->Stdout(); > __tmp%d__%s", stage,file.Data()));
               Bool_t hasGrep = (gSystem->Exec("grep --version 2>/dev/null > /dev/null")==0)?kTRUE:kFALSE;
               Bool_t nullFile = kFALSE;
               if (!hasGrep) {
                  Warning("CreateDataset", "'grep' command not available on this system - cannot validate the result of the grid 'find' command");
               } else {
                  nullFile = (gSystem->Exec(Form("grep -c /event __tmp%d__%s 2>/dev/null > __tmp__",stage,file.Data()))==0)?kFALSE:kTRUE;
                  if (nullFile) {
                     Warning("CreateDataset","Dataset %s produced by: <%s> is empty !", file.Data(), command.Data());
                     gSystem->Exec("rm -f __tmp*");
                     fRunNumbers.ReplaceAll(os->GetString().Data(), "");
                     break;
                  }   
                  TString line;
                  ifstream in;
                  in.open("__tmp__");
                  in >> line;
                  in.close();
                  gSystem->Exec("rm -f __tmp__");   
                  ncount = line.Atoi();
               }
               nullResult = kFALSE;         
            }
            if (ncount == gMaxEntries) {
               Info("CreateDataset", "Dataset %s has more than 15K entries. Trying to merge...", file.Data());
               if (fNrunsPerMaster > 1) {
                  Error("CreateDataset", "File %s has more than %d entries. Please set the number of runs per master to 1 !", 
                          file.Data(),gMaxEntries);
                  return kFALSE;
               }           
               cadd = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"__tmp%d__%s\", 1000000);",stage,file.Data()));
               if (!cbase) cbase = cadd;
               else {
                  cbase->Add(cadd);
                  delete cadd;
               }   
               nstart += ncount;
            } else {
               if (cbase && fNrunsPerMaster<2) {
                  cadd = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"__tmp%d__%s\", 1000000);",stage,file.Data()));
                  printf("... please wait - TAlienCollection::Add() scales badly...\n");
                  cbase->Add(cadd);
                  delete cadd;
                  cbase->ExportXML(Form("file://%s", file.Data()),kFALSE,kFALSE, file, "Merged entries for a run");
                  delete cbase; cbase = 0;               
               } else {
                  TFile::Cp(Form("__tmp%d__%s",stage, file.Data()), file.Data());
               }
               gSystem->Exec("rm -f __tmp*");   
               Info("CreateDataset", "Created dataset %s with %d files", file.Data(), nstart+ncount);
               break;
            }
         }   
         if (TestBit(AliAnalysisGrid::kTest)) break;
         // Check if there is one run per master job.
         if (fNrunsPerMaster<2) {
            if (FileExists(file)) {
               if (fOverwriteMode) gGrid->Rm(file);
               else {
                  Info("CreateDataset", "\n#####   Dataset %s exist. Skipping creation...", file.Data());
                  continue;
               }   
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
               printf("   Merging collection <%s> into masterjob input...\n", file.Data());
               cbase->Add(cadd);
               delete cadd;
            }
            if ((nruns%fNrunsPerMaster)!=0 && os!=arr->Last()) {
               continue;
            }   
            schunk += Form("_%s.xml", os->GetString().Data());
            if (FileExists(schunk)) {               
               if (fOverwriteMode) gGrid->Rm(file);
               else {
                  Info("CreateDataset", "\n#####   Dataset %s exist. Skipping creation...", schunk.Data());
                  continue;
               }   
            }        
            printf("Exporting merged collection <%s> and copying to AliEn\n", schunk.Data());
            cbase->ExportXML(Form("file://%s", schunk.Data()),kFALSE,kFALSE, schunk, "Merged runs");
            TFile::Cp(Form("file:%s",schunk.Data()), Form("alien://%s/%s",workdir.Data(), schunk.Data()));
            if (!FileExists(schunk)) {
               Error("CreateDataset", "Copy command did NOT succeed for %s", schunk.Data());
               delete arr;
               return kFALSE;
            }
         }
      }   
      delete arr;
      if (nullResult) {
         Error("CreateDataset", "No valid dataset corresponding to the query!");
         return kFALSE;
      }
   } else {
      // Process a full run range.
      for (Int_t irun=fRunRange[0]; irun<=fRunRange[1]; irun++) {
         format = Form("%%s/%s ", fRunPrefix.Data());
         nstart = 0;
         stage = 0;
         path = Form(format.Data(), fGridDataDir.Data(), irun);
         if (!DirectoryExists(path)) continue;
//         CdWork();
         format = Form("%s.xml", fRunPrefix.Data());
         if (TestBit(AliAnalysisGrid::kTest)) file = "wn.xml";
         else file = Form(format.Data(), irun);
         if (FileExists(file) && fNrunsPerMaster<2 && !TestBit(AliAnalysisGrid::kTest)) {         
            if (fOverwriteMode) gGrid->Rm(file);
            else {
               Info("CreateDataset", "\n#####   Dataset %s exist. Skipping creation...", file.Data());
               continue;
            }   
         }
         // If local collection file does not exist, create it via 'find' command.
         while (1) {
            ncount = 0;
            stage++;
            if (gSystem->AccessPathName(file) || TestBit(AliAnalysisGrid::kTest) || fOverwriteMode) {
               command = "find ";
               command +=  Form("%s -o %d ",options.Data(), nstart);
               command += path;
               command += pattern;
               command += conditions;
               TGridResult *res = gGrid->Command(command);
               if (res) delete res;
               // Write standard output to file
               gROOT->ProcessLine(Form("gGrid->Stdout(); > __tmp%d__%s", stage,file.Data()));
               Bool_t hasGrep = (gSystem->Exec("grep --version 2>/dev/null > /dev/null")==0)?kTRUE:kFALSE;
               Bool_t nullFile = kFALSE;
               if (!hasGrep) {
                  Warning("CreateDataset", "'grep' command not available on this system - cannot validate the result of the grid 'find' command");
               } else {
                  nullFile = (gSystem->Exec(Form("grep -c /event __tmp%d__%s 2>/dev/null > __tmp__",stage,file.Data()))==0)?kFALSE:kTRUE;
                  if (nullFile) {
                     Warning("CreateDataset","Dataset %s produced by: <%s> is empty !", file.Data(), command.Data());
                     gSystem->Exec("rm -f __tmp*");
                     break;
                  }   
                  TString line;
                  ifstream in;
                  in.open("__tmp__");
                  in >> line;
                  in.close();
                  gSystem->Exec("rm -f __tmp__");   
                  ncount = line.Atoi();
               }
               nullResult = kFALSE;         
            }   
            if (ncount == gMaxEntries) {
               Info("CreateDataset", "Dataset %s has more than 15K entries. Trying to merge...", file.Data());
               if (fNrunsPerMaster > 1) {
                  Error("CreateDataset", "File %s has more than %d entries. Please set the number of runs per master to 1 !", 
                          file.Data(),gMaxEntries);
                  return kFALSE;
               }           
               cadd = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"__tmp%d__%s\", 1000000);",stage,file.Data()));
               if (!cbase) cbase = cadd;
               else {
                  cbase->Add(cadd);
                  delete cadd;
               }   
               nstart += ncount;
            } else {
               if (cbase && fNrunsPerMaster<2) {
                  cadd = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"__tmp%d__%s\", 1000000);",stage,file.Data()));
                  printf("... please wait - TAlienCollection::Add() scales badly...\n");
                  cbase->Add(cadd);
                  delete cadd;
                  cbase->ExportXML(Form("file://%s", file.Data()),kFALSE,kFALSE, file, "Merged entries for a run");
                  delete cbase; cbase = 0;               
               } else {
                  TFile::Cp(Form("__tmp%d__%s",stage, file.Data()), file.Data());
               }
               Info("CreateDataset", "Created dataset %s with %d files", file.Data(), nstart+ncount);
               break;
            }
         }   
         if (TestBit(AliAnalysisGrid::kTest)) break;
         // Check if there is one run per master job.
         if (fNrunsPerMaster<2) {
            if (FileExists(file)) {
               if (fOverwriteMode) gGrid->Rm(file);
               else {
                  Info("CreateDataset", "\n#####   Dataset %s exist. Skipping creation...", file.Data());
                  continue;
               }   
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
            if (FileExists(fInputFiles->At(nchunk)->GetName())) {
               if (fOverwriteMode) gGrid->Rm(fInputFiles->At(nchunk)->GetName());
               else continue;
            }   
            printf("   Merging collection <%s> into %d runs chunk...\n",file.Data(),fNrunsPerMaster);
            if (((nruns-1)%fNrunsPerMaster) == 0) {
               schunk = Form(fRunPrefix.Data(), irun);
               cbase = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"%s\", 1000000);",file.Data()));
            } else {
               cadd = (TGridCollection*)gROOT->ProcessLine(Form("new TAlienCollection(\"%s\", 1000000);",file.Data()));
               cbase->Add(cadd);
               delete cadd;
            }
            format = Form("%%s_%s.xml", fRunPrefix.Data());
            schunk2 = Form(format.Data(), schunk.Data(), irun);
            if ((nruns%fNrunsPerMaster)!=0 && irun!=fRunRange[1] && schunk2 != fInputFiles->Last()->GetName()) {
               continue;
            }   
            schunk = schunk2;
            if (FileExists(schunk)) {
               if (fOverwriteMode) gGrid->Rm(schunk);
               else {
                  Info("CreateDataset", "\n#####   Dataset %s exist. Skipping creation...", schunk.Data());
                  continue;
               }   
            }        
            printf("Exporting merged collection <%s> and copying to AliEn.\n", schunk.Data());
            cbase->ExportXML(Form("file://%s", schunk.Data()),kFALSE,kFALSE, schunk, "Merged runs");
            if (FileExists(schunk)) {
               if (fOverwriteMode) gGrid->Rm(schunk);
               else {
                  Info("CreateDataset", "\n#####   Dataset %s exist. Skipping copy...", schunk.Data());
                  continue;
               }   
            }   
            TFile::Cp(Form("file:%s",schunk.Data()), Form("alien://%s/%s",workdir.Data(), schunk.Data()));
            if (!FileExists(schunk)) {
               Error("CreateDataset", "Copy command did NOT succeed for %s", schunk.Data());
               return kFALSE;
            }
         }   
      }
      if (nullResult) {
         Error("CreateDataset", "No valid dataset corresponding to the query!");
         return kFALSE;
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
   if (fProductionMode || TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   Bool_t generate = kTRUE;
   if (TestBit(AliAnalysisGrid::kTest) || TestBit(AliAnalysisGrid::kSubmit)) generate = kFALSE;
   if (!Connect()) {
      Error("CreateJDL", "Alien connection required");
      return kFALSE;
   }   
   // Check validity of alien workspace
   TString workdir;
   if (!fProductionMode && !fGridWorkingDir.BeginsWith("/alice")) workdir = gGrid->GetHomeDirectory();
   if (!fProductionMode &&  !TestBit(AliAnalysisGrid::kTest)) CdWork();
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
         if (!fProductionMode) {
            if (!fGridOutputDir.Contains("/")) fGridOutputDir = Form("%s/%s", workdir.Data(), fGridOutputDir.Data());
            if (!DirectoryExists(fGridOutputDir)) {
               if (gGrid->Mkdir(fGridOutputDir,"-p")) {
                  Info("CreateJDL", "\n#####   Created alien output directory %s", fGridOutputDir.Data());
               } else {
                  Error("CreateJDL", "Could not create alien output directory %s", fGridOutputDir.Data());
                  // error = kTRUE;
               }
            } else {
               Warning("CreateJDL", "#### Output directory %s exists! If this contains old data, jobs will fail with ERROR_SV !!! ###", fGridOutputDir.Data());
            }   
            gGrid->Cd(workdir);
         }   
      }   
      // Exit if any error up to now
      if (error) return kFALSE;   
      // Set JDL fields
      if (!fUser.IsNull()) {
         fGridJDL->SetValue("User", Form("\"%s\"", fUser.Data()));
         fMergingJDL->SetValue("User", Form("\"%s\"", fUser.Data()));
      }   
      fGridJDL->SetExecutable(fExecutable, "This is the startup script");
      TString mergeExec = fExecutable;
      mergeExec.ReplaceAll(".sh", "_merge.sh");
      fMergingJDL->SetExecutable(mergeExec, "This is the startup script");
      mergeExec.ReplaceAll(".sh", ".C");
      fMergingJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(),mergeExec.Data()), "List of input files to be uploaded to workers");
      if (!fArguments.IsNull())
         fGridJDL->SetArguments(fArguments, "Arguments for the executable command");
      if (IsOneStageMerging()) fMergingJDL->SetArguments(fGridOutputDir);
      else {
         if (fProductionMode)  fMergingJDL->SetArguments("wn.xml $4");    // xml, stage
         else                  fMergingJDL->SetArguments("wn.xml $2");    // xml, stage
     }               

      fGridJDL->SetValue("TTL", Form("\"%d\"",fTTL));
      fGridJDL->SetDescription("TTL", Form("Time after which the job is killed (%d min.)", fTTL/60));
      fMergingJDL->SetValue("TTL", Form("\"%d\"",fTTL));
      fMergingJDL->SetDescription("TTL", Form("Time after which the job is killed (%d min.)", fTTL/60));
        
      if (fMaxInitFailed > 0) {
         fGridJDL->SetValue("MaxInitFailed", Form("\"%d\"",fMaxInitFailed));
         fGridJDL->SetDescription("MaxInitFailed", "Maximum number of first failing jobs to abort the master job");
      }   
      if (fSplitMaxInputFileNumber > 0) {
         fGridJDL->SetValue("SplitMaxInputFileNumber", Form("\"%d\"", fSplitMaxInputFileNumber));
         fGridJDL->SetDescription("SplitMaxInputFileNumber", "Maximum number of input files to be processed per subjob");
      }
      if (!IsOneStageMerging()) {
         fMergingJDL->SetValue("SplitMaxInputFileNumber", Form("\"%d\"",fMaxMergeFiles));
         fMergingJDL->SetDescription("SplitMaxInputFileNumber", "Maximum number of input files to be merged in one go");
      }   
      if (fSplitMode.Length()) {
         fGridJDL->SetValue("Split", Form("\"%s\"", fSplitMode.Data()));
         fGridJDL->SetDescription("Split", "We split per SE or file");
      }
      fMergingJDL->SetValue("Split", "\"se\""); 
      fMergingJDL->SetDescription("Split", "We split per SE for merging in stages");
      if (!fAliROOTVersion.IsNull()) {
         fGridJDL->AddToPackages("AliRoot", fAliROOTVersion,"VO_ALICE", "List of requested packages");
         fMergingJDL->AddToPackages("AliRoot", fAliROOTVersion, "VO_ALICE", "List of requested packages");
      }   
      if (!fROOTVersion.IsNull()) {
         fGridJDL->AddToPackages("ROOT", fROOTVersion);
         fMergingJDL->AddToPackages("ROOT", fROOTVersion);
      }   
      if (!fAPIVersion.IsNull()) {
         fGridJDL->AddToPackages("APISCONFIG", fAPIVersion);
         fMergingJDL->AddToPackages("APISCONFIG", fAPIVersion);
      }   
      if (!fExternalPackages.IsNull()) {
         arr = fExternalPackages.Tokenize(" ");
         TIter next(arr);
         while ((os=(TObjString*)next())) {
            TString pkgname = os->GetString();
            Int_t index = pkgname.Index("::");
            TString pkgversion = pkgname(index+2, pkgname.Length());
            pkgname.Remove(index);
            fGridJDL->AddToPackages(pkgname, pkgversion);
            fMergingJDL->AddToPackages(pkgname, pkgversion);
         }   
         delete arr;   
      }   
      fGridJDL->SetInputDataListFormat(fInputFormat, "Format of input data");
      fGridJDL->SetInputDataList("wn.xml", "Collection name to be processed on each worker node");
      fMergingJDL->SetInputDataListFormat(fInputFormat, "Format of input data");
      fMergingJDL->SetInputDataList("wn.xml", "Collection name to be processed on each worker node");
      fGridJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(), fAnalysisMacro.Data()), "List of input files to be uploaded to workers");
      TString analysisFile = fExecutable;
      analysisFile.ReplaceAll(".sh", ".root");
      fGridJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(),analysisFile.Data()));
      fMergingJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(),analysisFile.Data()));
      if (fAdditionalLibs.Length()) {
         arr = fAdditionalLibs.Tokenize(" ");
         TIter next(arr);
         while ((os=(TObjString*)next())) {
            if (os->GetString().Contains(".so")) continue;
            fGridJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(), os->GetString().Data()));
            fMergingJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(), os->GetString().Data()));
         }   
         delete arr;   
      }
      if (fPackages) {
         TIter next(fPackages);
         TObject *obj;
         while ((obj=next())) {
            fGridJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(), obj->GetName()));
            fMergingJDL->AddToInputSandbox(Form("LF:%s/%s", workdir.Data(), obj->GetName()));
         }
      }
      if (fOutputArchive.Length()) {
         arr = fOutputArchive.Tokenize(" ");
         TIter next(arr);
         Bool_t first = kTRUE;
         const char *comment = "Files to be archived";
         const char *comment1 = comment;
         while ((os=(TObjString*)next())) {
            if (!first) comment = NULL;
            if (!os->GetString().Contains("@") && fCloseSE.Length())
               fGridJDL->AddToOutputArchive(Form("%s@%s",os->GetString().Data(), fCloseSE.Data()), comment); 
            else
               fGridJDL->AddToOutputArchive(os->GetString(), comment);
            first = kFALSE;   
         }      
         delete arr;
         // Output archive for the merging jdl
         TString outputArchive;
         if (TestBit(AliAnalysisGrid::kDefaultOutputs)) {
            outputArchive = "log_archive.zip:std*@disk=1 ";
            // Add normal output files, extra files + terminate files
            TString files = GetListOfFiles("outextter");
            // Do not register merge excludes
            if (!fMergeExcludes.IsNull()) {
               arr = fMergeExcludes.Tokenize(" ");
               TIter next1(arr);
               while ((os=(TObjString*)next1())) {
                  files.ReplaceAll(Form("%s,",os->GetString().Data()),"");
                  files.ReplaceAll(os->GetString(),"");
               }   
               delete arr;
            }
            files.ReplaceAll(".root", "*.root");
            outputArchive += Form("root_archive.zip:%s,*.stat@disk=%d",files.Data(),fNreplicas);
         } else {
            TString files = fOutputArchive;
            files.ReplaceAll(".root", "*.root"); // nreplicas etc should be already atttached by use
            outputArchive = files;
         }   
         arr = outputArchive.Tokenize(" ");
         TIter next2(arr);
         comment = comment1;
         first = kTRUE;
         while ((os=(TObjString*)next2())) {
            if (!first) comment = NULL;
            TString currentfile = os->GetString();
            if (!currentfile.Contains("@") && fCloseSE.Length())
               fMergingJDL->AddToOutputArchive(Form("%s@%s",currentfile.Data(), fCloseSE.Data()), comment);
            else
               fMergingJDL->AddToOutputArchive(currentfile, comment);
            first = kFALSE;   
         }      
         delete arr;         
      }      
      arr = fOutputFiles.Tokenize(",");
      TIter next(arr);
      Bool_t first = kTRUE;
      const char *comment = "Files to be saved";
      while ((os=(TObjString*)next())) {
         // Ignore ouputs in jdl that are also in outputarchive
         TString sout = os->GetString();
         sout.ReplaceAll("*", "");
         sout.ReplaceAll(".root", "");
         if (sout.Index("@")>0) sout.Remove(sout.Index("@"));
         if (fOutputArchive.Contains(sout)) continue;
         if (!first) comment = NULL;
         if (!os->GetString().Contains("@") && fCloseSE.Length())
            fGridJDL->AddToOutputSandbox(Form("%s@%s",os->GetString().Data(), fCloseSE.Data()), comment); 
         else
            fGridJDL->AddToOutputSandbox(os->GetString(), comment);
         first = kFALSE;
         if (fMergeExcludes.Contains(sout)) continue;   
         if (!os->GetString().Contains("@") && fCloseSE.Length())
            fMergingJDL->AddToOutputSandbox(Form("%s@%s",os->GetString().Data(), fCloseSE.Data()), comment); 
         else
            fMergingJDL->AddToOutputSandbox(os->GetString(), comment);
      }   
      delete arr;
      fGridJDL->SetPrice((UInt_t)fPrice, "AliEn price for this job");
      fMergingJDL->SetPrice((UInt_t)fPrice, "AliEn price for this job");
      TString validationScript = fValidationScript;
      fGridJDL->SetValidationCommand(Form("%s/%s", workdir.Data(),validationScript.Data()), "Validation script to be run for each subjob");
      validationScript.ReplaceAll(".sh", "_merge.sh");
      fMergingJDL->SetValidationCommand(Form("%s/%s", workdir.Data(),validationScript.Data()), "Validation script to be run for each subjob");
      if (fMasterResubmitThreshold) {
         fGridJDL->SetValue("MasterResubmitThreshold", Form("\"%d%%\"", fMasterResubmitThreshold));
         fGridJDL->SetDescription("MasterResubmitThreshold", "Resubmit failed jobs until DONE rate reaches this percentage");
      }   
      // Write a jdl with 2 input parameters: collection name and output dir name.
      WriteJDL(copy);
   }
   // Copy jdl to grid workspace   
   if (copy) {
      // Check if an output directory was defined and valid
      if (!fGridOutputDir.Length()) {
         Error("CreateJDL", "You must define AliEn output directory");
         return kFALSE;
      } else {
         if (!fGridOutputDir.Contains("/")) fGridOutputDir = Form("%s/%s", workdir.Data(), fGridOutputDir.Data());
         if (!fProductionMode && !DirectoryExists(fGridOutputDir)) {
            if (gGrid->Mkdir(fGridOutputDir,"-p")) {
               Info("CreateJDL", "\n#####   Created alien output directory %s", fGridOutputDir.Data());
            } else {
               Error("CreateJDL", "Could not create alien output directory %s", fGridOutputDir.Data());
               return kFALSE;
            }
         }
         gGrid->Cd(workdir);
      }   
      if (TestBit(AliAnalysisGrid::kSubmit)) {
         TString mergeJDLName = fExecutable;
         mergeJDLName.ReplaceAll(".sh", "_merge.jdl");
         TString locjdl = Form("%s/%s", fGridOutputDir.Data(),fJDLName.Data());
         TString locjdl1 = Form("%s/%s", fGridOutputDir.Data(),mergeJDLName.Data());
         if (fProductionMode) {
            locjdl = Form("%s/%s", workdir.Data(),fJDLName.Data());
            locjdl1 = Form("%s/%s", workdir.Data(),mergeJDLName.Data());
         }   
         if (FileExists(locjdl)) gGrid->Rm(locjdl);
         if (FileExists(locjdl1)) gGrid->Rm(locjdl1);
         Info("CreateJDL", "\n#####   Copying JDL file <%s> to your AliEn output directory", fJDLName.Data());
         if (!copyLocal2Alien("CreateJDL", fJDLName, locjdl)) 
            Fatal("","Terminating");
//         TFile::Cp(Form("file:%s",fJDLName.Data()), Form("alien://%s", locjdl.Data()));
         if (fMergeViaJDL) {
            Info("CreateJDL", "\n#####   Copying merging JDL file <%s> to your AliEn output directory", mergeJDLName.Data());
//            TFile::Cp(Form("file:%s",mergeJDLName.Data()), Form("alien://%s", locjdl1.Data()));
            if (!copyLocal2Alien("CreateJDL", mergeJDLName.Data(), locjdl1)) 
               Fatal("","Terminating");
         }   
      }
      if (fAdditionalLibs.Length()) {
         arr = fAdditionalLibs.Tokenize(" ");
         TObjString *os;
         TIter next(arr);
         while ((os=(TObjString*)next())) {
            if (os->GetString().Contains(".so")) continue;
            Info("CreateJDL", "\n#####   Copying dependency: <%s> to your alien workspace", os->GetString().Data());
            if (FileExists(os->GetString())) gGrid->Rm(os->GetString());
//            TFile::Cp(Form("file:%s",os->GetString().Data()), Form("alien://%s/%s", workdir.Data(), os->GetString().Data()));
            if (!copyLocal2Alien("CreateJDL", os->GetString().Data(), 
                Form("%s/%s", workdir.Data(), os->GetString().Data())))
              Fatal("","Terminating");
         }   
         delete arr;   
      }
      if (fPackages) {
         TIter next(fPackages);
         TObject *obj;
         while ((obj=next())) {
            if (FileExists(obj->GetName())) gGrid->Rm(obj->GetName());
            Info("CreateJDL", "\n#####   Copying dependency: <%s> to your alien workspace", obj->GetName());
//            TFile::Cp(Form("file:%s",obj->GetName()), Form("alien://%s/%s", workdir.Data(), obj->GetName()));
            if (!copyLocal2Alien("CreateJDL",obj->GetName(), 
                Form("%s/%s", workdir.Data(), obj->GetName()))) 
              Fatal("","Terminating"); 
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
   TObject *os;
   TString workdir;
   if (!fProductionMode && !fGridWorkingDir.BeginsWith("/alice")) workdir = gGrid->GetHomeDirectory();
   workdir += fGridWorkingDir;
   TString stageName = "$2";
   if (fProductionMode) stageName = "$4";
   if (!fMergeDirName.IsNull()) {
     fMergingJDL->AddToInputDataCollection(Form("LF:$1/%s/Stage_%s.xml,nodownload",fMergeDirName.Data(),stageName.Data()), "Collection of files to be merged for current stage");
     fMergingJDL->SetOutputDirectory(Form("$1/%s/Stage_%s/#alien_counter_03i#",fMergeDirName.Data(),stageName.Data()), "Output directory");
   } else {
     fMergingJDL->AddToInputDataCollection(Form("LF:$1/Stage_%s.xml,nodownload",stageName.Data()), "Collection of files to be merged for current stage");
     fMergingJDL->SetOutputDirectory(Form("$1/Stage_%s/#alien_counter_03i#",stageName.Data()), "Output directory");
   }
   if (fProductionMode) {
      TIter next(fInputFiles);
      while ((os=next())) {
         fGridJDL->AddToInputDataCollection(Form("LF:%s,nodownload", os->GetName()), "Input xml collections");
      }
      if (!fOutputToRunNo)
         fGridJDL->SetOutputDirectory(Form("%s/#alien_counter_04i#", fGridOutputDir.Data()));
      else  
         fGridJDL->SetOutputDirectory(fGridOutputDir);
   } else {            
      if (!fRunNumbers.Length() && !fRunRange[0]) {
         // One jdl with no parameters in case input data is specified by name.
         TIter next(fInputFiles);
         while ((os=next()))
            fGridJDL->AddToInputDataCollection(Form("LF:%s,nodownload", os->GetName()), "Input xml collections");
         if (!fOutputSingle.IsNull())
            fGridJDL->SetOutputDirectory(Form("#alienfulldir#/../%s",fOutputSingle.Data()), "Output directory");
         else {
            fGridJDL->SetOutputDirectory(Form("%s/#alien_counter_03i#", fGridOutputDir.Data()), "Output directory");
            fMergingJDL->SetOutputDirectory(fGridOutputDir);         
         }   
      } else {
         // One jdl to be submitted with 2 input parameters: data collection name and output dir prefix
         fGridJDL->AddToInputDataCollection(Form("LF:%s/$1,nodownload", workdir.Data()), "Input xml collections");
         if (!fOutputSingle.IsNull()) {
            if (!fOutputToRunNo) fGridJDL->SetOutputDirectory(Form("#alienfulldir#/%s",fOutputSingle.Data()), "Output directory");
            else fGridJDL->SetOutputDirectory(Form("%s/$2",fGridOutputDir.Data()), "Output directory");
         } else {   
            fGridJDL->SetOutputDirectory(Form("%s/$2/#alien_counter_03i#", fGridOutputDir.Data()), "Output directory");
         }   
      }
   }
      
   // Generate the JDL as a string
   TString sjdl = fGridJDL->Generate();
   TString sjdl1 = fMergingJDL->Generate();
   // Final merge jdl
   if (!fMergeDirName.IsNull()) {
     fMergingJDL->SetOutputDirectory(Form("$1/%s",fMergeDirName.Data()), "Output directory");
     fMergingJDL->AddToInputSandbox(Form("LF:$1/%s/Stage_%s.xml",fMergeDirName.Data(),stageName.Data()));
   } else {  
     fMergingJDL->SetOutputDirectory("$1", "Output directory");
     fMergingJDL->AddToInputSandbox(Form("LF:$1/Stage_%s.xml",stageName.Data()));
   }  
   TString sjdl2 = fMergingJDL->Generate();
   Int_t index, index1;
   sjdl.ReplaceAll("\"LF:", "\n   \"LF:");
   sjdl.ReplaceAll("(member", "\n   (member");
   sjdl.ReplaceAll("\",\"VO_", "\",\n   \"VO_");
   sjdl.ReplaceAll("{", "{\n   ");
   sjdl.ReplaceAll("};", "\n};");
   sjdl.ReplaceAll("{\n   \n", "{\n");
   sjdl.ReplaceAll("\n\n", "\n");
   sjdl.ReplaceAll("OutputDirectory", "OutputDir");
   sjdl1.ReplaceAll("\"LF:", "\n   \"LF:");
   sjdl1.ReplaceAll("(member", "\n   (member");
   sjdl1.ReplaceAll("\",\"VO_", "\",\n   \"VO_");
   sjdl1.ReplaceAll("{", "{\n   ");
   sjdl1.ReplaceAll("};", "\n};");
   sjdl1.ReplaceAll("{\n   \n", "{\n");
   sjdl1.ReplaceAll("\n\n", "\n");
   sjdl1.ReplaceAll("OutputDirectory", "OutputDir");
   sjdl2.ReplaceAll("\"LF:", "\n   \"LF:");
   sjdl2.ReplaceAll("(member", "\n   (member");
   sjdl2.ReplaceAll("\",\"VO_", "\",\n   \"VO_");
   sjdl2.ReplaceAll("{", "{\n   ");
   sjdl2.ReplaceAll("};", "\n};");
   sjdl2.ReplaceAll("{\n   \n", "{\n");
   sjdl2.ReplaceAll("\n\n", "\n");
   sjdl2.ReplaceAll("OutputDirectory", "OutputDir");
   sjdl += "JDLVariables = \n{\n   \"Packages\",\n   \"OutputDir\"\n};\n";
   sjdl.Prepend(Form("Jobtag = {\n   \"comment:%s\"\n};\n", fJobTag.Data()));
   index = sjdl.Index("JDLVariables");
   if (index >= 0) sjdl.Insert(index, "\n# JDL variables\n");
   sjdl += "Workdirectorysize = {\"5000MB\"};";
   sjdl1 += "Workdirectorysize = {\"5000MB\"};";
   sjdl1 += "JDLVariables = \n{\n   \"Packages\",\n   \"OutputDir\"\n};\n";
   index = fJobTag.Index(":");
   if (index < 0) index = fJobTag.Length();
   TString jobTag = fJobTag;
   if (fProductionMode) jobTag.Insert(index,"_Stage$4");
   sjdl1.Prepend(Form("Jobtag = {\n   \"comment:%s_Merging\"\n};\n", jobTag.Data()));
   if (fProductionMode) {   
     sjdl1.Prepend("# Generated merging jdl (production mode) \
                    \n# $1 = full alien path to output directory to be merged \
                    \n# $2 = train number \
                    \n# $3 = production (like LHC10b) \
                    \n# $4 = merging stage \
                    \n# Stage_<n>.xml made via: find <OutputDir> *Stage<n-1>/*root_archive.zip\n");
     sjdl2.Prepend(Form("Jobtag = {\n   \"comment:%s_FinalMerging\"\n};\n", jobTag.Data()));
     sjdl2.Prepend("# Generated merging jdl \
                    \n# $1 = full alien path to output directory to be merged \
                    \n# $2 = train number \
                    \n# $3 = production (like LHC10b) \
                    \n# $4 = merging stage \
                    \n# Stage_<n>.xml made via: find <OutputDir> *Stage<n-1>/*root_archive.zip\n");
   } else {
     sjdl1.Prepend("# Generated merging jdl \
                    \n# $1 = full alien path to output directory to be merged \
                    \n# $2 = merging stage \
                    \n# xml made via: find <OutputDir> *Stage<n-1>/*root_archive.zip\n");
     sjdl2.Prepend(Form("Jobtag = {\n   \"comment:%s_FinalMerging\"\n};\n", jobTag.Data()));
     sjdl2.Prepend("# Generated merging jdl \
                    \n# $1 = full alien path to output directory to be merged \
                    \n# $2 = merging stage \
                    \n# xml made via: find <OutputDir> *Stage<n-1>/*root_archive.zip\n");
   }
   index = sjdl1.Index("JDLVariables");
   if (index >= 0) sjdl1.Insert(index, "\n# JDL variables\n");
   index = sjdl2.Index("JDLVariables");
   if (index >= 0) sjdl2.Insert(index, "\n# JDL variables\n");
   sjdl1 += "Workdirectorysize = {\"5000MB\"};";
   sjdl2 += "Workdirectorysize = {\"5000MB\"};";
   index = sjdl2.Index("Split =");
   if (index>=0) {
      index1 = sjdl2.Index("\n", index);
      sjdl2.Remove(index, index1-index+1);
   }
   index = sjdl2.Index("SplitMaxInputFileNumber");
   if (index>=0) {
      index1 = sjdl2.Index("\n", index);
      sjdl2.Remove(index, index1-index+1);
   }
   index = sjdl2.Index("InputDataCollection");
   if (index>=0) {
      index1 = sjdl2.Index(";", index);
      sjdl2.Remove(index, index1-index+1);
   }
   index = sjdl2.Index("InputDataListFormat");
   if (index>=0) {
      index1 = sjdl2.Index("\n", index);
      sjdl2.Remove(index, index1-index+1);
   }
   index = sjdl2.Index("InputDataList");
   if (index>=0) {
      index1 = sjdl2.Index("\n", index);
      sjdl2.Remove(index, index1-index+1);
   }
   sjdl2.ReplaceAll("wn.xml", Form("Stage_%s.xml",stageName.Data()));
   // Write jdl to file
   ofstream out;
   out.open(fJDLName.Data(), ios::out);
   if (out.bad()) {
      Error("WriteJDL", "Bad file name: %s", fJDLName.Data());
      return kFALSE;
   }
   out << sjdl << endl;
   out.close();
   TString mergeJDLName = fExecutable;
   mergeJDLName.ReplaceAll(".sh", "_merge.jdl");
   if (fMergeViaJDL) {
      ofstream out1;
      out1.open(mergeJDLName.Data(), ios::out);
      if (out1.bad()) {
         Error("WriteJDL", "Bad file name: %s", mergeJDLName.Data());
         return kFALSE;
      }
      out1 << sjdl1 << endl;
      out1.close();
      ofstream out2;
      TString finalJDL = mergeJDLName;
      finalJDL.ReplaceAll(".jdl", "_final.jdl");
      out2.open(finalJDL.Data(), ios::out);
      if (out2.bad()) {
         Error("WriteJDL", "Bad file name: %s", finalJDL.Data());
         return kFALSE;
      }
      out2 << sjdl2 << endl;
      out2.close();
   }   

   // Copy jdl to grid workspace   
   if (!copy) {
      Info("WriteJDL", "\n#####   You may want to review jdl:%s and analysis macro:%s before running in <submit> mode", fJDLName.Data(), fAnalysisMacro.Data());
   } else {
      TString locjdl = Form("%s/%s", fGridOutputDir.Data(),fJDLName.Data());
      TString locjdl1 = Form("%s/%s", fGridOutputDir.Data(),mergeJDLName.Data());
      TString finalJDL = mergeJDLName;
      finalJDL.ReplaceAll(".jdl", "_final.jdl");
      TString locjdl2 = Form("%s/%s", fGridOutputDir.Data(),finalJDL.Data());
      if (fProductionMode) {
         locjdl = Form("%s/%s", workdir.Data(),fJDLName.Data());
         locjdl1 = Form("%s/%s", workdir.Data(),mergeJDLName.Data());
         locjdl2 = Form("%s/%s", workdir.Data(),finalJDL.Data());
      }   
      if (FileExists(locjdl)) gGrid->Rm(locjdl);
      if (FileExists(locjdl1)) gGrid->Rm(locjdl1);
      if (FileExists(locjdl2)) gGrid->Rm(locjdl2);
      Info("WriteJDL", "\n#####   Copying JDL file <%s> to your AliEn output directory", fJDLName.Data());
//      TFile::Cp(Form("file:%s",fJDLName.Data()), Form("alien://%s", locjdl.Data()));
      if (!copyLocal2Alien("WriteJDL",fJDLName.Data(),locjdl.Data())) 
         Fatal("","Terminating");
      if (fMergeViaJDL) {
         Info("WriteJDL", "\n#####   Copying merging JDL files <%s> to your AliEn output directory", mergeJDLName.Data());
//         TFile::Cp(Form("file:%s",mergeJDLName.Data()), Form("alien://%s", locjdl1.Data()));
//         TFile::Cp(Form("file:%s",finalJDL.Data()), Form("alien://%s", locjdl2.Data()));
         if (!copyLocal2Alien("WriteJDL",mergeJDLName.Data(),locjdl1.Data()))
            Fatal("","Terminating");
         if (!copyLocal2Alien("WriteJDL",finalJDL.Data(),locjdl2.Data()))
           Fatal("","Terminating");
      }   
   } 
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::FileExists(const char *lfn)
{
// Returns true if file exists.
   if (!gGrid) return kFALSE;
   TString slfn = lfn;
   slfn.ReplaceAll("alien://","");
   TGridResult *res = gGrid->Ls(slfn);
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
void AliAnalysisAlien::CheckDataType(const char *lfn, Bool_t &isCollection, Bool_t &isXml, Bool_t &useTags)
{
// Check input data type.
   isCollection = kFALSE;
   isXml = kFALSE;
   useTags = kFALSE;
   if (!gGrid) {
      Error("CheckDataType", "No connection to grid");
      return;
   }
   isCollection = IsCollection(lfn);
   TString msg = "\n#####   file: ";
   msg += lfn;
   if (isCollection) {
      msg += " type: raw_collection;";
   // special treatment for collections
      isXml = kFALSE;
      // check for tag files in the collection
      TGridResult *res = gGrid->Command(Form("listFilesFromCollection -z -v %s",lfn), kFALSE);
      if (!res) {
         msg += " using_tags: No (unknown)";
         Info("CheckDataType", "%s", msg.Data());
         return;
      }   
      const char* typeStr = res->GetKey(0, "origLFN");
      if (!typeStr || !strlen(typeStr)) {
         msg += " using_tags: No (unknown)";
         Info("CheckDataType", "%s", msg.Data());
         return;
      }   
      TString file = typeStr;
      useTags = file.Contains(".tag");
      if (useTags) msg += " using_tags: Yes";
      else          msg += " using_tags: No";
      Info("CheckDataType", "%s", msg.Data());
      return;
   }
   TString slfn(lfn);
   slfn.ToLower();
   isXml = slfn.Contains(".xml");
   if (isXml) {
   // Open xml collection and check if there are tag files inside
      msg += " type: xml_collection;";
      TGridCollection *coll = (TGridCollection*)gROOT->ProcessLine(Form("TAlienCollection::Open(\"alien://%s\",1);",lfn));
      if (!coll) {
         msg += " using_tags: No (unknown)";
         Info("CheckDataType", "%s", msg.Data());
         return;
      }   
      TMap *map = coll->Next();
      if (!map) {
         msg += " using_tags: No (unknown)";
         Info("CheckDataType", "%s", msg.Data());
         return;
      }   
      map = (TMap*)map->GetValue("");
      TString file;
      if (map && map->GetValue("name")) file = map->GetValue("name")->GetName();
      useTags = file.Contains(".tag");
      delete coll;
      if (useTags) msg += " using_tags: Yes";
      else          msg += " using_tags: No";
      Info("CheckDataType", "%s", msg.Data());
      return;
   }
   useTags = slfn.Contains(".tag");
   if (slfn.Contains(".root")) msg += " type: root file;";
   else                        msg += " type: unknown file;";
   if (useTags) msg += " using_tags: Yes";
   else          msg += " using_tags: No";
   Info("CheckDataType", "%s", msg.Data());
}

//______________________________________________________________________________
void AliAnalysisAlien::EnablePackage(const char *package)
{
// Enables a par file supposed to exist in the current directory.
   TString pkg(package);
   pkg.ReplaceAll(".par", "");
   pkg += ".par";
   if (gSystem->AccessPathName(pkg)) {
      Fatal("EnablePackage", "Package %s not found", pkg.Data());
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
TChain *AliAnalysisAlien::GetChainForTestMode(const char *treeName) const
{
// Make a tree from files having the location specified in fFileForTestMode. 
// Inspired from JF's CreateESDChain.
   if (fFileForTestMode.IsNull()) {
      Error("GetChainForTestMode", "For proof test mode please use SetFileForTestMode() pointing to a file that contains data file locations.");
      return NULL;
   }
   if (gSystem->AccessPathName(fFileForTestMode)) {
      Error("GetChainForTestMode", "File not found: %s", fFileForTestMode.Data());
      return NULL;
   }   
   // Open the file
   ifstream in;
   in.open(fFileForTestMode);
   Int_t count = 0;
    // Read the input list of files and add them to the chain
    TString line;
    TChain *chain = new TChain(treeName);
    while (in.good())
    {
      in >> line;
      if (line.IsNull()) continue;
      if (count++ == fNtestFiles) break;
      TString esdFile(line);
      TFile *file = TFile::Open(esdFile);
      if (file) {
         if (!file->IsZombie()) chain->Add(esdFile);
         file->Close();
      } else {
         Error("GetChainforTestMode", "Skipping un-openable file: %s", esdFile.Data());
      }   
    }
    in.close();
    if (!chain->GetListOfFiles()->GetEntries()) {
       Error("GetChainForTestMode", "No file from %s could be opened", fFileForTestMode.Data());
       delete chain;
       return NULL;
    }
//    chain->ls();
    return chain;
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
      pid = gROOT->ProcessLine(Form("atoi(((TAlienJobStatus*)%p)->GetKey(\"queueId\"));", status));
      if (pid<jobidstart) continue;
      if (pid == lastid) {
         gROOT->ProcessLine(Form("sprintf((char*)%p,((TAlienJobStatus*)%p)->GetKey(\"status\"));",mstatus, status));
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
Bool_t AliAnalysisAlien::IsSingleOutput() const
{
// Check if single-ouput option is on.
   return (!fOutputSingle.IsNull());
}
   
//______________________________________________________________________________
void AliAnalysisAlien::Print(Option_t *) const
{
// Print current plugin settings.
   printf("### AliEn analysis plugin current settings ###\n");
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (mgr && mgr->IsProofMode()) {
      TString proofType = "=   PLUGIN IN PROOF MODE ON CLUSTER:_________________";
      if (TestBit(AliAnalysisGrid::kTest))
         proofType = "=   PLUGIN IN PROOF LITE MODE ON CLUSTER:____________";
      printf("%s %s\n", proofType.Data(), fProofCluster.Data());
      if (!fProofDataSet.IsNull())
      printf("=   Requested data set:___________________________ %s\n", fProofDataSet.Data());
      if (fProofReset==1)
      printf("=   Soft reset signal will be send to master______ CHANGE BEHAVIOR AFTER COMPLETION\n");      
      if (fProofReset>1)   
      printf("=   Hard reset signal will be send to master______ CHANGE BEHAVIOR AFTER COMPLETION\n");      
      if (!fRootVersionForProof.IsNull())
      printf("=   ROOT version requested________________________ %s\n", fRootVersionForProof.Data());
      else
      printf("=   ROOT version requested________________________ default\n");
      printf("=   AliRoot version requested_____________________ %s\n", fAliROOTVersion.Data());
      if (!fAliRootMode.IsNull())
      printf("=   Requested AliRoot mode________________________ %s\n", fAliRootMode.Data());  
      if (fNproofWorkers)
      printf("=   Number of PROOF workers limited to____________ %d\n", fNproofWorkers);
      if  (fNproofWorkersPerSlave)
      printf("=   Maximum number of workers per slave___________ %d\n", fNproofWorkersPerSlave);
      if (TestSpecialBit(kClearPackages))
      printf("=   ClearPackages requested...\n");
      if (fIncludePath.Data())
      printf("=   Include path for runtime task compilation: ___ %s\n", fIncludePath.Data());
      printf("=   Additional libs to be loaded or souces to be compiled runtime: <%s>\n",fAdditionalLibs.Data());
      if (fPackages && fPackages->GetEntries()) {
         TIter next(fPackages);
         TObject *obj;
         TString list;
         while ((obj=next())) list += obj->GetName();
         printf("=   Par files to be used: ________________________ %s\n", list.Data());
      } 
      if (TestSpecialBit(kProofConnectGrid))
      printf("=   Requested PROOF connection to grid\n");
      return;
   }
   printf("=   OverwriteMode:________________________________ %d\n", fOverwriteMode);
   if (fOverwriteMode) {
      printf("***** NOTE: Overwrite mode will overwrite the input generated datasets and partial results from previous analysis. \
            \n*****       To disable, use: plugin->SetOverwriteMode(kFALSE);\n");
   }
   printf("=   Copy files to grid: __________________________ %s\n", (IsUseCopy())?"YES":"NO");
   printf("=   Check if files can be copied to grid: ________ %s\n", (IsCheckCopy())?"YES":"NO");
   printf("=   Production mode:______________________________ %d\n", fProductionMode);
   printf("=   Version of API requested: ____________________ %s\n", fAPIVersion.Data());
   printf("=   Version of ROOT requested: ___________________ %s\n", fROOTVersion.Data());
   printf("=   Version of AliRoot requested: ________________ %s\n", fAliROOTVersion.Data());
   if (fUser.Length()) 
   printf("=   User running the plugin: _____________________ %s\n", fUser.Data());
   printf("=   Grid workdir relative to user $HOME: _________ %s\n", fGridWorkingDir.Data());
   printf("=   Grid output directory relative to workdir: ___ %s\n", fGridOutputDir.Data());
   printf("=   Data base directory path requested: __________ %s\n", fGridDataDir.Data());
   printf("=   Data search pattern: _________________________ %s\n", fDataPattern.Data());
   printf("=   Input data format: ___________________________ %s\n", fInputFormat.Data());
   if (fRunNumbers.Length()) 
   printf("=   Run numbers to be processed: _________________ %s\n", fRunNumbers.Data());
   if (fRunRange[0])
   printf("=   Run range to be processed: ___________________ %d-%d\n", fRunRange[0], fRunRange[1]);
   if (!fRunRange[0] && !fRunNumbers.Length()) {
      TIter next(fInputFiles);
      TObject *obj;
      TString list;
      while ((obj=next())) list += obj->GetName();
      printf("=   Input files to be processed: _________________ %s\n", list.Data());
   }
   if (TestBit(AliAnalysisGrid::kTest))
   printf("=   Number of input files used in test mode: _____ %d\n", fNtestFiles);
   printf("=   List of output files to be registered: _______ %s\n", fOutputFiles.Data());
   printf("=   List of outputs going to be archived: ________ %s\n", fOutputArchive.Data());
   printf("=   List of outputs that should not be merged: ___ %s\n", fMergeExcludes.Data());
   printf("=   List of outputs produced during Terminate: ___ %s\n", fTerminateFiles.Data());
   printf("=====================================================================\n");
   printf("=   Job price: ___________________________________ %d\n", fPrice);
   printf("=   Time to live (TTL): __________________________ %d\n", fTTL);
   printf("=   Max files per subjob: ________________________ %d\n", fSplitMaxInputFileNumber);
   if (fMaxInitFailed>0) 
   printf("=   Max number of subjob fails to kill: __________ %d\n", fMaxInitFailed);
   if (fMasterResubmitThreshold>0) 
   printf("=   Resubmit master job if failed subjobs >_______ %d\n", fMasterResubmitThreshold);
   printf("=   Number of replicas for the output files_______ %d\n", fNreplicas);
   if (fNrunsPerMaster>0)
   printf("=   Number of runs per master job: _______________ %d\n", fNrunsPerMaster);
   printf("=   Number of files in one chunk to be merged: ___ %d\n", fMaxMergeFiles);
   printf("=   Name of the generated execution script: ______ %s\n", fExecutable.Data());
   printf("=   Executable command: __________________________ %s\n", fExecutableCommand.Data());
   if (fArguments.Length()) 
   printf("=   Arguments for the execution script: __________ %s\n",fArguments.Data());
   if (fExecutableArgs.Length()) 
   printf("=   Arguments after macro name in executable______ %s\n",fExecutableArgs.Data());
   printf("=   Name of the generated analysis macro: ________ %s\n",fAnalysisMacro.Data());
   printf("=   User analysis files to be deployed: __________ %s\n",fAnalysisSource.Data());
   printf("=   Additional libs to be loaded or souces to be compiled runtime: <%s>\n",fAdditionalLibs.Data());
   printf("=   Master jobs split mode: ______________________ %s\n",fSplitMode.Data());
   if (fDatasetName)
   printf("=   Custom name for the dataset to be created: ___ %s\n", fDatasetName.Data());
   printf("=   Name of the generated JDL: ___________________ %s\n", fJDLName.Data());
   if (fIncludePath.Data())
   printf("=   Include path for runtime task compilation: ___ %s\n", fIncludePath.Data());
   if (fCloseSE.Length())
   printf("=   Force job outputs to storage element: ________ %s\n", fCloseSE.Data());
   if (fFriendChainName.Length())
   printf("=   Open friend chain file on worker: ____________ %s\n", fFriendChainName.Data());
   if (fPackages && fPackages->GetEntries()) {
      TIter next(fPackages);
      TObject *obj;
      TString list;
      while ((obj=next())) list += obj->GetName();
      printf("=   Par files to be used: ________________________ %s\n", list.Data());
   }   
}

//______________________________________________________________________________
void AliAnalysisAlien::SetDefaults()
{
// Set default values for everything. What cannot be filled will be left empty.
   if (fGridJDL) delete fGridJDL;
   fGridJDL = (TGridJDL*)gROOT->ProcessLine("new TAlienJDL()");
   fMergingJDL = (TGridJDL*)gROOT->ProcessLine("new TAlienJDL()");
   fPrice                      = 1;
   fTTL                        = 30000;
   fSplitMaxInputFileNumber    = 100;
   fMaxInitFailed              = 0;
   fMasterResubmitThreshold    = 0;
   fNtestFiles                 = 10;
   fNreplicas                  = 2;
   fRunRange[0]                = 0;
   fRunRange[1]                = 0;
   fRunPrefix                  = "%d";
   fNrunsPerMaster             = 1;
   fMaxMergeFiles              = 100;
   fRunNumbers                 = "";
   fExecutable                 = "analysis.sh";
   fExecutableCommand          = "root -b -q";
   fArguments                  = "";
   fExecutableArgs             = "";
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
   fOutputArchive              = "log_archive.zip:std*@disk=1 root_archive.zip:*.root@disk=2";
   fOutputFiles                = "";  // Like "AliAODs.root histos.root"
   fInputFormat                = "xml-single";
   fJDLName                    = "analysis.jdl";
   fJobTag                     = "Automatically generated analysis JDL";
   fMergeExcludes              = "";
   fMergeViaJDL                = 0;
   SetUseCopy(kTRUE);
   SetCheckCopy(kTRUE);
   SetDefaultOutputs(kTRUE);
   fOverwriteMode              = 1;
}   

//______________________________________________________________________________
Bool_t AliAnalysisAlien::CheckMergedFiles(const char *filename, const char *aliendir, Int_t nperchunk, const char *jdl)
{
// Checks current merge stage, makes xml for the next stage, counts number of files, submits next stage.
   // First check if the result is already in the output directory.
   if (FileExists(Form("%s/%s",aliendir,filename))) {
      printf("Final merged results found. Not merging again.\n");
      return kFALSE;
   }
   // Now check the last stage done.
   Int_t stage = 0;
   while (1) {
      if (!FileExists(Form("%s/Stage_%d.xml",aliendir, stage+1))) break;
      stage++;
   }
   // Next stage of merging
   stage++;
   TString pattern = "*root_archive.zip";
   if (stage>1) pattern = Form("Stage_%d/*root_archive.zip", stage-1);
   TGridResult *res = gGrid->Command(Form("find -x Stage_%d %s %s", stage, aliendir, pattern.Data()));
   if (res) delete res;
   // Write standard output to file
   gROOT->ProcessLine(Form("gGrid->Stdout(); > %s", Form("Stage_%d.xml",stage)));
   // Count the number of files inside
   ifstream ifile;
   ifile.open(Form("Stage_%d.xml",stage));
   if (!ifile.good()) {
      ::Error("CheckMergedFiles", "Could not redirect result of the find command to file %s", Form("Stage_%d.xml",stage));
      return kFALSE;
   }   
   TString line;
   Int_t nfiles = 0;
   while (!ifile.eof()) {
      ifile >> line;
      if (line.Contains("/event")) nfiles++;
   }
   ifile.close();
   if (!nfiles) {
      ::Error("CheckMergedFiles", "Cannot start Stage_%d merging since Stage_%d did not produced yet output", stage, stage-1);
      return kFALSE;
   } else {
      printf("=== Stage_%d produced %d files\n", stage-1, nfiles);
   }   
   // Copy the file in the output directory
   printf("===> Copying collection %s in the output directory %s\n", Form("Stage_%d.xml",stage), aliendir);
   TFile::Cp(Form("Stage_%d.xml",stage), Form("alien://%s/Stage_%d.xml",aliendir,stage));
   if (!copyLocal2Alien("CheckMergedFiles", Form("Stage_%d.xml",stage), 
        Form("%s/Stage_%d.xml",aliendir,stage))) Fatal("","Terminating");
   // Check if this is the last stage to be done.
   Bool_t laststage = (nfiles<nperchunk);
   if (fMaxMergeStages && stage>=fMaxMergeStages) laststage = kTRUE;
   if (laststage) {
      printf("### Submiting final merging stage %d\n", stage);
      TString finalJDL = jdl;
      finalJDL.ReplaceAll(".jdl", "_final.jdl");
      TString query = Form("submit %s %s %d", finalJDL.Data(), aliendir, stage);
      Int_t jobId = SubmitSingleJob(query);
      if (!jobId) return kFALSE;      
   } else {
      printf("### Submiting merging stage %d\n", stage);
      TString query = Form("submit %s %s %d", jdl, aliendir, stage);
      Int_t jobId = SubmitSingleJob(query);
      if (!jobId) return kFALSE;           
   }
   return kTRUE;   
}        

//______________________________________________________________________________
AliAnalysisManager *AliAnalysisAlien::LoadAnalysisManager(const char *fname)
{
// Loat the analysis manager from a file.
   TFile *file = TFile::Open(fname);
   if (!file) {
      ::Error("LoadAnalysisManager", "Cannot open file %s", fname);
      return 0;
   }   
   TIter nextkey(file->GetListOfKeys());
   AliAnalysisManager *mgr = 0;
   TKey *key;
   while ((key=(TKey*)nextkey())) {
      if (!strcmp(key->GetClassName(), "AliAnalysisManager"))
         mgr = (AliAnalysisManager*)file->Get(key->GetName());
   }
   if (!mgr) 
      ::Error("LoadAnalysisManager", "No analysis manager found in file %s", fname);
   return mgr;
}      

//______________________________________________________________________________
Int_t AliAnalysisAlien::SubmitSingleJob(const char *query)
{
// Submits a single job corresponding to the query and returns job id. If 0 submission failed.
   if (!gGrid) return 0;
   printf("=> %s ------> ",query);
   TGridResult *res = gGrid->Command(query);
   if (!res) return 0;
   TString jobId = res->GetKey(0,"jobId");
   delete res;
   if (jobId.IsNull()) {
      printf("submission failed. Reason:\n");
      gGrid->Stdout();
      gGrid->Stderr();
      ::Error("SubmitSingleJob", "Your query %s could not be submitted", query);
      return 0;
   }
   printf(" Job id: %s\n", jobId.Data());
   return atoi(jobId);
}  

//______________________________________________________________________________
Bool_t AliAnalysisAlien::MergeOutput(const char *output, const char *basedir, Int_t nmaxmerge, Int_t stage)
{
// Merge given output files from basedir. Basedir can be an alien output directory
// but also an xml file with root_archive.zip locations. The file merger will merge nmaxmerge
// files in a group (ignored for xml input). Merging can be done in stages:
// stage=0 : will merge all existing files in a single stage, supporting resume if run locally
// stage=1 : works with an xml of all root_archive.zip in the output directory
// stage>1 : works with an xml of all root_archive.zip in the Stage_<n-1> directory
   TString outputFile = output;
   TString command;
   TString outputChunk;
   TString previousChunk = "";
   TObjArray *listoffiles = new TObjArray();
//   listoffiles->SetOwner();
   Int_t countChunk = 0;
   Int_t countZero = nmaxmerge;
   Bool_t merged = kTRUE;
   Int_t index = outputFile.Index("@");
   if (index > 0) outputFile.Remove(index);
   TString inputFile = outputFile;
   TString sbasedir = basedir;
   if (sbasedir.Contains(".xml")) {
      // Merge files pointed by the xml - ignore nmaxmerge and set ichunk to 0
      nmaxmerge = 9999999;
      TGridCollection *coll = (TGridCollection*)gROOT->ProcessLine(Form("TAlienCollection::Open(\"%s\");", basedir));
      if (!coll) {
         ::Error("MergeOutput", "Input XML collection empty.");
         return kFALSE;
      }
      // Iterate grid collection
      while (coll->Next()) {
         TString fname = gSystem->DirName(coll->GetTURL());
         fname += "/";
         fname += inputFile;      
         listoffiles->Add(new TNamed(fname.Data(),""));
      }   
   } else {   
      command = Form("find %s/ *%s", basedir, inputFile.Data());
      printf("command: %s\n", command.Data());
      TGridResult *res = gGrid->Command(command);
      if (!res) {
         ::Error("MergeOutput","No result for the find command\n");
         delete listoffiles;
         return kFALSE;
      }     
      TIter nextmap(res);
      TMap *map = 0;
      while ((map=(TMap*)nextmap())) {
         TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
         if (!objs || !objs->GetString().Length()) {
            // Nothing found - skip this output
            delete res;
            delete listoffiles;
            return kFALSE;
         }
         listoffiles->Add(new TNamed(objs->GetName(),""));
      }
      delete res;
   }
   if (!listoffiles->GetEntries()) {
      ::Error("MergeOutput","No result for the find command\n");
      delete listoffiles;
      return kFALSE;
   }     

   TFileMerger *fm = 0;
   TIter next0(listoffiles);
   TObjArray *listoffilestmp = new TObjArray();
   listoffilestmp->SetOwner();
   TObject *nextfile;
   TString snextfile;
   // Keep only the files at upper level
   Int_t countChar = 0;
   while ((nextfile=next0())) {
      snextfile = nextfile->GetName();
      Int_t crtCount = snextfile.CountChar('/');
      if (nextfile == listoffiles->First()) countChar = crtCount;
      if (crtCount < countChar) countChar = crtCount;
   }
   next0.Reset();
   while ((nextfile=next0())) {
      snextfile = nextfile->GetName();
      Int_t crtCount = snextfile.CountChar('/');
      if (crtCount > countChar) {
         delete nextfile;
         continue;
      }   
      listoffilestmp->Add(nextfile);
   }
   delete listoffiles;
   listoffiles = listoffilestmp;  // Now contains 'good' files
   listoffiles->Print();
   TIter next(listoffiles);   
   // Check if there is a merge operation to resume. Works only for stage 0 or 1.
   outputChunk = outputFile;
   outputChunk.ReplaceAll(".root", "_*.root");
   // Check for existent temporary merge files
   // Check overwrite mode and remove previous partial results if needed
   // Preserve old merging functionality for stage 0.
   if (stage==0) {
      if (!gSystem->Exec(Form("ls %s 2>/dev/null", outputChunk.Data()))) {
         while (1) {
            // Skip as many input files as in a chunk
            for (Int_t counter=0; counter<nmaxmerge; counter++) {
               nextfile = next();
               if (!nextfile) {
                  ::Error("MergeOutput", "Mismatch found. Please remove partial merged files from local dir.");
                  delete listoffiles;
                  return kFALSE;
               }   
               snextfile = nextfile->GetName();
            }
            outputChunk = outputFile;
            outputChunk.ReplaceAll(".root", Form("_%04d.root", countChunk));
            countChunk++;
            if (gSystem->AccessPathName(outputChunk)) continue;
            // Merged file with chunks up to <countChunk> found
            ::Info("MergeOutput", "Resume merging of <%s> from <%s>\n", outputFile.Data(), outputChunk.Data());
            previousChunk = outputChunk;
            break;
         }
      }   
      countZero = nmaxmerge;
   
      while ((nextfile=next())) {
         snextfile = nextfile->GetName();
         // Loop 'find' results and get next LFN
         if (countZero == nmaxmerge) {
            // First file in chunk - create file merger and add previous chunk if any.
            fm = new TFileMerger(kFALSE);
            fm->SetFastMethod(kTRUE);
            if (previousChunk.Length()) fm->AddFile(previousChunk.Data());
            outputChunk = outputFile;
            outputChunk.ReplaceAll(".root", Form("_%04d.root", countChunk));
         }
         // If last file found, put merged results in the output file
         if (nextfile == listoffiles->Last()) outputChunk = outputFile;
         // Add file to be merged and decrement chunk counter.
         fm->AddFile(snextfile);
         countZero--;
         if (countZero==0 || nextfile == listoffiles->Last()) {            
            if (!fm->GetMergeList() || !fm->GetMergeList()->GetSize()) {
            // Nothing found - skip this output
               ::Warning("MergeOutput", "No <%s> files found.", inputFile.Data());
               merged = kFALSE;
               break;
            }
            fm->OutputFile(outputChunk);
            // Merge the outputs, then go to next chunk      
            if (!fm->Merge()) {
               ::Error("MergeOutput", "Could not merge all <%s> files", outputFile.Data());
               merged = kFALSE;
               break;
            } else {
               ::Info("MergeOutputs", "\n#####   Merged %d output files to <%s>", fm->GetMergeList()->GetSize(), outputChunk.Data());
               gSystem->Unlink(previousChunk);
            }
            if (nextfile == listoffiles->Last()) break;
            countChunk++;
            countZero = nmaxmerge;
            previousChunk = outputChunk;
         }
      }
      delete listoffiles;
      delete fm;
      return merged;
   }
   // Merging stage different than 0.
   // Move to the begining of the requested chunk.
   fm = new TFileMerger(kFALSE);
   fm->SetFastMethod(kTRUE);
   while ((nextfile=next())) fm->AddFile(nextfile->GetName());
   delete listoffiles;
   if (!fm->GetMergeList() || !fm->GetMergeList()->GetSize()) {
      // Nothing found - skip this output
      ::Warning("MergeOutput", "No <%s> files found.", inputFile.Data());
      delete fm;
      return kFALSE;
   }
   fm->OutputFile(outputFile);
   // Merge the outputs
   if (!fm->Merge()) {
      ::Error("MergeOutput", "Could not merge all <%s> files", outputFile.Data());
      delete fm;
      return kFALSE;
   } else {
      ::Info("MergeOutput", "\n#####   Merged %d output files to <%s>", fm->GetMergeList()->GetSize(), outputFile.Data());
   }
   delete fm;
   return kTRUE;
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
   if (fMergeViaJDL) {
      if (!TestBit(AliAnalysisGrid::kMerge)) {
         Info("MergeOutputs", "### Re-run with <MergeViaJDL> option in terminate mode of the plugin to submit merging jobs ###");
         return kFALSE; 
      }     
      if (fProductionMode) {
         Info("MergeOutputs", "### Merging will be submitted by LPM manager... ###");
         return kFALSE;
      }
      Info("MergeOutputs", "Submitting merging JDL");
      if (!SubmitMerging()) return kFALSE;
      Info("MergeOutputs", "### Re-run with <MergeViaJDL> off to collect results after merging jobs are done ###");
      Info("MergeOutputs", "### The Terminate() method is executed by the merging jobs");
      return kFALSE;
   }   
   // Get the output path
   if (!fGridOutputDir.Contains("/")) fGridOutputDir = Form("%s/%s/%s", gGrid->GetHomeDirectory(), fGridWorkingDir.Data(), fGridOutputDir.Data());
   if (!DirectoryExists(fGridOutputDir)) {
      Error("MergeOutputs", "Grid output directory %s not found. Terminate() will NOT be executed", fGridOutputDir.Data());
      return kFALSE;
   }
   if (!fOutputFiles.Length()) {
      Error("MergeOutputs", "No output file names defined. Are you running the right AliAnalysisAlien configuration ?");
      return kFALSE;
   }
   // Check if fast read option was requested
   Info("MergeOutputs", "Started local merging of output files from: alien://%s \
        \n======= overwrite mode = %d", fGridOutputDir.Data(), (Int_t)fOverwriteMode);
   if (fFastReadOption) {
      Warning("MergeOutputs", "You requested FastRead option. Using xrootd flags to reduce timeouts. This may skip some files that could be accessed ! \
             \n+++ NOTE: To disable this option, use: plugin->SetFastReadOption(kFALSE)");
      gEnv->SetValue("XNet.ConnectTimeout",50);
      gEnv->SetValue("XNet.RequestTimeout",50);
      gEnv->SetValue("XNet.MaxRedirectCount",2);
      gEnv->SetValue("XNet.ReconnectTimeout",50);
      gEnv->SetValue("XNet.FirstConnectMaxCnt",1);
   }   
   // Make sure we change the temporary directory
   gSystem->Setenv("TMPDIR", gSystem->pwd());
   // Set temporary compilation directory to current one
   gSystem->SetBuildDir(gSystem->pwd(), kTRUE);   
   TObjArray *list = fOutputFiles.Tokenize(",");
   TIter next(list);
   TObjString *str;
   TString outputFile;
   Bool_t merged = kTRUE;
   while((str=(TObjString*)next())) {
      outputFile = str->GetString();
      Int_t index = outputFile.Index("@");
      if (index > 0) outputFile.Remove(index);
      TString outputChunk = outputFile;
      outputChunk.ReplaceAll(".root", "_*.root");
      // Skip already merged outputs
      if (!gSystem->AccessPathName(outputFile)) {
         if (fOverwriteMode) {
            Info("MergeOutputs", "Overwrite mode. Existing file %s was deleted.", outputFile.Data());
            gSystem->Unlink(outputFile);
            if (!gSystem->Exec(Form("ls %s 2>/dev/null", outputChunk.Data()))) {
               Info("MergeOutput", "Overwrite mode: partial merged files %s will removed",
                     outputChunk.Data());
               gSystem->Exec(Form("rm -f %s", outputChunk.Data()));
            }
         } else {   
            Info("MergeOutputs", "Output file <%s> found. Not merging again.", outputFile.Data());
            continue;
         }   
      } else {
         if (!gSystem->Exec(Form("ls %s 2>/dev/null", outputChunk.Data()))) {
            Info("MergeOutput", "Overwrite mode: partial merged files %s will removed",
                  outputChunk.Data());
            gSystem->Exec(Form("rm -f %s", outputChunk.Data()));
         }   
      }
      if (fMergeExcludes.Length() &&
          fMergeExcludes.Contains(outputFile.Data())) continue;
      // Perform a 'find' command in the output directory, looking for registered outputs    
      merged = MergeOutput(outputFile, fGridOutputDir, fMaxMergeFiles);
      if (!merged) {
         Error("MergeOutputs", "Terminate() will  NOT be executed");
         return kFALSE;
      }
      TFile *fileOpened = (TFile*)gROOT->GetListOfFiles()->FindObject(outputFile);
      if (fileOpened) fileOpened->Close();
   } 
   return kTRUE;
}   

//______________________________________________________________________________
void AliAnalysisAlien::SetDefaultOutputs(Bool_t flag)
{
// Use the output files connected to output containers from the analysis manager
// rather than the files defined by SetOutputFiles
   if (flag && !TObject::TestBit(AliAnalysisGrid::kDefaultOutputs))
      Info("SetDefaultOutputs", "Plugin will use the output files taken from analysis manager");
   TObject::SetBit(AliAnalysisGrid::kDefaultOutputs, flag);
}
      
//______________________________________________________________________________
void AliAnalysisAlien::SetOutputFiles(const char *list)
{
// Manually set the output files list.
// Removes duplicates. Not allowed if default outputs are not disabled.
   if (TObject::TestBit(AliAnalysisGrid::kDefaultOutputs)) {
      Fatal("SetOutputFiles", "You have to explicitly call SetDefaultOutputs(kFALSE) to manually set output files.");
      return;
   }
   Info("SetOutputFiles", "Output file list is set manually - you are on your own.");
   fOutputFiles = "";
   TString slist = list;
   if (slist.Contains("@")) Warning("SetOutputFiles","The plugin does not allow explicit SE's. Please use: SetNumberOfReplicas() instead.");
   TObjArray *arr = slist.Tokenize(" "); 
   TObjString *os;
   TIter next(arr);
   TString sout;
   while ((os=(TObjString*)next())) {
      sout = os->GetString();
      if (sout.Index("@")>0) sout.Remove(sout.Index("@"));
      if (fOutputFiles.Contains(sout)) continue;
      if (!fOutputFiles.IsNull()) fOutputFiles += ",";
      fOutputFiles += sout;
   }
   delete arr;   
}

//______________________________________________________________________________
void AliAnalysisAlien::SetOutputArchive(const char *list)
{
// Manually set the output archive list. Free text - you are on your own...
// Not allowed if default outputs are not disabled.
   if (TObject::TestBit(AliAnalysisGrid::kDefaultOutputs)) {
      Fatal("SetOutputArchive", "You have to explicitly call SetDefaultOutputs(kFALSE) to manually set the output archives.");
      return;
   }
   Info("SetOutputArchive", "Output archive is set manually - you are on your own.");
   fOutputArchive = list;
}

//______________________________________________________________________________
void AliAnalysisAlien::SetPreferedSE(const char */*se*/)
{
// Setting a prefered output SE is not allowed anymore.
   Warning("SetPreferedSE", "Setting a preferential SE is not allowed anymore via the plugin. Use SetNumberOfReplicas() and SetDefaultOutputs()");
}

//______________________________________________________________________________
void AliAnalysisAlien::SetProofParameter(const char *pname, const char *value)
{
// Set some PROOF special parameter.
   TPair *pair = dynamic_cast<TPair*>(fProofParam.FindObject(pname));
   if (pair) {
      TObject *old = pair->Key();
      TObject *val = pair->Value();
      fProofParam.Remove(old);
      delete old;
      delete val;
   }
   fProofParam.Add(new TObjString(pname), new TObjString(value));
}

//______________________________________________________________________________
const char *AliAnalysisAlien::GetProofParameter(const char *pname) const
{
// Returns a special PROOF parameter.
   TPair *pair = dynamic_cast<TPair*>(fProofParam.FindObject(pname));
   if (!pair) return 0;
   return pair->Value()->GetName();
}      

//______________________________________________________________________________
Bool_t AliAnalysisAlien::StartAnalysis(Long64_t /*nentries*/, Long64_t /*firstEntry*/)
{
// Start remote grid analysis.
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   Bool_t testMode = TestBit(AliAnalysisGrid::kTest);
   if (!mgr || !mgr->IsInitialized()) {
      Error("StartAnalysis", "You need an initialized analysis manager for this");
      return kFALSE;
   }
   // Are we in PROOF mode ?
   if (mgr->IsProofMode()) {
      if (testMode) Info("StartAnalysis", "##### Starting PROOF analysis with Proof Lite via the plugin #####");
      else Info("StartAnalysis", "##### Starting PROOF analysis on cluster <%s> via the plugin #####", fProofCluster.Data());
      if (fProofCluster.IsNull()) {
         Error("StartAnalysis", "You need to specify the proof cluster name via SetProofCluster");
         return kFALSE;
      }   
      if (fProofDataSet.IsNull() && !testMode) {
         Error("StartAnalysis", "You need to specify a dataset using SetProofDataSet()");
         return kFALSE;
      }   
      // Set the needed environment
      gEnv->SetValue("XSec.GSI.DelegProxy","2");
      // Do we need to reset PROOF ? The success of the Reset operation cannot be checked
      if (fProofReset && !testMode) {
         if (fProofReset==1) {
            Info("StartAnalysis", "Sending soft reset signal to proof cluster %s", fProofCluster.Data());
            gROOT->ProcessLine(Form("TProof::Reset(\"%s\", kFALSE);", fProofCluster.Data()));
         } else {         
            Info("StartAnalysis", "Sending hard reset signal to proof cluster %s", fProofCluster.Data());
            gROOT->ProcessLine(Form("TProof::Reset(\"%s\", kTRUE);", fProofCluster.Data()));
         }
         Info("StartAnalysis", "Stopping the analysis. Please use SetProofReset(0) to resume.");
         return kFALSE;
      }
      
      if (!testMode) {
        // Check if there is an old active session
        Long_t nsessions = gROOT->ProcessLine(Form("TProof::Mgr(\"%s\")->QuerySessions(\"\")->GetEntries();", fProofCluster.Data()));
        if (nsessions) {
          Error("StartAnalysis","You have to reset your old session first\n");
          return kFALSE;
        }
      }
      // Do we need to change the ROOT version ? The success of this cannot be checked.
      if (!fRootVersionForProof.IsNull() && !testMode) {
         gROOT->ProcessLine(Form("TProof::Mgr(\"%s\")->SetROOTVersion(\"%s\");", 
                            fProofCluster.Data(), fRootVersionForProof.Data()));
      }
      // Connect to PROOF and check the status
      Long_t proof = 0;
      TString sworkers;
      if (fNproofWorkersPerSlave) sworkers = Form("workers=%dx", fNproofWorkersPerSlave);
      else if (fNproofWorkers) sworkers = Form("workers=%d", fNproofWorkers);
      if (!testMode) {
         if (!sworkers.IsNull()) 
            proof = gROOT->ProcessLine(Form("TProof::Open(\"%s\", \"%s\");", fProofCluster.Data(), sworkers.Data()));
         else   
            proof = gROOT->ProcessLine(Form("TProof::Open(\"%s\");", fProofCluster.Data()));
      } else {
         proof = gROOT->ProcessLine("TProof::Open(\"\");");
         if (!proof) {
            Error("StartAnalysis", "Could not start PROOF in test mode");
            return kFALSE;
         }   
      }
      if (!proof) {
         Error("StartAnalysis", "Could not connect to PROOF cluster <%s>", fProofCluster.Data());
         return kFALSE;
      }   
      if (fNproofWorkersPerSlave*fNproofWorkers > 0)
         gROOT->ProcessLine(Form("gProof->SetParallel(%d);", fNproofWorkers));
      // Set proof special parameters if any
      TIter nextpp(&fProofParam);
      TObject *proofparam;
      while ((proofparam=nextpp())) {
         TString svalue = GetProofParameter(proofparam->GetName());
         gROOT->ProcessLine(Form("gProof->SetParameter(\"%s\",%s);", proofparam->GetName(), svalue.Data()));
      }   
      // Is dataset existing ?
      if (!testMode) {
         TString dataset = fProofDataSet;
         Int_t index = dataset.Index("#");
         if (index>=0) dataset.Remove(index);
//         if (!gROOT->ProcessLine(Form("gProof->ExistsDataSet(\"%s\");",fProofDataSet.Data()))) {
//            Error("StartAnalysis", "Dataset %s not existing", fProofDataSet.Data());
//            return kFALSE;
//         }
//         Info("StartAnalysis", "Dataset %s found", dataset.Data());
      }
      // Is ClearPackages() needed ?
      if (TestSpecialBit(kClearPackages)) {
         Info("StartAnalysis", "ClearPackages signal sent to PROOF. Use SetClearPackages(kFALSE) to reset this.");
         gROOT->ProcessLine("gProof->ClearPackages();");
      }
      // Is a given aliroot mode requested ?
      TList optionsList;
      TString parLibs;
      if (!fAliRootMode.IsNull()) {
         TString alirootMode = fAliRootMode;
         if (alirootMode == "default") alirootMode = "";
         Info("StartAnalysis", "You are requesting AliRoot mode: %s", fAliRootMode.Data());
         optionsList.SetOwner();
         optionsList.Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));
         // Check the additional libs to be loaded
         TString extraLibs;
         Bool_t parMode = kFALSE;
         if (!alirootMode.IsNull()) extraLibs = "ANALYSIS:OADB:ANALYSISalice";
         // Parse the extra libs for .so
         if (fAdditionalLibs.Length()) {
            TObjArray *list = fAdditionalLibs.Tokenize(" ");
            TIter next(list);
            TObjString *str;
            while((str=(TObjString*)next())) {
               if (str->GetString().Contains(".so")) {
                  if (parMode) {
                     Warning("StartAnalysis", "Plugin does not support loading libs after par files in PROOF mode. Library %s and following will not load on workers", str->GetName());
                     break;
                  }   
                  TString stmp = str->GetName();
                  if (stmp.BeginsWith("lib")) stmp.Remove(0,3);
                  stmp.ReplaceAll(".so","");
                  if (!extraLibs.IsNull()) extraLibs += ":";
                  extraLibs += stmp;
                  continue;
               }
               if (str->GetString().Contains(".par")) {
                  // The first par file found in the list will not allow any further .so
                  parMode = kTRUE;
                  if (!parLibs.IsNull()) parLibs += ":";
                  parLibs += str->GetName();
                  continue;
               }   
            }
            if (list) delete list;            
         }
        if (!extraLibs.IsNull()) {
          Info("StartAnalysis", "Adding extra libs: %s",extraLibs.Data());
          optionsList.Add(new TNamed("ALIROOT_EXTRA_LIBS",extraLibs.Data()));
        }
         // Check extra includes
         if (!fIncludePath.IsNull()) {
            TString includePath = fIncludePath;
            includePath.ReplaceAll(" ",":");
            includePath.ReplaceAll("$ALICE_ROOT/","");
            includePath.ReplaceAll("${ALICE_ROOT}/","");
            includePath.ReplaceAll("-I","");
            includePath.Remove(TString::kTrailing, ':');
            Info("StartAnalysis", "Adding extra includes: %s",includePath.Data()); 
            optionsList.Add(new TNamed("ALIROOT_EXTRA_INCLUDES",includePath.Data()));
         }
         // Check if connection to grid is requested
         if (TestSpecialBit(kProofConnectGrid)) 
            optionsList.Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));
         // Enable AliRoot par
         if (testMode) {
         // Enable proof lite package
            TString alirootLite = gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AliRootProofLite.par");
            for (Int_t i=0; i<optionsList.GetSize(); i++) {
               TNamed *obj = (TNamed*)optionsList.At(i);
               printf("%s  %s\n", obj->GetName(), obj->GetTitle());
            }   
            if (!gROOT->ProcessLine(Form("gProof->UploadPackage(\"%s\");",alirootLite.Data()))
              && !gROOT->ProcessLine(Form("gProof->EnablePackage(\"%s\", (TList*)%p);",alirootLite.Data(),&optionsList))) {
                  Info("StartAnalysis", "AliRootProofLite enabled");
            } else {                      
               Error("StartAnalysis", "There was an error trying to enable package AliRootProofLite.par");
               return kFALSE;
            }   
         } else {
           if ( ! fAliROOTVersion.IsNull() ) {
             if (gROOT->ProcessLine(Form("gProof->EnablePackage(\"VO_ALICE@AliRoot::%s\", (TList*)%p, kTRUE);", 
                                         fAliROOTVersion.Data(), &optionsList))) {
                Error("StartAnalysis", "There was an error trying to enable package VO_ALICE@AliRoot::%s", fAliROOTVersion.Data());
                return kFALSE;
             }
           }
         }
         // Enable first par files from fAdditionalLibs
         if (!parLibs.IsNull()) {
            TObjArray *list = parLibs.Tokenize(":");
            TIter next(list);
            TObjString *package;
            while((package=(TObjString*)next())) {
               TString spkg = package->GetName();
               spkg.ReplaceAll(".par", "");
               gSystem->Exec(TString::Format("rm -rf %s", spkg.Data()));
               if (!gROOT->ProcessLine(Form("gProof->UploadPackage(\"%s\");", package->GetName()))) {
                  TString enablePackage = (testMode)?Form("gProof->EnablePackage(\"%s\",kFALSE);", package->GetName()):Form("gProof->EnablePackage(\"%s\",kTRUE);", package->GetName());
                  if (gROOT->ProcessLine(enablePackage)) {
                     Error("StartAnalysis", "There was an error trying to enable package %s", package->GetName());
                     return kFALSE;
                  }
               } else {
                  Error("StartAnalysis", "There was an error trying to upload package %s", package->GetName());
                  return kFALSE;
               }
            }
            if (list) delete list; 
         }
      } else {
         if (fAdditionalLibs.Contains(".so") && !testMode) {
            Error("StartAnalysis", "You request additional libs to be loaded but did not enabled any AliRoot mode. Please refer to: \
                   \n http://aaf.cern.ch/node/83 and use a parameter for SetAliRootMode()");
            return kFALSE;       
         }
      }
      // Enable par files if requested
      if (fPackages && fPackages->GetEntries()) {
         TIter next(fPackages);
         TObject *package;
         while ((package=next())) {
            // Skip packages already enabled
            if (parLibs.Contains(package->GetName())) continue;
            TString spkg = package->GetName();
            spkg.ReplaceAll(".par", "");
            gSystem->Exec(TString::Format("rm -rf %s", spkg.Data()));
            if (!gROOT->ProcessLine(Form("gProof->UploadPackage(\"%s\");", package->GetName()))) {
               if (gROOT->ProcessLine(Form("gProof->EnablePackage(\"%s\",kTRUE);", package->GetName()))) {
                  Error("StartAnalysis", "There was an error trying to enable package %s", package->GetName());
                  return kFALSE;
               }
            } else {
               Error("StartAnalysis", "There was an error trying to upload package %s", package->GetName());
               return kFALSE;
            }
         }
      }
      // Do we need to load analysis source files ?
      // NOTE: don't load on client since this is anyway done by the user to attach his task.
      if (fAnalysisSource.Length()) {
         TObjArray *list = fAnalysisSource.Tokenize(" ");
         TIter next(list);
         TObjString *str;
         while((str=(TObjString*)next())) {
            gROOT->ProcessLine(Form("gProof->Load(\"%s+g\", kTRUE);", str->GetName()));
         }
         if (list) delete list;
      }
      if (testMode) {
      // Register dataset to proof lite.
         if (fFileForTestMode.IsNull()) {
            Error("GetChainForTestMode", "For proof test mode please use SetFileForTestMode() pointing to a file that contains data file locations.");
            return kFALSE;
         }
         if (gSystem->AccessPathName(fFileForTestMode)) {
            Error("GetChainForTestMode", "File not found: %s", fFileForTestMode.Data());
            return kFALSE;
         }   
         TFileCollection *coll = new TFileCollection();
         coll->AddFromFile(fFileForTestMode);
         gROOT->ProcessLine(Form("gProof->RegisterDataSet(\"test_collection\", (TFileCollection*)%p, \"OV\");", coll));
         gROOT->ProcessLine("gProof->ShowDataSets()");
      }
      return kTRUE;
   }
   
   // Check if output files have to be taken from the analysis manager
   if (TestBit(AliAnalysisGrid::kDefaultOutputs)) {
      // Add output files and AOD files
      fOutputFiles = GetListOfFiles("outaod");
      // Add extra files registered to the analysis manager
      TString extra = GetListOfFiles("ext");
      if (!extra.IsNull()) {
         extra.ReplaceAll(".root", "*.root");
         if (!fOutputFiles.IsNull()) fOutputFiles += ",";
         fOutputFiles += extra;
      }
      // Compose the output archive.
      fOutputArchive = "log_archive.zip:std*@disk=1 ";
      fOutputArchive += Form("root_archive.zip:%s,*.stat@disk=%d",fOutputFiles.Data(),fNreplicas);
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
      if (fMergeViaJDL) CheckInputData();
      return kTRUE;
   } else {
      Info("StartAnalysis","\n##### FULL ANALYSIS MODE ##### Producing needed files and submitting your analysis job...");   
   }   
      
   Print();   
   if (!Connect()) {
      Error("StartAnalysis", "Cannot start grid analysis without grid connection");
      return kFALSE;
   }
   if (IsCheckCopy() && gGrid) CheckFileCopy(gGrid->GetHomeDirectory());
   if (!CheckInputData()) {
      Error("StartAnalysis", "There was an error in preprocessing your requested input data");
      return kFALSE;
   }   
   if (!CreateDataset(fDataPattern)) {
      TString serror;
      if (!fRunNumbers.Length() && !fRunRange[0]) serror = Form("path to data directory: <%s>", fGridDataDir.Data());
      if (fRunNumbers.Length()) serror = "run numbers";
      if (fRunRange[0]) serror = Form("run range [%d, %d]", fRunRange[0], fRunRange[1]);
      serror += Form("\n   or data pattern <%s>", fDataPattern.Data());
      Error("StartAnalysis", "No data to process. Please fix %s in your plugin configuration.", serror.Data());
      return kFALSE;
   }   
   WriteAnalysisFile();
   WriteAnalysisMacro();
   WriteExecutable();
   WriteValidationScript();
   if (fMergeViaJDL) {
      WriteMergingMacro();
      WriteMergeExecutable();
      WriteValidationScript(kTRUE);
   }   
   if (!CreateJDL()) return kFALSE;
   if (TestBit(AliAnalysisGrid::kOffline)) return kFALSE;
   if (testMode) {
      // Locally testing the analysis
      Info("StartAnalysis", "\n_______________________________________________________________________ \
      \n   Running analysis script in a daughter shell as on a worker node \
      \n_______________________________________________________________________");
      TObjArray *list = fOutputFiles.Tokenize(",");
      TIter next(list);
      TObjString *str;
      TString outputFile;
      while((str=(TObjString*)next())) {
         outputFile = str->GetString();
         Int_t index = outputFile.Index("@");
         if (index > 0) outputFile.Remove(index);         
         if (!gSystem->AccessPathName(outputFile)) gSystem->Exec(Form("rm %s", outputFile.Data()));
      }
      delete list;
      gSystem->Exec(Form("bash %s 2>stderr", fExecutable.Data()));
      gSystem->Exec(Form("bash %s",fValidationScript.Data()));
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
      printf("*************************** %s\n",Form("submit %s", fJDLName.Data()));
      if (res) {
         const char *cjobId = res->GetKey(0,"jobId");
         if (!cjobId) {
            gGrid->Stdout();
            gGrid->Stderr();
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
      } else {
         Error("StartAnalysis", "No grid result after submission !!! Bailing out...");
         return kFALSE;      
      }   
   } else {
      // Submit for a range of enumeration of runs.
      if (!Submit()) return kFALSE;
   }   
         
   Info("StartAnalysis", "\n#### STARTING AN ALIEN SHELL FOR YOU. EXIT WHEN YOUR JOB %s HAS FINISHED. #### \
   \n You may exit at any time and terminate the job later using the option <terminate> \
   \n ##################################################################################", jobID.Data());
   gSystem->Exec("aliensh");
   return kTRUE;
}

//______________________________________________________________________________
const char *AliAnalysisAlien::GetListOfFiles(const char *type)
{
// Get a comma-separated list of output files of the requested type.
// Type can be (case unsensitive):
//    aod - list of aod files (std, extensions and filters)
//    out - list of output files connected to containers (but not aod's or extras)
//    ext - list of extra files registered to the manager
//    ter - list of files produced in terminate
   static TString files;
   files = "";
   TString stype = type;
   stype.ToLower();
   TString aodfiles, extra;
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("GetListOfFiles", "Cannot call this without analysis manager");
      return files.Data();
   }
   if (mgr->GetOutputEventHandler()) {
      aodfiles = mgr->GetOutputEventHandler()->GetOutputFileName();
      TString extraaod = mgr->GetOutputEventHandler()->GetExtraOutputs();
      if (!extraaod.IsNull()) {
         aodfiles += ",";
         aodfiles += extraaod;
      }
   }
   if (stype.Contains("aod")) {
      files = aodfiles;
      if (stype == "aod") return files.Data();
   }  
   // Add output files that are not in the list of AOD files 
   TString outputfiles = "";
   TIter next(mgr->GetOutputs());
   AliAnalysisDataContainer *output;
   const char *filename = 0;
   while ((output=(AliAnalysisDataContainer*)next())) {
      filename = output->GetFileName();
      if (!(strcmp(filename, "default"))) continue;
      if (outputfiles.Contains(filename)) continue;
      if (aodfiles.Contains(filename))    continue;
      if (!outputfiles.IsNull()) outputfiles += ",";
      outputfiles += filename;
   }
   if (stype.Contains("out")) {
      if (!files.IsNull()) files += ",";
      files += outputfiles;
      if (stype == "out") return files.Data();
   }   
   // Add extra files registered to the analysis manager
   TString sextra;
   extra = mgr->GetExtraFiles();
   if (!extra.IsNull()) {
      extra.Strip();
      extra.ReplaceAll(" ", ",");
      TObjArray *fextra = extra.Tokenize(",");
      TIter nextx(fextra);
      TObject *obj;
      while ((obj=nextx())) {
         if (aodfiles.Contains(obj->GetName())) continue;
         if (outputfiles.Contains(obj->GetName())) continue;
         if (sextra.Contains(obj->GetName())) continue;
         if (!sextra.IsNull()) sextra += ",";
         sextra += obj->GetName();
      }
      delete fextra;
      if (stype.Contains("ext")) {
         if (!files.IsNull()) files += ",";
         files += sextra;
      }
   }   
   if (stype == "ext") return files.Data();
   TString termfiles;
   if (!fTerminateFiles.IsNull()) {
      fTerminateFiles.Strip();
      fTerminateFiles.ReplaceAll(" ",",");
      TObjArray *fextra = fTerminateFiles.Tokenize(",");
      TIter nextx(fextra);
      TObject *obj;
      while ((obj=nextx())) {
         if (aodfiles.Contains(obj->GetName())) continue;
         if (outputfiles.Contains(obj->GetName())) continue;
         if (termfiles.Contains(obj->GetName())) continue;
         if (sextra.Contains(obj->GetName())) continue;
         if (!termfiles.IsNull()) termfiles += ",";
         termfiles += obj->GetName();
      }
      delete fextra;
   }   
   if (stype.Contains("ter")) {
      if (!files.IsNull() && !termfiles.IsNull()) {
         files += ",";
         files += termfiles;
      }   
   }   
   return files.Data();
}   

//______________________________________________________________________________
Bool_t AliAnalysisAlien::Submit()
{
// Submit all master jobs.
   Int_t nmasterjobs = fInputFiles->GetEntries();
   Long_t tshoot = gSystem->Now();
   if (!fNsubmitted && !SubmitNext()) return kFALSE;
   while (fNsubmitted < nmasterjobs) {
      Long_t now = gSystem->Now();
      if ((now-tshoot)>30000) {
         tshoot = now;
         if (!SubmitNext()) return kFALSE;
      }   
   }
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::SubmitMerging()
{
// Submit all merging jobs.
   if (!fGridOutputDir.Contains("/")) fGridOutputDir = Form("%s/%s/%s", gGrid->GetHomeDirectory(), fGridWorkingDir.Data(), fGridOutputDir.Data());
   gGrid->Cd(fGridOutputDir);
   TString mergeJDLName = fExecutable;
   mergeJDLName.ReplaceAll(".sh", "_merge.jdl");
   if (!fInputFiles) {
      Error("SubmitMerging", "You have to use explicit run numbers or run range to merge via JDL!");
      return kFALSE;
   }   
   Int_t ntosubmit = fInputFiles->GetEntries();
   for (Int_t i=0; i<ntosubmit; i++) {
      TString runOutDir = gSystem->BaseName(fInputFiles->At(i)->GetName());
      runOutDir.ReplaceAll(".xml", "");
      if (fOutputToRunNo) {
         // The output directory is the run number
         printf("### Submitting merging job for run <%s>\n", runOutDir.Data());
         runOutDir = Form("%s/%s", fGridOutputDir.Data(), runOutDir.Data());
      } else {
         if (!fRunNumbers.Length() && !fRunRange[0]) {
            // The output directory is the grid outdir
            printf("### Submitting merging job for the full output directory %s.\n", fGridOutputDir.Data());
            runOutDir = fGridOutputDir;
         } else {
            // The output directory is the master number in 3 digits format
            printf("### Submitting merging job for master <%03d>\n", i);
            runOutDir = Form("%s/%03d",fGridOutputDir.Data(), i);
         }   
      }
      // Check now the number of merging stages.
      TObjArray *list = fOutputFiles.Tokenize(",");
      TIter next(list);
      TObjString *str;
      TString outputFile;
      while((str=(TObjString*)next())) {
         outputFile = str->GetString();
         Int_t index = outputFile.Index("@");
         if (index > 0) outputFile.Remove(index);
         if (!fMergeExcludes.Contains(outputFile)) break;
      }
      delete list;
      Bool_t done = CheckMergedFiles(outputFile, runOutDir, fMaxMergeFiles, mergeJDLName);
      if (!done && (i==ntosubmit-1)) return kFALSE;
      if (!fRunNumbers.Length() && !fRunRange[0]) break;
   }
   if (!ntosubmit) return kTRUE;
   Info("StartAnalysis", "\n #### STARTING AN ALIEN SHELL FOR YOU. You can exit any time or inspect your jobs in a different shell.##########\
                          \n Make sure your jobs are in a final state (you can resubmit failed ones via 'masterjob <id> resubmit ERROR_ALL')\
                          \n Rerun in 'terminate' mode to submit all merging stages, each AFTER the previous one completed. The final merged \
                          \n output will be written to your alien output directory, while separate stages in <Stage_n>. \
                          \n ################################################################################################################");
   gSystem->Exec("aliensh");
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::SubmitNext()
{
// Submit next bunch of master jobs if the queue is free. The first master job is
// submitted right away, while the next will not be unless the previous was split.
// The plugin will not submit new master jobs if there are more that 500 jobs in
// waiting phase.
   static Bool_t iscalled = kFALSE;
   static Int_t firstmaster = 0;
   static Int_t lastmaster = 0;
   static Int_t npermaster  = 0;
   if (iscalled) return kTRUE;
   iscalled = kTRUE;
   Int_t nrunning=0, nwaiting=0, nerror=0, ndone=0;
   Int_t ntosubmit = 0;
   TGridResult *res;
   TString jobID = "";
   Int_t nmasterjobs = fInputFiles->GetEntries();
   if (!fNsubmitted) {
      ntosubmit = 1;
      if (!IsUseSubmitPolicy()) {
         if (nmasterjobs>5)
            Info("SubmitNext","### Warning submit policy not used ! Submitting too many jobs at a time may be prohibitted. \
                \n### You can use SetUseSubmitPolicy() to enable if you have problems.");
         ntosubmit = nmasterjobs;
      }   
   } else {
      TString status = GetJobStatus(firstmaster, lastmaster, nrunning, nwaiting, nerror, ndone);
      printf("=== master %d: %s\n", lastmaster, status.Data());
      // If last master not split, just return
      if (status != "SPLIT") {iscalled = kFALSE; return kTRUE;}
      // No more than 100 waiting jobs
      if (nwaiting>500) {iscalled = kFALSE; return kTRUE;}
      npermaster = (nrunning+nwaiting+nerror+ndone)/fNsubmitted;      
      if (npermaster) ntosubmit = (500-nwaiting)/npermaster;
      if (!ntosubmit) ntosubmit = 1;
      printf("=== WAITING(%d) RUNNING(%d) DONE(%d) OTHER(%d) NperMaster=%d => to submit %d jobs\n", 
             nwaiting, nrunning, ndone, nerror, npermaster, ntosubmit);
   }
   for (Int_t i=0; i<ntosubmit; i++) {
      // Submit for a range of enumeration of runs.
      if (fNsubmitted>=nmasterjobs) {iscalled = kFALSE; return kTRUE;}
      TString query;
      TString runOutDir = gSystem->BaseName(fInputFiles->At(fNsubmitted)->GetName());
      runOutDir.ReplaceAll(".xml", "");
      if (fOutputToRunNo)
         query = Form("submit %s %s %s", fJDLName.Data(), fInputFiles->At(fNsubmitted)->GetName(), runOutDir.Data());
      else
         query = Form("submit %s %s %03d", fJDLName.Data(), fInputFiles->At(fNsubmitted)->GetName(), fNsubmitted);
      printf("********* %s\n",query.Data());
      res = gGrid->Command(query);
      if (res) {
         TString cjobId1 = res->GetKey(0,"jobId");
         if (!cjobId1.Length()) {
            iscalled = kFALSE;
            gGrid->Stdout();
            gGrid->Stderr();
            Error("StartAnalysis", "Your JDL %s could not be submitted. The message was:", fJDLName.Data());
            return kFALSE;
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
      } else {
         Error("StartAnalysis", "No grid result after submission !!! Bailing out...");
         return kFALSE;
      }   
   }
   iscalled = kFALSE;
   return kTRUE;
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
         if (handler->InheritsFrom("AliMultiInputEventHandler")) {
            AliMultiInputEventHandler *multiIH = (AliMultiInputEventHandler*)handler;
            if (multiIH->GetFirstInputEventHandler()->InheritsFrom("AliESDInputHandler")) TObject::SetBit(AliAnalysisGrid::kUseESD);
            if (multiIH->GetFirstInputEventHandler()->InheritsFrom("AliAODInputHandler")) TObject::SetBit(AliAnalysisGrid::kUseAOD);
         } else {
            if (handler->InheritsFrom("AliESDInputHandler")) TObject::SetBit(AliAnalysisGrid::kUseESD);
            if (handler->InheritsFrom("AliAODInputHandler")) TObject::SetBit(AliAnalysisGrid::kUseAOD);
         }
      }
      TDirectory *cdir = gDirectory;
      TFile *file = TFile::Open(analysisFile, "RECREATE");
      if (file) {
         // Skip task Terminate calls for the grid job (but not in test mode, where we want to check also the terminate mode
         if (!TestBit(AliAnalysisGrid::kTest)) mgr->SetSkipTerminate(kTRUE);
         // Unless merging makes no sense
         if (IsSingleOutput()) mgr->SetSkipTerminate(kFALSE);
         mgr->Write();
         delete file;
         // Enable termination for local jobs
         mgr->SetSkipTerminate(kFALSE);
      }
      if (cdir) cdir->cd();
      Info("WriteAnalysisFile", "\n#####   Analysis manager: %s wrote to file <%s>\n", mgr->GetName(),analysisFile.Data());
   }   
   Bool_t copy = kTRUE;
   if (fProductionMode || TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      workdir += fGridWorkingDir;
      Info("WriteAnalysisFile", "\n#####   Copying file <%s> containing your initialized analysis manager to your alien workspace", analysisFile.Data());
      if (FileExists(analysisFile)) gGrid->Rm(analysisFile);
      TFile::Cp(Form("file:%s",analysisFile.Data()), Form("alien://%s/%s", workdir.Data(),analysisFile.Data()));
      if (!copyLocal2Alien("WriteAnalysisFile",analysisFile.Data(), 
          Form("%s/%s", workdir.Data(),analysisFile.Data()))) Fatal("","Terminating");
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
      Bool_t hasSTEERBase = kFALSE;
      Bool_t hasESD = kFALSE;
      Bool_t hasAOD = kFALSE;
      Bool_t hasANALYSIS = kFALSE;
      Bool_t hasOADB = kFALSE;
      Bool_t hasANALYSISalice = kFALSE;
      Bool_t hasCORRFW = kFALSE;
      TString func = fAnalysisMacro;
      TString type = "ESD";
      TString comment = "// Analysis using ";
      if (IsUseMCchain()) {
         type = "MC";
         comment += "MC";
      } else {   
         if (TObject::TestBit(AliAnalysisGrid::kUseESD)) comment += "ESD";
         if (TObject::TestBit(AliAnalysisGrid::kUseAOD)) {
            type = "AOD";
            comment += "AOD";
         }   
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
      // Change temp directory to current one
      out << "// Set temporary merging directory to current one" << endl;
      out << "   gSystem->Setenv(\"TMPDIR\", gSystem->pwd());" << endl << endl;   
      out << "// Set temporary compilation directory to current one" << endl;
      out << "   gSystem->SetBuildDir(gSystem->pwd(), kTRUE);" << endl << endl;   
      // Reset existing include path
      out << "// Reset existing include path and add current directory first in the search" << endl;
      out << "   gSystem->SetIncludePath(\"-I.\");" << endl;
      if (!fExecutableCommand.Contains("aliroot")) {
         out << "// load base root libraries" << endl;
         out << "   gSystem->Load(\"libTree\");" << endl;
         out << "   gSystem->Load(\"libGeom\");" << endl;
         out << "   gSystem->Load(\"libVMC\");" << endl;
         out << "   gSystem->Load(\"libPhysics\");" << endl << endl;
         out << "   gSystem->Load(\"libMinuit\");" << endl << endl;
      }   
      if (fAdditionalRootLibs.Length()) {
         // in principle libtree /lib geom libvmc etc. can go into this list, too
         out << "// Add aditional libraries" << endl;
         TObjArray *list = fAdditionalRootLibs.Tokenize(" ");
         TIter next(list);
         TObjString *str;
         while((str=(TObjString*)next())) {
            if (str->GetString().Contains(".so"))
            out << "   gSystem->Load(\"" << str->GetString().Data() << "\");" << endl;
         }
         if (list) delete list;
      }
      out << "// Load analysis framework libraries" << endl;
      TString setupPar = "AliAnalysisAlien::SetupPar";
      if (!fPackages) {
         if (!fExecutableCommand.Contains("aliroot")) {         
            out << "   gSystem->Load(\"libSTEERBase\");" << endl;
            out << "   gSystem->Load(\"libESD\");" << endl;
            out << "   gSystem->Load(\"libAOD\");" << endl;
         }   
         out << "   gSystem->Load(\"libANALYSIS\");" << endl;
         out << "   gSystem->Load(\"libOADB\");" << endl;
         out << "   gSystem->Load(\"libANALYSISalice\");" << endl;
         out << "   gSystem->Load(\"libCORRFW\");" << endl << endl;
      } else {
         TIter next(fPackages);
         TObject *obj;
         TString pkgname;
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
            if (pkgname == "OADB" ||
                pkgname == "OADB.par")      hasOADB = kTRUE;
            if (pkgname == "ANALYSISalice" ||
                pkgname == "ANALYSISalice.par") hasANALYSISalice = kTRUE;
            if (pkgname == "CORRFW" ||
                pkgname == "CORRFW.par")    hasCORRFW = kTRUE;
         }
         if (hasANALYSISalice) setupPar = "SetupPar";   
         if (!hasSTEERBase) out << "   gSystem->Load(\"libSTEERBase\");" << endl;
         else out << "   if (!" << setupPar << "(\"STEERBase\")) return;" << endl;
         if (!hasESD)       out << "   gSystem->Load(\"libESD\");" << endl;
         else out << "   if (!" << setupPar << "(\"ESD\")) return;" << endl;
         if (!hasAOD)       out << "   gSystem->Load(\"libAOD\");" << endl;
         else out << "   if (!" << setupPar << "(\"AOD\")) return;" << endl;
         if (!hasANALYSIS)  out << "   gSystem->Load(\"libANALYSIS\");" << endl;
         else out << "   if (!" << setupPar << "(\"ANALYSIS\")) return;" << endl;
         if (!hasOADB)  out << "   gSystem->Load(\"libOADB\");" << endl;
         else out << "   if (!" << setupPar << "(\"OADB\")) return;" << endl;
         if (!hasANALYSISalice)   out << "   gSystem->Load(\"libANALYSISalice\");" << endl;
         else out << "   if (!" << setupPar << "(\"ANALYSISalice\")) return;" << endl;
         if (!hasCORRFW)    out << "   gSystem->Load(\"libCORRFW\");" << endl << endl;
         else out << "   if (!" << setupPar << "(\"CORRFW\")) return;" << endl << endl;
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
                pkgname == "OADB" ||
                pkgname == "OADB.par" ||
                pkgname == "ANALYSISalice" ||
                pkgname == "ANALYSISalice.par" ||
                pkgname == "CORRFW" ||
                pkgname == "CORRFW.par") continue;
            out << "   if (!" << setupPar << "(\"" << obj->GetName() << "\")) return;" << endl;
         }   
      }   
      out << "// include path" << endl;
      // Get the include path from the interpreter and remove entries pointing to AliRoot
      out << "   TString intPath = gInterpreter->GetIncludePath();" << endl;
      out << "   TObjArray *listpaths = intPath.Tokenize(\" \");" << endl;
      out << "   TIter nextpath(listpaths);" << endl;
      out << "   TObjString *pname;" << endl;
      out << "   while ((pname=(TObjString*)nextpath())) {" << endl;
      out << "      TString current = pname->GetName();" << endl;
      out << "      if (current.Contains(\"AliRoot\") || current.Contains(\"ALICE_ROOT\")) continue;" << endl;
      out << "      gSystem->AddIncludePath(current);" << endl;
      out << "   }" << endl;
      out << "   if (listpaths) delete listpaths;" << endl;
      if (fIncludePath.Length()) out << "   gSystem->AddIncludePath(\"" << fIncludePath.Data() << "\");" << endl;
      out << "   gROOT->ProcessLine(\".include $ALICE_ROOT/include\");" << endl;
      out << "   printf(\"Include path: %s\\n\", gSystem->GetIncludePath());" << endl << endl;
      if (fAdditionalLibs.Length()) {
         out << "// Add aditional AliRoot libraries" << endl;
         TObjArray *list = fAdditionalLibs.Tokenize(" ");
         TIter next(list);
         TObjString *str;
         while((str=(TObjString*)next())) {
            if (str->GetString().Contains(".so"))
               out << "   gSystem->Load(\"" << str->GetString().Data() << "\");" << endl;
            if (str->GetString().Contains(".par"))
               out << "   if (!" << setupPar << "(\"" << str->GetString() << "\")) return;" << endl;
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
//      out << "   printf(\"Currently load libraries:\\n\");" << endl;
//      out << "   printf(\"%s\\n\", gSystem->GetLibraries());" << endl;
      if (fFastReadOption) {
         Warning("WriteAnalysisMacro", "!!! You requested FastRead option. Using xrootd flags to reduce timeouts in the grid jobs. This may skip some files that could be accessed !!! \
                \n+++ NOTE: To disable this option, use: plugin->SetFastReadOption(kFALSE)");
         out << "// fast xrootd reading enabled" << endl;
         out << "   printf(\"!!! You requested FastRead option. Using xrootd flags to reduce timeouts. Note that this may skip some files that could be accessed !!!\");" << endl;
         out << "   gEnv->SetValue(\"XNet.ConnectTimeout\",50);" << endl;
         out << "   gEnv->SetValue(\"XNet.RequestTimeout\",50);" << endl;
         out << "   gEnv->SetValue(\"XNet.MaxRedirectCount\",2);" << endl;
         out << "   gEnv->SetValue(\"XNet.ReconnectTimeout\",50);" << endl;
         out << "   gEnv->SetValue(\"XNet.FirstConnectMaxCnt\",1);" << endl << endl;
      } 
      if (!IsLocalTest()) {  
         out << "// connect to AliEn and make the chain" << endl;
         out << "   if (!TGrid::Connect(\"alien://\")) return;" << endl;
      }   
      out << "// read the analysis manager from file" << endl;
      TString analysisFile = fExecutable;
      analysisFile.ReplaceAll(".sh", ".root");
      out << "   AliAnalysisManager *mgr = AliAnalysisAlien::LoadAnalysisManager(\"" 
          << analysisFile << "\");" << endl;
      out << "   if (!mgr) return;" << endl;
      if (IsLocalTest()) {
         out << "   AliAnalysisAlien *plugin = new AliAnalysisAlien();" << endl;
         out << "   plugin->SetRunMode(\"test\");" << endl;
         out << "   plugin->SetFileForTestMode(\"data.txt\");" << endl;
         out << "   mgr->SetGridHandler(plugin);" << endl;
         out << "   mgr->SetDebugLevel(10);" << endl;
         out << "   mgr->SetNSysInfo(100);" << endl;
      }
      out << "   mgr->PrintStatus();" << endl;
      if (AliAnalysisManager::GetAnalysisManager()) {
         if (AliAnalysisManager::GetAnalysisManager()->GetDebugLevel()>3) {
            out << "   gEnv->SetValue(\"XNet.Debug\", \"1\");" << endl;
         } else {
            if (TestBit(AliAnalysisGrid::kTest))            
               out << "   AliLog::SetGlobalLogLevel(AliLog::kWarning);" << endl;
            else
               out << "   AliLog::SetGlobalLogLevel(AliLog::kError);" << endl;
         }
      }   
      if (!IsLocalTest()) {
         out << "   TChain *chain = CreateChain(\"wn.xml\", anatype);" << endl << endl;   
         out << "   mgr->StartAnalysis(\"localfile\", chain);" << endl;
      } else {
         out << "   mgr->StartAnalysis(\"localfile\");" << endl;
      }   
      out << "   timer.Stop();" << endl;
      out << "   timer.Print();" << endl;
      out << "}" << endl << endl;
      if (!IsLocalTest()) {
         out <<"//________________________________________________________________________________" << endl;
         out << "TChain* CreateChain(const char *xmlfile, const char *type=\"ESD\")" << endl;
         out << "{" << endl;
         out << "// Create a chain using url's from xml file" << endl;
         out << "   TString filename;" << endl;
         out << "   Int_t run = 0;" << endl;
         if (IsUseMCchain()) {
            out << "   TString treename = \"TE\";" << endl;
         } else {   
            out << "   TString treename = type;" << endl;
            out << "   treename.ToLower();" << endl;
            out << "   treename += \"Tree\";" << endl;
         }   
         out << "   printf(\"***************************************\\n\");" << endl;
         out << "   printf(\"    Getting chain of trees %s\\n\", treename.Data());" << endl;
         out << "   printf(\"***************************************\\n\");" << endl;
         out << "   TAlienCollection *coll = TAlienCollection::Open(xmlfile);" << endl;
         out << "   if (!coll) {" << endl;
         out << "      ::Error(\"CreateChain\", \"Cannot create an AliEn collection from %s\", xmlfile);" << endl;
         out << "      return NULL;" << endl;
         out << "   }" << endl;
         out << "   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();" << endl;
         out << "   TChain *chain = new TChain(treename);" << endl;
         if(fFriendChainName!="") {
            out << "   TChain *chainFriend = new TChain(treename);" << endl;
         }
         out << "   coll->Reset();" << endl;
         out << "   while (coll->Next()) {" << endl;
         out << "      filename = coll->GetTURL("");" << endl;
         out << "      if (mgr) {" << endl;
         out << "         Int_t nrun = AliAnalysisManager::GetRunFromAlienPath(filename);" << endl;
         out << "         if (nrun && nrun != run) {" << endl;
         out << "            printf(\"### Run number detected from chain: %d\\n\", nrun);" << endl;
         out << "            mgr->SetRunFromPath(nrun);" << endl;
         out << "            run = nrun;" << endl;
         out << "         }" << endl;
         out << "      }" << endl;
         out << "      chain->Add(filename);" << endl;
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
      if (hasANALYSISalice) {
         out <<"//________________________________________________________________________________" << endl;
         out << "Bool_t SetupPar(const char *package) {" << endl;
         out << "// Compile the package and set it up." << endl;
         out << "   TString pkgdir = package;" << endl;
         out << "   pkgdir.ReplaceAll(\".par\",\"\");" << endl;
         out << "   gSystem->Exec(TString::Format(\"tar xvzf %s.par\", pkgdir.Data()));" << endl;
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
   if (fProductionMode || TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      workdir += fGridWorkingDir;
      if (FileExists(fAnalysisMacro)) gGrid->Rm(fAnalysisMacro);
      Info("WriteAnalysisMacro", "\n#####   Copying analysis macro: <%s> to your alien workspace", fAnalysisMacro.Data());
//      TFile::Cp(Form("file:%s",fAnalysisMacro.Data()), Form("alien://%s/%s", workdir.Data(), fAnalysisMacro.Data()));
      if (!copyLocal2Alien("WriteAnalysisMacro",fAnalysisMacro.Data(), 
           Form("alien://%s/%s", workdir.Data(), 
           fAnalysisMacro.Data()))) Fatal("","Terminating");
   }
}

//______________________________________________________________________________
void AliAnalysisAlien::WriteMergingMacro()
{
// Write a macro to merge the outputs per master job.
   if (!fMergeViaJDL) return;
   if (!fOutputFiles.Length()) {
      Error("WriteMergingMacro", "No output file names defined. Are you running the right AliAnalysisAlien configuration ?");
      return;
   }   
   TString mergingMacro = fExecutable;
   mergingMacro.ReplaceAll(".sh","_merge.C");
   if (!fGridOutputDir.Contains("/")) fGridOutputDir = Form("%s/%s/%s", gGrid->GetHomeDirectory(), fGridWorkingDir.Data(), fGridOutputDir.Data());
   if (!TestBit(AliAnalysisGrid::kSubmit)) {  
      ofstream out;
      out.open(mergingMacro.Data(), ios::out);
      if (!out.good()) {
         Error("WriteMergingMacro", "could not open file %s for writing", fAnalysisMacro.Data());
         return;
      }
      Bool_t hasSTEERBase = kFALSE;
      Bool_t hasESD = kFALSE;
      Bool_t hasAOD = kFALSE;
      Bool_t hasANALYSIS = kFALSE;
      Bool_t hasOADB = kFALSE;
      Bool_t hasANALYSISalice = kFALSE;
      Bool_t hasCORRFW = kFALSE;
      TString func = mergingMacro;
      TString comment;
      func.ReplaceAll(".C", "");
      out << "void " << func.Data() << "(const char *dir, Int_t stage=0)" << endl;
      out << "{" << endl;
      out << "// Automatically generated merging macro executed in grid subjobs" << endl << endl;
      out << "   TStopwatch timer;" << endl;
      out << "   timer.Start();" << endl << endl;
      // Reset existing include path
      out << "// Reset existing include path and add current directory first in the search" << endl;
      out << "   gSystem->SetIncludePath(\"-I.\");" << endl;
      if (!fExecutableCommand.Contains("aliroot")) {
         out << "// load base root libraries" << endl;
         out << "   gSystem->Load(\"libTree\");" << endl;
         out << "   gSystem->Load(\"libGeom\");" << endl;
         out << "   gSystem->Load(\"libVMC\");" << endl;
         out << "   gSystem->Load(\"libPhysics\");" << endl << endl;
         out << "   gSystem->Load(\"libMinuit\");" << endl << endl;
      }   
      if (fAdditionalRootLibs.Length()) {
         // in principle libtree /lib geom libvmc etc. can go into this list, too
         out << "// Add aditional libraries" << endl;
         TObjArray *list = fAdditionalRootLibs.Tokenize(" ");
         TIter next(list);
         TObjString *str;
         while((str=(TObjString*)next())) {
            if (str->GetString().Contains(".so"))
            out << "   gSystem->Load(\"" << str->GetString().Data() << "\");" << endl;
         }
         if (list) delete list;
      }
      out << "// Load analysis framework libraries" << endl;
      if (!fPackages) {
         if (!fExecutableCommand.Contains("aliroot")) {
            out << "   gSystem->Load(\"libSTEERBase\");" << endl;
            out << "   gSystem->Load(\"libESD\");" << endl;
            out << "   gSystem->Load(\"libAOD\");" << endl;
         }
         out << "   gSystem->Load(\"libANALYSIS\");" << endl;
         out << "   gSystem->Load(\"libOADB\");" << endl;
         out << "   gSystem->Load(\"libANALYSISalice\");" << endl;
         out << "   gSystem->Load(\"libCORRFW\");" << endl << endl;
      } else {
         TIter next(fPackages);
         TObject *obj;
         TString pkgname;
         TString setupPar = "AliAnalysisAlien::SetupPar";
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
            if (pkgname == "OADB" ||
                pkgname == "OADB.par")      hasOADB = kTRUE;
            if (pkgname == "ANALYSISalice" ||
                pkgname == "ANALYSISalice.par") hasANALYSISalice = kTRUE;
            if (pkgname == "CORRFW" ||
                pkgname == "CORRFW.par")    hasCORRFW = kTRUE;
         }   
         if (hasANALYSISalice) setupPar = "SetupPar";   
         if (!hasSTEERBase) out << "   gSystem->Load(\"libSTEERBase\");" << endl;
         else out << "   if (!" << setupPar << "(\"STEERBase\")) return;" << endl;
         if (!hasESD)       out << "   gSystem->Load(\"libESD\");" << endl;
         else out << "   if (!" << setupPar << "(\"ESD\")) return;" << endl;
         if (!hasAOD)       out << "   gSystem->Load(\"libAOD\");" << endl;
         else out << "   if (!" << setupPar << "(\"AOD\")) return;" << endl;
         out << "   gSystem->Load(\"libOADB\");" << endl;
         if (!hasANALYSIS)  out << "   gSystem->Load(\"libANALYSIS\");" << endl;
         else out << "   if (!" << setupPar << "(\"ANALYSIS\")) return;" << endl;
         if (!hasOADB)  out << "   gSystem->Load(\"libOADB\");" << endl;
         else out << "   if (!" << setupPar << "(\"OADB\")) return;" << endl;
         if (!hasANALYSISalice)   out << "   gSystem->Load(\"libANALYSISalice\");" << endl;
         else out << "   if (!" << setupPar << "(\"ANALYSISalice\")) return;" << endl;
         if (!hasCORRFW)    out << "   gSystem->Load(\"libCORRFW\");" << endl << endl;
         else out << "   if (!" << setupPar << "(\"CORRFW\")) return;" << endl << endl;
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
                pkgname == "OADB" ||
                pkgname == "OADB.par" ||
                pkgname == "ANALYSISalice" ||
                pkgname == "ANALYSISalice.par" ||
                pkgname == "CORRFW" ||
                pkgname == "CORRFW.par") continue;
            out << "   if (!" << setupPar << "(\"" << obj->GetName() << "\")) return;" << endl;
         }   
      }   
      out << "// include path" << endl;
      // Get the include path from the interpreter and remove entries pointing to AliRoot
      out << "   TString intPath = gInterpreter->GetIncludePath();" << endl;
      out << "   TObjArray *listpaths = intPath.Tokenize(\" \");" << endl;
      out << "   TIter nextpath(listpaths);" << endl;
      out << "   TObjString *pname;" << endl;
      out << "   while ((pname=(TObjString*)nextpath())) {" << endl;
      out << "      TString current = pname->GetName();" << endl;
      out << "      if (current.Contains(\"AliRoot\") || current.Contains(\"ALICE_ROOT\")) continue;" << endl;
      out << "      gSystem->AddIncludePath(current);" << endl;
      out << "   }" << endl;
      out << "   if (listpaths) delete listpaths;" << endl;
      if (fIncludePath.Length()) out << "   gSystem->AddIncludePath(\"" << fIncludePath.Data() << "\");" << endl;
      out << "   gROOT->ProcessLine(\".include $ALICE_ROOT/include\");" << endl;
      out << "   printf(\"Include path: %s\\n\", gSystem->GetIncludePath());" << endl << endl;
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
      out << "// Analysis source to be compiled at runtime (if any)" << endl;
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

      if (fFastReadOption) {
         Warning("WriteMergingMacro", "!!! You requested FastRead option. Using xrootd flags to reduce timeouts in the grid merging jobs. Note that this may skip some files that could be accessed !!!");
         out << "// fast xrootd reading enabled" << endl;
         out << "   printf(\"!!! You requested FastRead option. Using xrootd flags to reduce timeouts. Note that this may skip some files that could be accessed !!!\");" << endl;
         out << "   gEnv->SetValue(\"XNet.ConnectTimeout\",50);" << endl;
         out << "   gEnv->SetValue(\"XNet.RequestTimeout\",50);" << endl;
         out << "   gEnv->SetValue(\"XNet.MaxRedirectCount\",2);" << endl;
         out << "   gEnv->SetValue(\"XNet.ReconnectTimeout\",50);" << endl;
         out << "   gEnv->SetValue(\"XNet.FirstConnectMaxCnt\",1);" << endl << endl;
      }
      // Change temp directory to current one
      out << "// Set temporary merging directory to current one" << endl;
      out << "   gSystem->Setenv(\"TMPDIR\", gSystem->pwd());" << endl << endl;   
      out << "// Set temporary compilation directory to current one" << endl;
      out << "   gSystem->SetBuildDir(gSystem->pwd(), kTRUE);" << endl << endl;   
      out << "// Connect to AliEn" << endl;
      out << "   if (!TGrid::Connect(\"alien://\")) return;" << endl;
      out << "   TString outputDir = dir;" << endl;  
      out << "   TString outputFiles = \"" << GetListOfFiles("out") << "\";" << endl;
      out << "   TString mergeExcludes = \"" << fMergeExcludes << "\";" << endl;
      out << "   TObjArray *list = outputFiles.Tokenize(\",\");" << endl;
      out << "   TIter *iter = new TIter(list);" << endl;
      out << "   TObjString *str;" << endl;
      out << "   TString outputFile;" << endl;
      out << "   Bool_t merged = kTRUE;" << endl;
      out << "   while((str=(TObjString*)iter->Next())) {" << endl;
      out << "      outputFile = str->GetString();" << endl;
      out << "      if (outputFile.Contains(\"*\")) continue;" << endl;
      out << "      Int_t index = outputFile.Index(\"@\");" << endl;
      out << "      if (index > 0) outputFile.Remove(index);" << endl;
      out << "      // Skip already merged outputs" << endl;
      out << "      if (!gSystem->AccessPathName(outputFile)) {" << endl;
      out << "         printf(\"Output file <%s> found. Not merging again.\",outputFile.Data());" << endl;
      out << "         continue;" << endl;
      out << "      }" << endl;
      out << "      if (mergeExcludes.Contains(outputFile.Data())) continue;" << endl;
      out << "      merged = AliAnalysisAlien::MergeOutput(outputFile, outputDir, " << fMaxMergeFiles << ", stage);" << endl;
      out << "      if (!merged) {" << endl;
      out << "         printf(\"ERROR: Cannot merge %s\\n\", outputFile.Data());" << endl;
      out << "         return;" << endl;
      out << "      }" << endl;
      out << "   }" << endl;
      out << "   // all outputs merged, validate" << endl;
      out << "   ofstream out;" << endl;
      out << "   out.open(\"outputs_valid\", ios::out);" << endl;
      out << "   out.close();" << endl;
      out << "   // read the analysis manager from file" << endl;
      TString analysisFile = fExecutable;
      analysisFile.ReplaceAll(".sh", ".root");
      out << "   if (!outputDir.Contains(\"Stage\")) return;" << endl;
      out << "   AliAnalysisManager *mgr = AliAnalysisAlien::LoadAnalysisManager(\"" 
          << analysisFile << "\");" << endl;
      out << "   if (!mgr) return;" << endl;
      out << "   mgr->SetRunFromPath(mgr->GetRunFromAlienPath(dir));" << endl;
      out << "   mgr->SetSkipTerminate(kFALSE);" << endl;
      out << "   mgr->PrintStatus();" << endl;
      if (AliAnalysisManager::GetAnalysisManager()) {
         if (AliAnalysisManager::GetAnalysisManager()->GetDebugLevel()>3) {
            out << "   gEnv->SetValue(\"XNet.Debug\", \"1\");" << endl;
         } else {
            if (TestBit(AliAnalysisGrid::kTest))            
               out << "   AliLog::SetGlobalLogLevel(AliLog::kWarning);" << endl;
            else
               out << "   AliLog::SetGlobalLogLevel(AliLog::kError);" << endl;
         }
      }   
      out << "   TTree *tree = NULL;" << endl;
      out << "   mgr->StartAnalysis(\"gridterminate\", tree);" << endl;
      out << "}" << endl << endl;
      if (hasANALYSISalice) {
         out <<"//________________________________________________________________________________" << endl;
         out << "Bool_t SetupPar(const char *package) {" << endl;
         out << "// Compile the package and set it up." << endl;
         out << "   TString pkgdir = package;" << endl;
         out << "   pkgdir.ReplaceAll(\".par\",\"\");" << endl;
         out << "   gSystem->Exec(TString::Format(\"tar xvzf %s.par\", pkgdir.Data()));" << endl;
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
   }   
   Bool_t copy = kTRUE;
   if (fProductionMode || TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      workdir += fGridWorkingDir;
      if (FileExists(mergingMacro)) gGrid->Rm(mergingMacro);
      Info("WriteMergingMacro", "\n#####   Copying merging macro: <%s> to your alien workspace", mergingMacro.Data());
//      TFile::Cp(Form("file:%s",mergingMacro.Data()), Form("alien://%s/%s", workdir.Data(), mergingMacro.Data()));
      if (!copyLocal2Alien("WriteMergeMacro",mergingMacro.Data(), 
           Form("%s/%s", workdir.Data(), mergingMacro.Data()))) Fatal("","Terminating");
   }
}

//______________________________________________________________________________
Bool_t AliAnalysisAlien::SetupPar(const char *package)
{
// Compile the par file archive pointed by <package>. This must be present in the current directory.
// Note that for loading the compiled library. The current directory should have precedence in
// LD_LIBRARY_PATH
   TString pkgdir = package;
   pkgdir.ReplaceAll(".par","");
   gSystem->Exec(TString::Format("tar xzf %s.par", pkgdir.Data()));
   TString cdir = gSystem->WorkingDirectory();
   gSystem->ChangeDirectory(pkgdir);
   // Check for BUILD.sh and execute
   if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("**************************************************\n");
      printf("*** Building PAR archive %s\n", package);
      printf("**************************************************\n");
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
         ::Error("SetupPar", "Cannot build par archive %s", pkgdir.Data());
         gSystem->ChangeDirectory(cdir);
         return kFALSE;
      }
   } else {
      ::Error("SetupPar","Cannot access PROOF-INF/BUILD.sh for package %s", pkgdir.Data());
      gSystem->ChangeDirectory(cdir);
      return kFALSE;
   }
   // Check for SETUP.C and execute
   if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("**************************************************\n");
      printf("*** Setup PAR archive %s\n", package);
      printf("**************************************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
      printf("*** Loaded library: %s\n", gSystem->GetLibraries(pkgdir,"",kFALSE));
   } else {
      ::Error("SetupPar","Cannot access PROOF-INF/SETUP.C for package %s", pkgdir.Data());
      gSystem->ChangeDirectory(cdir);
      return kFALSE;
   }   
   // Restore original workdir
   gSystem->ChangeDirectory(cdir);
   return kTRUE;
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
      // Make sure we can properly compile par files
      out << "export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH" << endl;
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
      out << "echo \"############## system limits : ##############\"" << endl;
      out << "ulimit -a" << endl;
      out << "echo \"############## memory : ##############\"" << endl;
      out << "free -m" << endl;
      out << "echo \"=========================================\"" << endl << endl;
      out << fExecutableCommand << " "; 
      out << fAnalysisMacro.Data() << " " << fExecutableArgs.Data() << endl << endl;
      out << "echo \"======== " << fAnalysisMacro.Data() << " finished with exit code: $? ========\"" << endl;
      out << "echo \"############## memory after: ##############\"" << endl;
      out << "free -m" << endl;
   }   
   Bool_t copy = kTRUE;
   if (fProductionMode || TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      TString bindir = Form("%s/bin", workdir.Data());
      if (!DirectoryExists(bindir)) gGrid->Mkdir(bindir,"-p");
      workdir += fGridWorkingDir;
      TString executable = Form("%s/bin/%s", gGrid->GetHomeDirectory(), fExecutable.Data());
      if (FileExists(executable)) gGrid->Rm(executable);
      Info("WriteExecutable", "\n#####   Copying executable file <%s> to your AliEn bin directory", fExecutable.Data());
//      TFile::Cp(Form("file:%s",fExecutable.Data()), Form("alien://%s", executable.Data()));
      if (!copyLocal2Alien("WriteExecutable",fExecutable.Data(), 
          executable.Data())) Fatal("","Terminating");
   } 
}

//______________________________________________________________________________
void AliAnalysisAlien::WriteMergeExecutable()
{
// Generate the alien executable script for the merging job.
   if (!fMergeViaJDL) return;
   TString mergeExec = fExecutable;
   mergeExec.ReplaceAll(".sh", "_merge.sh");
   if (!TestBit(AliAnalysisGrid::kSubmit)) {
      ofstream out;
      out.open(mergeExec.Data(), ios::out);
      if (out.bad()) {
         Error("WriteMergingExecutable", "Bad file name for executable: %s", mergeExec.Data());
         return;
      }
      out << "#!/bin/bash" << endl;
      // Make sure we can properly compile par files
      out << "export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH" << endl;
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
      out << "echo \"############## system limits : ##############\"" << endl;
      out << "ulimit -a" << endl;
      out << "echo \"############## memory : ##############\"" << endl;
      out << "free -m" << endl;
      out << "echo \"=========================================\"" << endl << endl;
      TString mergeMacro = fExecutable;
      mergeMacro.ReplaceAll(".sh", "_merge.C");
      if (IsOneStageMerging())
         out << "export ARG=\"" << mergeMacro << "(\\\"$1\\\")\"" << endl;
      else
         out << "export ARG=\"" << mergeMacro << "(\\\"$1\\\",$2)\"" << endl;
      out << fExecutableCommand << " " << "$ARG" << endl; 
      out << "echo \"======== " << mergeMacro.Data() << " finished with exit code: $? ========\"" << endl;
      out << "echo \"############## memory after: ##############\"" << endl;
      out << "free -m" << endl;
   }   
   Bool_t copy = kTRUE;
   if (fProductionMode || TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      TString bindir = Form("%s/bin", workdir.Data());
      if (!DirectoryExists(bindir)) gGrid->Mkdir(bindir,"-p");
      workdir += fGridWorkingDir;
      TString executable = Form("%s/bin/%s", gGrid->GetHomeDirectory(), mergeExec.Data());
      if (FileExists(executable)) gGrid->Rm(executable);
      Info("WriteMergeExecutable", "\n#####   Copying executable file <%s> to your AliEn bin directory", mergeExec.Data());
//      TFile::Cp(Form("file:%s",mergeExec.Data()), Form("alien://%s", executable.Data()));
      if (!copyLocal2Alien("WriteMergeExecutable",
          mergeExec.Data(), executable.Data())) Fatal("","Terminating");
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
   TString workdir;
   if (!fProductionMode && !fGridWorkingDir.BeginsWith("/alice"))
      workdir = gGrid->GetHomeDirectory();
   workdir += fGridWorkingDir;
   Int_t njobspermaster = 1000*fNrunsPerMaster/fSplitMaxInputFileNumber;
   TString locjdl = Form("%s/%s", workdir.Data(),fJDLName.Data());
   out << locjdl << " " << njobspermaster << endl;
   Int_t nmasterjobs = fInputFiles->GetEntries();
   for (Int_t i=0; i<nmasterjobs; i++) {
      TString runOutDir = gSystem->BaseName(fInputFiles->At(i)->GetName());
      runOutDir.ReplaceAll(".xml", "");
      if (fOutputToRunNo)
         out << Form("%s", fInputFiles->At(i)->GetName()) << " " << runOutDir << endl;
      else
         out << Form("%s", fInputFiles->At(i)->GetName()) << " " << Form("%03d", i) << endl;
   }
   if (gGrid) {
      Info("WriteProductionFile", "\n#####   Copying production file <%s> to your work directory", filename);
      if (FileExists(filename)) gGrid->Rm(filename);
//      TFile::Cp(Form("file:%s",filename), Form("alien://%s/%s", workdir.Data(),filename));
      if (!copyLocal2Alien("WriteProductionFile", filename, 
          Form("%s/%s", workdir.Data(),filename))) Fatal("","Terminating");
   }   
}

//______________________________________________________________________________
void AliAnalysisAlien::WriteValidationScript(Bool_t merge)
{
// Generate the alien validation script.
   // Generate the validation script
   TObjString *os;
   if (fValidationScript.IsNull()) {
      fValidationScript = fExecutable;
      fValidationScript.ReplaceAll(".sh", "_validation.sh");
   }   
   TString validationScript = fValidationScript;
   if (merge) validationScript.ReplaceAll(".sh", "_merge.sh");
   if (!Connect()) {
      Error("WriteValidationScript", "Alien connection required");
      return;
   }
   if (!fTerminateFiles.IsNull()) {
      fTerminateFiles.Strip();
      fTerminateFiles.ReplaceAll(" ",",");
   }   
   TString outStream = "";
   if (!TestBit(AliAnalysisGrid::kTest)) outStream = " >> stdout";
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
      out << "echo \"*******************************************************\"" << outStream << endl;
      out << "echo \"* Automatically generated validation script           *\""  << outStream << endl;
      out << "" << endl;
      out << "echo \"* Time:    $validatetime \""  << outStream << endl;
      out << "echo \"* Dir:     $validateout\""  << outStream << endl;
      out << "echo \"* Workdir: $validateworkdir\""  << outStream << endl;
      out << "echo \"* ----------------------------------------------------*\""  << outStream << endl;
      out << "ls -la ./"  << outStream << endl;
      out << "echo \"* ----------------------------------------------------*\""  << outStream << endl << endl;
      out << "##################################################" << endl;
      out << "" << endl;

      out << "if [ ! -f stderr ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"* ########## Job not validated - no stderr  ###\" " << outStream << endl;
      out << "   echo \"Error = $error\" " << outStream << endl;
      out << "fi" << endl;

      out << "parArch=`grep -Ei \"Cannot Build the PAR Archive\" stderr`" << endl;
      out << "segViol=`grep -Ei \"Segmentation violation\" stderr`" << endl;
      out << "segFault=`grep -Ei \"Segmentation fault\" stderr`" << endl;
      out << "glibcErr=`grep -Ei \"*** glibc detected ***\" stderr`" << endl;
      out << "" << endl;

      out << "if [ \"$parArch\" != \"\" ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"* ########## Job not validated - PAR archive not built  ###\" " << outStream << endl;
      out << "   echo \"$parArch\" " << outStream << endl;
      out << "   echo \"Error = $error\" " << outStream << endl;
      out << "fi" << endl;

      out << "if [ \"$segViol\" != \"\" ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"* ########## Job not validated - Segment. violation  ###\" " << outStream << endl;
      out << "   echo \"$segViol\" " << outStream << endl;
      out << "   echo \"Error = $error\" " << outStream << endl;
      out << "fi" << endl;

      out << "if [ \"$segFault\" != \"\" ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"* ########## Job not validated - Segment. fault  ###\" " << outStream << endl;
      out << "   echo \"$segFault\" " << outStream << endl;
      out << "   echo \"Error = $error\" " << outStream << endl;
      out << "fi" << endl;

      out << "if [ \"$glibcErr\" != \"\" ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"* ########## Job not validated - *** glibc detected ***  ###\" " << outStream << endl;
      out << "   echo \"$glibcErr\" " << outStream << endl;
      out << "   echo \"Error = $error\" " << outStream << endl;
      out << "fi" << endl;

      // Part dedicated to the specific analyses running into the train

      TString outputFiles = fOutputFiles;
      if (merge && !fTerminateFiles.IsNull()) {
         outputFiles += ",";
         outputFiles += fTerminateFiles;
      }
      TObjArray *arr = outputFiles.Tokenize(",");
      TIter next1(arr);
      TString outputFile;
      while (!merge && (os=(TObjString*)next1())) { 
         // No need to validate outputs produced by merging since the merging macro does this
         outputFile = os->GetString();
         Int_t index = outputFile.Index("@");
         if (index > 0) outputFile.Remove(index);
         if (fTerminateFiles.Contains(outputFile)) continue;
         if (outputFile.Contains("*")) continue;
         out << "if ! [ -f " << outputFile.Data() << " ] ; then" << endl;
         out << "   error=1" << endl;
         out << "   echo \"Output file " << outputFile << " not found. Job FAILED !\""  << outStream << endl;
         out << "   echo \"Output file " << outputFile << " not found. Job FAILED !\" >> stderr" << endl;
         out << "fi" << endl;
      }   
      delete arr;
      out << "if ! [ -f outputs_valid ] ; then" << endl;
      out << "   error=1" << endl;
      out << "   echo \"Output files were not validated by the analysis manager\" >> stdout" << endl;
      out << "   echo \"Output files were not validated by the analysis manager\" >> stderr" << endl;
      out << "fi" << endl;
      
      out << "if [ $error = 0 ] ; then" << endl;
      out << "   echo \"* ----------------   Job Validated  ------------------*\""  << outStream << endl;
      if (!IsKeepLogs()) {
         out << "   echo \"* === Logs std* will be deleted === \"" << endl;
         outStream = "";
         out << "   rm -f std*" << endl;
      }            
      out << "fi" << endl;

      out << "echo \"* ----------------------------------------------------*\""  << outStream << endl;
      out << "echo \"*******************************************************\""  << outStream << endl;
      out << "cd -" << endl;
      out << "exit $error" << endl;
   }    
   Bool_t copy = kTRUE;
   if (fProductionMode || TestBit(AliAnalysisGrid::kOffline) || TestBit(AliAnalysisGrid::kTest)) copy = kFALSE;
   if (copy) {
      CdWork();
      TString workdir = gGrid->GetHomeDirectory();
      workdir += fGridWorkingDir;
      Info("WriteValidationScript", "\n#####   Copying validation script <%s> to your AliEn working space", validationScript.Data());
      if (FileExists(validationScript)) gGrid->Rm(validationScript);
//      TFile::Cp(Form("file:%s",validationScript.Data()), Form("alien://%s/%s", workdir.Data(),validationScript.Data()));
      if (!copyLocal2Alien("WriteValidationScript", validationScript.Data(), 
          Form("%s/%s",workdir.Data(), validationScript.Data()))) Fatal("","Terminating");
   } 
}
