/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisTaskCfg.h"
#include "AliAnalysisManager.h"

#include "Riostream.h"
#include "TError.h"
#include "TMacro.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TList.h"

// Author: Andrei Gheata, 12/08/2011

//==============================================================================
//   AliAnalysysTaskCfg - Class embedding the configuration needed to run
// a given analysis task: libraries to be loaded, location and name of the macro
// used to add the task to the analysis manager, dependencies.
//==============================================================================

// This class is used to fully describe how to run a given analysis task. It
// requires that the user creates an AddTask macro for his task and defines:
//  - The needed libs separated by commas,
//  - The full path to the AddTask macro (starting with $ALICE_ROOT if needed)
//  - The list of arguments to be provided to the AddTask macro. One can use
//    here only constants that can be interpreted.
//  - The list of dependencies (other modules required to run this task). These
//    must be names of other AliAnalysisTaskCfg objects, separated by commas.
//  - Data types supported by the task (e.g. ESD, AOD, MC)
// The class has normal ROOT IO, but it can also read from and write to text files.
// An example:
// Content of file: QAsym.cfg
/*
# Lines that do not start with #Module are ignored, except those in embedded
  macro blocks
#Module.Begin        QAsym
#Module.Libs         PWG1
#Module.Deps         PhysicsSelection
#Module.DataTypes    ESD, AOD, MC
#Module.MacroName    $ALICE_ROOT/PWG1/PilotTrain/AddTaskQAsym.C
#Module.MacroArgs    0, AliVEvent::kAnyINT, AliVEvent::kHighMult, AliVEvent::kEMC7, AliVEvent::kMUU7
#Module.StartConfig  
__R_ADDTASK__->SelectCollisionCandidates();
#Module.EndConfig
*/
// The following special variable names can be used:
// __R_ADDTASK__ = the return value of the AddTask macro included
// __R_ESDH__    = pointer to ESD handler
// __R_AODH__    = pointer to AOD handler
// __R_MCH__     = pointer to MC handler
// The static method ExtractModulesFrom(const char *filename) allows reading
// several task configuration modules from the same text file and returning
// them in a TObjArray.
//
// A list of configuration modules representing a train should be injected in
// the right order in the grid handler to generate train macros.


ClassImp(AliAnalysisTaskCfg)

//______________________________________________________________________________
AliAnalysisTaskCfg::AliAnalysisTaskCfg()
                   :TNamed(),
                    fMacroName(),
                    fMacroArgs(),
                    fLibs(),
                    fDeps(),
                    fDataTypes(),
                    fMacro(0),
                    fConfigDeps(0),
                    fRAddTask(0)
{
// I/O constructor.
}

//______________________________________________________________________________
AliAnalysisTaskCfg::AliAnalysisTaskCfg(const char *name)
                   :TNamed(name,""),
                    fMacroName(),
                    fMacroArgs(),
                    fLibs(),
                    fDeps(),
                    fDataTypes(),
                    fMacro(0),
                    fConfigDeps(0),
                    fRAddTask(0)
{
// Constructor. All configuration objects need to be named since they are looked
// for by name.
}

//______________________________________________________________________________
AliAnalysisTaskCfg::AliAnalysisTaskCfg(const AliAnalysisTaskCfg &other)
                   :TNamed(other),
                    fMacroName(other.fMacroName),
                    fMacroArgs(other.fMacroArgs),
                    fLibs(other.fLibs),
                    fDeps(other.fDeps),
                    fDataTypes(other.fDataTypes),
                    fMacro(0),
                    fConfigDeps(0),
                    fRAddTask(0)
{
// Copy constructor.
   if (other.fMacro) fMacro = new TMacro(*other.fMacro);
   if (other.fConfigDeps) fConfigDeps = new TMacro(*other.fConfigDeps);
}   

//______________________________________________________________________________
AliAnalysisTaskCfg::~AliAnalysisTaskCfg()
{
// Destructor.
   delete fMacro;
   delete fConfigDeps;
}

//______________________________________________________________________________
AliAnalysisTaskCfg& AliAnalysisTaskCfg::operator=(const AliAnalysisTaskCfg &other)
{
// Assignment operator.
   if (&other == this) return *this;
   TNamed::operator=(other);
   fMacroName = other.fMacroName;
   fMacroArgs = other.fMacroArgs;
   fLibs      = other.fLibs;
   fDeps      = other.fDeps;
   fDataTypes = other.fDataTypes;
   if (other.fMacro) fMacro = new TMacro(*other.fMacro);
   if (other.fConfigDeps) fConfigDeps = new TMacro(*other.fConfigDeps);
   fRAddTask  = other.fRAddTask;
   return *this;
}
   
//______________________________________________________________________________
TMacro *AliAnalysisTaskCfg::OpenMacro(const char *name)
{
// Opens the specified macro if name is not empty. In case of success updates
// fMacroName, creates the maco object and returns its pointer.
   // Clean-up previous macro if any
   if (fMacro) {
      delete fMacro;
      fMacro = 0;
   }   
   TString expname;
   if (strlen(name)) expname = gSystem->ExpandPathName(name);
   else              expname = gSystem->ExpandPathName(fMacroName.Data());
   if (expname.IsNull()) {
      Error("OpenMacro", "Macro name not provided and not previously set");
      return 0;
   }   
   if (gSystem->AccessPathName(expname)) {
      Error("OpenMacro", "Macro: %s cannot be opened.", expname.Data());
      return 0;
   }
   if (strlen(name)) fMacroName = name;
   fMacro = new TMacro(expname);
   return fMacro;
}   

//______________________________________________________________________________
void AliAnalysisTaskCfg::SetMacro(TMacro *macro)
{
// Set the AddTask macro from outside. This will discard the existing macro if
// any. The provided macro will be owned by this object.
   if (fMacro) delete fMacro;
   fMacro = macro;
}   

//______________________________________________________________________________
Long64_t AliAnalysisTaskCfg::ExecuteMacro(const char *newargs)
{
// Execute AddTask macro. Opens first the macro pointed by fMacroName if not yet
// done. Checks if the requested libraries are loaded, else loads them. Executes 
// with stored fMacroArgs unless new arguments are provided. The flag IsLoaded
// is set once the macro was successfully executed.
   if (IsLoaded()) return kTRUE;
   if (!fMacro && !OpenMacro()) {
      Error("ExecuteMacro", "Cannot execute this macro");
      return -1;
   }
   if (!CheckLoadLibraries()) {
      Error("ExecuteMacro", "Cannot load requested libraries: %s", fLibs.Data());
      return -1;
   }
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("ExecuteMacro", "Analysis manager not defined yet");
      return -1;
   }
   Int_t ntasks0 = mgr->GetTasks()->GetEntriesFast();
   TString args = newargs;
   if (args.IsNull()) args = fMacroArgs;
   Int_t error = 0;
   Long64_t retval = fMacro->Exec(args, &error);
   if (error != TInterpreter::kNoError)
   {
      Error("ExecuteMacro", "Macro interpretation failed");
      return -1;
   }
   Int_t ntasks = mgr->GetTasks()->GetEntriesFast();
   if (ntasks<=ntasks0)
   {
      Error("ExecuteMacro", "The macro did not add any tasks to the manager");
      return -1;
   }
   Long64_t ptrTask = (Long64_t)mgr->GetTasks()->At(ntasks0);
   if (retval >= ptrTask) {
      TObject::SetBit(AliAnalysisTaskCfg::kLoaded, kTRUE);
      fRAddTask = reinterpret_cast<TObject*>(retval);
      if (fConfigDeps && dynamic_cast<TObject*>(fRAddTask)) {
         TString classname = fRAddTask->ClassName();
         classname += Form("* __R_ADDTASK__ = (%s*)0x%lx;", classname.Data(),(ULong_t)retval);
         classname.Prepend("  ");
         TObjString *line = fConfigDeps->GetLineWith("__R_ADDTASK__");
         if (line) {
            TList *lines = fConfigDeps->GetListOfLines();
            lines->AddBefore(line, new TObjString(classname));
         }
      }   
   }
   Info("ExecuteMacro", "Macro %s added %d tasks to the manager", fMacro->GetName(), ntasks-ntasks0);
   return retval;
}

//______________________________________________________________________________
Int_t AliAnalysisTaskCfg::GetNlibs() const
{
// Returns number of requested libraries.
   if (fLibs.IsNull()) return 0;
   TObjArray *list  = fLibs.Tokenize(",");
   Int_t nlibs = list->GetEntriesFast();
   delete list;
   return nlibs;
}

//______________________________________________________________________________
const char *AliAnalysisTaskCfg::GetLibrary(Int_t i) const
{
// Returns library name for the i-th library.
   Int_t nlibs = GetNlibs();
   if (i>=nlibs) return 0;
   static TString libname;
   TObjArray *list  = fLibs.Tokenize(",");
   libname = list->At(i)->GetName();
   libname.ReplaceAll(".so","");
   libname.ReplaceAll(" ","");
   if (libname.BeginsWith("lib")) libname.Remove(0, 3);
   delete list;
   return libname.Data();
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskCfg::CheckLoadLibraries() const
{
// Check if all requested libraries were loaded, otherwise load them. If some
// library cannot be loaded return false.
   TString library;
   Int_t nlibs = GetNlibs();
   for (Int_t i=0; i<nlibs; i++) {
      library = GetLibrary(i);
      library.Prepend("lib");
      Int_t loaded = strlen(gSystem->GetLibraries(library,"",kFALSE));
      if (!loaded) loaded = gSystem->Load(library);
      if (loaded < 0) {
         Error("CheckLoadLibraries", "Cannot load library %s", library.Data());
         return kFALSE;
      }
   }
   return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskCfg::NeedsLibrary(const char *lib) const
{
// Check if a given library is needed by the module.
   TString libname = lib;
   libname.ReplaceAll(".so","");
   if (libname.BeginsWith("lib")) libname.Remove(0, 3);
   return fLibs.Contains(libname);
}

//______________________________________________________________________________
Int_t AliAnalysisTaskCfg::GetNdeps() const
{
// Returns number of requested libraries.
   if (fDeps.IsNull()) return 0;
   Int_t ndeps = fDeps.CountChar(',')+1;
   return ndeps;
}

//______________________________________________________________________________
const char *AliAnalysisTaskCfg::GetDependency(Int_t i) const
{
// Returns library name for the i-th library.
   Int_t ndeps = GetNdeps();
   if (i>=ndeps) return 0;
   static TString depname;
   TObjArray *list  = fDeps.Tokenize(",");
   depname = list->At(i)->GetName();
   depname.ReplaceAll(" ","");
   delete list;
   return depname.Data();
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskCfg::NeedsDependency(const char *dep) const
{
// Check if a given library is needed by the module.
   return fDeps.Contains(dep);
}

//______________________________________________________________________________
TMacro *AliAnalysisTaskCfg::OpenConfigMacro(const char *name)
{
// Opens the specified macro if name is not empty.
   if (fConfigDeps) {
      delete fConfigDeps;
      fConfigDeps = 0;
   }
   
   TString expname = gSystem->ExpandPathName(name);
   if (expname.IsNull()) {
      Error("OpenConfigMacro", "Macro name not provided");
      return 0;
   }   
   if (gSystem->AccessPathName(expname)) {
      Error("OpenMacro", "Macro: %s cannot be opened.", expname.Data());
      return 0;
   }
   fConfigDeps = new TMacro(expname);
   return fConfigDeps;
}   

//______________________________________________________________________________
void AliAnalysisTaskCfg::SetConfigMacro(TMacro *macro)
{
// Set the macro for configuring deps from outside. This will discard the 
// existing macro if any. The provided macro will be owned by this object.
   if (fConfigDeps) delete fConfigDeps;
   fConfigDeps = macro;
}   
  
//______________________________________________________________________________
Long64_t AliAnalysisTaskCfg::ExecuteConfigMacro()
{
// Execute macro to configure dependencies. No arguments are supported.
   if (!fConfigDeps) {
      Error("ExecuteConfigMacro", "Call OpenConfigMacro() first");
      return -1;
   }
   if (!CheckLoadLibraries()) {
      Error("ExecuteConfigMacro", "Cannot load requested libraries: %s", fLibs.Data());
      return -1;
   }
   return fConfigDeps->Exec();
}

//______________________________________________________________________________
void AliAnalysisTaskCfg::SetDataTypes(const char *types)
{
// Sets the data types supported by the module. Stored in upper case.
   fDataTypes = types;
   fDataTypes.ToUpper();
}

//______________________________________________________________________________
Bool_t AliAnalysisTaskCfg::SupportsData(const char *type) const
{
// Checks if the given data type is supported.
   TString stype = type;
   stype.ToUpper();
   return fDataTypes.Contains(stype);
}

//______________________________________________________________________________
void AliAnalysisTaskCfg::Print(Option_t * option) const
{
// Print content of the module.
   TString opt(option);
   Bool_t full = (opt.Length())?kTRUE:kFALSE;
   printf("====================================================================\n");
   printf("# Analysis task:                %s\n", GetName());
   printf("# Supported data types:         %s\n", fDataTypes.Data());
   printf("# Extra libraries:              %s\n", fLibs.Data());
   printf("# Extra dependencies:           %s\n", fDeps.Data());
   if (fConfigDeps) {
      printf("# Macro to configure deps:      %s\n", fConfigDeps->GetTitle());
      if (full) fConfigDeps->Print();
   }   
   printf("# Macro connecting this task:   %s\n", fMacroName.Data());
   printf("# Arguments to run the macro:   %s\n", fMacroArgs.Data());
   if (full) {
      if (fMacro) fMacro->Print();
      else {
         TMacro macro(gSystem->ExpandPathName(fMacroName.Data()));
         macro.Print();
      }   
   }   
}
   
//______________________________________________________________________________
void AliAnalysisTaskCfg::SaveAs(const char *filename, Option_t *option) const
{
// Save the configuration module as text file in the form key:value. The
// option can be APPEND, otherwise the file gets overwritten.
   TString opt(option);
   opt.ToUpper();
   ios::openmode mode = ios::out;
   if (opt == "APPEND") mode = ios::app;
   ofstream out;
   out.open(filename, mode);
   if (out.bad()) {
      Error("SaveAs", "Bad file name: %s", filename);
      return;
   }
   out << "#Module.Begin " << GetName() << endl;
   out << "#Module.Libs " << fLibs << endl;
   out << "#Module.Deps " << fDeps << endl;
   out << "#Module.DataTypes " << fDataTypes << endl;
   out << "#Module.MacroName " << fMacroName << endl;
   out << "#Module.MacroArgs " << fMacroArgs << endl;
   if (fConfigDeps) {
      out << "#Config.Deps " << fConfigDeps->GetTitle() << endl;
   }      
}


//______________________________________________________________________________
const char *AliAnalysisTaskCfg::DecodeValue(TString &line)
{
// Decode the value string from the line
   static TString value;
   value = line(line.Index(' '),line.Length());
   value = value.Strip(TString::kLeading,' ');
   value = value.Strip(TString::kTrailing,' ');
   return value.Data();
}
   
//______________________________________________________________________________
TObjArray *AliAnalysisTaskCfg::ExtractModulesFrom(const char *filename)
{
// Read all modules from a text file and add them to an object array. The
// caller must delete the array at the end. Any module must start with a line
// containing: #Module.Begin
   TString expname = gSystem->ExpandPathName(filename);
   if (gSystem->AccessPathName(expname)) {
      ::Error("ExtractModulesFrom", "Cannot open file %s", filename);
      return 0;
   }   
   AliAnalysisTaskCfg *cfg = 0;
   TObjArray *array = 0;
   ifstream in;
   in.open(expname);
   char cline[1024];
   TString line;
   TString decode;
   TMacro *addMacro = 0;
   TMacro *addConfig = 0;
   while (in.good()) {
      in.getline(cline,1024);
      line = cline;
      if (line.BeginsWith("#Module.Begin")) {
      // New module found, save previous if any
         if (cfg) {
            if (addMacro) cfg->SetMacro(addMacro);
            if (addConfig) cfg->SetConfigMacro(addConfig);
            if (!array) array = new TObjArray();
            array->Add(cfg);
         }
         // Decode module name from the line
         decode = AliAnalysisTaskCfg::DecodeValue(line);
         cfg = new AliAnalysisTaskCfg(decode);
         addMacro = 0;
         addConfig = 0;
      } else if (cfg && line.BeginsWith("#Module.Libs")) {
         // Libraries section
         decode = AliAnalysisTaskCfg::DecodeValue(line);
         cfg->SetLibraries(decode);
      } else if (cfg && line.BeginsWith("#Module.Deps")) {
         // Dependencies section
         decode = AliAnalysisTaskCfg::DecodeValue(line);
         cfg->SetDependencies(decode);
      } else if (cfg && line.BeginsWith("#Module.DataTypes")) {
         // Data types
         decode = AliAnalysisTaskCfg::DecodeValue(line);
         cfg->SetDataTypes(decode);
      } else if (cfg && line.BeginsWith("#Module.MacroName")) {
         // Name of the add task macro (including path)
         decode = AliAnalysisTaskCfg::DecodeValue(line);
         cfg->SetMacroName(decode);
      } else if (cfg && line.BeginsWith("#Module.MacroArgs")) {
         // Arguments for the AddTask macro
         decode = AliAnalysisTaskCfg::DecodeValue(line);
         cfg->SetMacroArgs(decode);
      } else if (cfg && line.BeginsWith("#Module.StartMacro")) {
         // Marker for start of the AddTask macro
         addMacro = new TMacro();
         TString shortName = gSystem->BaseName(cfg->GetMacroName());
         shortName = shortName(0,shortName.Index("."));
         addMacro->SetName(shortName);
         addMacro->SetTitle(cfg->GetMacroName());
      } else if (cfg && line.BeginsWith("#Module.StartConfig")) {
         // Marker for start of the configuration macro
         addConfig = new TMacro();
         TString shortName = gSystem->BaseName(cfg->GetMacroName());
         shortName = shortName(0,shortName.Index("."));
         shortName += "Config";
         addConfig->SetName(shortName);
         shortName.Prepend("/");
         shortName.Prepend(gSystem->DirName(cfg->GetMacroName()));
         shortName += ".C";
         addConfig->SetTitle(shortName);
      } else if (cfg && line.BeginsWith("#Module.EndMacro")) {
         // Marker for the end of the embedded macro. EndMacro block mandatory.
         if (cfg && addMacro) {
            cfg->SetMacro(addMacro);
            addMacro = 0;
         } else {
            ::Error("ExtractModulesFrom", "Canot end un-started macro block");
         }
      } else if (cfg && line.BeginsWith("#Module.EndConfig")) {
         // Marker for the end of the config macro. EndConfig block is mandatory
         if (cfg && addConfig) {
            addConfig->GetListOfLines()->AddFirst(new TObjString(Form("%s() {",gSystem->BaseName(addConfig->GetName()))));
            addConfig->GetListOfLines()->AddLast(new TObjString("}"));
            cfg->SetConfigMacro(addConfig);
            addConfig = 0;
         } else {
            ::Error("ExtractModulesFrom", "Canot end un-started config macro block");
         }
      } else {   
         // Reading a block line
         if (addMacro) addMacro->AddLine(line);
         else if (addConfig) addConfig->AddLine(line);
      }   
   }
   // Add last found object to the list
   if (cfg) {
      if (addMacro) ::Error("ExtractModulesFrom", "#Module.EndMacro block not found");         
      if (addConfig) ::Error("ExtractModulesFrom", "#Module.EndConfig block not found");
      if (!array) array = new TObjArray();
      array->Add(cfg);
   }
   return array;  
}
