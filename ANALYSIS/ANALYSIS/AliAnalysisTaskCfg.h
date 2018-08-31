#ifndef ALIANALYSISTASKCFG_H
#define ALIANALYSISTASKCFG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliAnalysisTaskCfg
/// \brief AliAnalysisTaskCfg
/// ==============================================================================
///   AliAnalysysTaskCfg - Class embedding the configuration needed to run
/// a given analysis task: libraries to be loaded, location and name of the macro
/// used to add the task to the analysis manager, dependencies.
/// ==============================================================================
/// This class is used to fully describe how to run a given analysis task. It
/// requires that the user creates an AddTask macro for his task and defines:
///  - The needed libs separated by commas,
///  - The full path to the AddTask macro (starting with $ALICE_ROOT if needed)
///  - The list of arguments to be provided to the AddTask macro. One can use
///    here only constants that can be interpreted.
///  - The list of dependencies (other modules required to run this task). These
///    must be names of other AliAnalysisTaskCfg objects, separated by commas.
///  - Data types supported by the task (e.g. ESD, AOD, MC)
/// The class has normal ROOT IO, but it can also read from and write to text files.
/// An example:
/// Content of file: QAsym.cfg
/// The following special variable names can be used:
/// __R_ADDTASK__ = the return value of the AddTask macro included
/// __R_ESDH__    = pointer to ESD handler
/// __R_AODH__    = pointer to AOD handler
/// __R_MCH__     = pointer to MC handler
/// The static method ExtractModulesFrom(const char *filename) allows reading
/// several task configuration modules from the same text file and returning
/// them in a TObjArray.
///
/// A list of configuration modules representing a train should be injected in
/// the right order in the grid handler to generate train macros.
/// \author Andrei Gheata
/// \date 12/08/2011

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TMacro;
class TObjArray;

class AliAnalysisTaskCfg : public TNamed {
public:
enum ETaskCfgFlags {
  kLoaded      = BIT(14)
}; 
protected:
  TString                   fMacroName;     ///< Full path to AddTask macro
  TString                   fMacroArgs;     ///< Arguments to run the macro
  TString                   fLibs;          ///< List of custom libs needed to run the task (comma separated)
  TString                   fDeps;          ///< List of tasks this module depends on
  TString                   fDataTypes;     ///< List of supported data types (ESD, AOD, MC)
  TString                   fOutputFile;    ///< Desired output file name (via SetCommonFileName)
  TString                   fTerminateFile; ///< Custom output file written in Terminate
  TMacro                   *fMacro;         ///< Embedded AddTask macro
  TMacro                   *fConfigDeps;    ///< Macro used to configure the dependecies
                                            // (utility tasks or input handlers). The data type is passed as argument.
  TObject                  *fRAddTask;      ///< Object returned by AddTask method
public:  
  AliAnalysisTaskCfg();
  AliAnalysisTaskCfg(const char *name);
  AliAnalysisTaskCfg(const AliAnalysisTaskCfg &other); 
  virtual ~AliAnalysisTaskCfg();
  
  // Assignment
  AliAnalysisTaskCfg& operator=(const AliAnalysisTaskCfg &other);
  
  // AddTask macro handling
  const char               *GetMacroName() const {return fMacroName;}
  const char               *GetMacroArgs() const {return fMacroArgs;}
  void                      SetMacroName(const char *name) {fMacroName = name;}
  void                      SetMacroArgs(const char *args) {fMacroArgs = args;}
  TMacro                   *OpenMacro(const char *name="");
  void                      SetMacro(TMacro *macro);
  TMacro                   *GetMacro() const {return fMacro;}
  Long64_t                  ExecuteMacro(const char *newargs="");

  // Libraries
  const char               *GetLibs()      const {return fLibs;}
  Int_t                     GetNlibs()     const;
  const char *              GetLibrary(Int_t i) const;
  Bool_t                    NeedsLibrary(const char *lib) const;
  void                      SetLibraries(const char *libs) {fLibs = libs;}
  
  // Dependencies
  const char               *GetDeps()      const {return fDeps;}
  Int_t                     GetNdeps()     const;
  const char *              GetDependency(Int_t i) const;
  Bool_t                    NeedsDependency(const char *dep) const;
  void                      SetDependencies(const char *deps) {fDeps = deps;}
  
  // Customized macro to handle dependencies
  TMacro                   *OpenConfigMacro(const char *name);
  void                      SetConfigMacro(TMacro *macro);
  Long64_t                  ExecuteConfigMacro();
  TMacro                   *GetConfigMacro() const {return fConfigDeps;}
  
  // Supported data types
  const char               *GetDataTypes() const {return fDataTypes;}
  Bool_t                    SupportsData(const char *type) const;
  void                      SetDataTypes(const char *types);
  
  // Output files
  const char               *GetOutputFileName() const {return fOutputFile;}
  const char               *GetTerminateFileName() const {return fTerminateFile;}
  void                      SetOutputFileName(const char *name) {fOutputFile = name;}
  void                      SetTerminateFileName(const char *name) {fTerminateFile = name;}

  // Extra utilities  
  Bool_t                    CheckLoadLibraries() const;
  static const char        *DecodeValue(TString &line);
  TObject                  *GetRAddTask() const {return fRAddTask;}
  Bool_t                    IsLoaded() const {return TObject::TestBit(AliAnalysisTaskCfg::kLoaded);}
  void                      Print(Option_t *option="") const;
  void                      SaveAs(const char *filename, Option_t *option = "") const;
  static TObjArray         *ExtractModulesFrom(const char *filename);
    
  ClassDef(AliAnalysisTaskCfg,2)  // Class describing how to run a analysis task
};
#endif
