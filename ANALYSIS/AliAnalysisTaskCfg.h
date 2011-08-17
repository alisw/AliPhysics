#ifndef ALIANALYSISTASKCFG_H
#define ALIANALYSISTASKCFG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Andrei Gheata, 12/08/2011

//==============================================================================
//   AliAnalysysTaskCfg - Class embedding the configuration needed to run
// a given analysis task: libraries to be loaded, location and name of the macro
// used to add the task to the analysis manager, dependencies.
//==============================================================================

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TMacro;
class TObjArray;

class AliAnalysisTaskCfg : public TNamed {
protected:
  TString                   fMacroName;     // Full path to AddTask macro
  TString                   fMacroArgs;     // Arguments to run the macro
  TString                   fLibs;          // List of custom libs needed to run the task (comma separated)
  TString                   fDeps;          // List of tasks this module depends on
  TString                   fDataTypes;     // List of supported data types (ESD, AOD, MC)
  TMacro                   *fMacro;         // Embedded AddTask macro
  TMacro                   *fConfigDeps;    // Macro used to configure the dependecies
                                            // (utility tasks or input handlers). The data type is passed as argument.
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

  // Extra utilities  
  Bool_t                    CheckLoadLibraries() const;
  static const char        *DecodeValue(TString &line);
  void                      Print(Option_t *option="") const;
  void                      SaveAs(const char *filename, Option_t *option = "") const;
  static TObjArray         *ExtractModulesFrom(const char *filename);
    
  ClassDef(AliAnalysisTaskCfg,1)  // Class describing how to run a analysis task
};
#endif
