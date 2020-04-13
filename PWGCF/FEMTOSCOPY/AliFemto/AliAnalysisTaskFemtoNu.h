///
/// \file AliAnalysisTaskFemtoNu.h
///

#pragma once

// Used by the 'AddMultipleTasks.C' macro
#if defined(ROOT6_MACROLOAD_WORKAROUND) && !defined(ALIANALYSISTASKFEMTONU_H)
#define ALIANALYSISTASKFEMTONU_H

class AliAnalysisTask;
// class Macrocfg;


class AliAnalysisTaskFemtoNu {
public:
  static AliAnalysisTask* BuildForMacro(TString,TString,TString,TString);
};


#endif

#ifndef ALIANALYSISTASKFEMTONU_H
#define ALIANALYSISTASKFEMTONU_H

#include "AliAnalysisTaskFemto.h"
#include "AliFemtoResultStorage.h"


/// \class AliAnalysisTaskFemtoNu
/// \brief A class for handling and running AliFemtoAnalysis
///        objects, (writing output to a TDirectory)
///
/// Creates the necessary connection between the ESD or AOD input and
/// the femtoscopic code.
///
/// As with all [AliAnalysisTasks](@ref AliAnalysisTask), an
/// AliAnalysisTaskFemtoNu should be created from your macro and
/// added to the AliAnalysisManager.
/// The proper way to construct an AliAnalysisTaskFemto is through the
/// use of a configuration macro file which defines a function with
/// signature `AliFemtoManager* ConfigFemtoAnalysis()` whose name is
/// specified to the constructor.
/// The default name of this file is 'ConfigFemtoAnalysis.C'.
/// The macro is responsible for creating an AliFemtoAnalysisManager
/// and loading your custom AliFemtoAnalyses.
/// Upon the call to `::CreateOutputObjects`, the macro will be
/// loaded, and the analysis manager created and added to the
/// AliAnalysisTaskFemto.
/// Parameters can be specified in the constructor via the TString
/// paramter `aConfigParams`.
///
/// \author Andrew Kubera, Ohio State University <andrew.kubera@cern.ch>
///
class AliAnalysisTaskFemtoNu : public AliAnalysisTaskFemto {
public:

  /// Helper builder-pattern class used to construct tasks.
  ///
  /// Examples:
  ///
  /// ```cpp
  /// AliAnalysisTaskFemtoNu* task = AliAnalysisTaskFemtoNu::Builder()
  ///                                   .Verbose(false)
  ///                                   .Name("containername")
  ///                                   .Macro("/path/to/ConfigFemto.C")
  ///                                   .Verbose(false);
  /// ```
  ///
  struct Builder {
    TString name;
    TString macro;
    TString params;
    bool verbose;

    Builder()
      : name("femtoanalysis")
      , macro("ConfigFemtoAnalysis.C")
      , params("")
      , verbose(false)
    {
    }

    /// Construct with name
    Builder(const TString &name)
      : name(name)
      , macro("ConfigFemtoAnalysis.C")
      , params("")
      , verbose(false)
    {
    }

    Builder Name(TString n) const
      { Builder r(*this); r.name = n; return r; }

    Builder Macro(TString m) const
      { Builder r(*this); r.macro = m; return r; }

    Builder Params(TString m) const
      { Builder r(*this); r.params = m; return r; }

    Builder Verbose(bool v) const
      { Builder r(*this); r.verbose = v; return r; }

    AliAnalysisTaskFemtoNu* to_task() const
      { return params ? new AliAnalysisTaskFemtoNu(name, macro, verbose, params)
                      : new AliAnalysisTaskFemtoNu(name, macro, verbose); }

    operator AliAnalysisTaskFemtoNu*() const
      { return to_task(); }

  };

  struct MacroCfg : public TNamed {

    TString macro,
            auto_directory,
            output_filename,
            output_container,
            subwagon_array,
            subwagon_type;

    MacroCfg();
    virtual ~MacroCfg();


    TString GetFilename()
      {
        return (output_container == "")
             ? output_filename
             : TString::Format("%s:%s", output_filename.Data(), output_container.Data());
      }

  };

public:

  /// Default parameters initialize all values to zero (NULL) and empty strings
  ///
  /// Do not use this constructor - objects created with it have no way to set
  /// the configuration macro or parameters. Use an alternative constructor.
  AliAnalysisTaskFemtoNu();

  /// Full Constructor - Set the name of the task, configuration macro filename
  /// and paramters, and optional verbosity flag.
  AliAnalysisTaskFemtoNu(const TString &name,
                         const TString &aConfigMacro,
                         const TString &aConfigParams,
                         const Bool_t aVerbose=kFALSE);

  /// Construct with task name, configuration filename, and verbosity flag.
  ///
  /// The paramters are set to the empty string.
  AliAnalysisTaskFemtoNu(const TString &name,
                         const TString &aConfigMacro,
                         const Bool_t aVerbose=kFALSE);

  virtual ~AliAnalysisTaskFemtoNu();

  /// Automatically create and connect input/output containers
  ///
  /// The output container name will be the same as this task.
  ///
  /// \param filename The destination output filename.
  ///                 If empty this uses `mgr->GetCommonFileName()`
  ///
  virtual void SetupContainers(const TString &filename="");

  /// Run macro to generate the FemtoManager
  virtual void CreateOutputObjects();

  /// Connect Input Data
  virtual void ConnectInputData(Option_t *);

  /// Execute analysis
  virtual void Exec(Option_t *);

  void UseCustomExec()
    { fUseCustomExec = true; }

  /// Return manager
  AliFemtoManager *GetManager()
    { return fManager; }

  /// The output-slot number of the AliFemtoResultStorage object
  static const unsigned RESULT_STORAGE_OUTPUT_SLOT;

  // used to subvert ROOT6 macro loading bug
  static AliAnalysisTask* BuildForMacro(TString,TString,TString,TString);

protected:

  /// object used to write output objects to a TDirectory with this AnalysisTasks's name
  AliFemtoResultStorage *fStorage;

  bool fUseCustomExec;

private:
  /// non-copiable
#if __cplusplus < 201103L
  AliAnalysisTaskFemtoNu(const AliAnalysisTaskFemtoNu&);
  AliAnalysisTaskFemtoNu& operator=(const AliAnalysisTaskFemtoNu&);
#else
  AliAnalysisTaskFemtoNu(const AliAnalysisTaskFemtoNu&) = delete;
  AliAnalysisTaskFemtoNu& operator=(const AliAnalysisTaskFemtoNu&) = delete;
#endif

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskFemtoNu, 1);
  /// \endcond
};

#endif
