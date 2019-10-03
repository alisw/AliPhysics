///
/// \file AddTaskPiLam.C
/// \author Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///

///
/// \brief Adds an AliAnalysisTaskFemto analysis object to the global
///        AliAnalysisManager.
///
/// This macro creates and returns an AliAnalysisTaskFemto object. The task
/// is given the config macro "Train/PionLambdaFemto/ConfigFemtoAnalysis.C"
/// which is run to create the analysis objects. This is fixed (for now), and
/// if an alternative is required, you should use the general AddTaskFemto.C
/// macro.
///
/// Subwagons are supported, the name of which will be appended to the task
/// name. The task name by default is "TaskPiLamFemto" which cannot be changed.
///
/// The output container name is fixed at the standard "femtolist".
///
/// \param params A string forwarded to the ConfigFemtoAnalysis macro for
///               parsing. This string is wrapped in double-quotes so escaping
///               some to ensure results is a string is unneccessary.
///
/// \param maco_filename The path to the femto-confuration macro, passed to the
///                      created TaskFemto. If this parameter has no '/' characters,
///                      and the subwagon_suffix is empty, this is reinterpreted
///                      as the subwagon_suffix.
///
/// \param subwagon_suffix If this macro is run in a train with subwagons, this
///                        will be set with the identifier. NOTE: This
///                        parameter may be found in the macro_filename
///                        variable as the  suffix is simply the last argument
///                        passed. To keep from always having to explicitly set
///                        the macro_filename if using the default, this code
///                        will switch the two parameters if there is no '/'
///                        characters in the macro_filename. As such, keep '/'
///                        out of the subwagon_suffix!
///
AliAnalysisTaskFemto* AddTaskPiLam(TString params,
                                   TString macro_filename="",
                                   TString output_filename="",
                                   TString output_container="PWG2FEMTO",
                                   TString subwagon_suffix="")
{ // Adds a Pion-Lambda Femtoscopy task to the manager

  const TString DEFAULT_MACRO = "$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/Train/PionLambdaFemto/ConfigFemtoAnalysis.C";

  // Get the global manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPiLam", "Could not get the global AliAnalysisManager.");
    return NULL;
  }

  const TString default_name = "TaskPiLamFemto";

  // If the macro_filename was set, use it, else use the default
  if (macro_filename == "") {
    macro_filename = DEFAULT_MACRO;
  }
  // If subwagon_suffix was not set and there are no '/' characters in the
  // macro's path, interpret path as the subwagon_suffix and set default path
  if (subwagon_suffix == "" && !macro_filename.Contains("/")) {
    subwagon_suffix = macro_filename;
    macro_filename = DEFAULT_MACRO;
  }

  // build analysis name out of provided suffix
  const TString task_name = (subwagon_suffix == "")
                          ? default_name
                          : TString::Format("%s_%s", default_name, subwagon_suffix);

  cout << "[AddTaskPiLam]\n";
  cout << "   macro: '" << macro_filename << "'\n";
  cout << "   params: '" << params << "'\n";

  // The analysis config macro for PionLambdaFemto accepts a single string
  // argument, which it interprets.
  // This line escapes some escapable characters (backslash, newline, tab)
  // and wraps that string in double quotes, ensuring that the interpreter
  // reads a string when passing to the macro.
  const TString analysis_params = '"' + params.ReplaceAll("\\", "\\\\")
                                              .ReplaceAll("\n", "\\n")
                                              .ReplaceAll("\t", "\\t") + '"';

  AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto(
    task_name,
    macro_filename,
    analysis_params,
    kFALSE
  );

  mgr->AddTask(taskfemto);

  if (output_filename == "") {
    output_filename = AliAnalysisManager::GetCommonFileName();
  }

  const TString outputfile = TString::Format("%s:%s", output_filename.Data(), output_container.Data());

  AliAnalysisDataContainer *cout_femto = mgr->CreateContainer("femtolist",
                                                              TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              outputfile);

  // connect task to the containers
  mgr->ConnectInput(taskfemto, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskfemto, 0, cout_femto);

  // Return the task pointer
  return taskfemto;
}
