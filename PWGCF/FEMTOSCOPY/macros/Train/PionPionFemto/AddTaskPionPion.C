///
/// \file PWGCF/FEMTOSCOPY/macros/Train/PionPionFemto/AddTaskPionPion.C
/// \author Andrew Kubera, Ohio State University, andrew.kubera@cern.ch
///

///
/// \brief Adds an AliAnalysisTaskFemto analysis object, constructed
///        with parameters provided, to the global AliAnalysisManager.
///
/// This macro creates and returns an AliAnalysisTaskFemto object.
/// The task object is given a macro
/// is given the config macro "Train/PionPionFemto/ConfigFemtoAnalysis.C"
/// which is run to create the analysis objects. This is fixed (for now), and
/// if an alternative is required, you should use the general AddTaskFemto.C
/// macro.
///
/// Subwagons are supported, the name of which will be appended to the task
/// name. The task name by default is "TaskPionPionFemto" which cannot be changed.
///
/// The output container name is fixed at the standard "femtolist".
///
/// \param params
///     A string forwarded to the ConfigFemtoAnalysis macro for parsing.
///     This string is wrapped in double-quotes so escaping some to
///     ensure results is a string is unneccessary.
///
/// \param macro_path
///     The path to the femto-configuration macro, passed to the created
///     TaskFemto.
///     The string "%%" will be replaced by the location of this macro,
///     allowing short references to configuration files in this directory.
///     To set the name of the output_container, suffix this with a colon
///     followed by the name.
///     By default, the output_container name is PWG2FEMTO
///     If this parameter has no '/' characters, and the subwagon_suffix
///     is empty, this is reinterpreted as the subwagon_suffix.
///
/// \param output_filename
///     String with name of the output filename.
///     This string must contain ".root" or it will be interpreted as a
///     the subwagon_suffix.
///     If empty, this uses the default value returned by
///     AliAnalysisManager::GetCommonFileName().
///
/// \param subwagon_suffix
///     If this macro is run in a train with subwagons, this will be
///     set with the identifier.
///     NOTE: This parameter may be found in the macro_path
///           variable as the suffix is simply the last argument passed.
///           To keep from always having to explicitly set the
///           macro_path if using the default, this code will switch
///           the two parameters if there is no '/' characters in the
///           macro_path.
///           As such, keep '/' out of the subwagon_suffix!
///
AliAnalysisTaskFemto* AddTaskPionPion(TString container_name,
                                      TString params,
                                      TString macro_path="",
                                      TString output_filename="",
                                      TString subwagon_suffix="")
{ // Adds a Pion-Pion Femtoscopy task to the manager

  const TString AUTO_DIRECTORY = "$ALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/Train/PionPionFemto/"
              , DEFAULT_MACRO = "%%/ConfigFemtoAnalysis.C"
              , DEFAULT_CONTAINER_NAME = "pionpion_femtolist"
              , DEFAULT_TASK_NAME = "TaskPionPion"
              ;

  // Get the global manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPionPion", "Could not get the global AliAnalysisManager.");
    return NULL;
  }

  // If the macro_path was set, use it, else use the default
  if (macro_path == "") {
    macro_path = DEFAULT_MACRO;
  }

  // Replace %% with this directory for convenience
  macro_path.ReplaceAll("%%", AUTO_DIRECTORY);

  // If subwagon_suffix was not set and there are no '/' characters in the
  // macro's path, interpret path as the subwagon_suffix and set default path
  if (subwagon_suffix == "" && output_filename == "" && !macro_path.Contains("/")) {
    subwagon_suffix = macro_path;
    macro_path = DEFAULT_MACRO;
  }
  // subwagon_suffix not specified, and output_filename does not contain
  // ".root" - switch the two!
  else if (subwagon_suffix == "" && output_filename != "" && !output_filename.Contains(".root")) {
    subwagon_suffix = output_filename;
    output_filename = "";
  }

  // build analysis name out of provided suffix
  const TString task_name = (subwagon_suffix == "")
                          ? DEFAULT_TASK_NAME
                          : TString::Format("%s_%s", DEFAULT_TASK_NAME, subwagon_suffix);

  cout << "[AddTaskPionPion]\n"
          "   macro: '" << macro_path << "'\n"
          "   params: '" << params << "'\n";

  // The analysis config macro for PionPionFemto accepts a single string
  // argument, which it interprets.
  // This line escapes some escapable characters (backslash, newline, tab)
  // and wraps that string in double quotes, ensuring that the interpreter
  // reads a string when passing to the macro.
  const TString analysis_params = '"' + params.ReplaceAll("\\", "\\\\")
                                              .ReplaceAll("\n", "\\n")
                                              .ReplaceAll("\t", "\\t") + '"';

  AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto(
    task_name,
    macro_path,
    analysis_params,
    kFALSE
  );

  mgr->AddTask(taskfemto);

  if (output_filename == "") {
    output_filename = AliAnalysisManager::GetCommonFileName();
  }

  const TString outputfile = output_filename.Contains(":")
                           ? output_filename
                           : output_filename + ":PWG2FEMTO";

  container_name = (container_name == "")
                 ? DEFAULT_CONTAINER_NAME
                 : container_name;

  AliAnalysisDataContainer *out_container = mgr->CreateContainer(container_name,
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputfile);

  // connect task to the containers
  mgr->ConnectInput(taskfemto, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskfemto, 0, out_container);

  // Return the task pointer
  return taskfemto;
}
