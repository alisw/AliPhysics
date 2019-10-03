///
/// \file AddTaskFemtoNu.C
///


#if !defined(__CINT__) && !defined(__CLING__)

#include "AliAnalysisTaskFemtoNu.h"

#include <TString.h>
#include <TProofMgr.h>
#include <TProof.h>

#endif


/// Create a femtoscopic analysis task with results stored in an
/// AliFemtoResultStorage object instead of TList for more efficient
/// file merging.
///
/// Behavior is the same as standard AddTaskFemto:
///
/// - Output file is the manager's CommonFileName
/// - Results are stored in the `/PWG2FEMTO` TDirectory
/// - Config macro path is relative to `$ALICE_PHYSICS`
///
AliAnalysisTaskFemto*
AddTaskFemtoNu(TString containerName,
               TString configMacroName,
               TString configMacroParameters)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFemto", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFemto", "This task requires an input event handler");
    return NULL;
  }

  AliAnalysisTaskFemtoNu *task =
    new AliAnalysisTaskFemtoNu(containerName,
                               "$ALICE_PHYSICS/"+configMacroName,
                               configMacroParameters,
                               kFALSE);
  mgr->AddTask(task);

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG2FEMTO";

  task->SetupContainers(outputfile);

  return task;
}


/// Handle creating trains with subwagons:
///
/// - The subwagon name is joined to the containername with a '_'
/// - The subwagon name is appended to the configMacroParameters
///   as the last parameter, moving subwagon handling to the config
///   macro
///
AliAnalysisTaskFemto*
AddTaskFemtoNu(TString containerName,
               TString configMacroName,
               TString configMacroParameters,
               TString subwagon)
{
  containerName += "_" + subwagon;

  // add subwagon to list of config-macro parameters
  if (!configMacroParameters.IsWhitespace()) {
    configMacroParameters += ", ";
  }
  configMacroParameters += "\"" + subwagon + "\"";

  return AddTaskFemtoNu(configMacroName, containerName, configMacroParameters);
}
