///
/// \file AddTaskFemtoSubtrains.C
///


#if !defined(__CINT__) && !defined(__CLING__)

#include "AliAnalysisTaskFemtoNu.h"

#include <TString.h>
#include <TProofMgr.h>
#include <TProof.h>

#endif


/// Create Femtoscopy AliAnalysisTask
///
AliAnalysisTaskFemto*
AddTaskFemtoNu(TString configMacroName,
               TString containerName,
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
    new AliAnalysisTaskFemtoNu("TaskFemto",
                               "$ALICE_PHYSICS/"+configMacroName,
                               configMacroParameters,
                               kFALSE);
  mgr->AddTask(task);

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG2FEMTO";

  task->SetupContainers(outputfile);

  return task;
}


AliAnalysisTaskFemto*
AddTaskFemtoNu(TString configMacroName,
               TString containerName,
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
