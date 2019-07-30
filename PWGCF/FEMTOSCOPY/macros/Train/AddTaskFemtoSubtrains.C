///
/// \file AddTaskFemtoSubtrains.C
///


#if !defined(__CINT__) && !defined(__CLING__)

#include "AliAnalysisTaskFemto.h"

#include <TString.h>
#include <TProofMgr.h>
#include <TProof.h>

#endif


/// Create Femtoscopy AliAnalysisTask
///
AliAnalysisTaskFemto*
AddTaskFemtoSubtrains(TString configMacroName,
                      TString containerName,
                      TString configMacroParameters,
                      TString subtrain="")
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

  // Create the task, add it to manager.
  if (TProofMgr::GetListOfManagers()->GetEntries()) {
    gProof->Load(configMacroName);
  }

  // forward subtrain identifier to the macro
  if (!subtrain.IsWhitespace()) {
    configMacroParameters += ", \"" + subtrain + "\"";
  }

  AliAnalysisTaskFemto *taskfemto =
    new AliAnalysisTaskFemto("TaskFemto",
                             "$ALICE_PHYSICS/"+configMacroName,
                             configMacroParameters,
                             kFALSE);
  mgr->AddTask(taskfemto);

  // Get and connect other common input/output containers via the manager
  if (!subtrain.IsWhitespace()) {
    containerName += "_" + subtrain;
  }

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG2FEMTO";

  AliAnalysisDataContainer *cout_femto
    = mgr->CreateContainer(containerName,
                           TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           outputfile);

   mgr->ConnectInput(taskfemto, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskfemto, 0, cout_femto);

   return taskfemto;
}
