///
/// \file AddTaskFemtoCommands.C
///


#if !defined(__CINT__) && !defined(__CLING__)

#include "AliAnalysisTaskFemto.h"

#include <TString.h>
#include <TProofMgr.h>
#include <TProof.h>

#endif


/// Create Femtoscopy AliAnalysisTask
///
/// The commands are split on ';' and passed to `gROOT->ProcessLine()`
///
///
AliAnalysisTaskFemto*
AddTaskFemtoSubtrains(TString commands,
                      TString macro_params,
                      TString subwagon="")
{
  TString macro = "",
          container = "",
          output_container = "PWG2FEMTO",
          task_name = "TaskFemto";

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFemto", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFemto", "This task requires an input event handler");
    return NULL;
  }

  // loop thorugh and execute commands
  if (!commands.IsWhitespace()) {
    TObjArray *cmds = commands.Tokenize("\n;");
    TIter next_cmd(cmds);
    while (TObjString *cmd_obj = static_cast<TObjString*>(next_cmd())) {
      TString cmd = cmd_obj->String().Strip(TString::kBoth, ' ');
      if (cmd.IsWhitespace()) {
        continue;
      }

      cmd.ReplaceAll("'", '"');
      std::cout << "Running: `" << cmd << "`\n";
      gROOT->ProcessLineFast(cmd + ";");
    }

    delete cmds;
  }

  // Create the task, add it to manager.
  if (TProofMgr::GetListOfManagers()->GetEntries()) {
    gProof->Load(macro);
  }

  // forward subwagon identifier to the macro
  if (!subwagon.IsWhitespace()) {
    if (!macro_params.IsWhitespace()) {
      macro_params += ", ";
    }
    macro_params += "\"" + subwagon + "\"";
  }

  AliAnalysisTaskFemto *taskfemto =
    new AliAnalysisTaskFemto(task_name,
                             "$ALICE_PHYSICS/"+macro,
                             macro_params,
                             kFALSE);
  mgr->AddTask(taskfemto);

  // Get and connect other common input/output containers via the manager
  if (!subwagon.IsWhitespace()) {
    container += "_" + subwagon;
  }

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":" + output_container;

  AliAnalysisDataContainer *cout_femto
    = mgr->CreateContainer(container,
                           TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           outputfile);

   mgr->ConnectInput(taskfemto, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskfemto, 0, cout_femto);

   return taskfemto;
}
