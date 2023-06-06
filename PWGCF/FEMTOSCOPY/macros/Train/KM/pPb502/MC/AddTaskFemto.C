//May 2023, Konstantin.Mikhaylov@cern.ch
//AddTaskFemto.C to run MC train CF_pPb_MC
//~~~~~ local test passed ~~~~~~~~~~~~~~~~
#if !defined(__CINT__) && !defined(__CLING__)

#include <TROOT.h>
#include <TString.h>
#include <TList.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <AliAnalysisManager.h>
#include "AliAnalysisTaskFemto.h"

#endif

// Creates K+K- analysis task and adds it to the analysis manager.
AliAnalysisTaskFemto *AddTaskFemto(TString configMacroName, const char *containerName="femtolist", const char *configMacroParameters="" ) {

  //...Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFemto", "No analysis manager to connect to.");
    return NULL;
  }

  //...Check the analysis type using the event handlers connected to the analysis manager. The availability of MC handler cann also be checked here.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFemto", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  cout << "Found " <<type << " event handler" << endl;

  //...Create the task, add it to manager.
  if (TProofMgr::GetListOfManagers()->GetEntries()) {
    gProof->Load(configMacroName);
  }
  //local test:
  AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto("TaskFemto",configMacroName,configMacroParameters,kFALSE);
  //train:
  //  AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto("TaskFemto","$ALICE_PHYSICS/"+configMacroName,configMacroParameters,kFALSE);
  taskfemto->SelectCollisionCandidates(AliVEvent::kAny);//What does it mean?
  mgr->AddTask(taskfemto);

  //MC PID (not needed for data):
  //AliAnalysisTaskSE *pidresponse = AddTaskPIDResponse(kTRUE);

  //...Create ONLY the output containers for the data produced by the task. Get and connect other common input/output containers via the manager as below
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG2FEMTO";
  AliAnalysisDataContainer *cout_femto  = mgr->CreateContainer("Kpm_run2_MC",  TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   mgr->ConnectInput(taskfemto, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskfemto, 0, cout_femto);

   // Return task pointer at the end
   return taskfemto;
}
