#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliAnalysisDataContainer.h"
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskTOFSpectra.h>
#include <AliMCEventHandler.h>
#include <fstream>
#include <iostream>
#endif

AliAnalysisTaskTOFSpectra* AddTaskTOFSpectra(const Bool_t optTree = kTRUE, const Bool_t readMC = kFALSE, const Int_t system = 0 /*kPbPb*/, const Bool_t ChannelMismatch = kTRUE, const Bool_t CutVariation = kFALSE, const Int_t SimpleCutMode = -1, const TString prefix = "", const TString tname = "TOFSpectra", const Bool_t stdoutput = kFALSE)
{
  // Creates, configures and attaches to the train the task for pi, K , p spectra
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================

  Info("AddTaskTOFSpectra", "Adding a new task %s with this settings optTree = %i, readMC = %i, system = %i, ChannelMismatch = %i, CutVariation = %i, SimpleCutMode = %i", tname.Data(), optTree, readMC, system, ChannelMismatch, CutVariation, SimpleCutMode);

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTOFSpectra", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskTOFSpectra", "This task requires an input event handler");
    return NULL;
  }

  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type.Contains("AOD")) {
    Error("AddTaskTOFSpectra", "This task requires to run on ESD");
    return NULL;
  }

  // Add MC handler (for kinematics)
  if (readMC) {
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  // Create and configure the task
  AliAnalysisTaskTOFSpectra* taskTOF = new AliAnalysisTaskTOFSpectra(tname, system, readMC, optTree, ChannelMismatch, CutVariation, SimpleCutMode);
  mgr->AddTask(taskTOF);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  //   TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //   outputFileName += ":PWGLFSpectraTOF";
  //   Info("AddTaskTOFSpectra","The results of this task will be found in: %s",outputFileName.Data());

  //Create and attach input
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(taskTOF, 0, cinput);

  //Create and attach output
  TString ListFileName = "TListTOF";
  if (readMC)
    ListFileName.Append("_MC");
  //
  if (ChannelMismatch)
    ListFileName.Append("_Mismatch");
  //
  if (stdoutput)
    ListFileName = "AnalysisResults";
  //
  ListFileName.Append(".root");
  AliAnalysisDataContainer* cOutput1 = mgr->CreateContainer(Form("cOutputList%s", prefix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, ListFileName);
  mgr->ConnectOutput(taskTOF, 1, cOutput1);

  AliAnalysisDataContainer* cOutput2 = 0x0;
  if (optTree) {
    TString TreeFileName = "TreeTOF";
    if (readMC)
      TreeFileName.Append("_MC");
    //
    TreeFileName.Append(".root");
    cOutput2 = mgr->CreateContainer("cOutputTree", TTree::Class(), AliAnalysisManager::kOutputContainer, TreeFileName.Data());
    cOutput2->SetSpecialOutput();
    mgr->ConnectOutput(taskTOF, 2, cOutput2);
  }

  return taskTOF;
}
