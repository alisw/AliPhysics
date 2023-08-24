#include "TString.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskCVEVZEROCalib.h"

AliAnalysisTaskCVEVZEROCalib* AddTaskCVEVZEROCalib(
  TString trigger               = "kINT7+kCentral+kSemiCentral",
  TString period                = "LHC18q",
  TString uniqueID              = "")
{
  // Creates a pid task and adds it to the analysis manager
  // Get the pointer to the existing analysis manager via the static
  // access method
  //=========================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCVEVZEROCalib.C", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskCVEVZEROCalib.C", "This task requires an input event handler.");
    return NULL;
  }

  // --- instantiate analysis task
  AliAnalysisTaskCVEVZEROCalib* task = new AliAnalysisTaskCVEVZEROCalib("TaskCVEVZEROCalib");
  task->SetTrigger(trigger);
  task->SetPeriod(period);
  task->IfDebug(0);

  //=========================================================================
  // Read in Files
  TFile* fVZEROCalibFile = nullptr;
  TList* fVZEROCalibList = nullptr;

  //==============================Loading Calibration Files===========================================
  if (!gGrid) TGrid::Connect("alien://");

  if (period.EqualTo("LHC10h")) {
    fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/10hQnCalib.root", "READ");
    fVZEROCalibList = dynamic_cast<TList*>(fVZEROCalibFile->Get("10hlistqncalib"));
  }
  if (period.EqualTo("LHC15o")) {
    fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/VZEROCalibFile.root", "READ");
    fVZEROCalibList = dynamic_cast<TList*>(fVZEROCalibFile->Get("VZEROCalibList"));
  }
  if (period.EqualTo("LHC18q")) {
    fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/calibSpq2V0C18qP3.root", "READ");
    fVZEROCalibList = dynamic_cast<TList*>(fVZEROCalibFile->Get("fWgtsV0ZDC"));
  }
  if (period.EqualTo("LHC18r")) {
    fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/calibSpq2V0C18rP3.root", "READ");
    fVZEROCalibList = dynamic_cast<TList*>(fVZEROCalibFile->Get("fWgtsV0ZDC"));
  }
  if (fVZEROCalibList) {
    task->SetListForVZEROCalib(fVZEROCalibList);
    std::cout << "================  VZERO List Set =================" << std::endl;
  } else
    std::cout << "!!!!!!!!!!!!!!!VZERO List not Found!!!!!!!!!!!!!!!" << std::endl;

  //======================================================================
  mgr->AddTask(task);
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  TString outputFileName = mgr->GetCommonFileName();
  std::cout << "outputfileName::::==========:::" << outputFileName << std::endl;
  AliAnalysisDataContainer* coutput_QA = mgr->CreateContainer(Form("ListQA_%s", uniqueID.Data()), TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           Form("%s:%s", outputFileName.Data(), uniqueID.Data()));
  AliAnalysisDataContainer* coutput_result = mgr->CreateContainer(Form("ListResults_%s", uniqueID.Data()), TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           Form("%s:%s", outputFileName.Data(), uniqueID.Data()));
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput_QA);
  mgr->ConnectOutput(task, 2, coutput_result);

  //==============================================================
  // Return task pointer at the end

  std::cout << "================  Return task =================" << std::endl;
  return task;
}





