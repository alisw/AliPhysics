#include "TString.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskCVEPIDCMEDiff.h"

AliAnalysisTaskCVEPIDCMEDiff* AddTaskCVEPIDCMEDiff(
  TString trigger               = "kINT7+kSemiCentral",
  TString period                = "LHC18r",
  TString plane                 = "TPC",
  UInt_t filterBit              = 768,
  bool bCalculateLambdaHadron   = false,
  TString uniqueID              = "")
{
  // Creates a pid task and adds it to the analysis manager
  // Get the pointer to the existing analysis manager via the static
  // access method
  //=========================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCVEPIDCMEDiff.C", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskCVEPIDCMEDiff.C", "This task requires an input event handler.");
    return NULL;
  }

  bool bDoNUE                  = true;
  bool bDebug                  = false;
  bool bV0DaughterUseTOF       = false;
  bool bNarrowDcaCuts768       = true;
  bool bProtonCustomizedDCACut = true;
  bool bUsePionRejection       = false;
  // --- instantiate analysis task
  AliAnalysisTaskCVEPIDCMEDiff* task = new AliAnalysisTaskCVEPIDCMEDiff("TaskCVEPIDCMEDiff");
  task->SetTrigger(trigger);
  task->SetPeriod(period);
  task->SetFilterBit(filterBit);

  task->IfCalculateLambdaHadron(bCalculateLambdaHadron);
  task->IfDebug(bDebug);
  task->IfDoNUE(bDoNUE);
  task->IfV0DaughterUseTOF(bV0DaughterUseTOF);
  task->IfNarrowDcaCuts768(bNarrowDcaCuts768);
  task->IfProtonCustomizedDCACut(bProtonCustomizedDCACut);
  task->IfUsePionRejection(bUsePionRejection);
  task->SetPlaneEstimator(plane);

  //=========================================================================
  // Read in Files
  TFile* fNUEFile = nullptr;
  TFile* fNUAFile = nullptr;
  TFile* fVZEROCalibFile = nullptr;

  TList* fListNUE = nullptr;
  TList* fListNUA = nullptr;
  TList* fVZEROCalibList = nullptr;

  if (!gGrid) TGrid::Connect("alien://");
  if (bDoNUE) {
    if (period.EqualTo("LHC18q")) {
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/efficiency18q_FB768PL_dcazcut.root", "READ");
      fListNUE = dynamic_cast<TList*>(fNUEFile->Get("fListNUE"));
    }
    if (period.EqualTo("LHC18r")) {
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/efficiency18r_FB768PL_dcazcut.root", "READ");
      fListNUE = dynamic_cast<TList*>(fNUEFile->Get("fListNUE"));
    }
    if (fListNUE) {
      task->SetListForNUE(fListNUE);
      std::cout << "================  NUE List Set =================" << std::endl;
    } else
      std::cout << "!!!!!!!!!!!!!!!NUE List not Found!!!!!!!!!!!!!!!" << std::endl;
  }


  if (plane.EqualTo("TPC")) {
    TString fileNameNUA = "";
    if (period.EqualTo("LHC18q")) {
      if (uniqueID.EqualTo("Nhits60"))        fileNameNUA = TString("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_Nhits60.root");
      else if (uniqueID.EqualTo("Nhits80"))   fileNameNUA = TString("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_Nhits80.root");
      else if (uniqueID.EqualTo("ChiMax2"))   fileNameNUA = TString("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_ChiHg2.root");
      else if (uniqueID.EqualTo("ChiMax2p5")) fileNameNUA = TString("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_ChiHg2d5.root");
      else fileNameNUA = TString("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/WgtsNUAChargeAndPion_LHC18qPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root"); 
    } else if (period.EqualTo("LHC18r")) {
      if (uniqueID.EqualTo("Nhits60"))        fileNameNUA = TString("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_Nhits60.root");
      else if (uniqueID.EqualTo("Nhits80"))   fileNameNUA = TString("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_Nhits80.root");
      else if (uniqueID.EqualTo("ChiMax2"))   fileNameNUA = TString("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_ChiHg2.root");
      else if (uniqueID.EqualTo("ChiMax2p5")) fileNameNUA = TString("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_ChiHg2d5.root");
      else fileNameNUA = TString("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/WgtsNUAChargeAndPion_LHC18rPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root");
    }

    fNUAFile = TFile::Open(fileNameNUA, "READ");
    fListNUA = dynamic_cast<TList*>(fNUAFile->Get("fNUA_ChPosChNeg"));

    if (fListNUA) {
      task->SetListForNUA(fListNUA);
      std::cout << "================  NUA List Set =================" << std::endl;
    } else
      std::cout << "!!!!!!!!!!!!!!!NUA List not Found!!!!!!!!!!!!!!!" << std::endl;
  }

  if (plane.EqualTo("V0C")) {
    if (period.EqualTo("LHC18q")) {
      fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/calibq2V0C18qP3.root", "READ");
      fVZEROCalibList = dynamic_cast<TList*>(fVZEROCalibFile->Get("18qlistspPerc"));
    }
    if (period.EqualTo("LHC18r")) {
      fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/calibq2V0C18rP3.root", "READ");
      fVZEROCalibList = dynamic_cast<TList*>(fVZEROCalibFile->Get("18rlistspPerc"));
    }
    if (fVZEROCalibList) {
      task->SetListForVZEROCalib(fVZEROCalibList);
      std::cout << "================  VZERO List Set =================" << std::endl;
    } else
      std::cout << "!!!!!!!!!!!!!!!VZERO List not Found!!!!!!!!!!!!!!!" << std::endl;
  }

  //======================================================================
  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
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


