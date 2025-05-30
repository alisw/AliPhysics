#include "TString.h"
#include "TGrid.h"
#include "TFile.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskCVEPIDCMEDiff.h"
#include "TFile.h"
#include <TError.h>

static bool isLocalTest = 0;

TString GetCalibFilePath(const TString& fileName) {
    TString basePath;
    if (isLocalTest) {
        basePath = "file:./calibration_files/";
    } else {
        basePath = "alien:///alice/cern.ch/user/c/chunzhen/calibration_files/";
    }
    return basePath + fileName;
}

AliAnalysisTaskCVEPIDCMEDiff* AddTaskCVEPIDCMEDiff(
  TString period                = "LHC18r",
  TString plane                 = "TPC",
  TString pairs                 = "Proton",
  TString uniqueID              = "")
{
  // Creates a pid task and adds it to the analysis manager
  // Get the pointer to the existing analysis manager via the static
  // access method
  //=========================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Fatal("AddTaskCVEPIDCMEDiff.C", "No analysis manager to connect to.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Fatal("AddTaskCVEPIDCMEDiff.C", "This task requires an input event handler.");
    return nullptr;
  }

  // --- instantiate analysis task
  AliAnalysisTaskCVEPIDCMEDiff* task = new AliAnalysisTaskCVEPIDCMEDiff("TaskCVEPIDCMEDiff");
  task->SetPeriod(period);
  task->SetPlaneEstimator(plane);
  if (pairs.EqualTo("Proton")) {
      task->IfCalculateLambdaProton(true);
  } else if (pairs.EqualTo("Hadron")) {
      task->IfCalculateLambdaHadron(true);
  } else if (pairs.EqualTo("Pion")) {
      task->IfCalculateLambdaPion(true);
  } else if (pairs.EqualTo("Lambda")) {
     task->IfCalculateLambdaLambda(true);
  } else {
     Fatal("AddTaskCVEPIDCMEDiff.C", "Not support Lambda with this particle");
  }

  //=========================================================================
  // Read in Files
  TFile* fNUEFile = nullptr;
  TFile* fNUAFile = nullptr;
  TFile* fVZEROCalibFile = nullptr;

  TList* fListNUE = nullptr;
  TList* fListNUA = nullptr;
  TList* fVZEROCalibList = nullptr;

  if (!gGrid) TGrid::Connect("alien://");

  // NUE
  if (period.EqualTo("LHC18q") || period.EqualTo("LHC18r")) {
    TString nueFileName = "eff_pt_calib_cent.root";

    fNUEFile = TFile::Open(GetCalibFilePath(nueFileName), "READ");
    if (fNUEFile) {
      fListNUE = dynamic_cast<TList*>(fNUEFile->Get("fListNUE"));
      if (fListNUE) {
        task->SetListForNUE(fListNUE);
        std::cout << ">>>>================  NUE List Set =================<<<<" << std::endl;
      } else {
        Fatal("AddTaskCVEPIDCMEDiff.C", "!!!!!!!!!!!!!!!NUE List not Found!!!!!!!!!!!!!!!");
      }
    } else {
      Fatal("AddTaskCVEPIDCMEDiff.C", "!!!!!!!!!!!!!!!NUE File not Found!!!!!!!!!!!!!!!");
    }
  }

  // NUA
  if (plane.EqualTo("TPC")) {
    TString fileNameNUA = "";
    if (period.EqualTo("LHC18q")) {
      if (uniqueID.EqualTo("Nhits60"))        fileNameNUA = "LHC18q_pass3_NUA_Nhits60.root";
      else if (uniqueID.EqualTo("Nhits80"))   fileNameNUA = "LHC18q_pass3_NUA_Nhits80.root";
      else if (uniqueID.EqualTo("ChiMax2"))   fileNameNUA = "LHC18q_pass3_NUA_ChiHg2.root";
      else if (uniqueID.EqualTo("ChiMax2p5")) fileNameNUA = "LHC18q_pass3_NUA_ChiHg2d5.root";
      else fileNameNUA = "WgtsNUAChargeAndPion_LHC18qPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root";
    } else if (period.EqualTo("LHC18r")) {
      if (uniqueID.EqualTo("Nhits60"))        fileNameNUA = "LHC18r_pass3_NUA_Nhits60.root";
      else if (uniqueID.EqualTo("Nhits80"))   fileNameNUA = "LHC18r_pass3_NUA_Nhits80.root";
      else if (uniqueID.EqualTo("ChiMax2"))   fileNameNUA = "LHC18r_pass3_NUA_ChiHg2.root";
      else if (uniqueID.EqualTo("ChiMax2p5")) fileNameNUA = "LHC18r_pass3_NUA_ChiHg2d5.root";
      else fileNameNUA = "WgtsNUAChargeAndPion_LHC18rPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root";
    }

    fNUAFile = TFile::Open(GetCalibFilePath(fileNameNUA), "READ");
    if (fNUAFile) {
      fListNUA = dynamic_cast<TList*>(fNUAFile->Get("fNUA_ChPosChNeg"));
      if (fListNUA) {
        task->SetListForNUA(fListNUA);
        std::cout << ">>>>================  NUA List Set =================<<<<" << std::endl;
      } else {
        Fatal("AddTaskCVEPIDCMEDiff.C", "!!!!!!!!!!!!!!!NUA List not Found!!!!!!!!!!!!!!!");
      }
    } else {
      Fatal("AddTaskCVEPIDCMEDiff.C", "!!!!!!!!!!!!!!!NUA File not Found!!!!!!!!!!!!!!!");
    }
  }

  // VZERO
  if (plane.EqualTo("V0C")) {
    TString vzeroFileName;
    if      (period.EqualTo("LHC18q")) {
      vzeroFileName = "calibq2V0C18qP3.root";
    }
    else if (period.EqualTo("LHC18r")) {
      vzeroFileName = "calibq2V0C18rP3.root";
    } else {
      std::cout<<Form("Unsupported period for V0C calibration: %s", period.Data())<<std::endl;
      return 0;
    }

    fVZEROCalibFile = TFile::Open(GetCalibFilePath(vzeroFileName), "READ");
    if (fVZEROCalibFile) {
      fVZEROCalibList = dynamic_cast<TList*>(fVZEROCalibFile->Get(period.EqualTo("LHC18q") ? "18qlistspPerc" : "18rlistspPerc"));
      if (fVZEROCalibList) {
        task->SetListForVZEROCalib(fVZEROCalibList);
        std::cout << ">>>>================  VZERO List Set =================<<<<" << std::endl;
      }
    } else {
      Fatal("AddTaskCVEPIDCMEDiff.C", "!!!!!!!!!!!!!!!VZERO List not Found!!!!!!!!!!!!!!!");
    }
  }

  //======================================================================
  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
  //======================================================================
  mgr->AddTask(task);
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  TString outputFileName = mgr->GetCommonFileName();
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

  std::cout << ">>>>================  Return task =================<<<<" << std::endl;
  return task;
}
