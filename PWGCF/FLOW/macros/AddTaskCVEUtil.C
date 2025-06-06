#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCVEUtil.h"
#include "TError.h"
#include "TString.h"

static bool isLocalTest = false;

TString GetCalibFilePath(const TString& fileName) {
    TString basePath;
    if (isLocalTest) {
        basePath = "file:./calibration_files/";
    } else {
        basePath = "alien:///alice/cern.ch/user/c/chunzhen/calibration_files/";
    }
    return basePath + fileName;
}

AliAnalysisTaskCVEUtil* AddTaskCVEUtil(
    bool isMC = false,
    TString period = "LHC18q",
    TString uniqueID = "default"
)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Fatal("AddTaskCVEUtil.C", "No analysis manager to connect to.");
        return nullptr;
    }
    if (!mgr->GetInputEventHandler()) {
        Fatal("AddTaskCVEUtil.C", "No input event handler.");
        return nullptr;
    }

    AliAnalysisTaskCVEUtil* task = new AliAnalysisTaskCVEUtil("TaskCVEUtil");
    if(!task) {
      Fatal("AddTaskCVEUtil.C", "Failed to create task.");
      return nullptr;
    }
    task -> SetMC(isMC);
    task -> SetPeriod(period);

    if (!gGrid) TGrid::Connect("alien://");

    // NUE
    TList* fListNUE = nullptr;
    if (!isMC) {
      TString nueFileName = "eff_pt_calib_cent.root";
      TFile* fNUEFile = TFile::Open(GetCalibFilePath(nueFileName), "READ");
      if (fNUEFile) {
          fListNUE = dynamic_cast<TList*>(fNUEFile->Get("fListNUE"));
          if (!fListNUE) Fatal("AddTaskCVEPIDCMEDiff.C", "NUE List not found!");
      } else {
          Fatal("AddTaskCVEPIDCMEDiff.C", "NUE File not found!");
      }
    } else {
        fListNUE = new TList();
        fListNUE->SetName("dummyMCNUE");
    }

    // NUA
    TList* fListNUA = nullptr;
    if (!isMC) {
      TString fileNameNUA{""};
      if (period.EqualTo("LHC18q")) {
          fileNameNUA = "WgtsNUAChargeAndPion_LHC18qPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root";
      } else if (period.EqualTo("LHC18r")) {
          fileNameNUA = "WgtsNUAChargeAndPion_LHC18rPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root";
      }
      TFile* fNUAFile = TFile::Open(GetCalibFilePath(fileNameNUA), "READ");
      if (fNUAFile) {
          fListNUA = dynamic_cast<TList*>(fNUAFile->Get("fNUA_ChPosChNeg"));
          if (!fListNUA) Fatal("AddTaskCVEPIDCMEDiff.C", "NUA List not found!");
      } else {
          Fatal("AddTaskCVEPIDCMEDiff.C", "NUA File not found!");
      }
    } else {
        fListNUA = new TList();
        fListNUA->SetName("dummyMCNUA");
    }

    mgr->AddTask(task);

    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
    mgr->ConnectInput(task, 0, cinput);

    if (fListNUE) {
        AliAnalysisDataContainer* cinputNUE = mgr->CreateContainer(Form("NUEInput_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kInputContainer);
        cinputNUE->SetData(fListNUE);
        mgr->ConnectInput(task, 1, cinputNUE);
    } else {
        Fatal("AddTaskCVEPIDCMEDiff.C", "NUE List not found!");
    }


    if (fListNUA) {
        AliAnalysisDataContainer* cinputNUA = mgr->CreateContainer(Form("NUAInput_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kInputContainer);
        cinputNUA->SetData(fListNUA);
        mgr->ConnectInput(task, 2, cinputNUA);
    } else {
        Fatal("AddTaskCVEPIDCMEDiff.C", "NUA List not found!");
    }

    TString outputFileName = AliAnalysisManager::GetCommonFileName();
    AliAnalysisDataContainer* outputContainer = mgr->CreateContainer(Form("ResultsList_%s", uniqueID.Data()), TList::Class(),
                                                    AliAnalysisManager::kOutputContainer,
                                                    Form("%s:%s", outputFileName.Data(), uniqueID.Data()));
    if (!outputContainer) {
        Fatal("AddTaskCVEUtil.C", "Failed to create output container.");
        return nullptr;
    }
    mgr->ConnectOutput(task,1,outputContainer);
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
