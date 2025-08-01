#include "TString.h"
#include "TGrid.h"
#include "TFile.h"
#include "TError.h"
#include "TList.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskCVEPIDCMEDiff.h"

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

AliAnalysisTaskCVEPIDCMEDiff* AddTaskCVEPIDCMEDiff(
    TString period   = "LHC18r",
    TString plane    = "TPC",
    TString pairs    = "Proton",
    TString uniqueID = ""
) {
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Fatal("AddTaskCVEPIDCMEDiff.C", "No analysis manager to connect to.");
        return nullptr;
    }
    if (!mgr->GetInputEventHandler()) {
        Fatal("AddTaskCVEPIDCMEDiff.C", "This task requires an input event handler.");
        return nullptr;
    }

    AliAnalysisTaskCVEPIDCMEDiff* task = new AliAnalysisTaskCVEPIDCMEDiff("TaskCVEPIDCMEDiff");
    task->SetPeriod(period);
    task->SetPlaneEstimator(plane);
    if (pairs.EqualTo("Proton"))        task->IfCalculateLambdaProton(true);
    else if (pairs.EqualTo("Hadron"))   task->IfCalculateLambdaHadron(true);
    else if (pairs.EqualTo("Pion"))     task->IfCalculateLambdaPion(true);
    else if (pairs.EqualTo("Lambda"))   task->IfCalculateLambdaLambda(true);
    else {
        Fatal("AddTaskCVEPIDCMEDiff.C", "Not support Lambda with this particle");
    }

    if (!gGrid) TGrid::Connect("alien://");

    // NUE NUA
    TList* fListNUENUA = nullptr;
    {
        TString nueFileName = "eff_pt_calib_cent.root";
        TFile* fNUEFile = TFile::Open(GetCalibFilePath(nueFileName), "READ");
        if (fNUEFile) {
            fListNUENUA = dynamic_cast<TList*>(fNUEFile->Get("fListNUENUA"));
            if (!fListNUENUA) Fatal("AddTaskCVEPIDCMEDiff.C", "NUENUA List not found!");
        } else {
            Fatal("AddTaskCVEPIDCMEDiff.C", "NUENUA File not found!");
        }
    }

    // Plane NUA
    TList* fListPlaneNUA = nullptr;
    if (plane.EqualTo("TPC")) {
        TString fileNameNUA{""};
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
        TFile* fNUAFile = TFile::Open(GetCalibFilePath(fileNameNUA), "READ");
        if (fNUAFile) {
            fListPlaneNUA = dynamic_cast<TList*>(fNUAFile->Get("fNUA_ChPosChNeg"));
            if (!fListPlaneNUA) Fatal("AddTaskCVEPIDCMEDiff.C", "Plane NUA List not found!");
        } else {
            Fatal("AddTaskCVEPIDCMEDiff.C", "Plane NUA File not found!");
        }
    } else {
        fListPlaneNUA = new TList();
        fListPlaneNUA->SetName("dummyPlaneNUAList");
    }

    // TPC
    TList* fListTPC = nullptr;
    if (plane.EqualTo("TPC")) {
        TString tpcFileName;
        if      (period.EqualTo("LHC18q")) tpcFileName = "TPCQVectorMean18q.root";
        else if (period.EqualTo("LHC18r")) tpcFileName = "TPCQVectorMean18r.root";
        else {
            Fatal("AddTaskCVEPIDCMEDiff.C", "Unsupported period for TPC calibration.");
        }
        TFile* fTPCFile = TFile::Open(GetCalibFilePath(tpcFileName), "READ");
        if (fTPCFile) {
            fListTPC = dynamic_cast<TList*>(fTPCFile->Get("fListTPC"));
            if (!fListTPC) Fatal("AddTaskCVEPIDCMEDiff.C", "TPC List not found!");
        } else {
            Fatal("AddTaskCVEPIDCMEDiff.C", "TPC File not found!");
        }
    } else if (plane.EqualTo("V0C")) {
        fListTPC = new TList();
        fListTPC->SetName("dummyTPCList");
        Info("AddTaskCVEPIDCMEDiff.C", "No need for TPC Plane calibration, dummy list created.");
    } else {
        Fatal("AddTaskCVEPIDCMEDiff.C", "Unsupported plane for TPC calibration.");
    }

    // VZERO
    TList* fListVZERO = nullptr;
    if (plane.EqualTo("V0C")) {
        TString vzeroFileName;
        if      (period.EqualTo("LHC18q")) vzeroFileName = "calibq2V0C18qP3.root";
        else if (period.EqualTo("LHC18r")) vzeroFileName = "calibq2V0C18rP3.root";
        else {
            Fatal("AddTaskCVEPIDCMEDiff.C", "Unsupported period for V0C calibration.");
        }
        TFile* fVZEROFile = TFile::Open(GetCalibFilePath(vzeroFileName), "READ");
        if (fVZEROFile) {
            fListVZERO = dynamic_cast<TList*>(fVZEROFile->Get(period.EqualTo("LHC18q") ? "18qlistspPerc" : "18rlistspPerc"));
            if (!fListVZERO) Fatal("AddTaskCVEPIDCMEDiff.C", "VZERO List not found!");
        } else {
            Fatal("AddTaskCVEPIDCMEDiff.C", "VZERO File not found!");
        }
    } else if (plane.EqualTo("TPC")) {
        fListVZERO = new TList();
        fListVZERO->SetName("dummyVZEROList");
        Info("AddTaskCVEPIDCMEDiff.C","No need for VZERO calibration, dummy list created.");
    } else {
        Fatal("AddTaskCVEPIDCMEDiff.C", "Unsupported plane!");
    }

    mgr->AddTask(task);

    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
    mgr->ConnectInput(task, 0, cinput);

    if (fListNUENUA) {
        AliAnalysisDataContainer* cinputNUE = mgr->CreateContainer(Form("NUEInput_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kInputContainer);
        cinputNUE->SetData(fListNUENUA);
        mgr->ConnectInput(task, 1, cinputNUE);
    } else {
        Fatal("AddTaskCVEPIDCMEDiff.C", "NUE List not found!");
    }

    if (fListPlaneNUA) {
        AliAnalysisDataContainer* cinputNUA = mgr->CreateContainer(Form("NUAInput_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kInputContainer);
        cinputNUA->SetData(fListPlaneNUA);
        mgr->ConnectInput(task, 2, cinputNUA);
    } else {
        Fatal("AddTaskCVEPIDCMEDiff.C", "NUA List not found!");
    }

    if (fListTPC) {
        AliAnalysisDataContainer* cinputTPC = mgr->CreateContainer(Form("TPCInput_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kInputContainer);
        cinputTPC->SetData(fListTPC);
        mgr->ConnectInput(task, 3, cinputTPC);
    } else {
        Fatal("AddTaskCVEPIDCMEDiff.C", "TPC List not found!");
    }

    if (fListVZERO) {
        AliAnalysisDataContainer* cinputVZERO = mgr->CreateContainer(Form("VZEROInput_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kInputContainer);
        cinputVZERO->SetData(fListVZERO);
        mgr->ConnectInput(task, 4, cinputVZERO);
    }

    // Output containers (QA, results)
    TString outputFileName = mgr->GetCommonFileName();
    AliAnalysisDataContainer* coutput_QA = mgr->CreateContainer(Form("ListQA_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFileName.Data(), uniqueID.Data()));
    AliAnalysisDataContainer* coutput_result = mgr->CreateContainer(Form("ListResults_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", outputFileName.Data(), uniqueID.Data()));

    mgr->ConnectOutput(task, 1, coutput_QA);
    mgr->ConnectOutput(task, 2, coutput_result);

    std::cout << ">>>>================  Return task =================<<<<" << std::endl;
    return task;
}
