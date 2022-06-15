#include "TString.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskLambdaProtonCVE.h"

AliAnalysisTaskLambdaProtonCVE* AddTaskLambdaProtonCVE(
  int               debug=0, // debug level controls amount of output statements
  TString   trigger="kINT7",
  TString   period="LHC18r",
  int         filterBit=768, // AOD filter bit selection
  bool       v0calibOn=true,
  bool      zdccalibOn=true,
  bool         QAVZERO=true,
  bool           QAZDC=true,
  bool           QATPC=true,
  bool          doNUE=false,
  bool           doNUA=true,
  bool    checkPIDFlow=true,
  TString        uniqueID=""
  )
{  
  // Creates a pid task and adds it to the analysis manager
  // Get the pointer to the existing analysis manager via the static
  // access method
  //=========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskLambdaProtonCVE.C", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskLambdaProtonCVE.C", "This task requires an input event handler.");
    return NULL;
  }
  //TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // --- instantiate analysis task
  AliAnalysisTaskLambdaProtonCVE *task = new AliAnalysisTaskLambdaProtonCVE("TaskLambdaProtonCVE");
  task->SetDebug(debug);
  task->SetTrigger(trigger);
  task->SetPeriod(period);
  task->SetFilterBit(filterBit);
  task->SetNUEOn(doNUE);
  task->SetNUAOn(doNUA);  
  task->IfVZEROCalibOn(v0calibOn);
  task->IfZDCCalibOn(zdccalibOn);
  task->IfQAVZERO(QAVZERO);
  task->IfQAZDC(QAZDC);
  task->IfCheckPIDFlow(checkPIDFlow);

  //=========================================================================
  //Read in Files
  TFile* fNUEFile = nullptr;
  TFile* fNUAFile = nullptr;
  TFile* fVZEROCalibFile = nullptr;
  TFile* fZDCCalibFile = nullptr;

  TList* fListNUE = nullptr;
  TList* fListNUA = nullptr;
  TList* fVZEROCalibList = nullptr;
  TList* fZDCCalibList = nullptr;

  if (!gGrid) TGrid::Connect("alien://");
  if (doNUE) {
    if (period.EqualTo("LHC10h")) {
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/Run1NUE.root","READ");
      fListNUE = dynamic_cast <TList*> (fNUEFile->Get("listNUE"));
    }
    if (period.EqualTo("LHC15o")) {
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/efficiencyBothpol.root","READ");
      fListNUE = dynamic_cast <TList*> (fNUEFile->Get("fMcEffiHij"));
    }
    if (period.EqualTo("LHC18q")) {
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/efficiencyBothpol18qnew.root","READ");
      fListNUE = dynamic_cast <TList*> (fNUEFile->Get("fMcEffiHij"));
    }
    if (period.EqualTo("LHC18r")) {
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/efficiencyBothpol18qnew.root","READ");
      fListNUE = dynamic_cast <TList*> (fNUEFile->Get("fMcEffiHij"));
    }
    if(fListNUE) {
      task->SetListForNUE(fListNUE);
      std::cout<<"================  NUE List Set ================="<<std::endl;
    } else std::cout<<"!!!!!!!!!!!!!!!NUE List not Found!!!!!!!!!!!!!!!"<<std::endl;
  }

  if (doNUA) {
    if (period.EqualTo("LHC10h")) {
      if(filterBit == 1)   fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/10hNUAFB1.root","READ");
      if(filterBit == 768) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/10hNUAFB768.root","READ");
      fListNUA = dynamic_cast <TList*> (fNUAFile->Get("10hListNUA"));
    }
    if (period.EqualTo("LHC15o")) {
      fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/wgtPion_NUAFB768DeftwPUcut_LHC15op2_24Aug2021.root","READ");
      fListNUA = dynamic_cast <TList*> (fNUAFile->Get("15oListNUA"));
    }
    if (period.EqualTo("LHC18q")) {
      fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/WgtsNUAChargeAndPion_LHC18qPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root","READ");
      fListNUA = dynamic_cast <TList*> (fNUAFile->Get("fNUA_ChPosChNeg"));
    }
    if (period.EqualTo("LHC18r")) {
      fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/WgtsNUAChargeAndPion_LHC18rPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root","READ");
      fListNUA = dynamic_cast <TList*> (fNUAFile->Get("fNUA_ChPosChNeg"));
    }
    if(fListNUA) {
      task->SetListForNUA(fListNUA);
      std::cout<<"================  NUA List Set ================="<<std::endl;
    } else std::cout<<"!!!!!!!!!!!!!!!NUA List not Found!!!!!!!!!!!!!!!"<<std::endl;
  }

  if (v0calibOn) {
    if (period.EqualTo("LHC10h")) {
      fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/10hQnCalib.root","READ");
      fVZEROCalibList = dynamic_cast <TList*> (fVZEROCalibFile->Get("10hlistqncalib"));
    }
    if (period.EqualTo("LHC15o")) {
      fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/VZEROCalibFile.root","READ");
      fVZEROCalibList = dynamic_cast <TList*> (fVZEROCalibFile->Get("VZEROCalibList"));
    }
    if (period.EqualTo("LHC18q")) {
      fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/calibSpq2V0C18qP3.root","READ");
      fVZEROCalibList = dynamic_cast <TList*> (fVZEROCalibFile->Get("fWgtsV0ZDC"));
    }
    if (period.EqualTo("LHC18r")) {
      fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/calibSpq2V0C18rP3.root","READ");
      fVZEROCalibList = dynamic_cast <TList*> (fVZEROCalibFile->Get("fWgtsV0ZDC"));
    }
    if(fVZEROCalibList) {
      task->SetListForVZEROCalib(fVZEROCalibList);
      std::cout<<"================  VZERO List Set ================="<<std::endl;
    } else std::cout<<"!!!!!!!!!!!!!!!VZERO List not Found!!!!!!!!!!!!!!!"<<std::endl;
  }

  if (zdccalibOn) {
    if (period.EqualTo("LHC10h")) {
      fZDCCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/ZDCCalibFile.root","READ");
      fZDCCalibList = dynamic_cast <TList*> (fZDCCalibFile->Get("ZDCCalibList"));
    }
    // if (period.EqualTo("LHC15o")) {
    //   fZDCCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/10hQnCalib.root","READ");
    //   fZDCCalibList = dynamic_cast <TList*> (fZDCCalibFile->Get("10hlistqncalib"));
    // }
    if (period.EqualTo("LHC18q")) {
      fZDCCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/RecenteringResultFinal_2018q.root","READ");
      fZDCCalibList = dynamic_cast <TList*> (fZDCCalibFile->Get("fOutputRecenter"));
    }
    if (period.EqualTo("LHC18r")) {
      fZDCCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/RecenteringResultFinal_2018r.root","READ");
      fZDCCalibList = dynamic_cast <TList*> (fZDCCalibFile->Get("fOutputRecenter"));
    }
    if(fZDCCalibList) {
      task->SetListForZDCCalib(fZDCCalibList);
      std::cout<<"================  ZDC List Set ================="<<std::endl;
    } else std::cout<<"!!!!!!!!!!!!!!!ZDC List not Found!!!!!!!!!!!!!!!"<<std::endl;
  }

  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
  //======================================================================
    mgr->AddTask(task);
    AliAnalysisDataContainer* cinput  = mgr->GetCommonInputContainer();
    TString outputFileName = mgr->GetCommonFileName();
    std::cout<<"outputfileName::::==========:::"<<outputFileName<<std::endl;
    AliAnalysisDataContainer* coutput = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer,
                                                             Form("%s:%s", outputFileName.Data(), uniqueID.Data()));
    mgr->ConnectInput (task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);

  //==============================================================
  // Return task pointer at the end

  std::cout<<"================  Return task ================="<<std::endl;
  return task;
}  

