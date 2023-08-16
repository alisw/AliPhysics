#include "TString.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskCVEPIDCME.h"

AliAnalysisTaskCVEPIDCME* AddTaskCVEPIDCME(
  TString trigger               = "kINT7+kSemiCentral",
  TString period                = "LHC18r",
  int filterBit                 = 768,
  bool bCalculatePIDFlow        = false,
  bool bCalculateDiffResult     = false,
  bool bCalculateDeltaPhiSumPhi = false,
  bool bCalculateLambdaProton   = false,
  bool bCalculateLambdaHadron   = false,
  bool bCalculateLambdaLambda   = false,
  bool bCalculateProtonProton   = false,
  bool bCalculateProtonHadron   = false,
  bool bCalculateHadronHadron   = false,
  TString uniqueID              = "")
{
  // Creates a pid task and adds it to the analysis manager
  // Get the pointer to the existing analysis manager via the static
  // access method
  //=========================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskCVEPIDCME.C", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskCVEPIDCME.C", "This task requires an input event handler.");
    return NULL;
  }

  bool bDoNUE = false;
  bool bDoLambdaNUE = false;
  if (uniqueID.EqualTo("Effi")) bDoNUE = true;
  if (uniqueID.EqualTo("EffiLambda")) bDoLambdaNUE = true;
  if (uniqueID.EqualTo("EffiAll")) bDoNUE = true, bDoLambdaNUE = true;

  bool bDebug                  = false;
  bool bUseTPCPlane            = true;
  bool bUseVZEROPlane          = true;
  bool bUseZDCPlane            = false;
  bool bDoNUA                  = true;
  bool bV0DaughterUseTOF       = false;
  bool bQATPC                  = true;
  bool bQAVZERO                = true;
  bool bQAZDC                  = false;
  bool bNarrowDcaCuts768       = (filterBit == 768);
  bool bProtonCustomizedDCACut = (filterBit == 768);
  bool bUsePionRejection       = false;
  bool bCalculateLambdaProtonFromDecay = false;
  bool bCalculateLambdaProtonFromDecayFoundInTrackLoops = false;
  bool bUseOneSideTPCPlane     = false;
  // --- instantiate analysis task
  AliAnalysisTaskCVEPIDCME* task = new AliAnalysisTaskCVEPIDCME("TaskCVEPIDCME");
  task->SetTrigger(trigger);
  task->SetPeriod(period);
  task->SetFilterBit(filterBit);

  task->IfDebug(bDebug);
  task->IfUseTPCPlane(bUseTPCPlane);
  task->IfUseVZEROPlane(bUseVZEROPlane);
  task->IfUseZDCPlane(bUseZDCPlane);
  task->IfDoNUE(bDoNUE);
  task->IfDoLambdaNUE(bDoLambdaNUE);
  task->IfDoNUA(bDoNUA);
  task->IfV0DaughterUseTOF(bV0DaughterUseTOF);
  task->IfQATPC(bQATPC);
  task->IfQAVZERO(bQAVZERO);
  task->IfQAZDC(bQAZDC);
  task->IfNarrowDcaCuts768(bNarrowDcaCuts768);
  task->IfProtonCustomizedDCACut(bProtonCustomizedDCACut);
  task->IfUsePionRejection(bUsePionRejection);

  task->IfCalculatePIDFlow(bCalculatePIDFlow);
  task->IfCalculateDiffResult(bCalculateDiffResult);
  task->IfCalculateDeltaPhiSumPhi(bCalculateDeltaPhiSumPhi);

  task->IfCalculateLambdaProton(bCalculateLambdaProton);
  task->IfCalculateLambdaHadron(bCalculateLambdaHadron);
  task->IfCalculateLambdaLambda(bCalculateLambdaLambda);
  task->IfCalculateProtonProton(bCalculateProtonProton);
  task->IfCalculateProtonHadron(bCalculateProtonHadron);
  task->IfCalculateHadronHadron(bCalculateHadronHadron);
  
  task->IfCheckLambdaProtonFromDecay(bCalculateLambdaProtonFromDecay);
  task->IfCheckLambdaProtonFromDecayFoundInTrackLoops(bCalculateLambdaProtonFromDecayFoundInTrackLoops);

  task->IfUseOneSideTPCPlane(bUseOneSideTPCPlane);


  //=========================================================================
  // Read in Files
  TFile* fNUEFile = nullptr;
  TFile* fLambdaNUEFile = nullptr;
  TFile* fNUAFile = nullptr;
  TFile* fVZEROCalibFile = nullptr;
  TFile* fZDCCalibFile = nullptr;

  TList* fListNUE = nullptr;
  TList* fListLambdaNUE = nullptr;
  TList* fListNUA = nullptr;
  TList* fVZEROCalibList = nullptr;
  TList* fZDCCalibList = nullptr;

  bool isNeedLoadAlienFile = bDoNUA || bDoNUE || bUseVZEROPlane || bUseZDCPlane;

  if (isNeedLoadAlienFile && !gGrid)
    TGrid::Connect("alien://");
  if (bDoNUE) {
    if (period.EqualTo("LHC10h")) {
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/Run1NUE.root", "READ");
      fListNUE = dynamic_cast<TList*>(fNUEFile->Get("listNUE"));
    }
    if (period.EqualTo("LHC15o")) {
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/efficiencyBothpol.root", "READ");
      fListNUE = dynamic_cast<TList*>(fNUEFile->Get("fMcEffiHij"));
    }
    if (period.EqualTo("LHC18q")) {
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/efficiency18q.root", "READ");
      fListNUE = dynamic_cast<TList*>(fNUEFile->Get("fListNUE"));
    }
    if (period.EqualTo("LHC18r")) {
      fNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/efficiency18r.root", "READ");
      fListNUE = dynamic_cast<TList*>(fNUEFile->Get("fListNUE"));
    }
    if (fListNUE) {
      task->SetListForNUE(fListNUE);
      std::cout << "================  NUE List Set =================" << std::endl;
    } else
      std::cout << "!!!!!!!!!!!!!!!NUE List not Found!!!!!!!!!!!!!!!" << std::endl;
  }

  if(bDoLambdaNUE) {
    fLambdaNUEFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/LambdaEff.root", "READ");
    fListLambdaNUE = dynamic_cast<TList*>(fLambdaNUEFile->Get("LambdaEfficiency"));
    if (fListLambdaNUE) {
      task->SetListForLambdaNUE(fListLambdaNUE);
      std::cout << "================  Lambda NUE List Set =================" << std::endl;
    } else
      std::cout << "!!!!!!!!!!!!!!!Lambda NUE List not Found!!!!!!!!!!!!!!!" << std::endl;
  }

  if (bDoNUA) {
    if (period.EqualTo("LHC10h")) {
      if (filterBit == 1)
        fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/10hNUAFB1.root", "READ");
      if (filterBit == 768)
        fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/10hNUAFB768.root", "READ");
      fListNUA = dynamic_cast<TList*>(fNUAFile->Get("10hListNUA"));
    }
    if (period.EqualTo("LHC15o")) {
      fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/wgtPion_NUAFB768DeftwPUcut_LHC15op2_24Aug2021.root", "READ");
      fListNUA = dynamic_cast<TList*>(fNUAFile->Get("15oListNUA"));
    }
    if (period.EqualTo("LHC18q")) {
      if (uniqueID.EqualTo("FB96")) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/WgtsNUAChargeAndPion_LHC18qPass3_FB96_AlexPU_DeftMode_Oct2021.root", "READ");
      else if (uniqueID.EqualTo("Nhits60")) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_Nhits60.root");
      else if (uniqueID.EqualTo("Nhits80")) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_Nhits80.root");
      else if (uniqueID.EqualTo("ChiMax2")) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_ChiHg2.root");
	    else if (uniqueID.EqualTo("ChiMax2p5")) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18q/LHC18q_pass3_NUA_ChiHg2d5.root");
      else fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/WgtsNUAChargeAndPion_LHC18qPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root", "READ");
      fListNUA = dynamic_cast<TList*>(fNUAFile->Get("fNUA_ChPosChNeg"));
    }
    if (period.EqualTo("LHC18r")) {
      if (uniqueID.EqualTo("FB96")) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/WgtsNUAChargeAndPion_LHC18rPass3_FB96_AlexPU_DeftMode_Oct2021.root", "READ");
      else if (uniqueID.EqualTo("Nhits60")) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_Nhits60.root");
      else if (uniqueID.EqualTo("Nhits80")) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_Nhits80.root");
      else if (uniqueID.EqualTo("ChiMax2")) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_ChiHg2.root");
	    else if (uniqueID.EqualTo("ChiMax2p5")) fNUAFile = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc18r/LHC18r_pass3_NUA_ChiHg2d5.root");
      else fNUAFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/WgtsNUAChargeAndPion_LHC18rPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root", "READ");
      fListNUA = dynamic_cast<TList*>(fNUAFile->Get("fNUA_ChPosChNeg"));
    }
    if (fListNUA) {
      task->SetListForNUA(fListNUA);
      std::cout << "================  NUA List Set =================" << std::endl;
    } else
      std::cout << "!!!!!!!!!!!!!!!NUA List not Found!!!!!!!!!!!!!!!" << std::endl;
  }

  if (bUseVZEROPlane) {
    if (period.EqualTo("LHC10h")) {
      fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/10hQnCalib.root", "READ");
      fVZEROCalibList = dynamic_cast<TList*>(fVZEROCalibFile->Get("10hlistqncalib"));
    }
    if (period.EqualTo("LHC15o")) {
      fVZEROCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/VZEROCalibFile15o.root", "READ");
      fVZEROCalibList = dynamic_cast<TList*>(fVZEROCalibFile->Get("VZEROCalibList"));
    }
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

  if (bUseZDCPlane) {
    if (period.EqualTo("LHC10h")) {
      fZDCCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC10h/ZDCCalibFile.root", "READ");
      fZDCCalibList = dynamic_cast<TList*>(fZDCCalibFile->Get("ZDCCalibList"));
    }
    // if (period.EqualTo("LHC15o")) {
    //   fZDCCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC15o/10hQnCalib.root","READ");
    //   fZDCCalibList = dynamic_cast <TList*> (fZDCCalibFile->Get("10hlistqncalib"));
    // }
    if (period.EqualTo("LHC18q")) {
      fZDCCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18q/RecenteringResultFinal_2018q.root", "READ");
      fZDCCalibList = dynamic_cast<TList*>(fZDCCalibFile->Get("fOutputRecenter"));
    }
    if (period.EqualTo("LHC18r")) {
      fZDCCalibFile = TFile::Open("alien:///alice/cern.ch/user/c/chunzhen/CalibFiles/LHC18r/RecenteringResultFinal_2018r.root", "READ");
      fZDCCalibList = dynamic_cast<TList*>(fZDCCalibFile->Get("fOutputRecenter"));
    }
    if (fZDCCalibList) {
      task->SetListForZDCCalib(fZDCCalibList);
      std::cout << "================  ZDC List Set =================" << std::endl;
    } else
      std::cout << "!!!!!!!!!!!!!!!ZDC List not Found!!!!!!!!!!!!!!!" << std::endl;
  }

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


