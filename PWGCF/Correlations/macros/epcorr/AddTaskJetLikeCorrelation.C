#include "AliAnalysisTaskJetLikeCorrelation.h"

AliAnalysisTaskJetLikeCorrelation *AddTaskJetLikeCorrelation(int collision, int filterbit, float fTwoTrackEffCut, float fConversionsCut, float fResonancesCut)
{
  
  double centarrPbPb[] = {0, 5, 10, 20, 50, 80};
  double centarrpp[] = {0, 100};
  double zvertarr[] = {-9, -7, -5, -3, -1, 1, 3, 5, 7, 9};
  double ptarr[] = {0.8, 1, 1.5, 2, 3, 4, 6, 8,15, 25 };
  
  TArrayD dcentarr;
  if (collision == 0) {
    dcentarr.Set(2, centarrpp);
  } else {
    dcentarr.Set(6, centarrPbPb);
  }

  TArrayD dzvertarr(10, zvertarr);
  TArrayD dptarr(10, ptarr);
  float fEtaCut = 0.9;
  float fPhiCut = TMath::TwoPi();
//  float fTwoTrackEffCut = 0.02;
//  float fConversionsCut = 0.04;
//  float fResonancesCut = 0.01;
  int fMixingPoolSize = 200000;
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    AliError("AddTaskRunQA No Analysis Manager to connect!");
    return NULL;

  }

  if (!mgr->GetInputEventHandler()) {
    AliError("AddTaskRunQA : No Input Event Handler!");
    return NULL;
  }

  AliAnalysisTaskJetLikeCorrelation *taskjetlikecorr = new AliAnalysisTaskJetLikeCorrelation("taskjetlikecorr");
  taskjetlikecorr->UseSeparateMixingPool(kFALSE);
  taskjetlikecorr->SetMC(kFALSE);
  taskjetlikecorr->SetEtaCut(fEtaCut);
  taskjetlikecorr->SetPhiCut(fPhiCut);
  taskjetlikecorr->SetMinPtTrigCut(3.000);
  taskjetlikecorr->SetTwoTrackEffCut(fTwoTrackEffCut);
  taskjetlikecorr->SetCentArray(dcentarr);
  taskjetlikecorr->SetPtArray(dptarr);
  taskjetlikecorr->SetZVertexArray(dzvertarr);
//  taskjetlikecorr->SetEventMixingQueueSize(fEventMixingQueueSize);
  taskjetlikecorr->SetMixingPoolSize(fMixingPoolSize);
  taskjetlikecorr->SetMinNumTrack(5000);
  taskjetlikecorr->SetTrackDepth(50000);
  taskjetlikecorr->SetCollision(collision);
  taskjetlikecorr->SetFilterBit(filterbit); //tpconly
//  taskjetlikecorr->SetFilterBit(768); //hybrid
  taskjetlikecorr->SetResonancesCut(fResonancesCut);
  taskjetlikecorr->SetConversionsCut(fConversionsCut);
  taskjetlikecorr->SetDebugOption(1);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
    AliAnalysisDataContainer *coutputList = mgr->CreateContainer("InclList", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
    AliAnalysisDataContainer *coutputListIn = mgr->CreateContainer("InList", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
    AliAnalysisDataContainer *coutputListOut = mgr->CreateContainer("OutList", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
//    AliAnalysisDataContainer *coutputListM1 = mgr->CreateContainer("M1List", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
//    AliAnalysisDataContainer *coutputListM2 = mgr->CreateContainer("M2List", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  
  mgr->AddTask(taskjetlikecorr);

  mgr->ConnectInput(taskjetlikecorr, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskjetlikecorr, 1, coutputList);
    mgr->ConnectOutput(taskjetlikecorr, 2, coutputListIn);
    mgr->ConnectOutput(taskjetlikecorr, 3, coutputListOut);
//    mgr->ConnectOutput(taskjetlikecorr, 4, coutputListM1);
//    mgr->ConnectOutput(taskjetlikecorr, 5, coutputListM2);
  
  return taskjetlikecorr;

}
