#include "AliAnalysisTaskJetLikeCorrelation.h"

AliAnalysisTaskJetLikeCorrelation *AddTaskJetLikeCorrelation(int collision, int filterbit, float fTwoTrackEffCut, float fConversionsCut, float fResonancesCut, int fNumberOfPlanes)
{
  
  double centarrPbPb[] = {0, 5, 10, 20, 50, 80};
  double centarrpp[] = {0, 100};
  double zvertarr[] = {-9, -7, -5, -3, -1, 1, 3, 5, 7, 9};
  double zvertarrpp[] = {-7,-3,3, 7};
  double ptaarr[] = {0.8, 1, 2, 3, 4, 6, 8,15};
  double pttarr[] = {4, 6, 8,15};
  
  TArrayD dcentarr;
  TArrayD dzvertarr;
  if (collision == 0) {
    dcentarr.Set(2, centarrpp);
    dzvertarr.Set(4, zvertarr);
  } else {
    dcentarr.Set(6, centarrPbPb);
    dzvertarr.Set(10, zvertarr);
  }

  
  TArrayD dptaarr(8, ptaarr);
  TArrayD dpttarr(4, pttarr);
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
  taskjetlikecorr->SetMCCorrection(kFALSE);
  taskjetlikecorr->SetMCTruth(kFALSE);
  taskjetlikecorr->SetEtaCut(fEtaCut);
  taskjetlikecorr->SetPhiCut(fPhiCut);
  taskjetlikecorr->SetMinPtTrigCut(3.000);
  taskjetlikecorr->SetMinimumPtABinForMerging(5);
  taskjetlikecorr->SetTwoTrackEffCut(fTwoTrackEffCut);
  taskjetlikecorr->SetCentArray(dcentarr);
  taskjetlikecorr->SetPtaArray(dptaarr);
  taskjetlikecorr->SetPttArray(dpttarr);
  taskjetlikecorr->SetZVertexArray(dzvertarr);
  taskjetlikecorr->SetNumberOfPlanes(fNumberOfPlanes);
//  taskjetlikecorr->SetEventMixingQueueSize(fEventMixingQueueSize);
  taskjetlikecorr->SetMixingPoolSize(fMixingPoolSize);
  taskjetlikecorr->SetMinNumTrack(5000);
  taskjetlikecorr->SetTrackDepth(50000);
  taskjetlikecorr->SetCollision(collision);
  taskjetlikecorr->SetFilterBit(filterbit); //tpconly
//  taskjetlikecorr->SetFilterBit(768); //hybrid
  taskjetlikecorr->SetResonancesCut(fResonancesCut);
  taskjetlikecorr->SetConversionsCut(fConversionsCut);
  taskjetlikecorr->SetDebugOption(0);

  TString outputFileName = AliAnalysisManager::GetCommonFileName();
    AliAnalysisDataContainer *coutputList[7];
    for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
      string containername;
      coutputList[iplane] = mgr->CreateContainer(Form("List%02d", iplane), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
    }
  
  mgr->AddTask(taskjetlikecorr);

  mgr->ConnectInput(taskjetlikecorr, 0, mgr->GetCommonInputContainer());
  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    mgr->ConnectOutput(taskjetlikecorr, iplane+1, coutputList[iplane]);
  }
  
  return taskjetlikecorr;

}
