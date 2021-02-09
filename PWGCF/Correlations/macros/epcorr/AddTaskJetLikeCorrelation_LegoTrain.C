AliAnalysisTaskJetLikeCorrelation *AddTaskJetLikeCorrelation_LegoTrain(int collision, int filterbit, float fTwoTrackEffCut, float fConversionsCut, float fResonancesCut, int fNumberOfPlanes, TString dataset, TString taskname="JetLikeCorrelation")
{

 taskname = Form("%s_Fbit%d_Twotrackeff%.0f_fConv%.0f_fReso_%.0f", dataset.Data(), filterbit, fTwoTrackEffCut*100, fConversionsCut*100 , fResonancesCut*100);
  
  double centarrpbpb[] = {0, 5, 10, 20, 50, 80};
  double centarrpp[] = {0, 100};
  double zvertarr[] = {-9, -7, -5, -3, -1, 1, 3, 5, 7, 9};
  double zvertarrpp[] = {-7, -3,3, 7};
  double ptaarr[] = {0.8, 1, 2, 3, 4, 6, 8,15};
  double pttarr[] = { 4, 6, 8,15 };
  TArrayD dcentarr;
  TArrayD dzvertarr;
  if (collision == 0) {
    dcentarr.Set(2, centarrpp);
    dzvertarr.Set(4, zvertarrpp);
  } else {
    dcentarr.Set(6, centarrpbpb);
    dzvertarr.Set(10, zvertarr);
  }
  TArrayD dpttarr(4, pttarr);
  TArrayD dptaarr(8, ptaarr);
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

  AliAnalysisTaskJetLikeCorrelation *taskjetlikecorr = new AliAnalysisTaskJetLikeCorrelation(taskname);
  taskjetlikecorr->UseSeparateMixingPool(kFALSE);
  taskjetlikecorr->SetMCCorrection(kFALSE);
  taskjetlikecorr->SetMCTruth(kFALSE);
  taskjetlikecorr->SetEtaCut(fEtaCut);
  taskjetlikecorr->SetPhiCut(fPhiCut);
  taskjetlikecorr->SetMinPtTrigCut(3.000);
  taskjetlikecorr->SetTwoTrackEffCut(fTwoTrackEffCut);
  taskjetlikecorr->SetCentArray(dcentarr);
  taskjetlikecorr->SetPttArray(dpttarr);
  taskjetlikecorr->SetPtaArray(dptaarr);
  taskjetlikecorr->SetMinimumPtABinForMerging(4);
  taskjetlikecorr->SetZVertexArray(dzvertarr);
  taskjetlikecorr->SetNumberOfPlanes(fNumberOfPlanes);
//  taskjetlikecorr->SetEventMixingQueueSize(fEventMixingQueueSize);
  taskjetlikecorr->SetMixingPoolSize(fMixingPoolSize);
  taskjetlikecorr->SetMinNumTrack(3000);
  taskjetlikecorr->SetTrackDepth(30000);   //50000 for single grid, 30000 for train?
  taskjetlikecorr->SetCollision(collision);
  taskjetlikecorr->SetFilterBit(filterbit); //tpconly
//  taskjetlikecorr->SetFilterBit(768); //hybrid
  taskjetlikecorr->SetResonancesCut(fResonancesCut);
  taskjetlikecorr->SetConversionsCut(fConversionsCut);
  taskjetlikecorr->SetDebugOption(0);
  
  mgr->AddTask(taskjetlikecorr);

  TString outputFileName = TString::Format("%s:%s", AliAnalysisManager::GetCommonFileName(), taskname.Data());
  
  AliAnalysisDataContainer *coutputList[7];
  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    string containername;
    coutputList[iplane] = mgr->CreateContainer(Form("List%02d", iplane), TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  }

  mgr->ConnectInput(taskjetlikecorr, 0, mgr->GetCommonInputContainer());
  for (int iplane = 0; iplane < fNumberOfPlanes; iplane++) {
    mgr->ConnectOutput(taskjetlikecorr, iplane+1, coutputList[iplane]);
  }

  return taskjetlikecorr;

}
