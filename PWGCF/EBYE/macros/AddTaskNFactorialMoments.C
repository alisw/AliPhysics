#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNFactorialMoments.h"

AliAnalysisTaskNFactorialMoments* AddTaskNFactorialMoments(Double_t pt1 = 0.4, Double_t pt2 = 0.6, Double_t pt3 = 0.4, Double_t pt4 = 1.0, Double_t pt5 = 0.6, Double_t pt6 = 0.8, Double_t pt7 = 0.6, Double_t pt8 = 2.0, Int_t Mmax = 82, Bool_t fIsMC = kFALSE, const char* suffix = "")
{

  // year = "2010" or "2015"
  // fSelfAffAnalysis for self-affine analysis
  // cuts on detas and dphis work only if fTwoTrack = kTRUE
  // cuts on DCA, ITS/TPC clusters, ncrssed rows work only if they have non-zero

  // Initialize variables
  TString year = "2015";
  Bool_t fSelfAffAnalysis = kFALSE;
  Bool_t f2TrackQA = kFALSE;
  Bool_t fPileup = kTRUE;
  Double_t fEtaMin = -0.8;
  Double_t fEtaMax = 0.8;
  Int_t fCentMin = 0;
  Int_t fCentMax = 5;
  Int_t fBit = 768;
  Double_t fDeta = 0.02;
  Double_t fDphi = 0.02;
  Bool_t fTwoTrack = kFALSE;
  Double_t fSharedFraction = 0.05;
  Bool_t fSharity = kFALSE;
  Double_t fVxMax = 0.3;
  Double_t fVyMax = 0.4;
  Double_t fVzMax = 10.0;
  Bool_t fRejectElectrons = kFALSE;
  Double_t fDCAXYRangeMax = 0.0;
  Double_t fDCAZRangeMax = 0.0;
  Double_t fITSClusterCut = 0.0;
  Double_t fTPCClusterCut = 0.0;
  Double_t fnTPCrossedrows = 0.0;
  Double_t fSharedCls = 0.0;
  Double_t fSharedRows = 0.0;
  Double_t fFindableCls = 0.0;
  Bool_t defSharedFrac = kTRUE;
  Int_t bField = 0;
  Double_t nSigmaCutPr = 3.0;

  Double_t ptarr[8] = { pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8 };
  TArrayD ptarrD(8, ptarr);
  Int_t fNumberofPtBins = ptarrD.GetSize() / 2;

  Int_t Mbinsarr[40];
  Int_t Nbinsarr[40]; // used in case of selfaffine analysis
  for (size_t nMB = 0; nMB < 40; nMB++) {
    if (Mmax == 123)
      Mbinsarr[nMB] = 3 * (nMB + 2);
    if (Mmax == 82)
      Mbinsarr[nMB] = 2 * (nMB + 2);
  }
  if (fSelfAffAnalysis) {
    for (Int_t nMB = 0; nMB < 40; nMB++) {
      Mbinsarr[nMB] = nMB + 1;
      Nbinsarr[nMB] = 4 * (nMB + 1);
    }
  }

  TArrayI Mbins(40, Mbinsarr);
  TArrayI Nbins(40, Nbinsarr);

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }

  // creating an instance for tasknfms
  AliAnalysisTaskNFactorialMoments* tasknfms =
    new AliAnalysisTaskNFactorialMoments(suffix);
  if (!tasknfms)
    return 0x0;

  tasknfms->SetIsMC(fIsMC, "Hijing");
  tasknfms->SetYear(year);
  tasknfms->SetPSbinning(Mbins, Nbins, Mmax);
  tasknfms->SetRejectPileup(fPileup);
  tasknfms->SetEta(fEtaMin, fEtaMax);
  tasknfms->SetCentLim(fCentMin, fCentMax);
  tasknfms->SetPtArray(ptarrD);
  tasknfms->Setfbit(fBit);
  tasknfms->SetBfield(bField);
  tasknfms->SetTwoTrackCuts(fDeta, fDphi, fTwoTrack, f2TrackQA);
  tasknfms->SetSharingFraction(fSharedFraction, fSharity);
  tasknfms->SetVtxCut(fVxMax, fVyMax, fVzMax);
  tasknfms->SetRejectElectrons(fRejectElectrons);
  tasknfms->SetDCAXYRangeMax(fDCAXYRangeMax);   // 0.1
  tasknfms->SetDCAZRangeMax(fDCAZRangeMax);     // 1
  tasknfms->SetITSClusterCut(fITSClusterCut);   // chi2 per ITS 36
  tasknfms->SetTPCClusterCut(fTPCClusterCut);   // chi2 per TPC 4
  tasknfms->SetnTPCrossedrows(fnTPCrossedrows); // 70
  tasknfms->SetSelfAffine(fSelfAffAnalysis);
  tasknfms->SetSharedCls(fSharedCls, fSharedRows, fFindableCls, defSharedFrac);

  // add your tasknfms to the manager
  mgr->AddTask(tasknfms);

  // input for tasknfmss and connecting containers to it
  Int_t nout = (fNumberofPtBins) + 3 + fIsMC * ((fNumberofPtBins)-1);
  AliAnalysisDataContainer* coutputlist[nout];

  coutputlist[0] = mgr->CreateContainer(
    Form("outputlist%s", suffix), TList::Class(), AliAnalysisManager::kOutputContainer,
    AliAnalysisManager::GetCommonFileName());

  char* name = "ntplist";
  if (fIsMC)
    name = "ntplistRe";
  for (Int_t nbins = 0; nbins < fNumberofPtBins; nbins++) {
    coutputlist[nbins + 1] =
      mgr->CreateContainer(Form("%s%d%s", name, nbins + 1, suffix), TList::Class(),
                           AliAnalysisManager::kOutputContainer,
                           AliAnalysisManager::GetCommonFileName());
    if (fIsMC)
      coutputlist[nbins + fNumberofPtBins + 1] =
        mgr->CreateContainer(Form("ntplistGen%d%s", nbins + 1, suffix), TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             AliAnalysisManager::GetCommonFileName());
  }
  if (!fIsMC)
    coutputlist[nout - 2] = mgr->CreateContainer(
      Form("QAlist%s", suffix), TList::Class(), AliAnalysisManager::kOutputContainer,
      AliAnalysisManager::GetCommonFileName());
  coutputlist[nout - 1] = mgr->CreateContainer(
    Form("QAlist2%s", suffix), TList::Class(), AliAnalysisManager::kOutputContainer,
    AliAnalysisManager::GetCommonFileName());

  mgr->ConnectInput(tasknfms, 0, mgr->GetCommonInputContainer());
  for (Int_t i = 0; i < nout; i++) {
    mgr->ConnectOutput(tasknfms, i + 1, coutputlist[i]);
  }
  return tasknfms;
}
