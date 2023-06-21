#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNFactorialMoments.h"

AliAnalysisTaskNFactorialMoments *AddTaskNFactorialMoments(
    TString year = "2015", bool fSelfAffAnalysis = kFALSE, double pt1 = 0.4,
    double pt2 = 0.6, double pt3 = 0.4, double pt4 = 1.0, double pt5 = 0.6,
    double pt6 = 0.8, double pt7 = 0.6, double pt8 = 0.9, bool fPileup = kFALSE,
    double fEtaMin = -0.8, double fEtaMax = 0.8, int fCentMin = 0,
    int fCentMax = 10, int fBit = 128, double fDeta = 0.02, double fDphi = 0.02,
    bool fTwoTrack = kFALSE, double fSharedFraction = 0.05,
    bool fSharity = kFALSE, double fVxMax = 0.3, double fVyMax = 0.4,
    double fVzMax = 10.0, bool fRejectElectrons = kFALSE,
    double fDCAXYRangeMax = 0.0, double fDCAZRangeMax = 0.0,
    double fITSClusterCut = 0.0, double fTPCClusterCut = 0.0,
    double fnTPCrossedrows = 0.0, double fSharedCls = 0.0, double fSharedRows = 0.0, double fFindableCls = 0.0, int Mmax = 123, bool fIsMC = kTRUE) {

  // year = "2010" or "2015"
  // fSelfAffAnalysis for self-affine analysis
  // cuts on detas and dphis work only if fTwoTrack = kTRUE
  // cuts on DCA, ITS/TPC clusters, ncrssed rows work only if they have non-zero
  // values Particle Indentification uses only TPC for now TODO: add TOF

  double ptarr[8] = {pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8};
  TArrayD ptarrD(8, ptarr);
  int fNumberofPtBins = ptarrD.GetSize() / 2;

  int Mbinsarr[40];
  int Nbinsarr[40]; // used in case of selfaffine analysis
  for (size_t nMB = 0; nMB < 40; nMB++) {
    if (Mmax == 123)
    Mbinsarr[nMB] = 3 * (nMB + 2);
    if (Mmax == 82)
    Mbinsarr[nMB] = 2 * (nMB + 2);
  }
  if (fSelfAffAnalysis) {
    for (int nMB = 0; nMB < 40; nMB++) {
      Mbinsarr[nMB] = nMB + 1;
      Nbinsarr[nMB] = 4 * (nMB + 1);
    }
  }

  TArrayI Mbins(40, Mbinsarr);
  TArrayI Nbins(40, Nbinsarr);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }

  // creating an instance for tasknfms
  AliAnalysisTaskNFactorialMoments *tasknfms =
      new AliAnalysisTaskNFactorialMoments("tasknfms");
  if (!tasknfms)
    return 0x0;

  tasknfms->SetIsMC(fIsMC, "Hijing");
  tasknfms->SetYear(year);
  tasknfms->SetPSbinning(Mbins, Nbins, Mmax);
  tasknfms->SetRejectPileup(kTRUE);
  tasknfms->SetEta(fEtaMin, fEtaMax);
  tasknfms->SetCentLim(fCentMin, fCentMax);
  tasknfms->SetPtArray(ptarrD);
  tasknfms->Setfbit(fBit);
  tasknfms->SetTwoTrackCuts(fDeta, fDphi, fTwoTrack);
  tasknfms->SetSharingFraction(fSharedFraction, fSharity);
  tasknfms->SetVtxCut(fVxMax, fVyMax, fVzMax);
  tasknfms->SetRejectElectrons(fRejectElectrons);
  tasknfms->SetDCAXYRangeMax(fDCAXYRangeMax);   // 0.1
  tasknfms->SetDCAZRangeMax(fDCAZRangeMax);     // 1
  tasknfms->SetITSClusterCut(fITSClusterCut);   // chi2 per ITS 36
  tasknfms->SetTPCClusterCut(fTPCClusterCut);   // chi2 per TPC 4
  tasknfms->SetnTPCrossedrows(fnTPCrossedrows); // 70
  tasknfms->SetSelfAffine(fSelfAffAnalysis);
  tasknfms->SetSharedCls(fSharedCls, fSharedRows, fFindableCls);

  // add your tasknfms to the manager
  mgr->AddTask(tasknfms);

  // input for tasknfmss and connecting containers to it
  int nout = (fNumberofPtBins) + 3 + fIsMC * ((fNumberofPtBins)-1);
  AliAnalysisDataContainer *coutputlist[nout];

  coutputlist[0] = mgr->CreateContainer(
      "outputlist", TList::Class(), AliAnalysisManager::kOutputContainer,
      AliAnalysisManager::GetCommonFileName());

  char *name = "ntplist";
  if (fIsMC)
    name = "ntplistRe";
  for (int nbins = 0; nbins < fNumberofPtBins; nbins++) {
    coutputlist[nbins + 1] =
        mgr->CreateContainer(Form("%s%d", name, nbins + 1), TList::Class(),
                             AliAnalysisManager::kOutputContainer,
                             AliAnalysisManager::GetCommonFileName());
    if (fIsMC)
      coutputlist[nbins + fNumberofPtBins + 1] =
          mgr->CreateContainer(Form("ntplistGen%d", nbins + 1), TList::Class(),
                               AliAnalysisManager::kOutputContainer,
                               AliAnalysisManager::GetCommonFileName());
  }
  if (!fIsMC)
    coutputlist[nout - 2] = mgr->CreateContainer(
        "QAlist", TList::Class(), AliAnalysisManager::kOutputContainer,
        AliAnalysisManager::GetCommonFileName());
  coutputlist[nout - 1] = mgr->CreateContainer(
      "QAlist2", TList::Class(), AliAnalysisManager::kOutputContainer,
      AliAnalysisManager::GetCommonFileName());

  mgr->ConnectInput(tasknfms, 0, mgr->GetCommonInputContainer());
  for (int i = 0; i < nout; i++) {
    mgr->ConnectOutput(tasknfms, i + 1, coutputlist[i]);
  }
  return tasknfms;
}
