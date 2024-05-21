#include "AliAnalysisManager.h"
#include "AliAnalysisTaskNFactorialMomentsPID.h"

AliAnalysisTaskNFactorialMomentsPID* AddTaskNFactorialMomentsPID(Int_t nPt = 4, Double_t pt1 = 0.2, Double_t pt2 = 2.0, Double_t pt3 = 0.4, Double_t pt4 = 1.0, Double_t pt5 = 0.4, Double_t pt6 = 2.0, Double_t pt7 = 0.4, Double_t pt8 = 5.0, Double_t pt9 = 0.0, Double_t pt10 = 0.0, Int_t Mmax = 82, Bool_t fIsMC = kFALSE, const char* suffix = "")
{

  // year = "2010" or "2015"
  // Initialize variables
  TString year = "2015";
  Bool_t fPileup = kTRUE;
  Double_t fEtaMin = -0.8;
  Double_t fEtaMax = 0.8;
  Int_t fCentMin = 0;
  Int_t fCentMax = 5;
  Int_t fBit = 768;
  Double_t fVxMax = 0.3;
  Double_t fVyMax = 0.4;
  Double_t fVzMax = 10.0;
  Double_t fDCAXYRangeMax = 0.0;
  Double_t fDCAZRangeMax = 0.0;
  Double_t fITSClusterCut = 0.0;
  Double_t fTPCClusterCut = 0.0;
  Double_t fnTPCrossedrows = 0.0;
  Double_t fSharedCls = 0.0;
  Double_t fSharedRows = 0.0;
  Double_t fFindableCls = 0.0;
  Int_t bField = 0;
  Bool_t fSelfAffAnalysis = kFALSE;
  Double_t nSigmaCutPr = 3.0;

  std::vector<Double_t> ptVector;
  if (pt1 != 0.0)
    ptVector.push_back(pt1);
  if (pt2 != 0.0)
    ptVector.push_back(pt2);
  if (pt3 != 0.0)
    ptVector.push_back(pt3);
  if (pt4 != 0.0)
    ptVector.push_back(pt4);
  if (pt5 != 0.0)
    ptVector.push_back(pt5);
  if (pt6 != 0.0)
    ptVector.push_back(pt6);
  if (pt7 != 0.0)
    ptVector.push_back(pt7);
  if (pt8 != 0.0)
    ptVector.push_back(pt8);
  if (pt9 != 0.0)
    ptVector.push_back(pt9);
  if (pt10 != 0.0)
    ptVector.push_back(pt10);

  Int_t fNumberofPtBins = ptVector.size() / 2;
  std::cout << "Number of pt bins: " << fNumberofPtBins << std::endl;

  TGrid::Connect("alien://");

  if (fNumberofPtBins != nPt) {
    std::cout << "Number of pt bins is not equal to the number of pt values" << std::endl;
    return 0x0;
  }

  Double_t ptarr[nPt * 2];
  std::copy(ptVector.begin(), ptVector.end(), ptarr);
  TArrayD ptarrD(nPt * 2, ptarr);

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

  // creating an instance for taskFactorialMoments
  AliAnalysisTaskNFactorialMomentsPID* taskFactorialMoments =
    new AliAnalysisTaskNFactorialMomentsPID(suffix);
  if (!taskFactorialMoments)
    return 0x0;

  TString gen = "Hijing";

  taskFactorialMoments->SetIsMC(fIsMC, gen);
  taskFactorialMoments->SetYear(year);
  taskFactorialMoments->SetPSbinning(Mbins, Nbins, Mmax);
  taskFactorialMoments->SetRejectPileup(fPileup);
  taskFactorialMoments->SetEta(fEtaMin, fEtaMax);
  taskFactorialMoments->SetCentLim(fCentMin, fCentMax);
  taskFactorialMoments->SetPtArray(ptarrD, nPt);
  taskFactorialMoments->Setfbit(fBit);
  taskFactorialMoments->SetBfield(bField);
  taskFactorialMoments->SetVtxCut(fVxMax, fVyMax, fVzMax);
  taskFactorialMoments->SetDCAXYRangeMax(fDCAXYRangeMax);   // 0.1
  taskFactorialMoments->SetDCAZRangeMax(fDCAZRangeMax);     // 1
  taskFactorialMoments->SetITSClusterCut(fITSClusterCut);   // chi2 per ITS 36
  taskFactorialMoments->SetTPCClusterCut(fTPCClusterCut);   // chi2 per TPC 4
  taskFactorialMoments->SetnTPCrossedrows(fnTPCrossedrows); // 70
  taskFactorialMoments->SetSelfAffine(fSelfAffAnalysis);
  taskFactorialMoments->SetSharedCls(fSharedCls, fSharedRows, fFindableCls);

  // add your taskFactorialMoments to the manager
  mgr->AddTask(taskFactorialMoments);

  // input for taskFactorialMoments and connecting containers to it
  Int_t nout = (fNumberofPtBins) + 1 + fIsMC * (fNumberofPtBins);
  AliAnalysisDataContainer* coutputlist[nout];

  TString str;
  if (fIsMC)
    str = Form("list_MC_%dpt_%d-%dcent_FB%d_%s_%s", fNumberofPtBins, fCentMin, fCentMax, fBit, year.Data(), suffix);
  else
    str = Form("list_DATA_%dpt_%d-%dcent_FB%d_%s_%s", fNumberofPtBins, fCentMin, fCentMax, fBit, year.Data(), suffix);
  coutputlist[0] = mgr->CreateContainer(
    str, TList::Class(), AliAnalysisManager::kOutputContainer,
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

  mgr->ConnectInput(taskFactorialMoments, 0, mgr->GetCommonInputContainer());
  for (Int_t i = 0; i < nout; i++) {
    mgr->ConnectOutput(taskFactorialMoments, i + 1, coutputlist[i]);
  }
  return taskFactorialMoments;
}
