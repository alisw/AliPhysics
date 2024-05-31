/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Jonghan Park (jonghan@cern.ch)                                 *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "TChain.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliHFEextraCuts.h"
// #include "AliHFEmcQATest.h"
#include "AliHFEtools.h"
#include "AliAnalysisUtils.h"
#include "TClonesArray.h"
#include "AliAODMCParticle.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "AliHFEV0taginfo.h"
#include "AliAnalysisTaskBEpp13TeV.h"

class AliAnalysisTaskBEpp13TeV; // your analysis class

using namespace std; // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskBEpp13TeV) // classimp: necessary for root

AliAnalysisTaskBEpp13TeV::AliAnalysisTaskBEpp13TeV():
  AliAnalysisTaskSE(),
  fAOD(0),
  fOutputList(0),
  fTrackQA(0),
  fPIDResponse(0),
  fExtraCuts(0),
  fAODArrayMCInfo(0),
  fAODMCParticle(0),
  fV0Tagger(0),
  fIsMC(false),
  fMinTPCnCrossedRow(70),
  fMinTPCNclsPID(80),
  fMaxTPCchi2(4),
  fMinTPCclsRatio(0.6),
  fMinITSNcls(3),
  fITSlayer(2),
  fTPCnsigmaLow(-1),
  fTPCnsigmaHigh(3),
  fTOFnsigma(3),
  fMultRef(11.7),
  hVtxZbeforCut(0),
  hVtxZafterCut(0),
  hVtxZ(0),
  hNrEvents(0),
  hNrEventsMult(0),
  hNrEventsMult2(0),
  hSPDtracklet(0),
  hNtrklet_vtxZ(0),
  hMultEstimatorAvg(0),
  hSPDtracklet_Corr(0),
  hNtrklet_vtxZ_Corr(0),
  hMultEstimatorAvg_Corr(0),
  fMultEstimatorAvg(0),
  histNchCorr(0),
  funcNchCorr(0),
  hNtrkletCorr(0),
  hNtrkletCorr2(0),
  hSPDtrklet_Nch(0),
  hSPDtrklet_Nch_Corr(0),
  hSPDtrklet_Nch_Corr2(0),
  hNch_vtxZ(0),
  hFilterMask(0),
  hTPCnCrossedRow(0),
  hTPCclsPID(0),
  hTPCchi2(0),
  hTPCclsRatio(0),
  hITSNcls(0),
  hITSlayer(0),
  hDCAxy(0),
  hDCAz(0),
  hPt(0),
  hEta(0),
  hPhi(0),
  hBhadronPt(0),
  hBhadronPtCorr(0),
  hBhadronPtMult(0),
  hBhadronPtMultCorr(0),
  hD0Pt(0),
  hD0PtCorr(0),
  hD0PtMult(0),
  hD0PtMultCorr(0),
  hLcPt(0),
  hD0PtMultBin(0),
  hDsPtMultBin(0),
  hLcPtMultBin(0),
  hGenBePt(0),
  hGenBePtMult(0),
  hGenBePtMult2(0),
  hRecBePt_track(0),
  hRecBePt_tof(0),
  hRecBePt_tofMult(0),
  hRecBePt_tofMult2(0),
  hRecBePt_tpc(0),
  hITSnsigma(0),
  hITSnsigmaTOFcut(0),
  hITSnsigmaTOFTPCcut(0),
  hTPCnsigma(0),
  hTPCnsigmaTOFcut(0),
  hTPCnsigmaTOFcutPt(0),
  hTPCnsigmaQA(0),
  hTPCnsigmaPiQA(0),
  hTOFnsigma(0),
  hTOFnsigmaQA(0),
  hV0ElecTPCnsigma(0),
  hV0ElecTPCnsigmaTOFcut(0),
  hV0ElecTOFnsigmaDeno(0),
  hV0ElecTOFnsigmaNume(0),
  dcaTrack(0),
  dcaTrackMult(0),
  dcaTrackMult2(0),
  dcaPion(0),
  dcaBeauty(0),
  dcaBeautyCorr(0),
  dcaBeautyCorrVar1(0),
  dcaBeautyCorrVar2(0),
  dcaBeautyMult(0),
  dcaBeautyMultCorr(0),
  dcaBeautyMultCorr2(0),
  dcaBeautyMultCorrVar1(0),
  dcaBeautyMultCorrVar2(0),
  DelecVsDmother(0),
  dcaCharm(0),
  dcaDmeson(0),
  dcaDmesonCorr(0),
  dcaDmesonCorrVar1(0),
  dcaDmesonCorrVar2(0),
  dcaDmesonMult(0),
  dcaDmesonMultCorr(0),
  dcaDmesonMultCorrVar1(0),
  dcaDmesonMultCorrVar2(0),
  dcaDzero(0),
  dcaDzeroMult(0),
  dcaDzeroMult2(0),
  dcaDplus(0),
  dcaDplusMult(0),
  dcaDplusMult2(0),
  dcaDsplus(0),
  dcaDsplusMult(0),
  dcaDsplusMult2(0),
  dcaLc(0),
  dcaLcMult(0),
  dcaLcMult2(0),
  dcaDalitz(0),
  dcaDalitzMult(0),
  dcaDalitzMult2(0),
  dcaConv(0),
  dcaConvMult(0),
  dcaConvMult2(0),
  dcaPionSelected(0),
  dcaDalitzSelected(0),
  fBmesonCorrCentLow(0),
  fBmesonCorrCentHigh(0),
  fBmesonCorrMinLow(0),
  fBmesonCorrMinHigh(0),
  fBmesonCorrMaxLow(0),
  fBmesonCorrMaxHigh(0),
  fDmesonCorr(0),
  fDmesonCorrVar1(0),
  fDmesonCorrVar2(0),
  fDmesonCorrMultBin1(0),
  fDmesonCorrMultBin2(0),
  fDmesonCorrMultBin3(0),
  fDmesonCorrMultBin4(0),
  fDmesonCorrMultBin5(0),
  fDmesonCorrMultBin6(0),
  fDmesonCorrMultBin1_var1(0),
  fDmesonCorrMultBin2_var1(0),
  fDmesonCorrMultBin3_var1(0),
  fDmesonCorrMultBin4_var1(0),
  fDmesonCorrMultBin5_var1(0),
  fDmesonCorrMultBin6_var1(0),
  fDmesonCorrMultBin1_var2(0),
  fDmesonCorrMultBin2_var2(0),
  fDmesonCorrMultBin3_var2(0),
  fDmesonCorrMultBin4_var2(0),
  fDmesonCorrMultBin5_var2(0),
  fDmesonCorrMultBin6_var2(0),
  fLcCorr(0),
  fRnd(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
  fParentSelect[0][0] = 411;  // D+
  fParentSelect[0][1] = 421;  // D0
  fParentSelect[0][2] = 431;  // Ds+
  fParentSelect[0][3] = 4122; // Lambdac+
  fParentSelect[0][4] = 4132; // Ksic0
  fParentSelect[0][5] = 4232; // Ksic+
  fParentSelect[0][6] = 4332; // OmegaC0

  fParentSelect[1][0] = 511;  // B0
  fParentSelect[1][1] = 521;  // B+
  fParentSelect[1][2] = 531;  // Bs0
  fParentSelect[1][3] = 5122; // Lambdab0
  fParentSelect[1][4] = 5132; // Ksib-
  fParentSelect[1][5] = 5232; // Ksib0
  fParentSelect[1][6] = 5332; // Omegab-
}
//_____________________________________________________________________________
AliAnalysisTaskBEpp13TeV::AliAnalysisTaskBEpp13TeV(const char *name):
  AliAnalysisTaskSE(name),
  fAOD(0),
  fOutputList(0),
  fTrackQA(0),
  fPIDResponse(0),
  fExtraCuts(0),
  fAODArrayMCInfo(0),
  fAODMCParticle(0),
  fV0Tagger(0),
  fIsMC(false),
  fMinTPCnCrossedRow(70),
  fMinTPCNclsPID(80),
  fMaxTPCchi2(4),
  fMinTPCclsRatio(0.6),
  fMinITSNcls(3),
  fITSlayer(2),
  fTPCnsigmaLow(-1),
  fTPCnsigmaHigh(3),
  fTOFnsigma(3),
  fMultRef(11.7),
  hVtxZbeforCut(0),
  hVtxZafterCut(0),
  hVtxZ(0),
  hNrEvents(0),
  hNrEventsMult(0),
  hNrEventsMult2(0),
  hSPDtracklet(0),
  hNtrklet_vtxZ(0),
  hMultEstimatorAvg(0),
  hSPDtracklet_Corr(0),
  hNtrklet_vtxZ_Corr(0),
  hMultEstimatorAvg_Corr(0),
  fMultEstimatorAvg(0),
  histNchCorr(0),
  funcNchCorr(0),
  hNtrkletCorr(0),
  hNtrkletCorr2(0),
  hSPDtrklet_Nch(0),
  hSPDtrklet_Nch_Corr(0),
  hSPDtrklet_Nch_Corr2(0),
  hNch_vtxZ(0),
  hFilterMask(0),
  hTPCnCrossedRow(0),
  hTPCclsPID(0),
  hTPCchi2(0),
  hTPCclsRatio(0),
  hITSNcls(0),
  hITSlayer(0),
  hDCAxy(0),
  hDCAz(0),
  hPt(0),
  hEta(0),
  hPhi(0),
  hBhadronPt(0),
  hBhadronPtCorr(0),
  hBhadronPtMult(0),
  hBhadronPtMultCorr(0),
  hD0Pt(0),
  hD0PtCorr(0),
  hD0PtMult(0),
  hD0PtMultCorr(0),
  hLcPt(0),
  hD0PtMultBin(0),
  hDsPtMultBin(0),
  hLcPtMultBin(0),
  hGenBePt(0),
  hGenBePtMult(0),
  hGenBePtMult2(0),
  hRecBePt_track(0),
  hRecBePt_tof(0),
  hRecBePt_tofMult(0),
  hRecBePt_tofMult2(0),
  hRecBePt_tpc(0),
  hITSnsigma(0),
  hITSnsigmaTOFcut(0),
  hITSnsigmaTOFTPCcut(0),
  hTPCnsigma(0),
  hTPCnsigmaTOFcut(0),
  hTPCnsigmaTOFcutPt(0),
  hTPCnsigmaQA(0),
  hTPCnsigmaPiQA(0),
  hTOFnsigma(0),
  hTOFnsigmaQA(0),
  hV0ElecTPCnsigma(0),
  hV0ElecTPCnsigmaTOFcut(0),
  hV0ElecTOFnsigmaDeno(0),
  hV0ElecTOFnsigmaNume(0),
  dcaTrack(0),
  dcaTrackMult(0),
  dcaTrackMult2(0),
  dcaPion(0),
  dcaBeauty(0),
  dcaBeautyCorr(0),
  dcaBeautyCorrVar1(0),
  dcaBeautyCorrVar2(0),
  dcaBeautyMult(0),
  dcaBeautyMultCorr(0),
  dcaBeautyMultCorr2(0),
  dcaBeautyMultCorrVar1(0),
  dcaBeautyMultCorrVar2(0),
  DelecVsDmother(0),
  dcaCharm(0),
  dcaDmeson(0),
  dcaDmesonCorr(0),
  dcaDmesonCorrVar1(0),
  dcaDmesonCorrVar2(0),
  dcaDmesonMult(0),
  dcaDmesonMultCorr(0),
  dcaDmesonMultCorrVar1(0),
  dcaDmesonMultCorrVar2(0),
  dcaDzero(0),
  dcaDzeroMult(0),
  dcaDzeroMult2(0),
  dcaDplus(0),
  dcaDplusMult(0),
  dcaDplusMult2(0),
  dcaDsplus(0),
  dcaDsplusMult(0),
  dcaDsplusMult2(0),
  dcaLc(0),
  dcaLcMult(0),
  dcaLcMult2(0),
  dcaDalitz(0),
  dcaDalitzMult(0),
  dcaDalitzMult2(0),
  dcaConv(0),
  dcaConvMult(0),
  dcaConvMult2(0),
  dcaPionSelected(0),
  dcaDalitzSelected(0),
  fBmesonCorrCentLow(0),
  fBmesonCorrCentHigh(0),
  fBmesonCorrMinLow(0),
  fBmesonCorrMinHigh(0),
  fBmesonCorrMaxLow(0),
  fBmesonCorrMaxHigh(0),
  fDmesonCorr(0),
  fDmesonCorrVar1(0),
  fDmesonCorrVar2(0),
  fDmesonCorrMultBin1(0),
  fDmesonCorrMultBin2(0),
  fDmesonCorrMultBin3(0),
  fDmesonCorrMultBin4(0),
  fDmesonCorrMultBin5(0),
  fDmesonCorrMultBin6(0),
  fDmesonCorrMultBin1_var1(0),
  fDmesonCorrMultBin2_var1(0),
  fDmesonCorrMultBin3_var1(0),
  fDmesonCorrMultBin4_var1(0),
  fDmesonCorrMultBin5_var1(0),
  fDmesonCorrMultBin6_var1(0),
  fDmesonCorrMultBin1_var2(0),
  fDmesonCorrMultBin2_var2(0),
  fDmesonCorrMultBin3_var2(0),
  fDmesonCorrMultBin4_var2(0),
  fDmesonCorrMultBin5_var2(0),
  fDmesonCorrMultBin6_var2(0),
  fLcCorr(0),
  fRnd(0)
{
  fParentSelect[0][0] = 411;  // D+
  fParentSelect[0][1] = 421;  // D0
  fParentSelect[0][2] = 431;  // Ds+
  fParentSelect[0][3] = 4122; // Lambdac+
  fParentSelect[0][4] = 4132; // Ksic0
  fParentSelect[0][5] = 4232; // Ksic+
  fParentSelect[0][6] = 4332; // OmegaC0

  fParentSelect[1][0] = 511;  // B0
  fParentSelect[1][1] = 521;  // B+
  fParentSelect[1][2] = 531;  // Bs0
  fParentSelect[1][3] = 5122; // Lambdab0
  fParentSelect[1][4] = 5132; // Ksib-
  fParentSelect[1][5] = 5232; // Ksib0
  fParentSelect[1][6] = 5332; // Omegab-
  // constructor
  DefineInput(0, TChain::Class()); // define the input of the analysis: in this case we take a 'chain' of events
                                   // this chain is created by the analysis manager, so no need to worry about it,
                                   // it does its work automatically
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this case it's a list of histograms
                                   // you can add more output objects by calling DefineOutput(2, classname::Class())
                                   // if you add more output objects, make sure to call PostData for all of them, and to
                                   // make changes to your AddTask macro!
  DefineOutput(2, TList::Class());
  fV0Tagger = new AliHFEV0taginfo("Tagger");
}
//_____________________________________________________________________________
AliAnalysisTaskBEpp13TeV::~AliAnalysisTaskBEpp13TeV(){
  // destructor
  if (fOutputList) delete fOutputList;
  if (fTrackQA) delete fTrackQA;
  if (fExtraCuts) delete fExtraCuts;
  if (fV0Tagger) delete fV0Tagger;
}
//_____________________________________________________________________________
void AliAnalysisTaskBEpp13TeV::UserCreateOutputObjects(){
  fParentSelect[0][0] = 411;  // D+
  fParentSelect[0][1] = 421;  // D0
  fParentSelect[0][2] = 431;  // Ds+
  fParentSelect[0][3] = 4122; // Lambdac+
  fParentSelect[0][4] = 4132; // Ksic0
  fParentSelect[0][5] = 4232; // Ksic+
  fParentSelect[0][6] = 4332; // OmegaC0

  fParentSelect[1][0] = 511;  // B0
  fParentSelect[1][1] = 521;  // B+
  fParentSelect[1][2] = 531;  // Bs0
  fParentSelect[1][3] = 5122; // Lambdab0
  fParentSelect[1][4] = 5132; // Ksib-
  fParentSelect[1][5] = 5232; // Ksib0
  fParentSelect[1][6] = 5332; // Omegab-

  fExtraCuts = new AliHFEextraCuts("hfeExtraCuts", "HFE Extra Cuts");
  // fMCQA = new AliHFEmcQATest;
  // if(fMCQA) fMCQA->Init();

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fTrackQA = new TList();
  fTrackQA->SetOwner(kTRUE);

  int nPtBins = 11;
  double ptbinningX[12] = {1., 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 10.};
  double ptbinningX2[5] = {1.3, 2., 4., 6., 8.};
  // double ptbinningX[22] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 10., 12. };
  // double ptbinningD0[13] = { 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36. };
  double ptbinningD0[11] = {1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 24.};
  double ptbinningD0mult[7] = {1., 2., 4., 6., 8., 12., 24.};
  double ptbinningLc[7] = {1., 2., 4., 6., 8., 12., 24.};
  double ptbinningH[18] = {0.3, 0.5, 0.75, 1., 1.25, 1.5, 2., 2.5, 3., 4., 5., 6., 7., 8., 10., 12., 16., 20.};
  double multbinning[7] = {1., 9., 14., 20., 31., 60., 100.};
  double multbinning2[8] = {1., 9., 14., 20., 31., 50., 70., 100.};
  double evtbinning[2] = {0., 1.};

  int nBinsB = 100;
  double minB = 0.;
  double maxB = 100.;
  double binLimB[nBinsB + 1];
  for (int i = 0; i <= nBinsB; i++) binLimB[i] = minB + (maxB - minB) / nBinsB * (double)i;

  int nBinsPID = 400;
  double minPID = -10;
  double maxPID = 10;
  double binLimPID[nBinsPID + 1];
  for (int i = 0; i <= nBinsPID; i++) binLimPID[i] = minPID + (maxPID - minPID) / nBinsPID * (double)i;

  int nBinsIP = 4000;
  double minIP = -0.2;
  double maxIP = 0.2;
  double binLimIP[nBinsIP + 1];
  for (int i = 0; i <= nBinsIP; i++) binLimIP[i] = minIP + (maxIP - minIP) / nBinsIP * (double)i;

  // example of a histogram
  hVtxZbeforCut = new TH1F("hVtxZbeforCut", "vertex z", 600, -30, 30);
  fTrackQA->Add(hVtxZbeforCut);

  hVtxZafterCut = new TH1F("hVtxZafterCut", "vertex z", 600, -30, 30);
  fTrackQA->Add(hVtxZafterCut);

  hVtxZ = new TH2F("hVtxZ", "vertex z", 600, -30, 30, 2, 0., 2.);
  fTrackQA->Add(hVtxZ);

  hNrEvents = new TH1F("hNrEvents", "number of events", 1, 0., 1.);
  fOutputList->Add(hNrEvents);

  hNrEventsMult = new TH2F("hNrEventsMult", "number of events", 1, evtbinning, 6, multbinning);
  fOutputList->Add(hNrEventsMult);

  hNrEventsMult2 = new TH2F("hNrEventsMult2", "number of events", 1, evtbinning, 7, multbinning2);
  fOutputList->Add(hNrEventsMult2);

  hSPDtracklet = new TH1F("hSPDtracklet", "SPD tracklet distribution", 180, 0, 180);
  fOutputList->Add(hSPDtracklet);

  hNtrklet_vtxZ = new TH2F("hNtrklet_vtxZ", "", 300, -15., 15., 300, 0, 300);
  fOutputList->Add(hNtrklet_vtxZ);

  hMultEstimatorAvg = new TProfile("hMultEstimatorAvg", "", 300, -15., 15.);
  fOutputList->Add(hMultEstimatorAvg);

  hSPDtracklet_Corr = new TH1F("hSPDtracklet_Corr", "SPD tracklet distribution", 180, 0, 180);
  fOutputList->Add(hSPDtracklet_Corr);

  hNtrklet_vtxZ_Corr = new TH2F("hNtrklet_vtxZ_Corr", "", 300, -15., 15., 300, 0, 300);
  fOutputList->Add(hNtrklet_vtxZ_Corr);

  hMultEstimatorAvg_Corr = new TProfile("hMultEstimatorAvg_Corr", "", 300, -15., 15.);
  fOutputList->Add(hMultEstimatorAvg_Corr);

  hNtrkletCorr = new TH1F("hNtrkletCorr", "", 180, 0, 180);
  fOutputList->Add(hNtrkletCorr);

  hNtrkletCorr2 = new TH1F("hNtrkletCorr2", "", 180, 0, 180);
  fOutputList->Add(hNtrkletCorr2);

  hSPDtrklet_Nch = new TH2F("hSPDtrklet_Nch", "", 180, 0, 180, 180, 0, 180);
  fOutputList->Add(hSPDtrklet_Nch);

  hSPDtrklet_Nch_Corr = new TH2F("hSPDtrklet_Nch_Corr", "", 180, 0, 180, 180, 0, 180);
  fOutputList->Add(hSPDtrklet_Nch_Corr);

  hSPDtrklet_Nch_Corr2 = new TH2F("hSPDtrklet_Nch_Corr2", "", 180, 0, 180, 180, 0, 180);
  fOutputList->Add(hSPDtrklet_Nch_Corr2);

  hNch_vtxZ = new TH2F("hNch_vtxZ", "", 300, -15, 15, 300, 0, 300);
  fOutputList->Add(hNch_vtxZ);

  hFilterMask = new TH1F("hFilterMask", "", 2, 0., 2.);
  fTrackQA->Add(hFilterMask);

  hTPCnCrossedRow = new TH1F("hTPCnCrossedRow", "", 200, 0., 200.);
  fTrackQA->Add(hTPCnCrossedRow);

  hTPCclsPID = new TH1F("hTPCclsPID", "", 200, 0., 200.);
  fTrackQA->Add(hTPCclsPID);

  hTPCchi2 = new TH1F("hTPCchi2", "", 100, 0., 10.);
  fTrackQA->Add(hTPCchi2);

  hTPCclsRatio = new TH1F("hTPCclsRatio", "", 15, 0., 1.5);
  fTrackQA->Add(hTPCclsRatio);

  hITSNcls = new TH1F("hITSNcls", "", 10, 0., 10.);
  fTrackQA->Add(hITSNcls);

  hITSlayer = new TH1F("hITSlayer", "", 3, 0.5, 3.5);
  fTrackQA->Add(hITSlayer);

  hDCAxy = new TH1F("hDCAxy", "", 600, -3., 3.);
  fTrackQA->Add(hDCAxy);

  hDCAz = new TH1F("hDCAz", "", 600, -3., 3.);
  fTrackQA->Add(hDCAz);

  hPt = new TH1F("hPt", "pt; (GeV/c)", 300, 0., 30.);
  fTrackQA->Add(hPt);

  hEta = new TH1F("hEta", "", 200, -1., 1.);
  fTrackQA->Add(hEta);

  hPhi = new TH1F("hPhi", "", 700, -0.5, 6.5);
  fTrackQA->Add(hPhi);

  hBhadronPt = new TH1F("hBhadronPt", "", nBinsB, binLimB);
  fOutputList->Add(hBhadronPt);

  hBhadronPtCorr = new TH1F("hBhadronPtCorr", "", nBinsB, binLimB);
  fOutputList->Add(hBhadronPtCorr);

  hBhadronPtMult = new TH2F("hBhadronPtMult", "", nBinsB, binLimB, 6, multbinning);
  fOutputList->Add(hBhadronPtMult);

  hBhadronPtMultCorr = new TH2F("hBhadronPtMultCorr", "", nBinsB, binLimB, 6, multbinning);
  fOutputList->Add(hBhadronPtMultCorr);

  hD0Pt = new TH1F("hD0Pt", "", 10, ptbinningD0);
  fOutputList->Add(hD0Pt);

  hD0PtCorr = new TH1F("hD0PtCorr", "", 10, ptbinningD0);
  fOutputList->Add(hD0PtCorr);

  hD0PtMult = new TH2F("hD0PtMult", "", 6, ptbinningD0mult, 6, multbinning);
  fOutputList->Add(hD0PtMult);

  hD0PtMultCorr = new TH2F("hD0PtMultCorr", "", 6, ptbinningD0mult, 6, multbinning);
  fOutputList->Add(hD0PtMultCorr);

  hLcPt = new TH1F("hLcPt", "", 6, ptbinningLc);
  fOutputList->Add(hLcPt);

  hD0PtMultBin = new TH1D("hD0PtMultBin", "", 6, ptbinningLc);
  fOutputList->Add(hD0PtMultBin);

  hDsPtMultBin = new TH1D("hDsPtMultBin", "", 6, ptbinningLc);
  fOutputList->Add(hDsPtMultBin);

  hLcPtMultBin = new TH1D("hLcPtMultBin", "", 6, ptbinningLc);
  fOutputList->Add(hLcPtMultBin);

  hGenBePt = new TH1F("hGenBePt", "", nPtBins, ptbinningX);
  fOutputList->Add(hGenBePt);

  hGenBePtMult = new TH2F("hGenBePtMult", "", nPtBins, ptbinningX, 6, multbinning);
  fOutputList->Add(hGenBePtMult);

  hGenBePtMult2 = new TH2F("hGenBePtMult2", "", nPtBins, ptbinningX, 7, multbinning2);
  fOutputList->Add(hGenBePtMult2);

  hRecBePt_track = new TH1F("hRecBePt_track", "", nPtBins, ptbinningX);
  fOutputList->Add(hRecBePt_track);

  hRecBePt_tof = new TH1F("hRecBePt_tof", "", nPtBins, ptbinningX);
  fOutputList->Add(hRecBePt_tof);

  hRecBePt_tofMult = new TH2F("hRecBePt_tofMult", "", nPtBins, ptbinningX, 6, multbinning);
  fOutputList->Add(hRecBePt_tofMult);

  hRecBePt_tofMult2 = new TH2F("hRecBePt_tofMult2", "", nPtBins, ptbinningX, 7, multbinning2);
  fOutputList->Add(hRecBePt_tofMult2);

  hRecBePt_tpc = new TH1F("hRecBePt_tpc", "", nPtBins, ptbinningX);
  fOutputList->Add(hRecBePt_tpc);

  hITSnsigma = new TH2F("hITSnsigma", "n#sigma_{ITS} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{ITS}", 700, 0., 14., 800, -20., 20.);
  fOutputList->Add(hITSnsigma);

  hITSnsigmaTOFcut = new TH2F("hITSnsigmaTOFcut", "n#sigma_{ITS} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{ITS}", 700, 0., 14., 800, -20., 20.);
  fOutputList->Add(hITSnsigmaTOFcut);

  hITSnsigmaTOFTPCcut = new TH2F("hITSnsigmaTOFTPCcut", "n#sigma_{ITS} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{ITS}", 700, 0., 14., 800, -20., 20.);
  fOutputList->Add(hITSnsigmaTOFTPCcut);

  hTPCnsigma = new TH2F("hTPCnsigma", "n#sigma_{TPC} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TPC}", 700, 0., 14., 800, -20., 20.);
  fOutputList->Add(hTPCnsigma);

  hTPCnsigmaTOFcut = new TH2F("hTPCnsigmaTOFcut", "n#sigma_{TPC |TOFcut} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TPC}", 700, 0., 14., 800, -20., 20.);
  fOutputList->Add(hTPCnsigmaTOFcut);

  hTPCnsigmaTOFcutPt = new TH2F("hTPCnsigmaTOFcutPt", "n#sigma_{TPC |TOFcut} vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); n#sigma_{TPC}", nPtBins, ptbinningX, nBinsPID, binLimPID);
  fOutputList->Add(hTPCnsigmaTOFcutPt);

  hTPCnsigmaQA = new TH2F("hTPCnsigmaQA", "n#sigma_{TPC} QA vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TPC}", 700, 0., 14., 800, -20., 20.);
  fOutputList->Add(hTPCnsigmaQA);

  hTPCnsigmaPiQA = new TH2F("hTPCnsigmaPiQA", "n#sigma_{TPC |TOFcut} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TPC}", 700, 0., 14., 800, -20., 20.);
  fOutputList->Add(hTPCnsigmaPiQA);

  hTOFnsigma = new TH2F("hTOFnsigma", "n#sigma_{TOF} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TOF}", 700, 0., 14., 800, -20., 20.);
  fOutputList->Add(hTOFnsigma);

  hTOFnsigmaQA = new TH2F("hTOFnsigmaQA", "n#sigma_{TOF} vs #it{p}; #it{p} (GeV/#it{c}); n#sigma_{TOF}", 700, 0., 14., 800, -20., 20.);
  fOutputList->Add(hTOFnsigmaQA);

  hV0ElecTPCnsigma = new TH2F("hV0ElecTPCnsigma", "V0 elec TPC nsigma", nPtBins, ptbinningX, nBinsPID, binLimPID);
  fOutputList->Add(hV0ElecTPCnsigma);

  hV0ElecTPCnsigmaTOFcut = new TH2F("hV0ElecTPCnsigmaTOFcut", "V0 elec TPC nsigma with TOF cut", nPtBins, ptbinningX, nBinsPID, binLimPID);
  fOutputList->Add(hV0ElecTPCnsigmaTOFcut);

  hV0ElecTOFnsigmaDeno = new TH2F("hV0ElecTOFnsigmaDeno", "", nPtBins, ptbinningX, nBinsPID, binLimPID);
  fOutputList->Add(hV0ElecTOFnsigmaDeno);

  hV0ElecTOFnsigmaNume = new TH2F("hV0ElecTOFnsigmaNume", "", nPtBins, ptbinningX, nBinsPID, binLimPID);
  fOutputList->Add(hV0ElecTOFnsigmaNume);

  dcaTrack = new TH2F("dcaTrack", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaTrack);

  dcaTrackMult = new TH3F("dcaTrackMult", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaTrackMult);

  dcaTrackMult2 = new TH3F("dcaTrackMult2", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 7, multbinning2);
  fOutputList->Add(dcaTrackMult2);

  dcaPion = new TH2F("dcaPion", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaPion);

  dcaBeauty = new TH2F("dcaBeauty", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaBeauty);

  dcaBeautyCorr = new TH2F("dcaBeautyCorr", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaBeautyCorr);

  dcaBeautyCorrVar1 = new TH2F("dcaBeautyCorrVar1", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaBeautyCorrVar1);

  dcaBeautyCorrVar2 = new TH2F("dcaBeautyCorrVar2", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaBeautyCorrVar2);

  dcaBeautyMult = new TH3F("dcaBeautyMult", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaBeautyMult);

  dcaBeautyMultCorr = new TH3F("dcaBeautyMultCorr", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaBeautyMultCorr);

  dcaBeautyMultCorr2 = new TH3F("dcaBeautyMultCorr2", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 7, multbinning2);
  fOutputList->Add(dcaBeautyMultCorr2);

  dcaBeautyMultCorrVar1 = new TH3F("dcaBeautyMultCorrVar1", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaBeautyMultCorrVar1);

  dcaBeautyMultCorrVar2 = new TH3F("dcaBeautyMultCorrVar2", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaBeautyMultCorrVar2);

  DelecVsDmother = new TH2F("DelecVsDmother", "", 10, ptbinningD0, nPtBins, ptbinningX);
  fOutputList->Add(DelecVsDmother);

  dcaCharm = new TH2F("dcaCharm", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaCharm);

  dcaDmeson = new TH2F("dcaDmeson", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDmeson);

  dcaDmesonCorr = new TH2F("dcaDmesonCorr", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDmesonCorr);

  dcaDmesonCorrVar1 = new TH2F("dcaDmesonCorrVar1", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDmesonCorrVar1);

  dcaDmesonCorrVar2 = new TH2F("dcaDmesonCorrVar2", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDmesonCorrVar2);

  dcaDmesonMult = new TH3F("dcaDmesonMult", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaDmesonMult);

  dcaDmesonMultCorr = new TH3F("dcaDmesonMultCorr", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaDmesonMultCorr);

  dcaDmesonMultCorrVar1 = new TH3F("dcaDmesonMultCorrVar1", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaDmesonMultCorrVar1);

  dcaDmesonMultCorrVar2 = new TH3F("dcaDmesonMultCorrVar2", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaDmesonMultCorrVar2);

  dcaDzero = new TH2F("dcaDzero", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDzero);

  dcaDzeroMult = new TH3F("dcaDzeroMult", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaDzeroMult);

  dcaDzeroMult2 = new TH3F("dcaDzeroMult2", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 7, multbinning2);
  fOutputList->Add(dcaDzeroMult2);

  dcaDplus = new TH2F("dcaDplus", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDplus);

  dcaDplusMult = new TH3F("dcaDplusMult", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaDplusMult);

  dcaDplusMult2 = new TH3F("dcaDplusMult2", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 7, multbinning2);
  fOutputList->Add(dcaDplusMult2);

  dcaDsplus = new TH2F("dcaDsplus", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDsplus);

  dcaDsplusMult = new TH3F("dcaDsplusMult", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaDsplusMult);

  dcaDsplusMult2 = new TH3F("dcaDsplusMult2", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 7, multbinning2);
  fOutputList->Add(dcaDsplusMult2);

  dcaLc = new TH2F("dcaLc", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaLc);

  dcaLcMult = new TH3F("dcaLcMult", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaLcMult);

  dcaLcMult2 = new TH3F("dcaLcMult2", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 7, multbinning2);
  fOutputList->Add(dcaLcMult2);

  dcaDalitz = new TH2F("dcaDalitz", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDalitz);

  dcaDalitzMult = new TH3F("dcaDalitzMult", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaDalitzMult);

  dcaDalitzMult2 = new TH3F("dcaDalitzMult2", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 7, multbinning2);
  fOutputList->Add(dcaDalitzMult2);

  dcaConv = new TH2F("dcaConv", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaConv);

  dcaConvMult = new TH3F("dcaConvMult", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 6, multbinning);
  fOutputList->Add(dcaConvMult);

  dcaConvMult2 = new TH3F("dcaConvMult2", "", nPtBins, ptbinningX, nBinsIP, binLimIP, 7, multbinning2);
  fOutputList->Add(dcaConvMult2);

  dcaPionSelected = new TH2D("dcaPionSelected", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaPionSelected);

  dcaDalitzSelected = new TH2D("dcaDalitzSelected", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDalitzSelected);

  fRnd = new TRandom3(0);

  
  fDmesonCorrMultBin1 = new TF1("fDmesonCorrMultBin1", "(1/1.038302)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin1->FixParameter(0, 6.80225e+00);
  fDmesonCorrMultBin1->FixParameter(1, 1.44558e+01);
  fDmesonCorrMultBin1->FixParameter(2, 1.20305e+00);
  fDmesonCorrMultBin1->FixParameter(3, 2.23080e+00);
  fDmesonCorrMultBin1->FixParameter(4, 1.63620e+00);
  fDmesonCorrMultBin1->FixParameter(5, 5.37149e+00);
  fDmesonCorrMultBin1->FixParameter(6, 3.89111e-01);
  fDmesonCorrMultBin1->FixParameter(7, 1.56282e+01);


  fDmesonCorrMultBin2 = new TF1("fDmesonCorrMultBin2", "(1/1.146876)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin2->FixParameter(0, 6.71029e+00);
  fDmesonCorrMultBin2->FixParameter(1, 2.37157e+00);
  fDmesonCorrMultBin2->FixParameter(2, 2.64314e-02);
  fDmesonCorrMultBin2->FixParameter(3, 4.57065e-01);
  fDmesonCorrMultBin2->FixParameter(4, 2.26636e+00);
  fDmesonCorrMultBin2->FixParameter(5, 2.33583e+00);
  fDmesonCorrMultBin2->FixParameter(6, 2.61640e-01);
  fDmesonCorrMultBin2->FixParameter(7, 6.22554e-01);

  fDmesonCorrMultBin3 = new TF1("fDmesonCorrMultBin3", "(1/1.146876)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin3->FixParameter(0, 5.47737e+00);
  fDmesonCorrMultBin3->FixParameter(1, 2.40699e+00);
  fDmesonCorrMultBin3->FixParameter(2, 4.90741e-02);
  fDmesonCorrMultBin3->FixParameter(3, 8.71897e-01);
  fDmesonCorrMultBin3->FixParameter(4, 5.67444e+00);
  fDmesonCorrMultBin3->FixParameter(5, 2.44888e+00);
  fDmesonCorrMultBin3->FixParameter(6, 1.29810e-01);
  fDmesonCorrMultBin3->FixParameter(7, 3.10050e-01);

  fDmesonCorrMultBin4 = new TF1("fDmesonCorrMultBin4", "(1/1.109676)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin4->FixParameter(0, 5.29029e+00);
  fDmesonCorrMultBin4->FixParameter(1, 2.37934e+00);
  fDmesonCorrMultBin4->FixParameter(2, 4.68272e-02);
  fDmesonCorrMultBin4->FixParameter(3, 1.09655e+00);
  fDmesonCorrMultBin4->FixParameter(4, 6.34543e+00);
  fDmesonCorrMultBin4->FixParameter(5, 2.38059e+00);
  fDmesonCorrMultBin4->FixParameter(6, 1.01869e-01);
  fDmesonCorrMultBin4->FixParameter(7, 4.71452e-01);

  fDmesonCorrMultBin5 = new TF1("fDmesonCorrMultBin5", "(1/1.054812)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin5->FixParameter(0, 4.79937e+00);
  fDmesonCorrMultBin5->FixParameter(1, 2.99546e+00);
  fDmesonCorrMultBin5->FixParameter(2, 2.03118e-01);
  fDmesonCorrMultBin5->FixParameter(3, 4.60345e+00);
  fDmesonCorrMultBin5->FixParameter(4, 6.29361e+00);
  fDmesonCorrMultBin5->FixParameter(5, 2.86359e+00);
  fDmesonCorrMultBin5->FixParameter(6, 3.11677e-01);
  fDmesonCorrMultBin5->FixParameter(7, 4.18799e+00);

  fDmesonCorrMultBin6 = new TF1("fDmesonCorrMultBin6", "(1/1.268888)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin6->FixParameter(0, 6.05818e+00);
  fDmesonCorrMultBin6->FixParameter(1, 3.07025e+00);
  fDmesonCorrMultBin6->FixParameter(2, 2.42087e-02);
  fDmesonCorrMultBin6->FixParameter(3, 4.78095e-01);
  fDmesonCorrMultBin6->FixParameter(4, 2.54327e+00);
  fDmesonCorrMultBin6->FixParameter(5, 2.03081e+00);
  fDmesonCorrMultBin6->FixParameter(6, 2.16818e-01);
  fDmesonCorrMultBin6->FixParameter(7, 1.47484e+00);

  fDmesonCorrMultBin1_var1 = new TF1("fDmesonCorrMultBin1_var1", "(1/1.096594)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin1_var1->FixParameter(0, 6.80225e+00);
  fDmesonCorrMultBin1_var1->FixParameter(1, 1.44558e+01);
  fDmesonCorrMultBin1_var1->FixParameter(2, 1.20305e+00);
  fDmesonCorrMultBin1_var1->FixParameter(3, 2.23080e+00);
  fDmesonCorrMultBin1_var1->FixParameter(4, 1.63620e+00);
  fDmesonCorrMultBin1_var1->FixParameter(5, 5.37149e+00);
  fDmesonCorrMultBin1_var1->FixParameter(6, 4.09111e-01);
  fDmesonCorrMultBin1_var1->FixParameter(7, 1.56282e+01);

  fDmesonCorrMultBin2_var1 = new TF1("fDmesonCorrMultBin2_var1", "(1/1.306064)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin2_var1->FixParameter(0, 6.71029e+00);
  fDmesonCorrMultBin2_var1->FixParameter(1, 2.37157e+00);
  fDmesonCorrMultBin2_var1->FixParameter(2, 2.64314e-02);
  fDmesonCorrMultBin2_var1->FixParameter(3, 4.57065e-01);
  fDmesonCorrMultBin2_var1->FixParameter(4, 2.26636e+00);
  fDmesonCorrMultBin2_var1->FixParameter(5, 2.33583e+00);
  fDmesonCorrMultBin2_var1->FixParameter(6, 3.31640e-01);
  fDmesonCorrMultBin2_var1->FixParameter(7, 6.22554e-01);

  fDmesonCorrMultBin3_var1 = new TF1("fDmesonCorrMultBin3_var1", "(1/1.100170)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin3_var1->FixParameter(0, 5.47737e+00);
  fDmesonCorrMultBin3_var1->FixParameter(1, 2.40699e+00);
  fDmesonCorrMultBin3_var1->FixParameter(2, 4.90741e-02);
  fDmesonCorrMultBin3_var1->FixParameter(3, 8.71897e-01);
  fDmesonCorrMultBin3_var1->FixParameter(4, 5.67444e+00);
  fDmesonCorrMultBin3_var1->FixParameter(5, 2.44888e+00);
  fDmesonCorrMultBin3_var1->FixParameter(6, 1.39810e-01);
  fDmesonCorrMultBin3_var1->FixParameter(7, 3.10050e-01);

  fDmesonCorrMultBin4_var1 = new TF1("fDmesonCorrMultBin4_var1", "(1/1.060610)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin4_var1->FixParameter(0, 5.29029e+00);
  fDmesonCorrMultBin4_var1->FixParameter(1, 2.37934e+00);
  fDmesonCorrMultBin4_var1->FixParameter(2, 4.68272e-02);
  fDmesonCorrMultBin4_var1->FixParameter(3, 1.09655e+00);
  fDmesonCorrMultBin4_var1->FixParameter(4, 6.34543e+00);
  fDmesonCorrMultBin4_var1->FixParameter(5, 2.38059e+00);
  fDmesonCorrMultBin4_var1->FixParameter(6, 1.11869e-01);
  fDmesonCorrMultBin4_var1->FixParameter(7, 4.71452e-01);

  fDmesonCorrMultBin5_var1 = new TF1("fDmesonCorrMultBin5_var1", "(1/1.100586)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin5_var1->FixParameter(0, 4.79937e+00);
  fDmesonCorrMultBin5_var1->FixParameter(1, 2.99546e+00);
  fDmesonCorrMultBin5_var1->FixParameter(2, 2.03118e-01);
  fDmesonCorrMultBin5_var1->FixParameter(3, 4.60345e+00);
  fDmesonCorrMultBin5_var1->FixParameter(4, 6.29361e+00);
  fDmesonCorrMultBin5_var1->FixParameter(5, 2.86359e+00);
  fDmesonCorrMultBin5_var1->FixParameter(6, 3.31677e-01);
  fDmesonCorrMultBin5_var1->FixParameter(7, 4.18799e+00);

  fDmesonCorrMultBin6_var1 = new TF1("fDmesonCorrMultBin6_var1", "(1/1.852835)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin6_var1->FixParameter(0, 6.05818e+00);
  fDmesonCorrMultBin6_var1->FixParameter(1, 3.07025e+00);
  fDmesonCorrMultBin6_var1->FixParameter(2, 2.42087e-02);
  fDmesonCorrMultBin6_var1->FixParameter(3, 4.78095e-01);
  fDmesonCorrMultBin6_var1->FixParameter(4, 2.54327e+00);
  fDmesonCorrMultBin6_var1->FixParameter(5, 2.03081e+00);
  fDmesonCorrMultBin6_var1->FixParameter(6, 3.16818e-01);
  fDmesonCorrMultBin6_var1->FixParameter(7, 1.47484e+00);

  fDmesonCorrMultBin1_var2 = new TF1("fDmesonCorrMultBin1_var2", "(1/0.981082)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin1_var2->FixParameter(0, 6.80225e+00);
  fDmesonCorrMultBin1_var2->FixParameter(1, 1.44558e+01);
  fDmesonCorrMultBin1_var2->FixParameter(2, 1.20305e+00);
  fDmesonCorrMultBin1_var2->FixParameter(3, 2.23080e+00);
  fDmesonCorrMultBin1_var2->FixParameter(4, 1.63620e+00);
  fDmesonCorrMultBin1_var2->FixParameter(5, 5.37149e+00);
  fDmesonCorrMultBin1_var2->FixParameter(6, 3.69111e-01);
  fDmesonCorrMultBin1_var2->FixParameter(7, 1.56282e+01);

  fDmesonCorrMultBin2_var2 = new TF1("fDmesonCorrMultBin2_var2", "(1/1.057155)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin2_var2->FixParameter(0, 6.71029e+00);
  fDmesonCorrMultBin2_var2->FixParameter(1, 2.37157e+00);
  fDmesonCorrMultBin2_var2->FixParameter(2, 2.64314e-02);
  fDmesonCorrMultBin2_var2->FixParameter(3, 4.57065e-01);
  fDmesonCorrMultBin2_var2->FixParameter(4, 2.26636e+00);
  fDmesonCorrMultBin2_var2->FixParameter(5, 2.33583e+00);
  fDmesonCorrMultBin2_var2->FixParameter(6, 2.11640e-01);
  fDmesonCorrMultBin2_var2->FixParameter(7, 6.22554e-01);

  fDmesonCorrMultBin3_var2 = new TF1("fDmesonCorrMultBin3_var2", "(1/1.174481)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin3_var2->FixParameter(0, 5.47737e+00);
  fDmesonCorrMultBin3_var2->FixParameter(1, 2.40699e+00);
  fDmesonCorrMultBin3_var2->FixParameter(2, 4.90741e-02);
  fDmesonCorrMultBin3_var2->FixParameter(3, 8.71897e-01);
  fDmesonCorrMultBin3_var2->FixParameter(4, 5.67444e+00);
  fDmesonCorrMultBin3_var2->FixParameter(5, 2.44888e+00);
  fDmesonCorrMultBin3_var2->FixParameter(6, 1.19810e-01);
  fDmesonCorrMultBin3_var2->FixParameter(7, 3.10050e-01);

  fDmesonCorrMultBin4_var2 = new TF1("fDmesonCorrMultBin4_var2", "(1/1.173915)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin4_var2->FixParameter(0, 5.29029e+00);
  fDmesonCorrMultBin4_var2->FixParameter(1, 2.37934e+00);
  fDmesonCorrMultBin4_var2->FixParameter(2, 4.68272e-02);
  fDmesonCorrMultBin4_var2->FixParameter(3, 1.09655e+00);
  fDmesonCorrMultBin4_var2->FixParameter(4, 6.34543e+00);
  fDmesonCorrMultBin4_var2->FixParameter(5, 2.38059e+00);
  fDmesonCorrMultBin4_var2->FixParameter(6, 0.91869e-01);
  fDmesonCorrMultBin4_var2->FixParameter(7, 4.71452e-01);

  fDmesonCorrMultBin5_var2 = new TF1("fDmesonCorrMultBin5_var2", "(1/1.010567)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin5_var2->FixParameter(0, 4.79937e+00);
  fDmesonCorrMultBin5_var2->FixParameter(1, 2.99546e+00);
  fDmesonCorrMultBin5_var2->FixParameter(2, 2.03118e-01);
  fDmesonCorrMultBin5_var2->FixParameter(3, 4.60345e+00);
  fDmesonCorrMultBin5_var2->FixParameter(4, 6.29361e+00);
  fDmesonCorrMultBin5_var2->FixParameter(5, 2.86359e+00);
  fDmesonCorrMultBin5_var2->FixParameter(6, 2.91677e-01);
  fDmesonCorrMultBin5_var2->FixParameter(7, 4.18799e+00);

  fDmesonCorrMultBin6_var2 = new TF1("fDmesonCorrMultBin6_var2", "(1/0.806620)*(([0]*(([1]-1.)*([1]-2.))/([1]*[2]*([1]*[2]+[3]*([1]-2.)))*pow(1.+(sqrt([3]*[3] + x*x)-[3])/([1]*[2]),-[1])) / ([4]*(([5]-1.)*([5]-2.))/([5]*[6]*([5]*[6]+[3]*([5]-2.)))*pow(1.+(sqrt([7]*[7] + x*x)-[7])/([5]*[6]),-[5])))", 1., 100.);
  fDmesonCorrMultBin6_var2->FixParameter(0, 6.05818e+00);
  fDmesonCorrMultBin6_var2->FixParameter(1, 3.07025e+00);
  fDmesonCorrMultBin6_var2->FixParameter(2, 2.42087e-02);
  fDmesonCorrMultBin6_var2->FixParameter(3, 4.78095e-01);
  fDmesonCorrMultBin6_var2->FixParameter(4, 2.54327e+00);
  fDmesonCorrMultBin6_var2->FixParameter(5, 2.03081e+00);
  fDmesonCorrMultBin6_var2->FixParameter(6, 1.16818e-01);
  fDmesonCorrMultBin6_var2->FixParameter(7, 1.47484e+00);


  PostData(1, fOutputList); // postdata will notify the analysis manager of changes / updates to the
                            // fOutputList object. the manager will in the end take care of writing your output to file
                            // so it needs to know what's in the output
  PostData(2, fTrackQA);
}
//_____________________________________________________________________________
void AliAnalysisTaskBEpp13TeV::UserExec(Option_t *){

  fAOD = dynamic_cast<AliAODEvent *>(InputEvent());
  if(!fAOD) return;

  if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts", "HFE Extra Cuts");
  }
  fExtraCuts->SetRecEventInfo(fAOD);

  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse){
    AliDebug(1, "Using default PID Response");
    fPIDResponse = AliHFEtools::GetDefaultPID(false, fInputEvent->IsA() == AliAODEvent::Class());
  }

  // Initialize V0 electron tagger
  if(fV0Tagger){
    fV0Tagger->Reset();
    fV0Tagger->TagV0Tracks(fAOD);
  }

  bool isSelected = (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  if(!isSelected) return;

  // Event selection
  if(!PassEventCuts(fAOD)) return;

  // Pile-up removal
  if(PassPileUpEvent(fAOD)) return;

  // Look for kink mother
  double *fListOfMotherKink = 0;
  int fNumberOfVertices = 0;
  int fNumberOfMotherKink = 0;

  fNumberOfVertices = fAOD->GetNumberOfVertices();
  fListOfMotherKink = new double[fNumberOfVertices];

  for(int iVertex = 0; iVertex < fNumberOfVertices; iVertex++){
    AliAODVertex *aodvtx = fAOD->GetVertex(iVertex);
    if(!aodvtx) continue;
    if(aodvtx->GetType() == AliAODVertex::kKink){
      AliAODTrack *mother = (AliAODTrack *)aodvtx->GetParent();
      if(!mother) continue;
      int idmother = mother->GetID();
      fListOfMotherKink[fNumberOfMotherKink] = idmother;
      fNumberOfMotherKink++;
    }
  }

  // Multiplicity
  const AliAODVertex *vtx = fAOD->GetPrimaryVertex();
  double vtxZ = vtx->GetZ();

  int nAcceta = 0;
  AliAODTracklets *tracklets = ((AliAODEvent *)fAOD)->GetTracklets();
  // int nTracklets = tracklets->GetNumberOfTracklets();
  for(int iTrklet = 0; iTrklet < tracklets->GetNumberOfTracklets(); iTrklet++){
    double fTrkletTheta = tracklets->GetTheta(iTrklet);
    double fTrkletEta = -TMath::Log(TMath::Tan(fTrkletTheta / 2.));
    if(TMath::Abs(fTrkletEta) < 1.)
      nAcceta++;
  }
  hSPDtracklet->Fill(nAcceta);
  hNtrklet_vtxZ->Fill(vtxZ, nAcceta);
  hMultEstimatorAvg->Fill(vtxZ, nAcceta);

  int Corrected_Ntr = nAcceta;

  TProfile *estimatorAvg = GetEstimatorHistogram();
  if(estimatorAvg) Corrected_Ntr = static_cast<int>(GetCorrectedNtracklets(estimatorAvg, nAcceta, vtxZ, fMultRef));

  hSPDtracklet_Corr->Fill(Corrected_Ntr);
  hNtrklet_vtxZ_Corr->Fill(vtxZ, Corrected_Ntr);
  hMultEstimatorAvg_Corr->Fill(vtxZ, Corrected_Ntr);

  // generated MC loop
  if(fIsMC){
    fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fAODArrayMCInfo){
      AliError("No AOD MC particles");
      return;
    }
    // if(fMCQA) fMCQA->SetMCArray(fAODArrayMCInfo);

    for(int iMC = 0; iMC < fAODArrayMCInfo->GetEntries(); iMC++){
      fAODMCParticle = (AliAODMCParticle *)fAODArrayMCInfo->At(iMC);

      int hf = -999;
      double hfpt = -999., hfeta = -999., wghtD = -999., wghtmultD = -999., wghtB = -999.;
      hf = GetHeavyFlavours(fAODMCParticle, hfpt, hfeta);

      if(TMath::Abs(hfeta) < 0.5){
        if(hf==kPromptD0){
          hD0Pt->Fill(hfpt);
          hD0PtMult->Fill(hfpt, Corrected_Ntr);
          wghtD = fDmesonCorr->Eval(hfpt);
          hD0PtCorr->Fill(hfpt, wghtD);
		  if(Corrected_Ntr>=1 && Corrected_Ntr<9) wghtmultD = fDmesonCorrMultBin1->Eval(hfpt);
		  else if(Corrected_Ntr>=9 && Corrected_Ntr<14) wghtmultD = fDmesonCorrMultBin2->Eval(hfpt);
		  else if(Corrected_Ntr>=14 && Corrected_Ntr<20) wghtmultD = fDmesonCorrMultBin3->Eval(hfpt);
		  else if(Corrected_Ntr>=20 && Corrected_Ntr<31) wghtmultD = fDmesonCorrMultBin4->Eval(hfpt);
		  else if(Corrected_Ntr>=31 && Corrected_Ntr<60) wghtmultD = fDmesonCorrMultBin5->Eval(hfpt);
		  else if(Corrected_Ntr>=60 && Corrected_Ntr<100) wghtmultD = fDmesonCorrMultBin6->Eval(hfpt);
		  else wghtmultD = 1.;
          hD0PtMultCorr->Fill(hfpt, Corrected_Ntr, wghtmultD);
        }
        if(hf==kPromptLc) hLcPt->Fill(hfpt);
      }
      if(TMath::Abs(hfeta<0.8)){
        if(hf==kPromptB || hf==kNonPromptD){
          hBhadronPtMult->Fill(hfpt, Corrected_Ntr);
          hBhadronPt->Fill(hfpt);
          if(hfpt<=3.5) wghtB = fBmesonCorrCentLow->Eval(hfpt);
          if(hfpt>3.5) wghtB = fBmesonCorrCentHigh->Eval(hfpt);
          hBhadronPtCorr->Fill(hfpt, wghtB);
          hBhadronPtMultCorr->Fill(hfpt, Corrected_Ntr, wghtB);
        }
      }

      int src = -999, srcPdg = -999;
      double srcPt = -999.;
      src = GetElecSource(fAODMCParticle, true, srcPt, srcPdg);

      if(TMath::Abs(fAODMCParticle->Eta())<0.8){
        if(src==kDirectBeauty || src==kBeautyCharm){
          hGenBePt->Fill(fAODMCParticle->Pt());
          hGenBePtMult->Fill(fAODMCParticle->Pt(), Corrected_Ntr);
          hGenBePtMult2->Fill(fAODMCParticle->Pt(), Corrected_Ntr);
        }
      }
    }
  }

  if(fIsMC){
    double nchCorr = 1.;
    double nchCorr2 = 1.;
    if(Corrected_Ntr>50 && Corrected_Ntr<=100) nchCorr = funcNchCorr->Eval(Corrected_Ntr);
    else if(Corrected_Ntr<=50) nchCorr = histNchCorr->GetBinContent(histNchCorr->FindBin(Corrected_Ntr));
    nchCorr2 = histNchCorr->GetBinContent(histNchCorr->FindBin(Corrected_Ntr));
    hNtrkletCorr->Fill(Corrected_Ntr, nchCorr);
    hNtrkletCorr2->Fill(Corrected_Ntr, nchCorr2);
    int nch = GetNcharged();

    hSPDtrklet_Nch->Fill(Corrected_Ntr, nch);
    hSPDtrklet_Nch_Corr->Fill(Corrected_Ntr, nch, nchCorr);
    hSPDtrklet_Nch_Corr2->Fill(Corrected_Ntr, nch, nchCorr2);
    hNch_vtxZ->Fill(vtxZ, nch);
	
	//printf("%d---%d\n", Corrected_Ntr, nch);
  }

  double fBz = -999.;
  if(fAOD->GetMagneticField()<0) fBz = -1.;
  else if(fAOD->GetMagneticField()>0) fBz = 1.;
  else return;
  
  hNrEvents->Fill(0.5);
  hNrEventsMult->Fill(0.5, Corrected_Ntr);
  hNrEventsMult2->Fill(0.5, Corrected_Ntr);
  for(int iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++){
    AliAODTrack *aodTrack = static_cast<AliAODTrack *>(fAOD->GetTrack(iTracks));
    if (!aodTrack) continue;

    // Reject kink
    bool kinkmotherpass = true;
    for(int kinkmother=0; kinkmother<fNumberOfMotherKink; kinkmother++){
      if(aodTrack->GetID()==fListOfMotherKink[kinkmother]){
        kinkmotherpass = false;
        continue;
      }
    }
    if(!kinkmotherpass) continue;

    // Track selection
    if(!PassTrackCuts(aodTrack)) continue;

    double pt = aodTrack->Pt();
    double hfeImpactParam = -999., hfeImpactParamResol = -999.;
    fExtraCuts->GetHFEImpactParameters((AliVTrack *)aodTrack, hfeImpactParam, hfeImpactParamResol);
    double IP = hfeImpactParam * fBz * aodTrack->Charge();

    int mcelectronSource = -999, mcelectronSourcePDG = -999;
    double mcelectronSourcePt = -999.;
    if(fIsMC){
      fAODMCParticle = NULL;

      int label = TMath::Abs(aodTrack->GetLabel());
      if(label<fAODArrayMCInfo->GetEntriesFast()) fAODMCParticle = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(label));
      if(fAODMCParticle){
        AliDebug(2, "Associated MC particle found");
        mcelectronSource = GetElecSource(fAODMCParticle, true, mcelectronSourcePt, mcelectronSourcePDG);
      }
      // Fill beauty dca information
      if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm){
        // if(mcelectronSourcePt>70.) continue;
        hRecBePt_track->Fill(pt);
        dcaBeauty->Fill(pt, IP);
        dcaBeautyMult->Fill(pt, IP, Corrected_Ntr);

        double wghtB = -99., wghtBvar1 = -99., wghtBvar2 = -99.;
        if(pt>mcelectronSourcePt){
          wghtB = 1.;
          wghtBvar1 = 1.;
          wghtBvar2 = 1.;
        }
        else{
          if(mcelectronSourcePt<=3.5){
            wghtB = fBmesonCorrCentLow->Eval(mcelectronSourcePt);
            wghtBvar1 = fBmesonCorrMinLow->Eval(mcelectronSourcePt);
            wghtBvar2 = fBmesonCorrMaxLow->Eval(mcelectronSourcePt);
          }
          else if(mcelectronSourcePt>3.5){
            wghtB = fBmesonCorrCentHigh->Eval(mcelectronSourcePt);
            wghtBvar1 = fBmesonCorrMinHigh->Eval(mcelectronSourcePt);
            wghtBvar2 = fBmesonCorrMaxHigh->Eval(mcelectronSourcePt);
          }
        }

        // B hadron dca correction
        double rndmB = fRnd->Rndm();
        if(rndmB<wghtB){
		  dcaBeautyCorr->Fill(pt, IP);
		  dcaBeautyMultCorr->Fill(pt, IP, Corrected_Ntr);
		  dcaBeautyMultCorr2->Fill(pt, IP, Corrected_Ntr);
		}
        if(rndmB<wghtBvar1){
		  dcaBeautyCorrVar1->Fill(pt, IP);
		  dcaBeautyMultCorrVar1->Fill(pt, IP, Corrected_Ntr);
		}
        if(rndmB<wghtBvar2){
		  dcaBeautyCorrVar2->Fill(pt, IP);
		  dcaBeautyMultCorrVar2->Fill(pt, IP, Corrected_Ntr);
		}
      }

      // Fill charm dca information
      if(mcelectronSource==kDirectCharm){
        dcaCharm->Fill(pt, IP);
        if(TMath::Abs(mcelectronSourcePDG)==421 || TMath::Abs(mcelectronSourcePDG)==411 || TMath::Abs(mcelectronSourcePDG)==431){
          DelecVsDmother->Fill(mcelectronSourcePt, pt);
          dcaDmeson->Fill(pt, IP);
          dcaDmesonMult->Fill(pt, IP, Corrected_Ntr);
          double wghtD = -99., wghtDvar1 = -99., wghtDvar2 = -99.;
          double wghtmultD = -99., wghtmultDvar1 = -99., wghtmultDvar2 = -99.;
          if(pt>mcelectronSourcePt){
            wghtD = 1.;
            wghtDvar1 = 1.;
            wghtDvar2 = 1.;
            wghtmultD = 1.;
            wghtmultDvar1 = 1.;
            wghtmultDvar2 = 1.;
          }else{
            wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            wghtDvar1 = fDmesonCorrVar1->Eval(mcelectronSourcePt);
            wghtDvar2 = fDmesonCorrVar2->Eval(mcelectronSourcePt);
            // if(pt>=1. && pt<1.1)  wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=1.1 && pt<1.3) wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=1.3 && pt<1.5) wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=1.5 && pt<2.)  wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=2. && pt<2.5)  wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=2.5 && pt<3.)  wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=3. && pt<4.)   wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=4 && pt<5.)    wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=5. && pt<6.)   wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=6. && pt<8.)   wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=8. && pt<10.)  wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
            // if(pt>=10.)           wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			if(Corrected_Ntr>=1 && Corrected_Ntr<9){
			  wghtmultD = fDmesonCorrMultBin1->Eval(mcelectronSourcePt);
			  wghtmultDvar1 = fDmesonCorrMultBin1_var1->Eval(mcelectronSourcePt);
			  wghtmultDvar2 = fDmesonCorrMultBin1_var2->Eval(mcelectronSourcePt);
			}else if(Corrected_Ntr>=9 && Corrected_Ntr<14){
			  wghtmultD = fDmesonCorrMultBin2->Eval(mcelectronSourcePt);
			  wghtmultDvar1 = fDmesonCorrMultBin2_var1->Eval(mcelectronSourcePt);
			  wghtmultDvar2 = fDmesonCorrMultBin2_var2->Eval(mcelectronSourcePt);
			}else if(Corrected_Ntr>=14 && Corrected_Ntr<20){
			  wghtmultD = fDmesonCorrMultBin3->Eval(mcelectronSourcePt);
			  wghtmultDvar1 = fDmesonCorrMultBin3_var1->Eval(mcelectronSourcePt);
			  wghtmultDvar2 = fDmesonCorrMultBin3_var2->Eval(mcelectronSourcePt);
			}else if(Corrected_Ntr>=20 && Corrected_Ntr<31){
			  wghtmultD = fDmesonCorrMultBin4->Eval(mcelectronSourcePt);
			  wghtmultDvar1 = fDmesonCorrMultBin4_var1->Eval(mcelectronSourcePt);
			  wghtmultDvar2 = fDmesonCorrMultBin4_var2->Eval(mcelectronSourcePt);
			}else if(Corrected_Ntr>=31 && Corrected_Ntr<60){
			  wghtmultD = fDmesonCorrMultBin5->Eval(mcelectronSourcePt);
			  wghtmultDvar1 = fDmesonCorrMultBin5_var1->Eval(mcelectronSourcePt);
			  wghtmultDvar2 = fDmesonCorrMultBin5_var2->Eval(mcelectronSourcePt);
			}else if(Corrected_Ntr>=60 && Corrected_Ntr<100){
			  wghtmultD = fDmesonCorrMultBin6->Eval(mcelectronSourcePt);
			  wghtmultDvar1 = fDmesonCorrMultBin6_var1->Eval(mcelectronSourcePt);
			  wghtmultDvar2 = fDmesonCorrMultBin6_var2->Eval(mcelectronSourcePt);
			}else{
			  wghtmultD = 1.;
			  wghtmultDvar1 = 1.;
			  wghtmultDvar2 = 1.;
			}
		  }
          double rndmD = fRnd->Rndm();
          if(rndmD<wghtD){
            dcaDmesonCorr->Fill(pt, IP);
            if(TMath::Abs(mcelectronSourcePDG)==421) dcaDzero->Fill(pt, IP);
            if(TMath::Abs(mcelectronSourcePDG)==411) dcaDplus->Fill(pt, IP);
            if(TMath::Abs(mcelectronSourcePDG)==431) dcaDsplus->Fill(pt, IP);
          }
		  if(rndmD<wghtmultD){
            dcaDmesonMultCorr->Fill(pt, IP, Corrected_Ntr);
            if(TMath::Abs(mcelectronSourcePDG)==421){
			  dcaDzeroMult->Fill(pt, IP, Corrected_Ntr);
			  dcaDzeroMult2->Fill(pt, IP, Corrected_Ntr);
			  hD0PtMultBin->Fill(mcelectronSourcePt);
			}
            if(TMath::Abs(mcelectronSourcePDG)==411){
			  dcaDplusMult->Fill(pt, IP, Corrected_Ntr);
			  dcaDplusMult2->Fill(pt, IP, Corrected_Ntr);
			}
            if(TMath::Abs(mcelectronSourcePDG)==431){
			  dcaDsplusMult->Fill(pt, IP, Corrected_Ntr);
			  dcaDsplusMult2->Fill(pt, IP, Corrected_Ntr);
			  hDsPtMultBin->Fill(mcelectronSourcePt);
			}
		  }
          if(rndmD<wghtDvar1) dcaDmesonCorrVar1->Fill(pt, IP);
          if(rndmD<wghtDvar2) dcaDmesonCorrVar2->Fill(pt, IP);
		  if(rndmD<wghtmultDvar1) dcaDmesonMultCorrVar1->Fill(pt, IP, Corrected_Ntr);
		  if(rndmD<wghtmultDvar2) dcaDmesonMultCorrVar2->Fill(pt, IP, Corrected_Ntr);
        }
        if(TMath::Abs(mcelectronSourcePDG)==4122){
		  dcaLc->Fill(pt, IP);
		  dcaLcMult->Fill(pt, IP, Corrected_Ntr);
		  dcaLcMult2->Fill(pt, IP, Corrected_Ntr);
		  hLcPtMultBin->Fill(mcelectronSourcePt);
		}
      }
      // Fill Dalitz dca information
      if(mcelectronSource>=kPi0 && mcelectronSource<=kK2P){
		dcaDalitz->Fill(pt, IP);
		dcaDalitzMult->Fill(pt, IP, Corrected_Ntr);
		dcaDalitzMult2->Fill(pt, IP, Corrected_Ntr);
	  }
      // Fill conversion dca information
      if(mcelectronSource>=kGammaPi0 && mcelectronSource<=kGammaK2P){
		dcaConv->Fill(pt, IP);
		dcaConvMult->Fill(pt, IP, Corrected_Ntr);
		dcaConvMult2->Fill(pt, IP, Corrected_Ntr);
	  }
    }

    // electron identification
    double fITSnSigma = fPIDResponse->NumberOfSigmasITS(aodTrack, AliPID::kElectron);
    double fTPCnSigma = fPIDResponse->NumberOfSigmasTPC(aodTrack, AliPID::kElectron);
    double fTOFnSigma = fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kElectron);

    hITSnsigma->Fill(aodTrack->P(), fITSnSigma);
    hTPCnsigma->Fill(aodTrack->P(), fTPCnSigma);
    hTOFnsigma->Fill(aodTrack->P(), fTOFnSigma);

    // V0 electrons from systematic studies of TOF eID
    AliPID::EParticleType myv0pid = fV0Tagger->GetV0Info(aodTrack->GetID()); /// enum EParticleType: kElectron = 0, kMuon = 1, kPion = 2, etc
    if(myv0pid == AliPID::kElectron){
      hV0ElecTPCnsigma->Fill(aodTrack->Pt(), fTPCnSigma);
      if(TMath::Abs(fTOFnSigma)<=fTOFnsigma)
        hV0ElecTPCnsigmaTOFcut->Fill(aodTrack->Pt(), fTPCnSigma);
      // TOF eID systematics
      if(fTPCnSigma>=fTPCnsigmaLow && fTPCnSigma<=fTPCnsigmaHigh){
        hV0ElecTOFnsigmaDeno->Fill(aodTrack->Pt(), fTOFnSigma);
        if(TMath::Abs(fTOFnSigma)<=fTOFnsigma)
          hV0ElecTOFnsigmaNume->Fill(aodTrack->Pt(), fTOFnSigma);
      }
    }

    if(TMath::Abs(fTOFnSigma)>fTOFnsigma) continue;
    if(fIsMC){
      if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm){
        hRecBePt_tof->Fill(pt);
        hRecBePt_tofMult->Fill(pt, Corrected_Ntr);
        hRecBePt_tofMult2->Fill(pt, Corrected_Ntr);
	  }
    }

    hITSnsigmaTOFcut->Fill(aodTrack->P(), fITSnSigma);
    hTPCnsigmaTOFcut->Fill(aodTrack->P(), fTPCnSigma);
    hTPCnsigmaTOFcutPt->Fill(pt, fTPCnSigma);

    if(fTPCnSigma>-5 && fTPCnSigma<-3){
      hTPCnsigmaPiQA->Fill(aodTrack->P(), fTPCnSigma);
      dcaPion->Fill(pt, hfeImpactParam * fBz * aodTrack->Charge());
    }

    if(fTPCnSigma<fTPCnsigmaLow || fTPCnSigma>fTPCnsigmaHigh) continue;
    if(fIsMC){
      if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm){
        hRecBePt_tpc->Fill(pt);
	  }
      if(mcelectronSource>=kPi0 && mcelectronSource<=kK2P){
		dcaDalitzSelected->Fill(pt, IP);
	  }
	  if(TMath::Abs(fAODMCParticle->GetPdgCode()==211)){
		dcaPionSelected->Fill(pt, IP);
	  }
    }

    hITSnsigmaTOFTPCcut->Fill(aodTrack->P(), fITSnSigma);
    hTPCnsigmaQA->Fill(aodTrack->P(), fTPCnSigma);
    hTOFnsigmaQA->Fill(aodTrack->P(), fTOFnSigma);

    dcaTrack->Fill(pt, IP);
    dcaTrackMult->Fill(pt, IP, Corrected_Ntr);
    dcaTrackMult2->Fill(pt, IP, Corrected_Ntr);
  }

  PostData(1, fOutputList); // stream the results the analysis of this event to
                            // the output manager which will take care of writing
                            // it to a file
  PostData(2, fTrackQA);
}
//_____________________________________________________________________________
void AliAnalysisTaskBEpp13TeV::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}
//__________________________________________________________________
bool AliAnalysisTaskBEpp13TeV::PassEventCuts(AliAODEvent *event){

  // event selection cuts
  const AliAODVertex *vtx = event->GetPrimaryVertex();
  if(!vtx)
    return false;
  if(vtx->GetNContributors()<2)
    return false;

  // cut on the primary vertex on z position
  hVtxZbeforCut->Fill(vtx->GetZ());
  hVtxZ->Fill(vtx->GetZ(), 0);
  if(TMath::Abs(vtx->GetZ())>10)
    return false;
  hVtxZafterCut->Fill(vtx->GetZ());
  hVtxZ->Fill(vtx->GetZ(), 1);

  return true;
}
//_________________________________________________________________
bool AliAnalysisTaskBEpp13TeV::PassPileUpEvent(AliAODEvent *event){
  // This function checks if there was a pile up reconstructed with SPD
  bool isPileupfromSPDmulbins = event->IsPileupFromSPDInMultBins();
  if(isPileupfromSPDmulbins)
    return true;

  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(5); // Multi Vertex pileup selection
  utils.SetMaxPlpChi2MV(5);    // max value of Chi2perNDF of the pileup and multi-vertex
  utils.SetMinWDistMV(15);     // min of the sqrt of weighted distance between the primary and the pileup vertex, multi-vertex
  utils.SetCheckPlpFromDifferentBCMV(false);
  bool isPileupFromMV = utils.IsPileUpMV(event);
  return isPileupFromMV;
}
//_________________________________________________________________
double AliAnalysisTaskBEpp13TeV::GetCorrectedNtracklets(TProfile *estimatorAvg, double rawNtr, double vtxz, double refmult){
  if(TMath::Abs(vtxz) > 10) return rawNtr;
  if (!estimatorAvg) return rawNtr;
  double Ntr_mean = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxz));
  double deltaN = rawNtr * ((refmult/Ntr_mean) - 1);
  double correctedNtr = rawNtr + (deltaN > 0 ? 1 : -1) * gRandom->Poisson(TMath::Abs(deltaN));
  return correctedNtr;
}
//_________________________________________________________________
bool AliAnalysisTaskBEpp13TeV::PassTrackCuts(AliAODTrack *track){

  if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;
  if(TMath::Abs(track->Eta())>0.8) return false;
  if(track->Pt()<0.5 || track->Pt()>15.) return false;

  // basic tracking
  ULong_t status = track->GetStatus();
  if(!((status & AliVTrack::kITSrefit) && (status & AliVTrack::kTPCrefit))) return false;

  // TPC cut
  unsigned short TPCnCrossedRow = track->GetTPCNCrossedRows();
  unsigned short TPCclsPID = track->GetTPCsignalN();
  double TPCchi2 = track->Chi2perNDF();
  unsigned short findableTPC = track->GetTPCNclsF();
  unsigned short TPCsignalN = track->GetTPCsignalN();

  double FoundOverFindable = (findableTPC ? static_cast<float>(TPCnCrossedRow) / static_cast<float>(findableTPC) : 0);
  if(TPCnCrossedRow<fMinTPCnCrossedRow || TPCclsPID<fMinTPCNclsPID || TPCchi2>fMaxTPCchi2 || FoundOverFindable<fMinTPCclsRatio) return false;

  // ITS cut
  int ITSnCls = track->GetITSNcls();
  if (ITSnCls<fMinITSNcls) return false;
  if(fITSlayer==AliHFEextraCuts::kFirst){
    if(!(track->HasPointOnITSLayer(0))) return false;
    hITSlayer->Fill(1);
  }
  else if(fITSlayer==AliHFEextraCuts::kAny){
    if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return false;
    hITSlayer->Fill(2);
  }
  else if(fITSlayer==AliHFEextraCuts::kBoth){
    if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return false;
    hITSlayer->Fill(3);
  }
  else return false;

  // dca cut
  float dcaxy = -999.;
  float dcaz = -999.;
  fExtraCuts->GetImpactParameters((AliVTrack *)track, dcaxy, dcaz);
  if(TMath::Abs(dcaxy)>1. || TMath::Abs(dcaz)>2.) return false;

  hFilterMask->Fill(track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA));
  hEta->Fill(track->Eta());
  hPhi->Fill(track->Phi());
  hPt->Fill(track->Pt());
  hDCAxy->Fill(dcaxy);
  hDCAz->Fill(dcaz);
  hTPCnCrossedRow->Fill(TPCnCrossedRow);
  hTPCclsPID->Fill(TPCclsPID);
  hTPCchi2->Fill(TPCchi2);
  hTPCclsRatio->Fill(FoundOverFindable);
  hITSNcls->Fill(ITSnCls);
  return true;
}
//_________________________________________________________________________________________________________
int AliAnalysisTaskBEpp13TeV::GetElecSource(const AliAODMCParticle *const mcpart, double &mpt, int &mpdg){
  // printf("Modified version of source selection \n");
  if(!mcpart) return kMisID;
  if(!fAODArrayMCInfo) return -1;

  if(TMath::Abs(mcpart->GetPdgCode())!=11) return kElse;

  int origin = -1;
  bool isFinalOpenCharm = false;

  int iLabel = mcpart->GetMother();
  if((iLabel<0) || (iLabel>=fAODArrayMCInfo->GetEntriesFast())){
    AliDebug(1, "label is out of range, return\n");
    return -1;
  }

  AliAODMCParticle *mctrack = NULL; // will change all the time
  int tmpMomLabel = 0;
  if(!(mctrack=dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(iLabel))))) return -1;
  AliAODMCParticle *partMother = mctrack;     // mtrack
  AliAODMCParticle *partMotherCopy = mctrack; // mtrack
  int maPdgcode = TMath::Abs(partMother->GetPdgCode());   // mpdg
  mpt = partMother->Pt();                     // mpt
  mpdg = partMother->GetPdgCode();            // mpdg
  int gmaPdgcode, ggmaPdgcode;
  double gmpt, ggmpt;
  int gmpdg, ggmpdg;

  // if the mother is charmed hadron
  if((int(maPdgcode/100.)%10)==4 || (int(maPdgcode/1000.)%10)==4){
    if(maPdgcode==411 || maPdgcode==421 || maPdgcode==431 || maPdgcode==4122 || maPdgcode==4132 || maPdgcode==4232 || maPdgcode==4332){
      mpt = partMother->Pt();
      mpdg = partMother->GetPdgCode();
      isFinalOpenCharm = kTRUE;
    }
    if(!isFinalOpenCharm){
      return -1;
    }

    // iterate until find B hadron as a  mother
    for(int i=1; i<100; i++){
      int jLabel = partMother->GetMother();
      if(jLabel==-1){
        return kDirectCharm;
      }
      if(jLabel<0 || jLabel>=fAODArrayMCInfo->GetEntriesFast()){
        AliDebug(1, "Stack label is negative, return\n");
        return -1;
      }

      if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel))))){
        return -1;
      }
      int grandMaPDG = TMath::Abs(mctrack->GetPdgCode());
      if(grandMaPDG==511 || grandMaPDG == 521 || grandMaPDG == 531 || grandMaPDG == 5122 || grandMaPDG == 5132 || grandMaPDG == 5232 || grandMaPDG == 5332){
        mpt = mctrack->Pt();
        mpdg = mctrack->GetPdgCode();
        return kBeautyCharm;
      }
      partMother = mctrack;
    } // end of iteration
  }

  // if the mother is beauty hadron
  else if((int(maPdgcode/100.)%10)==5 || (int(maPdgcode/1000.)%10)==5){
    if(maPdgcode==511 || maPdgcode==521 || maPdgcode==531 || maPdgcode==5122 || maPdgcode==5132 || maPdgcode==5232 || maPdgcode==5332){
      mpt = partMotherCopy->Pt();
      mpdg = partMotherCopy->GetPdgCode();
      return kDirectBeauty;
    }
  }

  // if the mother is gamma
  else if(maPdgcode == 22){
    tmpMomLabel = partMotherCopy->GetMother(); // mother of photon
    mpt = partMotherCopy->Pt();                // pT of photon
    mpdg = partMotherCopy->GetPdgCode();
    if (tmpMomLabel == -1) return kGamma; // no grandmother
    if((tmpMomLabel<0) || (tmpMomLabel>=fAODArrayMCInfo->GetEntriesFast())){
      return -1;
    }
    if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))){
      return -1;
    }
    partMother = mctrack;                 // gmtrack
    partMotherCopy = mctrack;             // gmtrack
    mpt = partMother->Pt();               // grand mother pT
    mpdg = partMother->GetPdgCode();      // grand mother PDG
    maPdgcode = partMother->GetPdgCode(); // grand mother PDG

    // check if the ligth meson is the decay product of heavy mesons
    tmpMomLabel = partMother->GetMother(); // grand grand mother of photon
    if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())){ // grand grand mother
      if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))){
        partMother = mctrack;                  // ggmtrack
        gmaPdgcode = partMother->GetPdgCode(); // grand grand mother PDG
        mpt = partMother->Pt();                // grand grand mother pT
        mpdg = partMother->GetPdgCode();       // grand grand mother pT
        gmpt = partMother->Pt();               // grand grand mother pt
        gmpdg = partMother->GetPdgCode();      // grand grand mother pT

        if(maPdgcode==111){
          mpt = gmpt;
          mpdg = gmpdg;
          if(gmaPdgcode==310) return kGammaK0s2P;
          else if(gmaPdgcode==130) return kGammaK0l2P;
          else if(TMath::Abs(gmaPdgcode)==321) return kGammaK2P;
          else if(TMath::Abs(gmaPdgcode)==3122) return kGammaLamda2P;
          else if(gmaPdgcode==3222) return kGammaSigma2P;
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kGammaPi0;
        }
        else if(maPdgcode==221){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kGammaEta;
        }
        else if(maPdgcode==223){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kGammaOmega;
        }
        else if(maPdgcode==333){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kGammaPhi;
        }
        else if(maPdgcode==331){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kGammaEtaPrime;
        }
        else if(maPdgcode==113){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kGammaRho0;
        }
        else origin = kElse; // grand grand mother but nothing we identify
      }                   // mctrack grandgrandmother
    }
    else{
      // grandmother is primary
      if(maPdgcode==111){
        return kGammaPi0;
      }
      else if(maPdgcode==221){
        return kGammaEta;
      }
      else if(maPdgcode==223){
        return kGammaOmega;
      }
      else if(maPdgcode==333){
        return kGammaPhi;
      }
      else if(maPdgcode==331){
        return kGammaEtaPrime;
      }
      else if(maPdgcode==113){
        return kGammaRho0;
      }
      else origin = kElse; // grandmother is primary but nothing we identify
    }
    return origin;
  }

  // if the mother is light meson
  else{

    tmpMomLabel = partMotherCopy->GetMother(); // grand mother
    mpt = partMotherCopy->Pt();                // mother pT
    mpdg = partMotherCopy->GetPdgCode();       // mother PDG
    if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())){ // grand mother
      if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))){
        partMother = mctrack;                  // grand mother
        gmaPdgcode = partMother->GetPdgCode(); // grand mother PDG
        mpt = partMother->Pt();                // grand mother pT
        mpdg = partMother->GetPdgCode();       // grand mother PDG
        gmpt = partMother->Pt();               // grand mother pT
        gmpdg = partMother->GetPdgCode();      // grand mother PDG

        if(maPdgcode==111){
          mpt = gmpt;
          mpdg = gmpdg;
          if(gmaPdgcode == 310) return kK0s2P;
          else if(gmaPdgcode==130) return kK0l2P;
          else if(TMath::Abs(gmaPdgcode)==321) return kK2P;
          else if(TMath::Abs(gmaPdgcode)==3122) return kLamda2P;
          else if(gmaPdgcode==3222) return kSigma2P;
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kPi0;
        }
        else if(maPdgcode==221){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kEta;
        }
        else if(maPdgcode==223){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kOmega;
        }
        else if(maPdgcode==333){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kPhi;
        }
        else if(maPdgcode==331){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kEtaPrime;
        }
        else if(maPdgcode==113){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kRho0;
        }
        else if(maPdgcode==321){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kKe3;
        }
        else if(maPdgcode==130){
          mpt = partMotherCopy->Pt();
          mpdg = partMotherCopy->GetPdgCode();
          return kK0L;
        }
        else origin = kElse; // grandmother but nothing we identidy
      }                   // mctrack grandmother
    }
    else{
      // no grandmother
      if(maPdgcode==111) return kPi0;
      else if(maPdgcode==221) return kEta;
      else if(maPdgcode==223) return kOmega;
      else if(maPdgcode==333) return kPhi;
      else if(maPdgcode==331) return kEtaPrime;
      else if(maPdgcode==113) return kRho0;
      else if(maPdgcode==321) return kKe3;
      else if(maPdgcode==130) return kK0L;
      else origin = kElse; // mother but nothing we identify
    }
  } // mother is something different from J/psi,charm,beauty or gamma

  return origin;
}
//_______________________________________________________________________________________________________________
int AliAnalysisTaskBEpp13TeV::GetHeavyFlavours(const AliAODMCParticle *const mcpart, double &hfpt, double &hfeta){

  if (!mcpart) return -1;
  if (!fAODArrayMCInfo) return -1;

  int pdgHF = TMath::Abs(mcpart->GetPdgCode());
  hfpt = mcpart->Pt();
  hfeta = mcpart->Eta();
  if (!(pdgHF/100==4 || pdgHF/100==5 || pdgHF/1000==4 || pdgHF/1000==5)) return -1;

  AliAODMCParticle *mctrack = NULL;
  AliAODMCParticle *partMother = NULL;

  if(pdgHF == 411 || pdgHF == 421 || pdgHF == 431 || pdgHF == 4122 || pdgHF == 4132 || pdgHF == 4232 || pdgHF == 4332){
    // iterate until find B hadron as a mother
    int jLabel = -999;
    int maPdgcode = -999;
    for(int i = 1; i < 100; i++){
      if(i==1) jLabel = mcpart->GetMother();
      else jLabel = partMother->GetMother();

      if(jLabel==-1){
        if(pdgHF==421) return kPromptD0;
        if(pdgHF==4122) return kPromptLc;
      }
      if(jLabel<0 || jLabel>=fAODArrayMCInfo->GetEntriesFast()){
        AliDebug(1, "Stack label is negative, return\n");
        return -1;
      }
      if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel))))){
        return -1;
      }
      maPdgcode = TMath::Abs(mctrack->GetPdgCode());
      if(maPdgcode==511 || maPdgcode==521 || maPdgcode==531 || maPdgcode==5122 || maPdgcode==5132 || maPdgcode==5232 || maPdgcode==5332){
        hfpt = mctrack->Pt();
        hfeta = mctrack->Eta();
        return kNonPromptD;
      }
      partMother = mctrack;
    } // end of iteration
  }

  // prompt B mesons
  else if(pdgHF==511 || pdgHF==521 || pdgHF==531 || pdgHF==5122 || pdgHF==5132 || pdgHF==5232 || pdgHF==5332){
    return kPromptB;
  }

  return -1;
}
//_______________________________________________________________________________________________________________
int AliAnalysisTaskBEpp13TeV::GetElecSource(const AliAODMCParticle *const mcpart, bool isElec, double &mpt, int &mpdg) const
{
  // printf("Original version of source selection \n");

  if (!mcpart)
    return -1;
  if (!fAODArrayMCInfo)
    return -1;

  // if(isElec) if(TMath::Abs(mcpart->GetPdgCode()) != AliHFEmcQATest::kElectronPDG) return kMisID;
  if (isElec)
    if (TMath::Abs(mcpart->GetPdgCode()) != 11)
      return kMisID;

  int origin = -1;
  bool isFinalOpenCharm = kFALSE;

  int iLabel = mcpart->GetMother();
  if ((iLabel < 0) || (iLabel >= fAODArrayMCInfo->GetEntriesFast()))
  {
    AliDebug(1, "label is out of range, return\n");
    return -1;
  }

  AliAODMCParticle *mctrack = NULL; // will change all the time
  int tmpMomLabel = 0;
  if (!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(iLabel)))))
    return -1;
  AliAODMCParticle *partMother = mctrack;
  AliAODMCParticle *partMotherCopy = mctrack;
  int maPdgcode = partMother->GetPdgCode();
  mpt = partMother->Pt();
  mpdg = maPdgcode;
  int grmaPdgcode;
  int ggrmaPdgcode;
  double gmpt, ggmpt;

  // if the mother is charmed hadron
  if (maPdgcode == 443)
  {
    // J/spi
    int jLabel = partMother->GetMother();
    if ((jLabel >= 0) && (jLabel < fAODArrayMCInfo->GetEntriesFast()))
    {
      if ((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel)))))
      {
        int grandMaPDG = mctrack->GetPdgCode();
        mpt = mctrack->Pt();
        mpdg = grandMaPDG;
        if ((int(TMath::Abs(grandMaPDG) / 100.) % 10) == kBeauty || (int(TMath::Abs(grandMaPDG) / 1000.) % 10) == kBeauty)
          return kB2Jpsi;
      }
    }
    return kJpsi;
  }
  else if ((int(maPdgcode / 100.) % 10) == kCharm || (int(maPdgcode / 1000.) % 10) == kCharm)
  {
    // charm
    for (int i = 0; i < fNparents; i++)
    {
      if (maPdgcode == fParentSelect[0][i])
      {
        mpt = partMother->Pt();
        mpdg = maPdgcode;
        isFinalOpenCharm = kTRUE;
      }
    }
    if (!isFinalOpenCharm)
      return -1;

    // iterate until you find B hadron as a mother or become top ancester
    for (int i = 1; i < fgkMaxIter; i++)
    {

      int jLabel = partMother->GetMother();
      if (jLabel == -1)
        return kDirectCharm;

      if ((jLabel < 0) || (jLabel >= fAODArrayMCInfo->GetEntriesFast()))
      {
        AliDebug(1, "Stack label is negative, return\n");
        return -1;
      }

      // if there is an ancester
      if (!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel)))))
        return -1;
      int grandMaPDG = mctrack->GetPdgCode();
      for (int j = 0; j < fNparents; j++)
      {
        if (TMath::Abs(grandMaPDG) == fParentSelect[1][j])
        {
          mpt = mctrack->Pt();
          mpdg = grandMaPDG;
          return kBeautyCharm;
        }
      }
      partMother = mctrack;
    } // end of iteration

  } // end of charm
  else if ((int(maPdgcode / 100.) % 10) == kBeauty || (int(maPdgcode / 1000.) % 10) == kBeauty)
  {
    // beauty
    for (int i = 0; i < fNparents; i++)
    {
      if (maPdgcode == fParentSelect[1][i])
      {
        mpt = partMotherCopy->Pt();
        mpdg = maPdgcode;
        return kDirectBeauty;
      }
    }
  } // end of beauty
  else if (maPdgcode == 22)
  {
    // conversion
    tmpMomLabel = partMotherCopy->GetMother();
    mpt = partMotherCopy->Pt(); // pT of photon
    mpdg = maPdgcode;
    if (tmpMomLabel == -1)
      return kGamma;
    if ((tmpMomLabel < 0) || (tmpMomLabel >= fAODArrayMCInfo->GetEntriesFast()))
      return -1;
    if (!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel)))))
      return -1;

    partMother = mctrack;
    partMotherCopy = mctrack;
    mpt = partMother->Pt();
    maPdgcode = partMother->GetPdgCode();
    mpdg = maPdgcode;

    // check if the ligth meson is the decay product of heavy mesons
    tmpMomLabel = partMother->GetMother();
    if ((tmpMomLabel >= 0) && (tmpMomLabel < fAODArrayMCInfo->GetEntriesFast()))
    { // grandgrandmother
      if ((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel)))))
      {
        partMother = mctrack;
        grmaPdgcode = partMother->GetPdgCode();
        mpt = partMother->Pt();
        gmpt = partMother->Pt();
        mpdg = grmaPdgcode;

        // if((int(TMath::Abs(grmaPdgcode)/100.)%10)==kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10)==kBeauty) return kGammaB2M;
        // if((int(TMath::Abs(grmaPdgcode)/100.)%10)==kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10)==kCharm) return kGammaD2M;

        if ((int(TMath::Abs(grmaPdgcode) / 100.) % 10) == kCharm || (int(TMath::Abs(grmaPdgcode) / 1000.) % 10) == kCharm)
        {
          for (Int_t i = 1; i < fgkMaxIter; i++)
          {
            int kLabel = partMother->GetMother();
            if (kLabel == -1)
              return kGammaD2M;
            if ((kLabel < 0) || (kLabel >= fAODArrayMCInfo->GetEntriesFast()))
            {
              AliDebug(1, "Stack label is negative, return\n");
              return -1;
            }
            // if there is an ancester
            if (!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(kLabel)))))
              return -1;
            int ggrandMaPDG = mctrack->GetPdgCode();
            if ((int(TMath::Abs(ggrandMaPDG) / 100.) % 10) == kBeauty || (int(TMath::Abs(ggrandMaPDG) / 1000.) % 10) == kBeauty)
            {
              mpt = mctrack->Pt();
              mpdg = ggrandMaPDG;
              return kGammaB2M;
            }
            partMother = mctrack;
          }
        }
        else if ((int(TMath::Abs(grmaPdgcode) / 100.) % 10) == kBeauty || (int(TMath::Abs(grmaPdgcode) / 1000.) % 10) == kBeauty)
        {
          return kGammaB2M;
        }

        /*tmpMomLabel = partMother->GetMother();
        if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())) {//grandgrandgrandmother
          if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
            partMother = mctrack;
            ggrmaPdgcode = partMother->GetPdgCode();
            mpt = partMother->Pt();
            ggmpt = partMother->Pt();
            mpdg = ggrmaPdgcode;

            if((int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kBeauty) return kGammaB2M;
            if((int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kCharm) return kGammaD2M;
          }
        }//grandgrandgrandmother*/

        if (maPdgcode == 111)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113)
            return kGammaM2M;
          else if (grmaPdgcode == 310)
            return kGammaK0s2P;
          else if (grmaPdgcode == 130)
            return kGammaK0l2P;
          else if (TMath::Abs(grmaPdgcode) == 321)
            return kGammaK2P;
          else if (TMath::Abs(grmaPdgcode) == 3122)
            return kGammaLamda2P;
          else if (grmaPdgcode == 3222)
            return kGammaSigma2P;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kGammaPi0;
        }
        else if (maPdgcode == 221)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113)
            return kGammaM2M;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kGammaEta;
        }
        else if (maPdgcode == 223)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113)
            return kGammaM2M;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kGammaOmega;
        }
        else if (maPdgcode == 333)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113)
            return kGammaM2M;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kGammaPhi;
        }
        else if (maPdgcode == 331)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113)
            return kGammaM2M;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kGammaEtaPrime;
        }
        else if (maPdgcode == 113)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331)
            return kGammaM2M;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kGammaRho0;
        }
        else
          origin = kElse; // grandgrandmother but nothing we identify
      }                   // mctrack grandgrandmother
    }
    else
    {
      // grandmother is primary
      if (maPdgcode == 111)
        return kGammaPi0;
      else if (maPdgcode == 221)
        return kGammaEta;
      else if (maPdgcode == 223)
        return kGammaOmega;
      else if (maPdgcode == 333)
        return kGammaPhi;
      else if (maPdgcode == 331)
        return kGammaEtaPrime;
      else if (maPdgcode == 113)
        return kGammaRho0;
      else
        origin = kElse; // grandmother is primary but nothing we identify
    }
    return origin;
  }
  else
  {
    // check if the ligth meson is the decay product of heavy mesons
    tmpMomLabel = partMotherCopy->GetMother();
    mpt = partMotherCopy->Pt();
    mpdg = partMotherCopy->GetPdgCode();
    if ((tmpMomLabel >= 0) && (tmpMomLabel < fAODArrayMCInfo->GetEntriesFast()))
    { // grandmother
      if ((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel)))))
      {
        partMother = mctrack;
        grmaPdgcode = partMother->GetPdgCode();
        mpt = partMother->Pt();
        gmpt = partMother->Pt();
        mpdg = grmaPdgcode;
        // if((int(TMath::Abs(grmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kBeauty) return kB2M;
        // if((int(TMath::Abs(grmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kCharm) return kD2M;

        if ((int(TMath::Abs(grmaPdgcode) / 100.) % 10) == kCharm || (int(TMath::Abs(grmaPdgcode) / 1000.) % 10) == kCharm)
        {
          for (Int_t i = 1; i < fgkMaxIter; i++)
          {
            int lLabel = partMother->GetMother();
            if (lLabel == -1)
              return kD2M;
            if ((lLabel < 0) || (lLabel >= fAODArrayMCInfo->GetEntriesFast()))
            {
              AliDebug(1, "Stack label is negative, return\n");
              return -1;
            }
            // if there is an ancester
            if (!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(lLabel)))))
              return -1;
            int ggrandMaPDG = mctrack->GetPdgCode();
            if ((int(TMath::Abs(ggrandMaPDG) / 100.) % 10) == kBeauty || (int(TMath::Abs(ggrandMaPDG) / 1000.) % 10) == kBeauty)
            {
              mpt = mctrack->Pt();
              mpdg = ggrandMaPDG;
              return kB2M;
            }
            partMother = mctrack;
          }
        }
        else if ((int(TMath::Abs(grmaPdgcode) / 100.) % 10) == kBeauty || (int(TMath::Abs(grmaPdgcode) / 1000.) % 10) == kBeauty)
        {
          return kB2M;
        }

        /*tmpMomLabel = partMother->GetMother();
        if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())){//grandgrandmother
          if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))){
            partMother = mctrack;
            ggrmaPdgcode = partMother->GetPdgCode();
            mpt = partMother->Pt();
            ggmpt = partMother->Pt();
            mpdg = ggrmaPdgcode;

            if((int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kBeauty) return kB2M;
            if((int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kCharm) return kD2M;
          }
        }//grandgrandmother*/

        if (maPdgcode == 111)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113)
            return kM2M;
          else if (grmaPdgcode == 310)
            return kK0s2P;
          else if (grmaPdgcode == 130)
            return kK0l2P;
          else if (TMath::Abs(grmaPdgcode) == 321)
            return kK2P;
          else if (TMath::Abs(grmaPdgcode) == 3122)
            return kLamda2P;
          else if (grmaPdgcode == 3222)
            return kSigma2P;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kPi0;
        }
        else if (maPdgcode == 221)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113)
            return kM2M;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kEta;
        }
        else if (maPdgcode == 223)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113)
            return kM2M;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kOmega;
        }
        else if (maPdgcode == 333)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113)
            return kM2M;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kPhi;
        }
        else if (maPdgcode == 331)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113)
            return kM2M;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kEtaPrime;
        }
        else if (maPdgcode == 113)
        {
          mpt = gmpt;
          mpdg = grmaPdgcode;
          if (grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331)
            return kM2M;
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kRho0;
        }
        else if (maPdgcode == 321)
        {
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kKe3;
        }
        else if (maPdgcode == 130)
        {
          mpt = partMotherCopy->Pt();
          mpdg = maPdgcode;
          return kK0L;
        }
        else
          origin = kElse; // grandmother but nothing we identidy
      }                   // mctrack grandmother
    }
    else
    {
      // no grandmother
      if (maPdgcode == 111)
        return kPi0;
      else if (maPdgcode == 221)
        return kEta;
      else if (maPdgcode == 223)
        return kOmega;
      else if (maPdgcode == 333)
        return kPhi;
      else if (maPdgcode == 331)
        return kEtaPrime;
      else if (maPdgcode == 113)
        return kRho0;
      else if (maPdgcode == 321)
        return kKe3;
      else if (maPdgcode == 130)
        return kK0L;
      else
        origin = kElse; // mother but nothing we identify
    }
  } // mother is something different from J/psi,charm,beauty or gamma
  return origin;
}

TProfile *AliAnalysisTaskBEpp13TeV::GetEstimatorHistogram()
{
  return fMultEstimatorAvg;
}

int AliAnalysisTaskBEpp13TeV::GetNcharged()
{

  int Nch = 0;

  if (!fIsMC)
    return Nch; // if no MC info return 0

  // loop over all tracks
  for(int igen = 0; igen < fAODArrayMCInfo->GetEntriesFast(); igen++)
  {
    AliAODMCParticle *mctrack = (AliAODMCParticle *)fAODArrayMCInfo->UncheckedAt(igen);
    int charge = mctrack->Charge();
    double eta = mctrack->Eta();
    bool isPhysPrim = mctrack->IsPhysicalPrimary();
    if(charge != 0)
    {
      if(eta > -1.0 && eta < 1.0)
      {
        if(isPhysPrim)
        {
          Nch++;
        }
      }
    }
  }
  return Nch;
}
