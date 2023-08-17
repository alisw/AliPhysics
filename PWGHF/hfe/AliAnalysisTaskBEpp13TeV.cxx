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
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliHFEextraCuts.h"
//#include "AliHFEmcQATest.h"
#include "AliHFEtools.h"
#include "AliAnalysisUtils.h"
#include "TClonesArray.h"
#include "AliAODMCParticle.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "AliHFEV0taginfo.h"
#include "AliAnalysisTaskBEpp13TeV.h"

class AliAnalysisTaskBEpp13TeV;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskBEpp13TeV) // classimp: necessary for root

AliAnalysisTaskBEpp13TeV::AliAnalysisTaskBEpp13TeV()
:AliAnalysisTaskSE()
,fAOD(0)
,fOutputList(0)
,fPIDResponse(0)
,fExtraCuts(0)
//,fMCQA(0)
//,fAODMCHeader(0)
,fAODArrayMCInfo(0)
,fAODMCParticle(0)
,fV0Tagger(0)

,fIsMC(false)
,fMinTPCnCrossedRow(100)
,fMinTPCNclsPID(80)
,fMaxTPCchi2(4)
,fMinTPCclsRatio(0.6)
,fMinITSNcls(3)
,fITSlayer(2)
,fTPCnsigmaLow(-1)
,fTPCnsigmaHigh(3)
,fTOFnsigma(3)
,fMultRef(11.7)

,hVtxZbeforCut(0)
,hVtxZafterCut(0)
,hNrEvents(0)

,hSPDtracklet(0)
,hNtr_vtxZ(0)
,hMultEstimatorAvg(0)
,hSPDtracklet_Corr(0)
,hNtr_vtxZ_Corr(0)
,hMultEstimatorAvg_Corr(0)
,fMultEstimatorAvg(0)

,hFilterMask(0)
,hTPCnCrossedRow(0)
,hTPCclsPID(0)
,hTPCchi2(0)
,hTPCclsRatio(0)
,hITSNcls(0)
,hITSlayer(0)
,hDCAxy(0)
,hDCAz(0)
,hPt(0)
,hEta(0)
,hPhi(0)

,hBhadronPt(0)
,hBhadronPtCorr(0)
,hD0Pt(0)
,hD0PtCorr(0)
,hLcPt(0)
,hPtB(0)
,hPtD(0)
,hPtB2M(0)
,hPtD2M(0)
,hPtGammaB2M(0)
,hPtGammaD2M(0)
,hPtBe(0)
,hPtDe(0)
,hPtB2Me(0)
,hPtD2Me(0)
,hPtGammaB2Me(0)
,hPtGammaD2Me(0)

,hGenBePt(0)
,hRecBePt_track(0)
,hRecBePt_tof(0)
,hRecBePt_tpc(0)

,hITSnsigma(0)
,hITSnsigmaTOFcut(0)
,hITSnsigmaTOFTPCcut(0)
,hTPCnsigma(0)
,hTPCnsigmaTOFcut(0)
,hTPCnsigmaTOFcutPt(0)
,hTPCnsigmaQA(0)
,hTPCnsigmaPiQA(0)
,hTOFnsigma(0)
,hTOFnsigmaQA(0)
,hV0ElecTPCnsigma(0)
,hV0ElecTPCnsigmaTOFcut(0)
,hV0ElecTOFnsigmaDeno(0)
,hV0ElecTOFnsigmaNume(0)

,dcaTrack(0)
,dcaPion(0)
,dcaBeauty(0)
,dcaBeautyCorr(0)
,dcaBeautyCorrVar1(0)
,dcaBeautyCorrVar2(0)
,dcaBzero(0)
,dcaBplus(0)
,dcaBszero(0)
,dcaLb(0)
,DelecVsDmother(0)
,dcaCharm(0)
,dcaDmeson(0)
,dcaDmesonCorr(0)
,dcaDmesonCorrVar1(0)
,dcaDmesonCorrVar2(0)
,dcaDzero(0)
,dcaDplus(0)
,dcaDsplus(0)
,dcaLc(0)
,dcaDalitz(0)
,dcaConv(0)
,dcaB2M(0)
,dcaD2M(0)
,dcaGammaB2M(0)
,dcaGammaD2M(0)

,fBmesonCorrCentLow(0)
,fBmesonCorrCentHigh(0)
,fBmesonCorrMinLow(0)
,fBmesonCorrMinHigh(0)
,fBmesonCorrMaxLow(0)
,fBmesonCorrMaxHigh(0)
,fDmesonCorr(0)
,fDmesonCorrVar1(0)
,fDmesonCorrVar2(0)
,fLcCorr(0)
,fRnd(0)
{
	// default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
	fParentSelect[0][0] =  411; //D+  
  fParentSelect[0][1] =  421; //D0
  fParentSelect[0][2] =  431; //Ds+
  fParentSelect[0][3] = 4122; //Lambdac+
  fParentSelect[0][4] = 4132; //Ksic0
  fParentSelect[0][5] = 4232; //Ksic+
  fParentSelect[0][6] = 4332; //OmegaC0

  fParentSelect[1][0] =  511; //B0
  fParentSelect[1][1] =  521; //B+
  fParentSelect[1][2] =  531; //Bs0
  fParentSelect[1][3] = 5122; //Lambdab0
  fParentSelect[1][4] = 5132; //Ksib-
  fParentSelect[1][5] = 5232; //Ksib0
  fParentSelect[1][6] = 5332; //Omegab-
}
//_____________________________________________________________________________
AliAnalysisTaskBEpp13TeV::AliAnalysisTaskBEpp13TeV(const char* name)
:AliAnalysisTaskSE(name)
,fAOD(0)
,fOutputList(0)
,fPIDResponse(0)
,fExtraCuts(0)
//,fMCQA(0)
//,fAODMCHeader(0)
,fAODArrayMCInfo(0)
,fAODMCParticle(0)
,fV0Tagger(0)

,fIsMC(false)
,fMinTPCnCrossedRow(100)
,fMinTPCNclsPID(80)
,fMaxTPCchi2(4)
,fMinTPCclsRatio(0.6)
,fMinITSNcls(3)
,fITSlayer(2)
,fTPCnsigmaLow(-1)
,fTPCnsigmaHigh(3)
,fTOFnsigma(3)
,fMultRef(11.7)

,hVtxZbeforCut(0)
,hVtxZafterCut(0)
,hNrEvents(0)

,hSPDtracklet(0)
,hNtr_vtxZ(0)
,hMultEstimatorAvg(0)
,hSPDtracklet_Corr(0)
,hNtr_vtxZ_Corr(0)
,hMultEstimatorAvg_Corr(0)
,fMultEstimatorAvg(0)

,hFilterMask(0)
,hTPCnCrossedRow(0)
,hTPCclsPID(0)
,hTPCchi2(0)
,hTPCclsRatio(0)
,hITSNcls(0)
,hITSlayer(0)
,hDCAxy(0)
,hDCAz(0)
,hPt(0)
,hEta(0)
,hPhi(0)

,hBhadronPt(0)
,hBhadronPtCorr(0)
,hD0Pt(0)
,hD0PtCorr(0)
,hLcPt(0)
,hPtB(0)
,hPtD(0)
,hPtB2M(0)
,hPtD2M(0)
,hPtGammaB2M(0)
,hPtGammaD2M(0)
,hPtBe(0)
,hPtDe(0)
,hPtB2Me(0)
,hPtD2Me(0)
,hPtGammaB2Me(0)
,hPtGammaD2Me(0)

,hGenBePt(0)
,hRecBePt_track(0)
,hRecBePt_tof(0)
,hRecBePt_tpc(0)

,hITSnsigma(0)
,hITSnsigmaTOFcut(0)
,hITSnsigmaTOFTPCcut(0)
,hTPCnsigma(0)
,hTPCnsigmaTOFcut(0)
,hTPCnsigmaTOFcutPt(0)
,hTPCnsigmaQA(0)
,hTPCnsigmaPiQA(0)
,hTOFnsigma(0)
,hTOFnsigmaQA(0)
,hV0ElecTPCnsigma(0)
,hV0ElecTPCnsigmaTOFcut(0)
,hV0ElecTOFnsigmaDeno(0)
,hV0ElecTOFnsigmaNume(0)

,dcaTrack(0)
,dcaPion(0)
,dcaBeauty(0)
,dcaBeautyCorr(0)
,dcaBeautyCorrVar1(0)
,dcaBeautyCorrVar2(0)
,dcaBzero(0)
,dcaBplus(0)
,dcaBszero(0)
,dcaLb(0)
,DelecVsDmother(0)
,dcaCharm(0)
,dcaDmeson(0)
,dcaDmesonCorr(0)
,dcaDmesonCorrVar1(0)
,dcaDmesonCorrVar2(0)
,dcaDzero(0)
,dcaDplus(0)
,dcaDsplus(0)
,dcaLc(0)
,dcaDalitz(0)
,dcaConv(0)
,dcaB2M(0)
,dcaD2M(0)
,dcaGammaB2M(0)
,dcaGammaD2M(0)

,fBmesonCorrCentLow(0)
,fBmesonCorrCentHigh(0)
,fBmesonCorrMinLow(0)
,fBmesonCorrMinHigh(0)
,fBmesonCorrMaxLow(0)
,fBmesonCorrMaxHigh(0)
,fDmesonCorr(0)
,fDmesonCorrVar1(0)
,fDmesonCorrVar2(0)
,fLcCorr(0)
,fRnd(0)
{
  fParentSelect[0][0] =  411; //D+  
  fParentSelect[0][1] =  421; //D0
  fParentSelect[0][2] =  431; //Ds+
  fParentSelect[0][3] = 4122; //Lambdac+
  fParentSelect[0][4] = 4132; //Ksic0
  fParentSelect[0][5] = 4232; //Ksic+
  fParentSelect[0][6] = 4332; //OmegaC0

  fParentSelect[1][0] =  511; //B0
  fParentSelect[1][1] =  521; //B+
  fParentSelect[1][2] =  531; //Bs0
  fParentSelect[1][3] = 5122; //Lambdab0
  fParentSelect[1][4] = 5132; //Ksib-
  fParentSelect[1][5] = 5232; //Ksib0
  fParentSelect[1][6] = 5332; //Omegab-
  // constructor
  DefineInput(0, TChain::Class());                    // define the input of the analysis: in this case we take a 'chain' of events
                                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                                        // it does its work automatically
  DefineOutput(1, TList::Class());                    // define the ouptut of the analysis: in this case it's a list of histograms 
  																										// you can add more output objects by calling DefineOutput(2, classname::Class())
                                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                                        // make changes to your AddTask macro!
	fV0Tagger = new AliHFEV0taginfo("Tagger");
}
//_____________________________________________________________________________
AliAnalysisTaskBEpp13TeV::~AliAnalysisTaskBEpp13TeV()
{
  // destructor
  if(fOutputList) delete fOutputList;
  if(fExtraCuts) delete fExtraCuts;
	//if(fMCQA) delete fMCQA;
  if(fV0Tagger) delete fV0Tagger;
}
//_____________________________________________________________________________
void AliAnalysisTaskBEpp13TeV::UserCreateOutputObjects()
{
  fParentSelect[0][0] =  411; //D+  
  fParentSelect[0][1] =  421; //D0
  fParentSelect[0][2] =  431; //Ds+
  fParentSelect[0][3] = 4122; //Lambdac+
  fParentSelect[0][4] = 4132; //Ksic0
  fParentSelect[0][5] = 4232; //Ksic+
  fParentSelect[0][6] = 4332; //OmegaC0

  fParentSelect[1][0] =  511; //B0
  fParentSelect[1][1] =  521; //B+
  fParentSelect[1][2] =  531; //Bs0
  fParentSelect[1][3] = 5122; //Lambdab0
  fParentSelect[1][4] = 5132; //Ksib-
  fParentSelect[1][5] = 5232; //Ksib0
  fParentSelect[1][6] = 5332; //Omegab-
  
	fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
	//fMCQA = new AliHFEmcQATest;
	//if(fMCQA) fMCQA->Init();

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  int nPtBins = 11;
  double ptbinningX[12] = { 1., 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 10. };
  //double ptbinningX[22] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 10., 12. };
  //double ptbinningD0[13] = { 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36. };
  double ptbinningD0[11] = { 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 24. };
  double ptbinningLc[7] = { 1., 2., 4., 6., 8., 12., 24. };
  double ptbinningH[18] = { 0.3, 0.5, 0.75, 1., 1.25, 1.5, 2., 2.5, 3., 4., 5., 6., 7., 8., 10., 12., 16., 20. };

  int nBinsPID = 400;
  double minPID = -10;
  double maxPID = 10;
  double binLimPID[nBinsPID+1];
  for(int i=0; i<=nBinsPID; i++) binLimPID[i] = minPID + (maxPID-minPID)/nBinsPID*(double)i;

  int nBinsIP = 4000;
  double minIP = -0.2;
  double maxIP = 0.2;
  double binLimIP[nBinsIP+1];
  for(int i=0; i<=nBinsIP; i++) binLimIP[i] = minIP + (maxIP-minIP)/nBinsIP*(double)i;
	
  // example of a histogram
  hVtxZbeforCut = new TH1F("hVtxZbeforCut","vertex z",600,-30,30);
  fOutputList->Add(hVtxZbeforCut);
	
  hVtxZafterCut = new TH1F("hVtxZafterCut","vertex z",600,-30,30);
  fOutputList->Add(hVtxZafterCut);
	
  hNrEvents = new TH1F("hNrEvents","number of events",10,0.,10.);
  fOutputList->Add(hNrEvents);

  hSPDtracklet = new TH1F("hSPDtracklet", "SPD tracklet distribution", 300, 0, 300);
  fOutputList->Add(hSPDtracklet);

  hNtr_vtxZ = new TH2F("hNtr_vtxZ", "", 300, -15., 15., 300, 0, 300);
  fOutputList->Add(hNtr_vtxZ);

  hMultEstimatorAvg = new TProfile("hMultEstimatorAvg", "", 300, -15., 15.);
  fOutputList->Add(hMultEstimatorAvg);
  
  hSPDtracklet_Corr = new TH1F("hSPDtracklet_Corr", "SPD tracklet distribution", 300, 0, 300);
  fOutputList->Add(hSPDtracklet_Corr);

  hNtr_vtxZ_Corr = new TH2F("hNtr_vtxZ_Corr", "", 300, -15., 15., 300, 0, 300);
  fOutputList->Add(hNtr_vtxZ_Corr);

  hMultEstimatorAvg_Corr = new TProfile("hMultEstimatorAvg_Corr", "", 300, -15., 15.);
  fOutputList->Add(hMultEstimatorAvg_Corr);
  
  hFilterMask = new TH1F("hFilterMask", "", 2, 0., 2.);
  fOutputList->Add(hFilterMask);
  
  hTPCnCrossedRow = new TH1F("hTPCnCrossedRow", "", 200, 0., 200.);
  fOutputList->Add(hTPCnCrossedRow);
  
  hTPCclsPID = new TH1F("hTPCclsPID", "", 200, 0., 200.);
  fOutputList->Add(hTPCclsPID);
  
  hTPCchi2 = new TH1F("hTPCchi2", "", 100, 0., 10.);
  fOutputList->Add(hTPCchi2);
  
  hTPCclsRatio = new TH1F("hTPCclsRatio", "", 15, 0., 1.5);
  fOutputList->Add(hTPCclsRatio);
  
  hITSNcls= new TH1F("hITSNcls", "", 10, 0., 10.);
  fOutputList->Add(hITSNcls);
  
  hITSlayer = new TH1F("hITSlayer", "", 3, 0.5, 3.5);
  fOutputList->Add(hITSlayer);
  
  hDCAxy = new TH1F("hDCAxy", "", 600, -3., 3.);
  fOutputList->Add(hDCAxy);
  
  hDCAz = new TH1F("hDCAz", "", 600, -3., 3.);
  fOutputList->Add(hDCAz);
  
  hPt = new TH1F("hPt", "pt; (GeV/c)", 300, 0., 30.);
  fOutputList->Add(hPt);
  
  hEta = new TH1F("hEta", "", 200, -1., 1.);
  fOutputList->Add(hEta);
  
  hPhi = new TH1F("hPhi", "", 700, -0.5, 6.5);
  fOutputList->Add(hPhi);
  
  hBhadronPt = new TH1F("hBhadronPt", "", 100, 0, 100);
  fOutputList->Add(hBhadronPt);
  
  hBhadronPtCorr = new TH1F("hBhadronPtCorr", "", 100, 0, 100);
  fOutputList->Add(hBhadronPtCorr);
  
  hD0Pt = new TH1F("hD0Pt", "", 10, ptbinningD0);
  fOutputList->Add(hD0Pt);
  
  hD0PtCorr = new TH1F("hD0PtCorr", "", 10, ptbinningD0);
  fOutputList->Add(hD0PtCorr);
  
  hLcPt = new TH1F("hLcPt", "", 6, ptbinningLc);
  fOutputList->Add(hLcPt);
  
  hPtB = new TH1F("hPtB", "", 17, ptbinningH);
  fOutputList->Add(hPtB);
  
  hPtD = new TH1F("hPtD", "", 17, ptbinningH);
  fOutputList->Add(hPtD);
  
  hPtB2M = new TH1F("hPtB2M", "", 17, ptbinningH);
  fOutputList->Add(hPtB2M);
  
  hPtD2M = new TH1F("hPtD2M", "", 17, ptbinningH);
  fOutputList->Add(hPtD2M);
  
  hPtGammaB2M = new TH1F("hPtGammaB2M", "", 17, ptbinningH);
  fOutputList->Add(hPtGammaB2M);
  
  hPtGammaD2M = new TH1F("hPtGammaD2M", "", 17, ptbinningH);
  fOutputList->Add(hPtGammaD2M);
  
  hPtBe = new TH1F("hPtBe", "", nPtBins, ptbinningX);
  fOutputList->Add(hPtBe);
  
  hPtDe = new TH1F("hPtDe", "", nPtBins, ptbinningX);
  fOutputList->Add(hPtDe);
  
  hPtB2Me = new TH1F("hPtB2Me", "", nPtBins, ptbinningX);
  fOutputList->Add(hPtB2Me);
  
  hPtD2Me = new TH1F("hPtD2Me", "", nPtBins, ptbinningX);
  fOutputList->Add(hPtD2Me);
  
  hPtGammaB2Me = new TH1F("hPtGammaB2Me", "", nPtBins, ptbinningX);
  fOutputList->Add(hPtGammaB2Me);
  
  hPtGammaD2Me = new TH1F("hPtGammaD2Me", "", nPtBins, ptbinningX);
  fOutputList->Add(hPtGammaD2Me);
  
  hGenBePt = new TH1F("hGenBePt", "", nPtBins, ptbinningX);
  fOutputList->Add(hGenBePt);
  
  hRecBePt_track = new TH1F("hRecBePt_track", "", nPtBins, ptbinningX);
  fOutputList->Add(hRecBePt_track);
  
  hRecBePt_tof = new TH1F("hRecBePt_tof", "", nPtBins, ptbinningX);
  fOutputList->Add(hRecBePt_tof);
  
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

  dcaBzero = new TH2F("dcaBzero", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaBzero);

  dcaBplus = new TH2F("dcaBplus", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaBplus);

  dcaBszero = new TH2F("dcaBszero", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaBszero);

  dcaLb = new TH2F("dcaLb", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaLb);
  
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
  
  dcaDzero = new TH2F("dcaDzero", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDzero);
  
  dcaDplus = new TH2F("dcaDplus", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDplus);
  
  dcaDsplus = new TH2F("dcaDsplus", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDsplus);
  
  dcaLc = new TH2F("dcaLc", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaLc);
  
  dcaDalitz = new TH2F("dcaDalitz", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaDalitz);
  
  dcaConv = new TH2F("dcaConv", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaConv);

  dcaB2M = new TH2F("dcaB2M", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaB2M);

  dcaD2M = new TH2F("dcaD2M", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaD2M);

  dcaGammaB2M = new TH2F("dcaGammaB2M", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaGammaB2M);

  dcaGammaD2M = new TH2F("dcaGammaD2M", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
  fOutputList->Add(dcaGammaD2M);

  fRnd = new TRandom3(0);

  PostData(1, fOutputList);				// postdata will notify the analysis manager of changes / updates to the 
  										// fOutputList object. the manager will in the end take care of writing your output to file
										// so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskBEpp13TeV::UserExec(Option_t *){

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  hNrEvents->Fill(0);	// Analyzed events

  if(!fExtraCuts){
    fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  }
  fExtraCuts->SetRecEventInfo(fAOD);

  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse){
    AliDebug(1,"Using default PID Response");
    fPIDResponse = AliHFEtools::GetDefaultPID(false, fInputEvent->IsA()==AliAODEvent::Class());
  }
  
  //Initialize V0 electron tagger
  if(fV0Tagger){
    fV0Tagger->Reset();
    fV0Tagger->TagV0Tracks(fAOD);
  }

  bool isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  if(!isSelected) return;
  hNrEvents->Fill(1);

  // Event selection
  if(!PassEventCuts(fAOD)) return;
  hNrEvents->Fill(2);

  // Pile-up removal
  if(PassPileUpEvent(fAOD)) return;
  hNrEvents->Fill(3);

  //Look for kink mother
  double *fListOfMotherKink = 0;
  int fNumberOfVertices = 0;
  int fNumberOfMotherKink = 0;

  fNumberOfVertices = fAOD->GetNumberOfVertices();
  fListOfMotherKink = new double[fNumberOfVertices];

  for(int iVertex=0; iVertex<fNumberOfVertices; iVertex++){
	AliAODVertex *aodvtx = fAOD->GetVertex(iVertex);
	if(!aodvtx) continue;
	if(aodvtx->GetType()==AliAODVertex::kKink){
	  AliAODTrack *mother = (AliAODTrack*)aodvtx->GetParent();
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
  AliAODTracklets *tracklets = ((AliAODEvent*)fAOD)->GetTracklets();
  int nTracklets = tracklets->GetNumberOfTracklets();
  for(int nn=0; nn<nTracklets;nn++){
    double theta = tracklets->GetTheta(nn);
    double eta = -TMath::Log(TMath::Tan(theta/2.));
    if(TMath::Abs(eta)<1.) nAcceta++;
  }
  hSPDtracklet->Fill(nAcceta);
  hNtr_vtxZ->Fill(vtxZ, nAcceta);
  hMultEstimatorAvg->Fill(vtxZ, nAcceta);

  double Corrected_Ntr = nAcceta;

  TProfile *estimatorAvg = GetEstimatorHistogram();
  if(estimatorAvg)
    Corrected_Ntr = static_cast<int>(GetCorrectedNtracklets(estimatorAvg, nAcceta, vtxZ, fMultRef));

  hSPDtracklet_Corr->Fill(Corrected_Ntr);
  hNtr_vtxZ_Corr->Fill(vtxZ, Corrected_Ntr);
  hMultEstimatorAvg_Corr->Fill(vtxZ, Corrected_Ntr);


  // generated MC loop
  if(fIsMC){
    fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fAODArrayMCInfo){
      AliError("No AOD MC particles");
      return;
    }
    //if(fMCQA) fMCQA->SetMCArray(fAODArrayMCInfo);
	
    for(int iMC = 0; iMC<fAODArrayMCInfo->GetEntries(); iMC++){
      fAODMCParticle = (AliAODMCParticle*) fAODArrayMCInfo->At(iMC);

      int hf = -999;
      double hfpt = -999., hfeta = -999., wghtD =  -999., wghtB = -999.;
      hf = GetHeavyFlavours(fAODMCParticle, hfpt, hfeta);

      if(TMath::Abs(hfeta)<0.5){
        if(hf==kPromptD0){
          hD0Pt->Fill(hfpt);
          wghtD = fDmesonCorr->Eval(hfpt);
          hD0PtCorr->Fill(hfpt, wghtD);
        }
        if(hf==kPromptLc) hLcPt->Fill(hfpt);
      }
      if(TMath::Abs(hfeta<0.8)){
        if(hf==kPromptB || hf==kNonPromptD){
          hBhadronPt->Fill(hfpt);
          if(hfpt<=3.5) wghtB = fBmesonCorrCentLow->Eval(hfpt);
          if(hfpt>3.5) wghtB = fBmesonCorrCentHigh->Eval(hfpt);
          hBhadronPtCorr->Fill(hfpt, wghtB);
        }
      }

      int src = -999, srcPdg = -999;
      double srcPt = -999.;
      src = GetElecSource(fAODMCParticle, true, srcPt, srcPdg);

      if(TMath::Abs(fAODMCParticle->Eta()) < 0.8){
        if(src==kDirectBeauty || src==kBeautyCharm){
          hGenBePt->Fill(fAODMCParticle->Pt());
          hPtB->Fill(srcPt);
          hPtBe->Fill(fAODMCParticle->Pt());
        }
        if(src==kDirectCharm){
          hPtD->Fill(srcPt);
          hPtDe->Fill(fAODMCParticle->Pt());
        }
        if(src==kB2M){
          hPtB2M->Fill(srcPt);
          hPtB2Me->Fill(fAODMCParticle->Pt());
        }
        if(src==kD2M){
          hPtD2M->Fill(srcPt);
          hPtD2Me->Fill(fAODMCParticle->Pt());
        }
        if(src==kGammaB2M){
          hPtGammaB2M->Fill(srcPt);
          hPtGammaB2Me->Fill(fAODMCParticle->Pt());
        }
        if(src==kGammaD2M){
          hPtGammaD2M->Fill(srcPt);
          hPtGammaD2Me->Fill(fAODMCParticle->Pt());
        }
      }
    }
  }

  double fBz=-999.;
  //double fBz2 = fAOD->GetMagneticField() ? 1 : -1;
  if(fAOD->GetMagneticField()<0) fBz = -1.;
  else if(fAOD->GetMagneticField()>0) fBz = 1.;
  else return;

  hNrEvents->Fill(4);
  
  for(int iTracks = 0; iTracks<fAOD->GetNumberOfTracks(); iTracks++){
    AliAODTrack *aodTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
    if(!aodTrack) continue;

	// Reject kink
	bool kinkmotherpass = true;
	for(int kinkmother=0; kinkmother<fNumberOfMotherKink; kinkmother++){
	  if(aodTrack->GetID()==fListOfMotherKink[kinkmother]){
		kinkmotherpass=false;
		continue;
	  }
	}
	if(!kinkmotherpass) continue;
	
    // Track selection
	if(!PassTrackCuts(aodTrack)) continue;
    
	double pt = aodTrack->Pt();
	double hfeImpactParam = -999., hfeImpactParamResol = -999.;
	fExtraCuts->GetHFEImpactParameters((AliVTrack *)aodTrack, hfeImpactParam, hfeImpactParamResol);
	double IP = hfeImpactParam*fBz*aodTrack->Charge();
  
	int mcelectronSource=-999, mcelectronSourcePDG=-999;
	double mcelectronSourcePt=-999.;
	if(fIsMC){
	  fAODMCParticle = NULL;

	  int label = TMath::Abs(aodTrack->GetLabel());
	  if(label < fAODArrayMCInfo->GetEntriesFast())
		fAODMCParticle = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(label));
	  if(fAODMCParticle){
		AliDebug(2, "Associated MC particle found");
		mcelectronSource = GetElecSource(fAODMCParticle, true, mcelectronSourcePt, mcelectronSourcePDG);
	  }
	  // Fill beauty dca information
	  if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm){
        //if(mcelectronSourcePt>70.) continue;
		hRecBePt_track->Fill(pt);
		dcaBeauty->Fill(pt, IP);
		if(TMath::Abs(mcelectronSourcePDG)==511) dcaBzero->Fill(pt, IP);
		
		double wghtB = -99., wghtBvar1 = -99., wghtBvar2 = -99.;
		if(pt>mcelectronSourcePt){
		  wghtB = 1.;
		  wghtBvar1 = 1.;
		  wghtBvar2 = 1.;
		}else{
		  if(mcelectronSourcePt<=3.5){
		    wghtB = fBmesonCorrCentLow->Eval(mcelectronSourcePt);
			wghtBvar1 = fBmesonCorrMinLow->Eval(mcelectronSourcePt);
			wghtBvar2 = fBmesonCorrMaxLow->Eval(mcelectronSourcePt);
		  }else if(mcelectronSourcePt>3.5){
			wghtB = fBmesonCorrCentHigh->Eval(mcelectronSourcePt);
			wghtBvar1 = fBmesonCorrMinHigh->Eval(mcelectronSourcePt);
			wghtBvar2 = fBmesonCorrMaxHigh->Eval(mcelectronSourcePt);
		  }
        }
		
        // B hadron dca correction
		double rndmB = fRnd->Rndm();
    	if(rndmB<wghtB){
		  dcaBeautyCorr->Fill(pt, IP);
		  //if(TMath::Abs(mcelectronSourcePDG)==511) //dcaBzero->Fill(pt, IP);
		  //else if(TMath::Abs(mcelectronSourcePDG)==521) dcaBplus->Fill(pt, IP);
		  //else if(TMath::Abs(mcelectronSourcePDG)==531) dcaBszero->Fill(pt, IP);
		  //else dcaLb->Fill(pt, IP);
		}
		if(rndmB<wghtBvar1) dcaBeautyCorrVar1->Fill(pt, IP);
    	if(rndmB<wghtBvar2) dcaBeautyCorrVar2->Fill(pt, IP);
	  }
	  
	  // Fill charm dca information
	  if(mcelectronSource==kDirectCharm){
		dcaCharm->Fill(pt, IP);
		if(TMath::Abs(mcelectronSourcePDG)==421 || TMath::Abs(mcelectronSourcePDG)==411 || TMath::Abs(mcelectronSourcePDG)==431){
		  DelecVsDmother->Fill(mcelectronSourcePt, pt);
		  dcaDmeson->Fill(pt, IP);
		  double wghtD = -99., wghtDvar1 = -99., wghtDvar2 = -99.;
		  if(pt>mcelectronSourcePt){
			wghtD = 1.;
			wghtDvar1 = 1.;
			wghtDvar2 = 1.;
		  }else{
			wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			wghtDvar1 = fDmesonCorrVar1->Eval(mcelectronSourcePt);
			wghtDvar2 = fDmesonCorrVar2->Eval(mcelectronSourcePt);
			//if(pt>=1. && pt<1.1)  wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=1.1 && pt<1.3) wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=1.3 && pt<1.5) wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=1.5 && pt<2.)  wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=2. && pt<2.5)  wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=2.5 && pt<3.)  wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=3. && pt<4.)   wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=4 && pt<5.)    wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=5. && pt<6.)   wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=6. && pt<8.)   wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=8. && pt<10.)  wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
			//if(pt>=10.)           wghtD = fDmesonCorr->Eval(mcelectronSourcePt);
		  }
		  double rndmD = fRnd->Rndm();
		  if(rndmD<wghtD){
			dcaDmesonCorr->Fill(pt, IP);
			if(TMath::Abs(mcelectronSourcePDG)==421) dcaDzero->Fill(pt, IP);
			if(TMath::Abs(mcelectronSourcePDG)==411) dcaDplus->Fill(pt, IP);
			if(TMath::Abs(mcelectronSourcePDG)==431) dcaDsplus->Fill(pt, IP);
		  }
		  if(rndmD<wghtDvar1) dcaDmesonCorrVar1->Fill(pt, IP);
		  if(rndmD<wghtDvar2) dcaDmesonCorrVar2->Fill(pt, IP);
		}
		if(TMath::Abs(mcelectronSourcePDG)==4122) dcaLc->Fill(pt, IP);
	  }
	  // Fill Dalitz dca information
	  if(mcelectronSource>=kPi0 && mcelectronSource<=kK2P) dcaDalitz->Fill(pt, IP);
	  // Fill conversion dca information
	  if(mcelectronSource>=kGammaPi0 && mcelectronSource<=kGammaK2P) dcaConv->Fill(pt, IP);
	  if(mcelectronSource==kB2M) dcaB2M->Fill(pt, IP);
	  if(mcelectronSource==kD2M) dcaD2M->Fill(pt, IP);
	  if(mcelectronSource==kGammaB2M) dcaGammaB2M->Fill(pt, IP);
	  if(mcelectronSource==kGammaD2M) dcaGammaD2M->Fill(pt, IP);
	}

	// electron identification
	double fITSnSigma = fPIDResponse->NumberOfSigmasITS(aodTrack, AliPID::kElectron);
	double fTPCnSigma = fPIDResponse->NumberOfSigmasTPC(aodTrack, AliPID::kElectron);
	double fTOFnSigma = fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kElectron);
	
	hITSnsigma->Fill(aodTrack->P(), fITSnSigma);
	hTPCnsigma->Fill(aodTrack->P(), fTPCnSigma);
	hTOFnsigma->Fill(aodTrack->P(), fTOFnSigma);
	
	//V0 electrons from systematic studies of TOF eID
	AliPID::EParticleType myv0pid = fV0Tagger->GetV0Info(aodTrack->GetID()); /// enum EParticleType: kElectron = 0, kMuon = 1, kPion = 2, etc
	if(myv0pid == AliPID::kElectron){
	  hV0ElecTPCnsigma->Fill(aodTrack->Pt(), fTPCnSigma);
	  if(TMath::Abs(fTOFnSigma) <= fTOFnsigma) hV0ElecTPCnsigmaTOFcut->Fill(aodTrack->Pt(), fTPCnSigma);
	  // TOF eID systematics
	  if(fTPCnSigma >= fTPCnsigmaLow && fTPCnSigma <= fTPCnsigmaHigh){
		hV0ElecTOFnsigmaDeno->Fill(aodTrack->Pt(), fTOFnSigma);
		if(TMath::Abs(fTOFnSigma) <= fTOFnsigma) hV0ElecTOFnsigmaNume->Fill(aodTrack->Pt(), fTOFnSigma);
	  }
	}
	
	if(TMath::Abs(fTOFnSigma)>fTOFnsigma) continue;
	if(fIsMC){
	  if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm) hRecBePt_tof->Fill(pt);
	}


	hITSnsigmaTOFcut->Fill(aodTrack->P(), fITSnSigma);
	hTPCnsigmaTOFcut->Fill(aodTrack->P(), fTPCnSigma);
	hTPCnsigmaTOFcutPt->Fill(pt, fTPCnSigma);

	if(fTPCnSigma>-5 && fTPCnSigma<-3){
	  hTPCnsigmaPiQA->Fill(aodTrack->P(), fTPCnSigma);
	  dcaPion->Fill(pt, hfeImpactParam*fBz*aodTrack->Charge());
	}

	if(fTPCnSigma<fTPCnsigmaLow || fTPCnSigma>fTPCnsigmaHigh) continue;
	if(fIsMC){
	  if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm) hRecBePt_tpc->Fill(pt);
	}
		
	hITSnsigmaTOFTPCcut->Fill(aodTrack->P(), fITSnSigma);
	hTPCnsigmaQA->Fill(aodTrack->P(), fTPCnSigma);
	hTOFnsigmaQA->Fill(aodTrack->P(), fTOFnSigma);

	dcaTrack->Fill(pt, hfeImpactParam*fBz*aodTrack->Charge());

  }
  
  PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                     // the output manager which will take care of writing
                                                        // it to a file
}
//_____________________________________________________________________________
void AliAnalysisTaskBEpp13TeV::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//__________________________________________________________________
bool AliAnalysisTaskBEpp13TeV::PassEventCuts(AliAODEvent *event){

  //event selection cuts
  const AliAODVertex *vtx = event->GetPrimaryVertex();
  if(!vtx) return false;
  if(vtx->GetNContributors()<2) return false;

  // cut on the primary vertex on z position
  hVtxZbeforCut->Fill(vtx->GetZ());
  if(TMath::Abs(vtx->GetZ()) > 10) return false;
  hVtxZafterCut->Fill(vtx->GetZ());

  return true;
}
//_________________________________________________________________
bool AliAnalysisTaskBEpp13TeV::PassPileUpEvent(AliAODEvent *event){
  //This function checks if there was a pile up reconstructed with SPD
  bool isPileupfromSPDmulbins = event->IsPileupFromSPDInMultBins();
  if(isPileupfromSPDmulbins) return true;
  
  AliAnalysisUtils utils;
  utils.SetMinPlpContribMV(5); //Multi Vertex pileup selection
  utils.SetMaxPlpChi2MV(5); // max value of Chi2perNDF of the pileup and multi-vertex
  utils.SetMinWDistMV(15); // min of the sqrt of weighted distance between the primary and the pileup vertex, multi-vertex
  utils.SetCheckPlpFromDifferentBCMV(false);
  bool isPileupFromMV = utils.IsPileUpMV(event);
  return isPileupFromMV;
}
//_________________________________________________________________
double AliAnalysisTaskBEpp13TeV::GetCorrectedNtracklets(TProfile *estimatorAvg, double rawNtr, double vtxz, double refmult){
  if(TMath::Abs(vtxz)>10) return rawNtr;
  if(!estimatorAvg) return rawNtr;
  double Ntr_mean = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxz));
  double deltaN = rawNtr*(refmult/Ntr_mean - 1);
  double correctedNtr = rawNtr + (deltaN>0 ? 1: -1)*gRandom->Poisson(TMath::Abs(deltaN));

  return correctedNtr;
}
//_________________________________________________________________
bool AliAnalysisTaskBEpp13TeV::PassTrackCuts(AliAODTrack *track){

  if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;
  if(TMath::Abs(track->Eta()) > 0.8) return false;
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

  double FoundOverFindable = (findableTPC ? static_cast<float>(TPCnCrossedRow)/static_cast<float>(findableTPC) : 0);
  if(TPCnCrossedRow < fMinTPCnCrossedRow || TPCclsPID < fMinTPCNclsPID || TPCchi2 > fMaxTPCchi2 || FoundOverFindable < fMinTPCclsRatio) return false;

  // ITS cut
  int ITSnCls = track->GetITSNcls();
  if(ITSnCls < fMinITSNcls) return false;
  if(fITSlayer==AliHFEextraCuts::kFirst){
    if(!(track->HasPointOnITSLayer(0))) return false;
    hITSlayer->Fill(1);
  }else if(fITSlayer==AliHFEextraCuts::kAny){
    if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return false;
    hITSlayer->Fill(2);
  }else if(fITSlayer==AliHFEextraCuts::kBoth){
    if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return false;
    hITSlayer->Fill(3);
  }else return false;
	
  // dca cut
  float dcaxy = -999.; float dcaz = -999.;
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
int AliAnalysisTaskBEpp13TeV::GetElecSource(const AliAODMCParticle * const mcpart, double &mpt, int &mpdg){
	//printf("Modified version of source selection \n");
  if(!mcpart) return kMisID;
  if(!fAODArrayMCInfo) return -1;
	
  if(TMath::Abs(mcpart->GetPdgCode()) != 11 ) return kElse;

  int origin = -1;
  bool isFinalOpenCharm = kFALSE;

  int iLabel = mcpart->GetMother();
  if((iLabel<0) || (iLabel>=fAODArrayMCInfo->GetEntriesFast())){
	AliDebug(1, "label is out of range, return\n");
	return -1;
  }
	
  AliAODMCParticle *mctrack = NULL; // will change all the time
  int tmpMomLabel=0;
  if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(iLabel))))) return -1;
  AliAODMCParticle *partMother = mctrack;	//mtrack 
  AliAODMCParticle *partMotherCopy = mctrack;	//mtrack
  int maPdgcode = partMother->GetPdgCode();	//mpdg
  mpt = partMother->Pt();	//mpt
  mpdg = partMother->GetPdgCode();	//mpdg
  int gmaPdgcode, ggmaPdgcode;
  double gmpt, ggmpt;
  int gmpdg, ggmpdg;

  // if the mother is charmed hadron
  if((int(TMath::Abs(maPdgcode)/100.)%10)==4 || (int(TMath::Abs(maPdgcode)/1000.)%10)==4){
	if(TMath::Abs(maPdgcode)==411 || TMath::Abs(maPdgcode)==421 || TMath::Abs(maPdgcode)==431 || TMath::Abs(maPdgcode)==4122 || TMath::Abs(maPdgcode)==4132 || TMath::Abs(maPdgcode)==4232 || TMath::Abs(maPdgcode)==4332){
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
	  if(jLabel == -1){
		return kDirectCharm;
	  }
	  if(jLabel<0 || jLabel>=fAODArrayMCInfo->GetEntriesFast()){
		AliDebug(1, "Stack label is negative, return\n");
		return -1;
	  }
			
	  if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel))))){
		return -1;
	  }
	  int grandMaPDG = mctrack->GetPdgCode();
	  if(TMath::Abs(grandMaPDG)==511 || TMath::Abs(grandMaPDG)==521 || TMath::Abs(grandMaPDG)==531 || TMath::Abs(grandMaPDG)==5122 || TMath::Abs(grandMaPDG)==5132 || TMath::Abs(grandMaPDG)==5232 || TMath::Abs(grandMaPDG)==5332){
		mpt = mctrack->Pt();
		mpdg = mctrack->GetPdgCode();
		return kBeautyCharm;
	  }
	  partMother = mctrack;
	} // end of iteration 
  }
  
  // if the mother is beauty hadron
  else if((int(TMath::Abs(maPdgcode)/100.)%10)==5 || (int(TMath::Abs(maPdgcode)/1000.)%10)==5){
	if(TMath::Abs(maPdgcode)==511 || TMath::Abs(maPdgcode)==521 || TMath::Abs(maPdgcode)==531 || TMath::Abs(maPdgcode)==5122 || TMath::Abs(maPdgcode)==5132 || TMath::Abs(maPdgcode)==5232 || TMath::Abs(maPdgcode)==5332){
	  mpt = partMotherCopy->Pt();
	  mpdg = partMotherCopy->GetPdgCode();
	  return kDirectBeauty;
	}
  }
	
  // if the mother is gamma
  else if(TMath::Abs(maPdgcode)==22){
	tmpMomLabel = partMotherCopy->GetMother();  // mother of photon
	mpt = partMotherCopy->Pt(); // pT of photon
	mpdg = partMotherCopy->GetPdgCode();
	if(tmpMomLabel==-1) return kGamma;  // no grandmother
	if((tmpMomLabel<0) || (tmpMomLabel>=fAODArrayMCInfo->GetEntriesFast())) {
	  return -1;
	}
	if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
	  return -1;
	}
	partMother = mctrack; // gmtrack
	partMotherCopy = mctrack; // gmtrack
	mpt = partMother->Pt(); // grand mother pT
	mpdg = partMother->GetPdgCode(); // grand mother PDG
	maPdgcode = partMother->GetPdgCode(); // grand mother PDG
		
	// check if the ligth meson is the decay product of heavy mesons
	tmpMomLabel = partMother->GetMother(); // grand grand mother of photon
	if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())){//grand grand mother
	  if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))){
		partMother = mctrack; //ggmtrack
        gmaPdgcode = partMother->GetPdgCode(); //grand grand mother PDG
		mpt = partMother->Pt(); // grand grand mother pT
		mpdg = partMother->GetPdgCode(); // grand grand mother pT
		gmpt = partMother->Pt(); // grand grand mother pt
		gmpdg = partMother->GetPdgCode(); // grand grand mother pT

		if(TMath::Abs(maPdgcode)==111){
		  mpt = gmpt;
		  mpdg = gmpdg;
		  if(gmaPdgcode == 310) return kGammaK0s2P;
		  else if(gmaPdgcode == 130) return kGammaK0l2P;
		  else if(TMath::Abs(gmaPdgcode) == 321) return kGammaK2P;
		  else if(TMath::Abs(gmaPdgcode) == 3122) return kGammaLamda2P;
		  else if(gmaPdgcode == 3222) return kGammaSigma2P;
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaPi0;
		}
		else if(TMath::Abs(maPdgcode)==221){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaEta;
		}
		else if(TMath::Abs(maPdgcode)==223){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaOmega;
		}
		else if(TMath::Abs(maPdgcode)==333){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaPhi;
		}
		else if(TMath::Abs(maPdgcode)==331){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaEtaPrime;
		}
		else if(TMath::Abs(maPdgcode)==113){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kGammaRho0;
		}
		else origin = kElse;//grand grand mother but nothing we identify
	  }//mctrack grandgrandmother
	}
	else{
	  // grandmother is primary
	  if(TMath::Abs(maPdgcode)==111){
		return kGammaPi0;
	  }
	  else if(TMath::Abs(maPdgcode)==221){
		return kGammaEta;
	  }
	  else if(TMath::Abs(maPdgcode)==223){
		return kGammaOmega;
	  }
	  else if(TMath::Abs(maPdgcode)==333){
		return kGammaPhi;
	  }
	  else if(TMath::Abs(maPdgcode)==331){
		return kGammaEtaPrime;
	  }
	  else if(TMath::Abs(maPdgcode)==113){
		return kGammaRho0;
	  }
	  else origin = kElse;//grandmother is primary but nothing we identify
	}
	return origin;
  }

  // if the mother is light meson
  else{
	
	tmpMomLabel = partMotherCopy->GetMother(); // grand mother
	mpt = partMotherCopy->Pt(); // mother pT
	mpdg = partMotherCopy->GetPdgCode(); // mother PDG
	if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())){// grand mother
	  if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))){
		partMother = mctrack; // grand mother
		gmaPdgcode = partMother->GetPdgCode(); // grand mother PDG
		mpt = partMother->Pt(); // grand mother pT
		mpdg = partMother->GetPdgCode(); // grand mother PDG
		gmpt = partMother->Pt(); // grand mother pT
		gmpdg = partMother->GetPdgCode(); // grand mother PDG

		if(TMath::Abs(maPdgcode)==111){
		  mpt = gmpt;
		  mpdg = gmpdg;
		  if(gmaPdgcode == 310) return kK0s2P;
		  else if(gmaPdgcode == 130) return kK0l2P;
		  else if(TMath::Abs(gmaPdgcode) == 321) return kK2P;
		  else if(TMath::Abs(gmaPdgcode) == 3122) return kLamda2P;
		  else if(gmaPdgcode == 3222) return kSigma2P;
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kPi0;
		}
		else if(TMath::Abs(maPdgcode)==221){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kEta;
		}
		else if(TMath::Abs(maPdgcode)==223){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kOmega;
		}
		else if(TMath::Abs(maPdgcode)==333){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kPhi;
		}
		else if(TMath::Abs(maPdgcode)==331){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kEtaPrime;
		}
		else if(TMath::Abs(maPdgcode)==113){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kRho0;
		}
		else if(TMath::Abs(maPdgcode)==321){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kKe3;
		}
		else if(TMath::Abs(maPdgcode)==130){
		  mpt = partMotherCopy->Pt();
		  mpdg = partMotherCopy->GetPdgCode();
		  return kK0L;
		}
		else origin = kElse;//grandmother but nothing we identidy
	  }//mctrack grandmother
	}
	else {
	  // no grandmother
	  if(TMath::Abs(maPdgcode)==111) return kPi0;
	  else if(TMath::Abs(maPdgcode)==221) return kEta;
	  else if(TMath::Abs(maPdgcode)==223) return kOmega;
	  else if(TMath::Abs(maPdgcode)==333) return kPhi;
	  else if(TMath::Abs(maPdgcode)==331) return kEtaPrime;
	  else if(TMath::Abs(maPdgcode)==113) return kRho0;
	  else if(TMath::Abs(maPdgcode)==321) return kKe3;
	  else if(TMath::Abs(maPdgcode)==130) return kK0L;
	  else origin = kElse;//mother but nothing we identify
	}
  }//mother is something different from J/psi,charm,beauty or gamma
	
  return origin;
}
//_______________________________________________________________________________________________________________
int AliAnalysisTaskBEpp13TeV::GetHeavyFlavours(const AliAODMCParticle * const mcpart, double &hfpt, double &hfeta){

  if(!mcpart) return -1;
  if(!fAODArrayMCInfo) return -1;
  
  int pdgHF = TMath::Abs(mcpart->GetPdgCode());
  hfpt = mcpart->Pt();
  hfeta = mcpart->Eta();
  if(!(pdgHF/100==4 || pdgHF/100==5 || pdgHF/1000==4 || pdgHF/1000==5)) return -1;

  AliAODMCParticle *mctrack = NULL;
  AliAODMCParticle *partMother = NULL;
  
  if(pdgHF==411 || pdgHF==421 || pdgHF==431 || pdgHF==4122 || pdgHF==4132 || pdgHF==4232 || pdgHF==4332){
    // iterate until find B hadron as a mother
    int jLabel = -999;
    int maPdgcode = -999;
    for(int i=1; i<100; i++){
      if(i==1) jLabel = mcpart->GetMother();
      if(i!=1) jLabel = partMother->GetMother();
     
      if(jLabel==-1){
        if(pdgHF==421) return kPromptD0;
        if(pdgHF==4122) return kPromptLc;
      }    
      if(jLabel<0 || jLabel>=fAODArrayMCInfo->GetEntriesFast()){
        AliDebug(1, "Stack label is negative, return\n");
        return -1;
      }    
      if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel))))) {
        return -1;
      }    
      maPdgcode = TMath::Abs(mctrack->GetPdgCode());
      if(maPdgcode==511 || maPdgcode==521 || maPdgcode==531 || maPdgcode==5122 || maPdgcode==5132 || maPdgcode==5232 || maPdgcode==5332){
        hfpt = mctrack->Pt();
        hfeta = mctrack->Eta();
        return kNonPromptD;
      }    
      partMother = mctrack;
    }// end of iteration 
  }
  
  // prompt B mesons
  else if(pdgHF==511 || pdgHF==521 || pdgHF==531 || pdgHF==5122 || pdgHF==5132 || pdgHF==5232 || pdgHF==5332){
    return kPromptB;
  }

  return -1;
}
//_______________________________________________________________________________________________________________
int AliAnalysisTaskBEpp13TeV::GetElecSource(const AliAODMCParticle * const mcpart, bool isElec, double &mpt, int &mpdg) const
{
	//printf("Original version of source selection \n");

  if (!mcpart) return -1;
  if (!fAODArrayMCInfo) return -1;

  //if(isElec) if(TMath::Abs(mcpart->GetPdgCode()) != AliHFEmcQATest::kElectronPDG) return kMisID;
  if(isElec) if(TMath::Abs(mcpart->GetPdgCode()) != 11) return kMisID;

  int origin = -1;
  bool isFinalOpenCharm = kFALSE;

  int iLabel = mcpart->GetMother();
  if((iLabel<0) || (iLabel>=fAODArrayMCInfo->GetEntriesFast())){
	AliDebug(1, "label is out of range, return\n");
    return -1;
  }

  AliAODMCParticle *mctrack = NULL; // will change all the time
  int tmpMomLabel=0;
  if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(iLabel))))) return -1; 
  AliAODMCParticle *partMother = mctrack; 
  AliAODMCParticle *partMotherCopy = mctrack;
  int maPdgcode = partMother->GetPdgCode();
  mpt = partMother->Pt();
  mpdg = maPdgcode;
  int grmaPdgcode;
  int ggrmaPdgcode;
  double gmpt, ggmpt;

  // if the mother is charmed hadron  
  if(TMath::Abs(maPdgcode)==443){ 
    // J/spi
    int jLabel = partMother->GetMother();
    if((jLabel>=0) && (jLabel<fAODArrayMCInfo->GetEntriesFast())){
      if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel))))){
				int grandMaPDG = mctrack->GetPdgCode();
				mpt = mctrack->Pt();
				mpdg = grandMaPDG;
				if((int(TMath::Abs(grandMaPDG)/100.)%10) == kBeauty || (int(TMath::Abs(grandMaPDG)/1000.)%10) == kBeauty) return kB2Jpsi;
      }
    }
    return kJpsi;   
  }else if((int(TMath::Abs(maPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(maPdgcode)/1000.)%10) == kCharm){
    // charm
    for(int i=0; i<fNparents; i++){
      if(TMath::Abs(maPdgcode)==fParentSelect[0][i]){
				mpt = partMother->Pt();
				mpdg = maPdgcode;
				isFinalOpenCharm = kTRUE;
      }
    }
    if(!isFinalOpenCharm) return -1;
    
    // iterate until you find B hadron as a mother or become top ancester 
    for(int i=1; i<fgkMaxIter; i++){
      
      int jLabel = partMother->GetMother();
      if(jLabel == -1) return kDirectCharm;

      if((jLabel<0) || (jLabel>=fAODArrayMCInfo->GetEntriesFast())){
				AliDebug(1, "Stack label is negative, return\n");
				return -1;
      }
      
      // if there is an ancester
      if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel))))) return -1;
      int grandMaPDG = mctrack->GetPdgCode();
      for(int j=0; j<fNparents; j++){
				if(TMath::Abs(grandMaPDG)==fParentSelect[1][j]){
		  		mpt = mctrack->Pt();
		  		mpdg = grandMaPDG;
		  		return kBeautyCharm;
				}
      }
      partMother = mctrack;
    } // end of iteration 

  } // end of charm
  else if((int(TMath::Abs(maPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(maPdgcode)/1000.)%10) == kBeauty){
    // beauty
    for(int i=0; i<fNparents; i++){
      if(TMath::Abs(maPdgcode)==fParentSelect[1][i]){
				mpt = partMotherCopy->Pt();
				mpdg = maPdgcode;
				return kDirectBeauty;
      }
    }
  } // end of beauty
  else if(TMath::Abs(maPdgcode) == 22){ 
    //conversion
    tmpMomLabel = partMotherCopy->GetMother();
		mpt = partMotherCopy->Pt(); // pT of photon
    mpdg = maPdgcode;
		if(tmpMomLabel==-1) return kGamma;
    if((tmpMomLabel<0) || (tmpMomLabel>=fAODArrayMCInfo->GetEntriesFast())) return -1;
    if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) return -1;

    partMother = mctrack;
    partMotherCopy = mctrack;
		mpt = partMother->Pt();
    maPdgcode = partMother->GetPdgCode();
		mpdg = maPdgcode;
    
    // check if the ligth meson is the decay product of heavy mesons
    tmpMomLabel = partMother->GetMother();
    if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())){//grandgrandmother
      if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))){
				partMother = mctrack;
				grmaPdgcode = partMother->GetPdgCode();
				mpt = partMother->Pt();
				gmpt = partMother->Pt();
				mpdg = grmaPdgcode;

				//if((int(TMath::Abs(grmaPdgcode)/100.)%10)==kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10)==kBeauty) return kGammaB2M;
				//if((int(TMath::Abs(grmaPdgcode)/100.)%10)==kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10)==kCharm) return kGammaD2M;
	
				if((int(TMath::Abs(grmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kCharm){
					for (Int_t i=1; i<fgkMaxIter; i++){
      			int kLabel = partMother->GetMother();
      			if(kLabel == -1) return kGammaD2M;
      			if((kLabel<0) || (kLabel>=fAODArrayMCInfo->GetEntriesFast())){
							AliDebug(1, "Stack label is negative, return\n");
							return -1;
      			}
      			// if there is an ancester
      			if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(kLabel))))) return -1;
      			int ggrandMaPDG = mctrack->GetPdgCode();
						if((int(TMath::Abs(ggrandMaPDG)/100.)%10) == kBeauty || (int(TMath::Abs(ggrandMaPDG)/1000.)%10) == kBeauty){
		  				mpt = mctrack->Pt();
		  				mpdg = ggrandMaPDG;
							return kGammaB2M;
						}
      			partMother = mctrack;
					}
				}else if((int(TMath::Abs(grmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kBeauty){
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
	
				if(TMath::Abs(maPdgcode) == 111){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
		  		else if(grmaPdgcode == 310) return kGammaK0s2P;
		  		else if(grmaPdgcode == 130) return kGammaK0l2P;
		  		else if(TMath::Abs(grmaPdgcode) == 321) return kGammaK2P;
		  		else if(TMath::Abs(grmaPdgcode) == 3122) return kGammaLamda2P;
		  		else if(grmaPdgcode == 3222) return kGammaSigma2P;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kGammaPi0;
				}else if(TMath::Abs(maPdgcode) == 221){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kGammaEta;
				}else if(TMath::Abs(maPdgcode) == 223){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kGammaOmega;
				}else if(TMath::Abs(maPdgcode) == 333){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kGammaPhi;
				}else if(TMath::Abs(maPdgcode) == 331){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kGammaM2M;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kGammaEtaPrime; 
				}else if(TMath::Abs(maPdgcode) == 113){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kGammaM2M;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kGammaRho0;
				}else origin = kElse;//grandgrandmother but nothing we identify
	  	}//mctrack grandgrandmother
    }else{
      // grandmother is primary
      if(TMath::Abs(maPdgcode) == 111) return kGammaPi0;
      else if ( TMath::Abs(maPdgcode) == 221 ) return kGammaEta;
      else if ( TMath::Abs(maPdgcode) == 223 ) return kGammaOmega;
      else if ( TMath::Abs(maPdgcode) == 333 ) return kGammaPhi;
      else if ( TMath::Abs(maPdgcode) == 331 ) return kGammaEtaPrime; 
      else if ( TMath::Abs(maPdgcode) == 113 ) return kGammaRho0;
      else origin = kElse;//grandmother is primary but nothing we identify
		}
    return origin;
  }else{
    // check if the ligth meson is the decay product of heavy mesons
    tmpMomLabel = partMotherCopy->GetMother();
		mpt = partMotherCopy->Pt();
		mpdg = partMotherCopy->GetPdgCode();
    if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())){//grandmother
      if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))){
				partMother = mctrack;
				grmaPdgcode = partMother->GetPdgCode();
				mpt = partMother->Pt();
				gmpt = partMother->Pt();
				mpdg = grmaPdgcode;
				//if((int(TMath::Abs(grmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kBeauty) return kB2M;
				//if((int(TMath::Abs(grmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kCharm) return kD2M;
	
				if((int(TMath::Abs(grmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kCharm){
					for (Int_t i=1; i<fgkMaxIter; i++){
      			int lLabel = partMother->GetMother();
      			if(lLabel == -1) return kD2M;
      			if((lLabel<0) || (lLabel>=fAODArrayMCInfo->GetEntriesFast())){
							AliDebug(1, "Stack label is negative, return\n");
							return -1;
      			}
      			// if there is an ancester
      			if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(lLabel))))) return -1;
      			int ggrandMaPDG = mctrack->GetPdgCode();
						if((int(TMath::Abs(ggrandMaPDG)/100.)%10) == kBeauty || (int(TMath::Abs(ggrandMaPDG)/1000.)%10) == kBeauty){
		  				mpt = mctrack->Pt();
		  				mpdg = ggrandMaPDG;
							return kB2M;
						}
      			partMother = mctrack;
					}
				}else if((int(TMath::Abs(grmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kBeauty){
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
	
				if(TMath::Abs(maPdgcode) == 111){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
		  		else if(grmaPdgcode == 310) return kK0s2P;
		  		else if(grmaPdgcode == 130) return kK0l2P;
		  		else if(TMath::Abs(grmaPdgcode) == 321) return kK2P;
		  		else if(TMath::Abs(grmaPdgcode) == 3122) return kLamda2P;
		  		else if(grmaPdgcode == 3222) return kSigma2P;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kPi0;
				}else if(TMath::Abs(maPdgcode) == 221){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kEta;
				}else if(TMath::Abs(maPdgcode) == 223){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kOmega;
				}else if(TMath::Abs(maPdgcode) == 333){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kPhi;
				}else if(TMath::Abs(maPdgcode) == 331){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kM2M;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kEtaPrime;
				}else if(TMath::Abs(maPdgcode) == 113){
		  		mpt = gmpt;
		  		mpdg = grmaPdgcode;
		  		if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kM2M;
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kRho0;
				}else if(TMath::Abs(maPdgcode) == 321){
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kKe3;
				}else if(TMath::Abs(maPdgcode) == 130){
		  		mpt = partMotherCopy->Pt();
		  		mpdg = maPdgcode;
		  		return kK0L;
				}else origin = kElse;//grandmother but nothing we identidy
	  	}//mctrack grandmother
    }else{
      // no grandmother
      if ( TMath::Abs(maPdgcode) == 111 ) return kPi0;
      else if ( TMath::Abs(maPdgcode) == 221 ) return kEta;
      else if ( TMath::Abs(maPdgcode) == 223 ) return kOmega;
      else if ( TMath::Abs(maPdgcode) == 333 ) return kPhi;
      else if ( TMath::Abs(maPdgcode) == 331 ) return kEtaPrime;
      else if ( TMath::Abs(maPdgcode) == 113 ) return kRho0;
      else if ( TMath::Abs(maPdgcode) == 321 ) return kKe3;
      else if ( TMath::Abs(maPdgcode) == 130 ) return kK0L;
      else origin = kElse;//mother but nothing we identify
    }
  }//mother is something different from J/psi,charm,beauty or gamma
  return origin;
}

TProfile *AliAnalysisTaskBEpp13TeV::GetEstimatorHistogram(){
  return fMultEstimatorAvg;
}
