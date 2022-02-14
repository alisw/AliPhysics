/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *               
 *                                                                        *               
 * Author: The ALICE Off-line Project.                                    *               
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

//////////////////////////////////////////////////////////////////
//                                         			                //
//    Task for beauty electron analysis in Pb-Pb collisions			//
//																															//
//	  Authors																										//
//		Jonghan Park (jonghan@cern.ch)														//
//																															//
//////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskBtoElecPbPbTPCTOF.h"
#include "AliAODMCParticle.h"
#include "AliPIDResponse.h"
#include "AliHFEextraCuts.h"
#include "AliHFEtools.h"
#include "AliMultSelection.h"
#include "AliGenEventHeader.h"
#include "TAxis.h"
#include "AliKFParticle.h"
#include "TRandom3.h"
#include "THnSparse.h"
#include "AliAODv0.h"
#include "AliAODv0KineCuts.h"
#include "AliESDtrack.h"
#include "AliHFEV0taginfo.h"

//______________________________________________________________________
ClassImp(AliAnalysisTaskBtoElecPbPbTPCTOF)

//______________________________________________________________________
AliAnalysisTaskBtoElecPbPbTPCTOF::AliAnalysisTaskBtoElecPbPbTPCTOF(const char *name) 
  : AliAnalysisTaskSE(name)
,fIsMC(0)
,fCentrality(-999.)
,fCentralityCalib(-999.)
,fCentMin(30.)
,fCentMax(50.)
,fMinTPCNcls(100)
,fMinTPCNclsPID(80)
,fMaxTPCchi2(4.)
,fMinTPCclsRatio(0.6)
,fMinITSNcls(4)
,fMaxITSclsFrac(0.3)
,fMaxITSchi2(10.)
,fITSlayer("kBoth")
,fEta(0.8)
,fMinPt(0.5)
,fMaxPt(30.)
,fDCAxy(0.1)
,fDCAz(0.2)
,fTPCnsigmaLow(0.)
,fTPCnsigmaHigh(3.)
,fTOFnsigma(3.)

,fAOD(0)
,fOutputList(0)
,fPidResponse(0)
,fAODMCHeader(0)
,fAODArrayMCInfo(0)
,fExtraCuts(0)
,fAODMCParticle(0)
,fAODv0(0)
,fAODV0Cuts(0)
,fV0Tagger(0)

,hCent_nocut(0)
,hCent_nocut2(0)
,hCent_cut(0)
,hNrEvents(0)

,hFilterMask(0)
,hTPCNcls(0)
,hTPCclsPID(0)
,hTPCchi2(0)
,hTPCclsRatio(0)
,hITSNcls(0)
,hITSclsFrac(0)
,hITSchi2(0)
,hITSlayer(0)
,hDCAxy(0)
,hDCAz(0)
,hPt(0)
,hEta(0)
,hPhi(0)

,hITSnsigma(0)
,hITSnsigmaTOFcut(0)
,hITSnsigmaQA(0)
,hTPCnsigma(0)
,hTPCnsigmaTOFcut(0)
,hTPCnsigmaITSTOFcut(0)
,hTPCnsigmaQA(0)
,hTPCnsigmaPiQA(0)
,hTOFnsigma(0)
,hTOFnsigmaQA(0)

,hV0ElecTOFnsigmaDeno(0)
,hV0ElecTOFnsigmaNume(0)
,hV0ElecTPCnsigmaDeno(0)
,hV0ElecTPCnsigmaNume(0)

,hGenBtoElecPt(0)
,hRecBtoElecPt_nocut(0)
,hRecBtoElecPt_track(0)
,hRecBtoElecPt_tof(0)
,hRecBtoElecPt_tpc(0)

,hD0Pt(0)
,hD0PtCorr(0)

,hBhadronPt(0)
,hBhadronPtCorr(0)

,hLcPt(0)
,hLcPtCorr(0)

,dcaTrack(0)
,dcaDmeson(0)
,dcaDmesonCorr(0)
,dcaDzero(0)
,dcaDplus(0)
,dcaDsplus(0)
,dcaLc(0)
,dcaBeauty(0)
,dcaBeautyCorr(0)
,dcaDalitz(0)
,dcaConv(0)
,dcaPion1(0)
,dcaPion2(0)

,hProdR_Pt(0)
,hV0PionMinCut(0)
,hV0PionMaxCut(0)
,hV0PionMult(0)

,fBmesonCorr(0)
,fBmesonCorr1(0)
,fDmesonCorr(0)
,fDmesonCorr1(0)
,fDmesonCorr2(0)
,fDmesonCorr3(0)
,fDmesonCorr4(0)
,fDmesonCorr5(0)
,fDmesonCorr6(0)
,fDmesonCorr7(0)
,fDmesonCorr8(0)
,fDmesonCorr9(0)
,fDmesonCorr10(0)
,fDmesonCorr11(0)
,fDmesonCorr12(0)
,fLcCorr(0)
,fD0TauWeight(0)
,fDpTauWeight(0)
,fDsTauWeight(0)
,fB0TauWeight(0)
,fBpTauWeight(0)
,fBsTauWeight(0)
,fRnd(0)
{
  //Named constructor
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
	fV0Tagger = new AliHFEV0taginfo("Tagger");
}

//________________________________________________________________________
AliAnalysisTaskBtoElecPbPbTPCTOF::AliAnalysisTaskBtoElecPbPbTPCTOF() 
  : AliAnalysisTaskSE()

,fIsMC(0)
,fCentrality(-999.)
,fCentralityCalib(-999.)
,fCentMin(30.)
,fCentMax(50.)
,fMinTPCNcls(100)
,fMinTPCNclsPID(80)
,fMaxTPCchi2(4.)
,fMinTPCclsRatio(0.6)
,fMinITSNcls(4)
,fMaxITSclsFrac(0.3)
,fMaxITSchi2(10.)
,fITSlayer("kBoth")
,fEta(0.8)
,fMinPt(0.5)
,fMaxPt(30.)
,fDCAxy(0.1)
,fDCAz(0.2)
,fTPCnsigmaLow(0.)
,fTPCnsigmaHigh(3.)
,fTOFnsigma(3.)

,fAOD(0)
,fOutputList(0)
,fPidResponse(0)
,fAODMCHeader(0)
,fAODArrayMCInfo(0)
,fExtraCuts(0)
,fAODMCParticle(0)
,fAODv0(0)
,fAODV0Cuts(0)
,fV0Tagger(0)

,hCent_nocut(0)
,hCent_nocut2(0)
,hCent_cut(0)
,hNrEvents(0)

,hFilterMask(0)
,hTPCNcls(0)
,hTPCclsPID(0)
,hTPCchi2(0)
,hTPCclsRatio(0)
,hITSNcls(0)
,hITSclsFrac(0)
,hITSchi2(0)
,hITSlayer(0)
,hDCAxy(0)
,hDCAz(0)
,hPt(0)
,hEta(0)
,hPhi(0)

,hITSnsigma(0)
,hITSnsigmaTOFcut(0)
,hITSnsigmaQA(0)
,hTPCnsigma(0)
,hTPCnsigmaTOFcut(0)
,hTPCnsigmaITSTOFcut(0)
,hTPCnsigmaQA(0)
,hTPCnsigmaPiQA(0)
,hTOFnsigma(0)
,hTOFnsigmaQA(0)

,hV0ElecTOFnsigmaDeno(0)
,hV0ElecTOFnsigmaNume(0)
,hV0ElecTPCnsigmaDeno(0)
,hV0ElecTPCnsigmaNume(0)

,hGenBtoElecPt(0)
,hRecBtoElecPt_nocut(0)
,hRecBtoElecPt_track(0)
,hRecBtoElecPt_tof(0)
,hRecBtoElecPt_tpc(0)

,hD0Pt(0)
,hD0PtCorr(0)

,hBhadronPt(0)
,hBhadronPtCorr(0)

,hLcPt(0)
,hLcPtCorr(0)

,dcaTrack(0)
,dcaDmeson(0)
,dcaDmesonCorr(0)
,dcaDzero(0)
,dcaDplus(0)
,dcaDsplus(0)
,dcaLc(0)
,dcaBeauty(0)
,dcaBeautyCorr(0)
,dcaDalitz(0)
,dcaConv(0)
,dcaPion1(0)
,dcaPion2(0)

,hProdR_Pt(0)
,hV0PionMinCut(0)
,hV0PionMaxCut(0)
,hV0PionMult(0)

,fBmesonCorr(0)
,fBmesonCorr1(0)
,fDmesonCorr(0)
,fDmesonCorr1(0)
,fDmesonCorr2(0)
,fDmesonCorr3(0)
,fDmesonCorr4(0)
,fDmesonCorr5(0)
,fDmesonCorr6(0)
,fDmesonCorr7(0)
,fDmesonCorr8(0)
,fDmesonCorr9(0)
,fDmesonCorr10(0)
,fDmesonCorr11(0)
,fDmesonCorr12(0)
,fLcCorr(0)
,fD0TauWeight(0)
,fDpTauWeight(0)
,fDsTauWeight(0)
,fB0TauWeight(0)
,fBpTauWeight(0)
,fBsTauWeight(0)
,fRnd(0)
{
	// default constructor
}

//______________________________________________________________________
AliAnalysisTaskBtoElecPbPbTPCTOF::~AliAnalysisTaskBtoElecPbPbTPCTOF()
{
	//Destructor 
	if(fOutputList) delete fOutputList;
	if(fExtraCuts) delete fExtraCuts;
	if(fV0Tagger) delete fV0Tagger;
}

//______________________________________________________________________
//Create Output Objects
//Here we can define the histograms and others output files
//Called once
void AliAnalysisTaskBtoElecPbPbTPCTOF::UserCreateOutputObjects()
{
  
	fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
	fAODV0Cuts = new AliAODv0KineCuts();

	fOutputList = new TList();
	fOutputList->SetOwner();	

	int nPtBins = 11;
	double ptbinning[29] = { 1., 1.1, 1.3, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8., 9., 10., 11., 12., 13., 14., 16., 18., 20. };
	double ptbinningX[12] = { 1., 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 10. };
	double ptbinningD0[14] = { 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24., 36., 50. };
	double ptbinningLc[6] = {2., 4., 6., 8., 12., 24.};
	
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
	
	int nBinsInvMassK0s = 100;
	double minInvMassK0s = 0;
	double maxInvMassK0s = 1;
	double binLimInvMassK0s[nBinsInvMassK0s+1];
	for(int i=0; i<=nBinsInvMassK0s; i++) binLimInvMassK0s[i] = minInvMassK0s + double(i)*(maxInvMassK0s/nBinsInvMassK0s);
	
	int nBinsProdRadi = 200;
	double minProdRadi = 0.;
	double maxProdRadi = 20.;
	double binLimProdRadi[nBinsProdRadi+1];
	for(int i=0; i<=nBinsProdRadi; i++) binLimProdRadi[i] = minProdRadi + double(i)*(maxProdRadi/nBinsProdRadi);
	
	int nBinsCent = 20;
	double minCent = 0.;
	double maxCent = 100.;
	double binLimCent[nBinsCent+1];
	for(int i=0; i<=nBinsCent; i++) binLimCent[i] = minCent + double(i)*(maxCent/nBinsCent);
	
	int nBinsMult = 60;
	double minMult = 0.;
	double maxMult = 30000.;
	double binLimMult[nBinsMult+1];
	for(int i=0; i<=nBinsMult; i++) binLimMult[i] = minMult + double(i)*(maxMult/nBinsMult);
	
	// event qa
	hCent_nocut = new TH1F("hCent_nocut", "centrality without cut", 100, 0, 100);
	fOutputList->Add(hCent_nocut);

	hCent_nocut2 = new TH1F("hCent_nocut2", "centrality without cut", 100, 0, 100);
	fOutputList->Add(hCent_nocut2);

	hCent_cut = new TH1F("hCent_cut", "centrality with cut", 100, 0, 100);
	fOutputList->Add(hCent_cut);

	hNrEvents = new TH1F("hNrEvents", "Number of Events", 1, 0, 1);
	fOutputList->Add(hNrEvents);

	hFilterMask = new TH1F("hFilterMask", "", 2, 0., 2.);
	fOutputList->Add(hFilterMask);

	hTPCNcls = new TH1F("hTPCNcls", "", 200, 0., 200.);
	fOutputList->Add(hTPCNcls);

	hTPCclsPID = new TH1F("hTPCclsPID", "", 200, 0., 200.);
	fOutputList->Add(hTPCclsPID);

	hTPCchi2 = new TH1F("hTPCchi2", "", 100, 0., 10.);
	fOutputList->Add(hTPCchi2);

	hTPCclsRatio = new TH1F("hTPCclsRatio", "", 15, 0., 1.5);
	fOutputList->Add(hTPCclsRatio);

	hITSNcls= new TH1F("hITSNcls", "", 10, 0., 10.);
	fOutputList->Add(hITSNcls);

	hITSclsFrac = new TH1F("hITSclsFrac", "", 11, 0., 1.1);
	fOutputList->Add(hITSclsFrac);

	hITSchi2 = new TH1F("hITSchi2", "", 30, 0., 30.);
	fOutputList->Add(hITSchi2);

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

	hITSnsigma = new TH2F("hITSnsigma", "ITS n#sigma", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hITSnsigma);

	hITSnsigmaTOFcut = new TH2F("hITSnsigmaTOFcut", "ITS n#sigma after TOF cut", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hITSnsigmaTOFcut);

	hITSnsigmaQA = new TH2F("hITSnsigmaQA", "ITS n#sigma QA", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hITSnsigmaQA);

	hTPCnsigma = new TH2F("hTPCnsigma", "TPC n#sigma", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTPCnsigma);

	hTPCnsigmaTOFcut = new TH2F("hTPCnsigmaTOFcut", "TPC n#sigma after TOF cut", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTPCnsigmaTOFcut);

	hTPCnsigmaITSTOFcut = new TH2F("hTPCnsigmaITSTOFcut", "TPC n#sigma after TOF cut", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTPCnsigmaITSTOFcut);

	hTPCnsigmaQA = new TH2F("hTPCnsigmaQA", "TPC pid cut QA", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTPCnsigmaQA);

	hTPCnsigmaPiQA = new TH2F("hTPCnsigmaPiQA", "TPC pid cut QA", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTPCnsigmaPiQA);

	hTOFnsigma = new TH2F("hTOFnsigma", "TOF n#sigma", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTOFnsigma);

	hTOFnsigmaQA = new TH2F("hTOFnsigmaQA", "TOF pid cut QA", 500, 0., 10., 400, -10., 10.);
	fOutputList->Add(hTOFnsigmaQA);

	hV0ElecTOFnsigmaDeno = new TH2F("hV0ElecTOFnsigmaDeno", "", 28, ptbinning, nBinsPID, binLimPID);
	fOutputList->Add(hV0ElecTOFnsigmaDeno);

	hV0ElecTOFnsigmaNume = new TH2F("hV0ElecTOFnsigmaNume", "", 28, ptbinning, nBinsPID, binLimPID);
	fOutputList->Add(hV0ElecTOFnsigmaNume);

	hV0ElecTPCnsigmaDeno = new TH2F("hV0ElecTPCnsigmaDeno", "", 28, ptbinning, nBinsPID, binLimPID);
	fOutputList->Add(hV0ElecTPCnsigmaDeno);

	hV0ElecTPCnsigmaNume = new TH2F("hV0ElecTPCnsigmaNume", "", 28, ptbinning, nBinsPID, binLimPID);
	fOutputList->Add(hV0ElecTPCnsigmaNume);

	hProdR_Pt = new TH2F("hProdR_Pt","",nPtBins, ptbinningX, nBinsProdRadi, binLimProdRadi);
	fOutputList->Add(hProdR_Pt);

	hV0PionMinCut = new TH3F("hV0PionMinCut","", nPtBins, ptbinningX, nBinsProdRadi, binLimProdRadi, nBinsCent, binLimCent);
	fOutputList->Add(hV0PionMinCut);

	hV0PionMaxCut = new TH3F("hV0PionMaxCut","", nPtBins, ptbinningX, nBinsProdRadi, binLimProdRadi, nBinsCent, binLimCent);
	fOutputList->Add(hV0PionMaxCut);

	hV0PionMult = new TH3F("hV0PionMult","", nPtBins, ptbinningX, nBinsProdRadi, binLimProdRadi, nBinsMult, binLimMult);
	fOutputList->Add(hV0PionMult);

	hGenBtoElecPt = new TH1F("hGenBtoElecPt","Gen B hadron pt wo any cuts", nPtBins, ptbinningX);
	fOutputList->Add(hGenBtoElecPt);
	
	hRecBtoElecPt_nocut = new TH1F("hRecBtoElecPt_nocut","Rec B hadron pt wo any cuts", nPtBins, ptbinningX);
	fOutputList->Add(hRecBtoElecPt_nocut);
	
	hRecBtoElecPt_track = new TH1F("hRecBtoElecPt_track","Rec B hadron pt wo any cuts", nPtBins, ptbinningX);
	fOutputList->Add(hRecBtoElecPt_track);
	
	hRecBtoElecPt_tof = new TH1F("hRecBtoElecPt_tof","Rec B hadron pt wo any cuts", nPtBins, ptbinningX);
	fOutputList->Add(hRecBtoElecPt_tof);
	
	hRecBtoElecPt_tpc = new TH1F("hRecBtoElecPt_tpc","Rec B hadron pt wo any cuts", nPtBins, ptbinningX);
	fOutputList->Add(hRecBtoElecPt_tpc);
	
	hD0Pt = new TH1F("hD0Pt","D0 pt wo corr (outside track loop)", 13, ptbinningD0);
	fOutputList->Add(hD0Pt);
	
	hD0PtCorr = new TH1F("hD0PtCorr","D0 pt after corr", 13, ptbinningD0);
	fOutputList->Add(hD0PtCorr);
	
	hBhadronPt = new TH1F("hBhadronPt","B hadron pt wo corr (outside track loop)", 100, 0., 100.);
	fOutputList->Add(hBhadronPt);
	
	hBhadronPtCorr = new TH1F("hBhadronPtCorr","B hadron pt after corr", 100, 0., 100.);
	fOutputList->Add(hBhadronPtCorr);

	hLcPt = new TH1F("hLcPt","#Lambda_{c} pt before corr", 5, ptbinningLc);
	fOutputList->Add(hLcPt);
	
	hLcPtCorr = new TH1F("hLcPtCorr","#Lambda_{c} pt before corr", 5, ptbinningLc);
	fOutputList->Add(hLcPtCorr);
	
	dcaTrack = new TH2F("dcaTrack", "inclusive electron's dca", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaTrack);

	dcaDmeson = new TH2F("dcaDmeson", "charm temp wo corr", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaDmeson);

	dcaDmesonCorr = new TH2F("dcaDmesonCorr", "D meson temp w corr", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaDmesonCorr);

	dcaDzero = new TH2F("dcaDzero", "D0 template", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaDzero);

	dcaDplus = new TH2F("dcaDplus", "D+ template", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaDplus);

	dcaDsplus = new TH2F("dcaDsplus", "Ds+ template", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaDsplus);

	dcaLc = new TH2F("dcaLc", "Lc+ template", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaLc);

	dcaBeauty = new TH2F("dcaBeauty", "beauty temp wo corr", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaBeauty);

	dcaBeautyCorr = new TH2F("dcaBeautyCorr", "beauty temp w corr", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaBeautyCorr);

	dcaDalitz = new TH2F("dcaDalitz", "Dalitz template", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaDalitz);

	dcaConv = new TH2F("dcaConv", "gamma template", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaConv);

	dcaPion1 = new TH2F("dcaPion1", "", 28, ptbinning, nBinsIP, binLimIP);
	fOutputList->Add(dcaPion1);

	dcaPion2 = new TH2F("dcaPion2", "", nPtBins, ptbinningX, nBinsIP, binLimIP);
	fOutputList->Add(dcaPion2);

	fD0TauWeight = new TF1("fD0TauWeight", "exp(-x/124.4)/exp(-x/122.9)", 0., 100000.);
	fDpTauWeight = new TF1("fDpTauWeight", "exp(-x/317)/exp(-x/311.8)", 0., 100000.);
	fDsTauWeight = new TF1("fDsTauWeight", "exp(-x/140)/exp(-x/149.9)", 0., 100000.);
	fB0TauWeight = new TF1("fB0TauWeight", "exp(-x/468)/exp(-x/458.7)", 0., 100000000.);
	fBpTauWeight = new TF1("fBpTauWeight", "exp(-x/462)/exp(-x/491.1)", 0., 100000000.);
	fBsTauWeight = new TF1("fBsTauWeight", "exp(-x/483)/exp(-x/439)", 0., 100000000.);

	fRnd = new TRandom3(0);

	PostData(1, fOutputList);
	
}

//Main loop
//Called for each event
void AliAnalysisTaskBtoElecPbPbTPCTOF::UserExec(Option_t *) 
{
	//Check Event
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
	if(!fAOD)
	{
		printf("ERROR: fAOD not available\n");
		return;
	}

	if(fIsMC){
		fAODMCHeader = dynamic_cast<AliAODMCHeader *>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
		if(!fAODMCHeader){
			AliError("No AliAODMCHeader");
			return;
		}
		
		fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
		if(!fAODArrayMCInfo){
			AliError("No AOD MC particles");
			return;
		}
		if(fAODArrayMCInfo->GetEntries() < 1) return;
	}
	
	//Check HFEextraCut
	if(!fExtraCuts)
	{
		fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
	}
	fExtraCuts->SetRecEventInfo(fAOD);

	//PID response
	fPidResponse = fInputHandler->GetPIDResponse();
	if(!fPidResponse)
	{
		AliDebug(1, "Using default PID Response");
		fPidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class()); 
	}
	
	//Initialize V0 electron tagger
	if(fV0Tagger){
		fV0Tagger->Reset();
		fV0Tagger->TagV0Tracks(fAOD);
	}

	//======================
	// Centrality selection
	//======================	
	fCentrality = -999.;
	fCentralityCalib = -999.;

	AliMultSelection *fMultSelection = 0x0; 
	fMultSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
	if(!fMultSelection) {
		AliWarning("AliMultSelection object not found!");
	}else{
		fCentrality = fMultSelection->GetMultiplicityPercentile("V0M", false);
		fCentralityCalib = fMultSelection->GetMultiplicityPercentile("V0M", true);
	}

	hCent_nocut->Fill(fCentralityCalib);
	hCent_nocut2->Fill(fCentrality);
	if(fCentralityCalib<fCentMin || fCentralityCalib>fCentMax) return;

	//=================
	// Event selection
	//=================
	if(!PassEventCuts(fAOD)) return;		
	
	//==============================
	// Look for kink mother for AOD
	//==============================
	double *fListOfmotherkink = 0;
	int fNumberOfVertices = 0; 
	int fNumberOfMotherkink = 0;

	fNumberOfVertices = fAOD->GetNumberOfVertices();

	fListOfmotherkink = new double[fNumberOfVertices];

	for(int ivertex=0; ivertex < fNumberOfVertices; ivertex++) 
	{
		AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
		if(!aodvertex) continue;
		if(aodvertex->GetType()==AliAODVertex::kKink) 
		{
			AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
			if(!mother) continue;
			int idmother = mother->GetID();
			fListOfmotherkink[fNumberOfMotherkink] = idmother;
			fNumberOfMotherkink++;
		}
	}

	//=================
	// MC generated pt
	//=================
	if(fIsMC){
		
		TList *lh = fAODMCHeader->GetCocktailHeaders();
		if(!lh){
			AliError("no MC header");
			return;
		}
		
		int fNTotMCpart = 0;
		for(int igen=0; igen<lh->GetEntries(); igen++){
			AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igen);
			if(!gh) continue;

			fNTotMCpart += gh->NProduced();
		}

		int hf, src, srcPdg;
		double hfpt, hfeta, srcPt, srcTau, wghtD, wghtB, wghtLc;
		hf = src = srcPdg = srcTau = -999;
		hfpt = hfeta = srcPt = wghtD = wghtB = wghtLc = -999.;
		
		for(int iMC = 0; iMC < fNTotMCpart; iMC++){
			fAODMCParticle = (AliAODMCParticle*) fAODArrayMCInfo->At(iMC);

			hf = GetHeavyFlavours(fAODMCParticle, hfpt, hfeta);
			if(TMath::Abs(hfeta)<0.5){
				if(hf==kPromptD0){
					hD0Pt->Fill(hfpt);
				}
				if(hf==kPromptLc){
					hLcPt->Fill(hfpt);
				}
			}

			if(hf==kPromptB || hf==kNonPromptD){
				if(TMath::Abs(hfeta)<0.8){
					hBhadronPt->Fill(hfpt);
				}
			}

			src = GetElecSource(fAODMCParticle, srcPt, srcPdg, srcTau);
			
			//Pseudo-rapidity cut
			if(TMath::Abs(fAODMCParticle->Eta()) > fEta) continue;

			if(TMath::Abs(hfeta)<0.5){
				if(srcPdg==421){
					wghtD = fDmesonCorr->Eval(srcPt);
					hD0PtCorr->Fill(srcPt, wghtD);
				}
				if(srcPdg==4122){
					wghtLc = fLcCorr->Eval(srcPt);
					hLcPtCorr->Fill(srcPt, wghtLc);
				}
			}

			if(src==kDirectBeauty || src==kBeautyCharm){
				hGenBtoElecPt->Fill(fAODMCParticle->Pt());
				wghtB = fBmesonCorr->Eval(srcPt);
				hBhadronPtCorr->Fill(srcPt, wghtB);
			}
		}
	}

	//--------------------
	double fSignB = -999.;
	if(fAOD->GetMagneticField()<0) fSignB = -1;
	if(fAOD->GetMagneticField()>0) fSignB = 1;
	
	hNrEvents->Fill(0);
	hCent_cut->Fill(fCentralityCalib);

	double pt, mcelectronSourcePt, fVx, fVy, fProdR, hfeImpactParam, hfeImpactParamResol, wghtB, wghtD, rndmB, rndmD;
	double fITSnSigma, fTPCnSigma, fTOFnSigma, fTPCnSigmaPi, fTOFnSigmaPi;
	int mcelectronSource, mcelectronSourcePDG, label;
	bool kinkmotherpass;
	double momTime;
	//=======================================================================
	///Track loop
	for(int iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) 
	{
		pt = mcelectronSourcePt = fVx = fVy = fProdR = hfeImpactParam = hfeImpactParamResol = wghtB = wghtD = rndmB = rndmD = momTime = -999.;
		mcelectronSource = mcelectronSourcePDG = label = -999;

		AliAODTrack *aodTrack = static_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
		
		pt = aodTrack->Pt();
		
		if(fIsMC){
			fAODMCParticle = NULL;
			
			label = TMath::Abs(aodTrack->GetLabel());
			if(label < fAODArrayMCInfo->GetEntriesFast())
				fAODMCParticle = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(label));
			if(fAODMCParticle){
				AliDebug(2, "Associated MC particle found");
				mcelectronSource = GetElecSource(fAODMCParticle, mcelectronSourcePt, mcelectronSourcePDG, momTime);
				fVx = fAODMCParticle->Xv();
				fVy = fAODMCParticle->Yv();
				fProdR = TMath::Sqrt(fVx*fVx+fVy*fVy);
			}
		}
	
		//=========    
		// RecKink
		//=========
		kinkmotherpass = true;
		for(int kinkmother = 0; kinkmother < fNumberOfMotherkink; kinkmother++) 
		{
			if(aodTrack->GetID() == fListOfmotherkink[kinkmother]) 
			{
				kinkmotherpass = kFALSE;
				continue;
			}
		}
		if(!kinkmotherpass) continue;
		
		fExtraCuts->GetHFEImpactParameters((AliVTrack *)aodTrack, hfeImpactParam, hfeImpactParamResol);

		// Reconstructed beauty pt wo cuts
		if(fIsMC)
			if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm) hRecBtoElecPt_nocut->Fill(pt);
		
		if(!PassTrackCuts(aodTrack)) continue;
		
		if(fIsMC){
			if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm){
				dcaBeauty->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
				hRecBtoElecPt_track->Fill(pt);
				if(pt>mcelectronSourcePt){
					wghtB = 1.;
				}else{
					wghtB = fBmesonCorr1->Eval(mcelectronSourcePt);
				}
				rndmB = fRnd->Rndm();
				if(TMath::Abs(mcelectronSourcePDG)==511){
					wghtB *= fB0TauWeight->Eval(momTime);
					if(rndmB<wghtB) dcaBeautyCorr->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
				}else if(TMath::Abs(mcelectronSourcePDG)==521){
					wghtB *= fBpTauWeight->Eval(momTime);
					if(rndmB<wghtB) dcaBeautyCorr->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
				}else if(TMath::Abs(mcelectronSourcePDG)==531){
					wghtB *= fBsTauWeight->Eval(momTime);
					if(rndmB<wghtB) dcaBeautyCorr->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
				}else{
					if(rndmB<wghtB) dcaBeautyCorr->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
				}
			}
			
			if(mcelectronSource==kDirectCharm){
				if(TMath::Abs(mcelectronSourcePDG)==421 || TMath::Abs(mcelectronSourcePDG)==411 || TMath::Abs(mcelectronSourcePDG)==431){
					dcaDmeson->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
					if(pt>mcelectronSourcePt){
						wghtD = 1.;
					}else{
						if(pt>=1. && pt<1.1) wghtD = fDmesonCorr1->Eval(mcelectronSourcePt);
						if(pt>=1.1 && pt<1.3) wghtD = fDmesonCorr2->Eval(mcelectronSourcePt);
						if(pt>=1.3 && pt<1.5) wghtD = fDmesonCorr3->Eval(mcelectronSourcePt);
						if(pt>=1.5 && pt<2.) wghtD = fDmesonCorr4->Eval(mcelectronSourcePt);
						if(pt>=2. && pt<2.5) wghtD = fDmesonCorr5->Eval(mcelectronSourcePt);
						if(pt>=2.5 && pt<3.) wghtD = fDmesonCorr6->Eval(mcelectronSourcePt);
						if(pt>=3. && pt<4.) wghtD = fDmesonCorr7->Eval(mcelectronSourcePt);
						if(pt>=4 && pt<5.) wghtD = fDmesonCorr8->Eval(mcelectronSourcePt);
						if(pt>=5. && pt<6.) wghtD = fDmesonCorr9->Eval(mcelectronSourcePt);
						if(pt>=6. && pt<8.) wghtD = fDmesonCorr10->Eval(mcelectronSourcePt);
						if(pt>=8. && pt<10.) wghtD = fDmesonCorr11->Eval(mcelectronSourcePt);
						if(pt>=10.) wghtD = fDmesonCorr12->Eval(mcelectronSourcePt);
					}
					rndmD = fRnd->Rndm();
					if(rndmD<wghtD) dcaDmesonCorr->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
					
					if(TMath::Abs(mcelectronSourcePDG)==421){
						wghtD *= fD0TauWeight->Eval(momTime);
						if(rndmD<wghtD) dcaDzero->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
					}
					if(TMath::Abs(mcelectronSourcePDG)==411){
						wghtD *= fDpTauWeight->Eval(momTime);
						if(rndmD<wghtD) dcaDplus->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
					}
					if(TMath::Abs(mcelectronSourcePDG)==431){
						wghtD *= fDsTauWeight->Eval(momTime);
						if(rndmD<wghtD) dcaDsplus->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
					}
				}
				if(TMath::Abs(mcelectronSourcePDG)==4122) dcaLc->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
			}
			if(mcelectronSource>=5 && mcelectronSource<=15) dcaDalitz->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
			if(mcelectronSource>=18 && mcelectronSource<=28){
				dcaConv->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
				hProdR_Pt->Fill(pt, fProdR);
			}
		}

		//=======================================================================
		// QA plots after track selection
		//=======================================================================
	
		fITSnSigma = fPidResponse->NumberOfSigmasITS(aodTrack, AliPID::kElectron);
		fTPCnSigma = fPidResponse->NumberOfSigmasTPC(aodTrack, AliPID::kElectron);
		fTOFnSigma = fPidResponse->NumberOfSigmasTOF(aodTrack, AliPID::kElectron);
		fTPCnSigmaPi = fPidResponse->NumberOfSigmasTPC(aodTrack, AliPID::kPion);
		fTOFnSigmaPi = fPidResponse->NumberOfSigmasTOF(aodTrack, AliPID::kPion);
		
		hITSnsigma->Fill(aodTrack->P(), fITSnSigma);
		hTPCnsigma->Fill(aodTrack->P(), fTPCnSigma);
		hTOFnsigma->Fill(aodTrack->P(), fTOFnSigma);

		//V0 electrons from systematic studies of TOF PID cut
		AliPID::EParticleType myv0pid = fV0Tagger->GetV0Info(aodTrack->GetID()); /// enum EParticleType: kElectron = 0, kMuon = 1, kPion = 2, etc
		if(myv0pid == AliPID::kElectron){
			// TOF eID systematics
			if(fTPCnSigma >= fTPCnsigmaLow && fTPCnSigma <= fTPCnsigmaHigh){
				hV0ElecTOFnsigmaDeno->Fill(aodTrack->Pt(), fTOFnSigma);
				if(TMath::Abs(fTOFnSigma) <= fTOFnsigma) hV0ElecTOFnsigmaNume->Fill(aodTrack->Pt(), fTOFnSigma);
			}
			// TPC eID systematics
			if(TMath::Abs(fTOFnSigma) <= fTOFnsigma){
				hV0ElecTPCnsigmaDeno->Fill(aodTrack->Pt(), fTPCnSigma);
				if(fTPCnSigma >= fTPCnsigmaLow && fTPCnSigma <= fTPCnsigmaHigh) hV0ElecTPCnsigmaNume->Fill(aodTrack->Pt(), fTPCnSigma);
			}
		}


		if(TMath::Abs(fTOFnSigma) > fTOFnsigma) continue;

		hITSnsigmaTOFcut->Fill(aodTrack->P(), fITSnSigma);
		hTPCnsigmaTOFcut->Fill(aodTrack->P(), fTPCnSigma);
		if(fIsMC)
			if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm) hRecBtoElecPt_tof->Fill(pt);

		if(fTPCnSigma > -5 && fTPCnSigma < -3){
			hTPCnsigmaPiQA->Fill(aodTrack->P(), fTPCnSigma);
			dcaPion1->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
			dcaPion2->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
		}

		//if(fITSnSigma < -4 || fITSnSigma > 2) continue;
		hTPCnsigmaITSTOFcut->Fill(aodTrack->P(), fTPCnSigma);
		
		if(fTPCnSigma < fTPCnsigmaLow || fTPCnSigma > fTPCnsigmaHigh) continue;
		
		hITSnsigmaQA->Fill(aodTrack->P(), fITSnSigma);
		hTPCnsigmaQA->Fill(aodTrack->P(), fTPCnSigma);
		hTOFnsigmaQA->Fill(aodTrack->P(), fTOFnSigma);

		if(fIsMC)
			if(mcelectronSource==kDirectBeauty || mcelectronSource==kBeautyCharm) hRecBtoElecPt_tpc->Fill(pt);
		
		if(!fIsMC){
			dcaTrack->Fill(pt, hfeImpactParam*fSignB*aodTrack->Charge());
		}

	}//End of track loop
	
	SelectV0Pions(fAOD);

	//=======================================================================
	delete fListOfmotherkink;
	PostData(1, fOutputList);
}      

//=======================================================================
void AliAnalysisTaskBtoElecPbPbTPCTOF::Terminate(Option_t *) 
{
//Draw result to the screen
//Called once at the end of the query

	fOutputList = dynamic_cast<TList*> (GetOutputData(1));
	
	if(!fOutputList) 
	{
		printf("ERROR: Output list not available\n");
		return;
	}
}

//=======================================================================

//_________________________________________
bool AliAnalysisTaskBtoElecPbPbTPCTOF::PassEventCuts(AliAODEvent *event){

	//event selection cuts
	int ntracks = event->GetNumberOfTracks();
	if(ntracks < 2) return false;

	AliAODVertex* vtTrc = event->GetPrimaryVertex();
	if(!vtTrc) return false;
	int NcontV = vtTrc->GetNContributors();

	AliAODVertex* vtSPD = event->GetPrimaryVertexSPD();
	if(!vtSPD) return false;
	int NcontSPD = vtSPD->GetNContributors();
	
	if(NcontV<2 || NcontSPD<2)return false;

	double covTrc[6],covSPD[6];
	vtTrc->GetCovarianceMatrix(covTrc);
	vtSPD->GetCovarianceMatrix(covSPD);
	double dz = vtTrc->GetZ() - vtSPD->GetZ();
	double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
	double errTrc = TMath::Sqrt(covTrc[5]);
	double nsigTot = TMath::Abs(dz)/errTot;
	double nsigTrc = TMath::Abs(dz)/errTrc;
	if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return false;// bad vertexing
	
	double zvtx = vtTrc->GetZ();
	if(TMath::Abs(zvtx) > 10) return false;
	
	return true;
}

bool AliAnalysisTaskBtoElecPbPbTPCTOF::PassTrackCuts(AliAODTrack *track){

	if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;
	if(TMath::Abs(track->Eta()) > fEta) return false;

	// basic tracking
	ULong_t status = track->GetStatus();
	if(!((status & AliVTrack::kITSrefit) && (status & AliVTrack::kTPCrefit))) return false;
	
	// dca cut
	float dcaxy = -999.; float dcaz = -999.;
	fExtraCuts->GetImpactParameters((AliVTrack *)track, dcaxy, dcaz);
	if(TMath::Abs(dcaxy) > fDCAxy || TMath::Abs(dcaz) > fDCAz) return false;
	
	// TPC cut
	unsigned short findableTPC = track->GetTPCNclsF();
	unsigned short nclustersTPC = track->GetTPCNcls();
	double FoundOverFindable = (findableTPC ? static_cast<float>(nclustersTPC)/static_cast<float>(findableTPC) : 0);
	if(track->GetTPCNcls() < fMinTPCNcls || track->GetTPCsignalN() < fMinTPCNclsPID || track->Chi2perNDF() > fMaxTPCchi2 || FoundOverFindable < fMinTPCclsRatio) return false;
	
	// ITS cut
	if(track->GetITSNcls() < fMinITSNcls) return false;
	if(fITSlayer.Contains("kFirst")){
		if(!(track->HasPointOnITSLayer(0))) return false;
		hITSlayer->Fill(1);
	}
	if(fITSlayer.Contains("kAny")){
		if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return false;
		hITSlayer->Fill(2);
	}
	if(fITSlayer.Contains("kBoth")){
		if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return false;
		hITSlayer->Fill(3);
	}
	
	//ITS quality cut
	bool HasSharedCls = kFALSE;
	double fITSNSharedcls = 0.;
	for(int itsL = 0; itsL < 6; itsL++){
		HasSharedCls = track->HasSharedPointOnITSLayer(itsL);
		if(HasSharedCls) fITSNSharedcls++;
	}
	double fITSshaClsPerNcls = fITSNSharedcls/track->GetITSNcls();
	double fITSchi2perNcls = track->GetITSchi2()/track->GetITSNcls();
	if(fITSshaClsPerNcls > fMaxITSclsFrac || fITSchi2perNcls > fMaxITSchi2) return false;

	hFilterMask->Fill(track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA));
	hEta->Fill(track->Eta());
	hPhi->Fill(track->Phi());
	hPt->Fill(track->Pt());
	hDCAxy->Fill(dcaxy);
	hDCAz->Fill(dcaz);
	hTPCNcls->Fill(track->GetTPCNcls());
	hTPCclsPID->Fill(track->GetTPCsignalN());
	hTPCchi2->Fill(track->Chi2perNDF());
	hTPCclsRatio->Fill(FoundOverFindable);
	hITSNcls->Fill(track->GetITSNcls());
	hITSclsFrac->Fill(fITSshaClsPerNcls);
	hITSchi2->Fill(fITSchi2perNcls);
	return true;
}

bool AliAnalysisTaskBtoElecPbPbTPCTOF::PassV0PionMinCuts(AliAODTrack *track){
	
	if(TMath::Abs(track->Eta())>0.8) return false;
	
	ULong_t status = track->GetStatus();
	if(!((status & AliVTrack::kITSrefit) && (status & AliVTrack::kTPCrefit))) return false;
	
	return true;
}

bool AliAnalysisTaskBtoElecPbPbTPCTOF::PassV0PionMaxCuts(AliAODTrack *track){
	
	if(TMath::Abs(track->Eta())>0.8) return false;
	
	ULong_t status = track->GetStatus();
	if(!((status & AliVTrack::kITSrefit) && (status & AliVTrack::kTPCrefit))) return false;

	if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1) && track->GetITSNcls()>=fMinITSNcls)) return false;
	
	bool HasSharedCls = false;
	double fITSNSharedcls = 0.;
	for(int itsL = 0; itsL < 6; itsL++){
		HasSharedCls = track->HasSharedPointOnITSLayer(itsL);
		if(HasSharedCls) fITSNSharedcls++;
	}
	double fITSshaClsPerNcls = fITSNSharedcls/track->GetITSNcls();
	double fITSchi2perNcls = track->GetITSchi2()/track->GetITSNcls();
	if(fITSshaClsPerNcls > fMaxITSclsFrac || fITSchi2perNcls > fMaxITSchi2) return false;

	return true;
}
//_________________________________________

int AliAnalysisTaskBtoElecPbPbTPCTOF::GetElecSource(const AliAODMCParticle * const mcpart, double &mpt, int &mpdg, double &momTime){

	if(!mcpart) return -1;
	if(!fAODArrayMCInfo) return -1;
	
	if(TMath::Abs(mcpart->GetPdgCode()) != 11 ) return kMisID;

	int origin = -1;
	Bool_t isFinalOpenCharm = kFALSE;

	int iLabel = mcpart->GetMother();
	if ((iLabel<0) || (iLabel>=fAODArrayMCInfo->GetEntriesFast())){
		AliDebug(1, "label is out of range, return\n");
		return -1;
	}
	
	AliAODMCParticle *mctrack = NULL; // will change all the time
	int tmpMomLabel=0;
	if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(iLabel))))) return -1;
	AliAODMCParticle *partMother = mctrack;	//mtrack 
	AliAODMCParticle *partMotherCopy = mctrack;	//mtrack
	int maPdgcode = mctrack->GetPdgCode();	//mpdg
	mpt = partMother->Pt();	//mpt
	mpdg = partMother->GetPdgCode();	//mpdg
	int grmaPdgcode;
	int ggrmaPdgcode;
	double gmpt, ggmpt;
	int gmpdg, ggmpdg;
	double tau;
	double eVx, eVy, eVz, cVx, cVy, cVz, bVx, bVy, bVz;
	eVx = eVy = eVz = cVx = cVy = cVz = bVx = bVy = bVz = -999.;
	eVx = mcpart->Xv();
	eVy = mcpart->Yv();
	eVz = mcpart->Zv();

	// if the mother is charmed hadron
	if( (int(TMath::Abs(maPdgcode)/100.)%10) == 4 || (int(TMath::Abs(maPdgcode)/1000.)%10) == 4 ) {
		if(TMath::Abs(maPdgcode)==411 || TMath::Abs(maPdgcode)==421 || TMath::Abs(maPdgcode)==431 || TMath::Abs(maPdgcode)==4122 || TMath::Abs(maPdgcode)==4132 || TMath::Abs(maPdgcode)==4232 || TMath::Abs(maPdgcode)==4332){
			mpt = partMother->Pt();
			mpdg = partMother->GetPdgCode();
			cVx = partMother->Xv();
			cVy = partMother->Yv();
			cVz = partMother->Zv();
			tau = TMath::Sqrt(TMath::Power(cVx-eVx,2)+TMath::Power(cVy-eVy,2)+TMath::Power(cVz-eVz,2));
			momTime = (10000*tau*partMother->M())/partMother->P(); //The factor of 10000 is to convert the units to um
			isFinalOpenCharm = kTRUE;
		}
		if (!isFinalOpenCharm) {
			return -1;
		}
		
		// iterate until find B hadron as a  mother
		for (int i=1; i<100; i++){
			int jLabel = partMother->GetMother();
			if (jLabel == -1) {
				return kDirectCharm;
			}
			if(jLabel<0 || jLabel>=fAODArrayMCInfo->GetEntriesFast()){
				AliDebug(1, "Stack label is negative, return\n");
				return -1;
			}
			
			if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(jLabel))))) {
				return -1;
			}
			int grandMaPDG = mctrack->GetPdgCode();
			if(TMath::Abs(grandMaPDG)==511 || TMath::Abs(grandMaPDG)==521 || TMath::Abs(grandMaPDG)==531 || TMath::Abs(grandMaPDG)==5122 || TMath::Abs(grandMaPDG)==5132 || TMath::Abs(grandMaPDG)==5232 || TMath::Abs(grandMaPDG)==5332){
				mpt = mctrack->Pt();
				mpdg = mctrack->GetPdgCode();
				bVx = mctrack->Xv();
				bVy = mctrack->Yv();
				bVz = mctrack->Zv();
				tau = TMath::Sqrt(TMath::Power(bVx-cVx,2)+TMath::Power(bVy-cVy,2)+TMath::Power(bVz-cVz,2));
				momTime = (10000*tau*mctrack->M())/mctrack->P(); //The factor of 10000 is to convert the units to um
				return kBeautyCharm;
			}
			partMother = mctrack;
		} // end of iteration 
	}
	
	// if the mother is beauty hadron
	else if( (int(TMath::Abs(maPdgcode)/100.)%10) == 5 || (int(TMath::Abs(maPdgcode)/1000.)%10) == 5 ) {
		if (TMath::Abs(maPdgcode)==511 || TMath::Abs(maPdgcode)==521 || TMath::Abs(maPdgcode)==531 || TMath::Abs(maPdgcode)==5122 || TMath::Abs(maPdgcode)==5132 || TMath::Abs(maPdgcode)==5232 || TMath::Abs(maPdgcode)==5332){
			mpt = partMotherCopy->Pt();
			mpdg = partMotherCopy->GetPdgCode();
			bVx = partMotherCopy->Xv();
			bVy = partMotherCopy->Yv();
			bVz = partMotherCopy->Zv();
			tau = TMath::Sqrt(TMath::Power(bVx-eVx,2)+TMath::Power(bVy-eVy,2)+TMath::Power(bVz-eVz,2));
			momTime = (10000*tau*partMotherCopy->M())/partMotherCopy->P(); //The factor of 10000 is to convert the units to um
			return kDirectBeauty;
		}
	}
	
	// if the mother is gamma
	else if ( TMath::Abs(maPdgcode) == 22 ) {
		
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
		partMother = mctrack; // grand mother
		partMotherCopy = mctrack; // grand mother
		mpt = partMother->Pt(); // grand mother pT
		mpdg = partMother->GetPdgCode();
		maPdgcode = partMother->GetPdgCode(); // grand mother PDG
		
		// check if the ligth meson is the decay product of heavy mesons
		tmpMomLabel = partMother->GetMother();
		if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())) {//grandgrandmother
			if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
				partMother = mctrack; //grand grand mother
        grmaPdgcode = partMother->GetPdgCode(); //grand grand mother PDG
				mpt = partMother->Pt(); // grand grand mother PDG
				mpdg = partMother->GetPdgCode();
				gmpt = partMother->Pt(); // grand grand mother PDG
				gmpdg = partMother->GetPdgCode();

				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == 5 || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == 5 ) {
					return kGammaB2M;
				}
				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == 4 || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == 4 ) {
					return kGammaD2M;
				}
				
				tmpMomLabel = partMother->GetMother();
				if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())) {//grandgrandgrandmother
					if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
						partMother = mctrack; // grand grand grand mother
						ggrmaPdgcode = partMother->GetPdgCode(); // grand grand grand mother PDG
						mpt = partMother->Pt(); // grand grand grand mother PDG
						mpdg = partMother->GetPdgCode();
						ggmpt = partMother->Pt(); // grand grand grand mother PDG
						ggmpdg = partMother->GetPdgCode();

						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == 5 || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == 5 ) {
							return kGammaB2M;
						}
						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == 4 || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == 4 ) {
							return kGammaD2M;
						}
					}
				}
				
				if ( TMath::Abs(maPdgcode) == 111 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					else if(grmaPdgcode == 310) return kGammaK0s2P;
					else if(grmaPdgcode == 130) return kGammaK0l2P;
					else if(TMath::Abs(grmaPdgcode) == 321) return kGammaK2P;
					else if(TMath::Abs(grmaPdgcode) == 3122) return kGammaLamda2P;
					else if(grmaPdgcode == 3222) return kGammaSigma2P;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaPi0;
				}
				else if ( TMath::Abs(maPdgcode) == 221 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaEta;
				}
				else if ( TMath::Abs(maPdgcode) == 223 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaOmega;
				}
				else if ( TMath::Abs(maPdgcode) == 333 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaPhi;
				}
				else if ( TMath::Abs(maPdgcode) == 331 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kGammaM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaEtaPrime;
				}
				else if ( TMath::Abs(maPdgcode) == 113 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kGammaM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kGammaRho0;
				}
				else origin = kElse;//grandgrandmother but nothing we identify
			}//mctrack grandgrandmother
		}
		else {
			// grandmother is primary
			if ( TMath::Abs(maPdgcode) == 111 ) {
				return kGammaPi0;
			}
			else if ( TMath::Abs(maPdgcode) == 221 ) {
				return kGammaEta;
			}
			else if ( TMath::Abs(maPdgcode) == 223 ) {
				return kGammaOmega;
			}
			else if ( TMath::Abs(maPdgcode) == 333 ) {
				return kGammaPhi;
			}
			else if ( TMath::Abs(maPdgcode) == 331 ) {
				return kGammaEtaPrime;
			}
			else if ( TMath::Abs(maPdgcode) == 113 ) {
				return kGammaRho0;
			}
			else origin = kElse;//grandmother is primary but nothing we identify
		}
		return origin;
	}

	// if the mother is light meson
	else {
		
		tmpMomLabel = partMotherCopy->GetMother();
		mpt = partMotherCopy->Pt(); // mother pT
		mpdg = partMotherCopy->GetPdgCode();
		if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())) {//grandmother
			if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
				partMother = mctrack; // grand mother
				grmaPdgcode = partMother->GetPdgCode(); // grand mother PDG
				mpt = partMother->Pt(); // grand mother pT
				mpdg = partMother->GetPdgCode();
				gmpt = partMother->Pt(); // grand mother pT
				gmpdg = partMother->GetPdgCode();

				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == 5 || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == 5 ) {
					return kB2M;
				}
				if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == 4 || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == 4 ) {
					return kD2M;
				}
				
				tmpMomLabel = partMother->GetMother();
				if((tmpMomLabel>=0) && (tmpMomLabel<fAODArrayMCInfo->GetEntriesFast())) {//grandgrandmother
					if((mctrack = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(tmpMomLabel))))) {
						partMother = mctrack; // grand grand mother
						ggrmaPdgcode = partMother->GetPdgCode(); // grand grand mother PDG
						mpt = partMother->Pt(); // grand grand mother pT
						mpdg = partMother->GetPdgCode();
						ggmpt = partMother->Pt(); // grand grand mother pT
						ggmpdg = partMother->GetPdgCode();

						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == 5 || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == 5 ) {
							return kB2M;
						}
						if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == 4 || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == 4 ) {
							return kD2M;
						}
					}
				}
				
				if ( TMath::Abs(maPdgcode) == 111 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					else if(grmaPdgcode == 310) return kK0s2P;
					else if(grmaPdgcode == 130) return kK0l2P;
					else if(TMath::Abs(grmaPdgcode) == 321) return kK2P;
					else if(TMath::Abs(grmaPdgcode) == 3122) return kLamda2P;
					else if(grmaPdgcode == 3222) return kSigma2P;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kPi0;
				}
				else if ( TMath::Abs(maPdgcode) == 221 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kEta;
				}
				else if ( TMath::Abs(maPdgcode) == 223 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kOmega;
				}
				else if ( TMath::Abs(maPdgcode) == 333 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kPhi;
				}
				else if ( TMath::Abs(maPdgcode) == 331 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kEtaPrime;
				}
				else if ( TMath::Abs(maPdgcode) == 113 ) {
					mpt = gmpt;
					mpdg = gmpdg;
					if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kM2M;
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kRho0;
				}
				else if ( TMath::Abs(maPdgcode) == 321 ) {
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kKe3;
				}
				else if ( TMath::Abs(maPdgcode) == 130 ) {
					mpt = partMotherCopy->Pt();
					mpdg = partMotherCopy->GetPdgCode();
					return kK0L;
				}
				else origin = kElse;//grandmother but nothing we identidy
			}//mctrack grandmother
		}
		else {
			// no grandmother
			if ( TMath::Abs(maPdgcode) == 111 ) {
				return kPi0;
			}
			else if ( TMath::Abs(maPdgcode) == 221 ) {
				return kEta;
			}
			else if ( TMath::Abs(maPdgcode) == 223 ) {
				return kOmega;
			}
			else if ( TMath::Abs(maPdgcode) == 333 ) {
				return kPhi;
			}
			else if ( TMath::Abs(maPdgcode) == 331 ) {
				return kEtaPrime;
			}
			else if ( TMath::Abs(maPdgcode) == 113 ) {
				return kRho0;
			}
			else if ( TMath::Abs(maPdgcode) == 321 ) {
				return kKe3;
			}
			else if ( TMath::Abs(maPdgcode) == 130 ) {
				return kK0L;
			}
			else origin = kElse;//mother but nothing we identify
		}
	}//mother is something different from J/psi,charm,beauty or gamma
	
	return origin;

}

//_________________________________________

int AliAnalysisTaskBtoElecPbPbTPCTOF::GetHeavyFlavours(const AliAODMCParticle * const mcpart, double &mpt, double &meta){

	if(!mcpart) return -1;
	if(!fAODArrayMCInfo) return -1;
	
	int pdgHF = TMath::Abs(mcpart->GetPdgCode());
	mpt = mcpart->Pt();
	meta = mcpart->Eta();
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
				mpt = mctrack->Pt();
				meta = mctrack->Eta();
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

int AliAnalysisTaskBtoElecPbPbTPCTOF::GetGammaPt(const AliAODMCParticle * const mcpart, double &mpt){

	if(!mcpart) return -1;
	if(!fAODArrayMCInfo) return -1;
	
	if(TMath::Abs(mcpart->GetPdgCode())!=11) return -1;

	int moLabel = mcpart->GetMother();
	if ((moLabel<0) || (moLabel>=fAODArrayMCInfo->GetEntriesFast())){
		AliDebug(1, "label is out of range, return\n");
		return -1;
	}
	
	AliAODMCParticle *mother = NULL;
	if(!(mother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(moLabel))))) return -1;
	int moPdg = mother->GetPdgCode();	//mpdg
	
	// if the mother is gamma
	if(TMath::Abs(moPdg)==22){
		mpt = mother->Pt();	//mpt
		
		int gmoLabel=0;
		gmoLabel = mother->GetMother();  // mother of photon
		if(gmoLabel==-1) return kDirectGamma;  // no grandmother
		if((gmoLabel<0) || (gmoLabel>=fAODArrayMCInfo->GetEntriesFast())){
			return -1;
		}
		AliAODMCParticle *gmother = NULL;
		if(!(gmother=dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(gmoLabel))))){
			return -1;
		}
		int gmoPdg = gmother->GetPdgCode(); // grand mother PDG
		if(gmoPdg/100==4||gmoPdg/1000==4) return -1;
		if(gmoPdg/100==5||gmoPdg/1000==5) return -1;
		int ggmoLabel=0;
		ggmoLabel = gmother->GetMother();
		if((ggmoLabel>=0) && (ggmoLabel<fAODArrayMCInfo->GetEntriesFast())){//grandgrandmother
			AliAODMCParticle *ggmother = NULL;
			if((ggmother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(ggmoLabel))))){
        int ggmoPdg = ggmother->GetPdgCode(); //grand grand mother PDG
				if(ggmoPdg/100==4||ggmoPdg/1000==4) return -1;
				if(ggmoPdg/100==5||ggmoPdg/1000==5) return -1;

				if(TMath::Abs(gmoPdg)==111){
					if(ggmoPdg==221 || ggmoPdg==223 || ggmoPdg==333 || ggmoPdg==331 || ggmoPdg==113) return -1;
					else if(ggmoPdg == 310 || ggmoPdg == 130 || TMath::Abs(ggmoPdg) == 321 || TMath::Abs(ggmoPdg) == 3122 || ggmoPdg == 3222) return kBaryonGamma;
					return kDalitzGamma;
				}
				else if(TMath::Abs(gmoPdg)==221){
					if(ggmoPdg==111 || ggmoPdg==223 || ggmoPdg==333 || ggmoPdg==331 || ggmoPdg==113) return -1;
					return kDalitzGamma;
				}
			}
			int gggmoLabel=0;
			gggmoLabel = ggmother->GetMother();
			if((gggmoLabel>=0) && (gggmoLabel<fAODArrayMCInfo->GetEntriesFast())){
				AliAODMCParticle *gggmother = NULL;
				if((gggmother = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(TMath::Abs(gggmoLabel))))){
					int gggmoPdg = gggmother->GetPdgCode();
					if(gggmoPdg/100==4||gggmoPdg/1000==4) return -1;
					if(gggmoPdg/100==5||gggmoPdg/1000==5) return -1;
				}
			}
		}else{
			if(TMath::Abs(gmoPdg)==111) return kDalitzGamma;
			else if(TMath::Abs(gmoPdg)==221) return kDalitzGamma;
		}
		return -1;
	}
	return -1;
}

void AliAnalysisTaskBtoElecPbPbTPCTOF::SelectV0Pions(AliAODEvent *evt){
	
	fAODV0Cuts->SetEvent(evt);
	AliAODv0 *v0 = NULL;

	int V0MotherPdg, V0Daughter1Pdg, V0Daughter2Pdg;
	double recoRadius;
	AliAODTrack *pTrack, *nTrack;
	for(int iV0s = 0; iV0s<evt->GetNumberOfV0s(); iV0s++){
		v0 = evt->GetV0(iV0s);
		if(!v0) continue;

		if(fAODV0Cuts->ProcessV0(v0, V0MotherPdg, V0Daughter1Pdg, V0Daughter2Pdg)){
		
			recoRadius = v0->RadiusSecVtx();
			if(TMath::Abs(V0MotherPdg)==310 && recoRadius>0.5){

				pTrack=(AliAODTrack *)v0->GetDaughter(0); //0->Positive Daughter
				nTrack=(AliAODTrack *)v0->GetDaughter(1); //1->Negative Daughter
				if (!pTrack || !nTrack) {
					Printf("ERROR: Could not retreive one of the daughter tracks");
					continue;
				}
				if(pTrack->GetSign()==nTrack->GetSign()) continue;
				
				// minimal cut
				if(PassV0PionMinCuts(pTrack)) hV0PionMinCut->Fill(pTrack->Pt(), recoRadius, fCentrality);
				if(PassV0PionMinCuts(nTrack)) hV0PionMinCut->Fill(nTrack->Pt(), recoRadius, fCentrality);

				// maximal cut
				if(PassV0PionMaxCuts(pTrack)) hV0PionMaxCut->Fill(pTrack->Pt(), recoRadius, fCentrality);
				if(PassV0PionMaxCuts(nTrack)) hV0PionMaxCut->Fill(nTrack->Pt(), recoRadius, fCentrality);
			
				if(PassV0PionMaxCuts(pTrack)) hV0PionMult->Fill(pTrack->Pt(), recoRadius, evt->GetNumberOfV0s());
				if(PassV0PionMaxCuts(nTrack)) hV0PionMult->Fill(nTrack->Pt(), recoRadius, evt->GetNumberOfV0s());
			}
		}
	}
}







