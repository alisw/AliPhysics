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

/* $Id: AliAnalysisTaskCMEV0PID.cxx  Rihan Haque 14/02/2018 $ */

//-- general include---
#include "TChain.h"
#include "TTree.h"
#include "TGrid.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TMatrixDSym.h"

#include "TMath.h"
#include "stdio.h"
#include "Riostream.h"

//---- manager and handler---
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

//---V0 and ZDC info---
#include "AliAODZDC.h"
#include "AliAODVZERO.h"
#include "AliAODVertex.h"

//---AOD,ESD event--
#include "AliESDEvent.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"

//----for PID-----
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"

//----- Vevent and tracks
#include "AliVEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliCentrality.h"

//----- must include-------
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliPhysicsSelection.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskCMEV0PID.h"

//using namespace std;

using std::cout;
using std::endl;
using std::vector;


ClassImp(AliAnalysisTaskCMEV0PID)


AliAnalysisTaskCMEV0PID::AliAnalysisTaskCMEV0PID(const char *name): AliAnalysisTaskSE(name),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  mfileFBHijing(NULL),
  fListFBHijing(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),
  fHistTaskConfigParameters(NULL),
  fHistPileUpCount(NULL),
  fHistMultSelPUCount(NULL),
  fHistEtaPtBefore(NULL),
  fHistEtaPtAfter(NULL),
  fHistTPCvsGlobalMultBefore(NULL),
  fHistTPCvsGlobalMultAfter(NULL),
  fHistTPCdEdxvsPBefore(NULL),
  fHistTPCdEdxvsPAfter(NULL),
  fHistTOFBetavsPBefore(NULL),
  fHistTOFBetavsPAfter(NULL),
  fHistTOFMassvsPtBefore(NULL),
  fHistTOFMatchCount(NULL),
  fHistTPCVsESDTrkBefore(NULL),
  fHistTPCVsESDTrkAfter(NULL),
  fHistTPConlyVsCL1Before(NULL),
  fHistTPConlyVsV0MBefore(NULL),
  fHistTPConlyVsCL1After(NULL),
  fHistTPConlyVsV0MAfter(NULL),
  fHistGlobalVsV0MBefore(NULL),
  fHistGlobalVsV0MAfter(NULL),
  fHistRawVsCorrMultFB(NULL),
  hCentvsTPCmultCuts(NULL),
  fHV0AEventPlaneVsCent(NULL),
  fHV0CEventPlaneVsCent(NULL),
  fHTPCAEventPlaneVsCent(NULL),
  fHTPCCEventPlaneVsCent(NULL),
  fHTPCEventPlaneVsCent(NULL),
  fV0MultChVsRun(NULL),
  fCentDistBefore(NULL),
  fCentDistAfter(NULL),
  fHCorrectV0M(NULL),
  fHAvgerageQnV0A(NULL),  
  fHAvgerageQnV0C(NULL),
  fHCentWeightForRun(NULL),
  fQAEtaPhiAfterNUA(NULL),
  fQAEtaPhiAfterNUAPion(NULL),
  fQAEtaPhiAfterNUAKaon(NULL),
  fQAEtaPhiAfterNUAProton(NULL),
  fV0AQ2xVsCentRun(NULL),
  fV0AQ2yVsCentRun(NULL),
  fV0CQ2xVsCentRun(NULL),
  fV0CQ2yVsCentRun(NULL),
  fV0AQ3xVsCentRun(NULL),
  fV0AQ3yVsCentRun(NULL),
  fV0CQ3xVsCentRun(NULL),
  fV0CQ3yVsCentRun(NULL),
  fTPCAQ2xVsCentRun(NULL),
  fTPCAQ2yVsCentRun(NULL),
  fTPCCQ2xVsCentRun(NULL),
  fTPCCQ2yVsCentRun(NULL),
  fTPCAQ3xVsCentRun(NULL),
  fTPCAQ3yVsCentRun(NULL),
  fTPCCQ3xVsCentRun(NULL),
  fTPCCQ3yVsCentRun(NULL),
  fTPCAQ4xVsCentRun(NULL),
  fTPCAQ4yVsCentRun(NULL),
  fTPCCQ4xVsCentRun(NULL),
  fTPCCQ4yVsCentRun(NULL),
  fTPCFQ2xVsCentRun(NULL),
  fTPCFQ2yVsCentRun(NULL),
  fZPASignalPerChVsCent(NULL),
  fZPCSignalPerChVsCent(NULL),
  fZNASignalPerChVsCent(NULL),
  fZNCSignalPerChVsCent(NULL),
  fCentDistVsVz(NULL),
  fFilterBit(1),
  gN(1),
  gM(1),
  gPsiN(2),
  fOldRunNum(111),
  fEventCount(0),
  fTPCclustMin(70),
  fNSigmaCut(2),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fTrkChi2Min(0.1),
  fDCAxyMax(2.4),
  fDCAzMax(3.2),
  fdEdxMin(10.),
  fCentralityPercentMin(0),
  fCentralityPercentMax(90),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  fMinVzCut(-10.0),
  fMaxVzCut(10.0),
  bApplyMCcorr(kFALSE),
  bV0MGainCorr(kFALSE),
  bSkipPileUpCut(kFALSE), 
  bFillNUAHistPID(kFALSE),
  bUseKinkTracks(kFALSE),
  sPathOfMCFile("/alien"),
  sNucleiTP("PbPb"),
  sCentrEstimator("V0M"),
  fpX1X2PosT1(NULL),
  fpX1X3PosT1(NULL),
  fpX1X2PosT2(NULL),
  fpX1X3PosT2(NULL),
  fpX1X2PosT3(NULL),
  fpX1X3PosT3(NULL),
  fpX1X2PosT4(NULL),
  fpX1X3PosT4(NULL),
  fpX1X2NegT1(NULL),
  fpX1X3NegT1(NULL),
  fpX1X2NegT2(NULL),
  fpX1X3NegT2(NULL),
  fpX1X2NegT3(NULL),
  fpX1X3NegT3(NULL),
  fpX1X2NegT4(NULL),
  fpX1X3NegT4(NULL),
  fpX1X2OppT1(NULL),
  fpX1X2OppT2(NULL),
  fpX1X2OppT3(NULL),
  fpX1X2OppT4(NULL),
  fpX1X3PosT1EP2(NULL),
  fpX1X3PosT2EP2(NULL),
  fpX1X3PosT3EP2(NULL),
  fpX1X3PosT4EP2(NULL),
  fpX1X3NegT1EP2(NULL),
  fpX1X3NegT2EP2(NULL),
  fpX1X3NegT3EP2(NULL),
  fpX1X3NegT4EP2(NULL),
  fpX1PosCosT1(NULL),
  fpX1NegCosT1(NULL),
  fpX1PosSinT1(NULL),
  fpX1NegSinT1(NULL),
  fpX1EventT1EP1(NULL),
  fpX1EventT1EP2(NULL),
  fpY1EventT1EP1(NULL),
  fpY1EventT1EP2(NULL),
  fHistEventCount(NULL)
{
  for(int i=0;i<3;i++){
    fHistPtwithTPCNsigma[i]=NULL;
    fHistPtwithTOFmasscut[i]=NULL;
    fHistPtwithTOFSignal[i]=NULL;
    fHistTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCnSigmavsPtAfter[i]=NULL;
    fHistTPCTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCdEdxvsPtPIDAfter[i]=NULL;
  }
  for(int i=0;i<5;i++){
    fHCorrectNUApos[i] = NULL;
    fHCorrectNUAneg[i] = NULL;
  }
  for(int i=0;i<5;i++){  // for PID 
    fHCorrectNUAposPion[i] = NULL;
    fHCorrectNUAnegPion[i] = NULL;
    fHCorrectNUAposKaon[i] = NULL;
    fHCorrectNUAnegKaon[i] = NULL;
    fHCorrectNUAposProton[i] = NULL;
    fHCorrectNUAnegProton[i] = NULL;
  }
  //3p vs centrality
  for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
      fHist_Corr3p_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_NN[i][j]  =  NULL;
    }
    for(int j=0;j<4;j++) {
      fHist_Reso2n_EP_Norm_Det[i][j] =  NULL;
    }
  }

  //Differential 3p Charge:
  for(int i=0;i<2;i++){
    for(int j=0;j<6;j++){
      fHist_Corr3p_pTSum_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_pTSum_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_pTSum_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_pTSum_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_pTSum_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_pTSum_EP_V0C_NN[i][j] = NULL;
    
      fHist_Corr3p_pTDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_pTDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_pTDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_pTDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_pTDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_pTDiff_EP_V0C_NN[i][j] = NULL;

      fHist_Corr3p_EtaDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_EtaDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_EtaDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_EtaDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_EtaDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_EtaDiff_EP_V0C_NN[i][j] = NULL;
    }
  }  
 
  for(int i=0;i<4;i++){
    for(int j=0;j<5;j++){
      fHist3DEtaPhiVz_Pos_Run[i][j]=NULL;
      fHist3DEtaPhiVz_Neg_Run[i][j]=NULL;
    }
  }
  for(int i=0;i<10;i++){
    fFB_Efficiency_Cent[i] = NULL;
    //fFB_Efficiency_Pion_Cent[i] = NULL;
    //fFB_Efficiency_Kaon_Cent[i] = NULL;
    //fFB_Efficiency_Proton_Pos_Cent[i] = NULL;
    //fFB_Efficiency_Proton_Neg_Cent[i] = NULL;
  }


  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}
//______________________________empty constructor_______________________
AliAnalysisTaskCMEV0PID::AliAnalysisTaskCMEV0PID():
  AliAnalysisTaskSE(),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  mfileFBHijing(NULL),
  fListFBHijing(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),
  fHistTaskConfigParameters(NULL),
  fHistPileUpCount(NULL),
  fHistMultSelPUCount(NULL),
  fHistEtaPtBefore(NULL),
  fHistEtaPtAfter(NULL),
  fHistTPCvsGlobalMultBefore(NULL),
  fHistTPCvsGlobalMultAfter(NULL),
  fHistTPCdEdxvsPBefore(NULL),
  fHistTPCdEdxvsPAfter(NULL),
  fHistTOFBetavsPBefore(NULL),
  fHistTOFBetavsPAfter(NULL),
  fHistTOFMassvsPtBefore(NULL),
  fHistTOFMatchCount(NULL),
  fHistTPCVsESDTrkBefore(NULL),
  fHistTPCVsESDTrkAfter(NULL),
  fHistTPConlyVsCL1Before(NULL),
  fHistTPConlyVsV0MBefore(NULL),
  fHistTPConlyVsCL1After(NULL),
  fHistTPConlyVsV0MAfter(NULL),
  fHistGlobalVsV0MBefore(NULL),
  fHistGlobalVsV0MAfter(NULL),
  fHistRawVsCorrMultFB(NULL),
  hCentvsTPCmultCuts(NULL),
  fHV0AEventPlaneVsCent(NULL),
  fHV0CEventPlaneVsCent(NULL),
  fHTPCAEventPlaneVsCent(NULL),
  fHTPCCEventPlaneVsCent(NULL),
  fHTPCEventPlaneVsCent(NULL),
  fV0MultChVsRun(NULL),
  fCentDistBefore(NULL),
  fCentDistAfter(NULL),
  fHCorrectV0M(NULL),
  fHAvgerageQnV0A(NULL),  
  fHAvgerageQnV0C(NULL),
  fHCentWeightForRun(NULL),
  fQAEtaPhiAfterNUA(NULL),
  fQAEtaPhiAfterNUAPion(NULL),
  fQAEtaPhiAfterNUAKaon(NULL),
  fQAEtaPhiAfterNUAProton(NULL),
  fV0AQ2xVsCentRun(NULL),
  fV0AQ2yVsCentRun(NULL),
  fV0CQ2xVsCentRun(NULL),
  fV0CQ2yVsCentRun(NULL),
  fV0AQ3xVsCentRun(NULL),
  fV0AQ3yVsCentRun(NULL),
  fV0CQ3xVsCentRun(NULL),
  fV0CQ3yVsCentRun(NULL),  
  fTPCAQ2xVsCentRun(NULL),
  fTPCAQ2yVsCentRun(NULL),
  fTPCCQ2xVsCentRun(NULL),
  fTPCCQ2yVsCentRun(NULL),
  fTPCAQ3xVsCentRun(NULL),
  fTPCAQ3yVsCentRun(NULL),
  fTPCCQ3xVsCentRun(NULL),
  fTPCCQ3yVsCentRun(NULL),
  fTPCAQ4xVsCentRun(NULL),
  fTPCAQ4yVsCentRun(NULL),
  fTPCCQ4xVsCentRun(NULL),
  fTPCCQ4yVsCentRun(NULL),
  fTPCFQ2xVsCentRun(NULL),
  fTPCFQ2yVsCentRun(NULL),
  fZPASignalPerChVsCent(NULL),
  fZPCSignalPerChVsCent(NULL),
  fZNASignalPerChVsCent(NULL),
  fZNCSignalPerChVsCent(NULL),
  fCentDistVsVz(NULL),
  fFilterBit(1),
  gN(1),
  gM(1),
  gPsiN(2),
  fOldRunNum(111),
  fEventCount(0),
  fTPCclustMin(70),
  fNSigmaCut(2),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fTrkChi2Min(0.1),
  fDCAxyMax(2.4),
  fDCAzMax(3.2),
  fdEdxMin(10.),
  fCentralityPercentMin(0),
  fCentralityPercentMax(90),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  fMinVzCut(-10.0),
  fMaxVzCut(10.0),
  bApplyMCcorr(kFALSE),
  bV0MGainCorr(kFALSE),
  bSkipPileUpCut(kFALSE), 
  bFillNUAHistPID(kFALSE),
  bUseKinkTracks(kFALSE),
  sPathOfMCFile("/alien"),
  sNucleiTP("PbPb"),
  sCentrEstimator("V0M"),
  fpX1X2PosT1(NULL),
  fpX1X3PosT1(NULL),
  fpX1X2PosT2(NULL),
  fpX1X3PosT2(NULL),
  fpX1X2PosT3(NULL),
  fpX1X3PosT3(NULL),
  fpX1X2PosT4(NULL),
  fpX1X3PosT4(NULL),
  fpX1X2NegT1(NULL),
  fpX1X3NegT1(NULL),
  fpX1X2NegT2(NULL),
  fpX1X3NegT2(NULL),
  fpX1X2NegT3(NULL),
  fpX1X3NegT3(NULL),
  fpX1X2NegT4(NULL),
  fpX1X3NegT4(NULL),
  fpX1X2OppT1(NULL),
  fpX1X2OppT2(NULL),
  fpX1X2OppT3(NULL),
  fpX1X2OppT4(NULL),
  fpX1X3PosT1EP2(NULL),
  fpX1X3PosT2EP2(NULL),
  fpX1X3PosT3EP2(NULL),
  fpX1X3PosT4EP2(NULL),
  fpX1X3NegT1EP2(NULL),
  fpX1X3NegT2EP2(NULL),
  fpX1X3NegT3EP2(NULL),
  fpX1X3NegT4EP2(NULL),
  fpX1PosCosT1(NULL),
  fpX1NegCosT1(NULL),
  fpX1PosSinT1(NULL),
  fpX1NegSinT1(NULL),
  fpX1EventT1EP1(NULL),
  fpX1EventT1EP2(NULL),
  fpY1EventT1EP1(NULL),
  fpY1EventT1EP2(NULL),
  fHistEventCount(NULL)
{
  for(int i=0;i<3;i++){
    fHistPtwithTPCNsigma[i]=NULL;
    fHistPtwithTOFmasscut[i]=NULL;
    fHistPtwithTOFSignal[i]=NULL;
    fHistTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCnSigmavsPtAfter[i]=NULL;
    fHistTPCTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCdEdxvsPtPIDAfter[i]=NULL;
  }
  for(int i=0;i<5;i++){
    fHCorrectNUApos[i] = NULL;
    fHCorrectNUAneg[i] = NULL;
  }
  for(int i=0;i<5;i++){  // for PID NUA
    fHCorrectNUAposPion[i] = NULL;
    fHCorrectNUAnegPion[i] = NULL;
    fHCorrectNUAposKaon[i] = NULL;
    fHCorrectNUAnegKaon[i] = NULL;
    fHCorrectNUAposProton[i] = NULL;
    fHCorrectNUAnegProton[i] = NULL;
  }
  //3p vs Centrality
  for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
      fHist_Corr3p_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_NN[i][j]  =  NULL;
    }
    for(int j=0;j<4;j++) {
      fHist_Reso2n_EP_Norm_Det[i][j] =  NULL;
    }
  }

 //2p vs Centrality:

 //2p vs Refm:

  //Differential Charge:
  for(int i=0;i<2;i++){
    for(int j=0;j<6;j++){
      fHist_Corr3p_pTSum_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_pTSum_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_pTSum_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_pTSum_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_pTSum_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_pTSum_EP_V0C_NN[i][j] = NULL;
    
      fHist_Corr3p_pTDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_pTDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_pTDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_pTDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_pTDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_pTDiff_EP_V0C_NN[i][j] = NULL;

      fHist_Corr3p_EtaDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_EtaDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_EtaDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_EtaDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_EtaDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_EtaDiff_EP_V0C_NN[i][j] = NULL;
    }
  }  
  for(int i=0;i<4;i++){
    for(int j=0;j<5;j++){
      fHist3DEtaPhiVz_Pos_Run[i][j]=NULL;
      fHist3DEtaPhiVz_Neg_Run[i][j]=NULL;
    }
  }
  for(int i=0;i<10;i++){
    fFB_Efficiency_Cent[i] = NULL;
    //fFB_Efficiency_Pion_Cent[i] = NULL;
    //fFB_Efficiency_Kaon_Cent[i] = NULL;
    //fFB_Efficiency_Proton_Pos_Cent[i] = NULL;
    //fFB_Efficiency_Proton_Neg_Cent[i] = NULL;
  }
}

//___________________________ destructor ___________________________
AliAnalysisTaskCMEV0PID::~AliAnalysisTaskCMEV0PID()
{
  //Destructor
  //if(fPIDResponse)   delete fPIDResponse;
  //if(fMultSelection) delete fMultSelection;

  /* 
  if(fHCorrectV0M)    delete fHCorrectV0M;
  if(fHAvgerageQnV0A) delete fHAvgerageQnV0A;
  if(fHAvgerageQnV0C) delete fHAvgerageQnV0C;

  if(mfileFBHijing->IsOpen()){
     mfileFBHijing->Close();
     if(fListFBHijing) delete fListFBHijing;
  }
  for(int i=0;i<10;i++){
    if(fFB_Efficiency_Cent[i])
      delete fFB_Efficiency_Cent[i];
  }
  for(int i=0;i<5;i++){
    if(fHCorrectNUApos[i]) delete fHCorrectNUApos[i];
    if(fHCorrectNUAneg[i]) delete fHCorrectNUAneg[i];
  }
  for(int i=0;i<5;i++){  // for PID 
    if(fHCorrectNUAposPion[i])   delete fHCorrectNUAposPion[i];
    if(fHCorrectNUAnegPion[i])   delete fHCorrectNUAnegPion[i];
    if(fHCorrectNUAposKaon[i])   delete fHCorrectNUAposKaon[i];
    if(fHCorrectNUAnegKaon[i])   delete fHCorrectNUAnegKaon[i];
    if(fHCorrectNUAposProton[i]) delete fHCorrectNUAposProton[i];
    if(fHCorrectNUAnegProton[i]) delete fHCorrectNUAnegProton[i];
  }
*/

  //Delete the clones
  if(fListFBHijing) delete fListFBHijing;
  if(fListNUACorr)  delete fListNUACorr;
  if(fListV0MCorr)  delete fListV0MCorr;
      
  if(fListHist)      delete fListHist;  
  if(fAnalysisUtil)  delete fAnalysisUtil; // its 'new' !!

}//---------------- sanity ------------------------



void AliAnalysisTaskCMEV0PID::UserCreateOutputObjects()
{

  //std::cout<<"\n..UserCreateOutputObject called.. with isCorr = "<<isCorr<<"\n ....check if succeeded...\n"<<endl; 

  //input hander
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
  if (!inputHandler) {  printf("\n ...Input handler missing!!!...\n");    return; }

  //PileUp Multi-Vertex
  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);

  //pid response object
  fPIDResponse=inputHandler->GetPIDResponse();
    
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);


  SetupEventAndTaskConfigInfo();

  if(!gGrid){
    TGrid::Connect("alien://");
  }

  SetupMCcorrectionMap();
 


  hCentvsTPCmultCuts = new TH2F("hCentvsTPCmultCuts","TPCmult Low,high",100,0,100,5,0,5);
  hCentvsTPCmultCuts->GetYaxis()->SetBinLabel(1,"mean"); 
  hCentvsTPCmultCuts->GetYaxis()->SetBinLabel(2,"sigma");   
  fListHist->Add(hCentvsTPCmultCuts);

  SetUpCentralityOutlierCut();

  fHistTOFMatchCount = new TH2F("fHistTOFMatchCount","TofMatchFlag vs mismatch Prob",10,0,10,200,-5,5);
  fListHist->Add(fHistTOFMatchCount);

  fHistEtaPtBefore = new TH2F("fHistEtaPtBefore","Eta vs pT",100,-1.25,1.25,100,0,10);
  fListHist->Add(fHistEtaPtBefore);

  fHistEtaPtAfter  = new TH2F("fHistEtaPtAfter","Eta vs pT",100,-1.25,1.25,100,0,10);
  fListHist->Add(fHistEtaPtAfter);


  Int_t gMaxTPCFB1mult  = 0;
  Int_t gMaxGlobalmult  = 0;
  Int_t gMaxTPCcorrmult = 0;
  Int_t gMaxESDtracks   = 0;
  Int_t nBinRefMult  =  200;
  Int_t nRefMultMax  = 4000;
 


  if(sNucleiTP=="pp"||sNucleiTP=="PP"){  
    gMaxGlobalmult  = 200;
    gMaxTPCFB1mult  = 200;
    gMaxTPCcorrmult = 500;
    gMaxESDtracks   = 1000;
    nBinRefMult = 100;        //change binning for pp
    nRefMultMax = 500;
    //fSkipOutlierCut = 1;
  }
  else if(sNucleiTP=="pPb"||sNucleiTP=="Pbp"||sNucleiTP=="PbP"||sNucleiTP=="PPb"){  
    gMaxGlobalmult  = 400;
    gMaxTPCFB1mult  = 400;
    gMaxTPCcorrmult = 500;
    gMaxESDtracks   = 2000;
    nBinRefMult =  100;        //change binning for pPb
    nRefMultMax = 1000;
    //fSkipOutlierCut = 1;
  }
  else{
    gMaxGlobalmult  = 4000;
    gMaxTPCFB1mult  = 4000;
    gMaxTPCcorrmult = 5000;
    gMaxESDtracks   = 20000;
    //fSkipOutlierCut =  0;
  }



  //if(bSkipPileUpCut) { fSkipOutlierCut = 1;}


  fHistTPCVsESDTrkBefore = new TH2F("fHistTPCVsESDTrkBefore","Before; TPC1; ESD trk",500,0,gMaxTPCcorrmult,200,0,gMaxESDtracks);
  fListHist->Add(fHistTPCVsESDTrkBefore);
  fHistTPCVsESDTrkAfter  = new TH2F("fHistTPCVsESDTrkAfter"," After;  TPC1; ESD trk",500,0,gMaxTPCcorrmult,200,0,gMaxESDtracks);
  fListHist->Add(fHistTPCVsESDTrkAfter);

  fHistTPCvsGlobalMultBefore = new TH2F("fHistTPCvsGlobalMultBefore","Before; Global; TPC(fb1) ",200,0,gMaxGlobalmult,200,0,gMaxTPCFB1mult);
  fListHist->Add(fHistTPCvsGlobalMultBefore);
  fHistTPCvsGlobalMultAfter  = new TH2F("fHistTPCvsGlobalMultAfter"," After;  Global; TPC(fb1) ",200,0,gMaxGlobalmult,200,0,gMaxTPCFB1mult);
  fListHist->Add(fHistTPCvsGlobalMultAfter);

  fHistGlobalVsV0MBefore = new TH2F("fHistGlobalVsV0MBefore","Before;Cent(V0M);Global",100,0,100,500,0,gMaxGlobalmult);
  fListHist->Add(fHistGlobalVsV0MBefore);
  fHistGlobalVsV0MAfter  = new TH2F("fHistGlobalVsV0MAfter"," After; Cent(V0M);Global",100,0,100,500,0,gMaxGlobalmult);
  fListHist->Add(fHistGlobalVsV0MAfter);

  fHistTPConlyVsCL1Before = new TH2F("fHistTPConlyVsCL1Before","Before;Cent(CL1); TPC(FB1)",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsCL1Before);
  fHistTPConlyVsCL1After  = new TH2F("fHistTPConlyVsCL1After","After; Cent(CL1); TPC(FB1) ",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsCL1After);

  fHistTPConlyVsV0MBefore = new TH2F("fHistTPConlyVsV0MBefore","Before;Cent(V0M); TPC(FB1)",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsV0MBefore);
  fHistTPConlyVsV0MAfter  = new TH2F("fHistTPConlyVsV0MAfter","After; Cent(V0M); TPC(FB1) ",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsV0MAfter);

  fHistRawVsCorrMultFB = new TH2F("fHistRawVsCorrMultFB",Form("FB%d;Mult_{raw};Mult_{corr}",fFilterBit),gMaxTPCFB1mult,0,gMaxTPCFB1mult,gMaxTPCcorrmult,0,gMaxTPCcorrmult);
  fListHist->Add(fHistRawVsCorrMultFB);






  // Turn off PID QA histograms
  //---------------- PID QA Histograms --------------------- 
  /*
  fHistTPCdEdxvsPBefore = new TH2F("fHistTPCdEdxvsPBefore","Before; p (GeV/c); dEdx (arb)",200,-5,5,200,0,250);
  fListHist->Add(fHistTPCdEdxvsPBefore);
  fHistTPCdEdxvsPAfter  = new TH2F("fHistTPCdEdxvsPAfter"," After;  p (GeV/c); dEdx (arb)",200,-5,5, 200,0,250);
  fListHist->Add(fHistTPCdEdxvsPAfter);

  fHistTOFBetavsPBefore = new TH2F("fHistTOFBetavsPBefore","Before; p (GeV/c); beta ",200,-5,5,100,0.0,1.2);
  fListHist->Add(fHistTOFBetavsPBefore);
  fHistTOFBetavsPAfter  = new TH2F("fHistTOFBetavsPAfter"," After;  p (GeV/c); beta ",200,-5,5,100,0.0,1.2);
  fListHist->Add(fHistTOFBetavsPAfter);

  fHistTOFMassvsPtBefore = new TH2F("fHistTOFMassvsPtBefore","Before; p_{T}(GeV/c); m^{2}(GeV^{2}/c^{4})",200,-5,5,500,-0.5,4.5);
  fListHist->Add(fHistTOFMassvsPtBefore);

  //const char *gSpecies[4] = {"Pion","Kaon","proton","Charge"};

  for(int i=0;i<3;i++){
    fHistPtwithTPCNsigma[i]  = new TH1F(Form("fHistPtwithTPCNsigma_%d",i), Form("%d;p_{T}(GeV/c))",i),200,-5,5); //i
    fListHist->Add(fHistPtwithTPCNsigma[i]);
    fHistPtwithTOFmasscut[i] = new TH1F(Form("fHistPtwithTOFmasscut_%d",i),Form("%d;p_{T}(GeV/c))",i),200,-5,5);
    fListHist->Add(fHistPtwithTOFmasscut[i]);
    fHistPtwithTOFSignal[i]  = new TH1F(Form("fHistPtwithTOFSignal_%d", i),Form("%d;p_{T}(GeV/c))",i),200,-5,5);
    fListHist->Add(fHistPtwithTOFSignal[i]);

    fHistTOFnSigmavsPtAfter[i] = new TH2F(Form("fHistTOFnSigmavsPtAfter_%d",i),Form("%d;p_{T}(GeV/c);n#sigma_{TOF}",i),200,-5,5,400,-10.0,10.0);
    fListHist->Add(fHistTOFnSigmavsPtAfter[i]);

    fHistTPCnSigmavsPtAfter[i] = new TH2F(Form("fHistTPCnSigmavsPtAfter_%d",i),Form("%d;p_{T}(GeV/c);n#sigma_{TPC}",i),200,-5,5,400,-10.0,10.0);
    fListHist->Add(fHistTPCnSigmavsPtAfter[i]);

    fHistTPCTOFnSigmavsPtAfter[i] = new TH3F(Form("fHistTPCTOFnSigmavsPtAfter_%d",i),Form("%d; p_{T}(GeV/c); n#sigma_{TPC}; n#sigma_{TOF}",i),100,0,5,400,-10,10,400,-10,10);
    fListHist->Add(fHistTPCTOFnSigmavsPtAfter[i]);

    fHistTPCdEdxvsPtPIDAfter[i] = new TH2F(Form("fHistTPCdEdxvsPtAfter_%d",i),"AfterCut; p_{T} (GeV/c); dEdx (arb)",400,0,10,200,0,250);
    fListHist->Add(fHistTPCdEdxvsPtPIDAfter[i]);
  }// PID histograms done
  */




  fCentDistBefore = new TH1F("fCentDistBefore","no Cut; Cent (%); Events ",100,0,100);
  fListHist->Add(fCentDistBefore);

  fCentDistAfter = new TH1F("fCentDistAfter","with Cut; Cent (%); Events ",100,0,100);
  fListHist->Add(fCentDistAfter);













  Double_t centRange[11]   = {0,5,10,20,30,40,50,60,70,80,90};
  //const char *gDetForEP[4] = {"V0A","V0C","TPC-A","TPC-C"};
  // 10,centRange
 //------------------- 3p correlator vs Centrality (EP method) ------------------
  for(int i=0;i<2;i++){
    //Charged:
    for(int j=0;j<4;j++){ 
     //Detector: 0 = V0A, 1 = V0C, 3 = TPCA, 4 = TPCC 
      fHist_Corr3p_EP_Norm_PN[i][j] = new TProfile(Form("fHist_Corr3p_EP_Norm_PosNeg_Mag%d_Det%d",i,j+1),Form("US, #Psi_{2} %d",j),10,centRange,"");
      fHist_Corr3p_EP_Norm_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EP_Norm_PN[i][j]);
      fHist_Corr3p_EP_Norm_PP[i][j] = new TProfile(Form("fHist_Corr3p_EP_Norm_PosPos_Mag%d_Det%d",i,j+1),Form("P-P, #Psi_{2} %d",j),10,centRange,"");
      fHist_Corr3p_EP_Norm_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EP_Norm_PP[i][j]);
      fHist_Corr3p_EP_Norm_NN[i][j] = new TProfile(Form("fHist_Corr3p_EP_Norm_NegNeg_Mag%d_Det%d",i,j+1),Form("N-N, #Psi_{2}, %d",j),10,centRange,"");
      fHist_Corr3p_EP_Norm_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EP_Norm_NN[i][j]);
    }
    //EP Resolution:
    for(int j=0;j<4;j++){
    //Det: 0 = v0c-v0a, 1 = v0a-TPC, 2 = v0c-TPC, 3 =TPC-A TPC-C
      fHist_Reso2n_EP_Norm_Det[i][j]  = new TProfile(Form("fHist_Reso2n_EP_Norm_Mag%d_DetComb%d",i,j+1),"Event plane Resolution",10,centRange,"");
      fHist_Reso2n_EP_Norm_Det[i][j]->Sumw2();
      fListHist->Add(fHist_Reso2n_EP_Norm_Det[i][j]);
    }
  }//magfield loop




  Char_t  name[100];
  Char_t title[100];
 
  Double_t pTRange[21] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0}; 
//Double_t pTRange[11] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};

  //Charge:
  for(Int_t i=0;i<2;i++){ 
    for(Int_t j=0;j<6;j++){
      sprintf(name,"fHist_Corr3p_pTSum_EP_V0A_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0A_PN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0A_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTSum_EP_V0A_PN[i][j]);

      sprintf(name,"fHist_Corr3p_pTSum_EP_V0A_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0A_PP[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0A_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTSum_EP_V0A_PP[i][j]);

      sprintf(name,"fHist_Corr3p_pTSum_EP_V0A_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0A_NN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0A_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTSum_EP_V0A_NN[i][j]);
      //-----v0c----
      sprintf(name,"fHist_Corr3p_pTSum_EP_V0C_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0C_PN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0C_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTSum_EP_V0C_PN[i][j]);

      sprintf(name,"fHist_Corr3p_pTSum_EP_V0C_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0C_PP[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0C_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTSum_EP_V0C_PP[i][j]);

      sprintf(name,"fHist_Corr3p_pTSum_EP_V0C_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0C_NN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0C_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTSum_EP_V0C_NN[i][j]);
    }
  }

  for(Int_t i=0;i<2;i++){ 
    for(Int_t j=0;j<6;j++){
      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0A_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0A_PN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0A_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTDiff_EP_V0A_PN[i][j]);

      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0A_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0A_PP[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0A_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTDiff_EP_V0A_PP[i][j]);

      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0A_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0A_NN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0A_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTDiff_EP_V0A_NN[i][j]);
      //-----v0c----
      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0C_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0C_PN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0C_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTDiff_EP_V0C_PN[i][j]);

      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0C_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0C_PP[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0C_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTDiff_EP_V0C_PP[i][j]);

      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0C_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0C_NN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0C_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_pTDiff_EP_V0C_NN[i][j]);
    }
  }
 
  //Double_t EtaRange[9] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6}; //Use this after tests done
  //Now Eta binning: 16,0,1.6 for test

  for(Int_t i=0;i<2;i++){ 
    for(Int_t j=0;j<6;j++){
      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0A_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0A_PN[i][j] = new TProfile(name,title,8,0,1.6,"");
      fHist_Corr3p_EtaDiff_EP_V0A_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EtaDiff_EP_V0A_PN[i][j]);

      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0A_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0A_PP[i][j] = new TProfile(name,title,8,0,1.6,"");
      fHist_Corr3p_EtaDiff_EP_V0A_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EtaDiff_EP_V0A_PP[i][j]);

      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0A_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0A_NN[i][j] = new TProfile(name,title,8,0,1.6,"");
      fHist_Corr3p_EtaDiff_EP_V0A_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EtaDiff_EP_V0A_NN[i][j]);
      //-----v0c----
      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0C_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0C_PN[i][j] = new TProfile(name,title,8,0,1.6,""); 
      fHist_Corr3p_EtaDiff_EP_V0C_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EtaDiff_EP_V0C_PN[i][j]);

      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0C_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0C_PP[i][j] = new TProfile(name,title,8,0,1.6,"");
      fHist_Corr3p_EtaDiff_EP_V0C_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EtaDiff_EP_V0C_PP[i][j]);

      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0C_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0C_NN[i][j] = new TProfile(name,title,8,0,1.6,"");
      fHist_Corr3p_EtaDiff_EP_V0C_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EtaDiff_EP_V0C_NN[i][j]);
    }
  }




  // Alex's correction histograms: 
    
  fpX1X2PosT1 =  new TProfile("fpX1X2PosT1", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2PosT1);
  fpX1X3PosT1 =  new TProfile("fpX1X3PosT1", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3PosT1);
  fpX1X2PosT2 =  new TProfile("fpX1X2PosT2", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2PosT2);
  fpX1X3PosT2 =  new TProfile("fpX1X3PosT2", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3PosT2);
  fpX1X2PosT3 =  new TProfile("fpX1X2PosT3", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2PosT3);
  fpX1X3PosT3 =  new TProfile("fpX1X3PosT3", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3PosT3);
  fpX1X2PosT4 =  new TProfile("fpX1X2PosT4", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2PosT4);
  fpX1X3PosT4 =  new TProfile("fpX1X3PosT4", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3PosT4);
     
     
  fpX1X2NegT1 =  new TProfile("fpX1X2NegT1", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2NegT1);
  fpX1X3NegT1 =  new TProfile("fpX1X3NegT1", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3NegT1);
  fpX1X2NegT2 =  new TProfile("fpX1X2NegT2", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2NegT2);
  fpX1X3NegT2 =  new TProfile("fpX1X3NegT2", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3NegT2);
  fpX1X2NegT3 =  new TProfile("fpX1X2NegT3", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2NegT3);
  fpX1X3NegT3 =  new TProfile("fpX1X3NegT3", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3NegT3);
  fpX1X2NegT4 =  new TProfile("fpX1X2NegT4", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2NegT4);
  fpX1X3NegT4 =  new TProfile("fpX1X3NegT4", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3NegT4);
     
  fpX1X2OppT1 =  new TProfile("fpX1X2OppT1", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2OppT1);
  fpX1X2OppT2 =  new TProfile("fpX1X2OppT2", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2OppT2);
  fpX1X2OppT3 =  new TProfile("fpX1X2OppT3", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2OppT3);
  fpX1X2OppT4 =  new TProfile("fpX1X2OppT4", "; centrality percentile; #LT X1X2 #GT", 10,centRange);
  fListHist->Add(fpX1X2OppT4);


  fpX1X3PosT1EP2 =  new TProfile("fpX1X3PosT1EP2", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3PosT1EP2);
  fpX1X3PosT2EP2 =  new TProfile("fpX1X3PosT2EP2", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3PosT2EP2);
  fpX1X3PosT3EP2 =  new TProfile("fpX1X3PosT3EP2", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3PosT3EP2);
  fpX1X3PosT4EP2 =  new TProfile("fpX1X3PosT4EP2", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3PosT4EP2);

  fpX1X3NegT1EP2 =  new TProfile("fpX1X3NegT1EP2", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3NegT1EP2);
  fpX1X3NegT2EP2 =  new TProfile("fpX1X3NegT2EP2", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3NegT2EP2);
  fpX1X3NegT3EP2 =  new TProfile("fpX1X3NegT3EP2", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3NegT3EP2);
  fpX1X3NegT4EP2 =  new TProfile("fpX1X3NegT4EP2", "; centrality percentile; #LT X1X3 #GT", 10,centRange);
  fListHist->Add(fpX1X3NegT4EP2);
    
  /// Single terms:
  fpX1PosCosT1 =  new TProfile("fpX1PosCosT1", "; centrality percentile; #LT <cos(#varphi_{+})> #GT", 10, centRange);
  fListHist->Add(fpX1PosCosT1);
  fpX1NegCosT1 =  new TProfile("fpX1NegCosT1", "; centrality percentile; #LT <cos(#varphi_{-})> #GT", 10, centRange);
  fListHist->Add(fpX1NegCosT1);
  fpX1PosSinT1 =  new TProfile("fpX1PosSinT1", "; centrality percentile; #LT <sin(#varphi_{+})> #GT", 10, centRange);
  fListHist->Add(fpX1PosSinT1);
  fpX1NegSinT1 =  new TProfile("fpX1NegSinT1", "; centrality percentile; #LT <sin(#varphi_{-})> #GT", 10, centRange);
  fListHist->Add(fpX1NegSinT1);

  fpX1EventT1EP1 =  new TProfile("fpX1EventT1EP1", "; centrality percentile; #LT <cos(N#Psi_N)> #GT", 10, centRange);
  fListHist->Add(fpX1EventT1EP1);
  fpX1EventT1EP2 =  new TProfile("fpX1EventT1EP2", "; centrality percentile; #LT <cos(N#Psi_N)> #GT", 10, centRange);
  fListHist->Add(fpX1EventT1EP2);
  fpY1EventT1EP1 =  new TProfile("fpY1EventT1EP1", "; centrality percentile; #LT <sin(N#Psi_N)> #GT", 10, centRange);
  fListHist->Add(fpY1EventT1EP1);
  fpY1EventT1EP2 =  new TProfile("fpY1EventT1EP2", "; centrality percentile; #LT <sin(N#Psi_N)> #GT", 10, centRange);
  fListHist->Add(fpY1EventT1EP2);
  























  

  //---- to store NUA and calib histograms -----
  TList *fListNUACalib = new TList();
  fListNUACalib->SetName("fListNUACalib");
  fListNUACalib->SetOwner(kTRUE);


  Double_t truncPi = 3.1416;

  //-------------- Event plane distributions  --------------
  fHV0AEventPlaneVsCent = new TH2F("fHV0AEventPlaneVsCent",Form("Psi %d from V0A", gPsiN), 10,centRange,50,-0.0,truncPi);
  fListNUACalib->Add(fHV0AEventPlaneVsCent);
  
  fHV0CEventPlaneVsCent = new TH2F("fHV0CEventPlaneVsCent",Form("Psi %d from V0C", gPsiN), 10,centRange,50,-0.0,truncPi);
  fListNUACalib->Add(fHV0CEventPlaneVsCent);

  fHTPCAEventPlaneVsCent = new TH2F("fHTPCAEventPlaneVsCent",Form("Psi %d from pos eta",gPsiN),10,centRange,50,-0.0,truncPi);
  fListNUACalib->Add(fHTPCAEventPlaneVsCent);

  fHTPCCEventPlaneVsCent = new TH2F("fHTPCCEventPlaneVsCent",Form("Psi %d from neg eta",gPsiN),10,centRange,50,-0.0,truncPi);
  fListNUACalib->Add(fHTPCCEventPlaneVsCent);

  fHTPCEventPlaneVsCent = new TH2F("fHTPCEventPlaneVsCent",Form("Psi %d from Full TPC",gPsiN),10,centRange,50,-0.0,truncPi);
  fListNUACalib->Add(fHTPCEventPlaneVsCent);



  //-------------------- QA: how does the eta-phi looks after NUA correction ----------------
  fQAEtaPhiAfterNUA = new TH2F("fQAEtaPhiAfterNUA","eta vs phi with NUA corr",50,0,6.283185,16,-0.8,0.8);
  fListNUACalib->Add(fQAEtaPhiAfterNUA);

  fQAEtaPhiAfterNUAPion = new TH2F("fQAEtaPhiAfterNUAPion","Pion eta vs phi with NUA corr",50,0,6.283185,16,-0.8,0.8);
  fListNUACalib->Add(fQAEtaPhiAfterNUAPion);

  fQAEtaPhiAfterNUAKaon = new TH2F("fQAEtaPhiAfterNUAKaon","Kaon eta vs phi with NUA corr",50,0,6.283185,16,-0.8,0.8);
  fListNUACalib->Add(fQAEtaPhiAfterNUAKaon);

  fQAEtaPhiAfterNUAProton = new TH2F("fQAEtaPhiAfterNUAProton","Proton eta vs phi with NUA corr",50,0,6.283185,16,-0.8,0.8);
  fListNUACalib->Add(fQAEtaPhiAfterNUAProton);



  //----------------- ZPA,ZNA Calibration hist: ---------------------
  fZPASignalPerChVsCent = new TH2F("fZPASignalPerChVsCent","En-ZPA, Cent, Vz, ch",91,0,91,5,0,5); 
  fZPCSignalPerChVsCent = new TH2F("fZPCSignalPerChVsCent","En-ZPC, Cent, Vz, ch",91,0,91,5,0,5); 
  fZNASignalPerChVsCent = new TH2F("fZNASignalPerChVsCent","En-ZNA, Cent, Vz, ch",91,0,91,5,0,5); 
  fZNCSignalPerChVsCent = new TH2F("fZNCSignalPerChVsCent","En-ZNC, Cent, Vz, ch",91,0,91,5,0,5);

  fListNUACalib->Add(fZPASignalPerChVsCent);
  fListNUACalib->Add(fZPCSignalPerChVsCent);
  fListNUACalib->Add(fZNASignalPerChVsCent);
  fListNUACalib->Add(fZNCSignalPerChVsCent);

  fCentDistVsVz = new TProfile("fCentDistVsVz","<Cent> vs Vz.",80,-10,10);

  fListNUACalib->Add(fCentDistVsVz);




  //----------------- V0 Calibration hist: ---------------------
  fV0MultChVsRun = new TH2F("fV0MultChVsRun","1-32 V0C, 33-64 V0A",64,0,64,90,0,90);
  fListNUACalib->Add(fV0MultChVsRun);

  //const char *sCorrect[2]={"wo Corr","w Corr"};
  Int_t isCorr = 0;

  if(fListV0MCorr){
    isCorr = 1;
  }


  fV0AQ2xVsCentRun = new TProfile("fV0ACos2nVsCentRun",Form("<Cos2> vs cent (%d)",isCorr),90,0,90,""); //sCorrect[isCorr]
  fListNUACalib->Add(fV0AQ2xVsCentRun);
  fV0AQ2yVsCentRun = new TProfile("fV0ASin2nVsCentRun",Form("<Sin2> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fV0AQ2yVsCentRun);
  fV0CQ2xVsCentRun = new TProfile("fV0CCos2nVsCentRun",Form("<Cos2> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fV0CQ2xVsCentRun);
  fV0CQ2yVsCentRun = new TProfile("fV0CSin2nVsCentRun",Form("<Sin2> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fV0CQ2yVsCentRun);

  fV0AQ3xVsCentRun = new TProfile("fV0ACos3nVsCentRun",Form("<Cos3> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fV0AQ3xVsCentRun);
  fV0AQ3yVsCentRun = new TProfile("fV0ASin3nVsCentRun",Form("<Sin3> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fV0AQ3yVsCentRun);
  fV0CQ3xVsCentRun = new TProfile("fV0CCos3nVsCentRun",Form("<Cos3> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fV0CQ3xVsCentRun);
  fV0CQ3yVsCentRun = new TProfile("fV0CSin3nVsCentRun",Form("<Sin3> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fV0CQ3yVsCentRun);

  isCorr = 1;
  if(fListNUACorr){
    cout<<"\n =========> NUA file found for NUA correction <========== \n";
    isCorr = 0;
  }
  //------------------- TPC Qvector Recentering Histograms --------------
  fTPCAQ2xVsCentRun = new TProfile("fTPCACos2nVsCentRun",Form("<Cos2> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCAQ2xVsCentRun);
  fTPCAQ2yVsCentRun = new TProfile("fTPCASin2nVsCentRun",Form("<Sin2> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCAQ2yVsCentRun);
  fTPCCQ2xVsCentRun = new TProfile("fTPCCCos2nVsCentRun",Form("<Cos2> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCCQ2xVsCentRun);
  fTPCCQ2yVsCentRun = new TProfile("fTPCCSin2nVsCentRun",Form("<Sin2> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCCQ2yVsCentRun);

  fTPCAQ3xVsCentRun = new TProfile("fTPCACos3nVsCentRun",Form("<Cos3> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCAQ3xVsCentRun);
  fTPCAQ3yVsCentRun = new TProfile("fTPCASin3nVsCentRun",Form("<Sin3> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCAQ3yVsCentRun);
  fTPCCQ3xVsCentRun = new TProfile("fTPCCCos3nVsCentRun",Form("<Cos3> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCCQ3xVsCentRun);
  fTPCCQ3yVsCentRun = new TProfile("fTPCCSin3nVsCentRun",Form("<Sin3> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCCQ3yVsCentRun);

  fTPCAQ4xVsCentRun = new TProfile("fTPCACos4nVsCentRun",Form("<Cos4> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCAQ4xVsCentRun);
  fTPCAQ4yVsCentRun = new TProfile("fTPCASin4nVsCentRun",Form("<Sin4> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCAQ4yVsCentRun);
  fTPCCQ4xVsCentRun = new TProfile("fTPCCCos4nVsCentRun",Form("<Cos4> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCCQ4xVsCentRun);
  fTPCCQ4yVsCentRun = new TProfile("fTPCCSin4nVsCentRun",Form("<Sin4> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCCQ4yVsCentRun);


  fTPCFQ2xVsCentRun = new TProfile("fTPCFQ2xVsCentRun",Form("<Cos2> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCFQ2xVsCentRun);
  fTPCFQ2yVsCentRun = new TProfile("fTPCFQ2yVsCentRun",Form("<Sin2> vs cent (%d)",isCorr),90,0,90,"");
  fListNUACalib->Add(fTPCFQ2yVsCentRun);

  //-------------------------------------------------------------------------------

  
  
  
 //-------------------------- Define NUA Hist for PID -----------------------------
  Int_t gCentForNUA[6] = {0,5,10,20,40,90};
  //Char_t  name[100];
  //Char_t title[100];

  for(int i=0;i<4;i++){
    for(int j=0;j<5;j++){
      sprintf(name,"fHistEtaPhiVz_%d_Pos_Cent%d_Run%d",i,j,1); //gSpecies[i]
      sprintf(title,"eta,phi,Vz %dPos, Cent%d-%d, FB %d",i,gCentForNUA[j],gCentForNUA[j+1],fFilterBit);
      fHist3DEtaPhiVz_Pos_Run[i][j] = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
      fListNUACalib->Add(fHist3DEtaPhiVz_Pos_Run[i][j]);

      sprintf(name,"fHistEtaPhiVz_%d_Neg_Cent%d_Run%d",i,j,1); //gSpecies[i]
      sprintf(title,"eta,phi,Vz %dNeg, Cent%d-%d, FB %d",i,gCentForNUA[j],gCentForNUA[j+1],fFilterBit);
      fHist3DEtaPhiVz_Neg_Run[i][j] = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
      fListNUACalib->Add(fHist3DEtaPhiVz_Neg_Run[i][j]);
    }
  }
 //---------------------------------------------------------------------------------


  fListHist->Add(fListNUACalib);

  PostData(1,fListHist);
  cout<<"\n.........UserCreateOutputObject called.........\n fFilterBit = "<<fFilterBit<<" CentMax = "<<fCentralityPercentMax;
  cout<<" PU C = "<<fPileUpConstParm<<" gN = "<<gN<<" gM = "<<gM<<" PsiN = "<<gPsiN<<"\n\n"<<endl;
  //cout<<" ***********  checking my commit...!! ***************** \n\n ";
}

















//______________________________________________________________________
void AliAnalysisTaskCMEV0PID::UserExec(Option_t*) {
  //debug only
  //cout<<"\n Info:UserExec() called ..!!!\n";
  //watch.Start(kTRUE);  
  //if(fEventCount==501)  return;


  Float_t stepCount = 0.5;

  fHistEventCount->Fill(stepCount); //1
  stepCount++;

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    
  if(!(fESD || fAOD)){ printf("ERROR: fESD & fAOD not available\n"); return;  }
    
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());

  if (!fVevent) { printf("ERROR: fVevent not available\n");  return;  }
   
  fHistEventCount->Fill(stepCount); //2
  stepCount++;






  //--------- Check if I have PID response object --------
  
  if(!fPIDResponse){
     AliAnalysisManager   *mgr=AliAnalysisManager::GetAnalysisManager();
     AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
     if(inputHandler) fPIDResponse=inputHandler->GetPIDResponse();
     if(!fPIDResponse){
       printf("\n\n...... PIDResponse object not found..... \n\n");
       return;
     }
  }
 




  //-------------- Vtx cuts ---------------
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    
  Double_t pVtxZ = -999;
  pVtxZ = pVtx->GetZ();
    
//if(TMath::Abs(pVtxZ)>10.) return;
  //User defined cut:
  if(pVtxZ<fMinVzCut || pVtxZ>fMaxVzCut ) return;


  fHistEventCount->Fill(stepCount); //3
  stepCount++;



  Float_t centrality = -99.0;
  Float_t centrV0M   = -99.0;
  Float_t centrCL1   = -99.0;



  //---------- Centrality Estimators  -------------
  AliCentrality* Alicentr = ((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP();   // For Run1, 2010 data

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");  // Must never comment this
  if(!fMultSelection) { printf("\n\n **WARNING** \n::UserExec() AliMultSelection object not found.\n\n"); exit(1); }

  
  if(sNucleiTP=="PbPb2010" || sNucleiTP=="2010") {
    if(Alicentr){

      centrV0M  = Alicentr->GetCentralityPercentile("V0M"); 
      centrCL1  = Alicentr->GetCentralityPercentile("CL1");

      if(sCentrEstimator=="V0M" || sCentrEstimator=="V0"){
	centrality = centrV0M;
      }
      else if(sCentrEstimator=="CL1"){
	centrality = centrCL1;
      }
      else if(sCentrEstimator=="V0C"){
	centrality = Alicentr->GetCentralityPercentile("V0C");
      }
      else if(sCentrEstimator=="V0A"){
	centrality = Alicentr->GetCentralityPercentile("V0A");
      }
      else if(sCentrEstimator=="TRK"){
	centrality = Alicentr->GetCentralityPercentile("TRK");
      }
    }
  }
  else{  // fall back to MultSelection if other than 2010 data

    centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
    centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");

    if(sCentrEstimator=="V0M" || sCentrEstimator=="V0"){
      centrality = centrV0M;
    }
    else if(sCentrEstimator=="CL1"){
      centrality = centrCL1;
    }
    else if(sCentrEstimator=="V0C"){
      centrality = fMultSelection->GetMultiplicityPercentile("V0C");
    }
    else if(sCentrEstimator=="V0A"){
      centrality = fMultSelection->GetMultiplicityPercentile("V0A");
    }
    else if(sCentrEstimator=="TRK"){
      centrality = fMultSelection->GetMultiplicityPercentile("TRK");
    }
  }
 
  /* 
  //--------- Centrality Estimators -------------
  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");

  if(sCentrEstimator=="V0M" || sCentrEstimator=="V0"){
    centrality = centrV0M;
  }
  else if(sCentrEstimator=="CL1"){
    centrality = centrCL1;
  }
  else if(sCentrEstimator=="V0C"){
    centrality = fMultSelection->GetMultiplicityPercentile("V0C");
  }
  else if(sCentrEstimator=="V0A"){
    centrality = fMultSelection->GetMultiplicityPercentile("V0A");
  }
  else if(sCentrEstimator=="TRK"){
    centrality = fMultSelection->GetMultiplicityPercentile("TRK");
  }  */


  

  if(centrality<fCentralityPercentMin || centrality>fCentralityPercentMax){ 
    return;
  }

  fHistEventCount->Fill(stepCount); //4
  stepCount++;

  fCentDistBefore->Fill(centrality);
  

  



  Int_t ntracks=fAOD->GetNumberOfTracks();
  if(ntracks<2) return;              // Check this cut....!!!
    
  fHistEventCount->Fill(stepCount); //5
  stepCount++;






  Int_t cent10bin = -1;
  Int_t cIndex = -1;
  //cent10bin = GetCentralityScaled0to10(centrality); //Centrality in 0-10 scale

  if(centrality<5.0) {
    cent10bin  = 0; 
  }
  else if(centrality>=5.0 && centrality<10){
    cent10bin  = 1;
  }
  else if(centrality>=10.0) {
    cent10bin = abs(centrality/10.0)+1;
  }

  cIndex = cent10bin;

//Centrality array index for NUA correcion
  Int_t cForNUA = 0;  

  if(centrality<5.0) {
    cForNUA = 0;
  }
  else if(centrality>=5.0 && centrality<10){
    cForNUA = 1; // 1=5-10,
  }
  else if(centrality>=10.0 && centrality<20) {
    cForNUA = 2; // 2 = 10-20,
  }
  else if(centrality>=20 && centrality<40){ 
    cForNUA = 3; // 3=20-40
  }
  else if(centrality>=40){
    cForNUA = 4; // 4=40-90
  }






  //---------- Magnetic field --------
  Double_t fMagField = fAOD->GetMagneticField();

  const Int_t QAindex = (fMagField > 0) ? 1 : 0;
  //---------------------------------






  //Load NUA and V0M correction map run by run:
  Int_t runNumber = fAOD->GetRunNumber();
 
  if(runNumber!=fOldRunNum) {
   
    GetNUACorrectionHist(runNumber);
    
    if(bV0MGainCorr) {
      GetV0MCorrectionHist(runNumber);
    }

    fOldRunNum = runNumber;
  } 
  //------------------------------------------
 






  //----- Event Plane variables:-------
  Double_t PsiNV0A = 0;
  Double_t PsiNV0C = 0;

  Double_t PsiNTPCA = 0; // eta <0 
  Double_t PsiNTPCC = 0; // eta >0
  Double_t PsiNTPCF = 0; // Full TPC

  Double_t sumTPCQn2x[3] = {0,0,0};  //[0]= eta<0; [1]= eta>0; [2]= -0.8 < eta < 0.8
  Double_t sumTPCQn2y[3] = {0,0,0};
  //Double_t sumTPCQn3x[3] = {0,0,0};  //[0]= eta<0; [1]= eta>0; [2]= -0.8 < eta < 0.8
  //Double_t sumTPCQn3y[3] = {0,0,0};
  //Double_t sumTPCQn4x[3] = {0,0,0};  //[0]= eta<0; [1]= eta>0; [2]= -0.8 < eta < 0.8
  //Double_t sumTPCQn4y[3] = {0,0,0};
  //------------------------------------





 

 //Variables for MC tracking correction 
  Int_t   ptBinMC = 1;
  Int_t   iBinNUA = 1;
  Double_t  ptWgtMC = 1.0;
  Double_t  WgtNUA  = 1.0;
  Double_t  ptTrk   = 0.1;
  Double_t  dEdx    = 0.0;
  Double_t  Chi2Trk = 0.0;


 //-------------- Track loop for outlier and PileUp cut -------------------

 //---------------- a dobrin --------------

  Bool_t bIsPileup=kFALSE;

  Int_t isPileup = fAOD->IsPileupFromSPD(3);

  if(isPileup != 0) {
    fHistPileUpCount->Fill(0.5);
    bIsPileup=kTRUE;          
  }
  else if(PileUpMultiVertex(fAOD)) {
    fHistPileUpCount->Fill(1.5);
    bIsPileup=kTRUE;
  }
  else if(((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) {
    fHistPileUpCount->Fill(2.5);
    bIsPileup=kTRUE;
  }
  else if(fAOD->IsIncompleteDAQ())  {
    fHistPileUpCount->Fill(3.5);
    bIsPileup=kTRUE;
  }
  else if(fabs(centrV0M-centrCL1)> 5.0)  {//default: 7.5
//else if(fabs(centrV0M-centrCL1)> 7.5)  {//default: 7.5
    fHistPileUpCount->Fill(4.5);
    bIsPileup=kTRUE;
  }

  // check vertex consistency
  const AliAODVertex* vtTrc = fAOD->GetPrimaryVertex();
  const AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();

  if(vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
    fHistPileUpCount->Fill(5.5);
    bIsPileup=kTRUE;
  }

  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);

  double dz = vtTrc->GetZ() - vtSPD->GetZ();

  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = dz/errTot;
  double nsigTrc = dz/errTrc;

  if(TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)  {
    fHistPileUpCount->Fill(6.5);
    bIsPileup=kTRUE;
  }


  /*
  Float_t nSigPionTPC[40000]   = {0.,};
  Float_t nSigKaonTPC[40000]   = {0.,};
  Float_t nSigProtonTPC[40000] = {0.,};    
  Float_t nSigPionTOF[40000]   = {0.,};
  Float_t nSigKaonTOF[40000]   = {0.,};
  Float_t nSigProtonTOF[40000] = {0.,};  
  */

  if(ntracks > 40000)           return;     //Dont break segment for higher tracks:   
  

  /* 
  std::vector<float> nSigPionTPC;   
  std::vector<float> nSigKaonTPC;   
  std::vector<float> nSigProtonTPC; 
  std::vector<float> nSigPionTOF;   
  std::vector<float> nSigKaonTOF;   
  std::vector<float> nSigProtonTOF;  

  nSigPionTPC.reserve(10000);
  nSigKaonTPC.reserve(10000);
  nSigProtonTPC.reserve(10000);
  nSigPionTOF.reserve(10000);
  nSigKaonTOF.reserve(10000);
  nSigProtonTOF.reserve(10000); */
 



  Float_t multTPC     = 0;    // tpc mult estimate
  //Float_t RefMultRaw  = 0;    // tpc mult estimate
  //Float_t RefMultCorr = 0;    // tpc mult estimate
  Float_t RefMultRawFB = 0;
  Float_t RefMultCorrFB= 0;

  Float_t multTPCAll  = 0;    // tpc mult estimate
  Float_t multGlobal  = 0; // global multiplicity

  Int_t    multEtaNeg, multEtaPos, multEtaFull;
  Double_t SumWEtaNeg, SumWEtaPos, SumWEtaFull;

  Int_t   ChTrk;
  Int_t nClustTPC;

  Double_t etaTrk, phiTrk;  //  Never define eta as float, always double.!!

  //Double_t fMaxPtEP  =  5.0;
  //Double_t fMinPtEP  =  0.2;
  Double_t fMaxEtaEP =  0.8;
  Double_t fMinEtaEP = -0.8;

  multEtaNeg = 0;
  multEtaPos = 0;
  multEtaFull= 0;
  SumWEtaNeg = 0;
  SumWEtaPos = 0;
  SumWEtaFull= 0;



  for(Int_t iTrack = 0; iTrack < ntracks; iTrack++) { //-------------------------

    AliAODTrack* AODtrack =dynamic_cast<AliAODTrack*>(fVevent->GetTrack(iTrack));
    if(!AODtrack) continue;


   //------ Store PID info in array to use later:---------
   //Vector method: // taking same time as arrays
   /*
    nSigPionTPC.push_back(fPIDResponse->NumberOfSigmasTPC(AODtrack,   AliPID::kPion));
    nSigKaonTPC.push_back(fPIDResponse->NumberOfSigmasTPC( AODtrack,  AliPID::kKaon));
    nSigProtonTPC.push_back(fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kProton));
    nSigPionTOF.push_back(fPIDResponse->NumberOfSigmasTOF(AODtrack,   AliPID::kPion));
    nSigKaonTOF.push_back(fPIDResponse->NumberOfSigmasTOF( AODtrack,  AliPID::kKaon));
    nSigProtonTOF.push_back(fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kProton)); */

    //I used these Arrays for PID, now commented as not enough statistics: 28/09/18
    /*
    nSigPionTPC[iTrack]   = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kPion);
    nSigKaonTPC[iTrack]   = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kKaon);
    nSigProtonTPC[iTrack] = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kProton);
    nSigPionTOF[iTrack]   = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kPion);
    nSigKaonTOF[iTrack]   = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kKaon);
    nSigProtonTOF[iTrack] = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kProton);  
    */
    //-----------------------------------------------------


    ptTrk  = AODtrack->Pt();
    etaTrk = AODtrack->Eta();
    Chi2Trk = AODtrack->Chi2perNDF();
    nClustTPC = AODtrack->GetTPCNcls();
  //dEdx   = AODtrack->GetDetPid()->GetTPCsignal(); // This one breaks the code if called before checking FilterBit.


   //Track cuts for POIs: Same to be used for EP calc.
  //if((dPt2 > fMaxPtCut) || (dPt2 < fMinPtCut) || (dEta2 > fMaxEtaCut) || (dEta2 < fMinEtaCut) || (dEdx2 < fdEdxMin) || (track2->GetTPCNcls() < fTPCclustMin) || (Chi2Trk2 < fTrkChi2Min) || (Chi2Trk2 > 4.0) || (track2->DCA() > fDCAxyMax)|| (track2->ZAtDCA() > fDCAzMax) || !(TMath::Abs(gCharge2)))

    //cuts for EP calculation:
    if(AODtrack->TestFilterBit(fFilterBit)){

      phiTrk = AODtrack->Phi();
      dEdx   = AODtrack->GetDetPid()->GetTPCsignal(); 
      ChTrk  = AODtrack->Charge();

      if(gPsiN > 0 && (ptTrk <= fMaxPtCut) && (ptTrk >= fMinPtCut) && (etaTrk <= fMaxEtaEP) && (etaTrk >= fMinEtaEP) && (dEdx >= fdEdxMin) && (nClustTPC >= fTPCclustMin) && (AODtrack->DCA() <= fDCAxyMax) && (AODtrack->ZAtDCA() <= fDCAzMax) &&  (Chi2Trk >= fTrkChi2Min) && (Chi2Trk <= 4.0) && TMath::Abs(ChTrk))
      {

	ptWgtMC = 1.0;

	if(fFB_Efficiency_Cent[cent10bin]){
	  ptBinMC = fFB_Efficiency_Cent[cent10bin]->FindBin(ptTrk);    //Charge independent MC correction atm.
	  ptWgtMC = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBinMC);
	}
	  
	RefMultRawFB++;
	RefMultCorrFB += ptWgtMC;


	//------> Get NUA weights for EP <----------
	WgtNUA = 1.0;

	if(ChTrk>0){
	  if(fHCorrectNUApos[cForNUA]){
	    iBinNUA = fHCorrectNUApos[cForNUA]->FindBin(pVtxZ,phiTrk,etaTrk);
	    WgtNUA  = fHCorrectNUApos[cForNUA]->GetBinContent(iBinNUA);
	  }
	  //else{ WgtNUA = 1.0; }
	}
	else{
	  if(fHCorrectNUAneg[cForNUA]){
	    iBinNUA = fHCorrectNUAneg[cForNUA]->FindBin(pVtxZ,phiTrk,etaTrk);
	    WgtNUA  = fHCorrectNUAneg[cForNUA]->GetBinContent(iBinNUA);  
	  }
	  //else{ WgtNUA = 1.0; }
	}
	    
	if(etaTrk < -0.05){
	  sumTPCQn2x[0] += WgtNUA*TMath::Cos(gPsiN*phiTrk);
	  sumTPCQn2y[0] += WgtNUA*TMath::Sin(gPsiN*phiTrk);
	  //sumTPCQn3x[0] += WgtNUA*TMath::Cos(3*phiTrk);
	  //sumTPCQn3y[0] += WgtNUA*TMath::Sin(3*phiTrk);
	  //sumTPCQn4x[0] += WgtNUA*TMath::Cos(4*phiTrk);
	  //sumTPCQn4y[0] += WgtNUA*TMath::Sin(4*phiTrk);
	  multEtaNeg++;
	  SumWEtaNeg += WgtNUA;
	}
	else if(etaTrk > 0.05){
	  sumTPCQn2x[1] += WgtNUA*TMath::Cos(gPsiN*phiTrk);
	  sumTPCQn2y[1] += WgtNUA*TMath::Sin(gPsiN*phiTrk);
	  //sumTPCQn3x[1] += WgtNUA*TMath::Cos(3*phiTrk);
	  //sumTPCQn3y[1] += WgtNUA*TMath::Sin(3*phiTrk);
	  //sumTPCQn4x[1] += WgtNUA*TMath::Cos(4*phiTrk);
	  //sumTPCQn4y[1] += WgtNUA*TMath::Sin(4*phiTrk);
	  multEtaPos++;
	  SumWEtaPos += WgtNUA;
	}
	sumTPCQn2x[2] += WgtNUA*TMath::Cos(gPsiN*phiTrk);
	sumTPCQn2y[2] += WgtNUA*TMath::Sin(gPsiN*phiTrk);
	//sumTPCQn3x[2] += WgtNUA*TMath::Cos(3*phiTrk);
	//sumTPCQn3y[2] += WgtNUA*TMath::Sin(3*phiTrk);
	//sumTPCQn4x[2] += WgtNUA*TMath::Cos(4*phiTrk);
	//sumTPCQn4y[2] += WgtNUA*TMath::Sin(4*phiTrk);
	multEtaFull++;
	SumWEtaFull += WgtNUA;
      }// track cuts      
    } // AOD fiter bit


    if(AODtrack->TestFilterBit(128)) multTPCAll++; // A. Dobrin TPC vs ESD PileUp Cut.

    if(!(AODtrack->TestFilterBit(1))) continue;    // cuts for Outlier as in FlowEvent Task

    dEdx   = AODtrack->GetDetPid()->GetTPCsignal(); 

    if((ptTrk < 0.2) || (ptTrk > 5.0) || (TMath::Abs(etaTrk) > 0.8) || (nClustTPC < 70) || (dEdx < 10.0) || (Chi2Trk < 0.1)) continue;

    if(AODtrack->GetDetPid() && Chi2Trk > 0.2) multTPC++;
   
    if(!AODtrack->TestFilterBit(16) || AODtrack->Chi2perNDF() < 0.1) continue;
    Double_t b[2] = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};
    AliAODTrack copy(*AODtrack);
    if(copy.PropagateToDCA(fVevent->GetPrimaryVertex(), fVevent->GetMagneticField(), 100., b, bCov) && TMath::Abs(b[0]) < 0.3 && TMath::Abs(b[1]) < 0.3) multGlobal++;

  }//--- track loop outlier/PileUp ----






  Int_t multEsd = ((AliAODHeader*)fAOD->GetHeader())->GetNumberOfESDTracks();
  Float_t multESDTPCDiff = (Float_t) multEsd - fPileUpSlopeParm*multTPCAll;

  //cout<<" Info:UserExec() called ... I am after PU cut...  event = "<<fEventCount<<" \n";

  if(multESDTPCDiff > fPileUpConstParm) { 
    fHistPileUpCount->Fill(7.5);
    bIsPileup=kTRUE;
  }
  else if(bIsPileup==kFALSE) {
    
   if(!fMultSelection->GetThisEventIsNotPileup()){
      fHistMultSelPUCount->Fill(0.5);
      bIsPileup=kTRUE;
    }
    if(!fMultSelection->GetThisEventIsNotPileupMV()){
      fHistMultSelPUCount->Fill(1.5);
      bIsPileup=kTRUE;
    }
    if(!fMultSelection->GetThisEventIsNotPileupInMultBins()){
      fHistMultSelPUCount->Fill(2.5);
      bIsPileup=kTRUE;
    }
    if(!fMultSelection->GetThisEventHasNoInconsistentVertices()){
      fHistMultSelPUCount->Fill(2.5);
      bIsPileup=kTRUE;
    }      
    if(!fMultSelection->GetThisEventPassesTrackletVsCluster()){
      fHistMultSelPUCount->Fill(2.5);
      bIsPileup=kTRUE;
    } 
    if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()){
      fHistMultSelPUCount->Fill(2.5);
      bIsPileup=kTRUE;
    }      
    if(!fMultSelection->GetThisEventHasGoodVertex2016()){
      fHistMultSelPUCount->Fill(2.5);
      bIsPileup=kTRUE;
    }      
    if(bIsPileup) fHistPileUpCount->Fill(9.5);
  }  
  //-----------------------------------------------------------------









  fHistTPCVsESDTrkBefore->Fill(multTPCAll,multEsd);   //A. Dobrin
  fHistTPCvsGlobalMultBefore->Fill(multGlobal,multTPC);
     

  Bool_t  bIsOutLier=kFALSE;
  
  if(multTPC < (-20.0+1.15*multGlobal) || multTPC > (200.+1.45*multGlobal)) { bIsOutLier = kTRUE;}
       
  fHistEventCount->Fill(stepCount); //6
  stepCount++;


  fHistTPConlyVsCL1Before->Fill(centrCL1,multTPCAll);
  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTPCAll);
  fHistGlobalVsV0MBefore->Fill(centrV0M, multGlobal);


  //if bSkipPileUpCut is kTRUE then don't apply PileUp removal.
  if(!bSkipPileUpCut && bIsOutLier) return; //outlier TPC vs Global

  fHistTPCvsGlobalMultAfter->Fill(multGlobal,multTPC);

  fHistEventCount->Fill(stepCount); //7
  stepCount++;


  if(!bSkipPileUpCut && bIsPileup) return;         //PileUp A. Dobrin

  fHistTPCVsESDTrkAfter->Fill(multTPCAll,multEsd);  

  fHistEventCount->Fill(stepCount); //8
  stepCount++;

  //cout<<"After PU cut multTPC = "<<multTPC<<" multGlobal = "<<multGlobal<<endl;





  /*
  Int_t icentBin = centrality;
  icentBin++;

  Float_t TPCmultLowLimit  =  hCentvsTPCmultCuts->GetBinContent(icentBin,1);
  Float_t TPCmultHighLimit =  hCentvsTPCmultCuts->GetBinContent(icentBin,1);

  TPCmultLowLimit  -=  5.0 * hCentvsTPCmultCuts->GetBinContent(icentBin,2); //mean - 5sigma
  TPCmultHighLimit +=  5.0 * hCentvsTPCmultCuts->GetBinContent(icentBin,2); //mean + 5sigma
  //std::cout<<" Cent = "<<centrality<<"\t icent = "<<icentBin<<" low = "<<TPCmultLowLimit<<"\t high = "<<TPCmultHighLimit<<std::endl;

  //centrality outlier
  if(!bSkipPileUpCut){  if(multTPC<TPCmultLowLimit || multTPC>TPCmultHighLimit) return;   }
  */





  fHistEventCount->Fill(stepCount); //9
  stepCount++;



  fHistTPConlyVsCL1After->Fill(centrCL1,multTPCAll);
  fHistTPConlyVsV0MAfter->Fill(centrV0M,multTPCAll);
  fHistGlobalVsV0MAfter->Fill(centrV0M, multGlobal);


  // MC corrected Refmult:
  fHistRawVsCorrMultFB->Fill(RefMultRawFB,RefMultCorrFB); // FB set by AddTask.. 


  Float_t EvtCent = centrality;




  if(gPsiN > 0){
   if(multEtaNeg<2 || multEtaPos<2)  return;  //Minimum 2 tracks in each eta
  }

  fHistEventCount->Fill(stepCount); //10
  stepCount++;

  //--------------------------------------------------------


















  //--------------- cent CL1 <= 90 cut ------------------
  Int_t icentV0Qn = centrCL1; // cent CL1 used for V0 calibration.
  icentV0Qn += 1;

  if(icentV0Qn>90)     return;



  fCentDistVsVz->Fill(pVtxZ,centrality);
 

  
  fHistEventCount->Fill(stepCount); //11
  stepCount++;




  //=============== Get the ZDC data ====================
  const Int_t fNumberOfZDCChannels = 5;

  if(gPsiN>0) {  // Fill ZDC info for gPsiN > 0

    AliAODZDC *aodZDC = fAOD->GetZDCData();

    if(!aodZDC) {
      cout<<"\n ********* Error: could not find ZDC data ************ \n "<<endl;
    }

    const Double_t *gZNATowerRawAOD = aodZDC->GetZNATowerEnergy();
    const Double_t *gZNCTowerRawAOD = aodZDC->GetZNCTowerEnergy();
    const Double_t *gZPATowerRawAOD = aodZDC->GetZPATowerEnergy();
    const Double_t *gZPCTowerRawAOD = aodZDC->GetZPCTowerEnergy();

    Double_t gZNATowerRaw[5] = {0.,0.,0.,0.,0.};
    Double_t gZNCTowerRaw[5] = {0.,0.,0.,0.,0.};
    Double_t gZPATowerRaw[5] = {0.,0.,0.,0.,0.};
    Double_t gZPCTowerRaw[5] = {0.,0.,0.,0.,0.};

    for(Int_t iChannel = 0; iChannel < fNumberOfZDCChannels; iChannel++) {
      gZNATowerRaw[iChannel] = gZNATowerRawAOD[iChannel];
      gZNCTowerRaw[iChannel] = gZNCTowerRawAOD[iChannel];
      gZPATowerRaw[iChannel] = gZPATowerRawAOD[iChannel];
      gZPCTowerRaw[iChannel] = gZPCTowerRawAOD[iChannel];

      fZPASignalPerChVsCent->Fill(centrCL1,iChannel,gZPATowerRaw[iChannel]);
      fZPCSignalPerChVsCent->Fill(centrCL1,iChannel,gZPCTowerRaw[iChannel]);
      fZNASignalPerChVsCent->Fill(centrCL1,iChannel,gZNATowerRaw[iChannel]);
      fZNCSignalPerChVsCent->Fill(centrCL1,iChannel,gZNCTowerRaw[iChannel]);

    }// ZDC channel loop

  }



  //-------- V0M info ---------------

  //do v0m recentering
  Double_t QxanCor = 0, QyanCor = 0;
  Double_t QxcnCor = 0, QycnCor = 0;
   
  Double_t Qxan3  = 0., Qyan3 = 0.;
  Double_t Qxcn3  = 0., Qycn3 = 0.;
  Double_t Qxan2  = 0., Qyan2 = 0.;
  Double_t Qxcn2  = 0., Qycn2 = 0.;

  Double_t phiV0;
  Float_t fMultv0 = 0;
  Float_t sumMa = 0;
  Float_t sumMc = 0;

  
  if(gPsiN>0) {  // Calculate EP only for gPsiN > 0

    const AliAODVZERO *fAODV0 = fAOD->GetVZEROData();

    for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA

      fMultv0 = fAODV0->GetMultiplicity(iV0);

      fV0MultChVsRun->Fill(iV0+0.5,centrCL1,fMultv0);

      if(fHCorrectV0M){
	fMultv0 = fMultv0 * fHCorrectV0M->GetBinContent(iV0+1);  // Gain Correction
	//cout<<"info: run = "<<runNumber<<" cent = "<<centrCL1<<"\t channel = "<<iV0<<" gain = "<<fHCorrectV0M->GetBinContent(iV0+1)<<endl;
      }
     

      //fV0MultChVsRun->Fill(iV0+0.5,runindex,fMultv0);

      phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);
 
      if(iV0 < 32){
	Qxcn2  += TMath::Cos(2*phiV0) * fMultv0;
	Qycn2  += TMath::Sin(2*phiV0) * fMultv0;
	Qxcn3  += TMath::Cos(3*phiV0) * fMultv0;
	Qycn3  += TMath::Sin(3*phiV0) * fMultv0;
	sumMc += fMultv0;
      }
      else if(iV0 >= 32){
	Qxan2  += TMath::Cos(2*phiV0) * fMultv0;
	Qyan2  += TMath::Sin(2*phiV0) * fMultv0;
	Qxan3  += TMath::Cos(3*phiV0) * fMultv0;
	Qyan3  += TMath::Sin(3*phiV0) * fMultv0;
	sumMa += fMultv0;
      }
    }//----- channel loop ----------

   
    if(gPsiN==3){
      QxanCor = Qxan3/sumMa;  //3rd order event plane
      QyanCor = Qyan3/sumMa;
      QxcnCor = Qxcn3/sumMc;
      QycnCor = Qycn3/sumMc;

      if(fHAvgerageQnV0C && fHAvgerageQnV0A && icentV0Qn < 91){
	QxanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,3);  //x = Cos
	QxcnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,3);  //x = Cos
	QyanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,4);  //y = Sin
	QycnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,4);  //y = Sin

	fV0AQ3xVsCentRun->Fill(centrCL1,QxanCor); 
	fV0AQ3yVsCentRun->Fill(centrCL1,QyanCor); 
	fV0CQ3xVsCentRun->Fill(centrCL1,QxcnCor);
	fV0CQ3yVsCentRun->Fill(centrCL1,QycnCor);
      }
      //printf("\n .... I am using my own V0 gain correction for Psi3...\n");
    }
    else{
      QxanCor = Qxan2/sumMa;  //2nd order Event plane 
      QyanCor = Qyan2/sumMa;
      QxcnCor = Qxcn2/sumMc;
      QycnCor = Qycn2/sumMc;

      if(fHAvgerageQnV0C && fHAvgerageQnV0A && icentV0Qn < 91){
	QxanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,1);  //x = Cos
	QxcnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,1);  //x = Cos
	QyanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,2);  //y = Sin
	QycnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,2);  //y = Sin  

	fV0AQ2xVsCentRun->Fill(centrCL1,QxanCor); 
	fV0AQ2yVsCentRun->Fill(centrCL1,QyanCor); 
	fV0CQ2xVsCentRun->Fill(centrCL1,QxcnCor);
	fV0CQ2yVsCentRun->Fill(centrCL1,QycnCor);
      }
      //printf("\n .... I am using my own V0 gain correction for Psi2...\n");
    }
   
    //------ For V0-Qn Recenter and Event plane: Uncorrectd ----------
    if(!fHAvgerageQnV0C && !fHAvgerageQnV0A){
      fV0CQ2xVsCentRun->Fill(centrCL1,Qxcn2/sumMc);
      fV0CQ2yVsCentRun->Fill(centrCL1,Qycn2/sumMc);
      fV0AQ2xVsCentRun->Fill(centrCL1,Qxan2/sumMa); 
      fV0AQ2yVsCentRun->Fill(centrCL1,Qyan2/sumMa); 

      fV0CQ3xVsCentRun->Fill(centrCL1,Qxcn3/sumMc);
      fV0CQ3yVsCentRun->Fill(centrCL1,Qycn3/sumMc);
      fV0AQ3xVsCentRun->Fill(centrCL1,Qxan3/sumMa); 
      fV0AQ3yVsCentRun->Fill(centrCL1,Qyan3/sumMa);
    }

    if(gPsiN>2){
      PsiNV0C = 1.0/gPsiN*( TMath::ATan2(QycnCor,QxcnCor) + TMath::Pi() );
      PsiNV0A = 1.0/gPsiN*( TMath::ATan2(QyanCor,QxanCor) + TMath::Pi() );
    }
    else{
      PsiNV0C = 1.0/gPsiN*TMath::ATan2(QycnCor,QxcnCor) ;
      if(PsiNV0C<0.) PsiNV0C += 2*TMath::Pi()/gPsiN;

      PsiNV0A = 1.0/gPsiN*TMath::ATan2(QyanCor,QxanCor) ;
      if(PsiNV0A<0.) PsiNV0A += 2*TMath::Pi()/gPsiN;
    }


    fHV0CEventPlaneVsCent->Fill(EvtCent,PsiNV0C);
    fHV0AEventPlaneVsCent->Fill(EvtCent,PsiNV0A);

  } //---------------------  if(gPsiN>0)  --------------------------







  fEventCount++;
 








 //---- Copies of TPC-Q vectors to remove track -by- track auto-correlation -----

  Double_t sumQxTPCneg;
  Double_t sumQyTPCneg;
  Double_t sumQxTPCpos;
  Double_t sumQyTPCpos;
  Double_t sumQxTPCneg2;
  Double_t sumQyTPCneg2;
  Double_t sumQxTPCpos2;
  Double_t sumQyTPCpos2;
  Double_t SumWgtNeg, SumWgtPos;
  Double_t SumWgtNeg2, SumWgtPos2;



  //--------- Track variable for PID/Charge studies ----------------
  Double_t   PDGmassPion   = 0.13957;
  Double_t   PDGmassKaon   = 0.49368;
  Double_t   PDGmassProton = 0.93827;

  PDGmassProton *= PDGmassProton;
  PDGmassPion   *= PDGmassPion;
  PDGmassKaon   *= PDGmassKaon;
 
  Double_t dEdx1, dEdx2, dPhi1, dPhi2, dPt1, dPt2;
  Double_t dEta1, dEta2, ptw1,  ptw2, deltaPhi;
  Double_t mom, w1NUA, w2NUA, WgtEP;

  
  Double_t nSigTOFpion,  nSigTPCpion;
  Double_t nSigTOFkaon,  nSigTPCkaon;
  Double_t nSigTOFproton,nSigTPCproton;
  //Double_t nSigTOFpion2,  nSigTPCpion2;
  //Double_t nSigTOFkaon2,  nSigTPCkaon2;
  //Double_t nSigTOFproton2,nSigTPCproton2;
 

  //Tof variables
  //Double_t length, tofTime, probMis, mass, beta;
  //Double_t c = TMath::C()*1.E-9;                //bright light m/ns 
  //Int_t    TOFmatch=0;



  Int_t    ptBin,gCharge1,gCharge2;

  Double_t dcaXY, dcaZ ;


  //----------- Set the desired Harmonic ------------
  Int_t n = gN;
  Int_t m = gM;
  Int_t p =n+m;
  //------------------------------------------------


  Int_t skipPairHBT = 0;





  Double_t Chi2Trk1=0, Chi2Trk2=0;

  Double_t ptwPion1,ptwKaon1,ptwProton1; 
  Double_t wNUAPion1,wNUAKaon1,wNUAProton1;

  //Double_t ptwPion2,ptwKaon2,ptwProton2;
  //Double_t wNUAPion2,wNUAKaon2,wNUAProton2;

  //Double_t    WgtEPPion = 1.0;
  //Double_t    WgtEPKaon = 1.0;
  //Double_t    WgtEPProton = 1.0;


  Bool_t isPion1 = kFALSE;
  Bool_t isKaon1  = kFALSE;
  Bool_t isProton1 = kFALSE;

  Bool_t isPion2 = kFALSE;
  Bool_t isKaon2  = kFALSE;
  Bool_t isProton2 = kFALSE;

  Int_t multPOI1st = 0;
  Int_t multPOI2nd = 0;


  //const Int_t maxTrack =  40000;




  //Calling fPIDResponse in nested loop is CPU expensive.
  //Store nSigma values in a array:


  /*
  for(int itrack = 0; itrack < ntracks; itrack++) {

    AliAODTrack *trackForPID=dynamic_cast<AliAODTrack*>(fVevent->GetTrack(itrack));
    if(!trackForPID) continue;

   // Array method:
    if(trackForPID){
      nSigPionTPC[itrack]   = fPIDResponse->NumberOfSigmasTPC(trackForPID, AliPID::kPion);
      nSigKaonTPC[itrack]   = fPIDResponse->NumberOfSigmasTPC(trackForPID, AliPID::kKaon);
      nSigProtonTPC[itrack] = fPIDResponse->NumberOfSigmasTPC(trackForPID, AliPID::kProton);
      nSigPionTOF[itrack]   = fPIDResponse->NumberOfSigmasTOF(trackForPID, AliPID::kPion);
      nSigKaonTOF[itrack]   = fPIDResponse->NumberOfSigmasTOF(trackForPID, AliPID::kKaon);
      nSigProtonTOF[itrack] = fPIDResponse->NumberOfSigmasTOF(trackForPID, AliPID::kProton);   
    }
    else{
      nSigPionTPC[itrack]   = -99; 
      nSigKaonTPC[itrack]   = -99; 
      nSigProtonTPC[itrack] = -99; 
      nSigPionTOF[itrack]   = -99; 
      nSigKaonTOF[itrack]   = -99; 
      nSigProtonTOF[itrack] = -99; 
    }  

    //if(itrack%10==0)
    //cout<< "nSig pi = " <<nSigPionTPC[itrack] << "nSig K = " << nSigKaonTPC[itrack]  << "nSig p = " << nSigProtonTOF[itrack] << endl;

    // Vector method:
    //if(trackForPID){
      //nSigPionTPC.push_back(fPIDResponse->NumberOfSigmasTPC(trackForPID, AliPID::kPion));
      //nSigKaonTPC.push_back(fPIDResponse->NumberOfSigmasTPC(trackForPID,  AliPID::kKaon));
      //nSigProtonTPC.push_back(fPIDResponse->NumberOfSigmasTPC(trackForPID, AliPID::kProton));

      //nSigPionTOF.push_back(fPIDResponse->NumberOfSigmasTOF(trackForPID, AliPID::kPion));
      //nSigKaonTOF.push_back(fPIDResponse->NumberOfSigmasTOF(trackForPID,  AliPID::kKaon));
      //nSigProtonTOF.push_back(fPIDResponse->NumberOfSigmasTOF(trackForPID, AliPID::kProton));
    //}
    //else{
      //nSigPionTPC.push_back(-100.);
      //nSigKaonTPC.push_back(-100.);
      //nSigProtonTPC.push_back(-100.);

      //nSigPionTOF.push_back(-100.);
      //nSigKaonTOF.push_back(-100.);
      //nSigProtonTOF.push_back(-100.);
    //}
  }*/
   





  //Correct with Centrality Wgts :   in second pass 
  Int_t iCentBinWgt = (Int_t) centrality;
  iCentBinWgt += 1;

  Double_t fWgtCent = 1.0;

  
  if(fHCentWeightForRun){
    fWgtCent = fHCentWeightForRun->GetBinContent(iCentBinWgt);
  }

  //Fill Centrality for run-by-run:  in first pass over data
  fCentDistAfter->Fill(centrality,fWgtCent); 
  





  //#########################################
  //          Track loop starting          //
  //#########################################




  for(Int_t itrack = 0; itrack < ntracks; itrack++) {

    AliAODTrack *track=dynamic_cast<AliAODTrack*>(fVevent->GetTrack(itrack));
    if(!track) continue;

    if(!track->TestFilterBit(fFilterBit)) continue;



    mom      = track->P();
    dPt1     = track->Pt();
    dPhi1    = track->Phi();
    dEta1    = track->Eta();
    gCharge1 = track->Charge();
    dEdx1    = track->GetDetPid()->GetTPCsignal();
    Chi2Trk1 = track->Chi2perNDF();

    //dcaXY  = track->DCA();
    //dcaZ   = track->ZAtDCA();
        
    /* Turn off QAs 
    //-------------- Check TOF status ------------------
    AliPIDResponse::EDetPidStatus status;
    status = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,track);
    TOFmatch = 0;
    if(status==AliPIDResponse::kDetPidOk){
      TOFmatch++;
    }
    probMis = fPIDResponse->GetTOFMismatchProbability(track);
    fHistTOFMatchCount->Fill(TOFmatch,probMis);

    mass = -9.9;
    beta = -0.9;

    if(TOFmatch>0 && probMis < 0.01) { 
      //This conditions are called when detector status is checked above :  
      //if((track->IsOn(AliAODTrack::kTOFin)) && (track->IsOn(AliAODTrack::kTIME)) &&  (track->IsOn(AliAODTrack::kTOFout))) {
      //if((track->IsOn(AliAODTrack::kITSin)) && (track->IsOn(AliAODTrack::kTOFpid))) { //Naghmeh used it        
	tofTime = track->GetTOFsignal();  // in pico seconds
	length  = track->GetIntegratedLength();   
	tofTime = tofTime*1.E-3; // ns
	if (tofTime <= 0) tofTime = 9999;            
	length = length*0.01; // in meters
	tofTime = tofTime*c;
	beta = length/tofTime;
	mass = mom*mom*(1./(beta*beta) - 1);             
    }//------------ TOF signal -------------------------


    //QA histograms  before applying track cuts:
    fHistEtaPtBefore->Fill(dEta1,dPt1);
    fHistTPCdEdxvsPBefore->Fill(mom*gCharge1,dEdx1);
    fHistTOFBetavsPBefore->Fill(mom*gCharge1,beta);
    fHistTOFMassvsPtBefore->Fill(dPt1*gCharge1,mass);  */




    //-------- Apply Default track cuts for analysis: ---------
    // if((dPt1 > fMaxPtCut) || (dPt1 < fMinPtCut) || (dEta1 > fMaxEtaCut) || (dEta1 < fMinEtaCut) || (dEdx1 < fdEdxMin) || (track->GetTPCNcls() < fTPCclustMin)  || (Chi2Trk1 < fTrkChi2Min) || (Chi2Trk1 > 4.0) || (track->DCA() > fDCAxyMax) || (track->ZAtDCA() > fDCAzMax)  || !(TMath::Abs(gCharge1)))  continue;
    //Rihan Dec05,2019: Added new cuts for new Analysis. !!! Cuts for CME run-2 is commented out above^
    if((dPt1 > fMaxPtCut) || (dPt1 < fMinPtCut) || (dEta1 > fMaxEtaCut) || (dEta1 < fMinEtaCut) || (dEdx1 < fdEdxMin) || (track->GetTPCNcls() < fTPCclustMin)  || (Chi2Trk1 < fTrkChi2Min) || (Chi2Trk1 > 4.0) || !(TMath::Abs(gCharge1)))  continue;
    
    multPOI1st++;

    dcaXY = track->DCA();
    dcaZ = track->ZAtDCA();

    /// Use strict DCA cut for sytematic check: ///
    if(fabs(dcaXY)>1.0) continue;




    
    //--------------------- PID signals 1st track-------------------------
    //Vector/Array both works same way:
    /*
    nSigTPCpion   = nSigPionTPC[itrack];
    nSigTPCkaon   = nSigKaonTPC[itrack];
    nSigTPCproton = nSigProtonTPC[itrack];
    nSigTOFpion   = nSigPionTOF[itrack];
    nSigTOFkaon   = nSigKaonTOF[itrack];
    nSigTOFproton = nSigProtonTOF[itrack];    
    */


    //if(itrack%100==0)
    //cout<<"Trk "<<itrack<<" pt1 = "<<dPt1<<"\tnSigPion = "<<nSigTPCpion<<"\tnSigKaon = "<<nSigTPCkaon <<"\tnSigprot = "<<nSigTPCproton<<endl;

    
    isPion1 = kFALSE;
    isKaon1  = kFALSE;
    isProton1 = kFALSE;

    //----> Rihan Dec05,2019:  Using TPC/TOF for PID
    nSigTPCpion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    nSigTPCkaon  = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    nSigTPCproton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    
    nSigTOFpion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    nSigTOFkaon  = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    nSigTOFproton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
      
    //----- Pion
    if(dPt1<=0.5 && TMath::Abs(nSigTPCpion)<=2.0){
      isPion1 = kTRUE;
    }/// Pion and Kaon max pT = 2.0 GeV
    //else if(dPt1>0.5 && dPt1<=2.0 && TMath::Abs(nSigTPCpion)<=3.0 && TMath::Abs(nSigTOFpion)<=2.0 ){
    //isPion1 = kTRUE;
    //}
    //----- Kaon
    if(dPt1<=0.45 && TMath::Abs(nSigTPCkaon)<=2.5){
      isKaon1 = kTRUE;
    }/// Pion and Kaon max pT = 2.0 GeV
    //else if(dPt1>0.45 && dPt1<=1.0 && TMath::Abs(nSigTPCkaon)<=3.0 && TMath::Abs(nSigTOFkaon)<=2.0){
    else if(dPt1>0.45 && dPt1<=1.0 && TMath::Abs(nSigTOFkaon)<=2.5){
      isKaon1 = kTRUE;
    }
    //----- Proton 
    if(dPt1<=0.8 && TMath::Abs(nSigTPCproton)<=3.0){
      isProton1 = kTRUE;
      if(gCharge1 > 0 && dPt1 < 0.4) isProton1 = kFALSE;   //only Protons below 0.4 GeV is garbage
    }
    else if(dPt1>0.8 && dPt1<=2.0 && TMath::Abs(nSigTOFproton)<=3.0){  
      isProton1 = kTRUE;
    }

    //<---------------------------------------------------------------------------------









    //=================  MC wgt and NUA wgt for PID =================
    ptwPion1 = 1.0;
    ptwKaon1  = 1.0;
    ptwProton1 = 1.0;
    wNUAPion1 = 1.0;
    wNUAKaon1  = 1.0;
    wNUAProton1 = 1.0;
    
    /*
    //------ get MC weight and NUA for Pion track1 --------------
    if(isPion1){
      if(fFB_Efficiency_Pion_Cent[cent10bin]){ // <-------------------- !!!! WARNING: use Pion Efficiency file when available.
	ptBin    = fFB_Efficiency_Pion_Cent[cent10bin]->FindBin(dPt1);
	ptwPion1 = 1.0/fFB_Efficiency_Pion_Cent[cent10bin]->GetBinContent(ptBin);
      }
      //else{ ptwPion1 = 1.0; }

      if(gCharge1>0){
	if(fHCorrectNUAposPion[cForNUA]){
	  iBinNUA   = fHCorrectNUAposPion[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAPion1 = fHCorrectNUAposPion[cForNUA]->GetBinContent(iBinNUA);
	}
	//else{ wNUAPion1 = 1.0; }
      }
      else{
	if(fHCorrectNUAnegPion[cForNUA]){
	  iBinNUA   = fHCorrectNUAnegPion[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAPion1 = fHCorrectNUAnegPion[cForNUA]->GetBinContent(iBinNUA);  
	}
	//else{ wNUAPion1 = 1.0; }
      }
    }

    //------ get MC weight and NUA for Kaon track1 --------------
    if(isKaon1){
      if(fFB_Efficiency_Kaon_Cent[cent10bin]){ // <-------------------- !!!! WARNING: use Kaon Efficiency file when available.
	ptBin    = fFB_Efficiency_Kaon_Cent[cent10bin]->FindBin(dPt1);
	ptwKaon1 = 1.0/fFB_Efficiency_Kaon_Cent[cent10bin]->GetBinContent(ptBin);
      }
      //else{ ptwKaon1 = 1.0; }

      if(gCharge1>0){
	if(fHCorrectNUAposKaon[cForNUA]){
	  iBinNUA   = fHCorrectNUAposKaon[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAKaon1 = fHCorrectNUAposKaon[cForNUA]->GetBinContent(iBinNUA);
	}
	//else{ wNUAKaon1 = 1.0; }
      }
      else{
	if(fHCorrectNUAnegKaon[cForNUA]){
	  iBinNUA   = fHCorrectNUAnegKaon[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAKaon1 = fHCorrectNUAnegKaon[cForNUA]->GetBinContent(iBinNUA);  
	}
	//else{ wNUAKaon1 = 1.0; }
      }
    }

    //------ get MC weight and NUA for Proton track1 --------------
    if(isProton1){
      if(gCharge1>0){
	if(fFB_Efficiency_Proton_Pos_Cent[cent10bin]){ // <-------------------- !!!! WARNING: use Proton Efficiency file when available.
	  ptBin    = fFB_Efficiency_Proton_Pos_Cent[cent10bin]->FindBin(dPt1);
	  ptwProton1 = 1.0/fFB_Efficiency_Proton_Pos_Cent[cent10bin]->GetBinContent(ptBin);
	}
	if(fHCorrectNUAposProton[cForNUA]){
	  iBinNUA   = fHCorrectNUAposProton[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAProton1 = fHCorrectNUAposProton[cForNUA]->GetBinContent(iBinNUA);
	}
      }
      else{
	if(fFB_Efficiency_Proton_Neg_Cent[cent10bin]){ // <-------------------- !!!! WARNING: use Proton Efficiency file when available.
	  ptBin    = fFB_Efficiency_Proton_Neg_Cent[cent10bin]->FindBin(dPt1);
	  ptwProton1 = 1.0/fFB_Efficiency_Proton_Neg_Cent[cent10bin]->GetBinContent(ptBin);
	}
	if(fHCorrectNUAnegProton[cForNUA]){
	  iBinNUA   = fHCorrectNUAnegProton[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAProton1 = fHCorrectNUAnegProton[cForNUA]->GetBinContent(iBinNUA);  
	}
      }
      //else{ ptwProton1 = 1.0; }
    }
    */
    //=========================== X ===============================


      
    //------ get MC weight and NUA for Charged  track 1--------------
    ptw1  = 1.0;
    w1NUA = 1.0;

    if(fFB_Efficiency_Cent[cent10bin]){
      ptBin = fFB_Efficiency_Cent[cent10bin]->FindBin(dPt1);
      ptw1  = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBin);
    }
    //else{ ptw1 = 1.0; }

    if(gCharge1>0){
      if(fHCorrectNUApos[cForNUA]){
	iBinNUA = fHCorrectNUApos[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	w1NUA = fHCorrectNUApos[cForNUA]->GetBinContent(iBinNUA);
      }
      //else{ w1NUA = 1.0; }
    }
    else{
      if(fHCorrectNUAneg[cForNUA]){
        iBinNUA = fHCorrectNUAneg[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
        w1NUA = fHCorrectNUAneg[cForNUA]->GetBinContent(iBinNUA);  
      }
      //else{ w1NUA = 1.0; }
    }

    if(w1NUA==0) w1NUA = 1.0;



   //============= charged hadron analysis: ============
    sumQxTPCneg = sumTPCQn2x[0];   // first copy to remove 1st track fron Q vector (AutoCorrelation)
    sumQyTPCneg = sumTPCQn2y[0];
    sumQxTPCpos = sumTPCQn2x[1];
    sumQyTPCpos = sumTPCQn2y[1];
    SumWgtNeg   = SumWEtaNeg;
    SumWgtPos   = SumWEtaPos;

    //-------- Remove phi1 from EP calculation ----------

    if(dEta1 < -0.05){
      sumQxTPCneg -= w1NUA*TMath::Cos(gPsiN*dPhi1);
      sumQyTPCneg -= w1NUA*TMath::Sin(gPsiN*dPhi1); // [0] = eta <-0.05
      SumWgtNeg   -= w1NUA;
    }
    else if(dEta1 > 0.05){
      sumQxTPCpos -= w1NUA*TMath::Cos(gPsiN*dPhi1);
      sumQyTPCpos -= w1NUA*TMath::Sin(gPsiN*dPhi1); // [1] = eta > 0.05
      SumWgtPos   -= w1NUA;
    }
    //-----------------------------------------------------





    //Alexandru Corrections:
    /*
    if (gCharge1 > 0){
        
      fpX1X3PosT1->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Cos(gPsiN*PsiNTPCC),ptw1*w1NUA);
      fpX1X3PosT2->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Cos(gPsiN*PsiNTPCC),ptw1*w1NUA);
      fpX1X3PosT3->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Sin(gPsiN*PsiNTPCC),ptw1*w1NUA);
      fpX1X3PosT4->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Sin(gPsiN*PsiNTPCC),ptw1*w1NUA);

      fpX1PosCosT1->Fill(EvtCent,TMath::Cos(dPhi1),ptw1*w1NUA);
      fpX1PosSinT1->Fill(EvtCent,TMath::Sin(dPhi1),ptw1*w1NUA);
	  
      fpX1X3PosT1EP2->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Cos(gPsiN*PsiNTPCA),ptw1*w1NUA);
      fpX1X3PosT2EP2->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Cos(gPsiN*PsiNTPCA),ptw1*w1NUA);
      fpX1X3PosT3EP2->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Sin(gPsiN*PsiNTPCA),ptw1*w1NUA);
      fpX1X3PosT4EP2->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Sin(gPsiN*PsiNTPCA),ptw1*w1NUA);
          
    }  
    else if (gCharge1 < 0){
        
      fpX1X3NegT1->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Cos(gPsiN*PsiNTPCC),ptw1*w1NUA);
      fpX1X3NegT2->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Cos(gPsiN*PsiNTPCC),ptw1*w1NUA);
      fpX1X3NegT3->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Sin(gPsiN*PsiNTPCC),ptw1*w1NUA);
      fpX1X3NegT4->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Sin(gPsiN*PsiNTPCC),ptw1*w1NUA);

      fpX1NegCosT1->Fill(EvtCent,TMath::Cos(dPhi1),ptw1*w1NUA);
      fpX1NegSinT1->Fill(EvtCent,TMath::Sin(dPhi1),ptw1*w1NUA);
      
      fpX1X3NegT1EP2->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Cos(gPsiN*PsiNTPCA),ptw1*w1NUA);
      fpX1X3NegT2EP2->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Cos(gPsiN*PsiNTPCA),ptw1*w1NUA);
      fpX1X3NegT3EP2->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Sin(gPsiN*PsiNTPCA),ptw1*w1NUA);
      fpX1X3NegT4EP2->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Sin(gPsiN*PsiNTPCA),ptw1*w1NUA);          

    }

    */





    multPOI2nd = 0;   

    //===========> 2nd track loop (nested) <=================
    
    for(Int_t jtrack = 0; jtrack < ntracks; jtrack++) {

      if(jtrack==itrack) continue;

      if(n==0) continue;
    
      AliAODTrack *track2=dynamic_cast<AliAODTrack*>(fVevent->GetTrack(jtrack));
      if(!track2) continue;
      
      if(!(track2->TestFilterBit(fFilterBit)))    continue; 

      dPt2    = track2->Pt();
      dPhi2   = track2->Phi();
      dEta2   = track2->Eta();
      gCharge2= track2->Charge();
      dEdx2   = track2->GetDetPid()->GetTPCsignal();
      Chi2Trk2 = track2->Chi2perNDF();

      if((dPt2 > fMaxPtCut) || (dPt2 < fMinPtCut) || (dEta2 > fMaxEtaCut) || (dEta2 < fMinEtaCut) || (dEdx2 < fdEdxMin) || (track2->GetTPCNcls() < fTPCclustMin)  || (Chi2Trk2 < fTrkChi2Min) || (Chi2Trk2 > 4.0) || (track2->DCA() > fDCAxyMax) || (track2->ZAtDCA() > fDCAzMax) || !(TMath::Abs(gCharge2)))
        continue;

      multPOI2nd++;



      /*
      if(track2->GetDetPid()->GetTPCsignal() < 10.0) continue;
      if(track2->GetTPCNcls()  <  70)                continue;
      if(track2->Chi2perNDF()  < 0.2)                continue;
      if(dPt2 < fMinPtCut  ||  dPt2 > fMaxPtCut)     continue;
      if(dEta2< fMinEtaCut || dEta2 > fMaxEtaCut)    continue;
      if(!TMath::Abs(gCharge2))                      continue;
      //if(!(track2->TestFilterBit(fFilterBit)))     continue; */ 
      //-----------------------------------------------------
		  
      //calling fPIDResponse is too costly for CPU


      isPion2 = kFALSE;
      isKaon2  = kFALSE;
      isProton2 = kFALSE;

      //--------------------- PID signals 2nd track-------------------------
      //Vector/Array both works same way:
      /*
      nSigTPCpion2   = nSigPionTPC[jtrack];
      nSigTPCkaon2   = nSigKaonTPC[jtrack];
      nSigTPCproton2 = nSigProtonTPC[jtrack];
      nSigTOFpion2   = nSigPionTOF[jtrack];
      nSigTOFkaon2   = nSigKaonTOF[jtrack];
      nSigTOFproton2 = nSigProtonTOF[jtrack];    

      //------> Pion
      if(dPt2<=0.6 && TMath::Abs(nSigTPCpion2)<=3.0){
	isPion2 = kTRUE;
      }
      else if(dPt2>0.6 && dPt2<=2.0 && TMath::Abs(nSigTPCpion2)<=3.0 && TMath::Abs(nSigTOFpion2)<=2.0 ){
	isPion2 = kTRUE;
      }
      //------> Kaon
      if(dPt2<=0.45 && TMath::Abs(nSigTPCkaon2)<=3.0){
	isKaon2 = kTRUE;
      }
      else if(dPt2>0.45 && dPt2<=2.0 && TMath::Abs(nSigTPCkaon2)<=3.0 && TMath::Abs(nSigTOFkaon2)<=2.0){
	isKaon2 = kTRUE;
      }
      //------> Proton 
      if(dPt2<=0.8 && TMath::Abs(nSigTPCproton2)<=3.0){
	isProton2 = kTRUE;
	if(dPt2<0.4) isProton2 = kFALSE;      //Proton below 0.4 GeV is garbage
      }
      else if(dPt2>0.8 && dPt2<=2.0 && TMath::Abs(nSigTPCproton2)<=3.0 && TMath::Abs(nSigTOFproton2)<=3.0){  //0.8 to 3.5 was used earlier.
	isProton2 = kTRUE;
      }
      */

      //----------------------------------------------------------------


 


      //=================  MC wgt and NUA wgt for PID =================
      /*
      WgtEPPion = 1.0;
      WgtEPKaon = 1.0;
      WgtEPProton = 1.0;

      ptwPion2 = 1.0;
      ptwKaon2  = 1.0;
      ptwProton2 = 1.0;
      wNUAPion2 = 1.0;
      wNUAKaon2  = 1.0;
      wNUAProton2 = 1.0;

      //------ get MC weight and NUA for Pion track2 --------------
      if(isPion2){
	if(fFB_Efficiency_Pion_Cent[cent10bin]){ 
	  ptBin    = fFB_Efficiency_Pion_Cent[cent10bin]->FindBin(dPt2);
	  ptwPion2 = 1.0/fFB_Efficiency_Pion_Cent[cent10bin]->GetBinContent(ptBin);
	}
	//else{ ptwPion2 = 1.0; }

	if(gCharge2>0){
	  if(fHCorrectNUAposPion[cForNUA]){
	    iBinNUA   = fHCorrectNUAposPion[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAPion2 = fHCorrectNUAposPion[cForNUA]->GetBinContent(iBinNUA);
	  }
	  //else{ wNUAPion2 = 1.0; }
	}
	else{
	  if(fHCorrectNUAnegPion[cForNUA]){
	    iBinNUA   = fHCorrectNUAnegPion[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAPion2 = fHCorrectNUAnegPion[cForNUA]->GetBinContent(iBinNUA);  
	  }
	  //else{ wNUAPion2 = 1.0; }
	}

	WgtEPPion = ptwPion1*ptwPion2*wNUAPion1*wNUAPion2;
      }

      //------ get MC weight and NUA for Kaon track2 --------------
      if(isKaon2){
	if(fFB_Efficiency_Kaon_Cent[cent10bin]){
	  ptBin    = fFB_Efficiency_Kaon_Cent[cent10bin]->FindBin(dPt2);
	  ptwKaon2 = 1.0/fFB_Efficiency_Kaon_Cent[cent10bin]->GetBinContent(ptBin);
	}
	//else{ ptwKaon2 = 1.0; }

	if(gCharge2>0){
	  if(fHCorrectNUAposKaon[cForNUA]){
	    iBinNUA   = fHCorrectNUAposKaon[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAKaon2 = fHCorrectNUAposKaon[cForNUA]->GetBinContent(iBinNUA);
	  }
	  //else{ wNUAKaon2 = 1.0; }
	}
	else{
	  if(fHCorrectNUAnegKaon[cForNUA]){
	    iBinNUA   = fHCorrectNUAnegKaon[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAKaon2 = fHCorrectNUAnegKaon[cForNUA]->GetBinContent(iBinNUA);  
	  }
	  //else{ wNUAKaon2 = 1.0; }
	}

	WgtEPKaon = ptwKaon1*ptwKaon2*wNUAKaon1*wNUAKaon2;
      }

      //------ get MC weight and NUA for Proton track2 --------------
      if(isProton2){
	if(gCharge2>0){
	  if(fFB_Efficiency_Proton_Pos_Cent[cent10bin]){ 
	    ptBin    = fFB_Efficiency_Proton_Pos_Cent[cent10bin]->FindBin(dPt2);
	    ptwProton2 = 1.0/fFB_Efficiency_Proton_Pos_Cent[cent10bin]->GetBinContent(ptBin);
	  }
	  if(fHCorrectNUAposProton[cForNUA]){
	    iBinNUA   = fHCorrectNUAposProton[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAProton2 = fHCorrectNUAposProton[cForNUA]->GetBinContent(iBinNUA);
	  }
	}
	else{
	  if(fFB_Efficiency_Proton_Neg_Cent[cent10bin]){ 
	    ptBin    = fFB_Efficiency_Proton_Neg_Cent[cent10bin]->FindBin(dPt2);
	    ptwProton2 = 1.0/fFB_Efficiency_Proton_Neg_Cent[cent10bin]->GetBinContent(ptBin);
	  }
	  if(fHCorrectNUAnegProton[cForNUA]){
	    iBinNUA   = fHCorrectNUAnegProton[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAProton2 = fHCorrectNUAnegProton[cForNUA]->GetBinContent(iBinNUA);  
	  }
	}
	//else{ ptwProton2 = 1.0; }
        WgtEPProton =  ptwProton1*ptwProton2*wNUAProton1*wNUAProton2;
      }
      */
      //========================== X ================================




 
      //------ get MC weight and NUA (Charged) for track 2--------------
      WgtEP = 1.0;
      ptw2  = 1.0;
      w2NUA = 1.0;

      if(fFB_Efficiency_Cent[cent10bin]){
	ptBin = fFB_Efficiency_Cent[cent10bin]->FindBin(dPt2);
	ptw2  = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBin);
      }
      //else{ ptw2 = 1.0; }

      if(gCharge2>0){
	if(fHCorrectNUApos[cForNUA]){
	  iBinNUA = fHCorrectNUApos[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	  w2NUA   = fHCorrectNUApos[cForNUA]->GetBinContent(iBinNUA);
	}
	//else{ w2NUA = 1.0; }
      }
      else{
	if(fHCorrectNUAneg[cForNUA]){
	  iBinNUA = fHCorrectNUAneg[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	  w2NUA   = fHCorrectNUAneg[cForNUA]->GetBinContent(iBinNUA);  
	}
	//else{ w2NUA = 1.0; }
      }

      if(w2NUA==0) w2NUA = 1.0;


      //---------- Remove track2 from EP calculation ---------
      sumQxTPCneg2 = sumQxTPCneg; //second copy to remove 2nd track from EP
      sumQyTPCneg2 = sumQyTPCneg;
      sumQxTPCpos2 = sumQxTPCpos;
      sumQyTPCpos2 = sumQyTPCpos;
      SumWgtNeg2   = SumWgtNeg;
      SumWgtPos2   = SumWgtPos;
      //------------------------------------------------------



      if(dEta2 < -0.05){
	sumQxTPCneg2 -= w2NUA*TMath::Cos(gPsiN*dPhi2);
	sumQyTPCneg2 -= w2NUA*TMath::Sin(gPsiN*dPhi2); // [0] = eta <-0.05
	SumWgtNeg2   -= w2NUA;
      }
      else if(dEta2 > 0.05){
	sumQxTPCpos2 -= w2NUA*TMath::Cos(gPsiN*dPhi2);
	sumQyTPCpos2 -= w2NUA*TMath::Sin(gPsiN*dPhi2); // [1] = eta > 0.05
	SumWgtPos2   -= w2NUA;
      }

      if(SumWgtNeg2>0 && SumWgtPos2>0){
        sumQyTPCneg2 = sumQyTPCneg2/SumWgtNeg2;
        sumQxTPCneg2 = sumQxTPCneg2/SumWgtNeg2;

        sumQyTPCpos2 = sumQyTPCpos2/SumWgtPos2;
        sumQxTPCpos2 = sumQxTPCpos2/SumWgtPos2;
      }

      // track by track EP:
      if(gPsiN>0){
	PsiNTPCC = (1.0/gPsiN)*( TMath::ATan2(sumQyTPCneg2,sumQxTPCneg2) ) ; // negetive eta
	if(PsiNTPCC<0.) PsiNTPCC += 2*TMath::Pi()/gPsiN;

	PsiNTPCA = (1.0/gPsiN)*( TMath::ATan2(sumQyTPCpos2,sumQxTPCpos2) ) ; // positive eta
	if(PsiNTPCA<0.) PsiNTPCA += 2*TMath::Pi()/gPsiN;
      }

      //fHTPCAEventPlaneVsCent->Fill(EvtCent,PsiNTPCA);
      //fHTPCCEventPlaneVsCent->Fill(EvtCent,PsiNTPCC);
      //-----------------------------------------------------------




      // combined weight for EP:
      WgtEP = ptw1*ptw2*w1NUA*w2NUA;
 


      deltaPhi = dPhi1 - dPhi2;   //for 2p correlator

      /////   Unlike Sign   ///////

      if(gCharge1!=gCharge2) {
	fHist_Corr3p_EP_Norm_PN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEP*fWgtCent); 
	fHist_Corr3p_EP_Norm_PN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEP*fWgtCent);
	fHist_Corr3p_EP_Norm_PN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	fHist_Corr3p_EP_Norm_PN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);
	/* 
	fHist_Corr3p_EP_Refm_PN[QAindex][0]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEP*fWgtCent); 
	fHist_Corr3p_EP_Refm_PN[QAindex][1]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEP*fWgtCent);
	fHist_Corr3p_EP_Refm_PN[QAindex][2]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	fHist_Corr3p_EP_Refm_PN[QAindex][3]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);		
	*/
	
	fpX1X2OppT1->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Cos(dPhi2),WgtEP);            
	fpX1X2OppT2->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Sin(dPhi2),WgtEP);
	fpX1X2OppT3->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Cos(dPhi2),WgtEP);
	fpX1X2OppT4->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Sin(dPhi2),WgtEP);
              
	//Differential 	//Temporarily Turn off - Rihan 14Oct19
	/*
	if(cIndex<6){
	  fHist_Corr3p_pTSum_EP_V0A_PN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5,  TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	  fHist_Corr3p_pTSum_EP_V0C_PN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5,  TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);
	  fHist_Corr3p_pTDiff_EP_V0A_PN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	  fHist_Corr3p_pTDiff_EP_V0C_PN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);
	  fHist_Corr3p_EtaDiff_EP_V0A_PN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	  fHist_Corr3p_EtaDiff_EP_V0C_PN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);
	}
	*/
	//-------------> PID CME ---------------
	//Pion:
	if(isPion1 && isPion2){
	  /*
	  fHist_Corr3p_Pion_EP_Norm_PN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPPion*fWgtCent);
	  fHist_Corr3p_Pion_EP_Norm_PN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPPion*fWgtCent);
	  fHist_Corr3p_Pion_EP_Norm_PN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPPion*fWgtCent);
	  fHist_Corr3p_Pion_EP_Norm_PN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPPion*fWgtCent);
	  
	  //Differential:
	  if(cIndex<6){
	    fHist_Corr3p_Pion_pTSum_EP_V0A_PN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_pTSum_EP_V0C_PN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_pTDiff_EP_V0A_PN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_pTDiff_EP_V0C_PN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_EtaDiff_EP_V0A_PN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_EtaDiff_EP_V0C_PN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPPion*fWgtCent);
	  }
	  */
	}
	//Kaon:
	if(isKaon1 && isKaon2){
	  /*
	  fHist_Corr3p_Kaon_EP_Norm_PN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPKaon*fWgtCent);
	  fHist_Corr3p_Kaon_EP_Norm_PN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPKaon*fWgtCent);
	  fHist_Corr3p_Kaon_EP_Norm_PN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPKaon*fWgtCent);
	  fHist_Corr3p_Kaon_EP_Norm_PN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPKaon*fWgtCent);

	  //Differential:
	  if(cIndex<6){
	    fHist_Corr3p_Kaon_pTSum_EP_V0A_PN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_pTSum_EP_V0C_PN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_pTDiff_EP_V0A_PN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_pTDiff_EP_V0C_PN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_EtaDiff_EP_V0A_PN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_EtaDiff_EP_V0C_PN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPKaon*fWgtCent);
	  }
	  */
	}
	//Proton:
	if(isProton1 && isProton2){
          //cout<<"#pair pt1 = "<<dPt1<<"\tpt2 = "<<dPt2<<"\tnSigp1 = "<<nSigTPCproton<<"\tnSigp2 = "<<nSigTPCproton2<<endl;
	  /*
	  fHist_Corr3p_Proton_EP_Norm_PN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPProton*fWgtCent);
	  fHist_Corr3p_Proton_EP_Norm_PN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPProton*fWgtCent);
	  fHist_Corr3p_Proton_EP_Norm_PN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPProton*fWgtCent);
	  fHist_Corr3p_Proton_EP_Norm_PN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPProton*fWgtCent);

	  //Differential:
	  if(cIndex<6){
	    fHist_Corr3p_Proton_pTSum_EP_V0A_PN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_pTSum_EP_V0C_PN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_pTDiff_EP_V0A_PN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_pTDiff_EP_V0C_PN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_EtaDiff_EP_V0A_PN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_EtaDiff_EP_V0C_PN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPProton*fWgtCent);
	  }
	  */
	}
	//------------------------------------
      }



      /////// Pos-Pos //////////



      else if(gCharge1>0 && gCharge2>0 && skipPairHBT==0) {
	fHist_Corr3p_EP_Norm_PP[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEP*fWgtCent);
	fHist_Corr3p_EP_Norm_PP[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEP*fWgtCent);
	fHist_Corr3p_EP_Norm_PP[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);	
	fHist_Corr3p_EP_Norm_PP[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);	
	/*
	fHist_Corr3p_EP_Refm_PP[QAindex][0]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEP*fWgtCent); 
	fHist_Corr3p_EP_Refm_PP[QAindex][1]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEP*fWgtCent);
	fHist_Corr3p_EP_Refm_PP[QAindex][2]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	fHist_Corr3p_EP_Refm_PP[QAindex][3]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);	
	*/
	fpX1X2PosT1->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Cos(dPhi2),WgtEP);            
	fpX1X2PosT2->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Sin(dPhi2),WgtEP);
	fpX1X2PosT3->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Cos(dPhi2),WgtEP);
	fpX1X2PosT4->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Sin(dPhi2),WgtEP);
	//Differential 	//Temporarily Turn off - Rihan 14Oct19
	
	/*
	if(cIndex<6){
	  fHist_Corr3p_pTSum_EP_V0A_PP[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5,  TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	  fHist_Corr3p_pTSum_EP_V0C_PP[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5,  TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);
	  fHist_Corr3p_pTDiff_EP_V0A_PP[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	  fHist_Corr3p_pTDiff_EP_V0C_PP[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);
	  fHist_Corr3p_EtaDiff_EP_V0A_PP[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	  fHist_Corr3p_EtaDiff_EP_V0C_PP[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);
	}
	*/
	//-------------> PID CME ---------------
	//Pion:
	if(isPion1 && isPion2){
	  /*
	  fHist_Corr3p_Pion_EP_Norm_PP[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPPion*fWgtCent);
	  fHist_Corr3p_Pion_EP_Norm_PP[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPPion*fWgtCent);
	  fHist_Corr3p_Pion_EP_Norm_PP[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPPion*fWgtCent);
	  fHist_Corr3p_Pion_EP_Norm_PP[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPPion*fWgtCent);
	
	  //Differential:
	  if(cIndex<6){
	    fHist_Corr3p_Pion_pTSum_EP_V0A_PP[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_pTSum_EP_V0C_PP[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_pTDiff_EP_V0A_PP[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_pTDiff_EP_V0C_PP[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_EtaDiff_EP_V0A_PP[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_EtaDiff_EP_V0C_PP[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPPion*fWgtCent);
	  }
	  */
	}
	//Kaon:
	if(isKaon1 && isKaon2){
	  /*
	  fHist_Corr3p_Kaon_EP_Norm_PP[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPKaon*fWgtCent);
	  fHist_Corr3p_Kaon_EP_Norm_PP[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPKaon*fWgtCent);
	  fHist_Corr3p_Kaon_EP_Norm_PP[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPKaon*fWgtCent);
	  fHist_Corr3p_Kaon_EP_Norm_PP[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPKaon*fWgtCent);
	
	  //Differential:
	  if(cIndex<6){
	    fHist_Corr3p_Kaon_pTSum_EP_V0A_PP[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_pTSum_EP_V0C_PP[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_pTDiff_EP_V0A_PP[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_pTDiff_EP_V0C_PP[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_EtaDiff_EP_V0A_PP[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_EtaDiff_EP_V0C_PP[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPKaon*fWgtCent);
	  }
	  */
	}
	//Proton:
	if(isProton1 && isProton2){
	  /*
	  fHist_Corr3p_Proton_EP_Norm_PP[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPProton*fWgtCent);
	  fHist_Corr3p_Proton_EP_Norm_PP[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPProton*fWgtCent);
	  fHist_Corr3p_Proton_EP_Norm_PP[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPProton*fWgtCent);
	  fHist_Corr3p_Proton_EP_Norm_PP[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPProton*fWgtCent);
	
	  //Differential:
	  if(cIndex<6){
	    fHist_Corr3p_Proton_pTSum_EP_V0A_PP[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_pTSum_EP_V0C_PP[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_pTDiff_EP_V0A_PP[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_pTDiff_EP_V0C_PP[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_EtaDiff_EP_V0A_PP[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_EtaDiff_EP_V0C_PP[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPProton*fWgtCent);
	  }
	  */
	}
	//------------------------------------
      }





      //////  Neg-Neg ///////


      else if(gCharge1<0 && gCharge2<0 && skipPairHBT==0){
	fHist_Corr3p_EP_Norm_NN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEP*fWgtCent);
	fHist_Corr3p_EP_Norm_NN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEP*fWgtCent);
	fHist_Corr3p_EP_Norm_NN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);	
	fHist_Corr3p_EP_Norm_NN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);	
	/*
	fHist_Corr3p_EP_Refm_NN[QAindex][0]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEP*fWgtCent); 
	fHist_Corr3p_EP_Refm_NN[QAindex][1]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEP*fWgtCent);
	fHist_Corr3p_EP_Refm_NN[QAindex][2]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	fHist_Corr3p_EP_Refm_NN[QAindex][3]->Fill(RefMultCorrFB, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);	
	*/
	fpX1X2NegT1->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Cos(dPhi2),WgtEP);            
	fpX1X2NegT2->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Sin(dPhi2),WgtEP);
	fpX1X2NegT3->Fill(EvtCent, TMath::Sin(dPhi1)*TMath::Cos(dPhi2),WgtEP);
	fpX1X2NegT4->Fill(EvtCent, TMath::Cos(dPhi1)*TMath::Sin(dPhi2),WgtEP);
	
	//Differential 	//Temporarily Turn off - Rihan 14Oct19
	/*
	if(cIndex<6){
	  fHist_Corr3p_pTSum_EP_V0A_NN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5,   TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	  fHist_Corr3p_pTSum_EP_V0C_NN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5,   TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);
	  fHist_Corr3p_pTDiff_EP_V0A_NN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	  fHist_Corr3p_pTDiff_EP_V0C_NN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);
	  fHist_Corr3p_EtaDiff_EP_V0A_NN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP*fWgtCent);
	  fHist_Corr3p_EtaDiff_EP_V0C_NN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP*fWgtCent);
	}
	*/
	//-------------> PID CME ---------------
	//Pion:
	if(isPion1 && isPion2){
	  /*
	  fHist_Corr3p_Pion_EP_Norm_NN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPPion*fWgtCent);
	  fHist_Corr3p_Pion_EP_Norm_NN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPPion*fWgtCent);
	  fHist_Corr3p_Pion_EP_Norm_NN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPPion*fWgtCent);
	  fHist_Corr3p_Pion_EP_Norm_NN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPPion*fWgtCent);

	  //Differential:
	  if(cIndex<6){
	    fHist_Corr3p_Pion_pTSum_EP_V0A_NN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_pTSum_EP_V0C_NN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_pTDiff_EP_V0A_NN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_pTDiff_EP_V0C_NN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_EtaDiff_EP_V0A_NN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPPion*fWgtCent);
	    fHist_Corr3p_Pion_EtaDiff_EP_V0C_NN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPPion*fWgtCent);
	  }
	  */
	}
	//Kaon:
	if(isKaon1 && isKaon2){
	  /*
	  fHist_Corr3p_Kaon_EP_Norm_NN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPKaon*fWgtCent);
	  fHist_Corr3p_Kaon_EP_Norm_NN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPKaon*fWgtCent);
	  fHist_Corr3p_Kaon_EP_Norm_NN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPKaon*fWgtCent);
	  fHist_Corr3p_Kaon_EP_Norm_NN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPKaon*fWgtCent);

	  //Differential:
	  if(cIndex<6){
	    fHist_Corr3p_Kaon_pTSum_EP_V0A_NN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_pTSum_EP_V0C_NN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_pTDiff_EP_V0A_NN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_pTDiff_EP_V0C_NN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_EtaDiff_EP_V0A_NN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPKaon*fWgtCent);
	    fHist_Corr3p_Kaon_EtaDiff_EP_V0C_NN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPKaon*fWgtCent);
	  }
          */
	}
	//Proton:
	if(isProton1 && isProton2){
	  /*
	  fHist_Corr3p_Proton_EP_Norm_NN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPProton*fWgtCent);
	  fHist_Corr3p_Proton_EP_Norm_NN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPProton*fWgtCent);
	  fHist_Corr3p_Proton_EP_Norm_NN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPProton*fWgtCent);
	  fHist_Corr3p_Proton_EP_Norm_NN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPProton*fWgtCent);

	  //Differential:
	  if(cIndex<6){
	    fHist_Corr3p_Proton_pTSum_EP_V0A_NN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_pTSum_EP_V0C_NN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_pTDiff_EP_V0A_NN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_pTDiff_EP_V0C_NN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_EtaDiff_EP_V0A_NN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A),WgtEPProton*fWgtCent);
	    fHist_Corr3p_Proton_EtaDiff_EP_V0C_NN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C),WgtEPProton*fWgtCent);
	  }
          */
	}
	//------------------------------------
      } 

      //cout<<", passes 3 WgtEP = "<<WgtEP<<endl; 
    }//-------- nested track loop ends ------------------










    //------------- Fill QA histograms ------------
    if(TMath::Abs(pVtxZ) < 8.0 && gCharge1 > 0){
      
      fQAEtaPhiAfterNUA->Fill(dPhi1,dEta1,w1NUA);  
      //Rihan Dec05,2019: Filling PID QA:
      
      if(isPion1){
        fQAEtaPhiAfterNUAPion->Fill(dPhi1,dEta1,wNUAPion1);
      }
      if(isKaon1){
        fQAEtaPhiAfterNUAKaon->Fill(dPhi1,dEta1,wNUAKaon1);
      }
      if(isProton1){
	fQAEtaPhiAfterNUAProton->Fill(dPhi1,dEta1,wNUAProton1);
      } 
    }








    //-------------- Fill NUA for PID tracks ----------------
    if(bFillNUAHistPID) {

      // Pion 
      if(isPion1){
        if(gCharge1>0){
	  fHist3DEtaPhiVz_Pos_Run[0][cForNUA]->Fill(pVtxZ,dPhi1,dEta1);
	}
	else if(gCharge1<0){
	  fHist3DEtaPhiVz_Neg_Run[0][cForNUA]->Fill(pVtxZ,dPhi1,dEta1);
	}
      }
     
      // Kaon 
      if(isKaon1){
        if(gCharge1>0){
	  fHist3DEtaPhiVz_Pos_Run[1][cForNUA]->Fill(pVtxZ,dPhi1,dEta1);
	}
	else if(gCharge1<0){
	  fHist3DEtaPhiVz_Neg_Run[1][cForNUA]->Fill(pVtxZ,dPhi1,dEta1);
	}
      }

      // proton
      if(isProton1){
        if(gCharge1>0){
	  fHist3DEtaPhiVz_Pos_Run[2][cForNUA]->Fill(pVtxZ,dPhi1,dEta1);
	}
	else if(gCharge1<0){
	  fHist3DEtaPhiVz_Neg_Run[2][cForNUA]->Fill(pVtxZ,dPhi1,dEta1);
	}
      }

      // Charged
      if(gCharge1>0){
	fHist3DEtaPhiVz_Pos_Run[3][cForNUA]->Fill(pVtxZ,dPhi1,dEta1);
      }
      else if(gCharge1<0){
	fHist3DEtaPhiVz_Neg_Run[3][cForNUA]->Fill(pVtxZ,dPhi1,dEta1);
      }
    
    }// Fill NUA for PID or not?





    //============ PID business starts here =============

    /*
    fHistEtaPtAfter->Fill(dEta1,dPt1);

    fHistTPCdEdxvsPAfter->Fill(mom*gCharge1,dEdx1);

    if(TOFmatch>0 && probMis < 0.01){
      fHistTOFBetavsPAfter->Fill(mom*gCharge1,beta);
    }

    //nSigTOFpion=fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    //nSigTOFkaon=fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    //nSigTOFproton=fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

    fHistTOFMatchCount->Fill(TOFmatch+2,nSigTOFpion);
    fHistTOFMatchCount->Fill(TOFmatch+4,nSigTOFkaon);
    fHistTOFMatchCount->Fill(TOFmatch+6,nSigTOFproton);

    if(!TOFmatch || probMis > 0.01){  // I dont want mismatched track in my signal distribution
      nSigTOFpion   = -99;
      nSigTOFkaon   = -99;
      nSigTOFproton = -99;
    }

    //nSigTPCpion   = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    //nSigTPCkaon   = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    //nSigTPCproton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);


    // Fill only 3D nSigma TOF and TPC:
    //0=pi, 1=K, 2=Proton
    fHistTPCTOFnSigmavsPtAfter[0]->Fill(dPt1,nSigTPCpion,nSigTOFpion);
    fHistTPCTOFnSigmavsPtAfter[1]->Fill(dPt1,nSigTPCkaon,nSigTOFkaon);
    fHistTPCTOFnSigmavsPtAfter[2]->Fill(dPt1,nSigTPCproton,nSigTOFproton);

    
    //  Fill nSigTOF for fixed nSigmaTPC Cut
    if(TMath::Abs(nSigTPCpion)<=fNSigmaCut){
      fHistTOFnSigmavsPtAfter[0]->Fill(dPt1*gCharge1,nSigTOFpion);
      fHistPtwithTPCNsigma[0]->Fill(dPt1*gCharge1);
      fHistTPCdEdxvsPtPIDAfter[0]->Fill(dPt1,dEdx1);
      if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        fHistPtwithTOFSignal[0]->Fill(dPt1*gCharge1);
      }
    }
    if(TMath::Abs(nSigTPCkaon)<=fNSigmaCut){
      fHistPtwithTPCNsigma[1]->Fill(dPt1*gCharge1);
      fHistTOFnSigmavsPtAfter[1]->Fill(dPt1*gCharge1,nSigTOFkaon);
      fHistTPCdEdxvsPtPIDAfter[1]->Fill(dPt1,dEdx1);
      if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        fHistPtwithTOFSignal[1]->Fill(dPt1*gCharge1);
      }
    }
    if(TMath::Abs(nSigTPCproton)<=fNSigmaCut){
      fHistPtwithTPCNsigma[2]->Fill(dPt1*gCharge1);
      fHistTOFnSigmavsPtAfter[2]->Fill(dPt1*gCharge1,nSigTOFproton);
      fHistTPCdEdxvsPtPIDAfter[2]->Fill(dPt1,dEdx1);
      if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        fHistPtwithTOFSignal[2]->Fill(dPt1*gCharge1);
      }
    }
   

    //  Fill dEdx Histograms
    if(TMath::Abs(nSigTPCpion)<=fNSigmaCut){
      fHistPtwithTPCNsigma[0]->Fill(dPt1*gCharge1);
      fHistTPCdEdxvsPtPIDAfter[0]->Fill(dPt1,dEdx1);
    }
    if(TMath::Abs(nSigTPCkaon)<=fNSigmaCut){
      fHistPtwithTPCNsigma[1]->Fill(dPt1*gCharge1);
      fHistTPCdEdxvsPtPIDAfter[1]->Fill(dPt1,dEdx1);
    }
    if(TMath::Abs(nSigTPCproton)<=fNSigmaCut){
      fHistPtwithTPCNsigma[2]->Fill(dPt1*gCharge1);
      fHistTPCdEdxvsPtPIDAfter[2]->Fill(dPt1,dEdx1);
    }



    //========> nSigmaTOF distribution for circular cut <===============
    //if(TMath::Sqrt(nSigTPCpion*nSigTPCpion+nSigTOFpion*nSigTOFpion)<=fNSigmaCut){
      //fHistTOFnSigmavsPtAfter[0]->Fill(dPt1*gCharge1,nSigTOFpion);
      //if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        //fHistPtwithTOFSignal[0]->Fill(dPt1*gCharge1);
      //}
    //}
    //if(TMath::Sqrt(nSigTPCkaon*nSigTPCkaon+nSigTOFkaon*nSigTOFkaon)<=fNSigmaCut){
      //fHistTOFnSigmavsPtAfter[1]->Fill(dPt1*gCharge1,nSigTOFkaon);
      //if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        //fHistPtwithTOFSignal[1]->Fill(dPt1*gCharge1);
      //}
    //}
    //if(TMath::Sqrt(nSigTPCproton*nSigTPCproton+nSigTOFproton*nSigTOFproton)<=fNSigmaCut){
      //fHistTOFnSigmavsPtAfter[2]->Fill(dPt1*gCharge1,nSigTOFproton);
      //if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        //fHistPtwithTOFSignal[2]->Fill(dPt1*gCharge1);
      //}
    //} 
 

    //------------------------------------------------------------
   
    // TOF matching is not needed from data as done by Davide using MC.
   
    if(!TOFmatch || probMis > 0.01 || beta>0.2) continue;

    // nSigmaTPC distribution for Fixed nSigmaTOF cut
    if(TMath::Abs(nSigTOFpion)<=fNSigmaCut){
      fHistTPCnSigmavsPtAfter[0]->Fill(dPt1*gCharge1,nSigTPCpion);
      fHistPtwithTOFmasscut[0]->Fill(dPt1*gCharge1);
    }
    if(TMath::Abs(nSigTOFkaon)<=fNSigmaCut){
      fHistTPCnSigmavsPtAfter[1]->Fill(dPt1*gCharge1,nSigTPCkaon);
      fHistPtwithTOFmasscut[1]->Fill(dPt1*gCharge1);
    }
    if(TMath::Abs(nSigTOFproton)<=fNSigmaCut){
      fHistTPCnSigmavsPtAfter[2]->Fill(dPt1*gCharge1,nSigTPCproton);
      fHistPtwithTOFmasscut[2]->Fill(dPt1*gCharge1);
    }    */



 } 
//===================== track loop ends ============================



  //Alex Corrections EP terms:
  fpX1EventT1EP1->Fill(EvtCent, TMath::Cos(gPsiN*PsiNTPCC),multPOI1st);
  fpX1EventT1EP2->Fill(EvtCent, TMath::Cos(gPsiN*PsiNTPCA),multPOI1st);
  fpY1EventT1EP1->Fill(EvtCent, TMath::Sin(gPsiN*PsiNTPCC),multPOI1st);
  fpY1EventT1EP2->Fill(EvtCent, TMath::Sin(gPsiN*PsiNTPCA),multPOI1st);  




  //---------- TPC Event Planes -------
  //Sub events //
  sumTPCQn2x[0] = sumTPCQn2x[0]/SumWEtaNeg;
  sumTPCQn2y[0] = sumTPCQn2y[0]/SumWEtaNeg;
  sumTPCQn2x[1] = sumTPCQn2x[1]/SumWEtaPos;
  sumTPCQn2y[1] = sumTPCQn2y[1]/SumWEtaPos;

  sumTPCQn2x[2] = sumTPCQn2x[2]/SumWEtaFull;
  sumTPCQn2y[2] = sumTPCQn2y[2]/SumWEtaFull;

  /*
  sumTPCQn3x[0] = sumTPCQn3x[0]/SumWEtaNeg;
  sumTPCQn3y[0] = sumTPCQn3y[0]/SumWEtaNeg;
  sumTPCQn4x[0] = sumTPCQn4x[0]/SumWEtaNeg;
  sumTPCQn4y[0] = sumTPCQn4y[0]/SumWEtaNeg;

  sumTPCQn3x[1] = sumTPCQn3x[1]/SumWEtaPos;
  sumTPCQn3y[1] = sumTPCQn3y[1]/SumWEtaPos;
  sumTPCQn4x[1] = sumTPCQn4x[1]/SumWEtaPos;
  sumTPCQn4y[1] = sumTPCQn4y[1]/SumWEtaPos;

  sumTPCQn3x[2] = sumTPCQn3x[2]/SumWEtaFull;
  sumTPCQn3y[2] = sumTPCQn3y[2]/SumWEtaFull;
  sumTPCQn4x[2] = sumTPCQn4x[2]/SumWEtaFull;
  sumTPCQn4y[2] = sumTPCQn4y[2]/SumWEtaFull;
  */

  if(gPsiN==3){
    //Enable periodicity PsiN directly:
    PsiNTPCC = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[0],sumTPCQn2x[0]) + TMath::Pi() ) ; // negetive eta
    PsiNTPCA = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[1],sumTPCQn2x[1]) + TMath::Pi() ) ; // positive eta
    PsiNTPCF = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[2],sumTPCQn2x[2]) + TMath::Pi() ) ; // FUll TPC eta 
  }
  else{
    PsiNTPCC = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[0],sumTPCQn2x[0]) ) ; // negetive eta
    if(PsiNTPCC<0.) PsiNTPCC += 2*TMath::Pi()/gPsiN;
    PsiNTPCA = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[1],sumTPCQn2x[1]) ) ; // positive eta
    if(PsiNTPCA<0.) PsiNTPCA += 2*TMath::Pi()/gPsiN;
    PsiNTPCF = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[2],sumTPCQn2x[2]) ) ; // FUll TPC eta 
    if(PsiNTPCF<0.) PsiNTPCF += 2*TMath::Pi()/gPsiN;
  }

  fHTPCAEventPlaneVsCent->Fill(EvtCent,PsiNTPCA);
  fHTPCCEventPlaneVsCent->Fill(EvtCent,PsiNTPCC);
  fHTPCEventPlaneVsCent->Fill(EvtCent, PsiNTPCF);





  //-------------- vs Centrality -----------------
  //V0A-V0C 
  fHist_Reso2n_EP_Norm_Det[QAindex][0]->Fill(EvtCent, TMath::Cos(gPsiN*(PsiNV0A-PsiNV0C)), fWgtCent);
  //V0A-TPC 
  fHist_Reso2n_EP_Norm_Det[QAindex][1]->Fill(EvtCent, TMath::Cos(gPsiN*(PsiNV0A-PsiNTPCF)), fWgtCent);
  //V0C-TPC 
  fHist_Reso2n_EP_Norm_Det[QAindex][2]->Fill(EvtCent, TMath::Cos(gPsiN*(PsiNV0C-PsiNTPCF)), fWgtCent);
  //TPCa -TPCc 
  fHist_Reso2n_EP_Norm_Det[QAindex][3]->Fill(EvtCent, TMath::Cos(gPsiN*(PsiNTPCA-PsiNTPCC)), fWgtCent);


  //---------  Store TPC-Qn for Recenter ---------
  fTPCAQ2xVsCentRun->Fill(EvtCent,sumTPCQn2x[0]); 
  fTPCAQ2yVsCentRun->Fill(EvtCent,sumTPCQn2y[0]);
  fTPCCQ2xVsCentRun->Fill(EvtCent,sumTPCQn2x[1]);
  fTPCCQ2yVsCentRun->Fill(EvtCent,sumTPCQn2y[1]);
  fTPCFQ2xVsCentRun->Fill(EvtCent,sumTPCQn2x[2]);
  fTPCFQ2yVsCentRun->Fill(EvtCent,sumTPCQn2y[2]);
  
  /*
  fTPCAQ3xVsCentRun->Fill(EvtCent,sumTPCQn3x[0]); 
  fTPCAQ3yVsCentRun->Fill(EvtCent,sumTPCQn3y[0]);
  fTPCCQ3xVsCentRun->Fill(EvtCent,sumTPCQn3x[1]);
  fTPCCQ3yVsCentRun->Fill(EvtCent,sumTPCQn3y[1]);

  fTPCAQ4xVsCentRun->Fill(EvtCent,sumTPCQn4x[0]); 
  fTPCAQ4yVsCentRun->Fill(EvtCent,sumTPCQn4y[0]);
  fTPCCQ4xVsCentRun->Fill(EvtCent,sumTPCQn4x[1]);
  fTPCCQ4yVsCentRun->Fill(EvtCent,sumTPCQn4y[1]);
  */











  PostData(1,fListHist);

  fHistEventCount->Fill(14.5); //15th bin is last one
  stepCount++;






  //if(fEventCount%10==0) 
  //cout<<"Ev = "<<fEventCount<<"\tMult = "<<multEtaFull<<"\tPOIs1st = "<<multPOI1st<<"\t POI2nd = "<< multPOI2nd <<endl;

  //watch.Stop();

}//================ UserExec ==============






















//////////////////// FUNCTIONS //////////////////////




double AliAnalysisTaskCMEV0PID::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
  // calculate sqrt of weighted distance to other vertex
  if (!v0 || !v1) {
    AliDebug(2,"\n\n ::GetWDist => One of vertices is not valid\n\n");
    return 0;
  }
  static TMatrixDSym vVb(3);
  double dist = -1;
  double dx = v0->GetX()-v1->GetX();
  double dy = v0->GetY()-v1->GetY();
  double dz = v0->GetZ()-v1->GetZ();
  double cov0[6],cov1[6];
  v0->GetCovarianceMatrix(cov0);
  v1->GetCovarianceMatrix(cov1);
  vVb(0,0) = cov0[0]+cov1[0];
  vVb(1,1) = cov0[2]+cov1[2];
  vVb(2,2) = cov0[5]+cov1[5];
  vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
  vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
  vVb.InvertFast();
  if (!vVb.IsValid()) {
    AliDebug(2,"Singular Matrix\n");
    return dist;
  }
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
  +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1;
}


Bool_t AliAnalysisTaskCMEV0PID::PileUpMultiVertex(const AliAODEvent* faod)
 {  // check for multi-vertexer pile-up
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if(!(nPlp=faod->GetNumberOfPileupVerticesTracks()))
  return kFALSE;

  vtPrm = faod->GetPrimaryVertex();
  if(vtPrm == faod->GetPrimaryVertexSPD())
  return kTRUE;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for(int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)faod->GetPileupVertexTracks(ipl);
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF()    > kMaxPlpChi2)    continue;
    //int bcPlp = vtPlp->GetBC();
    //if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2)
    // return kTRUE; // pile-up from other BC

    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist)        continue;

    return kTRUE; // pile-up: well separated vertices
    }
  return kFALSE;
}



void AliAnalysisTaskCMEV0PID::SetupMCcorrectionMap(){
  /*if(bApplyMCcorr){
    if(!gGrid){
    TGrid::Connect("alien://");
    }
    if(!mfileFBHijing){
    mfileFBHijing = TFile::Open(sMCfilePath,"READ");
    fListFBHijing = dynamic_cast<TList*>(mfileFBHijing->FindObjectAny("fMcEffiHij")); */

  if(fListFBHijing) {
    cout<<"\n =========> Info: Using MC efficiency correction Map <=========== "<<endl;
    for(int i=0;i<10;i++) {
      fFB_Efficiency_Cent[i] = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_%d",i));
    }
    //PID :
    /*
    for(int i=0;i<10;i++) {
      fFB_Efficiency_Pion_Cent[i]   = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_Pion_%d",i));
      fFB_Efficiency_Kaon_Cent[i]   = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_Kaon_%d",i));
      fFB_Efficiency_Proton_Pos_Cent[i] = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_ProtonPos_%d",i));
      fFB_Efficiency_Proton_Neg_Cent[i] = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_ProtonNeg_%d",i));
    }
    */
    //------- Fill the flags: --------------
    if(fFB_Efficiency_Cent[0] && fFB_Efficiency_Cent[4]){
      fHistTaskConfigParameters->SetBinContent(19,1);
    }
    /*
    if(fFB_Efficiency_Pion_Cent[0] && fFB_Efficiency_Pion_Cent[4]){
      fHistTaskConfigParameters->SetBinContent(20,1);
    }
    if(fFB_Efficiency_Kaon_Cent[0] && fFB_Efficiency_Kaon_Cent[4]){
      fHistTaskConfigParameters->SetBinContent(21,1);
    }
    if(fFB_Efficiency_Proton_Pos_Cent[0] && fFB_Efficiency_Proton_Neg_Cent[0]){
      fHistTaskConfigParameters->SetBinContent(22,1);
    }    */
  }
  else if(!fListFBHijing){
    std::cout<<"\n\n !!!!**** Warning: FB Efficiency File/List not found Use Weight = 1.0 *****\n\n"<<std::endl;
    //exit(1);
  }
}




Int_t AliAnalysisTaskCMEV0PID::GetCentralityScaled0to10(Double_t fCent){

 Int_t cIndex = 0;

 if(fCent<5.0) {
   cIndex  = 0; 
 }
 else if(fCent>=5.0 && fCent<10){
   cIndex  = 1;
 }
 else if(fCent>=10.0) {
   cIndex = abs(fCent/10.0)+1;
 }
 return cIndex;
}

void AliAnalysisTaskCMEV0PID::SetUpCentralityOutlierCut(){
  //std::cout<<" centrality outlier function called "<<std::endl;
 Float_t fMeanTPC[100] = {2902.95,2758.33,2642.78,2536.67,2435.37,2340.06,2248.44,2163.71,2080.49,2001.54,1925.86,1852.64,1781.97,1715.56,1650.53,1587.23,1527.51,1468.19,1412.73,1357.86,1305.35,1254.33,1205.57,1157.28,1111.53,1066.42,1023.15,981.594,940.795,901.766,863.651,826.183,790.53,756.358,722.654,690.513,659.443,628.807,599.748,571.664,544.446,518.042,492.369,468.072,444.694,422.487,400.104,379.129,359.147,339.62,320.817,302.788,285.791,269.015,253.688,238.671,224.039,209.932,196.915,184.647,172.76,161.381,150.395,140.288,131.033,121.58,113.112,104.938,97.3078,90.2178,83.5974,77.2645,70.7126,65.4424,60.1404,55.5644,50.8314,46.3761,43.024,38.625,35.3435,32.2304,29.4192,26.821,24.3303,21.9332,19.4215,16.7163,14.9414,13.1092,0.};

 Float_t fSigmaTPC[100] = {122.209,107.901,103.452,100.498,97.7403,94.7845,93.2543,90.0548,88.1106,85.7382,84.0812,82.2978,80.3817,78.6002,77.3448,75.5086,73.6842,71.9733,70.3447,69.1999,67.878,66.3511,65.0406,63.4866,62.4409,60.7899,59.1328,58.426,56.8618,55.8871,54.1031,53.4959,52.0482,51.0441,49.6218,48.7646,47.5166,46.5247,45.0727,44.4311,43.4531,42.0404,41.0238,40.1384,39.2588,38.2461,36.5951,36.0552,35.3727,33.7883,32.7167,32.4486,31.3709,30.3444,29.505,28.5139,27.4471,26.5359,25.9506,25.127,24.3797,23.2985,22.279,21.4698,20.781,20.8193,19.9509,18.8036,17.9145,16.961,16.7375,15.852,14.9324,14.7663,13.5969,13.4533,12.3067,12.7835,11.7283,10.6758,10.6676,10.6492,9.04614,8.89065,8.66093,8.50997,7.98812,6.91087,7.12045,7.29593,0.};

 for(int i=0;i<90;i++) {
   hCentvsTPCmultCuts->SetBinContent(i+1,1,fMeanTPC[i]);
   hCentvsTPCmultCuts->SetBinContent(i+1,2,fSigmaTPC[i]);
 }
}






void AliAnalysisTaskCMEV0PID::SetupEventAndTaskConfigInfo(){

  fHistTaskConfigParameters = new TH1F("fHistTaskConfigParameters","Task parameters",25,0,25);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(1,"FilterBit");
  fHistTaskConfigParameters->SetBinContent(1,fFilterBit);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(2,"n#sigmaTPC");
  fHistTaskConfigParameters->SetBinContent(2,fNSigmaCut);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(3,"MinPt");
  fHistTaskConfigParameters->SetBinContent(3,fMinPtCut);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(4,"MaxPt");
  fHistTaskConfigParameters->SetBinContent(4,fMaxPtCut);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(5,"MinEta");
  fHistTaskConfigParameters->SetBinContent(5,fMinEtaCut);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(6,"MaxEta");
  fHistTaskConfigParameters->SetBinContent(6,fMaxEtaCut);

  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(11,"CentralityMin");
  fHistTaskConfigParameters->SetBinContent(11,fCentralityPercentMin);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(12,"CentralityMax");
  fHistTaskConfigParameters->SetBinContent(12,fCentralityPercentMax);

  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(13,"VertexMin(cm)");
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(14,"VertexMax(cm)");


  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(15,"NUA Charge");
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(16,"NUA Pion");
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(17,"NUA Kaon");
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(18,"NUA Proton");

  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(19,"MC Charge");
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(20,"MC Pion");
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(21,"MC Kaon");
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(22,"MC Proton");

  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(23,"V0 Gain Eq.");
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(24,"V0 Qn Recenter");


  fListHist->Add(fHistTaskConfigParameters);



  fHistPileUpCount = new TH1F("fHistPileUpCount", "fHistPileUpCount", 15, 0., 15.);
  fHistPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fHistPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fHistPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultComb08");
  fHistPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fHistPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>5.0");
  fHistPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fHistPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  Int_t puConst = fPileUpConstParm;
  fHistPileUpCount->GetXaxis()->SetBinLabel(8,Form("multESDTPCDif>%d",puConst));
  fHistPileUpCount->GetXaxis()->SetBinLabel(9,Form("multGlobTPCDif>%d",puConst));
  fHistPileUpCount->GetXaxis()->SetBinLabel(10,"PileUpMultSelTask");
  fListHist->Add(fHistPileUpCount);


  fHistMultSelPUCount = new TH1F("fHistMultSelPileUpCount", "no of PU Event from MultSelTask", 10, 0., 10);
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(1,"PileUp");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(2,"PileUpMV");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(3,"PileUpMultBins");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(4,"InconsistentVtx");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(5,"TrackletVsCluster");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(6,"IncompleteDAQ");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(7,"NotGoodVertex2016");
  fListHist->Add(fHistMultSelPUCount);


  fHistEventCount = new TH1F("fHistEventCount","Event counts",15,0,15);
  fHistEventCount->GetXaxis()->SetBinLabel(1,"Called UserExec()");
  fHistEventCount->GetXaxis()->SetBinLabel(2,"Called Exec()");
  fHistEventCount->GetXaxis()->SetBinLabel(3,"AOD Exist");
  fHistEventCount->GetXaxis()->SetBinLabel(4,"Vz < 10");
  fHistEventCount->GetXaxis()->SetBinLabel(5,Form("%2.0f<Cent<%2.0f",fCentralityPercentMin,fCentralityPercentMax));
  fHistEventCount->GetXaxis()->SetBinLabel(6,"noAODtrack > 2 ");
  fHistEventCount->GetXaxis()->SetBinLabel(7,"TPC vs Global");
  fHistEventCount->GetXaxis()->SetBinLabel(8,"TPC128 vs ESD");
  fHistEventCount->GetXaxis()->SetBinLabel(9,"Cent vs TPC");
  fHistEventCount->GetXaxis()->SetBinLabel(10,"mult eta+/- > 2");
  fHistEventCount->GetXaxis()->SetBinLabel(11,"centCL1 < 90");
  fHistEventCount->GetXaxis()->SetBinLabel(15,"Survived Events");
  fListHist->Add(fHistEventCount);

  //fHistEventCount->Fill(1);

}





void AliAnalysisTaskCMEV0PID::GetV0MCorrectionHist(Int_t run)
{
  /*
  TList *fListAppliedV0Corr = new TList();
  fListAppliedV0Corr->SetName("fListAppliedV0Corr");
  fListAppliedV0Corr->SetOwner(kTRUE);
  fListHist->Add(fListAppliedV0Corr);*/  // to be Tested later.

  if(fListV0MCorr){
    fHCorrectV0M    = (TH1D *) fListV0MCorr->FindObject(Form("fHistV0Gain_Run%d",run));
    fHAvgerageQnV0A = (TH2D *) fListV0MCorr->FindObject(Form("fHistAvgQnV0A_Run%d",run));
    fHAvgerageQnV0C = (TH2D *) fListV0MCorr->FindObject(Form("fHistAvgQnV0C_Run%d",run));

    //Now the wgts for centrality is added inside VOM gain file. 
    fHCentWeightForRun = (TH1D *) fListV0MCorr->FindObject(Form("fHistCentWeight_Run%d",run));

    if(fHCentWeightForRun){
      cout<<"\n =========== Info:: Using Centrality wgts for run = "<<run<<"============"<<endl;
    }
    else{
      cout<<"\n =========== Info:: No Centrality wgt. correction.!! for run = "<<run<<"============"<<endl;
    }

    if(fHCorrectV0M){
      cout<<"\n =========== Info:: Setting up V0 gain correction for run = "<<run<<"============"<<endl;
      fHistTaskConfigParameters->SetBinContent(23,1);
    }
    else{
      cout<<"\n =========== Info:: No V0 gain correction..!!! for run = "<<run<<"============"<<endl;
    }
    if(fHAvgerageQnV0A && fHAvgerageQnV0C){
      cout<<"\n =========== Info:: Setting up V0 <Qn> recentering for run = "<<run<<"============"<<endl;
      fHistTaskConfigParameters->SetBinContent(24,1);
    }
    else{
      cout<<"\n =========== Info:: No V0 <Qn> recentering..!!! for run = "<<run<<"============"<<endl;
    }
  }
  else{
    cout<<"\n ======== !!!! Error:: List For V0-Gain and Qn not found for run "<<run<<"============"<<endl;
  }

}







void AliAnalysisTaskCMEV0PID::GetNUACorrectionHist(Int_t run)
{

  if(fListNUACorr){
    for(int i=0;i<5;i++){
      fHCorrectNUApos[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pos_Cent%d_Run%d",i,run)); 
      fHCorrectNUAneg[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Neg_Cent%d_Run%d",i,run));
    }
    if(fHCorrectNUApos[3] && fHCorrectNUAneg[3]){
      cout<<"\n=========== Info:: Setting up NUA corrections for run "<<run<<"============"<<endl;
      fHistTaskConfigParameters->SetBinContent(15,1);
    }
  }
  else {
    printf("\n ******** Warning: No NUA Correction for Charge Particle in run %d, Use Wgt = 1.0 ********* \n",run);
  }

  //===================   PID:  ==========================
  if(fListNUACorr){
    /*
    for(int i=0;i<5;i++){
      fHCorrectNUAposPion[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Pos_Cent%d_Run%d",i,run)); //
      fHCorrectNUAnegPion[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Neg_Cent%d_Run%d",i,run));
      if(fHCorrectNUAposPion[3] && fHCorrectNUAnegPion[3])   fHistTaskConfigParameters->SetBinContent(16,1);

      fHCorrectNUAposKaon[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Pos_Cent%d_Run%d",i,run)); //
      fHCorrectNUAnegKaon[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Neg_Cent%d_Run%d",i,run));
      if(fHCorrectNUAposKaon[3] && fHCorrectNUAnegKaon[3])   fHistTaskConfigParameters->SetBinContent(17,1);

      fHCorrectNUAposProton[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Pos_Cent%d_Run%d",i,run)); 
      fHCorrectNUAnegProton[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Neg_Cent%d_Run%d",i,run));
      if(fHCorrectNUAposProton[3] && fHCorrectNUAnegProton[3])   fHistTaskConfigParameters->SetBinContent(18,1);
    }
    if(fHCorrectNUAposPion[3] && fHCorrectNUAposKaon[3] && fHCorrectNUAposProton[3]) {
      cout<<"\n=========== Info:: Setting up --> PID NUA corrections for run =  "<<run<<"============"<<endl;
    }
    else{
      cout<<"\n=========== WARNING :: PID NUA corrections NOT found for run =  "<<run<<"============"<<endl;
    }*/
  }
  else {
    printf("\n ******** Error/Warning: No NUA Correction found for PID for run %d, Use PID Wgt = 1.0 ********* \n",run);
  }

}

