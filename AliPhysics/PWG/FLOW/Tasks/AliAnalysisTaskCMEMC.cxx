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

/* $Id: AliAnalysisTaskCMEMC.cxx  Rihan Haque 14/02/2018 $ */

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

//---AOD,ESD,MC event--
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
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
#include "AliAnalysisTaskCMEMC.h"

//using namespace std;

using std::cout;
using std::endl;
using std::vector;


ClassImp(AliAnalysisTaskCMEMC)


AliAnalysisTaskCMEMC::AliAnalysisTaskCMEMC(const char *name): AliAnalysisTaskSE(name),
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
  fHistEventCount(NULL)
{
  for(int i=0;i<4;i++){
    for(int j=0;j<10;j++){
      fHistPtwithoutPIDcut[i][j]=NULL;
      fHistPtwithTPCNsigma[i][j]=NULL;
      fHistPtwithTPCTOFNsigma[i][j]=NULL;
      fHistPtwithNoPIDCuts[i][j]=NULL;
    }
  }
  for(int i=0;i<3;i++){
    fHistPtwithTOFSignal[i]=NULL;
    fHistTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCnSigmavsPtAfter[i]=NULL;
    fHistTPCTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCdEdxvsPtPIDAfter[i]=NULL;
    V2IntProQC[i]=NULL;
    fResolution_TPCEP[i]=NULL;
  }
  for(int i=0;i<8;i++){
    fDnnPro[i] = NULL;
    fCMEPro[i] = NULL;
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
    //PID:
    for(int j=0;j<4;j++){
      fHist_Corr3p_Pion_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_NN[i][j]  =  NULL;
    }
  }
  //3p vs RefMult
  for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
      fHist_Corr3p_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr3p_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr3p_EP_Refm_NN[i][j]  =  NULL;
    }
    for(int j=0;j<4;j++) {
      fHist_Reso2n_EP_Refm_Det[i][j] =  NULL;
    }
    //PID:
    for(int j=0;j<4;j++){
      fHist_Corr3p_Pion_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Refm_NN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Refm_NN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Refm_NN[i][j]  =  NULL;
    }
  }


//2p vs Centrality:
 for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
      fHist_Corr2p_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr2p_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr2p_EP_Norm_NN[i][j]  =  NULL;
    }
    //PID:
    for(int j=0;j<4;j++){
      fHist_Corr2p_Pion_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr2p_Pion_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr2p_Pion_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Norm_NN[i][j]  =  NULL;
    }
  }

//2p vs Refm:
 for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
      fHist_Corr2p_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr2p_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr2p_EP_Refm_NN[i][j]  =  NULL;
    }
    //PID:
    for(int j=0;j<4;j++){
      fHist_Corr2p_Pion_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr2p_Pion_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr2p_Pion_EP_Refm_NN[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Refm_NN[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Refm_NN[i][j]  =  NULL;
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
  //Differential 3p PID:
  for(int i=0;i<2;i++){
    for(int j=0;j<6;j++){
      //Pion
      fHist_Corr3p_Pion_pTSum_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Pion_pTSum_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Pion_pTSum_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Pion_pTSum_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Pion_pTSum_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Pion_pTSum_EP_V0C_NN[i][j] = NULL;
    
      fHist_Corr3p_Pion_pTDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Pion_pTDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Pion_pTDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Pion_pTDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Pion_pTDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Pion_pTDiff_EP_V0C_NN[i][j] = NULL;

      fHist_Corr3p_Pion_EtaDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Pion_EtaDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Pion_EtaDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Pion_EtaDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Pion_EtaDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Pion_EtaDiff_EP_V0C_NN[i][j] = NULL;
      //Kaon
      fHist_Corr3p_Kaon_pTSum_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTSum_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_pTSum_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTSum_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTSum_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_pTSum_EP_V0C_NN[i][j] = NULL;
    
      fHist_Corr3p_Kaon_pTDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_pTDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_pTDiff_EP_V0C_NN[i][j] = NULL;

      fHist_Corr3p_Kaon_EtaDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_EtaDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_EtaDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Kaon_EtaDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_EtaDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_EtaDiff_EP_V0C_NN[i][j] = NULL;
      //Proton
      fHist_Corr3p_Proton_pTSum_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Proton_pTSum_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Proton_pTSum_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Proton_pTSum_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Proton_pTSum_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Proton_pTSum_EP_V0C_NN[i][j] = NULL;
    
      fHist_Corr3p_Proton_pTDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Proton_pTDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Proton_pTDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Proton_pTDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Proton_pTDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Proton_pTDiff_EP_V0C_NN[i][j] = NULL;

      fHist_Corr3p_Proton_EtaDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Proton_EtaDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Proton_EtaDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Proton_EtaDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Proton_EtaDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Proton_EtaDiff_EP_V0C_NN[i][j] = NULL;
    }
  }  



  //Differential 2p Charge:
  for(int i=0;i<2;i++){
    for(int j=0;j<6;j++){
      fHist_Corr2p_pTSum_EP_Harm1_PN[i][j] = NULL;
      fHist_Corr2p_pTSum_EP_Harm1_PP[i][j] = NULL;
      fHist_Corr2p_pTSum_EP_Harm1_NN[i][j] = NULL;
      fHist_Corr2p_pTSum_EP_Harm2_PN[i][j] = NULL;
      fHist_Corr2p_pTSum_EP_Harm2_PP[i][j] = NULL;
      fHist_Corr2p_pTSum_EP_Harm2_NN[i][j] = NULL;
    
      fHist_Corr2p_pTDiff_EP_Harm1_PN[i][j] = NULL;
      fHist_Corr2p_pTDiff_EP_Harm1_PP[i][j] = NULL;
      fHist_Corr2p_pTDiff_EP_Harm1_NN[i][j] = NULL;
      fHist_Corr2p_pTDiff_EP_Harm2_PN[i][j] = NULL;
      fHist_Corr2p_pTDiff_EP_Harm2_PP[i][j] = NULL;
      fHist_Corr2p_pTDiff_EP_Harm2_NN[i][j] = NULL;

      fHist_Corr2p_EtaDiff_EP_Harm1_PN[i][j] = NULL;
      fHist_Corr2p_EtaDiff_EP_Harm1_PP[i][j] = NULL;
      fHist_Corr2p_EtaDiff_EP_Harm1_NN[i][j] = NULL;
      fHist_Corr2p_EtaDiff_EP_Harm2_PN[i][j] = NULL;
      fHist_Corr2p_EtaDiff_EP_Harm2_PP[i][j] = NULL;
      fHist_Corr2p_EtaDiff_EP_Harm2_NN[i][j] = NULL;
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
    fFB_Efficiency_Pion_Cent[i] = NULL;
    fFB_Efficiency_Kaon_Cent[i] = NULL;
    fFB_Efficiency_Proton_Cent[i] = NULL;
  }


  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}
//______________________________empty constructor_______________________
AliAnalysisTaskCMEMC::AliAnalysisTaskCMEMC():
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
  fHistEventCount(NULL)
{
  for(int i=0;i<4;i++){
   for(int j=0;j<10;j++){
      fHistPtwithoutPIDcut[i][j]=NULL;
      fHistPtwithTPCNsigma[i][j]=NULL;
      fHistPtwithTPCTOFNsigma[i][j]=NULL;
      fHistPtwithNoPIDCuts[i][j]=NULL;
    }
  }
  for(int i=0;i<3;i++){
    fHistPtwithTOFSignal[i]=NULL;
    fHistTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCnSigmavsPtAfter[i]=NULL;
    fHistTPCTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCdEdxvsPtPIDAfter[i]=NULL;
    V2IntProQC[i]=NULL;
    fResolution_TPCEP[i]=NULL;
  }
  for(int i=0;i<8;i++){
    fDnnPro[i] = NULL;
    fCMEPro[i] = NULL;
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
    //PID:
    for(int j=0;j<4;j++){
      fHist_Corr3p_Pion_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_NN[i][j]  =  NULL;
    }
  }
  //3p vs RefMult
  for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
      fHist_Corr3p_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr3p_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr3p_EP_Refm_NN[i][j]  =  NULL;
    }
    for(int j=0;j<4;j++) {
      fHist_Reso2n_EP_Refm_Det[i][j] =  NULL;
    }
    //PID:
    for(int j=0;j<4;j++){
      fHist_Corr3p_Pion_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Refm_NN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Refm_NN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Refm_NN[i][j]  =  NULL;
    }
  }



//2p vs Centrality:
 for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
      fHist_Corr2p_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr2p_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr2p_EP_Norm_NN[i][j]  =  NULL;
    }
    //PID:
    for(int j=0;j<4;j++){
      fHist_Corr2p_Pion_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr2p_Pion_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr2p_Pion_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Norm_NN[i][j]  =  NULL;
    }
  }
//2p vs Refm:
 for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
      fHist_Corr2p_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr2p_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr2p_EP_Refm_NN[i][j]  =  NULL;
    }
    //PID:
    for(int j=0;j<4;j++){
      fHist_Corr2p_Pion_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr2p_Pion_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr2p_Pion_EP_Refm_NN[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr2p_Kaon_EP_Refm_NN[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr2p_Proton_EP_Refm_NN[i][j]  =  NULL;
    }
  }















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
  //Differential PID:
  for(int i=0;i<2;i++){
    for(int j=0;j<6;j++){
      //Pion
      fHist_Corr3p_Pion_pTSum_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Pion_pTSum_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Pion_pTSum_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Pion_pTSum_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Pion_pTSum_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Pion_pTSum_EP_V0C_NN[i][j] = NULL;
    
      fHist_Corr3p_Pion_pTDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Pion_pTDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Pion_pTDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Pion_pTDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Pion_pTDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Pion_pTDiff_EP_V0C_NN[i][j] = NULL;

      fHist_Corr3p_Pion_EtaDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Pion_EtaDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Pion_EtaDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Pion_EtaDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Pion_EtaDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Pion_EtaDiff_EP_V0C_NN[i][j] = NULL;
      //Kaon
      fHist_Corr3p_Kaon_pTSum_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTSum_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_pTSum_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTSum_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTSum_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_pTSum_EP_V0C_NN[i][j] = NULL;
    
      fHist_Corr3p_Kaon_pTDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_pTDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_pTDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_pTDiff_EP_V0C_NN[i][j] = NULL;

      fHist_Corr3p_Kaon_EtaDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_EtaDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_EtaDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Kaon_EtaDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Kaon_EtaDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Kaon_EtaDiff_EP_V0C_NN[i][j] = NULL;
      //Proton
      fHist_Corr3p_Proton_pTSum_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Proton_pTSum_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Proton_pTSum_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Proton_pTSum_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Proton_pTSum_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Proton_pTSum_EP_V0C_NN[i][j] = NULL;
    
      fHist_Corr3p_Proton_pTDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Proton_pTDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Proton_pTDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Proton_pTDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Proton_pTDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Proton_pTDiff_EP_V0C_NN[i][j] = NULL;

      fHist_Corr3p_Proton_EtaDiff_EP_V0A_PN[i][j] = NULL;
      fHist_Corr3p_Proton_EtaDiff_EP_V0A_PP[i][j] = NULL;
      fHist_Corr3p_Proton_EtaDiff_EP_V0A_NN[i][j] = NULL;
      fHist_Corr3p_Proton_EtaDiff_EP_V0C_PN[i][j] = NULL;
      fHist_Corr3p_Proton_EtaDiff_EP_V0C_PP[i][j] = NULL;
      fHist_Corr3p_Proton_EtaDiff_EP_V0C_NN[i][j] = NULL;
    }
  }



  //Differential 2p Charge:
  for(int i=0;i<2;i++){
    for(int j=0;j<6;j++){
      fHist_Corr2p_pTSum_EP_Harm1_PN[i][j] = NULL;
      fHist_Corr2p_pTSum_EP_Harm1_PP[i][j] = NULL;
      fHist_Corr2p_pTSum_EP_Harm1_NN[i][j] = NULL;
      fHist_Corr2p_pTSum_EP_Harm2_PN[i][j] = NULL;
      fHist_Corr2p_pTSum_EP_Harm2_PP[i][j] = NULL;
      fHist_Corr2p_pTSum_EP_Harm2_NN[i][j] = NULL;
    
      fHist_Corr2p_pTDiff_EP_Harm1_PN[i][j] = NULL;
      fHist_Corr2p_pTDiff_EP_Harm1_PP[i][j] = NULL;
      fHist_Corr2p_pTDiff_EP_Harm1_NN[i][j] = NULL;
      fHist_Corr2p_pTDiff_EP_Harm2_PN[i][j] = NULL;
      fHist_Corr2p_pTDiff_EP_Harm2_PP[i][j] = NULL;
      fHist_Corr2p_pTDiff_EP_Harm2_NN[i][j] = NULL;

      fHist_Corr2p_EtaDiff_EP_Harm1_PN[i][j] = NULL;
      fHist_Corr2p_EtaDiff_EP_Harm1_PP[i][j] = NULL;
      fHist_Corr2p_EtaDiff_EP_Harm1_NN[i][j] = NULL;
      fHist_Corr2p_EtaDiff_EP_Harm2_PN[i][j] = NULL;
      fHist_Corr2p_EtaDiff_EP_Harm2_PP[i][j] = NULL;
      fHist_Corr2p_EtaDiff_EP_Harm2_NN[i][j] = NULL;
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
    fFB_Efficiency_Pion_Cent[i] = NULL;
    fFB_Efficiency_Kaon_Cent[i] = NULL;
    fFB_Efficiency_Proton_Cent[i] = NULL;
  }
}

//___________________________ destructor ___________________________
AliAnalysisTaskCMEMC::~AliAnalysisTaskCMEMC()
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



void AliAnalysisTaskCMEMC::UserCreateOutputObjects()
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
  */

  for(int i=0;i<4;i++){
    for(int j=0;j<10;j++){
      fHistPtwithoutPIDcut[i][j]  = new TH1F(Form("fHistPtwithoutPIDcut_Species%d_Cent%d",i,j), Form("Particle%d;p_{T}(GeV/c))",i),200,-5,5); //i
      fListHist->Add(fHistPtwithoutPIDcut[i][j]);
      fHistPtwithTPCNsigma[i][j]   = new TH1F(Form("fHistPtwithTPCNsigma_Species%d_Cent%d",i,j), Form("Particle%d;p_{T}(GeV/c))",i),200,-5,5); //i
      fListHist->Add(fHistPtwithTPCNsigma[i][j]);
      fHistPtwithTPCTOFNsigma[i][j] = new TH1F(Form("fHistPtwithTPCTOFNsigma_Species%d_Cent%d",i,j),Form("Particle%d;p_{T}(GeV/c))",i),200,-5,5);
      fListHist->Add(fHistPtwithTPCTOFNsigma[i][j]);
      fHistPtwithNoPIDCuts[i][j] = new TH1F(Form("fHistPtwithNoPIDCuts_Species%d_Cent%d",i,j),Form("Particle%d;p_{T}(GeV/c))",i),200,-5,5);
      fListHist->Add(fHistPtwithNoPIDCuts[i][j]);
    }
  }


  for(int i=0;i<3;i++){
    /*
    //These are not filled:
    fHistPtwithTOFSignal[i]  = new TH1F(Form("fHistPtwithTOFSignal_%d", i),Form("%d;p_{T}(GeV/c))",i),200,-10,10);
    fListHist->Add(fHistPtwithTOFSignal[i]);
    fHistTOFnSigmavsPtAfter[i] = new TH2F(Form("fHistTOFnSigmavsPtAfter_%d",i),Form("%d;p_{T}(GeV/c);n#sigma_{TOF}",i),200,-10,10,200,-10.0,10.0);
    fListHist->Add(fHistTOFnSigmavsPtAfter[i]);
    fHistTPCnSigmavsPtAfter[i] = new TH2F(Form("fHistTPCnSigmavsPtAfter_%d",i),Form("%d;p_{T}(GeV/c);n#sigma_{TPC}",i),200,-10,10,200,-10.0,10.0);
    fListHist->Add(fHistTPCnSigmavsPtAfter[i]); */

    //These are filled:
    fHistTPCTOFnSigmavsPtAfter[i] = new TH3F(Form("fHistTPCTOFnSigmavsPtAfter_%d",i),Form("%d; p_{T}(GeV/c); n#sigma_{TPC}; n#sigma_{TOF}",i),100,0,5,200,-10,10,200,-10,10);
    fListHist->Add(fHistTPCTOFnSigmavsPtAfter[i]);

    fHistTPCdEdxvsPtPIDAfter[i] = new TH2F(Form("fHistTPCdEdxvsPtAfter_%d",i),"AfterCut; p_{T} (GeV/c); dEdx (arb)",200,0,10,250,0,250);
    fListHist->Add(fHistTPCdEdxvsPtPIDAfter[i]);
  }// PID histograms done
 




  Double_t centRange[11]   = {0,5,10,20,30,40,50,60,70,80,90};
  //const char *gDetForEP[4] = {"V0A","V0C","TPC-A","TPC-C"};


  //--------------------- v2, delta and gamma ------------------
  for(int i=0;i<3;i++){
    V2IntProQC[i] = new TProfile(Form("V%dIntProQC",i+2),Form("V%dIntProQC",i+2),10,centRange,"s");
    V2IntProQC[i]->Sumw2();
    fListHist->Add(V2IntProQC[i]);
  }
 
  for(Int_t i=0; i<8; i++) {
    fDnnPro[i] = new TProfile(Form("fDnnPro_%d",i),Form("fDnnPro_%d",i),10,centRange,"s");
    fDnnPro[i]->Sumw2();
    fListHist->Add(fDnnPro[i]);
  }
  for(Int_t i=0; i<8; i++) {
    fCMEPro[i] = new TProfile(Form("CMEPro_%d",i),Form("CMEPro_%d",i),10,centRange,"");
    //fCMEPro[i]->Sumw2();
    fListHist->Add(fCMEPro[i]);
  }

  for(int i=0;i<3;i++){
    fResolution_TPCEP[i] = new TProfile(Form("fResolutionTPCPsi%d",i+2),Form("Resolution #Psi_{%d}",i+2),10,centRange,"");
    //V2IntProQC[i]->Sumw2();
    fListHist->Add(fResolution_TPCEP[i]);
  }
 
  //------------------------------------------------------------






  fCentDistBefore = new TH1F("fMultDistMC","Mult Dist; Events ",2000,0,4000);
  fListHist->Add(fCentDistBefore);

  //fCentDistBefore = new TH1F("fCentDistBefore","no Cut; Cent (%); Events ",100,0,100);
  //fListHist->Add(fCentDistBefore);

  fCentDistAfter = new TH1F("fCentDistAfter","with Cut; Cent (%); Events ",100,0,100);
  fListHist->Add(fCentDistAfter);


















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

  fHTPCAEventPlaneVsCent = new TH2F("fPsi2TPCVsCent",Form("Psi %d from TPC",2),10,centRange,50,-0.0,truncPi);
  fListNUACalib->Add(fHTPCAEventPlaneVsCent);

  fHTPCCEventPlaneVsCent = new TH2F("fPsi3TPCVsCent",Form("Psi %d from TPC",3),10,centRange,50,-0.0,truncPi);
  fListNUACalib->Add(fHTPCCEventPlaneVsCent);

  fHTPCEventPlaneVsCent = new TH2F("fPsi4TPCVsCent",Form("Psi %d from TPC",4),10,centRange,50,-0.0,truncPi);
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






  //----------------- V0 Calibration hist: ---------------------
  fV0MultChVsRun = new TH2F("fV0MultChVsRun","1-32 V0C, 33-64 V0A",64,0,64,20,0,100);
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
  Char_t  name[100];
  Char_t title[100];

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

}

















//______________________________________________________________________
void AliAnalysisTaskCMEMC::UserExec(Option_t*) {
  //debug only

  //cout<<"\n Info:UserExec() called ..!!!\n";


  //watch.Start(kTRUE);  
  //if(fEventCount==501)  return;


  Float_t stepCount = 0.5;

  fHistEventCount->Fill(stepCount); //1
  stepCount++;

  //fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  //fESD = dynamic_cast<AliESDEvent*>(InputEvent());   
  // if(!(fESD || fAOD)){ printf("ERROR: fESD & fAOD not available\n"); return;  }
  //fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  // if (!fVevent) { printf("ERROR: fVevent not available\n");  return;  }
   
  fHistEventCount->Fill(stepCount); //2
  stepCount++;


  //--------- Load MC Event ---------
  AliVEvent *eventMain = dynamic_cast<AliVEvent*>(MCEvent());
  AliMCEvent* mcEvent = dynamic_cast<AliMCEvent*>(eventMain);

  if (!mcEvent) {
    printf("ERROR: Could not retrieve MC event");
    return;
  }

  //return;



  //--------- Check if I have PID response object --------
  // if(!fPIDResponse){
  //    printf("\n\n...... PIDResponse object not found..... \n\n"); return;
  // }


  //-------------- Vtx cuts ---------------
  // const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
   Double_t pVtxZ = 0.0;
  // pVtxZ = pVtx->GetZ();
    
//if(TMath::Abs(pVtxZ)>10.) return;
  //User defined cut:
  // if(pVtxZ<fMinVzCut || pVtxZ>fMaxVzCut ) return;


  fHistEventCount->Fill(stepCount); //3
  stepCount++;



  Float_t centrality = -99.0;
  Float_t centrV0M   = -99.0;
  Float_t centrCL1   = -99.0;

  //AliCentrality* Alicentr = mcEvent->GetCentrality();   // For Run1, 2010 data
  //centrality  = Alicentr->GetCentralityPercentile("CL1");


  /*
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
    else{
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
 
*/

  //cout<<" vtx = "<<pVtxZ<<"\t cent = "<<centrality<<endl;


  //----------------- Cut on  Centrality -----------------------
  // if(centrality<fCentralityPercentMin || centrality>fCentralityPercentMax){ 
  //   return;
  // }




  fHistEventCount->Fill(stepCount); //4
  stepCount++;

  //fCentDistBefore->Fill(centrality);
  

  



  //Int_t ntracks=fAOD->GetNumberOfTracks();
  //if(ntracks<2) return;              // Check this cut....!!!
   
  fHistEventCount->Fill(stepCount); //5
  stepCount++;



  centrality = 4.0;




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

  cIndex = cent10bin; //0,1,2,

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
  Double_t fMagField = -1.0; //fAOD->GetMagneticField();

  const Int_t QAindex = (fMagField > 0) ? 1 : 0;
  //---------------------------------






  //Load NUA and V0M correction map run by run:
  Int_t runNumber = 246087; //fAOD->GetRunNumber();
 
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

  Double_t sumTPCQn2x[5] = {0,0,0};  //[0]= eta<0; [1]= eta>0; [2]= -0.8 < eta < 0.8
  Double_t sumTPCQn2y[5] = {0,0,0};
  Double_t sumTPCQn3x[5] = {0,0,0};  //[0]= eta<0; [1]= eta>0; [2]= -0.8 < eta < 0.8
  Double_t sumTPCQn3y[5] = {0,0,0};
  Double_t sumTPCQn4x[5] = {0,0,0};  //[0]= eta<0; [1]= eta>0; [2]= -0.8 < eta < 0.8
  Double_t sumTPCQn4y[5] = {0,0,0};
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
  /*
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

  */

  
  Float_t nSigPionTPC[40000]   = {0.,};
  Float_t nSigKaonTPC[40000]   = {0.,};
  Float_t nSigProtonTPC[40000] = {0.,};    
  Float_t nSigPionTOF[40000]   = {0.,};
  Float_t nSigKaonTOF[40000]   = {0.,};
  Float_t nSigProtonTOF[40000] = {0.,};  
  
  //if(ntracks > 40000)           return;     //Dont break segment for higher tracks:   
  

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


  /*
  for(Int_t iTrack = 0; iTrack < ntracks; iTrack++) { //-------------------------

    AliAODTrack* AODtrack =dynamic_cast<AliAODTrack*>(fVevent->GetTrack(iTrack));
    if(!AODtrack) continue;

    nSigPionTPC[iTrack]   = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kPion);
    nSigKaonTPC[iTrack]   = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kKaon);
    nSigProtonTPC[iTrack] = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kProton);

    nSigPionTOF[iTrack]   = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kPion);
    nSigKaonTOF[iTrack]   = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kKaon);
    nSigProtonTOF[iTrack] = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kProton);  

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
	  sumTPCQn3x[0] += WgtNUA*TMath::Cos(3*phiTrk);
	  sumTPCQn3y[0] += WgtNUA*TMath::Sin(3*phiTrk);
	  sumTPCQn4x[0] += WgtNUA*TMath::Cos(4*phiTrk);
	  sumTPCQn4y[0] += WgtNUA*TMath::Sin(4*phiTrk);
	  multEtaNeg++;
	  SumWEtaNeg += WgtNUA;
	}
	else if(etaTrk > 0.05){
	  sumTPCQn2x[1] += WgtNUA*TMath::Cos(gPsiN*phiTrk);
	  sumTPCQn2y[1] += WgtNUA*TMath::Sin(gPsiN*phiTrk);
	  sumTPCQn3x[1] += WgtNUA*TMath::Cos(3*phiTrk);
	  sumTPCQn3y[1] += WgtNUA*TMath::Sin(3*phiTrk);
	  sumTPCQn4x[1] += WgtNUA*TMath::Cos(4*phiTrk);
	  sumTPCQn4y[1] += WgtNUA*TMath::Sin(4*phiTrk);
	  multEtaPos++;
	  SumWEtaPos += WgtNUA;
	}
	sumTPCQn2x[3] += WgtNUA*TMath::Cos(gPsiN*phiTrk);
	sumTPCQn2y[3] += WgtNUA*TMath::Sin(gPsiN*phiTrk);
	sumTPCQn3x[3] += WgtNUA*TMath::Cos(3*phiTrk);
	sumTPCQn3y[3] += WgtNUA*TMath::Sin(3*phiTrk);
	sumTPCQn4x[3] += WgtNUA*TMath::Cos(4*phiTrk);
	sumTPCQn4y[3] += WgtNUA*TMath::Sin(4*phiTrk);
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

  */


  /*

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
  */



  //cout<<"After PU cut multTPC = "<<multTPC<<" multGlobal = "<<multGlobal<<endl;




  fHistEventCount->Fill(stepCount); //9
  stepCount++;



  fHistTPConlyVsCL1After->Fill(centrCL1,multTPCAll);
  fHistTPConlyVsV0MAfter->Fill(centrV0M,multTPCAll);
  fHistGlobalVsV0MAfter->Fill(centrV0M, multGlobal);


  // MC corrected Refmult:
  fHistRawVsCorrMultFB->Fill(RefMultRawFB,RefMultCorrFB); // FB set by AddTask.. 


  Float_t EvtCent = centrality;




  // if(gPsiN > 0 && (multEtaNeg<2 || multEtaPos<2)) return;  //Minimum 2 tracks in each eta

  fHistEventCount->Fill(stepCount); //10
  stepCount++;

  //--------------------------------------------------------


















  //--------------- cent CL1 <= 90 cut ------------------
  Int_t icentV0Qn = centrCL1; // cent CL1 used for V0 calibration.
  icentV0Qn += 1;

  // if(icentV0Qn>90)     return;
  
  fHistEventCount->Fill(stepCount); //11
  stepCount++;



  /*

  //-------- V0M info ---------------
  const AliAODVZERO *fAODV0 = fAOD->GetVZEROData();

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

    for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA

      fMultv0 = fAODV0->GetMultiplicity(iV0);
 
      if(fHCorrectV0M){
	fMultv0 = fMultv0 * fHCorrectV0M->GetBinContent(iV0+1);  // Gain Correction
	//cout<<"info: run = "<<runNumber<<" cent = "<<centrCL1<<"\t channel = "<<iV0<<" gain = "<<fHCorrectV0M->GetBinContent(iV0+1)<<endl;
      }
     
      fV0MultChVsRun->Fill(iV0+0.5,centrCL1,fMultv0);
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


  */












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
  Double_t nSigTOFpion2,  nSigTPCpion2;
  Double_t nSigTOFkaon2,  nSigTPCkaon2;
  Double_t nSigTOFproton2,nSigTPCproton2;

  //Tof variables
  //Double_t length, tofTime, probMis, mass, beta;
  //Double_t c = TMath::C()*1.E-9;                //bright light m/ns 
  //Int_t    TOFmatch=0;



  Int_t    ptBin,gCharge1,gCharge2;

  //Double_t dcaXY, dcaZ ;


  //----------- Set the desired Harmonic ------------
  Int_t n = gN;
  Int_t m = gM;
  Int_t p =n+m;
  //------------------------------------------------


  Int_t skipPairHBT = 0;





  Double_t Chi2Trk1,ptwPion1,ptwKaon1,ptwProton1; 
  Double_t Chi2Trk2,ptwPion2,ptwKaon2,ptwProton2;
  Double_t wNUAPion1,wNUAKaon1,wNUAProton1;
  Double_t wNUAPion2,wNUAKaon2,wNUAProton2;

  Double_t    WgtEPPion = 1.0;
  Double_t    WgtEPKaon = 1.0;
  Double_t    WgtEPProton = 1.0;


  Bool_t isPion1 = kFALSE;
  Bool_t isKaon1  = kFALSE;
  Bool_t isProton1 = kFALSE;
  Bool_t isPion2 = kFALSE;
  Bool_t isKaon2  = kFALSE;
  Bool_t isProton2 = kFALSE;

  Int_t multPOI1st = 0;
  Int_t multPOI2nd = 0;







  //Correct with Centrality Wgts :   in second pass 
  Int_t iCentBinWgt = (Int_t) centrality;
  iCentBinWgt += 1;

  Double_t fWgtCent = 1.0;

  if(fHCentWeightForRun){
    fWgtCent = fHCentWeightForRun->GetBinContent(iCentBinWgt);
  }

  //Fill Centrality for run-by-run:  in first pass over data
  fCentDistAfter->Fill(centrality); 
  






  // cout<<" Info:UserExec()  I am before ... AliAODMCHeader, cent = "<<centrality<<" \n";





  // For productions with injected signals, figure out above which label to skip particles/tracks
  Int_t skipParticlesAbove = 0;

  AliGenEventHeader* eventHeader = 0;
  //Int_t headers = 0;
    
  // AOD only
  // AliAODMCHeader* header = (AliAODMCHeader*) fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  // if(!header) {
  //   printf("fInjectedSignals set but no MC header found");
  //   return;
  // }
    
  //headers     = header->GetNCocktailHeaders();
  //eventHeader = header->GetCocktailHeader(0);
    
  // if(!eventHeader)
  // {
  //   return;
  // }
    
  //skipParticlesAbove = eventHeader->NProduced();
  //AliInfo(Form("Injected signals in this event (%d headers). Keeping particles/tracks of %s. Will skip particles/tracks above %d.", headers, eventHeader->ClassName(), skipParticlesAbove)); 
  




















  //++++++++++++++++++ MC tracks+++++++++++++++++++++//
  Int_t nMCLabelCounter  = 0;

  //cout<<mcEvent<<endl;

  //Int_t nMCParticles = mcEvent->GetNumberOfTracks();
  Int_t nMCParticles = mcEvent->GetNumberOfPrimaries();

  // cout<<" number of mc primaries = "<<nMCParticles<<endl;

  TArrayI labelMCArray(nMCParticles);
  Bool_t labelTPC = kTRUE;     
  Int_t pdgIndex = 0;
  Int_t LabelOfMother = 0;

  Double_t ReQn[3] = {0.,};
  Double_t ImQn[3] = {0.,};
 
  Double_t ReQ1P=0, ReQ1N=0, ImQ1P=0, ImQ1N=0;// = 0;
  Double_t ReQ2P=0, ImQ2P=0, ReQ2N=0, ImQ2N=0;// = 0;
  Double_t ReQ3P=0, ImQ3P=0, ReQ3N=0, ImQ3N=0;
  Double_t ReQ4P=0, ImQ4P=0, ReQ4N=0, ImQ4N=0;
  Double_t CMEres = 0, MQP=0,   MQN=0,   MQ=0;
  Double_t MEtaPos = 0, MEtaNeg = 0;
  Double_t QxEP[3] = {0.,};
  Double_t QyEP[3] = {0.,};
  Double_t QxEPPos[3] = {0.,};
  Double_t QyEPPos[3] = {0.,};
  Double_t QxEPNeg[3] = {0.,};
  Double_t QyEPNeg[3] = {0.,};

  //Array for Nested Loop:
  //Double_t pTMCtrack[5000] = {0,}; //dont need pT
  Double_t PhiMCtrack[5000] = {0,};
  Double_t EtaMCtrack[5000] = {0,};
  Int_t ChMCtrack[5000] = {0,};
  Int_t iLoop = 0;



  for (Int_t iTracks = 0; iTracks < nMCParticles; iTracks++) {
    //AliAODMCParticle *mcTrack = (AliAODMCParticle*) mcEvent->GetTrack(iTracks); 
    AliMCParticle* mcTrack = dynamic_cast<AliMCParticle *> (mcEvent->GetTrack(iTracks));
    if(!mcTrack) {
      AliError(Form("ERROR: Could not receive track %d (mc loop)", iTracks));
      continue;
    }

    if(!(mcEvent->IsPhysicalPrimary(iTracks))) continue;
	
    //exclude particles generated out of the acceptance
    Double_t vz = mcTrack->Zv();
    if(TMath::Abs(vz) > 50.0) continue;
    //acceptance
    //if(mcTrack->Eta() > fMaxEtaCut || mcTrack->Eta() < fMinEtaCut) 
    //continue;
    Double_t y = mcTrack->Y();
    if(TMath::Abs(y) > 0.8) continue;  //  why do I need this cut?
    //if(mcTrack->Pt() > 5.0 || mcTrack->Pt() < 0.1)   continue;

    //TString generatorName;
    //Bool_t hasGenerator = mcEvent->GetCocktailGenerator(iTracks,generatorName);
    //cout<<" generator Name : "<<generatorName.Data()<<endl;
    //if((!hasGenerator) || (!generatorName.Contains("AMPT"))) continue;
		
    LabelOfMother = mcTrack->GetMother();
    if(LabelOfMother>0) continue;

    //if(!mcTrack->IsPhysicalPrimary()) continue; works for AODmc


    //--------- Now load mcTrack parameters ------------
    dPt1     = mcTrack->Pt();
    gCharge1 = mcTrack->Charge();

    if(gCharge1>0) gCharge1 = 1;
    else if(gCharge1<0) gCharge1 = -1;

    Int_t pdgcode = TMath::Abs(mcTrack->PdgCode());

    //cout<<" pdgcode : "<<pdgcode<<endl;



    /*
    //check PID strategy
    switch (pdgcode) {
    case 11:
      pdgIndex = 0; //e
      //hTrue[pdgIndex][chargeBin]->Fill(pt);
      break;
    case 13: 
      pdgIndex = 1; //mu 
      //hTrue[pdgIndex][chargeBin]->Fill(pt);
      break;
    case 211:
      pdgIndex = 2; //pi 
      fHistPtwithoutPIDcut[0][cIndex]->Fill(dPt1*gCharge1); //Pion
      fHistPtwithoutPIDcut[3][cIndex]->Fill(dPt1*gCharge1); //charge
      break;	
    case 321:
      pdgIndex = 3; //K
      //hTrue[pdgIndex][chargeBin]->Fill(pt);
      fHistPtwithoutPIDcut[1][cIndex]->Fill(dPt1*gCharge1); //Kaon
      fHistPtwithoutPIDcut[3][cIndex]->Fill(dPt1*gCharge1); //charge
      break;	
    case 2212:
      pdgIndex = 4; //p
      //hTrue[pdgIndex][chargeBin]->Fill(pt);
      fHistPtwithoutPIDcut[2][cIndex]->Fill(dPt1*gCharge1); //Proton
      fHistPtwithoutPIDcut[3][cIndex]->Fill(dPt1*gCharge1); //charge
      break;
    } */

    if(pdgcode==211 || pdgcode==321 || pdgcode==2212){

      //labelMCArray.AddAt(iTracks,nMCLabelCounter);

      //------------- For v2 and gamma/delta calculation -----------
      if(dPt1>0.2 && dPt1<10.0){

	dEta1 = mcTrack->Eta();
 
	if(TMath::Abs(dEta1)<0.8){
	dPhi1 = mcTrack->Phi();
        //cout<<iTracks<<" pt = "<<dPt1<<"\t phi = "<<dPhi1<<"\t dEta1 = "<<dEta1<<"\tch = "<<gCharge1<<endl;

	//Storing for next loop
	PhiMCtrack[iLoop] = dPhi1;
	EtaMCtrack[iLoop] = dEta1;
	ChMCtrack[iLoop] = gCharge1;
	iLoop++;

	  for(int h=0; h<3; h++) {
	    ReQn[h] += TMath::Cos((h+2.)*dPhi1);
	    ImQn[h] += TMath::Sin((h+2.)*dPhi1);
	    //TPC Event plane for |n|<0.4
	    if(dEta1>0){
	      QxEPPos[h] += TMath::Cos((h+2.)*dPhi1);
	      QyEPPos[h] += TMath::Sin((h+2.)*dPhi1);
	      MEtaPos++;
	    }
	    else if(dEta1<0){
	      QxEPNeg[h] += TMath::Cos((h+2.)*dPhi1);
	      QyEPNeg[h] += TMath::Sin((h+2.)*dPhi1);
	      MEtaNeg++;
	    }
	  }
	  MQ++;



	  if(gCharge1 > 0) {
	    ReQ1P += TMath::Cos(dPhi1);
	    ImQ1P += TMath::Sin(dPhi1);
	    ReQ2P += TMath::Cos(2*dPhi1);
	    ImQ2P += TMath::Sin(2*dPhi1);
	    ReQ3P += TMath::Cos(3*dPhi1);
	    ImQ3P += TMath::Sin(3*dPhi1);
	    ReQ4P += TMath::Cos(4*dPhi1);
	    ImQ4P += TMath::Sin(4*dPhi1);
	    MQP++;
	  } else {
	    ReQ1N += TMath::Cos(dPhi1);
	    ImQ1N += TMath::Sin(dPhi1);
	    ReQ2N += TMath::Cos(2*dPhi1);
	    ImQ2N += TMath::Sin(2*dPhi1);
	    ReQ3N += TMath::Cos(3*dPhi1);
	    ImQ3N += TMath::Sin(3*dPhi1);
	    ReQ4N += TMath::Cos(4*dPhi1);
	    ImQ4N += TMath::Sin(4*dPhi1);
	    MQN++;
	  }	
	}//global EtaCut for tracks

      }
      //---------------------- upto this (v2,gamma) -----------------------
    }//PID check using PDG value

    if(nMCLabelCounter>40000) break;
    nMCLabelCounter++;
   
  }//------------ mc track loop

  // cout<<""<<" MQ = "<<MQ<<"\t MQP = "<<MQP<<"\t MQN = "<<MQN<<"\t centrality = "<<centrality<<endl;





  Double_t newQx2,newQx3,newQx4;
  Double_t newQy2,newQy3,newQy4;

  Double_t newQx2b,newQx3b,newQx4b;
  Double_t newQy2b,newQy3b,newQy4b;

  Double_t Psi2TPC=0,Psi3TPC=0,Psi4TPC=0;

  //Now nested loop for gamma Correlator:
  for(int i=0;i<iLoop;i++){

    dPhi1 = PhiMCtrack[i];
    dEta1 = EtaMCtrack[i];
    gCharge1 = ChMCtrack[i];

    //'newQn' should be restored to original value:
    newQx2 = QxEPPos[0];
    newQx3 = QxEPPos[1];
    newQx4 = QxEPPos[2];

    newQy2 = QyEPPos[0];
    newQy3 = QyEPPos[1];
    newQy4 = QyEPPos[2];

    if(dEta1>0){ //remove 1st track
      newQx2 -= TMath::Cos(2.*dPhi1);
      newQx3 -= TMath::Cos(3.*dPhi1);
      newQx4 -= TMath::Cos(4.*dPhi1);

      newQy2 -= TMath::Sin(2.*dPhi1);
      newQy3 -= TMath::Sin(3.*dPhi1);
      newQy4 -= TMath::Sin(4.*dPhi1);
    } 

    //---------- 2nd loop -----------
    for(int j=0;j<iLoop;j++){

      if(i==j) continue; 

      dPhi2 = PhiMCtrack[j];
      dEta2 = EtaMCtrack[j];
      gCharge2 = ChMCtrack[j];

      //'newQnb' should be restored to 1st track removed value;
      newQx2b = newQx2;
      newQx3b = newQx3;
      newQx4b = newQx4;

      newQy2b = newQy2;
      newQy3b = newQy3;
      newQy4b = newQy4;

      if(dEta2>0){ //remove 2nd track
	newQx2b -= TMath::Cos(2.*dPhi2);
	newQx3b -= TMath::Cos(3.*dPhi2);
	newQx4b -= TMath::Cos(4.*dPhi2);

	newQy2b -= TMath::Sin(2.*dPhi2);
	newQy3b -= TMath::Sin(3.*dPhi2);
	newQy4b -= TMath::Sin(4.*dPhi2);
      } 

      //now form the EPs:
      Psi2TPC = (1/2.)*(TMath::ATan2(newQy2b,newQx2b)); 
      if(Psi2TPC<0.) Psi2TPC += 2*TMath::Pi()/2.0;

      Psi3TPC = (1/3.)*(TMath::ATan2(newQy3b,newQx3b) + TMath::Pi()); 

      Psi4TPC = (1/4.)*(TMath::ATan2(newQy4b,newQx4b));
      if(Psi4TPC<0.) Psi4TPC += 2*TMath::Pi()/4.0;

      if(gCharge1!=gCharge2){//Unlike Sign
	fCMEPro[1]->Fill(centrality, TMath::Cos( dPhi1  + dPhi2  - 2*Psi2TPC));
	fCMEPro[3]->Fill(centrality, TMath::Cos( dPhi1  + 2*dPhi2 - 3*Psi3TPC));
	fCMEPro[5]->Fill(centrality, TMath::Cos( dPhi1  - 3*dPhi2 + 2*Psi2TPC));
	fCMEPro[7]->Fill(centrality, TMath::Cos(2*dPhi1 + 2*dPhi2 - 4*Psi4TPC));
      } 
      else{ //Like Sign
	fCMEPro[0]->Fill(centrality, TMath::Cos( dPhi1  + dPhi2 - 2*Psi2TPC));
	fCMEPro[2]->Fill(centrality, TMath::Cos( dPhi1  + 2*dPhi2 - 3*Psi3TPC));
	fCMEPro[4]->Fill(centrality, TMath::Cos( dPhi1  - 3*dPhi2 + 2*Psi2TPC));
	fCMEPro[6]->Fill(centrality, TMath::Cos(2*dPhi1 + 2*dPhi2 - 4*Psi4TPC));
      }
    }//j-loop
  }//i-loop


  fCentDistBefore->Fill(iLoop);

  // cout<<"track for nested loop: "<<iLoop<<"\t MQP = "<<MQP<<"\t MQN = "<<MQN<<"\t centrality = "<<centrality<<endl;
  






  Double_t Psi2TPCA = (1/2.)*(TMath::ATan2(QyEPPos[0],QxEPPos[0]));
  if(Psi2TPCA<0.) Psi2TPCA += 2*TMath::Pi()/2.0;

  Double_t Psi3TPCA = (1/3.)*(TMath::ATan2(QyEPPos[1],QxEPPos[1]) + TMath::Pi()); 

  Double_t Psi4TPCA = (1/4.)*(TMath::ATan2(QyEPPos[2],QxEPPos[2]));
  if(Psi4TPCA<0.) Psi4TPCA += 2*TMath::Pi()/4.0;
  //------------------
  Double_t Psi2TPCC = (1/2.)*(TMath::ATan2(QyEPNeg[0],QxEPNeg[0]));
  if(Psi2TPCC<0.) Psi2TPCC += 2*TMath::Pi()/2.0;

  Double_t Psi3TPCC = (1/3.)*(TMath::ATan2(QyEPNeg[1],QxEPNeg[1]) + TMath::Pi()); 

  Double_t Psi4TPCC = (1/4.)*(TMath::ATan2(QyEPNeg[2],QxEPNeg[2]));
  if(Psi4TPCC<0.) Psi4TPCC += 2*TMath::Pi()/4.0;



  fHTPCAEventPlaneVsCent->Fill(centrality,Psi2TPCA);
  fHTPCCEventPlaneVsCent->Fill(centrality,Psi3TPCA);
  fHTPCEventPlaneVsCent->Fill(centrality, Psi4TPCA);

  fResolution_TPCEP[0]->Fill(centrality,TMath::Cos(2.*(Psi2TPCA - Psi2TPCC)));
  fResolution_TPCEP[1]->Fill(centrality,TMath::Cos(3.*(Psi3TPCA - Psi3TPCC)));
  fResolution_TPCEP[2]->Fill(centrality,TMath::Cos(4.*(Psi4TPCA - Psi4TPCC)));


  //Fill v2 values: Rihan
  if(MQ>1){
    for(int h=0; h<3; h++) {
      V2IntProQC[h]->Fill(centrality,(ReQn[h]*ReQn[h]+ImQn[h]*ImQn[h]-MQ)/(MQ*(MQ-1.)),MQ*(MQ-1.));
    }
  }

  if(MQP>1 && MQN>1) {
  //Fill Delta Correlations:
    CMEres = (ReQ1P*ReQ1P+ImQ1P*ImQ1P-MQP)/(MQP*MQP-MQP);
    fDnnPro[0]->Fill(centrality,CMEres,MQP*MQP-MQP);
    CMEres = (ReQ1P*ReQ1N+ImQ1P*ImQ1N)/(MQP*MQN);
    fDnnPro[1]->Fill(centrality,CMEres,MQP*MQN);

    CMEres = (ReQ2P*ReQ2P+ImQ2P*ImQ2P-MQP)/(MQP*MQP-MQP);
    fDnnPro[2]->Fill(centrality,CMEres,MQP*MQP-MQP);
    CMEres = (ReQ2P*ReQ2N+ImQ2P*ImQ2N)/(MQP*MQN);
    fDnnPro[3]->Fill(centrality,CMEres,MQP*MQN);

    CMEres = (ReQ3P*ReQ3P+ImQ3P*ImQ3P-MQP)/(MQP*MQP-MQP);
    fDnnPro[4]->Fill(centrality,CMEres,MQP*MQP-MQP);
    CMEres = (ReQ3P*ReQ3N+ImQ3P*ImQ3N)/(MQP*MQN);
    fDnnPro[5]->Fill(centrality,CMEres,MQP*MQN);

    CMEres = (ReQ4P*ReQ4P+ImQ4P*ImQ4P-MQP)/(MQP*MQP-MQP);
    fDnnPro[6]->Fill(centrality,CMEres,MQP*MQP-MQP);
    CMEres = (ReQ4P*ReQ4N+ImQ4P*ImQ4N)/(MQP*MQN);
    fDnnPro[7]->Fill(centrality,CMEres,MQP*MQN);


   //----------- Gamma Correlators ------------
    /*
   //C112
    CMEres = ((ReQ1P*ReQ1P-ImQ1P*ImQ1P-ReQ2P)*cos(2.*Psi2TPC)+(2.*ReQ1P*ImQ1P-ImQ2P)*sin(2.*Psi2TPC))/(MQP*MQP-MQP);
    fCMEPro[0]->Fill(centrality,CMEres,MQP*MQP-MQP);
    CMEres = ((ReQ1P*ReQ1N-ImQ1P*ImQ1N)*cos(2.*Psi2TPC)+(ReQ1P*ImQ1N+ImQ1P*ReQ1N)*sin(2.*Psi2TPC))/(MQP*MQN);
    fCMEPro[1]->Fill(centrality,CMEres,MQP*MQN);
   
    //Rihan C123
    CMEres = ((ReQ1P*ReQ2P-ImQ1P*ImQ2P-ReQ3P)*cos(3.*Psi3TPC)+(ReQ1P*ImQ2P+ReQ2P*ImQ1P-ImQ3P)*sin(3.*Psi3TPC))/(MQP*(MQP-1.));
    fCMEPro[2]->Fill(centrality,CMEres,MQP*(MQP-1.));
    CMEres = ((ReQ1P*ReQ2N-ImQ1P*ImQ2N)*cos(3.*Psi3TPC)+(ReQ1P*ImQ2N+ImQ1P*ReQ2N)*sin(3.*Psi3TPC))/(MQP*MQN);
    fCMEPro[3]->Fill(centrality,CMEres,MQP*MQN);

    //Rihan C132
    CMEres = ((ReQ1P*ReQ3P+ImQ1P*ImQ3P-ReQ2P)*cos(2.*Psi2TPC)+(ReQ1P*ImQ3P-ImQ1P*ReQ3P-ImQ2P)*sin(2.*Psi2TPC))/(MQP*(MQP-1.));
    fCMEPro[4]->Fill(centrality,CMEres,MQP*(MQP-1.));
    CMEres = ((ReQ1P*ReQ3N+ImQ1P*ImQ3N)*cos(2.*Psi2TPC)+(ReQ1P*ImQ3N-ReQ3N*ImQ1P)*sin(2.*Psi2TPC))/(MQP*MQN);
    fCMEPro[5]->Fill(centrality,CMEres,MQP*MQN);
   
    //Rihan C224
    CMEres = ((ReQ2P*ReQ2P-ImQ2P*ImQ2P-ReQ4P)*cos(4.*Psi4TPC)+(2.*ReQ2P*ImQ2P-ImQ4P)*sin(4.*Psi4TPC))/(MQP*MQP-MQP);
    fCMEPro[6]->Fill(centrality,CMEres,MQP*MQP-MQP);
    CMEres = ((ReQ2P*ReQ2N-ImQ2P*ImQ2N)*cos(4.*Psi4TPC)+(ReQ2P*ImQ2N+ReQ2N*ImQ2P)*sin(4.*Psi4TPC))/(MQP*MQN);
    fCMEPro[7]->Fill(centrality,CMEres,MQP*MQN);
    */
  }





  PostData(1,fListHist);

  fHistEventCount->Fill(14.5); //15th bin is last one
  stepCount++;







  //if(fEventCount%5==0) 
  //cout<<"Ev = "<<fEventCount<<"\tMult = "<<multEtaFull<<"\t nMCLabelCounter = "<<nMCLabelCounter<<"\t MatchedReco = "<< nMatchedRecoTrk <<endl;


  fEventCount++;
 

}//================ UserExec ==============






















//////////////////// FUNCTIONS //////////////////////




double AliAnalysisTaskCMEMC::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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


Bool_t AliAnalysisTaskCMEMC::PileUpMultiVertex(const AliAODEvent* faod)
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



void AliAnalysisTaskCMEMC::SetupMCcorrectionMap(){
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
    for(int i=0;i<10;i++) {
      fFB_Efficiency_Pion_Cent[i]   = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_Pion_%d",i));
      fFB_Efficiency_Kaon_Cent[i]   = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_Kaon_%d",i));
      fFB_Efficiency_Proton_Cent[i] = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_Proton_%d",i));
    }

    //------- Fill the flags: --------------
    if(fFB_Efficiency_Cent[0] && fFB_Efficiency_Cent[4]){
      fHistTaskConfigParameters->SetBinContent(19,1);
    }
    if(fFB_Efficiency_Pion_Cent[0] && fFB_Efficiency_Pion_Cent[4]){
      fHistTaskConfigParameters->SetBinContent(20,1);
    }
    if(fFB_Efficiency_Kaon_Cent[0] && fFB_Efficiency_Kaon_Cent[4]){
      fHistTaskConfigParameters->SetBinContent(21,1);
    }
    if(fFB_Efficiency_Proton_Cent[0] && fFB_Efficiency_Proton_Cent[4]){
      fHistTaskConfigParameters->SetBinContent(22,1);
    }
  }
  else if(!fListFBHijing){
    std::cout<<"\n\n !!!!**** Warning: FB Efficiency File/List not found Use Weight = 1.0 *****\n\n"<<std::endl;
    //exit(1);
  }
}


//____________________________________________________________________//
Bool_t AliAnalysisTaskCMEMC::IsLabelUsed(TArrayI labelArray, Int_t label) {
  //Checks if the label is used already
  Bool_t status = kFALSE;
  for(Int_t i = 0; i < labelArray.GetSize(); i++) {
    if(labelArray.At(i) == label)
      status = kTRUE;
  }

  return status;
}



Int_t AliAnalysisTaskCMEMC::GetCentralityScaled0to10(Double_t fCent){

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

void AliAnalysisTaskCMEMC::SetUpCentralityOutlierCut(){
  //std::cout<<" centrality outlier function called "<<std::endl;
 Float_t fMeanTPC[100] = {2902.95,2758.33,2642.78,2536.67,2435.37,2340.06,2248.44,2163.71,2080.49,2001.54,1925.86,1852.64,1781.97,1715.56,1650.53,1587.23,1527.51,1468.19,1412.73,1357.86,1305.35,1254.33,1205.57,1157.28,1111.53,1066.42,1023.15,981.594,940.795,901.766,863.651,826.183,790.53,756.358,722.654,690.513,659.443,628.807,599.748,571.664,544.446,518.042,492.369,468.072,444.694,422.487,400.104,379.129,359.147,339.62,320.817,302.788,285.791,269.015,253.688,238.671,224.039,209.932,196.915,184.647,172.76,161.381,150.395,140.288,131.033,121.58,113.112,104.938,97.3078,90.2178,83.5974,77.2645,70.7126,65.4424,60.1404,55.5644,50.8314,46.3761,43.024,38.625,35.3435,32.2304,29.4192,26.821,24.3303,21.9332,19.4215,16.7163,14.9414,13.1092,0.};

 Float_t fSigmaTPC[100] = {122.209,107.901,103.452,100.498,97.7403,94.7845,93.2543,90.0548,88.1106,85.7382,84.0812,82.2978,80.3817,78.6002,77.3448,75.5086,73.6842,71.9733,70.3447,69.1999,67.878,66.3511,65.0406,63.4866,62.4409,60.7899,59.1328,58.426,56.8618,55.8871,54.1031,53.4959,52.0482,51.0441,49.6218,48.7646,47.5166,46.5247,45.0727,44.4311,43.4531,42.0404,41.0238,40.1384,39.2588,38.2461,36.5951,36.0552,35.3727,33.7883,32.7167,32.4486,31.3709,30.3444,29.505,28.5139,27.4471,26.5359,25.9506,25.127,24.3797,23.2985,22.279,21.4698,20.781,20.8193,19.9509,18.8036,17.9145,16.961,16.7375,15.852,14.9324,14.7663,13.5969,13.4533,12.3067,12.7835,11.7283,10.6758,10.6676,10.6492,9.04614,8.89065,8.66093,8.50997,7.98812,6.91087,7.12045,7.29593,0.};

 for(int i=0;i<90;i++) {
   hCentvsTPCmultCuts->SetBinContent(i+1,1,fMeanTPC[i]);
   hCentvsTPCmultCuts->SetBinContent(i+1,2,fSigmaTPC[i]);
 }
}






void AliAnalysisTaskCMEMC::SetupEventAndTaskConfigInfo(){

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





void AliAnalysisTaskCMEMC::GetV0MCorrectionHist(Int_t run)
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







void AliAnalysisTaskCMEMC::GetNUACorrectionHist(Int_t run)
{

  if(fListNUACorr){
    for(int i=0;i<5;i++){
      fHCorrectNUApos[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pos_Cent%d_Run%d",i,run)); 
      fHCorrectNUAneg[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Neg_Cent%d_Run%d",i,run));
    }
    if(fHCorrectNUApos[3] && fHCorrectNUAneg[3]){
      cout<<"\n=========== Info:: Setting up NUA corrections for run "<<run<<"============"<<endl;
      fHistTaskConfigParameters->SetBinContent(15,1);
    }
  }
  else {
    printf("\n ******** Warning: No NUA Correction for Charge Particle in run %d, Use Wgt = 1.0 ********* \n",run);
    //Do not allocate memory:
    /*
    for(int i=0;i<5;i++){
      fHCorrectNUApos[i] = new TH3D(Form("fHCorrectNUApos_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAneg[i] = new TH3D(Form("fHCorrectNUAneg_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUApos[i]->SetBinContent(1,1,1,1.0);
      fHCorrectNUAneg[i]->SetBinContent(1,1,1,1.0);
      //exit(1);
    }*/
  }

  //===================   PID:  ==========================
  if(fListNUACorr){
    for(int i=0;i<5;i++){
      fHCorrectNUAposPion[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Pos_Cent%d_Run%d",i,run)); //
      fHCorrectNUAnegPion[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Neg_Cent%d_Run%d",i,run));
      if(fHCorrectNUAposPion[3] && fHCorrectNUAnegPion[3])   fHistTaskConfigParameters->SetBinContent(16,1);

      fHCorrectNUAposKaon[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Pos_Cent%d_Run%d",i,run)); //
      fHCorrectNUAnegKaon[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Neg_Cent%d_Run%d",i,run));
      if(fHCorrectNUAposKaon[3] && fHCorrectNUAnegKaon[3])   fHistTaskConfigParameters->SetBinContent(17,1);

      fHCorrectNUAposProton[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Pos_Cent%d_Run%d",i,run)); 
      fHCorrectNUAnegProton[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Neg_Cent%d_Run%d",i,run));
      if(fHCorrectNUAposProton[3] && fHCorrectNUAnegProton[3])   fHistTaskConfigParameters->SetBinContent(18,1);
    }
    if(fHCorrectNUAposPion[3] && fHCorrectNUAposKaon[3] && fHCorrectNUAposProton[3]) {
      cout<<"\n=========== Info:: Setting up --> PID NUA corrections for run =  "<<run<<"============"<<endl;
    }
    else{
      cout<<"\n=========== WARNING :: PID NUA corrections NOT found for run =  "<<run<<"============"<<endl;
    }
  }
  else {
    printf("\n ******** Error/Warning: No NUA Correction found for PID for run %d, Use PID Wgt = 1.0 ********* \n",run);
    /* 
    for(int i=0;i<5;i++){
      fHCorrectNUAposPion[i] = new TH3D(Form("fHCorrectNUAPionpos_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposPion[i] = new TH3D(Form("fHCorrectNUAPionneg_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposPion[i]->SetBinContent(1,1,1,1.0);
      fHCorrectNUAposPion[i]->SetBinContent(1,1,1,1.0);

      fHCorrectNUAposKaon[i] = new TH3D(Form("fHCorrectNUAKaonpos_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposKaon[i] = new TH3D(Form("fHCorrectNUAKaonneg_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposKaon[i]->SetBinContent(1,1,1,1.0);
      fHCorrectNUAposKaon[i]->SetBinContent(1,1,1,1.0);

      fHCorrectNUAposProton[i] = new TH3D(Form("fHCorrectNUAProtonpos_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposProton[i] = new TH3D(Form("fHCorrectNUAProtonneg_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposProton[i]->SetBinContent(1,1,1,1.0);
      fHCorrectNUAposProton[i]->SetBinContent(1,1,1,1.0);
      //exit(1);
      }*/
  }

}

