/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////
// CME correlator study using VO.
// Author: Rihan Haque (mhaque@cern.ch)
///////////////////////////////////////////////



#include <TGrid.h>
#include <TFile.h>
#include <TList.h>
#include "TMatrixDSym.h"
#include "Riostream.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliVVertex.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliFlowEvent.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskCMEV0.h"


using std::endl;
using std::cout;

ClassImp(AliAnalysisTaskCMEV0)

AliAnalysisTaskCMEV0::AliAnalysisTaskCMEV0(const TString name):AliAnalysisTaskSE(name),
  fEvent(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHistos(NULL),
  fListCalibs(NULL),
  fListNUAHist(NULL),
  fListFBHijing(NULL),
  fListNUACorr(NULL),
  fListZDNCorr(NULL),
  fListV0MCorr(NULL),
  fRejectPileUp(kTRUE),
  fRejectPileUpTight(kTRUE),
  bFillAvgTPCQn(kFALSE),
  bFillEtaPhiNUA(kFALSE),
  bApplyNUACorr(kFALSE),
  bApplyZDCCorr(kFALSE),
  bApplyNUAforEP(kFALSE),
  bFillZDCinfo(kFALSE),
  bSkipNestedTrk(kFALSE),
  bRemNegTrkRndm(kFALSE),
  bApplyV0MCorr(kFALSE),
  sDataSet("2015"),
  sAnalysisSet("DoGainEq"),
  sCentEstimator("V0"),
  sFileNUA("NewR"),
  sMCdimension("1D"),
  fRunFlag(0),
  fOldRunNum(0),
  fievent(0),
  EvtCent(0.0),
  fHarmonicN(1),
  fHarmonicM(1),
  fHarmonicPsi(2),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  gFilterBit(96),
  fHBTCutValue(0.0),
  fRefMultCorr(0),
  fRefMultRaw(0),
  fHist_Event_count(NULL),
  fPileUpMultSelCount(NULL),
  fPileUpCount(NULL),
  fTaskConfigParm(NULL),
  fMultV0(NULL),
  fQxnmV0A(NULL), 
  fQynmV0A(NULL),     
  fQxnsV0A(NULL),       
  fQynsV0A(NULL),       
  fQxnmV0C(NULL),        
  fQynmV0C(NULL),    
  fQxnsV0C(NULL),     
  fQynsV0C(NULL),
  fHV0AEventPlaneVsCent(NULL),
  fHV0CEventPlaneVsCent(NULL),
  fHTPCEventPlaneVsCent(NULL),
  fHEnergyZNCvsCent(NULL),
  fHEnergyZNAvsCent(NULL),
  fHEnergyZPCvsCent(NULL),
  fHEnergyZPAvsCent(NULL),
  fHEnergyZNCvsCentRun(NULL),
  fHEnergyZNAvsCentRun(NULL),
  fHEnergyZPCvsCentRun(NULL),
  fHEnergyZPAvsCentRun(NULL),
  fHEnergyZPCvsZPA(NULL),
  fHEnergyZNCvsZNA(NULL),
  hUnderOverBinNUApos(NULL),
  hUnderOverBinNUAneg(NULL),
  fHCentBinTrkRecenter(NULL),
  fHCorrectZDNP(NULL),
  fV0AQ2xVsCentRun(NULL),
  fV0AQ2yVsCentRun(NULL),
  fV0CQ2xVsCentRun(NULL),
  fV0CQ2yVsCentRun(NULL),
  fV0AQ3xVsCentRun(NULL),
  fV0AQ3yVsCentRun(NULL),
  fV0CQ3xVsCentRun(NULL),
  fV0CQ3yVsCentRun(NULL),
  fTPCQ2xVsCentRun(NULL),
  fTPCQ2yVsCentRun(NULL),
  fTPCQ3xVsCentRun(NULL),
  fTPCQ3yVsCentRun(NULL),
  mListNUAPos(NULL),
  mListNUANeg(NULL),
  fileNUApos(NULL),
  fileNUAneg(NULL),
  fAvgMultCentRun(NULL),
  fAvgWgtMultCentRun(NULL),
  fAvgPOIposCentRun(NULL),
  fAvgPOInegCentRun(NULL),
  fAvgPOIPPCentRun(NULL),
  fAvgPOINNCentRun(NULL),
  fAvgPOIOSCentRun(NULL),
  fV0MultChVsRun(NULL),
  fEventStatvsRun(NULL),
  fEtaBinFinderForQA(NULL),
  fVzBinFinderForNUA(NULL),
  fHistVtxZvsRun(NULL),
  fHistVtxXvsRun(NULL),
  fHistVtxYvsRun(NULL),
  fRejectRatioVsCR(NULL),
  fCentDistvsRun(NULL),
  fHCorrectV0M(NULL),
  fHAvgerageQnV0A(NULL),
  fHAvgerageQnV0C(NULL),
  fCentV0MvsVzRun(NULL),
  fCent3pvsVzRun(NULL),
  fRefMultCorrvsRaw(NULL),
  fTPCvsGlobalTrk(NULL),
  fTPCuncutvsGlobal(NULL),
  fGlobalTracks(NULL),
  fTPCvsITSfb96(NULL),
  fTPCvsITSfb32(NULL),
  fTPCFEvsITSfb96(NULL),
  fCentCL1vsVzRun(NULL),
  fdPhiFemtoCut(NULL),
  fVzDistribuion(NULL)
{
  for(int i=0;i<5;i++){
    fHCorrectNUApos[i] = NULL;
    fHCorrectNUAneg[i] = NULL;
  }
  for(int i=0;i<90;i++){
    runNums[i] = 0;
    for(int j=0;j<5;j++){
     fHist3DEtaPhiVz_Pos_Run[j][i] = NULL;
     fHist3DEtaPhiVz_Neg_Run[j][i] = NULL;
    }
  }
  for(int i=0;i<10;i++){
    fFB_Efficiency_Cent[i] = NULL;
    fFB_Efficiency_Pos[i] = NULL;
    fFB_Efficiency_Neg[i] = NULL;
  }
  for(int i=0;i<2;i++){
    for(int j=0;j<3;j++){
      fHist_Corr3p_SP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_SP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_SP_Norm_NN[i][j]  =  NULL;
      fHist_Reso2n_SP_Norm_Det[i][j] =  NULL;
    }  
  }
  for(int i=0;i<2;i++){
    for(int j=0;j<3;j++){
      fHist_Corr3p_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_NN[i][j]  =  NULL;
      fHist_Reso2n_EP_Norm_Det[i][j] =  NULL;
    }  
    for(int j=0;j<3;j++){
      fHist_Corr3p_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr3p_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr3p_EP_Refm_NN[i][j]  =  NULL;
      fHist_Reso2n_EP_Refm_Det[i][j] =  NULL;
    }  
  }
  for(int i=0;i<2;i++){
    fHist_Corr3p_vsRun_EP_PN[i] = NULL;
    fHist_Corr3p_vsRun_EP_PP[i] = NULL;
    fHist_Corr3p_vsRun_EP_NN[i] = NULL;
  }
  for(int i=0;i<3;i++){
    fHist_Corr3p_ZDN_SP_PN[i]  =  NULL;
    fHist_Corr3p_ZDN_SP_PP[i]  =  NULL;
    fHist_Corr3p_ZDN_SP_NN[i]  =  NULL;
    fHist_Reso2n_ZDN_SP_Det[i] = NULL;
  }
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
  for(int i=0;i<10;i++){
    fHistChPosvsEtaPtRun[i] = NULL;
    fHistChNegvsEtaPtRun[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_Corr3p_QAEta_SP_V0A_PN[i] = NULL;
    fHist_Corr3p_QAEta_SP_V0A_PP[i] = NULL;
    fHist_Corr3p_QAEta_SP_V0A_NN[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_Corr3p_QAEta_SP_V0C_PN[i] = NULL;
    fHist_Corr3p_QAEta_SP_V0C_PP[i] = NULL;
    fHist_Corr3p_QAEta_SP_V0C_NN[i] = NULL;
  }
  for(int i=0;i<4;i++){
    fHCos1nPosChEtaVz[i] = NULL;
    fHCos2nPosChEtaVz[i] = NULL;
    fHCos3nPosChEtaVz[i] = NULL;
    fHCos4nPosChEtaVz[i] = NULL;
    fHSin1nPosChEtaVz[i] = NULL;
    fHSin2nPosChEtaVz[i] = NULL;
    fHSin3nPosChEtaVz[i] = NULL;
    fHSin4nPosChEtaVz[i] = NULL;

    fHCos1nNegChEtaVz[i] = NULL;
    fHCos2nNegChEtaVz[i] = NULL;
    fHCos3nNegChEtaVz[i] = NULL;
    fHCos4nNegChEtaVz[i] = NULL;
    fHSin1nNegChEtaVz[i] = NULL;
    fHSin2nNegChEtaVz[i] = NULL;
    fHSin3nNegChEtaVz[i] = NULL;
    fHSin4nNegChEtaVz[i] = NULL;

    fHCos2nDWPosChEtaVz[i] = NULL;
    fHSin2nDWPosChEtaVz[i] = NULL;
    fHCos2nDWNegChEtaVz[i] = NULL;
    fHSin2nDWNegChEtaVz[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_NonIso_SP_PP_Mag0[i] = NULL;
    fHist_NonIso_SP_NN_Mag0[i] = NULL;
    fHist_NonIso_SP_PP_Mag1[i] = NULL;
    fHist_NonIso_SP_NN_Mag1[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_Corr2p_EP_Norm_PN[i] = NULL;
    fHist_Corr2p_EP_Norm_PP[i] = NULL;
    fHist_Corr2p_EP_Norm_NN[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_Corr2p_EP_Refm_PN[i] = NULL;
    fHist_Corr2p_EP_Refm_PP[i] = NULL;
    fHist_Corr2p_EP_Refm_NN[i] = NULL;
  }

  DefineInput(1, AliFlowEventSimple::Class()); // Input slot #1: AliFlowEventSimple

  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
  DefineOutput(3,TList::Class());

}//-------- real ---------

AliAnalysisTaskCMEV0::AliAnalysisTaskCMEV0(): AliAnalysisTaskSE(),
  fEvent(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHistos(NULL),
  fListCalibs(NULL),
  fListNUAHist(NULL),
  fListFBHijing(NULL),
  fListNUACorr(NULL),
  fListZDNCorr(NULL),
  fListV0MCorr(NULL),
  fRejectPileUp(kTRUE),
  fRejectPileUpTight(kTRUE),
  bFillAvgTPCQn(kFALSE),
  bFillEtaPhiNUA(kFALSE),
  bApplyNUACorr(kFALSE),
  bApplyZDCCorr(kFALSE),
  bApplyNUAforEP(kFALSE),
  bFillZDCinfo(kFALSE),
  bSkipNestedTrk(kFALSE),
  bRemNegTrkRndm(kFALSE),
  bApplyV0MCorr(kFALSE),
  sDataSet("2015"),
  sAnalysisSet("DoGainEq"),
  sCentEstimator("V0"),
  sFileNUA("NewR"),
  sMCdimension("1D"),
  fRunFlag(0),
  fOldRunNum(0),
  fievent(0),
  EvtCent(0.0),
  fHarmonicN(1),
  fHarmonicM(1),
  fHarmonicPsi(2),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  gFilterBit(96),
  fHBTCutValue(0.0),
  fRefMultCorr(0),
  fRefMultRaw(0),
  fHist_Event_count(NULL),
  fPileUpMultSelCount(NULL),
  fPileUpCount(NULL),
  fTaskConfigParm(NULL),
  fMultV0(NULL),
  fQxnmV0A(NULL), 
  fQynmV0A(NULL),     
  fQxnsV0A(NULL),       
  fQynsV0A(NULL),       
  fQxnmV0C(NULL),        
  fQynmV0C(NULL),    
  fQxnsV0C(NULL),     
  fQynsV0C(NULL),
  fHV0AEventPlaneVsCent(NULL),
  fHV0CEventPlaneVsCent(NULL),
  fHTPCEventPlaneVsCent(NULL),
  fHEnergyZNCvsCent(NULL),
  fHEnergyZNAvsCent(NULL),
  fHEnergyZPCvsCent(NULL),
  fHEnergyZPAvsCent(NULL),
  fHEnergyZNCvsCentRun(NULL),
  fHEnergyZNAvsCentRun(NULL),
  fHEnergyZPCvsCentRun(NULL),
  fHEnergyZPAvsCentRun(NULL),
  fHEnergyZPCvsZPA(NULL),
  fHEnergyZNCvsZNA(NULL),
  hUnderOverBinNUApos(NULL),
  hUnderOverBinNUAneg(NULL),
  fHCentBinTrkRecenter(NULL),
  fHCorrectZDNP(NULL),
  fV0AQ2xVsCentRun(NULL),
  fV0AQ2yVsCentRun(NULL),
  fV0CQ2xVsCentRun(NULL),
  fV0CQ2yVsCentRun(NULL),
  fV0AQ3xVsCentRun(NULL),
  fV0AQ3yVsCentRun(NULL),
  fV0CQ3xVsCentRun(NULL),
  fV0CQ3yVsCentRun(NULL),
  fTPCQ2xVsCentRun(NULL),
  fTPCQ2yVsCentRun(NULL),
  fTPCQ3xVsCentRun(NULL),
  fTPCQ3yVsCentRun(NULL),
  mListNUAPos(NULL),
  mListNUANeg(NULL),
  fileNUApos(NULL),
  fileNUAneg(NULL),
  fAvgMultCentRun(NULL),
  fAvgWgtMultCentRun(NULL),
  fAvgPOIposCentRun(NULL),
  fAvgPOInegCentRun(NULL),
  fAvgPOIPPCentRun(NULL),
  fAvgPOINNCentRun(NULL),
  fAvgPOIOSCentRun(NULL),
  fV0MultChVsRun(NULL),
  fEventStatvsRun(NULL),
  fEtaBinFinderForQA(NULL),
  fVzBinFinderForNUA(NULL),
  fHistVtxZvsRun(NULL),
  fHistVtxXvsRun(NULL),
  fHistVtxYvsRun(NULL),
  fRejectRatioVsCR(NULL),
  fCentDistvsRun(NULL),
  fHCorrectV0M(NULL),
  fHAvgerageQnV0A(NULL),
  fHAvgerageQnV0C(NULL),
  fCentV0MvsVzRun(NULL),
  fCent3pvsVzRun(NULL),
  fRefMultCorrvsRaw(NULL),
  fTPCvsGlobalTrk(NULL),
  fTPCuncutvsGlobal(NULL),
  fGlobalTracks(NULL),
  fTPCvsITSfb96(NULL),
  fTPCvsITSfb32(NULL),
  fTPCFEvsITSfb96(NULL),
  fCentCL1vsVzRun(NULL),
  fdPhiFemtoCut(NULL),
  fVzDistribuion(NULL)
{
  for(int i=0;i<5;i++){
    fHCorrectNUApos[i] = NULL;
    fHCorrectNUAneg[i] = NULL;
  }
  for(int i=0;i<90;i++){
    runNums[i] = 0;
    for(int j=0;j<5;j++){
     fHist3DEtaPhiVz_Pos_Run[j][i] = NULL;
     fHist3DEtaPhiVz_Neg_Run[j][i] = NULL;
    }
  }
  for(int i=0;i<10;i++){
    fFB_Efficiency_Cent[i] = NULL;
    fFB_Efficiency_Pos[i] = NULL;
    fFB_Efficiency_Neg[i] = NULL;
  }
  for(int i=0;i<2;i++){
    for(int j=0;j<3;j++){
      fHist_Corr3p_SP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_SP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_SP_Norm_NN[i][j]  =  NULL;
      fHist_Reso2n_SP_Norm_Det[i][j] =  NULL;
    }  
  }
  for(int i=0;i<2;i++){
    for(int j=0;j<3;j++){
      fHist_Corr3p_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_NN[i][j]  =  NULL;
      fHist_Reso2n_EP_Norm_Det[i][j] =  NULL;
    }
    for(int j=0;j<3;j++){
      fHist_Corr3p_EP_Refm_PN[i][j]  =  NULL;
      fHist_Corr3p_EP_Refm_PP[i][j]  =  NULL;
      fHist_Corr3p_EP_Refm_NN[i][j]  =  NULL;
      fHist_Reso2n_EP_Refm_Det[i][j] =  NULL;
    }    
  }
  for(int i=0;i<2;i++){
    fHist_Corr3p_vsRun_EP_PN[i] = NULL;
    fHist_Corr3p_vsRun_EP_PP[i] = NULL;
    fHist_Corr3p_vsRun_EP_NN[i] = NULL;
  }
  for(int i=0;i<3;i++){
    fHist_Corr3p_ZDN_SP_PN[i]  =  NULL;
    fHist_Corr3p_ZDN_SP_PP[i]  =  NULL;
    fHist_Corr3p_ZDN_SP_NN[i]  =  NULL;
    fHist_Reso2n_ZDN_SP_Det[i] = NULL;
  }
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
  for(int i=0;i<10;i++){
    fHistChPosvsEtaPtRun[i] = NULL;
    fHistChNegvsEtaPtRun[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_Corr3p_QAEta_SP_V0A_PN[i] = NULL;
    fHist_Corr3p_QAEta_SP_V0A_PP[i] = NULL;
    fHist_Corr3p_QAEta_SP_V0A_NN[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_Corr3p_QAEta_SP_V0C_PN[i] = NULL;
    fHist_Corr3p_QAEta_SP_V0C_PP[i] = NULL;
    fHist_Corr3p_QAEta_SP_V0C_NN[i] = NULL;
  }
  for(int i=0;i<4;i++){
    fHCos1nPosChEtaVz[i] = NULL;
    fHCos2nPosChEtaVz[i] = NULL;
    fHCos3nPosChEtaVz[i] = NULL;
    fHCos4nPosChEtaVz[i] = NULL;
    fHSin1nPosChEtaVz[i] = NULL;
    fHSin2nPosChEtaVz[i] = NULL;
    fHSin3nPosChEtaVz[i] = NULL;
    fHSin4nPosChEtaVz[i] = NULL;

    fHCos1nNegChEtaVz[i] = NULL;
    fHCos2nNegChEtaVz[i] = NULL;
    fHCos3nNegChEtaVz[i] = NULL;
    fHCos4nNegChEtaVz[i] = NULL;
    fHSin1nNegChEtaVz[i] = NULL;
    fHSin2nNegChEtaVz[i] = NULL;
    fHSin3nNegChEtaVz[i] = NULL;
    fHSin4nNegChEtaVz[i] = NULL;

    fHCos2nDWPosChEtaVz[i] = NULL;
    fHSin2nDWPosChEtaVz[i] = NULL;
    fHCos2nDWNegChEtaVz[i] = NULL;
    fHSin2nDWNegChEtaVz[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_NonIso_SP_PP_Mag0[i] = NULL;
    fHist_NonIso_SP_NN_Mag0[i] = NULL;
    fHist_NonIso_SP_PP_Mag1[i] = NULL;
    fHist_NonIso_SP_NN_Mag1[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_Corr2p_EP_Norm_PN[i] = NULL;
    fHist_Corr2p_EP_Norm_PP[i] = NULL;
    fHist_Corr2p_EP_Norm_NN[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_Corr2p_EP_Refm_PN[i] = NULL;
    fHist_Corr2p_EP_Refm_PP[i] = NULL;
    fHist_Corr2p_EP_Refm_NN[i] = NULL;
  }

}//-----------




void AliAnalysisTaskCMEV0::UserCreateOutputObjects()
{

 this->InitializeRunArray(sDataSet);

 fAnalysisUtil = new AliAnalysisUtils();
 fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
 fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);


 if(fListFBHijing) {
  if(sMCdimension=="1D"){
    for(int i=0;i<10;i++) {
      fFB_Efficiency_Cent[i] = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_%d",i));
    }
  }
  else if(sMCdimension=="3D"){
    Int_t centFB[11] = {0,5,10,20,30,40,50,60,70,80,90};
    for(int i=0;i<10;i++) {
      fFB_Efficiency_Pos[i] = (TH3F *) fListFBHijing->FindObject(Form("fHistCorrectionPlus%d-%d",centFB[i],centFB[i+1]));
      fFB_Efficiency_Neg[i] = (TH3F *) fListFBHijing->FindObject(Form("fHistCorrectionMinus%d-%d",centFB[i],centFB[i+1]));
    }
  }
  else{ 
    printf("\n\n!!*****  Warning *****!!\n Enter correct Dimention !!\n\n"); exit(1);
  }
 }
 else{ // if MC efficiency not found then use weight = 1.
  printf("\n\n!!*****  Warning *****!!\n FilterBit efficiency not found, use = 1.0 !!\n\n");
   for(int i=0;i<10;i++){
     fFB_Efficiency_Cent[i] = new TH1D(Form("eff_unbiased_%d",i),"",1,0,50.); 
     fFB_Efficiency_Cent[i]->SetBinContent(1,1.0);
     fFB_Efficiency_Pos[i]  = new TH3F(Form("eff_unbiased_Pos_%d",i),"",1,-0.9,0.9,1,0,50,1,0,7);
     fFB_Efficiency_Pos[i]->SetBinContent(1,1,1,1.0);
     fFB_Efficiency_Neg[i]  = new TH3F(Form("eff_unbiased_Neg_%d",i),"",1,-0.9,0.9,1,0,50,1,0,7);
     fFB_Efficiency_Neg[i]->SetBinContent(1,1,1,1.0);
   } 
 }



 fListHistos = new TList();
 fListHistos->SetOwner(kTRUE);

 fListCalibs = new TList();
 fListCalibs->SetOwner(kTRUE);

 fListNUAHist = new TList();
 fListNUAHist->SetOwner(kTRUE);

 this->DefineHistograms();

  if(fListFBHijing){
    fTaskConfigParm->SetBinContent(11,1.0);
  }
  if(bApplyNUACorr){
    fTaskConfigParm->SetBinContent(12,1.0);
  }
  if(bApplyV0MCorr){
    fTaskConfigParm->SetBinContent(13,1.0);
  }
  if(bApplyZDCCorr){
    fTaskConfigParm->SetBinContent(14,1.0);
  }

 fTaskConfigParm->SetBinContent(15,fHarmonicN);
 fTaskConfigParm->SetBinContent(16,fHarmonicM);
 fTaskConfigParm->SetBinContent(17,fHarmonicPsi);

 PostData(1,fListHistos); 
 PostData(2,fListCalibs); 
 PostData(3,fListNUAHist);

 printf("\n ======================== Info::UserCreateOutPutObject Called ========================  \n");
 std::cout<<" sDataSet = "<<sDataSet.Data()<<" FB = "<<gFilterBit<<" PileUp C = "<<fPileUpConstParm<<" HBT cut = "<<fHBTCutValue<<"\n"<<std::endl;


}

AliAnalysisTaskCMEV0::~AliAnalysisTaskCMEV0()
{
  delete                     fListHistos;        
  //delete                   fListCalibs;    
  //delete                  fListNUAHist;     

  if(fAnalysisUtil) delete fAnalysisUtil; // it is '= new' !!!

  if(fHCentBinTrkRecenter)
    delete fHCentBinTrkRecenter;

  if(fHCorrectZDNP)   delete fHCorrectZDNP;
  if(fHCorrectV0M)    delete fHCorrectV0M;
  if(fHAvgerageQnV0A) delete fHAvgerageQnV0A;
  if(fHAvgerageQnV0C) delete fHAvgerageQnV0C;


  for(int i=0;i<5;i++){
    if(fHCorrectNUApos[i]) delete fHCorrectNUApos[i];
    if(fHCorrectNUAneg[i]) delete fHCorrectNUAneg[i];
  }
  for(int i=0;i<10;i++){
    if(fFB_Efficiency_Cent[i]) delete fFB_Efficiency_Cent[i];
    if(fFB_Efficiency_Pos[i])  delete fFB_Efficiency_Pos[i];
    if(fFB_Efficiency_Neg[i])  delete fFB_Efficiency_Neg[i];
  }
}







void AliAnalysisTaskCMEV0::UserExec(Option_t *)
{
//printf("\n ======================== Debug::UserExec Called ========================  \n");

 Float_t stepCount = 0.5;
 fHist_Event_count->Fill(stepCount); //1
 stepCount++;

 AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(1));

 if(!aod || !fEvent){
   printf("\n ... ::UserExec = no AOD or Flow Event, \n.... EXIT .....  \n");
   return;
 }

 fHist_Event_count->Fill(stepCount); //2
 stepCount++;


 //---------pileup rejection: --------

 Bool_t kPileupEvent = kFALSE;

 kPileupEvent = CheckEventIsPileUp(aod);

 if(kPileupEvent)    return;

 

 fHist_Event_count->Fill(stepCount); //3
 stepCount++;



 Double_t centrV0M = -1;
 Double_t centrCL1 = -1;
 Double_t centrCL0 = -1;
 Double_t centrTRK = -1;

 EvtCent = -1;

 if(sDataSet=="2010"||sDataSet=="2011"){
    centrV0M = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
    centrCL1 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1");
    centrCL0 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL0");
    centrTRK = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("TRK");
 }
 else{
  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
   if(!fMultSelection) {
     printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
     exit(1);
   }
    centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
    centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
    centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
    centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");
 }



 if(sCentEstimator=="V0" || sCentEstimator=="V0M"){
   EvtCent = centrV0M;
 }
 else{
   EvtCent = centrCL1;
 }


 Int_t iCentSPD = centrCL1;

 if(iCentSPD > 90) return;

 fHist_Event_count->Fill(stepCount); //4
 stepCount++;




 Int_t cIndex = 0;
 Int_t cForNUA = 0;
 Int_t cInputNUA = 0;

 if(EvtCent<5.0) {
   cIndex  = 0; 
   cForNUA = 0; //0=0-5,
   cInputNUA = 0;
 }
 else if(EvtCent>=5.0 && EvtCent<10){
   cIndex  = 1;
   cForNUA = 1; // 1=5-10,
   cInputNUA = 1;
 }
 else if(EvtCent>=10.0) {
   cIndex = abs(EvtCent/10.0)  +  1;
   cForNUA = cIndex;  //2 = 10-20,

   if(EvtCent>=10 && EvtCent<40.0)
     cInputNUA = 2;

   if(EvtCent>=20 && EvtCent<40.0){ 
     cForNUA = 3; // 3=20-40
   }
   if(EvtCent>=40.0){
     cForNUA = 4; // 4=40-90
     cInputNUA = 3; 
   }
 }



  //------- load v0 calib file of Alex --------
 Int_t runindex = -111;
 Int_t runNumber = aod->GetRunNumber();
 Float_t psiN = fHarmonicPsi;
 runindex = GetCurrentRunIndex(runNumber);

 if(runNumber!=fOldRunNum){ 

  
   if(sDataSet=="2015" || sDataSet=="2015HI"){
     OpenInfoCalbration(runNumber, psiN);
   }
   if(bApplyV0MCorr) {
     GetV0MCorrectionHist(runNumber);
   }
   if(bApplyNUACorr){
     GetNUACorrectionHist(runNumber,sFileNUA);
   }
   if(bApplyZDCCorr && fListZDNCorr){
     GetZDCCorrectionHist(runNumber);
   }

   fOldRunNum = runNumber;
 }




 Double_t fMagField = aod->GetMagneticField();

 Int_t QAindex = -1;

 //Femto:
 Int_t bSign = (fMagField > 0) ? 1 : -1;

 //Double_t fHBTCutValue = 0.02;


 if(fMagField<0)
   QAindex = 0;
 else if(fMagField>0)
   QAindex = 1;

//-------- V0 info ---------------
 const AliAODVZERO *fAODV0 = aod->GetVZEROData();


 Double_t QyanCor = 0., QycnCor = 0.;
 Double_t QxanCor = 0., QxcnCor = 0.;
 Double_t Qxan2  = 0., Qyan2  = 0.;
 Double_t Qxcn2  = 0., Qycn2  = 0.;
 Double_t sumMa = 0., sumMc = 0.;
 Double_t fMultv0 = 0.;
 Double_t phiV0 = 0.;


 if(sDataSet=="2015" || sDataSet=="2015HI") {
   GetV0QvectAndMult(fAODV0, psiN, Qxan2, Qyan2, sumMa, Qxcn2, Qycn2, sumMc);

   QyanCor = (Qyan2 - fQynmV0A->GetBinContent(iCentSPD+1))/fQynsV0A->GetBinContent(iCentSPD+1);
   QycnCor = (Qycn2 - fQynmV0C->GetBinContent(iCentSPD+1))/fQynsV0C->GetBinContent(iCentSPD+1);

   QxanCor = Qxan2;
   QxcnCor = Qxcn2;

   if(psiN != 4.){
     QxanCor = (Qxan2 - fQxnmV0A->GetBinContent(iCentSPD+1))/fQxnsV0A->GetBinContent(iCentSPD+1);
     QxcnCor = (Qxcn2 - fQxnmV0C->GetBinContent(iCentSPD+1))/fQxnsV0C->GetBinContent(iCentSPD+1);
   }
 }
 //else{

   Double_t Qxan3  = 0., Qyan3  = 0.;
   Double_t Qxcn3  = 0., Qycn3  = 0.;
   Qxan2  = 0., Qyan2 = 0.;
   Qxcn2  = 0., Qycn2 = 0.;

   for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA

     fMultv0 = fAODV0->GetMultiplicity(iV0);
 
     if(fHCorrectV0M){
       fMultv0 = fMultv0 * fHCorrectV0M->GetBinContent(iV0+1);  // Gain Correction
     }

     fV0MultChVsRun->Fill(iV0+0.5,runindex,fMultv0);

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
   }//channel loop

   //For V0-Qn Recenter
   fV0CQ2xVsCentRun->Fill(centrCL1,runindex,Qxcn2);
   fV0CQ2yVsCentRun->Fill(centrCL1,runindex,Qycn2);
   fV0AQ2xVsCentRun->Fill(centrCL1,runindex,Qxan2); 
   fV0AQ2yVsCentRun->Fill(centrCL1,runindex,Qyan2); 

   fV0CQ3xVsCentRun->Fill(centrCL1,runindex,Qxcn3);
   fV0CQ3yVsCentRun->Fill(centrCL1,runindex,Qycn3);
   fV0AQ3xVsCentRun->Fill(centrCL1,runindex,Qxan3); 
   fV0AQ3yVsCentRun->Fill(centrCL1,runindex,Qyan3); 

   Int_t icentV0Qn = centrCL1;
   icentV0Qn += 1;

   if(sDataSet=="2015pPb" || sDataSet=="pPb" || sDataSet=="2015LI") {
     if(psiN==3.0){
       QxanCor = Qxan3;  //3rd
       QxcnCor = Qxcn3;
       QyanCor = Qyan3;
       QycnCor = Qycn3;

       if(fHAvgerageQnV0C && fHAvgerageQnV0A && icentV0Qn < 91){
         QxanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,3);  //x = Cos
         QxcnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,3);  //x = Cos
         QyanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,4);  //y = Sin
         QycnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,4);  //y = Sin
       }
       //printf("\n .... I am using my own V0 gain correction for Psi3...\n");
     }
     else{
       QxanCor = Qxan2;  //2nd or higher
       QxcnCor = Qxcn2;
       QyanCor = Qyan2;
       QycnCor = Qycn2;
       if(fHAvgerageQnV0C && fHAvgerageQnV0A && icentV0Qn < 91){
         QxanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,1);  //x = Cos
         QxcnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,1);  //x = Cos
         QyanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,2);  //y = Sin
         QycnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,2);  //y = Sin  
       }
       //printf("\n .... I am using my own V0 gain correction for Psi2...\n");
     }
   }

 //} else (other than 2015o_pass1 HI data)



 //Fill V0M Mult for all System: filled above
 /*for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA
   fMultv0 = fAODV0->GetMultiplicity(iV0);
   fV0MultChVsRun->Fill(iV0+0.5,runindex,fMultv0);
   }*/


 


 if(sumMa <= 0 || sumMc <= 0)   return;
 fHist_Event_count->Fill(stepCount); //5
 stepCount++;


 AliAODVertex *pVertex = aod->GetPrimaryVertex();
 Double_t VtxZ    =   pVertex->GetZ();
 Double_t VtxX    =   pVertex->GetX();
 Double_t VtxY    =   pVertex->GetY();

 fVzDistribuion->Fill(VtxZ);

 if(sDataSet=="2015pPb" || sDataSet=="pPb"){

  if(VtxZ < 1.e-4 && VtxX < 1.e-4 && VtxY < 1.e-4)
    return;

  fHistVtxZvsRun->Fill(runindex,VtxZ);
  fHistVtxXvsRun->Fill(runindex,VtxX);
  fHistVtxYvsRun->Fill(runindex,VtxY);
 }

 
 //------------- Load ZDC info: -----------
 Double_t energyZNC=0.,energyZNA=0.,energyZPC=0.,energyZPA=0.;

 energyZNC = ((AliVAODHeader*)aod->GetHeader())->GetZDCN1Energy();
 energyZNA = ((AliVAODHeader*)aod->GetHeader())->GetZDCN2Energy();
 energyZPC = ((AliVAODHeader*)aod->GetHeader())->GetZDCP1Energy();
 energyZPA = ((AliVAODHeader*)aod->GetHeader())->GetZDCP2Energy();


 //Apply Gain Correction:
 if(fHCorrectZDNP){
   energyZNC = energyZNC * fHCorrectZDNP->GetBinContent((abs(EvtCent)+1),1);
   energyZNA = energyZNA * fHCorrectZDNP->GetBinContent((abs(EvtCent)+1),2);
   energyZPC = energyZPC * fHCorrectZDNP->GetBinContent((abs(EvtCent)+1),3);
   energyZPA = energyZPA * fHCorrectZDNP->GetBinContent((abs(EvtCent)+1),4);
 }

 //fill ZDC info: Uses 
 if(bFillZDCinfo){

   fHEnergyZNCvsCentRun->Fill(EvtCent,runindex,energyZNC);
   fHEnergyZNAvsCentRun->Fill(EvtCent,runindex,energyZNA);
   fHEnergyZPCvsCentRun->Fill(EvtCent,runindex,energyZPC);
   fHEnergyZPAvsCentRun->Fill(EvtCent,runindex,energyZPA);

   //if(EvtCent>30. && EvtCent<40.){
   //fHEnergyZPCvsZPA->Fill(energyZPC,energyZPA);
   //fHEnergyZNCvsZNA->Fill(energyZNC,energyZNA);
   //}
 }





 Int_t BadRuns[10] = {246871,246870,246867,246865,246864,246859,246858,246676,246675,246540};
 Int_t iBadrun = 10;
 Int_t isBadRun = 0;

 for(int ib=0;ib<iBadrun;ib++){
   if(runNumber==BadRuns[ib]){
     isBadRun = 1;
     break;
   }
 }

 

 //Store Qvector for tracks:
 Double_t QxPos[2] = {0.,};  //[0]=n*phi,[1]=2*n*phi
 Double_t QyPos[2] = {0.,};
 Double_t QxNeg[2] = {0.,};  
 Double_t QyNeg[2] = {0.,};

 //Store Qvector for Auto-correlations
 Double_t QxAutoPos[2] = {0.,};  
 Double_t QyAutoPos[2] = {0.,};
 Double_t QxAutoNeg[2] = {0.,};  
 Double_t QyAutoNeg[2] = {0.,};

 
 //temp: for QA only
 Double_t QxPosQAEta[80] = {0.,};  
 Double_t QyPosQAEta[80] = {0.,};
 Double_t QxNegQAEta[80] = {0.,};  
 Double_t QyNegQAEta[80] = {0.,};
 Double_t QxAutoPosQAEta[80] = {0.,};  
 Double_t QyAutoPosQAEta[80] = {0.,};
 Double_t QxAutoNegQAEta[80] = {0.,};  
 Double_t QyAutoNegQAEta[80] = {0.,};
 Double_t MPOIposQAEta[80]   = {0.,};
 Double_t MPOInegQAEta[80]   = {0.,};
 Int_t    iEtaQA = -1; 

 Double_t AvgCos1n=0.,AvgSin1n=0.;
 //Double_t AvgCos2n=0.,AvgSin2n=0.;
 Double_t AvgDWCos2n=0.,AvgDWSin2n=0.;


 //Multiplicity of POIs:
 Double_t MPOIpos = 0;
 Double_t MPOIneg = 0;
 Double_t McorrPos = 0.;
 Double_t McorrNeg = 0.;

 //TPC Event plane Qvectors:
 Double_t QxTPC[3] = {0.,};
 Double_t QyTPC[3] = {0.,};

 Int_t ptBin;
 Int_t nRefMult = 0;
 Int_t n = fHarmonicN; 
 Int_t m = fHarmonicM; 
 Int_t p = (fHarmonicN+fHarmonicM);  // cos(nphi1 + mphi2 - p*Psi_{2,3}) //not 4 

 Double_t nRefMultWgt = 0.;
 Int_t skipPairHBT = 0;




 //Int_t iCentRec = 0;  //centrality bin for track recenter 
 //iCentRec = fHCentBinTrkRecenter->FindBin(EvtCent);

 //Int_t nBinVz=1,nBinPhi=1,nBinEta=1;

 //CME using Event plane:
 Double_t pi =  TMath::Pi();
 Double_t Psi2V0C = 1./psiN*TMath::ATan2(QycnCor,QxcnCor);
 if(Psi2V0C<0.) Psi2V0C += 2*pi/psiN;

 Double_t Psi2V0A = 1./psiN*TMath::ATan2(QyanCor,QxanCor);
 if(Psi2V0A<0.) Psi2V0A += 2*pi/psiN;

 
 /*
 fV0AQ2xVsCentRun->Fill(centrCL1,runindex,QycnCor); // For V0-Q Recenter QA.
 fV0AQ2yVsCentRun->Fill(centrCL1,runindex,QyanCor); 
 fV0CQ2xVsCentRun->Fill(centrCL1,runindex,QxcnCor);
 fV0CQ2yVsCentRun->Fill(centrCL1,runindex,QycnCor);*/



 Int_t iBinNUA = 1;
 Int_t iVzNUA = -1;
 iVzNUA = fVzBinFinderForNUA->FindBin(VtxZ) - 1;



 Int_t iTracks = fEvent->NumberOfTracks();

 Int_t icReject = 0;


 AliFlowTrackSimple*   pTrack1 = NULL;
 AliFlowTrackSimple*   pTrack2 = NULL;

 Double_t dPhi1,dPt1,dEta1,dChrg1;
 Double_t dPhi2,dPt2,dEta2,dChrg2;
 Double_t WgtEP = 1.0;

 Double_t ptw1  = 1.0, w1NUA = 1.0;
 Double_t ptw2  = 1.0; 
 Double_t w2NUA = 1.0;

 //Femto:
 Double_t deltaPhi=0.,deltaEta=0.;
 


 for(int i=0; i<iTracks; i++) {
   pTrack1    =     fEvent->GetTrack(i);
   if(!pTrack1)                continue;
   dPhi1      =          pTrack1->Phi(); //0-2pi range
   dPt1       =          pTrack1-> Pt();
   dEta1      =          pTrack1->Eta();
   dChrg1     =       pTrack1->Charge();

   if(!pTrack1->InPOISelection())
   continue;

   //remove randomly -Ve tracks:
   //if(bRemNegTrkRndm && dChrg1<0 && fRand.Rndm()<0.15){ 
   //icReject++;
   //continue;
   //}


   if(bFillEtaPhiNUA){
     if(dChrg1>0){
       fHist3DEtaPhiVz_Pos_Run[cForNUA][runindex]->Fill(VtxZ,dPhi1,dEta1);
     }
     else if(dChrg1<0){
       fHist3DEtaPhiVz_Neg_Run[cForNUA][runindex]->Fill(VtxZ,dPhi1,dEta1);
     }
   }
   /*
   if(bFillEtaPhiNUA){
     if(dChrg1>0){
       fHist3DEtaPhiVz_Pos_Run[iVzNUA][runindex]->Fill(dPt1,dPhi1,dEta1);
     }
     else if(dChrg1<0){
       fHist3DEtaPhiVz_Neg_Run[iVzNUA][runindex]->Fill(dPt1,dPhi1,dEta1);
     }
   }*/

   if(bApplyNUACorr){
     //get NUA weights: 
     if(sFileNUA=="New2018"){
       if(dChrg1>0){
         iBinNUA = fHCorrectNUApos[cForNUA]->FindBin(VtxZ,dPhi1,dEta1);
         w1NUA = fHCorrectNUApos[cForNUA]->GetBinContent(iBinNUA);     
       }
       else if(dChrg1<0){
         iBinNUA = fHCorrectNUAneg[cForNUA]->FindBin(VtxZ,dPhi1,dEta1);
         w1NUA = fHCorrectNUAneg[cForNUA]->GetBinContent(iBinNUA);     
       }
     }
     else if(sFileNUA=="OldJ"||sFileNUA=="NewR"){
       if(dChrg1>0){
         iBinNUA = fHCorrectNUApos[cInputNUA]->FindBin(VtxZ,dPhi1,dEta1);
         w1NUA = fHCorrectNUApos[cInputNUA]->GetBinContent(iBinNUA);     
       }
       else if(dChrg1<0){
         iBinNUA = fHCorrectNUAneg[cInputNUA]->FindBin(VtxZ,dPhi1,dEta1);
         w1NUA = fHCorrectNUAneg[cInputNUA]->GetBinContent(iBinNUA);     
       }
       if(sFileNUA=="OldJ") w1NUA = 1./w1NUA;
     }
     else if(sFileNUA=="NewPt"||sFileNUA=="NewpT"){
       if(dChrg1>0){
         iBinNUA = fHCorrectNUApos[iVzNUA]->FindBin(dPt1,dPhi1,dEta1);
         w1NUA = fHCorrectNUApos[iVzNUA]->GetBinContent(iBinNUA);     
       }
       else if(dChrg1<0){
         iBinNUA = fHCorrectNUAneg[iVzNUA]->FindBin(dPt1,dPhi1,dEta1);
         w1NUA = fHCorrectNUAneg[iVzNUA]->GetBinContent(iBinNUA);     
       }
     }
     if(w1NUA > 1e3 ) w1NUA = 1.0;
   }



   if(dChrg1>0){
     fHistChPosvsEtaPtRun[cIndex]->Fill(dPt1,dEta1,runindex,w1NUA);
   }
   else if(dChrg1<0){
     fHistChNegvsEtaPtRun[cIndex]->Fill(dPt1,dEta1,runindex,w1NUA);
   }



   if(sMCdimension=="1D"){
     ptBin = fFB_Efficiency_Cent[cIndex]->FindBin(dPt1);
     ptw1   = 1.0/fFB_Efficiency_Cent[cIndex]->GetBinContent(ptBin);
   }
   else if(sMCdimension=="3D"){
     if(dChrg1>0){
       ptBin = fFB_Efficiency_Pos[cIndex]->FindBin(dEta1,dPt1,dPhi1); //left here
       ptw1  = fFB_Efficiency_Pos[cIndex]->GetBinContent(ptBin); 
     }
     else if(dChrg1<0){
       ptBin = fFB_Efficiency_Neg[cIndex]->FindBin(dEta1,dPt1,dPhi1); //left here
       ptw1  = fFB_Efficiency_Neg[cIndex]->GetBinContent(ptBin); 
     }
   }

   if (ptw1>1e3) ptw1 = 1.0;

   
   //if(bFillAvgTPCQn) { //fill TPC Q-vect: 
     if(VtxZ>0 && dEta1>=0){
       if(dChrg1>0){
	 //fHCos1nPosChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Cos(dPhi1));
	 //fHSin1nPosChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Sin(dPhi1));
	 fHCos2nPosChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Cos(2.*dPhi1));
	 fHSin2nPosChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Sin(2.*dPhi1));
	 fHCos3nPosChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Cos(3.*dPhi1));
	 fHSin3nPosChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Sin(3.*dPhi1));
	 //fHCos4nPosChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Cos(4.*dPhi1));
	 //fHSin4nPosChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Sin(4.*dPhi1));
	 //fHCos2nDWPosChEtaVz[0]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Cos(2.*dPhi1));
	 //fHSin2nDWPosChEtaVz[0]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Sin(2.*dPhi1));
       }
       else if(dChrg1<0){
	 //fHCos1nNegChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Cos(dPhi1));
	 //fHSin1nNegChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Sin(dPhi1));
	 fHCos2nNegChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Cos(2.*dPhi1));
	 fHSin2nNegChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Sin(2.*dPhi1));
	 fHCos3nNegChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Cos(3.*dPhi1));
	 fHSin3nNegChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Sin(3.*dPhi1));
	 //fHCos4nNegChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Cos(4.*dPhi1));
	 //fHSin4nNegChEtaVz[0]->Fill(EvtCent,runindex,ptw1*TMath::Sin(4.*dPhi1));
	 //fHCos2nDWNegChEtaVz[0]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Cos(2.*dPhi1));
	 //fHSin2nDWNegChEtaVz[0]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Sin(2.*dPhi1));
       }
     }
     else if(VtxZ>0 && dEta1<0){
       if(dChrg1>0){
	 //fHCos1nPosChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Cos(dPhi1));
	 //fHSin1nPosChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Sin(dPhi1));
	 fHCos2nPosChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Cos(2.*dPhi1));
	 fHSin2nPosChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Sin(2.*dPhi1));
	 fHCos3nPosChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Cos(3.*dPhi1));
	 fHSin3nPosChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Sin(3.*dPhi1));
	 //fHCos4nPosChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Cos(4.*dPhi1));
	 //fHSin4nPosChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Sin(4.*dPhi1));
	 //fHCos2nDWPosChEtaVz[1]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Cos(2.*dPhi1));
	 //fHSin2nDWPosChEtaVz[1]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Sin(2.*dPhi1));
       }
       else if(dChrg1<0){
	 // fHCos1nNegChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Cos(dPhi1));
	 // fHSin1nNegChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Sin(dPhi1));
	 fHCos2nNegChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Cos(2.*dPhi1));
	 fHSin2nNegChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Sin(2.*dPhi1));
	 fHCos3nNegChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Cos(3.*dPhi1));
	 fHSin3nNegChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Sin(3.*dPhi1));
	 // fHCos4nNegChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Cos(4.*dPhi1));
	 // fHSin4nNegChEtaVz[1]->Fill(EvtCent,runindex,ptw1*TMath::Sin(4.*dPhi1));
	 // fHCos2nDWNegChEtaVz[1]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Cos(2.*dPhi1));
	 // fHSin2nDWNegChEtaVz[1]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Sin(2.*dPhi1));
       }
     }
     else if(VtxZ<0 && dEta1>=0){
        if(dChrg1>0){
         // fHCos1nPosChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Cos(dPhi1));
         // fHSin1nPosChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Sin(dPhi1));
	 fHCos2nPosChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Cos(2.*dPhi1));
	 fHSin2nPosChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Sin(2.*dPhi1));
	 fHCos3nPosChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Cos(3.*dPhi1));
	 fHSin3nPosChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Sin(3.*dPhi1));
	 // fHCos4nPosChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Cos(4.*dPhi1));
	 // fHSin4nPosChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Sin(4.*dPhi1));
	 // fHCos2nDWPosChEtaVz[2]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Cos(2.*dPhi1));
	 // fHSin2nDWPosChEtaVz[2]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Sin(2.*dPhi1));
       }
       else if(dChrg1<0){
	 // fHCos1nNegChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Cos(dPhi1));
	 // fHSin1nNegChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Sin(dPhi1));
	 fHCos2nNegChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Cos(2.*dPhi1));
	 fHSin2nNegChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Sin(2.*dPhi1));
	 fHCos3nNegChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Cos(3.*dPhi1));
	 fHSin3nNegChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Sin(3.*dPhi1));
	 // fHCos4nNegChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Cos(4.*dPhi1));
	 // fHSin4nNegChEtaVz[2]->Fill(EvtCent,runindex,ptw1*TMath::Sin(4.*dPhi1));
	 // fHCos2nDWNegChEtaVz[2]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Cos(2.*dPhi1));
	 // fHSin2nDWNegChEtaVz[2]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Sin(2.*dPhi1));
       }

     }
     else if(VtxZ<0 && dEta1<0){
       if(dChrg1>0){
	 // fHCos1nPosChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Cos(dPhi1));
	 // fHSin1nPosChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Sin(dPhi1));
	 fHCos2nPosChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Cos(2.*dPhi1));
	 fHSin2nPosChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Sin(2.*dPhi1));
	 fHCos3nPosChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Cos(3.*dPhi1));
	 fHSin3nPosChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Sin(3.*dPhi1));
	 // fHCos4nPosChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Cos(4.*dPhi1));
	 // fHSin4nPosChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Sin(4.*dPhi1));
	 // fHCos2nDWPosChEtaVz[3]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Cos(2.*dPhi1));
	 // fHSin2nDWPosChEtaVz[3]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Sin(2.*dPhi1));
       }
       else if(dChrg1<0){
	 // fHCos1nNegChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Cos(dPhi1));
	 // fHSin1nNegChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Sin(dPhi1));
	 fHCos2nNegChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Cos(2.*dPhi1));
	 fHSin2nNegChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Sin(2.*dPhi1));
	 fHCos3nNegChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Cos(3.*dPhi1));
	 fHSin3nNegChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Sin(3.*dPhi1));
	 // fHCos4nNegChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Cos(4.*dPhi1));
	 // fHSin4nNegChEtaVz[3]->Fill(EvtCent,runindex,ptw1*TMath::Sin(4.*dPhi1));
	 // fHCos2nDWNegChEtaVz[3]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Cos(2.*dPhi1));
	 // fHSin2nDWNegChEtaVz[3]->Fill(EvtCent,runindex,ptw1*ptw1*TMath::Sin(2.*dPhi1));
       }
     }
     //} //if(fillFlag==kTRUE)






   iEtaQA = fEtaBinFinderForQA->FindBin(dEta1) - 1;

   if(dChrg1>0){
     QxPos[0] += (ptw1*w1NUA*TMath::Cos(n*dPhi1) - AvgCos1n);
     QyPos[0] += (ptw1*w1NUA*TMath::Sin(n*dPhi1) - AvgSin1n);
     QxAutoPos[0] += (ptw1*ptw1*w1NUA*w1NUA*TMath::Cos(2.*n*dPhi1) - AvgDWCos2n);
     QyAutoPos[0] += (ptw1*ptw1*w1NUA*w1NUA*TMath::Sin(2.*n*dPhi1) - AvgDWSin2n);
     MPOIpos  += w1NUA*ptw1;
     McorrPos += w1NUA*ptw1*w1NUA*ptw1;
     //QA: eta dependence
     if(EvtCent>=10 && EvtCent<20){
       QxPosQAEta[iEtaQA]     += (ptw1*w1NUA*TMath::Cos(n*dPhi1) - AvgCos1n);
       QyPosQAEta[iEtaQA]     += (ptw1*w1NUA*TMath::Sin(n*dPhi1) - AvgSin1n);
       QxAutoPosQAEta[iEtaQA] += (ptw1*ptw1*w1NUA*w1NUA*TMath::Cos(2.*n*dPhi1) - AvgDWCos2n);
       QyAutoPosQAEta[iEtaQA] += (ptw1*ptw1*w1NUA*w1NUA*TMath::Sin(2.*n*dPhi1) - AvgDWSin2n);
       MPOIposQAEta[iEtaQA]   +=  ptw1*w1NUA;
     }
   }
   else if(dChrg1<0){
     QxNeg[0] += (ptw1*w1NUA*TMath::Cos(n*dPhi1) - AvgCos1n);
     QyNeg[0] += (ptw1*w1NUA*TMath::Sin(n*dPhi1) - AvgSin1n);
     QxAutoNeg[0] += (ptw1*ptw1*w1NUA*w1NUA*TMath::Cos(2.*n*dPhi1) - AvgDWCos2n);
     QyAutoNeg[0] += (ptw1*ptw1*w1NUA*w1NUA*TMath::Sin(2.*n*dPhi1) - AvgDWSin2n);
     MPOIneg  += w1NUA*ptw1;
     McorrNeg += w1NUA*ptw1*w1NUA*ptw1;
     //QA: eta dependence
     if(EvtCent>=10 && EvtCent<20){
       QxNegQAEta[iEtaQA]     += (ptw1*w1NUA*TMath::Cos(n*dPhi1) - AvgCos1n);
       QyNegQAEta[iEtaQA]     += (ptw1*w1NUA*TMath::Sin(n*dPhi1) - AvgSin1n);
       QxAutoNegQAEta[iEtaQA] += (ptw1*ptw1*w1NUA*w1NUA*TMath::Cos(2.*n*dPhi1) - AvgDWCos2n);
       QyAutoNegQAEta[iEtaQA] += (ptw1*ptw1*w1NUA*w1NUA*TMath::Sin(2.*n*dPhi1) - AvgDWSin2n);
       MPOInegQAEta[iEtaQA]   +=  ptw1*w1NUA;
     }
   }

 //QxTPC[0] += ptw1*w1NUA*TMath::Cos(n*dPhi1);
 //QyTPC[0] += ptw1*w1NUA*TMath::Sin(n*dPhi1);
   QxTPC[1] += ptw1*w1NUA*TMath::Cos(2*dPhi1);
   QyTPC[1] += ptw1*w1NUA*TMath::Sin(2*dPhi1);

   QxTPC[2] += ptw1*w1NUA*TMath::Cos(3*dPhi1);
   QyTPC[2] += ptw1*w1NUA*TMath::Sin(3*dPhi1);

   nRefMult++;
   nRefMultWgt += ptw1*w1NUA;
  

   //if(i%10==0)
   //cout<<" track "<<i<<" eta = "<<dEta1<<"\tdPhi1 = "<<dPhi1<<"\tptw1 = "<<ptw1<<"\tw1NUA = "<<w1NUA<<"\tAvgCos1n = "<<AvgCos1n<<endl;

   if(bSkipNestedTrk) continue;

   //2nd track loop:
   for(int j=0; j<iTracks; j++) {

     if(j==i) continue;  //Auto-correlation removed.

     pTrack2    =     fEvent->GetTrack(j);
     if(!pTrack2)                continue;
     dPhi2      =          pTrack2->Phi(); //0-2pi range
     dPt2       =          pTrack2-> Pt();
     dEta2      =          pTrack2->Eta();
     dChrg2     =       pTrack2->Charge();

     if(!pTrack2->InPOISelection())
     continue;

     //remove randomly -Ve tracks:
     //if(bRemNegTrkRndm && dChrg1<0 && fRand.Rndm()<0.15){ 
     //continue;
     //}

   
     if(sMCdimension=="1D"){
       ptBin = fFB_Efficiency_Cent[cIndex]->FindBin(dPt2);
       ptw2   = 1.0/fFB_Efficiency_Cent[cIndex]->GetBinContent(ptBin);
     }
     else if(sMCdimension=="3D"){
       if(dChrg2>0){
         ptBin = fFB_Efficiency_Pos[cIndex]->FindBin(dEta2,dPt2,dPhi2); //left here
         ptw2  = fFB_Efficiency_Pos[cIndex]->GetBinContent(ptBin); 
       }
       else if(dChrg2<0){
         ptBin = fFB_Efficiency_Neg[cIndex]->FindBin(dEta2,dPt2,dPhi2); //left here
         ptw2  = fFB_Efficiency_Neg[cIndex]->GetBinContent(ptBin); 
       }
     }

     if (ptw2>1e3) ptw1 = 1.0;


     if(bApplyNUACorr) {
     //get NUA weights: 
      if(sFileNUA=="New2018"){
        if(dChrg2>0){
          iBinNUA = fHCorrectNUApos[cForNUA]->FindBin(VtxZ,dPhi1,dEta1);
          w1NUA = fHCorrectNUApos[cForNUA]->GetBinContent(iBinNUA);     
        }
        else if(dChrg2<0){
          iBinNUA = fHCorrectNUAneg[cForNUA]->FindBin(VtxZ,dPhi1,dEta1);
          w1NUA = fHCorrectNUAneg[cForNUA]->GetBinContent(iBinNUA);     
        }
      }
       else if(sFileNUA=="OldJ"||sFileNUA=="NewR"){
         if(dChrg2>0){
           iBinNUA = fHCorrectNUApos[cInputNUA]->FindBin(VtxZ,dPhi2,dEta2);
           w2NUA = fHCorrectNUApos[cInputNUA]->GetBinContent(iBinNUA);     
         }
         else if(dChrg2<0){
           iBinNUA = fHCorrectNUAneg[cInputNUA]->FindBin(VtxZ,dPhi2,dEta2);
           w2NUA = fHCorrectNUAneg[cInputNUA]->GetBinContent(iBinNUA);     
         }
         if(sFileNUA=="OldJ") w2NUA = 1./w2NUA;
       }
       else if(sFileNUA=="NewPt"||sFileNUA=="NewpT"){
         if(dChrg2>0){
           iBinNUA = fHCorrectNUApos[iVzNUA]->FindBin(dPt2,dPhi2,dEta2);
           w2NUA = fHCorrectNUApos[iVzNUA]->GetBinContent(iBinNUA);     
         }
         else if(dChrg2<0){
           iBinNUA = fHCorrectNUAneg[iVzNUA]->FindBin(dPt2,dPhi2,dEta2);
           w2NUA = fHCorrectNUAneg[iVzNUA]->GetBinContent(iBinNUA);     
         }
       }
       if(w2NUA > 1e3 ) w2NUA = 1.0;
     }
     if(bApplyNUAforEP){
       WgtEP = ptw1*ptw2*w1NUA*w2NUA;
     }
     else{
       WgtEP = 1.0;
     } 


     //------------- test    ------ before femto cuts:
     if(dChrg1!=dChrg2){
       fHist_Corr3p_SP_Norm_PN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       fHist_Corr3p_SP_Norm_PN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
     }
     else if(dChrg1>0 && dChrg2>0){
       fHist_Corr3p_SP_Norm_PP[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       fHist_Corr3p_SP_Norm_PP[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
     }
     else if(dChrg1<0 && dChrg2<0){
       fHist_Corr3p_SP_Norm_NN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       fHist_Corr3p_SP_Norm_NN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
     }

     if(dChrg1==dChrg2) {
       fdPhiFemtoCut->Fill(2.5,1e-8);
     }
     else{
       fdPhiFemtoCut->Fill(0.5,1e-8);
     }


     //---------- Femto Cuts ----------
     deltaEta = dEta1 - dEta2;
     deltaPhi = dPhi1 - dPhi2;
     skipPairHBT = 0;

     if(TMath::Abs(deltaEta) < fHBTCutValue * 2.5 * 3) //fHBTCutValue = 0.02 [default for dphicorrelations]
     {	
	 Float_t phi1rad = dPhi1;   // phi in rad
	 Float_t phi2rad = dPhi2;

	 // check first boundaries to see if is worth to loop and find the minimum
	 Float_t dphistar1 = GetDPhiStar(phi1rad, dPt1, dChrg1, phi2rad, dPt2, dChrg2, 0.8, bSign);
	 Float_t dphistar2 = GetDPhiStar(phi1rad, dPt1, dChrg1, phi2rad, dPt2, dChrg2, 2.5, bSign);

	 const Float_t kLimit = fHBTCutValue * 3;

	 Float_t dphistarminabs = 1e5;
	 Float_t dphistarmin = 1e5;

	 if(TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0 ) {

	   for(Double_t rad=0.8; rad<2.51; rad+=0.01) {
	     Float_t dphistar = GetDPhiStar(phi1rad, dPt1, dChrg1, phi2rad, dPt2, dChrg2, rad, bSign);
	     Float_t dphistarabs = TMath::Abs(dphistar);

	     if(dphistarabs < dphistarminabs) {
	       dphistarmin = dphistar;
	       dphistarminabs = dphistarabs;
	     }
	   }
	   //Float_t dPtDiff = fabs(dPt1-dPt2);
	   if(dphistarminabs < fHBTCutValue && TMath::Abs(deltaEta) < fHBTCutValue) {
	     //cout<<"dEta "<<deltaEta<<"  dphi "<<deltaPhi<<" ch1 = "<<dChrg1<<" ch2 = "<<dChrg2<<" fHBTCutValue = "<<fHBTCutValue<<endl;
	     if(dChrg1==dChrg2) {
	       fdPhiFemtoCut->Fill(3.5,1e-8);
	     }
	     else{
	       fdPhiFemtoCut->Fill(1.5,1e-8);
	     }
	     skipPairHBT = 1;
	     //continue;
	   }   
	 }
       }

   //---------------- Femto Cuts done --------------------












     if(dChrg1!=dChrg2) {
       //if(skipPairHBT > 0)
       //cout<<"All dEta "<<deltaEta<<"  dphi "<<deltaPhi<<" ch1 = "<<dChrg1<<" ch2 = "<<dChrg2<<" skipPairHBT = "<<skipPairHBT<<endl;
	
       fdPhiFemtoCut->Fill(1.5,1e-8);
       
       fHist_Corr3p_EP_Norm_PN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       fHist_Corr3p_EP_Norm_PN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
       //vs Refmult:
       fHist_Corr3p_EP_Refm_PN[QAindex][0]->Fill(fRefMultCorr, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       fHist_Corr3p_EP_Refm_PN[QAindex][1]->Fill(fRefMultCorr, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);

       if(EvtCent>=10. && EvtCent<20){
         fCent3pvsVzRun->Fill(VtxZ,runindex,TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A));
       }

       //fHist_Corr3p_vsRun_EP_PN[0]->Fill(EvtCent, runindex, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       //fHist_Corr3p_vsRun_EP_PN[1]->Fill(EvtCent, runindex, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);

       fHist_Corr2p_EP_Norm_PN[QAindex]->Fill(EvtCent,0.5,TMath::Cos((dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Norm_PN[QAindex]->Fill(EvtCent,1.5,TMath::Cos(2.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Norm_PN[QAindex]->Fill(EvtCent,2.5,TMath::Cos(3.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Norm_PN[QAindex]->Fill(EvtCent,3.5,TMath::Cos(4.*(dPhi1 - dPhi2)),WgtEP);

       //2p Vs Refmult:
       fHist_Corr2p_EP_Refm_PN[QAindex]->Fill(fRefMultCorr,0.5,TMath::Cos((dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Refm_PN[QAindex]->Fill(fRefMultCorr,1.5,TMath::Cos(2.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Refm_PN[QAindex]->Fill(fRefMultCorr,2.5,TMath::Cos(3.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Refm_PN[QAindex]->Fill(fRefMultCorr,3.5,TMath::Cos(4.*(dPhi1 - dPhi2)),WgtEP);




       if(cIndex<6){
         fHist_Corr3p_pTSum_EP_V0A_PN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
         fHist_Corr3p_pTSum_EP_V0C_PN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
         fHist_Corr3p_pTDiff_EP_V0A_PN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
         fHist_Corr3p_pTDiff_EP_V0C_PN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
         fHist_Corr3p_EtaDiff_EP_V0A_PN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
         fHist_Corr3p_EtaDiff_EP_V0C_PN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
       }
     }
     else if(dChrg1>0 && dChrg2>0 && skipPairHBT==0) {

       //if(deltaEta < fHBTCutValue)
       //cout<<"++ dEta "<<deltaEta<<"  dphi "<<deltaPhi<<" ch1 = "<<dChrg1<<" ch2 = "<<dChrg2<<" skipPairHBT = "<<skipPairHBT<<endl;
	
       fdPhiFemtoCut->Fill(3.5,1e-8);
	       
       fHist_Corr3p_EP_Norm_PP[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       fHist_Corr3p_EP_Norm_PP[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
       //vs Refmult:
       fHist_Corr3p_EP_Refm_PP[QAindex][0]->Fill(fRefMultCorr, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       fHist_Corr3p_EP_Refm_PP[QAindex][1]->Fill(fRefMultCorr, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);

       //fHist_Corr3p_vsRun_EP_PP[0]->Fill(EvtCent, runindex, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       //fHist_Corr3p_vsRun_EP_PP[1]->Fill(EvtCent, runindex, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);

       fHist_Corr2p_EP_Norm_PP[QAindex]->Fill(EvtCent,0.5,TMath::Cos((dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Norm_PP[QAindex]->Fill(EvtCent,1.5,TMath::Cos(2.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Norm_PP[QAindex]->Fill(EvtCent,2.5,TMath::Cos(3.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Norm_PP[QAindex]->Fill(EvtCent,3.5,TMath::Cos(4.*(dPhi1 - dPhi2)),WgtEP);

       //2p Vs Refmult:
       fHist_Corr2p_EP_Refm_PP[QAindex]->Fill(fRefMultCorr,0.5,TMath::Cos((dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Refm_PP[QAindex]->Fill(fRefMultCorr,1.5,TMath::Cos(2.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Refm_PP[QAindex]->Fill(fRefMultCorr,2.5,TMath::Cos(3.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Refm_PP[QAindex]->Fill(fRefMultCorr,3.5,TMath::Cos(4.*(dPhi1 - dPhi2)),WgtEP);



       if(cIndex<6){
         fHist_Corr3p_pTSum_EP_V0A_PP[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
         fHist_Corr3p_pTSum_EP_V0C_PP[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
         fHist_Corr3p_pTDiff_EP_V0A_PP[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
         fHist_Corr3p_pTDiff_EP_V0C_PP[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
         fHist_Corr3p_EtaDiff_EP_V0A_PP[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
         fHist_Corr3p_EtaDiff_EP_V0C_PP[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
       }
     }
     else if(dChrg1<0 && dChrg2<0 && skipPairHBT==0){
       
       fdPhiFemtoCut->Fill(3.5,1e-8);
       
       fHist_Corr3p_EP_Norm_NN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       fHist_Corr3p_EP_Norm_NN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
       //vs Refmult:
       fHist_Corr3p_EP_Refm_NN[QAindex][0]->Fill(fRefMultCorr, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       fHist_Corr3p_EP_Refm_NN[QAindex][1]->Fill(fRefMultCorr, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);

       //fHist_Corr3p_vsRun_EP_NN[0]->Fill(EvtCent, runindex, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
       //fHist_Corr3p_vsRun_EP_NN[1]->Fill(EvtCent, runindex, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);

       fHist_Corr2p_EP_Norm_NN[QAindex]->Fill(EvtCent,0.5,TMath::Cos((dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Norm_NN[QAindex]->Fill(EvtCent,1.5,TMath::Cos(2.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Norm_NN[QAindex]->Fill(EvtCent,2.5,TMath::Cos(3.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Norm_NN[QAindex]->Fill(EvtCent,3.5,TMath::Cos(4.*(dPhi1 - dPhi2)),WgtEP);
       //2p Vs Refmult:
       fHist_Corr2p_EP_Refm_NN[QAindex]->Fill(fRefMultCorr,0.5,TMath::Cos((dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Refm_NN[QAindex]->Fill(fRefMultCorr,1.5,TMath::Cos(2.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Refm_NN[QAindex]->Fill(fRefMultCorr,2.5,TMath::Cos(3.*(dPhi1 - dPhi2)),WgtEP);
       fHist_Corr2p_EP_Refm_NN[QAindex]->Fill(fRefMultCorr,3.5,TMath::Cos(4.*(dPhi1 - dPhi2)),WgtEP);


       if(cIndex<6){
         fHist_Corr3p_pTSum_EP_V0A_NN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
         fHist_Corr3p_pTSum_EP_V0C_NN[QAindex][cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
         fHist_Corr3p_pTDiff_EP_V0A_NN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
         fHist_Corr3p_pTDiff_EP_V0C_NN[QAindex][cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
         fHist_Corr3p_EtaDiff_EP_V0A_NN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0A),WgtEP);
         fHist_Corr3p_EtaDiff_EP_V0C_NN[QAindex][cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + m*dPhi2 - p*Psi2V0C),WgtEP);
       }
     }
   }//nested loop ends
 }//------ track loop ends ------

 fEventStatvsRun->Fill(runindex);

 fCentDistvsRun->Fill(EvtCent,runindex);


 fRefMultCorrvsRaw->Fill(fRefMultRaw,fRefMultCorr);



 if(EvtCent<=20){
   fCentV0MvsVzRun->Fill(VtxZ,runindex,EvtCent);
   fCentCL1vsVzRun->Fill(VtxZ,runindex,(centrCL1-EvtCent));
 }


 Double_t QTPCRe = QxTPC[1];
 Double_t QTPCIm = QyTPC[1];

 Double_t QTPCRe3 = QxTPC[2];
 Double_t QTPCIm3 = QyTPC[2];


 Double_t Psi2TPC = 0;

 if(psiN<=2){
   Psi2TPC = 1./2*(TMath::ATan2(QTPCIm,QTPCRe));
   if(Psi2TPC < 0.) Psi2TPC += TMath::Pi();

   fTPCQ2xVsCentRun->Fill(EvtCent,runindex,TMath::Cos(Psi2TPC));
   fTPCQ2yVsCentRun->Fill(EvtCent,runindex,TMath::Sin(Psi2TPC)); 
 }
 else if(psiN>=3){
   Psi2TPC = 1./psiN*(TMath::ATan2(QTPCIm3,QTPCRe3));
   if(Psi2TPC < 0.) Psi2TPC += 2./psiN*TMath::Pi();

   fTPCQ3xVsCentRun->Fill(EvtCent,runindex,TMath::Cos(Psi2TPC));
   fTPCQ3yVsCentRun->Fill(EvtCent,runindex,TMath::Sin(Psi2TPC));
 }




 PostData(3,fListNUAHist);   // Store NUA for all Runs including bad one. 


 if(QTPCRe==0 && QTPCIm==0) return;  // is this is one to one with bad run list?
 fHist_Event_count->Fill(stepCount); //6
 stepCount++;


 if(isBadRun) return;       
 fHist_Event_count->Fill(stepCount); //7
 stepCount++;



 Float_t fNegTracks = iTracks*0.5;
 Float_t ratioReject = (Float_t) icReject/fNegTracks;

 fRejectRatioVsCR->Fill(EvtCent,runindex,ratioReject);







//V0A-V0C EP resolution:
 fHist_Reso2n_EP_Norm_Det[QAindex][0]->Fill(EvtCent, TMath::Cos(psiN*(Psi2V0A-Psi2V0C)));
//V0A-TPC EP resolution:
 fHist_Reso2n_EP_Norm_Det[QAindex][1]->Fill(EvtCent, TMath::Cos(psiN*(Psi2V0A-Psi2TPC)));
//V0C-TPC EP resolution:
 fHist_Reso2n_EP_Norm_Det[QAindex][2]->Fill(EvtCent, TMath::Cos(psiN*(Psi2V0C-Psi2TPC)));

//V0A-V0C EP resolution:
 fHist_Reso2n_EP_Refm_Det[QAindex][0]->Fill(fRefMultCorr, TMath::Cos(psiN*(Psi2V0A-Psi2V0C)));
//V0A-TPC EP resolution:
 fHist_Reso2n_EP_Refm_Det[QAindex][1]->Fill(fRefMultCorr, TMath::Cos(psiN*(Psi2V0A-Psi2TPC)));
//V0C-TPC EP resolution:
 fHist_Reso2n_EP_Refm_Det[QAindex][2]->Fill(fRefMultCorr, TMath::Cos(psiN*(Psi2V0C-Psi2TPC)));


 fHV0AEventPlaneVsCent->Fill(EvtCent,Psi2V0A);
 fHV0CEventPlaneVsCent->Fill(EvtCent,Psi2V0C);
 fHTPCEventPlaneVsCent->Fill(EvtCent,Psi2TPC);






 //------- SP Method -------------


 Double_t uPRe=0.,uNRe=0.,uPIm=0.,uNIm=0.,uN2Re=0.,uN2Im=0.,uP2Re=0.,uP2Im=0.;
 Double_t uPM=0.,uNM=0.;
 //Double_t VCIm=0., VPRe=0.,VPIm=0.;


 uPM = MPOIpos; uNM = MPOIneg;

 uPRe =  QxPos[0];  uNRe =  QxNeg[0]; 
 uPIm =  QyPos[0];  uNIm =  QyNeg[0];

 uP2Re = QxAutoPos[0]; uP2Im =  QyAutoPos[0];
 uN2Re = QxAutoNeg[0]; uN2Im =  QyAutoNeg[0];

 Double_t TwoQpQnV = 0.,TwoQpQpV=0.,TwoQnQnV=0.;

 if(uPM > 1 && uNM > 1) {
   //CME w.r.t V0A EP:
   TwoQpQnV = ((uPRe*uNRe-uPIm*uNIm)*QxanCor + (uPRe*uNIm+uPIm*uNRe)*QyanCor) / (uPM*uNM) ;
   TwoQpQpV = ((uPRe*uPRe-uPIm*uPIm-uP2Re)*QxanCor + (2.*uPRe*uPIm-uP2Im)*QyanCor) / (uPM*(uPM-1)) ;
   TwoQnQnV = ((uNRe*uNRe-uNIm*uNIm-uN2Re)*QxanCor + (2.*uNRe*uNIm-uN2Im)*QyanCor) / (uNM*(uNM-1)) ;
   //Fill profiles:
   /*
   fHist_Corr3p_SP_Norm_PN[QAindex][0]->Fill(EvtCent, TwoQpQnV, uPM*uNM);
   fHist_Corr3p_SP_Norm_PP[QAindex][0]->Fill(EvtCent, TwoQpQpV, (uPM*(uPM-1)));
   fHist_Corr3p_SP_Norm_NN[QAindex][0]->Fill(EvtCent, TwoQnQnV, (uNM*(uNM-1)));
   */
  //-------- Fill NonIsotropic terms ------

   Double_t QnNonIsoRe = 0.,QnNonIsoIm=0.;

   if(QAindex==0){ //B < 0
    //charge pos:
     fHist_NonIso_SP_PP_Mag0[0]->Fill(EvtCent, 0.5,uPRe/uPM,uPM);  //<Cos(nPhi)> ChPos
     fHist_NonIso_SP_PP_Mag0[0]->Fill(EvtCent, 1.5,uPIm/uPM,uPM);  //<Sin(nPhi)> ChPos

     QnNonIsoRe = (uPRe*QxanCor + uPIm*QyanCor)/uPM;
     QnNonIsoIm = (uPIm*QxanCor - uPRe*QyanCor)/uPM; 
     fHist_NonIso_SP_PP_Mag0[0]->Fill(EvtCent, 2.5,QnNonIsoRe,uPM);  //<Cos(nPhi-mPsi)> ChPos
     fHist_NonIso_SP_PP_Mag0[0]->Fill(EvtCent, 3.5,QnNonIsoIm,uPM);  //<Sin(nPhi-mPsi)> ChPos

     QnNonIsoRe = (uPRe*uPRe-uPIm*uPIm-uP2Re)/(uPM*(uPM-1)); 
     QnNonIsoIm = (2.*uPRe*uPIm-uP2Im)/(uPM*(uPM-1));
     fHist_NonIso_SP_PP_Mag0[0]->Fill(EvtCent, 4.5,QnNonIsoRe,(uPM*(uPM-1))); //<Cos(nPhi1+nPhi2)> ChPos
     fHist_NonIso_SP_PP_Mag0[0]->Fill(EvtCent, 5.5,QnNonIsoIm,(uPM*(uPM-1))); //<Sin(nPhi1+nPhi2)> ChPos

     QnNonIsoRe = (uPRe*uNRe-uPIm*uNIm)/(uPM*uNM);
     QnNonIsoIm = (uPRe*uNIm+uPIm*uNRe)/(uPM*uNM);
     fHist_NonIso_SP_PP_Mag0[0]->Fill(EvtCent, 6.5,QnNonIsoRe,(uPM*uNM)); //<Cos(nPhi1+nPhi2)> phi1,phi2 opposite charge
     fHist_NonIso_SP_PP_Mag0[0]->Fill(EvtCent, 7.5,QnNonIsoIm,(uPM*uNM)); //<Sin(nPhi1+nPhi2)> phi1,phi2 opposite charge

     fHist_NonIso_SP_PP_Mag0[0]->Fill(EvtCent, 8.5,QxanCor,uPM); //<Cos(mPsiEP)> 
     fHist_NonIso_SP_PP_Mag0[0]->Fill(EvtCent, 9.5,QyanCor,uPM); //<Sin(mPsiEP)>

     //charge neg:
     fHist_NonIso_SP_NN_Mag0[0]->Fill(EvtCent, 0.5,uNRe/uNM,uNM);  //<Cos(nPhi)> ChNeg
     fHist_NonIso_SP_NN_Mag0[0]->Fill(EvtCent, 1.5,uNIm/uNM,uNM);  //<Sin(nPhi)> ChNeg

     QnNonIsoRe = (uNRe*QxanCor + uNIm*QyanCor)/uNM;
     QnNonIsoIm = (uNIm*QxanCor - uNRe*QyanCor)/uNM; 
     fHist_NonIso_SP_NN_Mag0[0]->Fill(EvtCent, 2.5,QnNonIsoRe,uNM); //<Cos(nPhi-mPsi)> ChNeg
     fHist_NonIso_SP_NN_Mag0[0]->Fill(EvtCent, 3.5,QnNonIsoIm,uNM); //<Cos(nPhi-mPsi)> ChNeg

     QnNonIsoRe = (uNRe*uNRe-uNIm*uNIm-uN2Re)/(uNM*(uNM-1)); 
     QnNonIsoIm = (2.*uNRe*uNIm-uN2Im)/(uNM*(uNM-1));
     fHist_NonIso_SP_NN_Mag0[0]->Fill(EvtCent, 4.5,QnNonIsoRe,(uNM*(uNM-1))); //<Cos(nPhi1+nPhi2)> ChNeg
     fHist_NonIso_SP_NN_Mag0[0]->Fill(EvtCent, 5.5,QnNonIsoIm,(uNM*(uNM-1))); //<Sin(nPhi1+nPhi2)> ChNeg

     QnNonIsoRe = (uPRe*uNRe-uPIm*uNIm)/(uPM*uNM);
     QnNonIsoIm = (uPRe*uNIm+uPIm*uNRe)/(uPM*uNM);
     fHist_NonIso_SP_NN_Mag0[0]->Fill(EvtCent, 6.5,QnNonIsoRe,(uPM*uNM)); //<Cos(nPhi1+nPhi2)> phi1,phi2 opposite charge
     fHist_NonIso_SP_NN_Mag0[0]->Fill(EvtCent, 7.5,QnNonIsoIm,(uPM*uNM)); //<Sin(nPhi1+nPhi2)> phi1,phi2 opposite charge

     fHist_NonIso_SP_NN_Mag0[0]->Fill(EvtCent, 8.5,QxanCor,uNM); //<Cos(mPsiEP)> 
     fHist_NonIso_SP_NN_Mag0[0]->Fill(EvtCent, 9.5,QyanCor,uNM); //<Sin(mPsiEP)>
   }
   else if(QAindex==1){ //B > 0
    //charge pos:
     fHist_NonIso_SP_PP_Mag1[0]->Fill(EvtCent, 0.5,uPRe/uPM,uPM);  //<Cos(nPhi)> ChPos
     fHist_NonIso_SP_PP_Mag1[0]->Fill(EvtCent, 1.5,uPIm/uPM,uPM);  //<Sin(nPhi)> ChPos

     QnNonIsoRe = (uPRe*QxanCor + uPIm*QyanCor)/uPM;
     QnNonIsoIm = (uPIm*QxanCor - uPRe*QyanCor)/uPM; 
     fHist_NonIso_SP_PP_Mag1[0]->Fill(EvtCent, 2.5,QnNonIsoRe,uPM);  //<Cos(nPhi-mPsi)> ChPos
     fHist_NonIso_SP_PP_Mag1[0]->Fill(EvtCent, 3.5,QnNonIsoIm,uPM);  //<Sin(nPhi-mPsi)> ChPos

     QnNonIsoRe = (uPRe*uPRe-uPIm*uPIm-uP2Re)/(uPM*(uPM-1)); 
     QnNonIsoIm = (2.*uPRe*uPIm-uP2Im)/(uPM*(uPM-1));
     fHist_NonIso_SP_PP_Mag1[0]->Fill(EvtCent, 4.5,QnNonIsoRe,(uPM*(uPM-1))); //<Cos(nPhi1+nPhi2)> ChPos
     fHist_NonIso_SP_PP_Mag1[0]->Fill(EvtCent, 5.5,QnNonIsoIm,(uPM*(uPM-1))); //<Sin(nPhi1+nPhi2)> ChPos

     QnNonIsoRe = (uPRe*uNRe-uPIm*uNIm)/(uPM*uNM);
     QnNonIsoIm = (uPRe*uNIm+uPIm*uNRe)/(uPM*uNM);
     fHist_NonIso_SP_PP_Mag1[0]->Fill(EvtCent, 6.5,QnNonIsoRe,(uPM*uNM)); //<Cos(nPhi1+nPhi2)> phi1,phi2 opposite charge
     fHist_NonIso_SP_PP_Mag1[0]->Fill(EvtCent, 7.5,QnNonIsoIm,(uPM*uNM)); //<Sin(nPhi1+nPhi2)> phi1,phi2 opposite charge

     fHist_NonIso_SP_PP_Mag1[0]->Fill(EvtCent, 8.5,QxanCor,uPM); //<Cos(mPsiEP)> 
     fHist_NonIso_SP_PP_Mag1[0]->Fill(EvtCent, 9.5,QyanCor,uPM); //<Sin(mPsiEP)>

     //charge neg:
     fHist_NonIso_SP_NN_Mag1[0]->Fill(EvtCent, 0.5,uNRe/uNM,uNM);  //<Cos(nPhi)> ChNeg
     fHist_NonIso_SP_NN_Mag1[0]->Fill(EvtCent, 1.5,uNIm/uNM,uNM);  //<Sin(nPhi)> ChNeg

     QnNonIsoRe = (uNRe*QxanCor + uNIm*QyanCor)/uNM;
     QnNonIsoIm = (uNIm*QxanCor - uNRe*QyanCor)/uNM; 
     fHist_NonIso_SP_NN_Mag1[0]->Fill(EvtCent, 2.5,QnNonIsoRe,uNM); //<Cos(nPhi-mPsi)> ChNeg
     fHist_NonIso_SP_NN_Mag1[0]->Fill(EvtCent, 3.5,QnNonIsoIm,uNM); //<Cos(nPhi-mPsi)> ChNeg

     QnNonIsoRe = (uNRe*uNRe-uNIm*uNIm-uN2Re)/(uNM*(uNM-1)); 
     QnNonIsoIm = (2.*uNRe*uNIm-uN2Im)/(uNM*(uNM-1));
     fHist_NonIso_SP_NN_Mag1[0]->Fill(EvtCent, 4.5,QnNonIsoRe,(uNM*(uNM-1))); //<Cos(nPhi1+nPhi2)> ChNeg
     fHist_NonIso_SP_NN_Mag1[0]->Fill(EvtCent, 5.5,QnNonIsoIm,(uNM*(uNM-1))); //<Sin(nPhi1+nPhi2)> ChNeg

     QnNonIsoRe = (uPRe*uNRe-uPIm*uNIm)/(uPM*uNM);
     QnNonIsoIm = (uPRe*uNIm+uPIm*uNRe)/(uPM*uNM);
     fHist_NonIso_SP_NN_Mag1[0]->Fill(EvtCent, 6.5,QnNonIsoRe,(uPM*uNM)); //<Cos(nPhi1+nPhi2)> phi1,phi2 opposite charge
     fHist_NonIso_SP_NN_Mag1[0]->Fill(EvtCent, 7.5,QnNonIsoIm,(uPM*uNM)); //<Sin(nPhi1+nPhi2)> phi1,phi2 opposite charge

     fHist_NonIso_SP_NN_Mag1[0]->Fill(EvtCent, 8.5,QxanCor,uNM); //<Cos(mPsiEP)> 
     fHist_NonIso_SP_NN_Mag1[0]->Fill(EvtCent, 9.5,QyanCor,uNM); //<Sin(mPsiEP)>
   }


   //-------- CME - ZDN correlators:---------
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  0., TwoQpQnV, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  1., energyZPA, 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  2., energyZPC, 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  4., energyZNA, 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  5., energyZNC, 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  7., TwoQpQnV*energyZPA, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  8., TwoQpQnV*energyZPC, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  9., TwoQpQnV*(energyZPC+energyZPA), uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent, 10., TwoQpQnV*energyZNA, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent, 11., TwoQpQnV*energyZNC, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent, 12., TwoQpQnV*(energyZNC+energyZNA), uPM*uNM);

   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  0., TwoQpQpV, (uPM*(uPM-1)));
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  1., energyZPA, 1.);
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  2., energyZPC, 1.);
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  4., energyZNA, 1.);
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  5., energyZNC, 1.);
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  7., TwoQpQpV*energyZPA, (uPM*(uPM-1)));
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  8., TwoQpQpV*energyZPC, (uPM*(uPM-1)));
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  9., TwoQpQpV*(energyZPC+energyZPA), (uPM*(uPM-1)));
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent, 10., TwoQpQpV*energyZNA, (uPM*(uPM-1)));
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent, 11., TwoQpQpV*energyZNC, (uPM*(uPM-1)));
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent, 12., TwoQpQpV*(energyZNC+energyZNA), (uPM*(uPM-1)));

   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  0., TwoQnQnV, (uNM*(uNM-1)));
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  1., energyZPA, 1.);
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  2., energyZPC, 1.);
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  4., energyZNA, 1.);
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  5., energyZNC, 1.);
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  7., TwoQnQnV*energyZPA, (uNM*(uNM-1)));
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  8., TwoQnQnV*energyZPC, (uNM*(uNM-1)));
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  9., TwoQnQnV*(energyZPC+energyZPA), (uNM*(uNM-1)));
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent, 10., TwoQnQnV*energyZNA, (uNM*(uNM-1)));
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent, 11., TwoQnQnV*energyZNC, (uNM*(uNM-1)));
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent, 12., TwoQnQnV*(energyZNC+energyZNA), (uNM*(uNM-1)));
   //-----------------------------------------




   //CME w.r.t V0C EP:
   TwoQpQnV = 0.; TwoQpQpV = 0.; TwoQnQnV = 0.;
   TwoQpQnV = ((uPRe*uNRe-uPIm*uNIm)*QxcnCor + (uPRe*uNIm+uPIm*uNRe)*QycnCor) / (uPM*uNM) ;
   TwoQpQpV = ((uPRe*uPRe-uPIm*uPIm-uP2Re)*QxcnCor + (2.*uPRe*uPIm-uP2Im)*QycnCor) / (uPM*(uPM-1)) ;
   TwoQnQnV = ((uNRe*uNRe-uNIm*uNIm-uN2Re)*QxcnCor + (2.*uNRe*uNIm-uN2Im)*QycnCor) / (uNM*(uNM-1)) ;

   //Fill profiles:
  /*
   fHist_Corr3p_SP_Norm_PN[QAindex][1]->Fill(EvtCent, TwoQpQnV, uPM*uNM);
   fHist_Corr3p_SP_Norm_PP[QAindex][1]->Fill(EvtCent, TwoQpQpV, (uPM*(uPM-1)));
   fHist_Corr3p_SP_Norm_NN[QAindex][1]->Fill(EvtCent, TwoQnQnV, (uNM*(uNM-1)));
  */

  if(QAindex==0){ //B < 0
    //charge pos:
     fHist_NonIso_SP_PP_Mag0[1]->Fill(EvtCent, 0.5,uPRe/uPM,uPM);  //<Cos(nPhi)> ChPos
     fHist_NonIso_SP_PP_Mag0[1]->Fill(EvtCent, 1.5,uPIm/uPM,uPM);  //<Sin(nPhi)> ChPos

     QnNonIsoRe = (uPRe*QxcnCor + uPIm*QycnCor)/uPM;
     QnNonIsoIm = (uPIm*QxcnCor - uPRe*QycnCor)/uPM; 
     fHist_NonIso_SP_PP_Mag0[1]->Fill(EvtCent, 2.5,QnNonIsoRe,uPM);  //<Cos(nPhi-mPsi)> ChPos
     fHist_NonIso_SP_PP_Mag0[1]->Fill(EvtCent, 3.5,QnNonIsoIm,uPM);  //<Sin(nPhi-mPsi)> ChPos

     QnNonIsoRe = (uPRe*uPRe-uPIm*uPIm-uP2Re)/(uPM*(uPM-1)); 
     QnNonIsoIm = (2.*uPRe*uPIm-uP2Im)/(uPM*(uPM-1));
     fHist_NonIso_SP_PP_Mag0[1]->Fill(EvtCent, 4.5,QnNonIsoRe,(uPM*(uPM-1))); //<Cos(nPhi1+nPhi2)> ChPos
     fHist_NonIso_SP_PP_Mag0[1]->Fill(EvtCent, 5.5,QnNonIsoIm,(uPM*(uPM-1))); //<Sin(nPhi1+nPhi2)> ChPos

     QnNonIsoRe = (uPRe*uNRe-uPIm*uNIm)/(uPM*uNM);
     QnNonIsoIm = (uPRe*uNIm+uPIm*uNRe)/(uPM*uNM);
     fHist_NonIso_SP_PP_Mag0[1]->Fill(EvtCent, 6.5,QnNonIsoRe,(uPM*uNM)); //<Cos(nPhi1+nPhi2)> phi1,phi2 opposite charge
     fHist_NonIso_SP_PP_Mag0[1]->Fill(EvtCent, 7.5,QnNonIsoIm,(uPM*uNM)); //<Sin(nPhi1+nPhi2)> phi1,phi2 opposite charge

     fHist_NonIso_SP_PP_Mag0[1]->Fill(EvtCent, 8.5,QxcnCor,uPM); //<Cos(mPsiEP)> 
     fHist_NonIso_SP_PP_Mag0[1]->Fill(EvtCent, 9.5,QycnCor,uPM); //<Sin(mPsiEP)>

     //charge neg:
     fHist_NonIso_SP_NN_Mag0[1]->Fill(EvtCent, 0.5,uNRe/uNM,uNM);  //<Cos(nPhi)> ChNeg
     fHist_NonIso_SP_NN_Mag0[1]->Fill(EvtCent, 1.5,uNIm/uNM,uNM);  //<Sin(nPhi)> ChNeg

     QnNonIsoRe = (uNRe*QxcnCor + uNIm*QycnCor)/uNM;
     QnNonIsoIm = (uNIm*QxcnCor - uNRe*QycnCor)/uNM; 
     fHist_NonIso_SP_NN_Mag0[1]->Fill(EvtCent, 2.5,QnNonIsoRe,uNM); //<Cos(nPhi-mPsi)> ChNeg
     fHist_NonIso_SP_NN_Mag0[1]->Fill(EvtCent, 3.5,QnNonIsoIm,uNM); //<Cos(nPhi-mPsi)> ChNeg

     QnNonIsoRe = (uNRe*uNRe-uNIm*uNIm-uN2Re)/(uNM*(uNM-1)); 
     QnNonIsoIm = (2.*uNRe*uNIm-uN2Im)/(uNM*(uNM-1));
     fHist_NonIso_SP_NN_Mag0[1]->Fill(EvtCent, 4.5,QnNonIsoRe,(uNM*(uNM-1))); //<Cos(nPhi1+nPhi2)> ChNeg
     fHist_NonIso_SP_NN_Mag0[1]->Fill(EvtCent, 5.5,QnNonIsoIm,(uNM*(uNM-1))); //<Sin(nPhi1+nPhi2)> ChNeg

     QnNonIsoRe = (uPRe*uNRe-uPIm*uNIm)/(uPM*uNM);
     QnNonIsoIm = (uPRe*uNIm+uPIm*uNRe)/(uPM*uNM);
     fHist_NonIso_SP_NN_Mag0[1]->Fill(EvtCent, 6.5,QnNonIsoRe,(uPM*uNM)); //<Cos(nPhi1+nPhi2)> phi1,phi2 opposite charge
     fHist_NonIso_SP_NN_Mag0[1]->Fill(EvtCent, 7.5,QnNonIsoIm,(uPM*uNM)); //<Sin(nPhi1+nPhi2)> phi1,phi2 opposite charge

     fHist_NonIso_SP_NN_Mag0[1]->Fill(EvtCent, 8.5,QxcnCor,uNM); //<Cos(mPsiEP)> 
     fHist_NonIso_SP_NN_Mag0[1]->Fill(EvtCent, 9.5,QycnCor,uNM); //<Sin(mPsiEP)>
   }
   else if(QAindex==1) { //B > 0
    //charge pos:
     fHist_NonIso_SP_PP_Mag1[1]->Fill(EvtCent, 0.5,uPRe/uPM,uPM);  //<Cos(nPhi)> ChPos
     fHist_NonIso_SP_PP_Mag1[1]->Fill(EvtCent, 1.5,uPIm/uPM,uPM);  //<Sin(nPhi)> ChPos

     QnNonIsoRe = (uPRe*QxcnCor + uPIm*QycnCor)/uPM;
     QnNonIsoIm = (uPIm*QxcnCor - uPRe*QycnCor)/uPM; 
     fHist_NonIso_SP_PP_Mag1[1]->Fill(EvtCent, 2.5,QnNonIsoRe,uPM);  //<Cos(nPhi-mPsi)> ChPos
     fHist_NonIso_SP_PP_Mag1[1]->Fill(EvtCent, 3.5,QnNonIsoIm,uPM);  //<Sin(nPhi-mPsi)> ChPos

     QnNonIsoRe = (uPRe*uPRe-uPIm*uPIm-uP2Re)/(uPM*(uPM-1)); 
     QnNonIsoIm = (2.*uPRe*uPIm-uP2Im)/(uPM*(uPM-1));
     fHist_NonIso_SP_PP_Mag1[1]->Fill(EvtCent, 4.5,QnNonIsoRe,(uPM*(uPM-1))); //<Cos(nPhi1+nPhi2)> ChPos
     fHist_NonIso_SP_PP_Mag1[1]->Fill(EvtCent, 5.5,QnNonIsoIm,(uPM*(uPM-1))); //<Sin(nPhi1+nPhi2)> ChPos

     QnNonIsoRe = (uPRe*uNRe-uPIm*uNIm)/(uPM*uNM);
     QnNonIsoIm = (uPRe*uNIm+uPIm*uNRe)/(uPM*uNM);
     fHist_NonIso_SP_PP_Mag1[1]->Fill(EvtCent, 6.5,QnNonIsoRe,(uPM*uNM)); //<Cos(nPhi1+nPhi2)> phi1,phi2 opposite charge
     fHist_NonIso_SP_PP_Mag1[1]->Fill(EvtCent, 7.5,QnNonIsoIm,(uPM*uNM)); //<Sin(nPhi1+nPhi2)> phi1,phi2 opposite charge

     fHist_NonIso_SP_PP_Mag1[1]->Fill(EvtCent, 8.5,QxcnCor,uPM); //<Cos(mPsiEP)> 
     fHist_NonIso_SP_PP_Mag1[1]->Fill(EvtCent, 9.5,QycnCor,uPM); //<Sin(mPsiEP)>

     //charge neg:
     fHist_NonIso_SP_NN_Mag1[1]->Fill(EvtCent, 0.5,uNRe/uNM,uNM);  //<Cos(nPhi)> ChNeg
     fHist_NonIso_SP_NN_Mag1[1]->Fill(EvtCent, 1.5,uNIm/uNM,uNM);  //<Sin(nPhi)> ChNeg

     QnNonIsoRe = (uNRe*QxcnCor + uNIm*QycnCor)/uNM;
     QnNonIsoIm = (uNIm*QxcnCor - uNRe*QycnCor)/uNM; 
     fHist_NonIso_SP_NN_Mag1[1]->Fill(EvtCent, 2.5,QnNonIsoRe,uNM); //<Cos(nPhi-mPsi)> ChNeg
     fHist_NonIso_SP_NN_Mag1[1]->Fill(EvtCent, 3.5,QnNonIsoIm,uNM); //<Cos(nPhi-mPsi)> ChNeg

     QnNonIsoRe = (uNRe*uNRe-uNIm*uNIm-uN2Re)/(uNM*(uNM-1)); 
     QnNonIsoIm = (2.*uNRe*uNIm-uN2Im)/(uNM*(uNM-1));
     fHist_NonIso_SP_NN_Mag1[1]->Fill(EvtCent, 4.5,QnNonIsoRe,(uNM*(uNM-1))); //<Cos(nPhi1+nPhi2)> ChNeg
     fHist_NonIso_SP_NN_Mag1[1]->Fill(EvtCent, 5.5,QnNonIsoIm,(uNM*(uNM-1))); //<Sin(nPhi1+nPhi2)> ChNeg

     QnNonIsoRe = (uPRe*uNRe-uPIm*uNIm)/(uPM*uNM);
     QnNonIsoIm = (uPRe*uNIm+uPIm*uNRe)/(uPM*uNM);
     fHist_NonIso_SP_NN_Mag1[1]->Fill(EvtCent, 6.5,QnNonIsoRe,(uPM*uNM)); //<Cos(nPhi1+nPhi2)> phi1,phi2 opposite charge
     fHist_NonIso_SP_NN_Mag1[1]->Fill(EvtCent, 7.5,QnNonIsoIm,(uPM*uNM)); //<Sin(nPhi1+nPhi2)> phi1,phi2 opposite charge

     fHist_NonIso_SP_NN_Mag1[1]->Fill(EvtCent, 8.5,QxcnCor,uNM); //<Cos(mPsiEP)> 
     fHist_NonIso_SP_NN_Mag1[1]->Fill(EvtCent, 9.5,QycnCor,uNM); //<Sin(mPsiEP)>
   }


   //-------- CME - ZDN correlators:---------
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  0., TwoQpQnV, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  1., energyZPA, 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  2., energyZPC, 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  4., energyZNA, 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  5., energyZNC, 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  7., TwoQpQnV*energyZPA, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  8., TwoQpQnV*energyZPC, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  9., TwoQpQnV*(energyZPC+energyZPA), uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent, 10., TwoQpQnV*energyZNA, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent, 11., TwoQpQnV*energyZNC, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent, 12., TwoQpQnV*(energyZNC+energyZNA), uPM*uNM);

   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  0., TwoQpQpV, (uPM*(uPM-1)));
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  1., energyZPA, 1.);
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  2., energyZPC, 1.);
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  4., energyZNA, 1.);
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  5., energyZNC, 1.);
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  7., TwoQpQpV*energyZPA, (uPM*(uPM-1)));
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  8., TwoQpQpV*energyZPC, (uPM*(uPM-1)));
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  9., TwoQpQpV*(energyZPC+energyZPA), (uPM*(uPM-1)));
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent, 10., TwoQpQpV*energyZNA, (uPM*(uPM-1)));
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent, 11., TwoQpQpV*energyZNC, (uPM*(uPM-1)));
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent, 12., TwoQpQpV*(energyZNC+energyZNA), (uPM*(uPM-1)));

   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  0., TwoQnQnV, (uNM*(uNM-1)));
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  1., energyZPA, 1.);
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  2., energyZPC, 1.);
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  4., energyZNA, 1.);
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  5., energyZNC, 1.);
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  7., TwoQnQnV*energyZPA, (uNM*(uNM-1)));
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  8., TwoQnQnV*energyZPC, (uNM*(uNM-1)));
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  9., TwoQnQnV*(energyZPC+energyZPA), (uNM*(uNM-1)));
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent, 10., TwoQnQnV*energyZNA, (uNM*(uNM-1)));
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent, 11., TwoQnQnV*energyZNC, (uNM*(uNM-1)));
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent, 12., TwoQnQnV*(energyZNC+energyZNA), (uNM*(uNM-1)));
   //-----------------------------------------


   //Fill Mult, POIs vs cent:

   fAvgMultCentRun->Fill(EvtCent,runindex,nRefMult); 
   fAvgWgtMultCentRun->Fill(EvtCent,runindex,nRefMultWgt);
   fAvgPOIposCentRun->Fill(EvtCent,runindex,uPM); 
   fAvgPOInegCentRun->Fill(EvtCent,runindex,uNM); 
   fAvgPOIPPCentRun->Fill(EvtCent,runindex,(uPM*(uPM-1)));
   fAvgPOINNCentRun->Fill(EvtCent,runindex,(uNM*(uNM-1)));
   fAvgPOIOSCentRun->Fill(EvtCent,runindex,uPM*uNM);


   //------------ eta depedence -----------------
   Double_t fEtaCent = -100.;
   Int_t nBinsX = fEtaBinFinderForQA->GetNbinsX();

   for(int ie=0; ie<nBinsX; ie++){
     
     uPM = MPOIposQAEta[ie];
     uNM = MPOInegQAEta[ie];

     if(uPM>1 && uNM>1){

       uPRe =  QxPosQAEta[ie];
       uNRe =  QxNegQAEta[ie]; 
       uPIm =  QyPosQAEta[ie];
       uNIm =  QyNegQAEta[ie];

       uP2Re = QxAutoPosQAEta[ie];
       uP2Im = QyAutoPosQAEta[ie];
       uN2Re = QxAutoNegQAEta[ie];
       uN2Im = QyAutoNegQAEta[ie];


       fEtaCent = fEtaBinFinderForQA->GetBinCenter(ie+1);

       //w.r.t V0A:
       TwoQpQnV = 0.; TwoQpQpV = 0.; TwoQnQnV = 0.;
       TwoQpQnV = ((uPRe*uNRe-uPIm*uNIm)*QxanCor + (uPRe*uNIm+uPIm*uNRe)*QyanCor) / (uPM*uNM) ;
       TwoQpQpV = ((uPRe*uPRe-uPIm*uPIm-uP2Re)*QxanCor + (2.*uPRe*uPIm-uP2Im)*QyanCor) / (uPM*(uPM-1.)) ;
       TwoQnQnV = ((uNRe*uNRe-uNIm*uNIm-uN2Re)*QxanCor + (2.*uNRe*uNIm-uN2Im)*QyanCor) / (uNM*(uNM-1.)) ;

       fHist_Corr3p_QAEta_SP_V0A_PN[QAindex]->Fill(fEtaCent,TwoQpQnV,uPM*uNM);
       fHist_Corr3p_QAEta_SP_V0A_PP[QAindex]->Fill(fEtaCent,TwoQpQpV,uPM*(uPM-1.));
       fHist_Corr3p_QAEta_SP_V0A_NN[QAindex]->Fill(fEtaCent,TwoQnQnV,uNM*(uNM-1.));

       //w.r.t V0C:
       TwoQpQnV = 0.; TwoQpQpV = 0.; TwoQnQnV = 0.;
       TwoQpQnV = ((uPRe*uNRe-uPIm*uNIm)*QxcnCor + (uPRe*uNIm+uPIm*uNRe)*QycnCor) / (uPM*uNM) ;
       TwoQpQpV = ((uPRe*uPRe-uPIm*uPIm-uP2Re)*QxcnCor + (2.*uPRe*uPIm-uP2Im)*QycnCor) / (uPM*(uPM-1.)) ;
       TwoQnQnV = ((uNRe*uNRe-uNIm*uNIm-uN2Re)*QxcnCor + (2.*uNRe*uNIm-uN2Im)*QycnCor) / (uNM*(uNM-1.)) ;

       fHist_Corr3p_QAEta_SP_V0C_PN[QAindex]->Fill(fEtaCent,TwoQpQnV,uPM*uNM);
       fHist_Corr3p_QAEta_SP_V0C_PP[QAindex]->Fill(fEtaCent,TwoQpQpV,uPM*(uPM-1.));
       fHist_Corr3p_QAEta_SP_V0C_NN[QAindex]->Fill(fEtaCent,TwoQnQnV,uNM*(uNM-1.));
     }
   }
   //---------------------------------------




   if(MPOIpos>1 && MPOIneg>1){
     //V0A-V0C SP resolution:
     fHist_Reso2n_SP_Norm_Det[QAindex][0]->Fill(EvtCent, (QxcnCor*QxanCor+QycnCor*QyanCor));
     //V0A-TPC SP resolution:
     fHist_Reso2n_SP_Norm_Det[QAindex][1]->Fill(EvtCent, (QTPCRe*QxanCor+QTPCIm*QyanCor));
     //V0C-TPC SP resolution:
     fHist_Reso2n_SP_Norm_Det[QAindex][2]->Fill(EvtCent, (QTPCRe*QxcnCor+QTPCIm*QycnCor));
  

   //V0A-V0C:
     fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 0., (QxcnCor*QxanCor+QycnCor*QyanCor)); 
     fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 1., (QxcnCor*QxanCor+QycnCor*QyanCor)*energyZPA); 
     fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 2., (QxcnCor*QxanCor+QycnCor*QyanCor)*energyZPC); 
     fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 3., (QxcnCor*QxanCor+QycnCor*QyanCor)*(energyZPA+energyZPC)); 
     fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 4., (QxcnCor*QxanCor+QycnCor*QyanCor)*energyZNA); 
     fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 5., (QxcnCor*QxanCor+QycnCor*QyanCor)*energyZNC); 
     fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 6., (QxcnCor*QxanCor+QycnCor*QyanCor)*(energyZNA+energyZNC)); 

   //V0A-TPC 
     fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 0., (QTPCRe*QxanCor+QTPCIm*QyanCor));
     fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 1., (QTPCRe*QxanCor+QTPCIm*QyanCor)*energyZPA);
     fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 2., (QTPCRe*QxanCor+QTPCIm*QyanCor)*energyZPC);
     fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 3., (QTPCRe*QxanCor+QTPCIm*QyanCor)*(energyZPA+energyZPC));
     fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 4., (QTPCRe*QxanCor+QTPCIm*QyanCor)*energyZNA);
     fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 5., (QTPCRe*QxanCor+QTPCIm*QyanCor)*energyZNC);
     fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 6., (QTPCRe*QxanCor+QTPCIm*QyanCor)*(energyZNA+energyZNC));

  //V0C-TPC
     fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 0., (QTPCRe*QxcnCor+QTPCIm*QycnCor));
     fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 1., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*energyZPA);
     fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 2., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*energyZPC);
     fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 3., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*(energyZPA+energyZPC));
     fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 1., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*energyZNA);
     fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 2., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*energyZNC);
     fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 3., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*(energyZNA+energyZNC));
   }

 }//---- SP method -------


 PostData(1,fListHistos);
 PostData(2,fListCalibs); 
 
 fHist_Event_count->Fill(9.5);

 //if(fievent%100==0) {
 //cout<<"irun = "<<runindex<<" n "<<n<<" m = "<<m<<" p = "<<p<<" cent= "<<EvtCent<<"\tiCentSPD = "<<iCentSPD<<"\tQRe= "<<QTPCRe<<"\tQIm= "<<QTPCIm<<"\tpsiN = "<<psiN<<endl;
 //cout<<" cent= "<<EvtCent<<"\teZNC= "<<energyZNC<<"\teZPC = "<<energyZPC<<"\teZNA= "<<energyZNA<<"\teZPA = "<<energyZPA<<endl;
 //}

 fievent++;

}
//======================= UserExec done =========================



void AliAnalysisTaskCMEV0::Terminate(Option_t *)
{
 AliDebug(2,"\n ... AliAnalysisTaskCMEV0::Terminate() is being called ...  \n");
}













void AliAnalysisTaskCMEV0::GetV0QvectAndMult(const AliAODVZERO *aodV0,Float_t fHarmonic, Double_t& Qxan,Double_t& Qyan,Double_t& sumMa,Double_t& Qxcn,Double_t& Qycn,Double_t& sumMc) 
{
  for(Int_t iV0 = 0; iV0 < 64; iV0++) {
  /*if(fRemChV0A){
     if(iV0 == 46)
        continue;
    }*/
    Double_t phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);
    Float_t multv0 = aodV0->GetMultiplicity(iV0);

    if(iV0 < 32) {
      Double_t multCorC = -10;

      if(iV0 < 8)
         multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(1);
      else if(iV0 >= 8 && iV0 < 16)
         multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(9);
      else if(iV0 >= 16 && iV0 < 24)
         multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(17);
      else if(iV0 >= 24 && iV0 < 32)
         multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(25);

      if(multCorC < 0){
         cout<<"Problem with multiplicity in V0C"<<endl;
         continue;
      }
      Qxcn += TMath::Cos(fHarmonic*phiV0) * multCorC;
      Qycn += TMath::Sin(fHarmonic*phiV0) * multCorC;

      sumMc = sumMc + multCorC;
    } 
    else{
      Double_t multCorA = -10;

      if(iV0 >= 32 && iV0 < 40)
         multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(33);
      else if(iV0 >= 40 && iV0 < 48)
         multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(41);
      else if(iV0 >= 48 && iV0 < 56)
         multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(49);
      else if(iV0 >= 56 && iV0 < 64)
         multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(57);

      if(multCorA < 0){
         cout<<"Problem with multiplicity in V0A"<<endl;
         continue;
      }
      Qxan += TMath::Cos(fHarmonic*phiV0) * multCorA;
      Qyan += TMath::Sin(fHarmonic*phiV0) * multCorA;

      sumMa = sumMa + multCorA;
    }
  }
}//-------- Get V0 QVect and Multiplicity -----

void AliAnalysisTaskCMEV0::GetZDCCorrectionHist(Int_t run)
{
  if(fListZDNCorr){
    fHCorrectZDNP = (TH2D *) fListZDNCorr->FindObject(Form("fHist_ZDC_All_Wgt_VsCent_Run%d",run));
  }
  else{
    fHCorrectZDNP = new TH2D("fHCorrectZDNP","",100,0,100,4,0,4) ;
    for(int i=1;i<=100;i++){
      for(int j=1;j<=4;j++){
        fHCorrectZDNP->SetBinContent(i,j,1.0);
      }
    }
    //printf("\n\n ********** ZDC Wgt Histograms NotFound ***************\n\n");
    //exit(1);
  }
}

void AliAnalysisTaskCMEV0::GetNUACorrectionHist(Int_t run, TString sfileNUA)
{

  //old method
  /*
  Int_t centBin = -1;
  if(cent<5.0){ centBin = 0;}
  else if(cent>=5.0 && cent<10.0){centBin = 1;}
  else if(cent>=10.0 && cent<40.0){centBin = 2;}
  else if(cent>=40.0) {centBin = 3;}
  if(fListNUACorr){
    fHCorrectNUApos = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pos_Cent%d_Run%d",centBin,run));
    fHCorrectNUAneg = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Neg_Cent%d_Run%d",centBin,run));
  } */

  if(sfileNUA=="New2018" && fListNUACorr){
    for(int i=0;i<5;i++){
      fHCorrectNUApos[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pos_Cent%d_Run%d",i,run));
      fHCorrectNUAneg[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Neg_Cent%d_Run%d",i,run));
    }
    //cout<<"\n sfileNUA = "<<sfileNUA<<" opening Rihan's NUA file.\n "<<endl;
  }
  else if(sfileNUA=="NewR" && fListNUACorr){
    for(int i=0;i<4;i++){
      fHCorrectNUApos[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pos_Cent%d_Run%d",i,run));
      fHCorrectNUAneg[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Neg_Cent%d_Run%d",i,run));
    }
    //cout<<"\n sfileNUA = "<<sfileNUA<<" opening Rihan's NUA file.\n "<<endl;
  }
  else if(sfileNUA=="NewPt"|| sfileNUA=="NewpT"){
    if(fListNUACorr){
      for(int i=0;i<4;i++){
        fHCorrectNUApos[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_pTPhiEta_Pos_Vz%d_Run%d",i,run));
        fHCorrectNUAneg[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_pTPhiEta_Neg_Vz%d_Run%d",i,run));
      }
    }
  }
  else if(sfileNUA=="OldJ") {//from Jacopo's NUA file:
   if(!gGrid){
     TGrid::Connect("alien://");
   }
   if(!fileNUApos){
     fileNUApos = TFile::Open("alien:///alice/cern.ch/user/m/mhaque/calib_files/15oHI_FB768_PosCh_CenPhiEtaWeights_VtxRbR.root");
     mListNUAPos = dynamic_cast<TList*> (fileNUApos->FindObjectAny("CenPhiEta Weights"));
   }
   if(!fileNUAneg){
     fileNUAneg = TFile::Open("alien:///alice/cern.ch/user/m/mhaque/calib_files/15oHI_FB768_NegCh_CenPhiEtaWeights_VtxRbR.root");
     mListNUANeg = dynamic_cast<TList*> (fileNUAneg->FindObjectAny("CenPhiEta Weights"));
   }
   if(mListNUAPos){
    for(int i=0;i<4;i++){
      fHCorrectNUApos[i] = (TH3D *) mListNUAPos->FindObject(Form("CRCQVecPhiHistVtx[%d][%d]",i,run));
    } 
  }
   if(mListNUANeg){
    for(int i=0;i<4;i++){
      fHCorrectNUAneg[i] = (TH3D *) mListNUANeg->FindObject(Form("CRCQVecPhiHistVtx[%d][%d]",i,run));
    }
   }   
  }


  if(!fHCorrectNUApos[0] || !fHCorrectNUAneg[0]){
    printf("\n\n ******** could not open NUA Histograms for run %d, Use Wgt = 1.0 *********\n\n",run);
    if(sfileNUA=="New2018"){
     for(int i=0;i<5;i++){
        fHCorrectNUApos[i] = new TH3D(Form("fHCorrectNUApos_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
        fHCorrectNUAneg[i] = new TH3D(Form("fHCorrectNUAneg_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
        fHCorrectNUApos[i]->SetBinContent(1,1,1,1.0);
        fHCorrectNUAneg[i]->SetBinContent(1,1,1,1.0);
      //exit(1);
     }
    }
    else if(sfileNUA=="OldJ"|| sfileNUA=="NewR"){
     for(int i=0;i<4;i++){
        fHCorrectNUApos[i] = new TH3D(Form("fHCorrectNUApos_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
        fHCorrectNUAneg[i] = new TH3D(Form("fHCorrectNUAneg_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
        fHCorrectNUApos[i]->SetBinContent(1,1,1,1.0);
        fHCorrectNUAneg[i]->SetBinContent(1,1,1,1.0);
      //exit(1);
     }
    }
    else if(sfileNUA=="NewPt"|| sfileNUA=="NewpT"){
     for(int i=0;i<4;i++){
        fHCorrectNUApos[i] = new TH3D(Form("fHCorrectNUApos_cent%d",i),"",1,0,10,1,0,6.284,1,-0.9,0.9); 
        fHCorrectNUAneg[i] = new TH3D(Form("fHCorrectNUAneg_cent%d",i),"",1,0,10,1,0,6.284,1,-0.9,0.9); 
        fHCorrectNUApos[i]->SetBinContent(1,1,1,1.0);
        fHCorrectNUAneg[i]->SetBinContent(1,1,1,1.0);
      //exit(1);
     }
    }
  }
}



void AliAnalysisTaskCMEV0::OpenInfoCalbration(Int_t run, Float_t fHarmonic)
{
  if(!gGrid){
    TGrid::Connect("alien://");
  }

  TFile* foadb = 0;
  //if(!fRemChV0A)
     foadb = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/calibV0HIR.root");
  //else
   //foadb = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/calibV0HIRNoCh46V0A.root");

    if(!foadb){
        printf("OADB V0 calibration file cannot be opened\n");
        return;
    }

    AliOADBContainer* cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorPfpx");
    if(!cont){
        printf("OADB object hMultV0BefCorr is not available in the file\n");
        return;
    }
    if(!(cont->GetObject(run))){
        printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", run);
        return;
    }
    fMultV0 = ((TH1D*) cont->GetObject(run));

    AliOADBContainer* contQxnam = 0;
    if (fHarmonic == 2.)
        contQxnam = (AliOADBContainer*) foadb->Get("fqxa2m");
    else
        contQxnam = (AliOADBContainer*) foadb->Get("fqxa3m");

    if(!contQxnam){
        printf("OADB object fqxanm is not available in the file\n");
        return;
    }
    if(!(contQxnam->GetObject(run))){
        printf("OADB object fqxanm is not available for run %i\n", run);
        return;
    }
    fQxnmV0A = ((TH1D*) contQxnam->GetObject(run));

    AliOADBContainer* contQynam = 0;
    if (fHarmonic == 2.)
        contQynam = (AliOADBContainer*) foadb->Get("fqya2m");
    else if (fHarmonic == 3.)
        contQynam = (AliOADBContainer*) foadb->Get("fqya3m");
    else if (fHarmonic == 4.)
        contQynam = (AliOADBContainer*) foadb->Get("fqya4m");

    if(!contQynam){
        printf("OADB object fqyanm is not available in the file\n");
        return;
    }
    if(!(contQynam->GetObject(run))){
        printf("OADB object fqyanm is not available for run %i\n", run);
        return;
    }
    fQynmV0A = ((TH1D*) contQynam->GetObject(run));

    AliOADBContainer* contQxnas = 0;
    if (fHarmonic == 2.)
        contQxnas = (AliOADBContainer*) foadb->Get("fqxa2s");
    else
        contQxnas = (AliOADBContainer*) foadb->Get("fqxa3s");

    if(!contQxnas){
        printf("OADB object fqxans is not available in the file\n");
        return;
    }
    if(!(contQxnas->GetObject(run))){
        printf("OADB object fqxans is not available for run %i\n", run);
        return;
    }
    fQxnsV0A = ((TH1D*) contQxnas->GetObject(run));

    AliOADBContainer* contQynas = 0;
    if (fHarmonic == 2.)
        contQynas = (AliOADBContainer*) foadb->Get("fqya2s");
    else if (fHarmonic == 3.)
        contQynas = (AliOADBContainer*) foadb->Get("fqya3s");
    else if (fHarmonic == 4.)
        contQynas = (AliOADBContainer*) foadb->Get("fqya4s");

    if(!contQynas){
        printf("OADB object fqyans is not available in the file\n");
        return;
    }
    if(!(contQynas->GetObject(run))){
        printf("OADB object fqyans is not available for run %i\n", run);
        return;
    }
    fQynsV0A = ((TH1D*) contQynas->GetObject(run));



    AliOADBContainer* contQxncm = 0;
    if (fHarmonic == 2.)
        contQxncm = (AliOADBContainer*) foadb->Get("fqxc2m");
    else
        contQxncm = (AliOADBContainer*) foadb->Get("fqxc3m");

    if(!contQxncm){
        printf("OADB object fqxcnm is not available in the file\n");
        return;
    }
    if(!(contQxncm->GetObject(run))){
        printf("OADB object fqxcnm is not available for run %i\n", run);
        return;
    }
    fQxnmV0C = ((TH1D*) contQxncm->GetObject(run));



    AliOADBContainer* contQyncm = 0;
    if (fHarmonic == 2.)
        contQyncm = (AliOADBContainer*) foadb->Get("fqyc2m");
    else if (fHarmonic == 3.)
        contQyncm = (AliOADBContainer*) foadb->Get("fqyc3m");
    else if (fHarmonic == 4.)
        contQyncm = (AliOADBContainer*) foadb->Get("fqyc4m");

    if(!contQyncm){
        printf("OADB object fqyc2m is not available in the file\n");
        return;
    }
    if(!(contQyncm->GetObject(run))){
        printf("OADB object fqyc2m is not available for run %i\n", run);
        return;
    }
    fQynmV0C = ((TH1D*) contQyncm->GetObject(run));


    AliOADBContainer* contQxncs = 0;
    if (fHarmonic == 2.)
        contQxncs = (AliOADBContainer*) foadb->Get("fqxc2s");
    else
        contQxncs = (AliOADBContainer*) foadb->Get("fqxc3s");

    if(!contQxncs){
        printf("OADB object fqxc2s is not available in the file\n");
        return;
    }
    if(!(contQxncs->GetObject(run))){
        printf("OADB object fqxc2s is not available for run %i\n", run);
        return;
    }
    fQxnsV0C = ((TH1D*) contQxncs->GetObject(run));


    AliOADBContainer* contQyncs = 0;
    if (fHarmonic == 2.)
        contQyncs = (AliOADBContainer*) foadb->Get("fqyc2s");
    else if (fHarmonic == 3.)
        contQyncs = (AliOADBContainer*) foadb->Get("fqyc3s");
    else if (fHarmonic == 4.)
        contQyncs = (AliOADBContainer*) foadb->Get("fqyc4s");

    if(!contQyncs){
        printf("OADB object fqycnm is not available in the file\n");
        return;
    }
    if(!(contQyncs->GetObject(run))){
        printf("OADB object fqycns is not available for run %i\n", run);
        return;
    }
    fQynsV0C = ((TH1D*) contQyncs->GetObject(run));
}//------- OADB container ------




double AliAnalysisTaskCMEV0::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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

Bool_t AliAnalysisTaskCMEV0::PileUpMultiVertex(const AliAODEvent* faod)
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



void AliAnalysisTaskCMEV0::GetV0MCorrectionHist(Int_t run)
{
/*
  if(!gGrid){
    TGrid::Connect("alien://");
  }

  TFile *fileV0M  = TFile::Open(fFileV0MCorr);
  if(!fileV0M){
    printf("\n\n ********** File for V0M gain eq. not found ***************\n\n");
  }

  TList *mListV0M = dynamic_cast<TList*> (fileV0M->FindObjectAny("fV0MChWgts"));
*/
  if(fListV0MCorr){
    fHCorrectV0M    = (TH1D *) fListV0MCorr->FindObject(Form("fHistV0Gain_Run%d",run));
    fHAvgerageQnV0A = (TH2D *) fListV0MCorr->FindObject(Form("fHistAvgQnV0A_Run%d",run));
    fHAvgerageQnV0C = (TH2D *) fListV0MCorr->FindObject(Form("fHistAvgQnV0C_Run%d",run));
  }
  else{
    fHCorrectV0M  = new TH1D("fHCorrectV0M","",64,0,64);
    for(int i=1;i<=64;i++){
      fHCorrectV0M->SetBinContent(i,1.0);
    }
    fHAvgerageQnV0A = new TH2D("fHAvgerageQnV0A_empty","<Cos2>,<Sin2>,<Cos3>,<Sin3> V0A",90,0,90,8,0,8);
    fHAvgerageQnV0C = new TH2D("fHAvgerageQnV0C_empty","<Cos2>,<Sin2>,<Cos3>,<Sin3> V0A",90,0,90,8,0,8);
    for(int i=1;i<=90;i++){
      for(int j=1;j<=8;j++){
        fHAvgerageQnV0A->SetBinContent(i,j,0.0);
        fHAvgerageQnV0C->SetBinContent(i,j,0.0);
      }
    }
  }
}




Bool_t AliAnalysisTaskCMEV0::CheckEventIsPileUp(AliAODEvent *faod) {

 Bool_t BisPileup=kFALSE;

 Double_t centrV0M=300;
 Double_t centrCL1=300;
 Double_t centrCL0=300;
 Double_t centrTRK=300;

 if(sDataSet=="2010"||sDataSet=="2011") {
   centrV0M = ((AliVAODHeader*)faod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
   centrCL1 = ((AliVAODHeader*)faod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1");
   centrCL0 = ((AliVAODHeader*)faod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL0");
   centrTRK = ((AliVAODHeader*)faod->GetHeader())->GetCentralityP()->GetCentralityPercentile("TRK");
 }
 else{
   fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
   if(!fMultSelection) {
     printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
     exit(1);
   }
   centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
   centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
   centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
   centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");
 }// 2015,pPb


 if(fRejectPileUp && InputEvent()) {

  //if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;

    if(sDataSet!="2015" && sDataSet!="2015LI" && sDataSet!="2015pPb" && sDataSet!="pPb") {
      if(PileUpMultiVertex(faod)) {
         fPileUpCount->Fill(0.5);
         BisPileup=kTRUE;
      }
      Int_t isPileup = faod->IsPileupFromSPD(3);
      if(isPileup != 0) {
         fPileUpCount->Fill(1.5);
       //BisPileup=kTRUE; //
      }
      if(((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
         fPileUpCount->Fill(2.5);
         BisPileup=kTRUE;
      }
      if(faod->IsIncompleteDAQ())  {
         fPileUpCount->Fill(3.5);
         BisPileup=kTRUE;
      }

    //check vertex consistency
      const AliAODVertex* vtTrc = faod->GetPrimaryVertex();
      const AliAODVertex* vtSPD = faod->GetPrimaryVertexSPD();

      if(vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
        fPileUpCount->Fill(5.5);
        BisPileup=kTRUE;
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
         fPileUpCount->Fill(6.5);
         BisPileup=kTRUE;
      }
      if(fAnalysisUtil->IsPileUpEvent(InputEvent())) {
         fPileUpCount->Fill(7.5);
         BisPileup=kTRUE;
      }
    } ////------ dataset 2010,2011 pile up ---------

    else { //------------ pileup for 2015HI,LI, 2016,pPb data ----------------- 

      if(!fMultSelection->GetThisEventIsNotPileup())
         fPileUpMultSelCount->Fill(0.5);
      if(!fMultSelection->GetThisEventIsNotPileupMV())
         fPileUpMultSelCount->Fill(1.5);
      if(!fMultSelection->GetThisEventIsNotPileupInMultBins())
         fPileUpMultSelCount->Fill(2.5);
      if(!fMultSelection->GetThisEventHasNoInconsistentVertices())
         fPileUpMultSelCount->Fill(3.5);
      if(!fMultSelection->GetThisEventPassesTrackletVsCluster())
         fPileUpMultSelCount->Fill(4.5);
      if(!fMultSelection->GetThisEventIsNotAsymmetricInVZERO())
         fPileUpMultSelCount->Fill(5.5);
      if(!fMultSelection->GetThisEventIsNotIncompleteDAQ())
         fPileUpMultSelCount->Fill(6.5);
      if(!fMultSelection->GetThisEventHasGoodVertex2016())
         fPileUpMultSelCount->Fill(7.5);

         BisPileup=kFALSE;

     //-- pile-up a la Dobrin for LHC15o -----
      if(PileUpMultiVertex(faod)) {
         fPileUpCount->Fill(0.5);
         BisPileup=kTRUE;
      }
      Int_t isPileup = faod->IsPileupFromSPD(3);
      if(isPileup != 0) {
         fPileUpCount->Fill(1.5);
         BisPileup=kTRUE;          
      }
      if(((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
         fPileUpCount->Fill(2.5);
         BisPileup=kTRUE;
      }
      if(faod->IsIncompleteDAQ())  {
         fPileUpCount->Fill(3.5);
         BisPileup=kTRUE;
      }
      if(fabs(centrV0M-centrCL1)> 5.0)  {//default: 7.5
         fPileUpCount->Fill(4.5);
         BisPileup=kTRUE;
      }

     // check vertex consistency
      const AliAODVertex* vtTrc = faod->GetPrimaryVertex();
      const AliAODVertex* vtSPD = faod->GetPrimaryVertexSPD();

      if(vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
        fPileUpCount->Fill(5.5);
        BisPileup=kTRUE;
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
        fPileUpCount->Fill(6.5);
        BisPileup=kTRUE;
      }

      //cuts on tracks
      //Int_t multTrk = 0;
      //Int_t multTrkBefC = 0;
      //Int_t multTrkTOFBefC = 0;
      


      Int_t multTPC = 0;
      Int_t multITSfb96 = 0;
      Int_t multITSfb32 = 0;

      Int_t multTPCFE = 0;
      Int_t multGlobal = 0;
      Int_t multTPCuncut = 0;
      Float_t dEta1,dPt1,dPhi1,ptw1;
      Int_t ptBin,dChrg1;

      Int_t multEsd = ((AliAODHeader*)faod->GetHeader())->GetNumberOfESDTracks();

      const Int_t nTracks = faod->GetNumberOfTracks();


      fRefMultCorr = 0;  //Global Variable
      fRefMultRaw  = 0;   //Global Variable
      Int_t cIndex = 0;

      EvtCent = centrV0M;

      if(EvtCent<5.0) {
	cIndex  = 0; 
      }
      else if(EvtCent>=5.0 && EvtCent<10){
	cIndex  = 1;
      }
      else if(EvtCent>=10.0) {
	cIndex = abs(EvtCent/10.0)  +  1;
      }


      for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
	//AliNanoAODTrack* track = dynamic_cast<AliNanoAODTrack*>(faod->GetTrack(iTracks));
         AliAODTrack* track = (AliAODTrack*)faod->GetTrack(iTracks);
         if(!track)  continue;
	 //---------- old method -----------
        if(track->TestFilterBit(128))
	   multTPC++;
        if(track->TestFilterBit(96))
           multITSfb96++;
	//----------------------------------
         if(track->TestFilterBit(1))  multTPCuncut++;
         if(track->TestFilterBit(32)) multITSfb32++;

	 dEta1  = track->Eta();
	 dPt1   = track->Pt();
	 dPhi1  = track->Phi();
         dChrg1 = track->Charge();

         if(dPt1<0.2 || dPt1>5.0 || TMath::Abs(dEta1)>0.8 || track->GetTPCNcls()<70 || track->GetTPCsignal()<10.0)
            continue;

	 //-------- GetRefMult corrected --------------
	 if(sMCdimension=="1D"){
	   ptBin = fFB_Efficiency_Cent[cIndex]->FindBin(dPt1);
	   ptw1   = 1.0/fFB_Efficiency_Cent[cIndex]->GetBinContent(ptBin);
	 }
	 else if(sMCdimension=="3D"){
	   if(dChrg1>0){
	     ptBin = fFB_Efficiency_Pos[cIndex]->FindBin(dEta1,dPt1,dPhi1); //left here
	     ptw1  = fFB_Efficiency_Pos[cIndex]->GetBinContent(ptBin); 
	   }
	   else if(dChrg1<0){
	     ptBin = fFB_Efficiency_Neg[cIndex]->FindBin(dEta1,dPt1,dPhi1); //left here
	     ptw1  = fFB_Efficiency_Neg[cIndex]->GetBinContent(ptBin); 
	   }
	 }

	 if (ptw1>1e3) ptw1 = 1.0;


	 //Get Corrected mulltiplicity:
	 if(track->TestFilterBit(gFilterBit) && dPt1>=0.2 && TMath::Abs(dEta1) < 0.8) {
	   fRefMultCorr +=  ptw1;
	   fRefMultRaw++;
	 }


         if(track->TestFilterBit(1) && track->Chi2perNDF()>0.2)  multTPCFE++;
         if(!track->TestFilterBit(16) || track->Chi2perNDF()<0.1)   continue;
                

         Double_t b[2]    = {-99., -99.};
         Double_t bCov[3] = {-99., -99., -99.};
                
         AliAODTrack copy(*track);
         Double_t magField = faod->GetMagneticField();
                
         if(magField!=0){     
           if(track->PropagateToDCA(faod->GetPrimaryVertex(), magField, 100., b, bCov) && TMath::Abs(b[0]) < 0.3 && TMath::Abs(b[1]) < 0.3) multGlobal++;    
         }
      }


      /*
      for(Int_t it = 0; it < nTracks; it++) {
        AliAODTrack* aodTrk = (AliAODTrack*)faod->GetTrack(it);
        if(!aodTrk) {
          delete aodTrk;
          continue;
        }
       //if(aodTrk->TestFilterBit(32)){
       //   multTrkBefC++;
       //   if(TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
       //     multTrkTOFBefC++;
       //     if((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2) && (aodTrk->Pt() < 20.))
       //       multTrk++;
       //}
        if(aodTrk->TestFilterBit(128))
           multTPC++;
        if(aodTrk->TestFilterBit(96))
           multITSfb96++;
      } // end of for AOD track loop
    
      */

      Double_t multTPCn      = multTPC;
      Double_t multEsdn      = multEsd;
 
      Double_t multESDTPCDif  = multEsdn  - fPileUpSlopeParm*multTPCn;
      Double_t multTPCGlobDif = multTPCFE - fPileUpSlopeParm*multGlobal;


      //cout<<" multTPC = "<<multTPCn<<"\tESD = "<<multEsdn<<"\tMultRaw = "<<fRefMultRaw<<"\tMultCorr="<<fRefMultCorr<<endl;

      fTPCvsITSfb96->Fill(multITSfb96,multTPC);     //multITSfb96 = FB 96, multTPC = FB 128
      fTPCvsITSfb32->Fill(multITSfb32,multTPC);     //multITSfb32 = FB 32
      fTPCFEvsITSfb96->Fill(multITSfb96,multTPCFE); //multTPCFE = FB 1

      fGlobalTracks->Fill(multGlobal);
      fTPCvsGlobalTrk->Fill(multGlobal,multTPCFE);
      fTPCuncutvsGlobal->Fill(multGlobal,multTPCuncut);



      /*if(multESDTPCDif > (fRejectPileUpTight?700.:15000.)) {
        fPileUpCount->Fill(7.5);
        BisPileup=kTRUE;
        }*/
      /*if(multESDTPCDif > 15000.){ //default: 15000
        fPileUpCount->Fill(7.5);
        BisPileup=kTRUE;
       }*/

       if(fRejectPileUpTight) {
       if(multESDTPCDif > fPileUpConstParm) { 
          fPileUpCount->Fill(7.5);
          BisPileup=kTRUE;
        }
        if(multTPCGlobDif > fPileUpConstParm) { 
          fPileUpCount->Fill(8.5);
          //BisPileup=kTRUE;
        }
        if(BisPileup==kFALSE) {
           if(!fMultSelection->GetThisEventIsNotPileup()) BisPileup=kTRUE;
           if(!fMultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
           if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
           if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
           if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
           if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
           if(!fMultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
           if(BisPileup) fPileUpCount->Fill(9.5);
        }  
      }
    }
  }

 return BisPileup; 
}//-------pile up function ------



Float_t AliAnalysisTaskCMEV0::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign) {

 static const Double_t kPi = TMath::Pi();
 
 Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);

 if(dphistar < -kPi)
   dphistar = -kPi * 2 - dphistar;

 if(dphistar > kPi) // might look funny but is needed
   dphistar = kPi * 2 - dphistar;

 return dphistar;
}







Int_t AliAnalysisTaskCMEV0::GetCurrentRunIndex(Int_t  run) {

 Int_t irun = -1;

 for(int i=0;i<fRunFlag;i++){
  if(run==runNums[i])
   {
    irun = i;
    break;
   }
 }
 if(irun<0) {
   printf("\n ... **WARNING** \n::UserExec() runnumber not listed.\n EXIT..\n");
 }

 return irun;
}//------ GetCurrentRunIndex --------


void AliAnalysisTaskCMEV0::InitializeRunArray(TString sPeriod){

 Int_t runArray_2010[89] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137243, 137236, 137235, 137232, 137231, 137162, 137161};

 Int_t runArray_2011[68] = {167915, 168115, 168460, 169035, 169238, 169859, 170228, 167920, 168310, 168464, 169091, 169411, 169923, 170230, 167985, 168311, 168467, 169094, 169415, 170027, 170268, 167987, 168322, 168511, 169138, 169417, 170081, 170269, 167988, 168325, 168512, 169144, 169835, 170155, 170270, 168069, 168341, 168514, 169145, 169837, 170159, 170306, 168076, 168342, 168777, 169148, 169838, 170163, 170308, 168105, 168361, 168826, 169156, 169846, 170193, 170309, 168107, 168362, 168988, 169160, 169855, 170203, 168108, 168458, 168992, 169167, 169858, 170204};

 Int_t runArray_2015[90] = {246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246676, 246675, 246540, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};

 Int_t runArray_2015LI[13] = {244917, 244918, 244975, 244980, 244982, 244983, 245061, 245064, 245066, 245068, 246390, 246391, 246392};

//Int_t runArray_pPb_13cpass2[14] = {195677, 195675, 195673, 195644, 195635, 195633, 195596, 195593, 195592, 195568, 195567, 195566, 195531, 195529};

 Int_t runArray_pPb_16q_pass1[32] = {265309, 265332, 265334, 265335, 265336, 265338, 265339, 265342, 265343, 265344, 265377, 265378, 265381, 265383, 265384, 265385, 265387, 265388, 265419, 265420, 265421, 265422, 265424, 265425, 265426, 265427, 265435, 265499, 265500, 265501, 265521, 265525}; 

 


 if(sPeriod=="2010"){
  fRunFlag = 89;
  for(int i=0;i<fRunFlag;i++)
    runNums[i] = runArray_2010[i];
 }
 else if(sPeriod=="2011"){
   fRunFlag = 68;                 //<-- 2011 check
  for(int i=0;i<fRunFlag;i++)
    runNums[i] = runArray_2011[i];
 }
 else if(sPeriod=="2015" || sPeriod=="2015PbPb"){
  fRunFlag = 90;
  for(int i=0;i<fRunFlag;i++)
    runNums[i] = runArray_2015[i];
 }
 else if(sPeriod=="2015pPb" || sPeriod=="pPb"){
  fRunFlag = 32;
  for(int i=0;i<fRunFlag;i++)
    runNums[i] = runArray_pPb_16q_pass1[i];
 }
 else if(sPeriod=="2015LI"){
  fRunFlag = 13;
  for(int i=0;i<fRunFlag;i++)
    runNums[i] = runArray_2015LI[i];
 }


 else{
   printf("\n\n ***** Run Number not defined for this data set. *******\n\n Please modify code..\n\n");
   exit(1);
 }

}//------- InitializeRunArray ---------



void AliAnalysisTaskCMEV0::DefineHistograms(){

  fHist_Event_count = new TH1F("fHist_Event_count"," ",20,0,20);
  fHist_Event_count->GetXaxis()->SetBinLabel(1,"Called Exec()");
  fHist_Event_count->GetXaxis()->SetBinLabel(2,"AOD Exist");
  fHist_Event_count->GetXaxis()->SetBinLabel(3,"PileUp");
  fHist_Event_count->GetXaxis()->SetBinLabel(4,"iCentSPD<90");
  fHist_Event_count->GetXaxis()->SetBinLabel(5,"V0 Mult>0");
  fHist_Event_count->GetXaxis()->SetBinLabel(6,"TPC Q!=0");
  fHist_Event_count->GetXaxis()->SetBinLabel(7,"Bad Runs");
  fHist_Event_count->GetXaxis()->SetBinLabel(8,"..TBA..");

  fHist_Event_count->GetXaxis()->SetBinLabel(10,"Final Event");
  fListHistos->Add(fHist_Event_count);

  fPileUpMultSelCount = new TH1F("fPileUpMultSelCount", "fPileUpMultSelCount", 10, 0., 10.);
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(1,"IsNotPileup");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(2,"IsNotPileupMV");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(3,"IsNotPileupInMultBins");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(4,"InconsistentVertices");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(5,"TrackletVsCluster");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(6,"AsymmetricInVZERO");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(7,"IncompleteDAQ");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(8,"GoodVertex2016");
  fListHistos->Add(fPileUpMultSelCount);

  fPileUpCount = new TH1F("fPileUpCount", "fPileUpCount", 12, 0., 12.);
  fPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultiplicityComb08");
  fPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>5.0");
  fPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  Int_t puConst = fPileUpConstParm;
  fPileUpCount->GetXaxis()->SetBinLabel(8,Form("multESDTPCDif>%d",puConst));
  fPileUpCount->GetXaxis()->SetBinLabel(9,Form("multGlobTPCDif>%d",puConst));
  fPileUpCount->GetXaxis()->SetBinLabel(10,"extraPileUpMultSel");
  fListHistos->Add(fPileUpCount);

  //To be implemented fully:
  fTaskConfigParm = new TH1F("fTaskConfigParm", "Connfig Values", 20, 0., 20.);
  fTaskConfigParm->GetXaxis()->SetBinLabel(1,"CentMin");
  fTaskConfigParm->GetXaxis()->SetBinLabel(2,"CentMax");
  fTaskConfigParm->GetXaxis()->SetBinLabel(3,"VzMin");
  fTaskConfigParm->GetXaxis()->SetBinLabel(4,"VzMax");
  fTaskConfigParm->GetXaxis()->SetBinLabel(5,"EtaMin");
  fTaskConfigParm->GetXaxis()->SetBinLabel(6,"EtaMax");
  fTaskConfigParm->GetXaxis()->SetBinLabel(7,"PtMin");
  fTaskConfigParm->GetXaxis()->SetBinLabel(8,"PtMax");
  fTaskConfigParm->GetXaxis()->SetBinLabel(9,"DCAxy");
  fTaskConfigParm->GetXaxis()->SetBinLabel(10,"DCAz");
  fTaskConfigParm->GetXaxis()->SetBinLabel(11,"MCcorrected?");
  fTaskConfigParm->GetXaxis()->SetBinLabel(12,"NUACorrected?");
  fTaskConfigParm->GetXaxis()->SetBinLabel(13,"V0gainCorr?");
  fTaskConfigParm->GetXaxis()->SetBinLabel(14,"ZDCgainCorr?");
  fTaskConfigParm->GetXaxis()->SetBinLabel(15,"harmonic N");
  fTaskConfigParm->GetXaxis()->SetBinLabel(16,"harmonic M");
  fTaskConfigParm->GetXaxis()->SetBinLabel(17,"harmonicPsi");
  fListHistos->Add(fTaskConfigParm);

  if(sDataSet=="2015pPb" || sDataSet=="pPb") {
    fHistVtxZvsRun = new TH2F("fHistVtxZvsRun", " Vz (cm) ",fRunFlag,0,fRunFlag,200, -10, 10.);
    fListHistos->Add(fHistVtxZvsRun);
    fHistVtxXvsRun = new TH2F("fHistVtxXvsRun", " Vx (cm) ",fRunFlag,0,fRunFlag,100, -1, 1.);
    fListHistos->Add(fHistVtxXvsRun);
    fHistVtxYvsRun = new TH2F("fHistVtxYvsRun", " Vy (cm) ",fRunFlag,0,fRunFlag,100, -1, 1.);
    fListHistos->Add(fHistVtxYvsRun);
  }


  //----- CME-ZDN correlator SP method---------
  for(int i=0;i<3;i++){ //ZDN_SP
   //Detector: 0 = V0A, 1 = V0C, 3 = Q-cumulant
    fHist_Corr3p_ZDN_SP_PN[i] = new TProfile2D(Form("fHist_Corr3p_ZDN_SP_PosNeg_Det%d",i+1),"opposit charge correlator",90,0,90,15,0,15,"");
    fHist_Corr3p_ZDN_SP_PN[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_ZDN_SP_PN[i]);
    fHist_Corr3p_ZDN_SP_PP[i] = new TProfile2D(Form("fHist_Corr3p_ZDN_SP_PosPos_Det%d",i+1),"pos-pos charge correlator",90,0,90,15,0,15,"");
    fHist_Corr3p_ZDN_SP_PP[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_ZDN_SP_PP[i]);
    fHist_Corr3p_ZDN_SP_NN[i] = new TProfile2D(Form("fHist_Corr3p_ZDN_SP_NegNeg_Det%d",i+1),"neg-neg charge correlator",90,0,90,15,0,15,"");
    fHist_Corr3p_ZDN_SP_NN[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_ZDN_SP_NN[i]);
  }
  //SP Resolution:
  for(int i=0;i<3;i++){
  //Det: 0 = v0c-v0a, 1 = v0a-TPC, 2 = v0c-TPC, 
    fHist_Reso2n_ZDN_SP_Det[i]  = new TProfile2D(Form("fHist_Reso2n_ZDN_SP_DetComb%d",i+1),"Event plane Resolution",90,0,90,8,0,8,"");
    fHist_Reso2n_ZDN_SP_Det[i]->Sumw2();
    fListHistos->Add(fHist_Reso2n_ZDN_SP_Det[i]);
  }
  //--------------------------------------------

  Double_t centRange[11]    = {0,5,10,20,30,40,50,60,70,80,90};


  fRejectRatioVsCR =  new TProfile2D("fRejectRatioVsCentRun","",10,centRange,fRunFlag,0,fRunFlag);
  fListHistos->Add(fRejectRatioVsCR);



  //----- CME SP method histograms ---------
  for(int i=0;i<2;i++){
    for(int j=0;j<3;j++){
     //Detector: 0 = V0A, 1 = V0C, 3 = Q-cumulant
      fHist_Corr3p_SP_Norm_PN[i][j] = new TProfile(Form("fHist_Corr3p_SP_Norm_PosNeg_Mag%d_Det%d",i,j+1),"opposit charge correlator",10,centRange,"");
      fHist_Corr3p_SP_Norm_PN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_SP_Norm_PN[i][j]);
      fHist_Corr3p_SP_Norm_PP[i][j] = new TProfile(Form("fHist_Corr3p_SP_Norm_PosPos_Mag%d_Det%d",i,j+1),"pos-pos charge correlator",10,centRange,"");
      fHist_Corr3p_SP_Norm_PP[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_SP_Norm_PP[i][j]);
      fHist_Corr3p_SP_Norm_NN[i][j] = new TProfile(Form("fHist_Corr3p_SP_Norm_NegNeg_Mag%d_Det%d",i,j+1),"neg-neg charge correlator",10,centRange,"");
      fHist_Corr3p_SP_Norm_NN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_SP_Norm_NN[i][j]);
    }
    //EP Resolution:
    for(int j=0;j<3;j++){
    //Det: 0 = v0c-v0a, 1 = v0a-TPC, 2 = v0c-TPC, 
      fHist_Reso2n_SP_Norm_Det[i][j]  = new TProfile(Form("fHist_Reso2n_SP_Norm_Mag%d_DetComb%d",i,j+1),"Event plane Resolution",10,centRange,"");
      fHist_Reso2n_SP_Norm_Det[i][j]->Sumw2();
      fListHistos->Add(fHist_Reso2n_SP_Norm_Det[i][j]);
    }
  }

  Int_t gMaxRefMult = 4000;

  if(sDataSet=="2015pPb" || sDataSet=="pPb"){
    gMaxRefMult = 800;
  }

  //----- CME EP method histograms ---------
  for(int i=0;i<2;i++){
    for(int j=0;j<3;j++){
     //Detector: 0 = V0A, 1 = V0C, 3 = Q-cumulant
      fHist_Corr3p_EP_Norm_PN[i][j] = new TProfile(Form("fHist_Corr3p_EP_Norm_PosNeg_Mag%d_Det%d",i,j+1),"opposit charge correlator",10,centRange,"");
      fHist_Corr3p_EP_Norm_PN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EP_Norm_PN[i][j]);
      fHist_Corr3p_EP_Norm_PP[i][j] = new TProfile(Form("fHist_Corr3p_EP_Norm_PosPos_Mag%d_Det%d",i,j+1),"pos-pos charge correlator",10,centRange,"");
      fHist_Corr3p_EP_Norm_PP[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EP_Norm_PP[i][j]);
      fHist_Corr3p_EP_Norm_NN[i][j] = new TProfile(Form("fHist_Corr3p_EP_Norm_NegNeg_Mag%d_Det%d",i,j+1),"neg-neg charge correlator",10,centRange,"");
      fHist_Corr3p_EP_Norm_NN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EP_Norm_NN[i][j]);
    }
    //EP Resolution:
    for(int j=0;j<3;j++){
    //Det: 0 = v0c-v0a, 1 = v0a-TPC, 2 = v0c-TPC, 
      fHist_Reso2n_EP_Norm_Det[i][j]  = new TProfile(Form("fHist_Reso2n_EP_Norm_Mag%d_DetComb%d",i,j+1),"Event plane Resolution",10,centRange,"");
      fHist_Reso2n_EP_Norm_Det[i][j]->Sumw2();
      fListHistos->Add(fHist_Reso2n_EP_Norm_Det[i][j]);
    }

    //Same obsevables in Refmult bins:
    for(int j=0;j<3;j++){
     //Detector: 0 = V0A, 1 = V0C, 3 = Q-cumulant
      fHist_Corr3p_EP_Refm_PN[i][j] = new TProfile(Form("fHist_Corr3p_EP_Refm_PosNeg_Mag%d_Det%d",i,j+1),"opposit charge correlator",800,0,gMaxRefMult,"");
      fHist_Corr3p_EP_Refm_PN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EP_Refm_PN[i][j]);
      fHist_Corr3p_EP_Refm_PP[i][j] = new TProfile(Form("fHist_Corr3p_EP_Refm_PosPos_Mag%d_Det%d",i,j+1),"pos-pos charge correlator",800,0,gMaxRefMult,"");
      fHist_Corr3p_EP_Refm_PP[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EP_Refm_PP[i][j]);
      fHist_Corr3p_EP_Refm_NN[i][j] = new TProfile(Form("fHist_Corr3p_EP_Refm_NegNeg_Mag%d_Det%d",i,j+1),"neg-neg charge correlator",800,0,gMaxRefMult,"");
      fHist_Corr3p_EP_Refm_NN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EP_Refm_NN[i][j]);
    }
    for(int j=0;j<3;j++){
    //Det: 0 = v0c-v0a, 1 = v0a-TPC, 2 = v0c-TPC, 
      fHist_Reso2n_EP_Refm_Det[i][j]  = new TProfile(Form("fHist_Reso2n_EP_Refm_Mag%d_DetComb%d",i,j+1),"Event plane Resolution",800,0,gMaxRefMult,"");
      fHist_Reso2n_EP_Refm_Det[i][j]->Sumw2();
      fListHistos->Add(fHist_Reso2n_EP_Refm_Det[i][j]);
    }
  }







  //two particle correlators:
  for(int i=0;i<2;i++){
    fHist_Corr2p_EP_Norm_PN[i] = new TProfile2D(Form("fHist_Corr2p_EP_Norm_PN_Mag%d",i),"<cos(n(phiA-phiB))>",10,centRange,4,0,4,"");
    fListHistos->Add(fHist_Corr2p_EP_Norm_PN[i]);
    fHist_Corr2p_EP_Norm_PP[i] = new TProfile2D(Form("fHist_Corr2p_EP_Norm_PP_Mag%d",i),"<cos(n(phiA-phiB))>",10,centRange,4,0,4,"");
    fListHistos->Add(fHist_Corr2p_EP_Norm_PP[i]);
    fHist_Corr2p_EP_Norm_NN[i] = new TProfile2D(Form("fHist_Corr2p_EP_Norm_NN_Mag%d",i),"<cos(n(phiA-phiB))>",10,centRange,4,0,4,"");
    fListHistos->Add(fHist_Corr2p_EP_Norm_NN[i]);
  }
  for(int i=0;i<2;i++){
    fHist_Corr2p_EP_Refm_PN[i] = new TProfile2D(Form("fHist_Corr2p_EP_Refm_PN_Mag%d",i),"<cos(n(phiA-phiB))>",800,0,gMaxRefMult,4,0,4,"");
    fListHistos->Add(fHist_Corr2p_EP_Refm_PN[i]);
    fHist_Corr2p_EP_Refm_PP[i] = new TProfile2D(Form("fHist_Corr2p_EP_Refm_PP_Mag%d",i),"<cos(n(phiA-phiB))>",800,0,gMaxRefMult,4,0,4,"");
    fListHistos->Add(fHist_Corr2p_EP_Refm_PP[i]);
    fHist_Corr2p_EP_Refm_NN[i] = new TProfile2D(Form("fHist_Corr2p_EP_Refm_NN_Mag%d",i),"<cos(n(phiA-phiB))>",800,0,gMaxRefMult,4,0,4,"");
    fListHistos->Add(fHist_Corr2p_EP_Refm_NN[i]);
  }






  fHV0AEventPlaneVsCent = new TH2F("fHV0AEventPlaneVsCent","Psi2 from V0A",10,centRange,50,0,3.1415);
  fListHistos->Add(fHV0AEventPlaneVsCent);
  fHV0CEventPlaneVsCent = new TH2F("fHV0CEventPlaneVsCent","Psi2 from V0C",10,centRange,50,0,3.1415);
  fListHistos->Add(fHV0CEventPlaneVsCent);
  fHTPCEventPlaneVsCent = new TH2F("fHTPCEventPlaneVsCent","Psi2 from TPC",10,centRange,50,0,3.1415);
  fListHistos->Add(fHTPCEventPlaneVsCent);




  //TH1::SetDefaultSumw2();

  Char_t name[100], title[100];
 


  //Differential in pT:
  //Double_t pTRange[24] = {0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.33,2.66,3.,3.5,4.,5.,6.,8.,10.};

  Double_t pTRange[21] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0}; 

  for(int i=0;i<2;i++){ 
    for(int j=0;j<6;j++){
      sprintf(name,"fHist_Corr3p_pTSum_EP_V0A_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0A_PN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0A_PN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTSum_EP_V0A_PN[i][j]);

      sprintf(name,"fHist_Corr3p_pTSum_EP_V0A_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0A_PP[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0A_PP[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTSum_EP_V0A_PP[i][j]);

      sprintf(name,"fHist_Corr3p_pTSum_EP_V0A_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0A_NN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0A_NN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTSum_EP_V0A_NN[i][j]);
      //-----v0c----
      sprintf(name,"fHist_Corr3p_pTSum_EP_V0C_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0C_PN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0C_PN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTSum_EP_V0C_PN[i][j]);

      sprintf(name,"fHist_Corr3p_pTSum_EP_V0C_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0C_PP[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0C_PP[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTSum_EP_V0C_PP[i][j]);

      sprintf(name,"fHist_Corr3p_pTSum_EP_V0C_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTSum_EP_V0C_NN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTSum_EP_V0C_NN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTSum_EP_V0C_NN[i][j]);
    }
  }

  for(int i=0;i<2;i++){ 
    for(int j=0;j<6;j++){
      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0A_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0A_PN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0A_PN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0A_PN[i][j]);

      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0A_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0A_PP[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0A_PP[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0A_PP[i][j]);

      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0A_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0A_NN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0A_NN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0A_NN[i][j]);
      //-----v0c----
      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0C_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0C_PN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0C_PN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0C_PN[i][j]);

      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0C_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0C_PP[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0C_PP[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0C_PP[i][j]);

      sprintf(name,"fHist_Corr3p_pTDiff_EP_V0C_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_pTDiff_EP_V0C_NN[i][j] = new TProfile(name,title,20,pTRange,"");
      fHist_Corr3p_pTDiff_EP_V0C_NN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0C_NN[i][j]);
    }
  }
 
  //Double_t EtaRange[9] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6}; //Use this after tests done
  //Now Eta binning: 160,0,1.6 for test

  for(int i=0;i<2;i++){ 
    for(int j=0;j<6;j++){
      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0A_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0A_PN[i][j] = new TProfile(name,title,160,0,1.6,"");
      fHist_Corr3p_EtaDiff_EP_V0A_PN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0A_PN[i][j]);

      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0A_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0A_PP[i][j] = new TProfile(name,title,160,0,1.6,"");
      fHist_Corr3p_EtaDiff_EP_V0A_PP[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0A_PP[i][j]);

      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0A_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0A_NN[i][j] = new TProfile(name,title,160,0,1.6,"");
      fHist_Corr3p_EtaDiff_EP_V0A_NN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0A_NN[i][j]);
      //-----v0c----
      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0C_PN_Mag%d_Cent%d",i,j);
      sprintf(title,"PN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0C_PN[i][j] = new TProfile(name,title,160,0,1.6,""); 
      fHist_Corr3p_EtaDiff_EP_V0C_PN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0C_PN[i][j]);

      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0C_PP_Mag%d_Cent%d",i,j);
      sprintf(title,"PP 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0C_PP[i][j] = new TProfile(name,title,160,0,1.6,"");
      fHist_Corr3p_EtaDiff_EP_V0C_PP[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0C_PP[i][j]);

      sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0C_NN_Mag%d_Cent%d",i,j);
      sprintf(title,"NN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
      fHist_Corr3p_EtaDiff_EP_V0C_NN[i][j] = new TProfile(name,title,160,0,1.6,"");
      fHist_Corr3p_EtaDiff_EP_V0C_NN[i][j]->Sumw2();
      fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0C_NN[i][j]);
    }
  }





  //non-Isotropic terms for correction:
  for(int i=0;i<2;i++){
    fHist_NonIso_SP_PP_Mag0[i] = new TProfile2D(Form("fHist_NonIso_SP_PP_Mag0_Det%d",i),"Non Isotropic terms",10,centRange,10,0,10,"");
    fListHistos->Add(fHist_NonIso_SP_PP_Mag0[i]);
    fHist_NonIso_SP_NN_Mag0[i] = new TProfile2D(Form("fHist_NonIso_SP_NN_Mag0_Det%d",i),"Non Isotropic terms",10,centRange,10,0,10,"");
    fListHistos->Add(fHist_NonIso_SP_NN_Mag0[i]);

    fHist_NonIso_SP_PP_Mag1[i] = new TProfile2D(Form("fHist_NonIso_SP_PP_Mag1_Det%d",i),"Non Isotropic terms",10,centRange,10,0,10,"");
    fListHistos->Add(fHist_NonIso_SP_PP_Mag1[i]);
    fHist_NonIso_SP_NN_Mag1[i] = new TProfile2D(Form("fHist_NonIso_SP_NN_Mag1_Det%d",i),"Non Isotropic terms",10,centRange,10,0,10,"");
    fListHistos->Add(fHist_NonIso_SP_NN_Mag1[i]);
  }










  //=========== Calibtation Histograms etc ==================



  fListCalibs->Add(fHist_Event_count);

  fdPhiFemtoCut = new TH1F("fdPhiFemtoCut"," total comb ",5,0,5);
  fdPhiFemtoCut->GetXaxis()->SetBinLabel(1,"US Total");
  fdPhiFemtoCut->GetXaxis()->SetBinLabel(2,"US Reject");
  fdPhiFemtoCut->GetXaxis()->SetBinLabel(3,"LS Total");
  fdPhiFemtoCut->GetXaxis()->SetBinLabel(4,"LS Reject");
  fListCalibs->Add(fdPhiFemtoCut);

  fVzDistribuion = new TH1F("fVzDistribuion","Vz (cm)",100,-10,10);
  fListCalibs->Add(fVzDistribuion);

  fEventStatvsRun = new TH1F("fEventStatvsRun","Event stat per run",fRunFlag,0,fRunFlag);
  fListCalibs->Add(fEventStatvsRun);

  fCentDistvsRun  = new TH2F("fCentStatvsRun","Cent Dist. vs  run",90,0,90,fRunFlag,0,fRunFlag);
  fListCalibs->Add(fCentDistvsRun);


  fCentV0MvsVzRun =  new TProfile2D("fCentV0MvsVzRun","",100,-10,10,fRunFlag,0,fRunFlag);
  fListCalibs->Add(fCentV0MvsVzRun);

  fCent3pvsVzRun =  new TProfile2D("fCent3pvsVzRun","",100,-10,10,fRunFlag,0,fRunFlag);
  fListCalibs->Add(fCent3pvsVzRun);

  fCentCL1vsVzRun =  new TH3F("fCentCL1vsVzRun","",80,-10,10,fRunFlag,0,fRunFlag,20,-10,10);
  fListCalibs->Add(fCentCL1vsVzRun);


  fRefMultCorrvsRaw = new TH2F("fRefMultCorrvsRaw",Form("Raw vs Corr Mult (FB %d)",gFilterBit),gMaxRefMult,0,gMaxRefMult,gMaxRefMult,0,gMaxRefMult);
  fListCalibs->Add(fRefMultCorrvsRaw);



  fTPCvsGlobalTrk = new TH2F("fTPCvsGlobalTrk","ESDTrk vs TPC(FB128)",400,0,gMaxRefMult,400,0,gMaxRefMult);
  fListCalibs->Add(fTPCvsGlobalTrk);

  fTPCuncutvsGlobal = new TH2F("fTPCuncutvsGlobal","ESDTrk vs ITS(FB96) ",400,0,gMaxRefMult,400,0,gMaxRefMult);
  fListCalibs->Add(fTPCuncutvsGlobal);

  fTPCvsITSfb96 = new TH2F("fTPCvsITSfb96","FB96 vs TPC(FB128)",400,0,gMaxRefMult,400,0,gMaxRefMult);
  fListCalibs->Add(fTPCvsITSfb96);

  fTPCvsITSfb32 = new TH2F("fTPCvsITSfb32","FB32 vs TPC(FB128)",400,0,gMaxRefMult,400,0,gMaxRefMult);
  fListCalibs->Add(fTPCvsITSfb32);

  fTPCFEvsITSfb96 = new TH2F("fTPCFEvsITSfb96","FB96 vs TPC(FB1)",400,0,gMaxRefMult,400,0,gMaxRefMult);
  fListCalibs->Add(fTPCFEvsITSfb96);

  fGlobalTracks  = new TH1F("fGlobalTracks","Global Multiplilcity",gMaxRefMult,0,gMaxRefMult);
  fListCalibs->Add(fGlobalTracks);


  //for debug only, remove after stable code
  hUnderOverBinNUApos = new TH1F("hUnderOverBinNUApos","",90,0,90);
  fListCalibs->Add(hUnderOverBinNUApos);
  hUnderOverBinNUAneg = new TH1F("hUnderOverBinNUAneg","",90,0,90); 
  fListCalibs->Add(hUnderOverBinNUAneg);

  Double_t fCentBinQvect[16] = {0.,2.5,5,10,15,20,25,30,35,40,45,50,60,70,80,90};

  fHEnergyZNCvsCent = new TH2F("fHEnergyZNCvsCent","ZNC Energy vs cent",15,fCentBinQvect,4000,0,200000);
  fListCalibs->Add(fHEnergyZNCvsCent);
  fHEnergyZNAvsCent = new TH2F("fHEnergyZNAvsCent","ZNA Energy vs cent",15,fCentBinQvect,4000,0,200000);
  fListCalibs->Add(fHEnergyZNAvsCent);
  fHEnergyZPCvsCent = new TH2F("fHEnergyZPCvsCent","ZPC Energy vs cent",15,fCentBinQvect,2000,0,50000);
  fListCalibs->Add(fHEnergyZPCvsCent);
  fHEnergyZPAvsCent = new TH2F("fHEnergyZPAvsCent","ZPA Energy vs cent",15,fCentBinQvect,2000,0,50000);
  fListCalibs->Add(fHEnergyZPAvsCent);

  fHEnergyZNCvsCentRun = new TProfile2D("fHEnergyZNCvsCentRun","",90,0,90,fRunFlag,0,fRunFlag);
  fListCalibs->Add(fHEnergyZNCvsCentRun);
  fHEnergyZNAvsCentRun = new TProfile2D("fHEnergyZNAvsCentRun","",90,0,90,fRunFlag,0,fRunFlag);
  fListCalibs->Add(fHEnergyZNAvsCentRun);
  fHEnergyZPCvsCentRun = new TProfile2D("fHEnergyZPCvsCentRun","",90,0,90,fRunFlag,0,fRunFlag);
  fListCalibs->Add(fHEnergyZPCvsCentRun);
  fHEnergyZPAvsCentRun = new TProfile2D("fHEnergyZPAvsCentRun","",90,0,90,fRunFlag,0,fRunFlag);
  fListCalibs->Add(fHEnergyZPAvsCentRun);

  fHEnergyZPCvsZPA = new TH2F("fHEnergyZPCvsZPA","ZNC Energy vs cent",500,0,50000,500,0,50000);
  fListCalibs->Add(fHEnergyZPCvsZPA);
  fHEnergyZNCvsZNA = new TH2F("fHEnergyZNCvsZNA","ZNC Energy vs cent",500,50000,150000,500,50000,150000);
  fListCalibs->Add(fHEnergyZNCvsZNA);

  fHCentBinTrkRecenter = new TH1F("fHCentBinTrkRecenter","centrality Binning",15,fCentBinQvect);




  Int_t magField[2] = {0,1}; //0 = Neg, 1 = Pos

  for(int i=0;i<2;i++){
    fHist_Corr3p_QAEta_SP_V0A_PN[i] = new TProfile(Form("fHist_Corr3p_QAEta_SP_V0A_PN_Mag%d",magField[i]),"PN, Cent 10-20%",80,-0.8,0.8,"");
    fHist_Corr3p_QAEta_SP_V0A_PN[i]->Sumw2();
    fListCalibs->Add(fHist_Corr3p_QAEta_SP_V0A_PN[i]);
    fHist_Corr3p_QAEta_SP_V0A_PP[i] = new TProfile(Form("fHist_Corr3p_QAEta_SP_V0A_PP_Mag%d",magField[i]),"PP, Cent 10-20%",80,-0.8,0.8,"");
    fHist_Corr3p_QAEta_SP_V0A_PP[i]->Sumw2();
    fListCalibs->Add(fHist_Corr3p_QAEta_SP_V0A_PP[i]);
    fHist_Corr3p_QAEta_SP_V0A_NN[i] = new TProfile(Form("fHist_Corr3p_QAEta_SP_V0A_NN_Mag%d",magField[i]),"NN, Cent 10-20%",80,-0.8,0.8,"");
    fHist_Corr3p_QAEta_SP_V0A_NN[i]->Sumw2();
    fListCalibs->Add(fHist_Corr3p_QAEta_SP_V0A_NN[i]);
  }
  for(int i=0;i<2;i++){
    fHist_Corr3p_QAEta_SP_V0C_PN[i] = new TProfile(Form("fHist_Corr3p_QAEta_SP_V0C_PN_Mag%d",magField[i]),"PN, Cent 10-20%",80,-0.8,0.8,"");
    fHist_Corr3p_QAEta_SP_V0C_PN[i]->Sumw2();
    fListCalibs->Add(fHist_Corr3p_QAEta_SP_V0C_PN[i]);
    fHist_Corr3p_QAEta_SP_V0C_PP[i] = new TProfile(Form("fHist_Corr3p_QAEta_SP_V0C_PP_Mag%d",magField[i]),"PP, Cent 10-20%",80,-0.8,0.8,"");
    fHist_Corr3p_QAEta_SP_V0C_PP[i]->Sumw2();
    fListCalibs->Add(fHist_Corr3p_QAEta_SP_V0C_PP[i]);
    fHist_Corr3p_QAEta_SP_V0C_NN[i] = new TProfile(Form("fHist_Corr3p_QAEta_SP_V0C_NN_Mag%d",magField[i]),"NN, Cent 10-20%",80,-0.8,0.8,"");
    fHist_Corr3p_QAEta_SP_V0C_NN[i]->Sumw2();
    fListCalibs->Add(fHist_Corr3p_QAEta_SP_V0C_NN[i]);
  }

  fEtaBinFinderForQA = new TH1F("fEtaBinFinderForQA","",80,-0.8,0.8);
  fListCalibs->Add(fEtaBinFinderForQA);






  Int_t suffixEta[4] = {1,0,1,0}; //{"Pos","Neg","Pos","Neg"};
  Int_t suffixVz[4]  = {1,1,0,0}; //{"Pos","Pos","Neg","Neg"};

  if(bFillAvgTPCQn){
    for(int i=0;i<4;i++){
      fHCos1nPosChEtaVz[i] = new TProfile2D(Form("fHCos1nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos1n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos1nPosChEtaVz[i]);
      fHSin1nPosChEtaVz[i] = new TProfile2D(Form("fHSin1nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin1n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin1nPosChEtaVz[i]);
      fHCos1nNegChEtaVz[i] = new TProfile2D(Form("fHCos1nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos1n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos1nNegChEtaVz[i]);
      fHSin1nNegChEtaVz[i] = new TProfile2D(Form("fHSin1nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin1n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin1nNegChEtaVz[i]);
      /*
      fHCos2nPosChEtaVz[i] = new TProfile2D(Form("fHCos2nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos2nPosChEtaVz[i]);
      fHSin2nPosChEtaVz[i] = new TProfile2D(Form("fHSin2nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin2nPosChEtaVz[i]);
      fHCos2nNegChEtaVz[i] = new TProfile2D(Form("fHCos2nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos2nNegChEtaVz[i]);
      fHSin2nNegChEtaVz[i] = new TProfile2D(Form("fHSin2nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin2nNegChEtaVz[i]);

      fHCos3nPosChEtaVz[i] = new TProfile2D(Form("fHCos3nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos3n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos3nPosChEtaVz[i]);
      fHSin3nPosChEtaVz[i] = new TProfile2D(Form("fHSin3nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin3n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin3nPosChEtaVz[i]);
      fHCos3nNegChEtaVz[i] = new TProfile2D(Form("fHCos3nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos3n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos3nNegChEtaVz[i]);
      fHSin3nNegChEtaVz[i] = new TProfile2D(Form("fHSin3nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin3n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin3nNegChEtaVz[i]);
       */
      fHCos4nPosChEtaVz[i] = new TProfile2D(Form("fHCos4nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos4n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos4nPosChEtaVz[i]);
      fHSin4nPosChEtaVz[i] = new TProfile2D(Form("fHSin4nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin4n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin4nPosChEtaVz[i]);
      fHCos4nNegChEtaVz[i] = new TProfile2D(Form("fHCos4nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos4n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos4nNegChEtaVz[i]);
      fHSin4nNegChEtaVz[i] = new TProfile2D(Form("fHSin4nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin4n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin4nNegChEtaVz[i]);

      //Double weighted <..>
      fHCos2nDWPosChEtaVz[i] = new TProfile2D(Form("fHCos2nDWPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*ptw*cos2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos2nDWPosChEtaVz[i]);
      fHSin2nDWPosChEtaVz[i] = new TProfile2D(Form("fHSin2nDWPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*ptw*cos2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin2nDWPosChEtaVz[i]);
      fHCos2nDWNegChEtaVz[i] = new TProfile2D(Form("fHCos2nDWNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*ptw*cos2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos2nDWNegChEtaVz[i]);
      fHSin2nDWNegChEtaVz[i] = new TProfile2D(Form("fHSin2nDWNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*ptw*cos2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin2nDWNegChEtaVz[i]);
    }
  }

  //fill <Q2x> and <Q2y> track by track: for Test
    for(int i=0;i<4;i++){
      fHCos2nPosChEtaVz[i] = new TProfile2D(Form("fHCos2nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos2nPosChEtaVz[i]);
      fHSin2nPosChEtaVz[i] = new TProfile2D(Form("fHSin2nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin2nPosChEtaVz[i]);
      fHCos2nNegChEtaVz[i] = new TProfile2D(Form("fHCos2nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos2nNegChEtaVz[i]);
      fHSin2nNegChEtaVz[i] = new TProfile2D(Form("fHSin2nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin2n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin2nNegChEtaVz[i]);

      fHCos3nPosChEtaVz[i] = new TProfile2D(Form("fHCos3nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos3n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos3nPosChEtaVz[i]);
      fHSin3nPosChEtaVz[i] = new TProfile2D(Form("fHSin3nPosChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin3n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin3nPosChEtaVz[i]);
      fHCos3nNegChEtaVz[i] = new TProfile2D(Form("fHCos3nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*cos3n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHCos3nNegChEtaVz[i]);
      fHSin3nNegChEtaVz[i] = new TProfile2D(Form("fHSin3nNegChEta%dVz%d",suffixEta[i],suffixVz[i]),"<ptw*sin3n>",90,0,90,fRunFlag,0,fRunFlag);
      fListCalibs->Add(fHSin3nNegChEtaVz[i]);
    }



  //Average Multiplicity, same-sign, opposite-sign pair vs Centrality 1% :

  fAvgMultCentRun = new TProfile2D("fAvgMultCentRun","<Mult> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListCalibs->Add(fAvgMultCentRun);
  fAvgWgtMultCentRun = new TProfile2D("fAvgWgtMultCentRun","<wgt*Mult> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListCalibs->Add(fAvgWgtMultCentRun);
  fAvgPOIposCentRun = new TProfile2D("fAvgPOIposCentRun","<wgt*POIs> ch-pos vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListCalibs->Add(fAvgPOIposCentRun);
  fAvgPOInegCentRun = new TProfile2D("fAvgPOInegCentRun","<wgt*POIs> ch-neg vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListCalibs->Add(fAvgPOInegCentRun);
  fAvgPOIPPCentRun = new TProfile2D("fAvgPOIPPCentRun","<wgt*POIs> Pos-Pos vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListCalibs->Add(fAvgPOIPPCentRun);
  fAvgPOINNCentRun = new TProfile2D("fAvgPOINNCentRun","<wgt*POIs> Neg-Neg vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListCalibs->Add(fAvgPOINNCentRun);
  fAvgPOIOSCentRun = new TProfile2D("fAvgPOIOSCentRun","<wgt*POIs> oppo.-sign vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListCalibs->Add(fAvgPOIOSCentRun);



  //Store CME run by Run: Not used for huge memory usage.
  /*
  for(int i=0;i<2;i++){//0=V0A,1=V0C
    fHist_Corr3p_vsRun_EP_PN[i] = new TProfile2D(Form("fHist_Corr3p_vsRun_EP_PN_Det%d",i+1),"opposit charge correlator",90,0,90,fRunFlag,0,fRunFlag,"");
    fListHistos->Add(fHist_Corr3p_vsRun_EP_PN[i]);
    fHist_Corr3p_vsRun_EP_PP[i] = new TProfile2D(Form("fHist_Corr3p_vsRun_EP_PP_Det%d",i+1),"PosPos charge correlator",90,0,90,fRunFlag,0,fRunFlag,"");
    fListHistos->Add(fHist_Corr3p_vsRun_EP_PP[i]);
    fHist_Corr3p_vsRun_EP_NN[i] = new TProfile2D(Form("fHist_Corr3p_vsRun_EP_NN_Det%d",i+1),"NegNeg charge correlator",90,0,90,fRunFlag,0,fRunFlag,"");
    fListHistos->Add(fHist_Corr3p_vsRun_EP_NN[i]);
  }*/








  //======== Following are all type of NUA written in List: fListNUAHist ================


  //---------- V0 histograms --------- 
  fV0MultChVsRun = new TProfile2D("fV0MultChVsRun","1-32 V0C, 33-64 V0A",64,0,64,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fV0MultChVsRun);

  fV0AQ2xVsCentRun = new TProfile2D("fV0ACos2nVsCentRun","<Cos2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fV0AQ2xVsCentRun);
  fV0AQ2yVsCentRun = new TProfile2D("fV0ASin2nVsCentRun","<Sin2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fV0AQ2yVsCentRun);
  fV0CQ2xVsCentRun = new TProfile2D("fV0CCos2nVsCentRun","<Cos2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fV0CQ2xVsCentRun);
  fV0CQ2yVsCentRun = new TProfile2D("fV0CSin2nVsCentRun","<Sin2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fV0CQ2yVsCentRun);

  fV0AQ3xVsCentRun = new TProfile2D("fV0ACos3nVsCentRun","<Cos3> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fV0AQ3xVsCentRun);
  fV0AQ3yVsCentRun = new TProfile2D("fV0ASin3nVsCentRun","<Sin3> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fV0AQ3yVsCentRun);
  fV0CQ3xVsCentRun = new TProfile2D("fV0CCos3nVsCentRun","<Cos3> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fV0CQ3xVsCentRun);
  fV0CQ3yVsCentRun = new TProfile2D("fV0CSin3nVsCentRun","<Sin3> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fV0CQ3yVsCentRun);
  //-----------------------------------


  fTPCQ2xVsCentRun = new TProfile2D("fTPCCos2nVsCentRun","<Cos2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fTPCQ2xVsCentRun);
  fTPCQ2yVsCentRun = new TProfile2D("fTPCSin2nVsCentRun","<Sin2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fTPCQ2yVsCentRun);

  fTPCQ3xVsCentRun = new TProfile2D("fTPCCos3nVsCentRun","<Cos3> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fTPCQ3xVsCentRun);
  fTPCQ3yVsCentRun = new TProfile2D("fTPCSin3nVsCentRun","<Sin3> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"");
  fListNUAHist->Add(fTPCQ3yVsCentRun);


  for(int i=0;i<10;i++){
    sprintf(name,"fHistChPosvsEtaPtRun_Cent%d",i);
    sprintf(title,"Pos Ch, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHistChPosvsEtaPtRun[i]  = new TH3F(name,title,25,0.2,5.2,16,-0.8,0.8,fRunFlag,0,fRunFlag);
    fListNUAHist->Add(fHistChPosvsEtaPtRun[i]);
    sprintf(name,"fHistChNegvsEtaPtRun_Cent%d",i);
    sprintf(title,"Neg Ch, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHistChNegvsEtaPtRun[i]  = new TH3F(name,title,25,0.2,5.2,16,-0.8,0.8,fRunFlag,0,fRunFlag);
    fListNUAHist->Add(fHistChNegvsEtaPtRun[i]);
  }

  
  Int_t gCentForNUA[6] = {0,5,10,20,40,90};

  if(bFillEtaPhiNUA) {
   for(int i=0;i<5;i++){
    for(int j=0;j<fRunFlag;j++){
      sprintf(name,"fHistEtaPhiVz_Pos_Cent%d_Run%d",i,runNums[j]);
      sprintf(title,"eta,phi,Vz Pos Cent%d-%d%%",gCentForNUA[i],gCentForNUA[i+1]);
      fHist3DEtaPhiVz_Pos_Run[i][j] = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
      fListNUAHist->Add(fHist3DEtaPhiVz_Pos_Run[i][j]);

      sprintf(name,"fHistEtaPhiVz_Neg_Cent%d_Run%d",i,runNums[j]);
      sprintf(title,"eta,phi,Vz Pos Cent%d-%d%%",gCentForNUA[i],gCentForNUA[i+1]);
      fHist3DEtaPhiVz_Neg_Run[i][j] = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
      fListNUAHist->Add(fHist3DEtaPhiVz_Neg_Run[i][j]);
    }
   }
  }

  fVzBinFinderForNUA = new TH1F("fVzBinFinderForNUA","",4,-10,10); //only needed for pT dependent NUA.
  fListNUAHist->Add(fVzBinFinderForNUA);

  /*
  //  pT Dependent NUA correction:
  Double_t pTbinNUA[16] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.6,3.0,3.5,4.0,5.0};
  Double_t etabinNUA[9] = {-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8}; 
  Double_t phibinNUA[51] = {0.,};
  Double_t fPhi = 6.2832/50;
  for(int i=0;i<51;i++){
    phibinNUA[i] = fPhi*i;
  }

  if(bFillEtaPhiNUA) {
   for(int i=0;i<5;i++){
    for(int j=0;j<fRunFlag;j++){
      sprintf(name,"fHistEtaPhiPt_Pos_Vz%d_Run%d",i,runNums[j]);
      sprintf(title,"eta,phi,Pt Pos Vz%d",i);
      fHist3DEtaPhiVz_Pos_Run[i][j] = new TH3F(name,title,15,pTbinNUA,50,phibinNUA,8,etabinNUA); 
      fListNUAHist->Add(fHist3DEtaPhiVz_Pos_Run[i][j]);

      sprintf(name,"fHistEtaPhiPt_Neg_Vz%d_Run%d",i,runNums[j]);
      sprintf(title,"eta,phi,Pt Pos Vz%d",i);
      fHist3DEtaPhiVz_Neg_Run[i][j] = new TH3F(name,title,15,pTbinNUA,50,phibinNUA,8,etabinNUA); 
      fListNUAHist->Add(fHist3DEtaPhiVz_Neg_Run[i][j]);
    }
   }
  }  */



}


