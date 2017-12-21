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
// 
// Author: Rihan Haque (mhaque@cern.ch)
///////////////////////////////////////////////


#include "Riostream.h" //needed as include

class TFile;
class TList;
class AliAnalysisTaskSE;

#include "TProfile.h"  //needed as include
#include "TProfile2D.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "AliMultSelection.h"
#include "AliVVertex.h"

#include "AliAnalysisManager.h"
#include "AliFlowEvent.h"
#include "AliFlowEventSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliLog.h"
#include "AliAnalysisTaskZDCGainEq.h"
#include "AliFlowVector.h"

using std::endl;
using std::cout;


ClassImp(AliAnalysisTaskZDCGainEq)

//________________________________________________________________________
AliAnalysisTaskZDCGainEq::AliAnalysisTaskZDCGainEq(const char *name) :
  AliAnalysisTaskSE(name),
  fEvent(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHistos(NULL),
  fListZDCQxy(NULL),
  fListZDCWgt(NULL),
  fListDummy1(NULL),
  fListHijing(NULL),
  fListSubRun(NULL),
  fRejectPileUpTight(kFALSE),
  fRejectPileUp(kFALSE),
  bFillCosSin(kFALSE),
  bFillZDCQAon(kFALSE),
  bRunAveragedQn(kFALSE),
  bApplyRecent(kFALSE),
  bCentCutShift(kFALSE),
  fHarmonic(2),
  frunflag(0),
  fievent(0),
  vxBin(0),
  vyBin(0),
  vzBin(0),
  fcheckOnce(0),
  fOldRunNum(0),
  fHist_Event_count(NULL),
  fPileUpMultSelCount(NULL),
  fPileUpCount(NULL),
  fHist_ChanWgt_ZDCC(NULL),
  fHist_ChanWgt_ZDCA(NULL),
  fHist_Vx_ArrayFinder(NULL),
  fHist_Vy_ArrayFinder(NULL),
  fHist_Vz_ArrayFinder(NULL),
  fHist_Task_config(NULL),
  fHist_Cent_woZDCcut(NULL),
  fHist_Cent_wiZDCcut(NULL),
  fHist_CutParameters(NULL),
  fHist_Psi1_ZDCC_wGainCorr(NULL),
  fHist_Psi1_ZDCA_wGainCorr(NULL),
  fHist_Psi1_ZDCC_wRectCorr(NULL),
  fHist_Psi1_ZDCA_wRectCorr(NULL),
  fHist_Psi1_ZDCC_wCorrFull(NULL),
  fHist_Psi1_ZDCA_wCorrFull(NULL),
  fHist_Psi1_ZDCC_RunByRun(NULL),
  fHist_Psi1_ZDCA_RunByRun(NULL),
  fHist_v2xV1_ZDN_Norm_All(NULL),
  fHist_v2xV1_ZDN_Refm_All(NULL),
  fHist_v2xV1_ZDN_Cent_All(NULL),
  fHist_v3xV1_ZDN_Norm_Comb1(NULL),
  fHist_v3xV1_ZDN_Norm_Comb2(NULL),
  fHist_v3xV1_ZDN_Cent_Comb1(NULL),
  fHist_v3xV1_ZDN_Cent_Comb2(NULL),
  fHist_v4xV1_ZDN_Norm_Comb1(NULL),
  fHist_v4xV1_ZDN_Cent_Comb1(NULL),
  fHist_ZDN_resol_Norm_All(NULL),     
  fHist_ZDN_resol_Cent_All(NULL),
  fHist_ZDN_resol_Refm_All(NULL),
  fHist_ZDN_resol_Norm_XX(NULL),
  fHist_ZDN_resol_Norm_YY(NULL),
  fHist_ZDN_resol_Cent_XX(NULL),
  fHist_ZDN_resol_Cent_YY(NULL),
  fHist_Vx_vs_runnum(NULL),
  fHist_Vy_vs_runnum(NULL),
  fHist_Vz_vs_runnum(NULL),
  fWeight_Cent(NULL),
  fHist_Vxy_RunAveraged(NULL),
  fHist_Event_counter_vRun(NULL),
  fHist_PsiSumAC_woCorr(NULL),
  fHist_PsiSumAC_wiCorr(NULL),
  fHist_ZeroQnXY_Conter(NULL),
  fHist_ZDCC_AvgQx_VsRun(NULL),
  fHist_ZDCC_AvgQy_VsRun(NULL),
  fHist_ZDCA_AvgQx_VsRun(NULL),
  fHist_ZDCA_AvgQy_VsRun(NULL),
  fHist_ZDC_AvgXXminusYY_VsRun(NULL),
  fHist_ZDC_AvgXYplusXY_VsRun(NULL),
  fHist_ZDCC_AvgQx_VsCR(NULL),
  fHist_ZDCC_AvgQy_VsCR(NULL),
  fHist_ZDCA_AvgQx_VsCR(NULL),
  fHist_ZDCA_AvgQy_VsCR(NULL),
  fHist_PsiSumAC_ZeroQn(NULL),
  fHist_PsiZDCA_ZeroQn(NULL),
  fHist_PsiZDCC_ZeroQn(NULL),
  fHist_ZDC_dTermXXYY_VsRun(NULL),
  fHist_ZDC_dTermXYXY_VsRun(NULL),
  fHist_Corr3p_ZDN_Norm_PN(NULL),
  fHist_Corr3p_ZDN_Norm_PP(NULL),
  fHist_Corr3p_ZDN_Norm_NN(NULL),
  fHist_Corr3p_ZDN_Cent_PN(NULL),
  fHist_Corr3p_ZDN_Cent_PP(NULL),
  fHist_Corr3p_ZDN_Cent_NN(NULL),
  fHist_Reso2EP_TPC_Norm(NULL),
  fHist_Reso2EP_TPC_Cent(NULL),
  fHist_ZDC_dTermXXYY_VsCR(NULL),
  fHist_ZDC_dTermXYXY_VsCR(NULL),
  fHist_NormalCentralityBins(NULL),
  fDataSet("2010"),
  fAnalysisSet("DoGainEq"),
  sCentEstimator("V0")
{
  for(int i=0;i<90;i++){
    runNums[i] = 0;
    fHist_ZDCA_En_Run[i]  = NULL;
    fHist_ZDCC_En_Run[i]  = NULL;

    for(int j=0;j<10;j++){
     fHist_znCx_V0_VxVy[i][j] = NULL;
     fHist_znCy_V0_VxVy[i][j] = NULL;
     fHist_znAx_V0_VxVy[i][j] = NULL;
     fHist_znAy_V0_VxVy[i][j] = NULL;
    }     
  }

  for(int i=0;i<4;i++){
    fHist_Qx_wiCorr_RunByRun[i] = NULL;
    fHist_Qx_Trig_woCorr[i] = NULL;
    fHist_XX_Trig_woCorr[i] = NULL;
    fHist_Qx_Trig_wiCorr[i] = NULL;
    fHist_XX_Trig_wiCorr[i] = NULL;
    for(int j=0;j<5;j++){
     fHist_Qx_vs_Obs_woCorr[i][j] = NULL;
     fHist_XX_vs_Obs_woCorr[i][j] = NULL;
     fHist_Qx_vs_Obs_wiCorr[i][j] = NULL;
     fHist_XX_vs_Obs_wiCorr[i][j] = NULL;
    }
  }
  for(int i=0;i<20;i++){
    fHist_Recenter_ZDCCx[i] = NULL;
    fHist_Recenter_ZDCCy[i] = NULL;
    fHist_Recenter_ZDCAx[i] = NULL;
    fHist_Recenter_ZDCAy[i] = NULL;
  }
  for(int i=0;i<20;i++){
    fHist_ZDCC_En_CommonCh[i] = NULL;
    fHist_ZDCA_En_CommonCh[i] = NULL;
  }

  for(int i=0;i<2;i++){
    VxCut[i] = 0;
    VyCut[i] = 0;
    VzCut[i] = 0;
  }

  for(int i=0;i<10;i++){
    fFB_Efficiency_Cent[i] = NULL;
    fHist_v2xV1_ZDN_pTDiff_All[i] = NULL;
    fHist_v4xV1_ZDN_pTDiff_All[i] = NULL;
    fHist_v3xV1_ZDN_EtaDiff_Comb1[i]= NULL;
    fHist_v3xV1_ZDN_EtaDiff_Comb2[i]= NULL;

    for(int j=0;j<4;j++){
      fHist_v1xV1_ZDN_EtaDiff[j][i] = NULL;
      fHist_v1xV1_ZDN_pTDiff[j][i]  = NULL;
    }
  }

  for(int i=0;i<4;i++){
    fHist_v2xV1_ZDN_Norm_Sep[i] = NULL;
    fHist_v2xV1_ZDN_Cent_Sep[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_ZDN_resol_Norm_Sep[i] = NULL;
    fHist_ZDN_resol_Cent_Sep[i] = NULL;
  }
  for(int i=0;i<4;i++){ 
    //here
    fHist_ZDCC_AvgCosNPsi[i]  = NULL;
    fHist_ZDCC_AvgSinNPsi[i]  = NULL;
    fHist_ZDCA_AvgCosNPsi[i]  = NULL;
    fHist_ZDCA_AvgSinNPsi[i]  = NULL;
    fHist_ZDC_AvgCosNPsiAC[i] = NULL;
    fHist_ZDC_AvgSinNPsiAC[i] = NULL;
  }
  for(int i=0;i<4;i++){
    fHist_ZDCC_AvgCos_VsRun[i]  = NULL;
    fHist_ZDCC_AvgSin_VsRun[i]  = NULL;
    fHist_ZDCA_AvgCos_VsRun[i]  = NULL;
    fHist_ZDCA_AvgSin_VsRun[i]  = NULL;
    fHist_ZDC_AvgCosPsiSum_VsRun[i] = NULL;
    fHist_ZDC_AvgSinPsiSum_VsRun[i] = NULL;
    fHist_ZDC_AvgCosPsiDif_VsRun[i] = NULL;
    fHist_ZDC_AvgSinPsiDif_VsRun[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_XXYY_vs_Cent_woCorr[i] = NULL;
    fHist_XXYY_vs_Cent_wiCorr[i] = NULL;
  }

  for(int i=0;i<4;i++){
   for(int j=0;j<6;j++){
    fHist_XX_vs_QnC_2DwoCorr_PosMag[i][j] = NULL;
    fHist_XX_vs_QnA_2DwoCorr_PosMag[i][j] = NULL;
    fHist_XX_vs_QnC_2DwoCorr_NegMag[i][j] = NULL;
    fHist_XX_vs_QnA_2DwoCorr_NegMag[i][j] = NULL;
   }
  }

  for(int i=0;i<90;i++){
   fHist_VZERO_Mult_vsRun[i] = NULL;
  }


  DefineInput(1, AliFlowEventSimple::Class()); // Input slot #1 works with an AliFlowEventSimple
  DefineInput(2, AliFlowEventSimple::Class()); // Input slot #2 for ZDC flow event


  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());

  //fDataSet="2010";
  //fAnalysisSet="DoGainEq";
  //sCentEstimator="V0";

 //fTotalQvector = new TString("QaQb");         // "QaQb" (means Qa+Qb), "Qa"  or "Qb"

}//-------------constructor-----------------

//________________________________________________
AliAnalysisTaskZDCGainEq::AliAnalysisTaskZDCGainEq() :
  AliAnalysisTaskSE(),
  fEvent(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHistos(NULL),
  fListZDCQxy(NULL),
  fListZDCWgt(NULL),
  fListDummy1(NULL),
  fListHijing(NULL),
  fListSubRun(NULL),
  fRejectPileUpTight(kFALSE),
  fRejectPileUp(kFALSE),
  bFillCosSin(kFALSE),
  bFillZDCQAon(kFALSE),
  bRunAveragedQn(kFALSE),
  bApplyRecent(kFALSE),
  bCentCutShift(kFALSE),
  fHarmonic(2),
  frunflag(0),
  fievent(0),
  vxBin(0),
  vyBin(0),
  vzBin(0),
  fcheckOnce(0),
  fOldRunNum(0),
  fHist_Event_count(NULL),
  fPileUpMultSelCount(NULL),
  fPileUpCount(NULL),
  fHist_ChanWgt_ZDCC(NULL),
  fHist_ChanWgt_ZDCA(NULL),
  fHist_Vx_ArrayFinder(NULL),
  fHist_Vy_ArrayFinder(NULL),
  fHist_Vz_ArrayFinder(NULL),
  fHist_Task_config(NULL),
  fHist_Cent_woZDCcut(NULL),
  fHist_Cent_wiZDCcut(NULL),
  fHist_CutParameters(NULL),
  fHist_Psi1_ZDCC_wGainCorr(NULL),
  fHist_Psi1_ZDCA_wGainCorr(NULL),
  fHist_Psi1_ZDCC_wRectCorr(NULL),
  fHist_Psi1_ZDCA_wRectCorr(NULL),
  fHist_Psi1_ZDCC_wCorrFull(NULL),
  fHist_Psi1_ZDCA_wCorrFull(NULL),
  fHist_Psi1_ZDCC_RunByRun(NULL),
  fHist_Psi1_ZDCA_RunByRun(NULL),
  fHist_v2xV1_ZDN_Norm_All(NULL),
  fHist_v2xV1_ZDN_Refm_All(NULL),
  fHist_v2xV1_ZDN_Cent_All(NULL),
  fHist_v3xV1_ZDN_Norm_Comb1(NULL),
  fHist_v3xV1_ZDN_Norm_Comb2(NULL),
  fHist_v3xV1_ZDN_Cent_Comb1(NULL),
  fHist_v3xV1_ZDN_Cent_Comb2(NULL),
  fHist_v4xV1_ZDN_Norm_Comb1(NULL),
  fHist_v4xV1_ZDN_Cent_Comb1(NULL),
  fHist_ZDN_resol_Norm_All(NULL),     
  fHist_ZDN_resol_Cent_All(NULL),
  fHist_ZDN_resol_Refm_All(NULL),
  fHist_ZDN_resol_Norm_XX(NULL),
  fHist_ZDN_resol_Norm_YY(NULL),
  fHist_ZDN_resol_Cent_XX(NULL),
  fHist_ZDN_resol_Cent_YY(NULL),
  fHist_Vx_vs_runnum(NULL),
  fHist_Vy_vs_runnum(NULL),
  fHist_Vz_vs_runnum(NULL),
  fWeight_Cent(NULL),
  fHist_Vxy_RunAveraged(NULL),
  fHist_Event_counter_vRun(NULL),
  fHist_PsiSumAC_woCorr(NULL),
  fHist_PsiSumAC_wiCorr(NULL),
  fHist_ZeroQnXY_Conter(NULL),
  fHist_ZDCC_AvgQx_VsRun(NULL),
  fHist_ZDCC_AvgQy_VsRun(NULL),
  fHist_ZDCA_AvgQx_VsRun(NULL),
  fHist_ZDCA_AvgQy_VsRun(NULL),
  fHist_ZDC_AvgXXminusYY_VsRun(NULL),
  fHist_ZDC_AvgXYplusXY_VsRun(NULL),
  fHist_ZDCC_AvgQx_VsCR(NULL),
  fHist_ZDCC_AvgQy_VsCR(NULL),
  fHist_ZDCA_AvgQx_VsCR(NULL),
  fHist_ZDCA_AvgQy_VsCR(NULL),
  fHist_PsiSumAC_ZeroQn(NULL),
  fHist_PsiZDCA_ZeroQn(NULL),
  fHist_PsiZDCC_ZeroQn(NULL),
  fHist_ZDC_dTermXXYY_VsRun(NULL),
  fHist_ZDC_dTermXYXY_VsRun(NULL),
  fHist_Corr3p_ZDN_Norm_PN(NULL),
  fHist_Corr3p_ZDN_Norm_PP(NULL),
  fHist_Corr3p_ZDN_Norm_NN(NULL),
  fHist_Corr3p_ZDN_Cent_PN(NULL),
  fHist_Corr3p_ZDN_Cent_PP(NULL),
  fHist_Corr3p_ZDN_Cent_NN(NULL),
  fHist_Reso2EP_TPC_Norm(NULL),
  fHist_Reso2EP_TPC_Cent(NULL),
  fHist_ZDC_dTermXXYY_VsCR(NULL),
  fHist_ZDC_dTermXYXY_VsCR(NULL),
  fHist_NormalCentralityBins(NULL),
  fDataSet("2010"),
  fAnalysisSet("DoGainEq"),
  sCentEstimator("V0")
{
  for(int i=0;i<90;i++){
    runNums[i] = 0;
    fHist_ZDCA_En_Run[i]  = NULL;
    fHist_ZDCC_En_Run[i]  = NULL;

    for(int j=0;j<10;j++){
     fHist_znCx_V0_VxVy[i][j] = NULL;
     fHist_znCy_V0_VxVy[i][j] = NULL;
     fHist_znAx_V0_VxVy[i][j] = NULL;
     fHist_znAy_V0_VxVy[i][j] = NULL;
    }     
  }

  for(int i=0;i<4;i++){
    fHist_Qx_wiCorr_RunByRun[i] = NULL;
    fHist_Qx_Trig_woCorr[i] = NULL;
    fHist_XX_Trig_woCorr[i] = NULL;
    fHist_Qx_Trig_wiCorr[i] = NULL;
    fHist_XX_Trig_wiCorr[i] = NULL;

    for(int j=0;j<5;j++){
     fHist_Qx_vs_Obs_woCorr[i][j] = NULL;
     fHist_XX_vs_Obs_woCorr[i][j] = NULL;
     fHist_Qx_vs_Obs_wiCorr[i][j] = NULL;
     fHist_XX_vs_Obs_wiCorr[i][j] = NULL;
    }
  }
  for(int i=0;i<20;i++){
    fHist_Recenter_ZDCCx[i] = NULL;
    fHist_Recenter_ZDCCy[i] = NULL;
    fHist_Recenter_ZDCAx[i] = NULL;
    fHist_Recenter_ZDCAy[i] = NULL;
  }
  for(int i=0;i<20;i++){
    fHist_ZDCC_En_CommonCh[i] = NULL;
    fHist_ZDCA_En_CommonCh[i] = NULL;
  }

  for(int i=0;i<2;i++){
    VxCut[i] = 0;
    VyCut[i] = 0;
    VzCut[i] = 0;
  }

  for(int i=0;i<10;i++){
    fFB_Efficiency_Cent[i] = NULL;
    fHist_v2xV1_ZDN_pTDiff_All[i] = NULL;
    fHist_v4xV1_ZDN_pTDiff_All[i] = NULL;
    fHist_v3xV1_ZDN_EtaDiff_Comb1[i]= NULL;
    fHist_v3xV1_ZDN_EtaDiff_Comb2[i]= NULL;

    for(int j=0;j<4;j++){
      fHist_v1xV1_ZDN_EtaDiff[j][i] = NULL;
      fHist_v1xV1_ZDN_pTDiff[j][i]  = NULL;
    }
  }

  for(int i=0;i<4;i++){
    fHist_v2xV1_ZDN_Norm_Sep[i] = NULL;
    fHist_v2xV1_ZDN_Cent_Sep[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_ZDN_resol_Norm_Sep[i] = NULL;
    fHist_ZDN_resol_Cent_Sep[i] = NULL;
  }
  for(int i=0;i<4;i++){
    fHist_ZDCC_AvgCosNPsi[i]  = NULL;
    fHist_ZDCC_AvgSinNPsi[i]  = NULL;
    fHist_ZDCA_AvgCosNPsi[i]  = NULL;
    fHist_ZDCA_AvgSinNPsi[i]  = NULL;
    fHist_ZDC_AvgCosNPsiAC[i] = NULL;
    fHist_ZDC_AvgSinNPsiAC[i] = NULL;
  }
  for(int i=0;i<4;i++){
    fHist_ZDCC_AvgCos_VsRun[i]  = NULL;
    fHist_ZDCC_AvgSin_VsRun[i]  = NULL;
    fHist_ZDCA_AvgCos_VsRun[i]  = NULL;
    fHist_ZDCA_AvgSin_VsRun[i]  = NULL;
    fHist_ZDC_AvgCosPsiSum_VsRun[i] = NULL;
    fHist_ZDC_AvgSinPsiSum_VsRun[i] = NULL;
    fHist_ZDC_AvgCosPsiDif_VsRun[i] = NULL;
    fHist_ZDC_AvgSinPsiDif_VsRun[i] = NULL;
  }
  for(int i=0;i<2;i++){
    fHist_XXYY_vs_Cent_woCorr[i] = NULL;
    fHist_XXYY_vs_Cent_wiCorr[i] = NULL;
  }
  for(int i=0;i<4;i++){
   for(int j=0;j<6;j++){
    fHist_XX_vs_QnC_2DwoCorr_PosMag[i][j] = NULL;
    fHist_XX_vs_QnA_2DwoCorr_PosMag[i][j] = NULL;
    fHist_XX_vs_QnC_2DwoCorr_NegMag[i][j] = NULL;
    fHist_XX_vs_QnA_2DwoCorr_NegMag[i][j] = NULL;
   }
  }

  for(int i=0;i<90;i++){
   fHist_VZERO_Mult_vsRun[i] = NULL;
  }


  //fDataSet="2010";
  //fAnalysisSet="DoGainEq";
  //sCentEstimator="V0";
}





//________________________________________________________________________
AliAnalysisTaskZDCGainEq::~AliAnalysisTaskZDCGainEq()
{
  delete                 fListHistos;        
  delete                 fListZDCQxy;         
  delete                 fListZDCWgt;         
  delete                 fListDummy1;        
  delete                 fListHijing;        

  delete              fMultSelection; 
  delete               fAnalysisUtil; // it is '= new' !!!


 //these histograms are not in any list:
 //therefore, deleted manually.
 //delete         *fHist_ChanWgt_ZDCC;  //can't delete, Execption thrown.!! why?
 //delete         *fHist_ChanWgt_ZDCA;  
 //delete       *fHist_Vx_ArrayFinder; 
 //delete       *fHist_Vy_ArrayFinder; 
 //delete       *fHist_Vz_ArrayFinder; 


 //printf("\n\n ~AliAnalysisTaskZDCGainEq::Destructor is called..!!\n\n");
}

//________________________________________________________________________
void AliAnalysisTaskZDCGainEq::UserCreateOutputObjects()
{

 Int_t runArray_2010[89] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137243, 137236, 137235, 137232, 137231, 137162, 137161};

 Int_t runArray_2011[68] = {167915, 168115, 168460, 169035, 169238, 169859, 170228, 167920, 168310, 168464, 169091, 169411, 169923, 170230, 167985, 168311, 168467, 169094, 169415, 170027, 170268, 167987, 168322, 168511, 169138, 169417, 170081, 170269, 167988, 168325, 168512, 169144, 169835, 170155, 170270, 168069, 168341, 168514, 169145, 169837, 170159, 170306, 168076, 168342, 168777, 169148, 169838, 170163, 170308, 168105, 168361, 168826, 169156, 169846, 170193, 170309, 168107, 168362, 168988, 169160, 169855, 170203, 168108, 168458, 168992, 169167, 169858, 170204};

 Int_t runArray_2015[90] = {246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246676, 246675, 246540, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};


 if(fDataSet=="2010"){
  frunflag = 89;
  for(int i=0;i<frunflag;i++)
    runNums[i] = runArray_2010[i];
 }
 if(fDataSet=="2011"){
   frunflag = 10;                 //<--------- 2011 runnumbers have to be checked and put...
  for(int i=0;i<frunflag;i++)
    runNums[i] = runArray_2011[i];
 }
 if(fDataSet=="2015"){
  frunflag = 90;
  for(int i=0;i<frunflag;i++)
    runNums[i] = runArray_2015[i];
 }

  fListHistos = new TList();
  fListHistos->SetOwner(kTRUE);

  fHist_Event_count    = new TH1F("fHist_Event_count"," ",20,0,20);
  fHist_Event_count->GetXaxis()->SetBinLabel(1,"Called Exec()");
  fHist_Event_count->GetXaxis()->SetBinLabel(2,"AOD Exist");
  fHist_Event_count->GetXaxis()->SetBinLabel(3,"Pass VzCut");
  fHist_Event_count->GetXaxis()->SetBinLabel(4,"Pass VxCut");
  fHist_Event_count->GetXaxis()->SetBinLabel(5,"Pass VyCut");
  fHist_Event_count->GetXaxis()->SetBinLabel(6,"NotPileUp");
  fHist_Event_count->GetXaxis()->SetBinLabel(7,"ZNC Ei>=0");
  fHist_Event_count->GetXaxis()->SetBinLabel(8,"ZNA Ei>=0");
  fHist_Event_count->GetXaxis()->SetBinLabel(9,"|QnC| < 1.5");
  fHist_Event_count->GetXaxis()->SetBinLabel(10,"|QnA| < 1.5");
  fHist_Event_count->GetXaxis()->SetBinLabel(11,"#Psi_{A}=0 && #Psi_{C}=0");
  fHist_Event_count->GetXaxis()->SetBinLabel(12,"QyA=0 && QyC=0");
  fListHistos->Add(fHist_Event_count);

  fPileUpCount = new TH1F("fPileUpCount", "fPileUpCount", 9, 0., 9.);
  fPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultiplicityComb08");
  fPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>7.5");
  fPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(8,"multESDTPCDif");
  fPileUpCount->GetXaxis()->SetBinLabel(9,"extraPileUpMultSel");
  fListHistos->Add(fPileUpCount);

  fPileUpMultSelCount = new TH1F("fPileUpMultSelCount", "fPileUpMultSelCount", 8, 0., 8.);
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(1,"IsNotPileup");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(2,"IsNotPileupMV");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(3,"IsNotPileupInMultBins");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(4,"InconsistentVertices");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(5,"TrackletVsCluster");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(6,"AsymmetricInVZERO");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(7,"IncompleteDAQ");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(8,"GoodVertex2016");
  fListHistos->Add(fPileUpMultSelCount);

  fHist_Task_config = new TH1F("fTask_Configuration", "Task Configuration", 20, 0., 20.);
  fHist_Task_config->GetXaxis()->SetBinLabel(1,"IsZDCQAon");
  fHist_Task_config->GetXaxis()->SetBinLabel(2,"IsCentV0M");
  fHist_Task_config->GetXaxis()->SetBinLabel(3,"IsCentTPC");
  fHist_Task_config->GetXaxis()->SetBinLabel(4,"IsCentCL1");
  fHist_Task_config->GetXaxis()->SetBinLabel(5,"IsPileUpOn");
  fHist_Task_config->GetXaxis()->SetBinLabel(6,"IsPileUpTightOn");
  fHist_Task_config->GetXaxis()->SetBinLabel(7,"IsFillGainEq");
  fHist_Task_config->GetXaxis()->SetBinLabel(8,"IsDoGainEq");
  fHist_Task_config->GetXaxis()->SetBinLabel(9,"IsAvailGainFile");
  fHist_Task_config->GetXaxis()->SetBinLabel(10,"IsFillCosSin");
  fHist_Task_config->GetXaxis()->SetBinLabel(11,"IsDoRecenter");
  fHist_Task_config->GetXaxis()->SetBinLabel(12,"IsAvailRecntFile");
  fHist_Task_config->GetXaxis()->SetBinLabel(13,"IsRunByRun");
  fHist_Task_config->GetXaxis()->SetBinLabel(14,"IsCentCutforShift");
  fHist_Task_config->GetXaxis()->SetBinLabel(15,"IsApplyShiftCorr");
  fHist_Task_config->GetXaxis()->SetBinLabel(16,"IsShiftCorrVsCent");
  fHist_Task_config->GetXaxis()->SetBinLabel(17,"IsTrigonoMetricQn");
  fListHistos->Add(fHist_Task_config);


  fHist_CutParameters = new TH1F("fTask_CutParameters", "Variable Cut Values", 10, 0., 10.);
  fHist_CutParameters->GetXaxis()->SetBinLabel(1,"VzCutLowValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(2,"VzCutHighValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(3,"NBins_Vz");
  fHist_CutParameters->GetXaxis()->SetBinLabel(4,"VxCutHighValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(5,"VxCutLowValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(6,"NBins_Vx");
  fHist_CutParameters->GetXaxis()->SetBinLabel(7,"VyCutHighValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(8,"VyCutLowValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(9,"NBins_Vy");
  fListHistos->Add(fHist_CutParameters);

  fHist_ZeroQnXY_Conter  =  new TH1F("fHist_ZeroQnXY_Conter","", 20, 0, 20);
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(1,"QxC = 0 (w/o)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(2,"QxC = 0 (wi)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(3,"QyC = 0 (w/o)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(4,"QyC = 0 (wi)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(5,"QxA = 0 (w/o)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(6,"QxA = 0 (wi)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(7,"QyA = 0 (w/o)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(8,"QyA = 0 (wi)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(9, "QxC = 0 && QyC = 0 (w/o)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(10,"QxC = 0 && QyC = 0 (wi)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(11,"QxA = 0 && QyA = 0 (w/o)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(12,"QxA = 0 && QyA = 0 (wi)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(13,"QxA,QyA,QxC,QyC= 0 (w/o)");
  fHist_ZeroQnXY_Conter->GetXaxis()->SetBinLabel(14,"QxA,QyA,QxC,QyC= 0 (wi)");
  fListHistos->Add(fHist_ZeroQnXY_Conter);



  fHist_Cent_woZDCcut =  new TH1F("fHist_Cent_before_ZDCcut"," ",100,0,100);
  fListHistos->Add(fHist_Cent_woZDCcut);
  fHist_Cent_wiZDCcut =  new TH1F("fHist_Cent_afterr_ZDCcut"," ",100,0,100);
  fListHistos->Add(fHist_Cent_wiZDCcut);





 //-------- define and add the recenter histograms ----------------
  if(fDataSet=="2010") {
    //vxBin = 9; 
    //vyBin = 9;
    //vzBin = 10;   

    VxCut[0] =  -0.030;
    VxCut[1] =   0.015;
    VyCut[0] =   0.150;
    VyCut[1] =   0.204;
    VzCut[0] =   -10.0;
    VzCut[1] =     9.4;
  }

  if(fDataSet=="2015"){
    //vxBin = 8;
    //vyBin = 8;
    //vzBin = 10;  

    VxCut[0] =   0.060;
    VxCut[1] =   0.086;
    VyCut[0] =   0.321;
    VyCut[1] =   0.345;
    VzCut[0] =   -10.0;
    VzCut[1] =    10.0;
  }

  if(fDataSet=="2011"){
     AliDebug(2,"\n\n !!** WARNING ***!!  \n cuts not defined for DATASET: %s\n ...EXIT...\n\n)");
     exit(1);
  }

  const int NbinVt =  vyBin*vxBin;

  fHist_Vx_ArrayFinder = new TH1F("fHist_Vx_ArrayFinder","",vxBin,VxCut[0],VxCut[1]);
  fHist_Vy_ArrayFinder = new TH1F("fHist_Vy_ArrayFinder","",vyBin,VyCut[0],VyCut[1]);
  fHist_Vz_ArrayFinder = new TH1F("fHist_Vz_ArrayFinder","",vzBin,VzCut[0],VzCut[1]);

  fHist_Vxy_RunAveraged   = new TH2F("fHist_Vxy_RunAveraged","Vy Vs Vx (RbR)",100,0.050,0.10,100,0.31,0.36);
  fListHistos->Add(fHist_Vxy_RunAveraged);


  fHist_ZDCC_AvgQx_VsRun = new TProfile2D(Form("fHist_ZDCC_AvgQx_VsCentRun"),"",60,0,60,90,0,90,""); //cent,run
  fListHistos->Add(fHist_ZDCC_AvgQx_VsRun);
  fHist_ZDCC_AvgQy_VsRun = new TProfile2D(Form("fHist_ZDCC_AvgQy_VsCentRun"),"",60,0,60,90,0,90,""); //cent,run
  fListHistos->Add(fHist_ZDCC_AvgQy_VsRun);

  fHist_ZDCA_AvgQx_VsRun = new TProfile2D(Form("fHist_ZDCA_AvgQx_VsCentRun"),"",60,0,60,90,0,90,""); //cent,run
  fListHistos->Add(fHist_ZDCA_AvgQx_VsRun);
  fHist_ZDCA_AvgQy_VsRun = new TProfile2D(Form("fHist_ZDCA_AvgQy_VsCentRun"),"",60,0,60,90,0,90,""); //cent,run
  fListHistos->Add(fHist_ZDCA_AvgQy_VsRun);

  fHist_ZDC_AvgXXminusYY_VsRun = new TProfile2D(Form("fHist_ZDCC_AvgXXminusYY_VsCentRun"),"",60,0,60,90,0,90,""); //cent,run
  fListHistos->Add(fHist_ZDC_AvgXXminusYY_VsRun);
  fHist_ZDC_AvgXYplusXY_VsRun = new TProfile2D(Form("fHist_ZDCC_AvgXYplusXY_VsCentRun"),"",60,0,60,90,0,90,""); //cent,run
  fListHistos->Add(fHist_ZDC_AvgXYplusXY_VsRun);


  fHist_ZDC_dTermXXYY_VsRun = new TProfile2D(Form("fHist_ZDC_dTermXXYY_VsCentRun"),"",60,0,60,90,0,90,""); //cent,run
  fListHistos->Add(fHist_ZDC_dTermXXYY_VsRun);
  fHist_ZDC_dTermXYXY_VsRun = new TProfile2D(Form("fHist_ZDC_dTermXYXY_VsCentRun"),"",60,0,60,90,0,90,""); //cent,run
  fListHistos->Add(fHist_ZDC_dTermXYXY_VsRun);


  //Shift histograms: Fill Cos(nPsi) Sin(nPsi) vs Cent,Run
  for(int i=0; i<4; i++){
    fHist_ZDCC_AvgCos_VsRun[i] = new TProfile2D(Form("fHist_ZDC_AvgCos%dPsiC_VsCentRun",i+1),"",60,0,60,90,0,90,""); //cent,run
    fListHistos->Add(fHist_ZDCC_AvgCos_VsRun[i]);
    fHist_ZDCA_AvgCos_VsRun[i] = new TProfile2D(Form("fHist_ZDC_AvgCos%dPsiA_VsCentRun",i+1),"",60,0,60,90,0,90,"");
    fListHistos->Add(fHist_ZDCA_AvgCos_VsRun[i]);

    fHist_ZDCC_AvgSin_VsRun[i] = new TProfile2D(Form("fHist_ZDC_AvgSin%dPsiC_VsCentRun",i+1),"",60,0,60,90,0,90,"");
    fListHistos->Add(fHist_ZDCC_AvgSin_VsRun[i]);
    fHist_ZDCA_AvgSin_VsRun[i] = new TProfile2D(Form("fHist_ZDC_AvgSin%dPsiA_VsCentRun",i+1),"",60,0,60,90,0,90,"");
    fListHistos->Add(fHist_ZDCA_AvgSin_VsRun[i]);

    fHist_ZDC_AvgCosPsiSum_VsRun[i] = new TProfile2D(Form("fHist_ZDC_AvgCos%dPsi_AplusC_VsCentRun",i+1),"",60,0,60,90,0,90,"");
    fListHistos->Add(fHist_ZDC_AvgCosPsiSum_VsRun[i]);
    fHist_ZDC_AvgSinPsiSum_VsRun[i] = new TProfile2D(Form("fHist_ZDC_AvgSin%dPsi_AplusC_VsCentRun",i+1),"",60,0,60,90,0,90,"");
    fListHistos->Add(fHist_ZDC_AvgSinPsiSum_VsRun[i]);

    fHist_ZDC_AvgCosPsiDif_VsRun[i] = new TProfile2D(Form("fHist_ZDC_AvgCos%dPsi_AminusC_VsCentRun",i+1),"",60,0,60,90,0,90,"");
    fListHistos->Add(fHist_ZDC_AvgCosPsiDif_VsRun[i]);
    fHist_ZDC_AvgSinPsiDif_VsRun[i] = new TProfile2D(Form("fHist_ZDC_AvgSin%dPsi_AminusC_VsCentRun",i+1),"",60,0,60,90,0,90,"");
    fListHistos->Add(fHist_ZDC_AvgSinPsiDif_VsRun[i]);
  }




  fHist_XXYY_vs_Cent_woCorr[0] = new TProfile(Form("fHist_XXminusYY_vs_Cent_woCorr"),"XX-YY",90,0,90,"");
  fListHistos->Add(fHist_XXYY_vs_Cent_woCorr[0]);
  fHist_XXYY_vs_Cent_woCorr[1] = new TProfile(Form("fHist_XYplusXY_vs_Cent_woCorr"), "XY+XY",90,0,90,"");
  fListHistos->Add(fHist_XXYY_vs_Cent_woCorr[1]);
  fHist_XXYY_vs_Cent_wiCorr[0] = new TProfile(Form("fHist_XXminusYY_vs_Cent_wiCorr"),"XX-YY",90,0,90,"");
  fListHistos->Add(fHist_XXYY_vs_Cent_wiCorr[0]);
  fHist_XXYY_vs_Cent_wiCorr[1] = new TProfile(Form("fHist_XYplusXY_vs_Cent_wiCorr"), "XY+XY",90,0,90,"");
  fListHistos->Add(fHist_XXYY_vs_Cent_wiCorr[1]);

  fHist_Event_counter_vRun    =  new TH1F("fHist_Event_counter_vRun","",95,0,95);
  fListHistos->Add(fHist_Event_counter_vRun);


  char name[100];


  if(fAnalysisSet=="FillGainEq") {

    fHist_Vx_vs_runnum = new TProfile("fHist_Vx_vs_runnum","<Vx>_vs_runnum",frunflag,0,frunflag,"s");
    fListHistos->Add(fHist_Vx_vs_runnum);
    fHist_Vy_vs_runnum = new TProfile("fHist_Vy_vs_runnum","<Vy>_vs_runnum",frunflag,0,frunflag,"s");
    fListHistos->Add(fHist_Vy_vs_runnum);
    fHist_Vz_vs_runnum = new TProfile("fHist_Vz_vs_runnum","<Vz>_vs_runnum",frunflag,0,frunflag,"s");
    fListHistos->Add(fHist_Vz_vs_runnum);

    for(int i=0;i<vzBin;i++) {
      fHist_ZDCC_En_CommonCh[i] = new TProfile2D(Form("fHist_ZDCC_En_CommonCh_Vz%d",i+1),"",100,0,100,NbinVt,0,NbinVt,"");
      fHist_ZDCA_En_CommonCh[i] = new TProfile2D(Form("fHist_ZDCA_En_CommonCh_Vz%d",i+1),"",100,0,100,NbinVt,0,NbinVt,"");
      fListHistos->Add(fHist_ZDCC_En_CommonCh[i]);
      fListHistos->Add(fHist_ZDCA_En_CommonCh[i]);
    }
    for(int i=0;i<frunflag;i++) {
      //store: ZDC energy for gain calibration:
      fHist_ZDCA_En_Run[i]  = new TProfile2D(Form("fHist_ZDCA_En_Run%d",runNums[i]),"",100,0,100,5,0,5,"");
      fListHistos->Add(fHist_ZDCA_En_Run[i]);
      fHist_ZDCC_En_Run[i]  = new TProfile2D(Form("fHist_ZDCC_En_Run%d",runNums[i]),"",100,0,100,5,0,5,"");
      fListHistos->Add(fHist_ZDCC_En_Run[i]);
    }
  }


  if(bFillZDCQAon){

   fHist_Task_config->Fill(0.5);

   if(fAnalysisSet=="FillGainEq"){
     fHist_Psi1_ZDCC_wGainCorr      = new TH1F("fHist_Psi1_ZDCC_woGainCorr","", 200, 0,2.*TMath::Pi());
     fHist_Psi1_ZDCA_wGainCorr      = new TH1F("fHist_Psi1_ZDCA_woGainCorr","", 200, 0,2.*TMath::Pi());
   }
   else{
     fHist_Psi1_ZDCC_wGainCorr      = new TH1F("fHist_Psi1_ZDCC_wiGainCorr","", 200, 0,2.*TMath::Pi());
     fHist_Psi1_ZDCA_wGainCorr      = new TH1F("fHist_Psi1_ZDCA_wiGainCorr","", 200, 0,2.*TMath::Pi());
   }

   fListHistos->Add(fHist_Psi1_ZDCC_wGainCorr);
   fListHistos->Add(fHist_Psi1_ZDCA_wGainCorr);

   fHist_Psi1_ZDCC_wRectCorr      = new TH1F("fHist_Psi1_ZDCC_wiRectCorr","", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_Psi1_ZDCC_wRectCorr);
   fHist_Psi1_ZDCA_wRectCorr      = new TH1F("fHist_Psi1_ZDCA_wiRectCorr","", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_Psi1_ZDCA_wRectCorr);

   fHist_Psi1_ZDCC_RunByRun      = new TH2F("fHist_Psi1_ZDCC_RunByRun","", 100, 0, 100, 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_Psi1_ZDCC_RunByRun);
   fHist_Psi1_ZDCA_RunByRun      = new TH2F("fHist_Psi1_ZDCA_RunByRun","", 100, 0, 100, 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_Psi1_ZDCA_RunByRun);

   fHist_Psi1_ZDCC_wCorrFull      = new TH1F("fHist_Psi1_ZDCC_wCorrFull","", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_Psi1_ZDCC_wCorrFull);
   fHist_Psi1_ZDCA_wCorrFull      = new TH1F("fHist_Psi1_ZDCA_wCorrFull","", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_Psi1_ZDCA_wCorrFull);


   fHist_PsiSumAC_woCorr    =  new TH1F("fHist_PsiSumAC_woCorr","", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_PsiSumAC_woCorr);
   fHist_PsiSumAC_wiCorr    =  new TH1F("fHist_PsiSumAC_wiCorr","", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_PsiSumAC_wiCorr);


   fHist_PsiSumAC_ZeroQn  =  new TH1F("fHist_PsiSumAC_ZeroQn","", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_PsiSumAC_ZeroQn);
   fHist_PsiZDCA_ZeroQn   =  new TH1F("fHist_PsiZDCA_ZeroQn"," ", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_PsiZDCA_ZeroQn);
   fHist_PsiZDCC_ZeroQn   =  new TH1F("fHist_PsiZDCC_ZeroQn"," ", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_PsiZDCC_ZeroQn);


   TString     sNameQn[4] = {"Xa","Xc","Ya","Yc"};
   TString    sNameQn2[4] = {"XaXc","YaYc","XcYa","YcXa"};
   TString     sNameVs[5] = {"Cent","Mult","Vx","Vy","Vz"};      
   Int_t     nBinNumVs[5] = {100,  100,       50,       50, 400}; // number of bins for "Cent", "Mult", "Vx","Vy","Vz"
   Double_t  lBinLowVs[5] = {  0,    0, VxCut[0], VyCut[0], -10}; // bin low  value: "Cent", "Mult", "Vx", "Vy", "Vz"
   Double_t lBinHighVs[5] = {100, 4000, VxCut[1], VyCut[1],  10}; // bin high value: "Cent", "Mult", "Vx", "Vy", "Vz"
   Double_t fNormCent[11] = {0,5,10,20,30,40,50,60,70,80,90};

   fHist_NormalCentralityBins = new TH1F("fHist_NormalCentralityBins","",10,fNormCent);
   fListHistos->Add(fHist_NormalCentralityBins);

   //here:
   for(int i=0;i<6;i++) {
     for(int j=0;j<4;j++) {
       sprintf(name,"fHist_%s_vs_Cx_Cy_woCorr_PosMag_Cent%d",static_cast<const char*>(sNameQn2[j]),i);
       fHist_XX_vs_QnC_2DwoCorr_PosMag[j][i] = new TProfile2D(name,"",100, -1.0, 1.0, 100, -1.0, 1.0, "");
       fListHistos->Add(fHist_XX_vs_QnC_2DwoCorr_PosMag[j][i]);

       sprintf(name,"fHist_%s_vs_Ax_Ay_woCorr_PosMag_Cent%d",static_cast<const char*>(sNameQn2[j]),i);
       fHist_XX_vs_QnA_2DwoCorr_PosMag[j][i] = new TProfile2D(name,"",100, -1.0, 1.0, 100, -1.0, 1.0, "");
       fListHistos->Add(fHist_XX_vs_QnA_2DwoCorr_PosMag[j][i]);

       sprintf(name,"fHist_%s_vs_Cx_Cy_woCorr_NegMag_Cent%d",static_cast<const char*>(sNameQn2[j]),i);
       fHist_XX_vs_QnC_2DwoCorr_NegMag[j][i] = new TProfile2D(name,"",100, -1.0, 1.0, 100, -1.0, 1.0, "");
       fListHistos->Add(fHist_XX_vs_QnC_2DwoCorr_NegMag[j][i]);

       sprintf(name,"fHist_%s_vs_Ax_Ay_woCorr_NegMag_Cent%d",static_cast<const char*>(sNameQn2[j]),i);
       fHist_XX_vs_QnA_2DwoCorr_NegMag[j][i] = new TProfile2D(name,"",100, -1.0, 1.0, 100, -1.0, 1.0, "");
       fListHistos->Add(fHist_XX_vs_QnA_2DwoCorr_NegMag[j][i]);
     }
   }



   for(int i=0;i<4;i++) {
     sprintf(name,"fHist_%s_wiCorr_RunByRun",static_cast<const char*>(sNameQn[i]));
     fHist_Qx_wiCorr_RunByRun[i] = new TProfile(name,"", 90, 0, 90,"");
     fListHistos->Add(fHist_Qx_wiCorr_RunByRun[i]);

     sprintf(name,"fHist_%s_TrigonoMetric_woCorr",static_cast<const char*>(sNameQn[i]));
     fHist_Qx_Trig_woCorr[i] = new TProfile(name,"", 90, 0, 90,""); //centrality range
     fListHistos->Add(fHist_Qx_Trig_woCorr[i]);
     sprintf(name,"fHist_%s_TrigonoMetric_woCorr",static_cast<const char*>(sNameQn2[i]));
     fHist_XX_Trig_woCorr[i] = new TProfile(name,"", 90, 0, 90,""); //centrality range
     fListHistos->Add(fHist_XX_Trig_woCorr[i]);

     sprintf(name,"fHist_%s_TrigonoMetric_wiCorr",static_cast<const char*>(sNameQn[i]));
     fHist_Qx_Trig_wiCorr[i] = new TProfile(name,"", 90, 0, 90,""); //centrality range
     fListHistos->Add(fHist_Qx_Trig_wiCorr[i]);
     sprintf(name,"fHist_%s_TrigonoMetric_wiCorr",static_cast<const char*>(sNameQn2[i]));
     fHist_XX_Trig_wiCorr[i] = new TProfile(name,"", 90, 0, 90,""); //centrality range
     fListHistos->Add(fHist_XX_Trig_wiCorr[i]);

     for(int j=0;j<5;j++) {//fHist_Qx_vs_Obs_woCorr
      //store: X,Y position for recenter:
       sprintf(name,"fHist_%s_vs_%s_woCorr",static_cast<const char*>(sNameQn[i]),static_cast<const char*>(sNameVs[j]));
       fHist_Qx_vs_Obs_woCorr[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
       fListHistos->Add(fHist_Qx_vs_Obs_woCorr[i][j]);

       sprintf(name,"fHist_%s_vs_%s_woCorr",static_cast<const char*>(sNameQn2[i]),static_cast<const char*>(sNameVs[j]));
       fHist_XX_vs_Obs_woCorr[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
       fListHistos->Add(fHist_XX_vs_Obs_woCorr[i][j]);

       if(bFillCosSin && fAnalysisSet=="FillGainEq") {
         sprintf(name,"fHist_%s_vs_%s_woCorrCosSin",static_cast<const char*>(sNameQn[i]),static_cast<const char*>(sNameVs[j]));
       }
       else{
         sprintf(name,"fHist_%s_vs_%s_wiCorr",static_cast<const char*>(sNameQn[i]),static_cast<const char*>(sNameVs[j]));
       }
       fHist_Qx_vs_Obs_wiCorr[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
       fListHistos->Add(fHist_Qx_vs_Obs_wiCorr[i][j]);

       if(bFillCosSin && fAnalysisSet=="FillGainEq") {
         sprintf(name,"fHist_%s_vs_%s_woCorrCosSin",static_cast<const char*>(sNameQn2[i]),static_cast<const char*>(sNameVs[j]));
       }
       else{
         sprintf(name,"fHist_%s_vs_%s_wiCorr",static_cast<const char*>(sNameQn2[i]),static_cast<const char*>(sNameVs[j]));
       }

       fHist_XX_vs_Obs_wiCorr[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
       fListHistos->Add(fHist_XX_vs_Obs_wiCorr[i][j]);
     }
   }
  }


  if(vzBin>20){
    printf("\n\n::UserCreateOutPutObject()\n Vz Binning more than 20 not allowed\n ==> Exit <== \n\n");
    exit(1);
  }

  const int VzBinIter = (const int) vzBin;


  //filling the Q vectors for recenter:
  if(fAnalysisSet=="DoGainEq") {
   if(!bApplyRecent) {

    fListZDCQxy = new TList();
    fListZDCQxy->SetOwner(kTRUE);

    if(!bRunAveragedQn) {
      for(int i=0;i<frunflag;i++) {
       for(int j=0;j<VzBinIter;j++) {
          fHist_znCx_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znCx_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90);
          fHist_znCy_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znCy_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90);
          fHist_znAx_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znAx_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90);
          fHist_znAy_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znAy_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90);

          fListZDCQxy->Add(fHist_znCx_V0_VxVy[i][j]);
          fListZDCQxy->Add(fHist_znCy_V0_VxVy[i][j]);
          fListZDCQxy->Add(fHist_znAx_V0_VxVy[i][j]);
          fListZDCQxy->Add(fHist_znAy_V0_VxVy[i][j]);
       }
      }
    }
    else{
      for(int j=0;j<VzBinIter;j++) {
          fHist_znCx_V0_VxVy[0][j] = new TProfile2D(Form("fHist_znCx_V0_Run%d_Vz%d",0,j+1),"",NbinVt,0,NbinVt,90,0,90); 
          fHist_znCy_V0_VxVy[0][j] = new TProfile2D(Form("fHist_znCy_V0_Run%d_Vz%d",0,j+1),"",NbinVt,0,NbinVt,90,0,90); 
          fHist_znAx_V0_VxVy[0][j] = new TProfile2D(Form("fHist_znAx_V0_Run%d_Vz%d",0,j+1),"",NbinVt,0,NbinVt,90,0,90);
          fHist_znAy_V0_VxVy[0][j] = new TProfile2D(Form("fHist_znAy_V0_Run%d_Vz%d",0,j+1),"",NbinVt,0,NbinVt,90,0,90);

          fListZDCQxy->Add(fHist_znCx_V0_VxVy[0][j]);
          fListZDCQxy->Add(fHist_znCy_V0_VxVy[0][j]);
          fListZDCQxy->Add(fHist_znAx_V0_VxVy[0][j]);
          fListZDCQxy->Add(fHist_znAy_V0_VxVy[0][j]);
      }
    }
   }
  }

  for(int i=0;i<90;i++){
    fHist_VZERO_Mult_vsRun[i]  = new TProfile2D(Form("fHist_V0Mult_vsCent_Run%d",runArray_2015[i]),"",90,0,90,64,0,64,""); //cent,channel
    //added to second list
  }


  Double_t centRange[12]    = {0,5,10,20,30,40,50,60,70,80,90,100};
  //Double_t pTRange[21] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.3,3.6,4.0,4.5,5.0};
  Double_t pTRange[28] = {0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.33,2.66,3.,3.5,4.,5.,6.,8.,10.,14.,20.,30.,50.};
  const int fPtDiffNBins = 27;
  TString  sNameComp[4] = {"uxCx","uyCy","uxAx","uyAy"};


  fHist_Corr3p_ZDN_Norm_PN = new TProfile("fHist_Corr3p_ZDN_Norm_PosNeg","opposit charge correlator",11,centRange,"");
  fHist_Corr3p_ZDN_Norm_PN->Sumw2();
  fListHistos->Add(fHist_Corr3p_ZDN_Norm_PN);
  fHist_Corr3p_ZDN_Norm_PP = new TProfile("fHist_Corr3p_ZDN_Norm_PosPos","pos-pos charge correlator",11,centRange,"");
  fHist_Corr3p_ZDN_Norm_PP->Sumw2();
  fListHistos->Add(fHist_Corr3p_ZDN_Norm_PP);
  fHist_Corr3p_ZDN_Norm_NN = new TProfile("fHist_Corr3p_ZDN_Norm_NegNeg","neg-neg charge correlator",11,centRange,"");
  fHist_Corr3p_ZDN_Norm_NN->Sumw2();
  fListHistos->Add(fHist_Corr3p_ZDN_Norm_NN);

  fHist_Corr3p_ZDN_Cent_PN = new TProfile("fHist_Corr3p_ZDN_Cent_PosNeg","opposit charge correlator",90,0,90,"");
  fHist_Corr3p_ZDN_Cent_PN->Sumw2();
  fListHistos->Add(fHist_Corr3p_ZDN_Cent_PN);
  fHist_Corr3p_ZDN_Cent_PP = new TProfile("fHist_Corr3p_ZDN_Cent_PosPos","pos-pos charge correlator",90,0,90,"");
  fHist_Corr3p_ZDN_Cent_PP->Sumw2();
  fListHistos->Add(fHist_Corr3p_ZDN_Cent_PP);
  fHist_Corr3p_ZDN_Cent_NN = new TProfile("fHist_Corr3p_ZDN_Cent_NegNeg","neg-neg charge correlator",90,0,90,"");
  fHist_Corr3p_ZDN_Cent_NN->Sumw2();
  fListHistos->Add(fHist_Corr3p_ZDN_Cent_NN);

  fHist_Reso2EP_TPC_Norm  = new TProfile("fHist_Reso2EP_TPC_Norm","SubEvent Resolution",11,centRange,"");
  fHist_Reso2EP_TPC_Norm->Sumw2();
  fListHistos->Add(fHist_Reso2EP_TPC_Norm);
  fHist_Reso2EP_TPC_Cent  = new TProfile("fHist_Reso2EP_TPC_Cent","SubEvent Resolution",90,0,90,"");
  fHist_Reso2EP_TPC_Cent->Sumw2();
  fListHistos->Add(fHist_Reso2EP_TPC_Cent);


  if(bApplyRecent) {
  //v2 integrated:
   fHist_v2xV1_ZDN_Norm_All  = new TProfile("fHist_v2xV1_ZDN_Norm_All","v2 X V1^2 (ZDC) all terms",11,centRange,"");
   fHist_v2xV1_ZDN_Norm_All->Sumw2();
   fListHistos->Add(fHist_v2xV1_ZDN_Norm_All);
   fHist_v2xV1_ZDN_Cent_All  = new TProfile("fHist_v2xV1_ZDN_Cent_All","v2 X V1^2 (ZDC) all terms",90, 0,  90,"");
   fHist_v2xV1_ZDN_Cent_All->Sumw2();
   fListHistos->Add(fHist_v2xV1_ZDN_Cent_All);
   fHist_v2xV1_ZDN_Refm_All  = new TProfile("fHist_v2xV1_ZDN_Refm_All","v2 X V1^2 (ZDC) all terms", 40, 0, 4000,"");
   fHist_v2xV1_ZDN_Refm_All->Sumw2();
   fListHistos->Add(fHist_v2xV1_ZDN_Refm_All);

   //v3 integrated:
   fHist_v3xV1_ZDN_Norm_Comb1  = new TProfile("fHist_v3xV1_ZDN_Norm_Comb1","v3 X V1^2 (ZDC) all terms",11,centRange,"");
   fHist_v3xV1_ZDN_Norm_Comb1->Sumw2();
   fListHistos->Add(fHist_v3xV1_ZDN_Norm_Comb1);
   fHist_v3xV1_ZDN_Cent_Comb1  = new TProfile("fHist_v3xV1_ZDN_Cent_Comb1","v3 X V1^2 (ZDC) all terms",90, 0,  90,"");
   fHist_v3xV1_ZDN_Cent_Comb1->Sumw2();
   fListHistos->Add(fHist_v3xV1_ZDN_Cent_Comb1);
   fHist_v3xV1_ZDN_Norm_Comb2  = new TProfile("fHist_v3xV1_ZDN_Norm_Comb2","v3 X V1^2 (ZDC) all terms",11,centRange,"");
   fHist_v3xV1_ZDN_Norm_Comb2->Sumw2();
   fListHistos->Add(fHist_v3xV1_ZDN_Norm_Comb2);
   fHist_v3xV1_ZDN_Cent_Comb2  = new TProfile("fHist_v3xV1_ZDN_Cent_Comb2","v3 X V1^2 (ZDC) all terms",90, 0,  90,"");
   fHist_v3xV1_ZDN_Cent_Comb2->Sumw2();
   fListHistos->Add(fHist_v3xV1_ZDN_Cent_Comb2);

   //v4 integrated:
   fHist_v4xV1_ZDN_Norm_Comb1  = new TProfile("fHist_v4xV1_ZDN_Norm_Comb1","v4 X V1^2 (ZDC) all terms",11,centRange,"");
   fHist_v4xV1_ZDN_Norm_Comb1->Sumw2();
   fListHistos->Add(fHist_v4xV1_ZDN_Norm_Comb1);
   fHist_v4xV1_ZDN_Cent_Comb1  = new TProfile("fHist_v4xV1_ZDN_Cent_Comb1","v4 X V1^2 (ZDC) all terms",90, 0,  90,"");
   fHist_v4xV1_ZDN_Cent_Comb1->Sumw2();
   fListHistos->Add(fHist_v4xV1_ZDN_Cent_Comb1);

   //Summed ZDC resolution:
   fHist_ZDN_resol_Norm_All  = new TProfile("fHist_ZDN_resol_Norm_All","Resol. V1^2(ZDC) All",11,centRange,"");
   fHist_ZDN_resol_Norm_All->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Norm_All);
   fHist_ZDN_resol_Cent_All  = new TProfile("fHist_ZDN_resol_Cent_All","Resol. V1^2(ZDC) All",90, 0,  90,"");
   fHist_ZDN_resol_Cent_All->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Cent_All);
   fHist_ZDN_resol_Refm_All  = new TProfile("fHist_ZDN_resol_Refm_All","Resol. V1^2(ZDC) All",40, 0, 4000,"");
   fHist_ZDN_resol_Refm_All->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Refm_All);  

   //Individual ZDC resolutions vs cent:
   fHist_ZDN_resol_Norm_XX  = new TProfile("fHist_ZDN_resol_Norm_XX","Resol. V1^2(ZDC) XX",11,centRange,"");
   fHist_ZDN_resol_Norm_XX->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Norm_XX);
   fHist_ZDN_resol_Norm_YY  = new TProfile("fHist_ZDN_resol_Norm_YY","Resol. V1^2(ZDC) YY",11,centRange,"");
   fHist_ZDN_resol_Norm_YY->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Norm_YY);
   fHist_ZDN_resol_Cent_XX  = new TProfile("fHist_ZDN_resol_Cent_XX","Resol. V1^2(ZDC) XX",90, 0,  90,"");
   fHist_ZDN_resol_Cent_XX->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Cent_XX);
   fHist_ZDN_resol_Cent_YY  = new TProfile("fHist_ZDN_resol_Cent_YY","Resol. V1^2(ZDC) YY",90, 0,  90,"");
   fHist_ZDN_resol_Cent_YY->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Cent_YY);



   //v1 integrated, 4 terms, normal bin, fine Centr bin
   for(int i=0; i<4; i++) {
     sprintf(name,"fHist_v1_%s_ZDN_Norm",static_cast<const char*>(sNameComp[i]));
     fHist_v2xV1_ZDN_Norm_Sep[i] = new TProfile(name,"v1 vs Cent",11,centRange,"");
     fHist_v2xV1_ZDN_Norm_Sep[i]->Sumw2();
     fListHistos->Add(fHist_v2xV1_ZDN_Norm_Sep[i]);

     sprintf(name,"fHist_v1_%s_ZDN_Cent",static_cast<const char*>(sNameComp[i]));
     fHist_v2xV1_ZDN_Cent_Sep[i] = new TProfile(name,"v1 vs Cent",100,0,100,"");
     fHist_v2xV1_ZDN_Cent_Sep[i]->Sumw2();
     fListHistos->Add(fHist_v2xV1_ZDN_Cent_Sep[i]);
   }
   //v1 resulution, 2 terms, normal bin, fine Centr bin
   fHist_ZDN_resol_Norm_Sep[0] = new TProfile("fHist_ZDN_resol_XaXc_Norm"," XaXc vs Cent",11,centRange,"");
   fHist_ZDN_resol_Norm_Sep[0]->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Norm_Sep[0]);
   fHist_ZDN_resol_Norm_Sep[1] = new TProfile("fHist_ZDN_resol_YaYc_Norm"," YaYc vs Cent",11,centRange,"");
   fHist_ZDN_resol_Norm_Sep[1]->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Norm_Sep[1]);

   fHist_ZDN_resol_Cent_Sep[0] = new TProfile("fHist_ZDN_resol_XaXc_Cent"," XaXc vs Cent",100,0,100,"");
   fHist_ZDN_resol_Cent_Sep[0]->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Cent_Sep[0]);
   fHist_ZDN_resol_Cent_Sep[1] = new TProfile("fHist_ZDN_resol_YaYc_Cent"," YaYc vs Cent",100,0,100,"");
   fHist_ZDN_resol_Cent_Sep[1]->Sumw2();
   fListHistos->Add(fHist_ZDN_resol_Cent_Sep[1]);



   //v2 differential pT:
   for(int i=0; i<10; i++) {
     sprintf(name,"fHist_v2xV1_ZDN_pTDiff_Cent%d",i);
     fHist_v2xV1_ZDN_pTDiff_All[i]  = new TProfile(name,"v2 X V1^2 pTdiff",fPtDiffNBins,pTRange);
     fHist_v2xV1_ZDN_pTDiff_All[i]->Sumw2();
     fListHistos->Add(fHist_v2xV1_ZDN_pTDiff_All[i]);
   }

  //v3 differential Eta:
   for(int i=0; i<10; i++) {
     sprintf(name,"fHist_v3xV1_ZDN_EtaDiff_Comb1_Cent%d",i);
     fHist_v3xV1_ZDN_EtaDiff_Comb1[i]  = new TProfile(name,"v3 X V1^2 Etadiff",10,-0.8,0.8);
     fHist_v3xV1_ZDN_EtaDiff_Comb1[i]->Sumw2();
     fListHistos->Add(fHist_v3xV1_ZDN_EtaDiff_Comb1[i]);
     sprintf(name,"fHist_v3xV1_ZDN_EtaDiff_Comb2_Cent%d",i);
     fHist_v3xV1_ZDN_EtaDiff_Comb2[i]  = new TProfile(name,"v3 X V1^2 Etadiff",10,-0.8,0.8);
     fHist_v3xV1_ZDN_EtaDiff_Comb2[i]->Sumw2();
     fListHistos->Add(fHist_v3xV1_ZDN_EtaDiff_Comb2[i]);
   }

  //v4 differential pT:
   for(int i=0; i<10; i++) {
     sprintf(name,"fHist_v4xV1_ZDN_pTDiff_Cent%d",i);
     fHist_v4xV1_ZDN_pTDiff_All[i]  = new TProfile(name,"v4 X V1^2 pTdiff",fPtDiffNBins,pTRange);
     fHist_v4xV1_ZDN_pTDiff_All[i]->Sumw2();
     fListHistos->Add(fHist_v4xV1_ZDN_pTDiff_All[i]);
   }


   //v1 differential Eta, pT: not needed now.
   for(int i=0; i<10; i++) {
     for(int j=0; j<4; j++) {
       sprintf(name,"fHist_%s_ZDN_EtaDiff_Cent%d",static_cast<const char*>(sNameComp[j]),i);
       fHist_v1xV1_ZDN_EtaDiff[j][i] = new TProfile(name,"v1 Eta-diff", 5, -0.8, 0.8);
       fHist_v1xV1_ZDN_EtaDiff[j][i]->Sumw2();
       fListHistos->Add(fHist_v1xV1_ZDN_EtaDiff[j][i]);

       sprintf(name,"fHist_%s_ZDN_pTDiff_Cent%d",static_cast<const char*>(sNameComp[j]),i);
       fHist_v1xV1_ZDN_pTDiff[j][i] = new TProfile(name,"v1 Eta-diff", fPtDiffNBins, pTRange);
       fHist_v1xV1_ZDN_pTDiff[j][i]->Sumw2();
       fListHistos->Add(fHist_v1xV1_ZDN_pTDiff[j][i]);
     }
   } 

   //------------ calculate centrality weight: ------------
   TH1F *fCent_fromDATA;

   if(fListZDCQxy) {
     fCent_fromDATA = (TH1F *) fListZDCQxy->FindObject("fHist_Cent_afterr_ZDCcut");

     for(int i=0;i<4;i++){                                        
       fHist_ZDCC_AvgCosNPsi[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDC_AvgCos%dPsiC_VsCR",i+1));
       fHist_ZDCA_AvgCosNPsi[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDC_AvgCos%dPsiA_VsCR",i+1));
       fHist_ZDCC_AvgSinNPsi[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDC_AvgSin%dPsiC_VsCR",i+1));
       fHist_ZDCA_AvgSinNPsi[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDC_AvgSin%dPsiA_VsCR",i+1));

       fHist_ZDC_AvgCosNPsiAC[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDC_AvgCos%dPsiAplusC_VsCR",i+1));
       fHist_ZDC_AvgSinNPsiAC[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDC_AvgSin%dPsiAplusC_VsCR",i+1)); 
     }

     fHist_ZDCC_AvgQx_VsCR = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDCC_AvgQx_VsCR"));
     fHist_ZDCC_AvgQy_VsCR = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDCC_AvgQy_VsCR"));
     fHist_ZDCA_AvgQx_VsCR = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDCA_AvgQx_VsCR"));
     fHist_ZDCA_AvgQy_VsCR = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDCA_AvgQy_VsCR"));

     fHist_ZDC_dTermXXYY_VsCR = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDC_dTermXXYY_VsCR"));
     fHist_ZDC_dTermXYXY_VsCR = (TH2F *) fListZDCQxy->FindObject(Form("fHist_ZDC_dTermXYXY_VsCR"));
   }
   else{
     printf("\n\n ******** Running without Centrality weight, No Shift Correction !!******** \n\n");
     fCent_fromDATA = new TH1F("fCent_DATA","Centrality distribution",100,0,100);
     for(int i=1;i<=fCent_fromDATA->GetNbinsX();i++){
       fCent_fromDATA->SetBinContent(i,100);
     }
     for(int i=0;i<4;i++){                     
       fHist_ZDCC_AvgCosNPsi[i] = new TH2F(Form("fHist_ZDC_AvgCos%dPsiC_VsCR_dummy",i+1),"",60,0,60,90,0,90); //cent,run
       fHist_ZDCA_AvgCosNPsi[i] = new TH2F(Form("fHist_ZDC_AvgCos%dPsiA_VsCR_dummy",i+1),"",60,0,60,90,0,90);   
       fHist_ZDCC_AvgSinNPsi[i] = new TH2F(Form("fHist_ZDC_AvgSin%dPsiC_VsCR_dummy",i+1),"",60,0,60,90,0,90);   
       fHist_ZDCA_AvgSinNPsi[i] = new TH2F(Form("fHist_ZDC_AvgSin%dPsiA_VsCR_dummy",i+1),"",60,0,60,90,0,90);   
 
       fHist_ZDC_AvgCosNPsiAC[i] = new TH2F(Form("fHist_ZDC_AvgCos%dPsiAplusC_VsCR_dummy",i+1),"",60,0,60,90,0,90);   
       fHist_ZDC_AvgSinNPsiAC[i] = new TH2F(Form("fHist_ZDC_AvgSin%dPsiAplusC_VsCR_dummy",i+1),"",60,0,60,90,0,90);   
     }
     fHist_ZDCC_AvgQx_VsCR = new TH2F(Form("fHist_ZDCC_AvgQx_VsCR_dummy"),"",60,0,60,90,0,90); //cent,run
     fHist_ZDCC_AvgQy_VsCR = new TH2F(Form("fHist_ZDCC_AvgQy_VsCR_dummy"),"",60,0,60,90,0,90);   
     fHist_ZDCA_AvgQx_VsCR = new TH2F(Form("fHist_ZDCA_AvgQx_VsCR_dummy"),"",60,0,60,90,0,90);   
     fHist_ZDCA_AvgQy_VsCR = new TH2F(Form("fHist_ZDCA_AvgQy_VsCR_dummy"),"",60,0,60,90,0,90);   
   }

   Double_t Content,error;
   Double_t maxCount = -1;
   Int_t    maxCent  = -1;
   Double_t weight   = 0.;

   for(int i=1;i<=fCent_fromDATA->GetNbinsX();i++){
    Content = fCent_fromDATA->GetBinContent(i);
    if(maxCount < Content){
       maxCount = Content;
       maxCent = i;
    }
   }

   fWeight_Cent = new TH1F("fWeight_Cent","Weight for centrality",100,0,100);

   for(int i=1;i<=fWeight_Cent->GetNbinsX();i++)
   {
     Content = fCent_fromDATA->GetBinContent(i);
     error   = fCent_fromDATA->GetBinError(i);
     if(Content){
       weight  = maxCount/Content;
     }
     else{ 
       weight = 0;
     }
     fWeight_Cent->SetBinContent(i,weight);
     fWeight_Cent->SetBinError(i,0.0);
     //cout<<"cent "<<i-1<<"-"<<i<<" \t wgt = "<<fWeight_Cent->GetBinContent(i)<<" error = "<<fWeight_Cent->GetBinError(i)<<endl;
   }
 
   fListHistos->Add(fWeight_Cent);
   //delete fCent_fromDATA;



 //---------Filter bit efficiency----------
   if(fListHijing) {
    for(int i=0;i<10;i++) {
      fFB_Efficiency_Cent[i] = (TH1D *) fListHijing->FindObject(Form("eff_unbiased_%d",i));
      //cout<<"cent = "<<i<<" name = "<<fFB_Efficiency_Cent[i]->GetName()<<endl;
    }
   }
   else{ // if MC efficiency not found then use weight = 1. Code won't crash.
     printf("\n\n !!*****  Warning *****!!  \n TList for FilterBit efficiency not found !!\n\n");
     for(int i=0;i<10;i++){
      fFB_Efficiency_Cent[i] = new TH1D(Form("eff_unbiased_%d",i),"",100,0,20); 
      for(int j=1;j<=fFB_Efficiency_Cent[i]->GetNbinsX();j++){
         fFB_Efficiency_Cent[i]->SetBinContent(j,1.000);
       }
     } 
   }

  } //if(bApplyRecent)









  if(sCentEstimator=="V0"){
    fHist_Task_config->Fill(1.5);
  }
  if(sCentEstimator=="TPC"){
    fHist_Task_config->Fill(2.5);
  }
  if(sCentEstimator=="CL1"){
    fHist_Task_config->Fill(3.5);
  }
  if(fRejectPileUp){
    fHist_Task_config->Fill(4.5);
  }
  if(fRejectPileUpTight){
    fHist_Task_config->Fill(5.5);
  }
  if(fAnalysisSet=="FillGainEq"){
    fHist_Task_config->Fill(6.5);
  }
  if(fAnalysisSet=="DoGainEq"){
    fHist_Task_config->Fill(7.5);
    if(fListZDCWgt){
      fHist_Task_config->Fill(8.5);
    }
  }
  if(bFillCosSin){
    fHist_Task_config->Fill(9.5);
  }
  if(bApplyRecent){
    fHist_Task_config->Fill(10.5);
    if(fListZDCQxy){
      fHist_Task_config->Fill(11.5);
    }
  }
  if(!bRunAveragedQn){
    fHist_Task_config->Fill(12.5);
  }
  if(bCentCutShift){
    fHist_Task_config->Fill(13.5);
  }
  if(bApplyShiftCorr){
    fHist_Task_config->Fill(14.5);
  }
  if(bShiftCorrOnCent){
    fHist_Task_config->Fill(15.5);
  }
  if(bUseTrigonQn){
    fHist_Task_config->Fill(16.5);
  }


  fHist_CutParameters->SetBinContent(1,VzCut[0]);
  fHist_CutParameters->SetBinContent(2,VzCut[1]);
  fHist_CutParameters->SetBinContent(3,vzBin);
  fHist_CutParameters->SetBinContent(4,VxCut[0]);
  fHist_CutParameters->SetBinContent(5,VxCut[1]);
  fHist_CutParameters->SetBinContent(6,vxBin);
  fHist_CutParameters->SetBinContent(7,VyCut[0]);
  fHist_CutParameters->SetBinContent(8,VyCut[1]);
  fHist_CutParameters->SetBinContent(9,vyBin);






  PostData(1,fListHistos); 

  if(!bApplyRecent && fAnalysisSet=="DoGainEq") {
    PostData(2,fListZDCQxy); 
  }
  else{
    fListDummy1 = new TList();
    fListDummy1->SetOwner(kTRUE);
    fListDummy1->Add(fHist_Event_count);
    for(int i=0;i<90;i++){
      fListDummy1->Add(fHist_VZERO_Mult_vsRun[i]);
    }

    PostData(2,fListDummy1); 
  }


  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);


  fcheckOnce = 0;
  fOldRunNum = 0;
  printf("\n\n::UserCreateOutPutObject()\n Runflag= %d, dataset: %s, VxyBin = %d, Analysis= %s\n\n",frunflag,fDataSet.Data(),NbinVt,fAnalysisSet.Data());

}

//________________________________________________________________________
void AliAnalysisTaskZDCGainEq::UserExec(Option_t *)
{

  float stepCount = 0.5;

  //printf("\n ... ::UserExec() is being called. 1 Step %d...  \n",stepCount);

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
  fEvent           = dynamic_cast<AliFlowEventSimple*>(GetInputData(1));

  AliFlowEvent* anEvent = dynamic_cast<AliFlowEvent*>(GetInputData(2));

  AliFlowVector vQarray[2];

  if(anEvent) {
  // Get Q vectors for the subevents
   anEvent->GetZDC2Qsub(vQarray);
  }




  if(!aod){
    printf("\n ... ::UserExec = no aod found.....  \n");
    return;
  }

  fHist_Event_count->Fill(stepCount);
  stepCount++;


  AliAODVertex *pVertex = aod->GetPrimaryVertex();
  Double_t Vxyz[3]    =         {0,0,0};
  Vxyz[0]             = pVertex->GetX();
  Vxyz[1]             = pVertex->GetY();
  Vxyz[2]             = pVertex->GetZ();

  fHist_Vxy_RunAveraged->Fill(Vxyz[0], Vxyz[1]);

  //------- Apply Necessary event cuts ---------
  if(Vxyz[2] >= VzCut[1]  || Vxyz[2] <= VzCut[0])  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(Vxyz[0] >= VxCut[1]  || Vxyz[0] <= VxCut[0])  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(Vxyz[1] >= VyCut[1]  || Vxyz[1] <= VyCut[0])  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;


  Double_t EvtCent = fEvent->GetCentrality();


  //---------- get runindex: --------------
  Int_t runNumber = aod->GetRunNumber();
  Int_t runindex = -111;

 for(int i=0;i<frunflag;i++){
   if(runNumber==runNums[i])
    {
      runindex = i;
      break;
    }
  }
 if(runindex<0) {
    printf("\n ... **WARNING** \n::UserExec() runnumber not listed.\n EXIT..\n");
    //exit(1);
  }
 //-----------------------------------------

 //Get MagneticField:
 Double_t fMagField = aod->GetMagneticField();









 //--------- starting pileup rejection work: --------
   Double_t centrV0M=300;
   Double_t centrCL1=300;
   Double_t centrCL0=300;
   Double_t centrTRK=300;

  if(fDataSet=="2010"||fDataSet=="2011"){
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
  }// 2015

  Bool_t BisPileup=kFALSE;

  if(fRejectPileUp && InputEvent()) {
    //if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;
    if(fDataSet!="2015") {
          if(plpMV(aod)) {
            fPileUpCount->Fill(0.5);
            BisPileup=kTRUE;
          }
          Int_t isPileup = aod->IsPileupFromSPD(3);
          if(isPileup != 0) {
            fPileUpCount->Fill(1.5);
            //BisPileup=kTRUE;       //  ----- Rihan: pileup from SPD not used for 2010
          }
          if(((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
            fPileUpCount->Fill(2.5);
            BisPileup=kTRUE;
          }
          if(aod->IsIncompleteDAQ())  {
            fPileUpCount->Fill(3.5);
            BisPileup=kTRUE;
          }

    //check vertex consistency
     const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
     const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();

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
        }

        else {
          // pileup from AliMultSelection for 2015
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

      // pile-up a la Dobrin for LHC15o
          if(plpMV(aod)) {
            fPileUpCount->Fill(0.5);
            BisPileup=kTRUE;
          }

          Int_t isPileup = aod->IsPileupFromSPD(3);
          if(isPileup != 0) {
            fPileUpCount->Fill(1.5);
            BisPileup=kTRUE;          
          }

          if(((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
            fPileUpCount->Fill(2.5);
            BisPileup=kTRUE;
          }

          if(aod->IsIncompleteDAQ())  {
            fPileUpCount->Fill(3.5);
            BisPileup=kTRUE;
          }

          if(fabs(centrV0M-centrCL1)>7.5)  {
            fPileUpCount->Fill(4.5);
            BisPileup=kTRUE;
          }

     // check vertex consistency
     const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
     const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();

     if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
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

       // cuts on tracks
          const Int_t nTracks = aod->GetNumberOfTracks();
          Int_t multEsd = ((AliAODHeader*)aod->GetHeader())->GetNumberOfESDTracks();

          Int_t multTrk = 0;
          Int_t multTrkBefC = 0;
          Int_t multTrkTOFBefC = 0;
          Int_t multTPC = 0;

    for(Int_t it = 0; it < nTracks; it++) {
       AliAODTrack* aodTrk = (AliAODTrack*)aod->GetTrack(it);
        if(!aodTrk) {
          delete aodTrk;
          continue;
         }
//      if(aodTrk->TestFilterBit(32)){
//         multTrkBefC++;
//       if(TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
//         multTrkTOFBefC++;
//         if((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2) && (aodTrk->Pt() < 20.))
//            multTrk++;
//       }
         if(aodTrk->TestFilterBit(128))
              multTPC++;
        } // end of for (Int_t it = 0; it < nTracks; it++)

      Double_t multTPCn = multTPC;
      Double_t multEsdn = multEsd;
      Double_t multESDTPCDif = multEsdn - multTPCn*3.38;

     if(multESDTPCDif > (fRejectPileUpTight?700.:15000.)) {
        fPileUpCount->Fill(7.5);
        BisPileup=kTRUE;
      }

     if(fRejectPileUpTight) {
       if(BisPileup==kFALSE) {
              if(!fMultSelection->GetThisEventIsNotPileup()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
              if(BisPileup) fPileUpCount->Fill(8.5);
            }
          }
	}
      }

 //----------- pile up rejection done ---------
  if(BisPileup) {
    return;
  }

  fHist_Event_count->Fill(stepCount);
  stepCount++;







  if(runNumber!=fOldRunNum){  //if runNumber changed, re-read new list.
    fcheckOnce = 0;
  }

  if(!fcheckOnce && fAnalysisSet=="DoGainEq") {
    fHist_ChanWgt_ZDCC = (TH1F *) fListZDCWgt->FindObject(Form("fHist1F_ZDCC_ChannelWgt_Run%d",runNumber));
    fHist_ChanWgt_ZDCA = (TH1F *) fListZDCWgt->FindObject(Form("fHist1F_ZDCA_ChannelWgt_Run%d",runNumber));

    fcheckOnce++;
    fOldRunNum = runNumber;
  }



  /*
  if(!fcheckOnce && fAnalysisSet=="DoGainEq") {
    fHist_ChanWgt_ZDCC = (TH1F *) fListZDCWgt->FindObject(Form("fHist1F_ZDCC_ChannelWgt_Run%d",runNumber));
    fHist_ChanWgt_ZDCA = (TH1F *) fListZDCWgt->FindObject(Form("fHist1F_ZDCA_ChannelWgt_Run%d",runNumber));

    if(bApplyRecent) {
     if(!bRunAveragedQn) {
       if(fListZDCQxy) {
        for(int i=0; i<vzBin; i++) {
          fHist_Recenter_ZDCCx[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist2F_znCx_V0_Run%d_Vz%d",runNumber,i+1));
          fHist_Recenter_ZDCCy[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist2F_znCy_V0_Run%d_Vz%d",runNumber,i+1));
          fHist_Recenter_ZDCAx[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist2F_znAx_V0_Run%d_Vz%d",runNumber,i+1));
          fHist_Recenter_ZDCAy[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist2F_znAy_V0_Run%d_Vz%d",runNumber,i+1));
        }
       }
      }
     else {
       if(fListZDCQxy) {
        for(int i=0; i<vzBin; i++) {
          fHist_Recenter_ZDCCx[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist2F_znCx_V0_Run%d_Vz%d",0,i+1));
          fHist_Recenter_ZDCCy[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist2F_znCy_V0_Run%d_Vz%d",0,i+1));
          fHist_Recenter_ZDCAx[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist2F_znAx_V0_Run%d_Vz%d",0,i+1));
          fHist_Recenter_ZDCAy[i] = (TH2F *) fListZDCQxy->FindObject(Form("fHist2F_znAy_V0_Run%d_Vz%d",0,i+1));
        }
       }
     }
    }

    fcheckOnce++;
    fOldRunNum = runNumber;
  } */




  Double_t ChanWgtZDCC[5] = {1.,1.,1.,1.,1.};
  Double_t ChanWgtZDCA[5] = {1.,1.,1.,1.,1.};

  Int_t iCentBin = abs(EvtCent) + 1;
  Int_t iWgtBin = -1;
  Int_t iCommon = -1;

  
  if(fAnalysisSet=="DoGainEq") {
    if(fHist_ChanWgt_ZDCC && fHist_ChanWgt_ZDCA){

     for(int ich=1; ich<=5;  ich++){
        iWgtBin = 5*(iCentBin-1) + ich;
        ChanWgtZDCC[ich-1] = fHist_ChanWgt_ZDCC->GetBinContent(iWgtBin);
        ChanWgtZDCA[ich-1] = fHist_ChanWgt_ZDCA->GetBinContent(iWgtBin);
      }
    }
    else{
      //printf("\n\n **WARNING**\n ZDC Channel Weights not found. Using weights = 1.0 \n\n");
      //exit(1);
    }
  }



  fHist_Cent_woZDCcut->Fill(EvtCent);


  //----------- Read ZDC information ----------
  AliAODZDC *aodZDC  = aod->GetZDCData();
  Float_t                            fZDCGainAlpha = 0.500;  // fZDCGainAlpha : Jacopo uses 0.35 ??
  Float_t energyZNC  = (Float_t)  (aodZDC->GetZNCEnergy());
  Float_t energyZPC  = (Float_t)  (aodZDC->GetZPCEnergy());
  Float_t energyZNA  = (Float_t)  (aodZDC->GetZNAEnergy());
  Float_t energyZPA  = (Float_t)  (aodZDC->GetZPAEnergy());
  Float_t energyZEM1 = (Float_t) (aodZDC->GetZEM1Energy());
  Float_t energyZEM2 = (Float_t) (aodZDC->GetZEM2Energy());

  const Double_t * towZNC   = aodZDC->GetZNCTowerEnergy();
  const Double_t * towZPC   = aodZDC->GetZPCTowerEnergy();
  const Double_t * towZNA   = aodZDC->GetZNATowerEnergy();
  const Double_t * towZPA   = aodZDC->GetZPATowerEnergy();

  const Double_t * towZNClg = aodZDC->GetZNCTowerEnergyLR(); // Low gain something, should not be used.
  const Double_t * towZNAlg = aodZDC->GetZNATowerEnergyLR();

  Double_t towZPClg[5] = {0.,};
  Double_t towZPAlg[5] = {0.,};

  for(Int_t it=0; it<5; it++) {
    towZPClg[it] = 8*towZPC[it];
    towZPAlg[it] = 8*towZNA[it];
  }

  Int_t BadChannel  = 0;

  //sanity: remove if any of ZDC_C_A has negetive Energy: 
  //This was a bad choice of cut. I need to see the Physics effect later..!! 
  //if(towZNC[1]<0 || towZNC[2]<0 || towZNC[3]<0 || towZNC[4]<0)  return; 

  for(int i=0; i<4; i++) {
    if(towZNC[i] < 0) {
       BadChannel++;
    }
  }

  if(BadChannel>=2)  return; // Remove Event if more than one channel has Energy < 0 

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  //if(towZNA[1]<0 || towZNA[2]<0 || towZNA[3]<0 || towZNA[4]<0)  return; 
  BadChannel  = 0;

  for(int i=0; i<4; i++) {
    if(towZNA[i] < 0) {
       BadChannel++;
    }
  }

  if(BadChannel>=2)  return; 

  fHist_Event_count->Fill(stepCount);
  stepCount++;












//********** Get centroid from ZDCs **************

  Double_t xyZNC[2]={0.,0.};
  Double_t xyZNA[2]={0.,0.};

  Float_t zncEnergy=0., znaEnergy=0.;

  

/*-----------------------Not used---------------------------------
  Int_t CenBin = GetCenBin(centrperc);
  Double_t zvtxpos[3]={0.,0.,0.};
  fFlowEvent->GetVertexPosition(zvtxpos);
  Int_t RunNum=fFlowEvent->GetRun();
  if(fTowerEqList) {
   if(RunNum!=fCachedRunNum) {
      for(Int_t i=0; i<8; i++) {
      fTowerGainEq[i] = (TH1D*)(fTowerEqList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fhnTowerGainEqFactor[%d][%d]",RunNum,i)));
     }
   }
 }
  Bool_t fUseMCCen = kFALSE; //rihan:hardcoded
  if (fUseMCCen) {
   if(aod->GetRunNumber() < 209122)
    aodZDC->GetZNCentroidInPbPb(1380., xyZNC, xyZNA);
    else
    aodZDC->GetZNCentroidInPbPb(2510., xyZNC, xyZNA);
  }
  else {
  //set tower gain equalization, if available
  if(fTowerEqList) {
   for(Int_t i=0; i<8; i++)
   {
    if(fTowerGainEq[i])
    AvTowerGain[i] = fTowerGainEq[i]->GetBinContent(fTowerGainEq[i]->FindBin(centrperc));
   }
 }//---------------------------------------------------------------  */



  Double_t  towCalibZNC[5] = {0,}; 
  Double_t  towCalibZNA[5] = {0,};

  // towZNC[] is constant; so Need to make a copy
  for(int i=0;i<5;i++){
    towCalibZNC[i] = towZNC[i];
    towCalibZNA[i] = towZNA[i];
  }

  // If Ei < 0 , Ei = 0 (Maxim's idea)
  for(int i=1;i<5;i++){
    if(towCalibZNC[i] < 0) towCalibZNC[i] = 0;
    if(towCalibZNA[i] < 0) towCalibZNA[i] = 0;    
  }



  Int_t indexVx = fHist_Vx_ArrayFinder->FindBin(Vxyz[0]);
  Int_t indexVy = fHist_Vy_ArrayFinder->FindBin(Vxyz[1]);
  Int_t indexVz = fHist_Vz_ArrayFinder->FindBin(Vxyz[2]);


  Double_t tVertexBin1 = 0;
  tVertexBin1 = (Double_t) (indexVy-1)*vxBin + (Double_t)indexVx - 0.5 ; 


  if(fAnalysisSet=="FillGainEq") {
    fHist_Vx_vs_runnum->Fill(runindex,Vxyz[0]);
    fHist_Vy_vs_runnum->Fill(runindex,Vxyz[1]);
    fHist_Vz_vs_runnum->Fill(runindex,Vxyz[2]);

    fHist_ZDCC_En_CommonCh[indexVz-1]->Fill(EvtCent,tVertexBin1,towCalibZNC[0]);
    fHist_ZDCA_En_CommonCh[indexVz-1]->Fill(EvtCent,tVertexBin1,towCalibZNA[0]);

    for(int ich=0; ich<5; ich++) {
       fHist_ZDCC_En_Run[runindex]->Fill(EvtCent,ich,towCalibZNC[ich]);
       fHist_ZDCA_En_Run[runindex]->Fill(EvtCent,ich,towCalibZNA[ich]);
    }
  }







  // Now calibrate the energy of channel [0-4]:
  for(int i=0;i<5;i++){
    towCalibZNC[i] = ChanWgtZDCC[i]*towCalibZNC[i];

    if(ChanWgtZDCA[i] < 4.){
      towCalibZNA[i] = ChanWgtZDCA[i]*towCalibZNA[i];
    }
  }

 //manually get Energy in ZDC-A channel [2]:
  if(ChanWgtZDCA[2] >= 4.0){
    towCalibZNA[2] = towCalibZNA[0] - towCalibZNA[1] - towCalibZNA[3] - towCalibZNA[4];
  }


  for(Int_t i=0; i<5; i++){
    zncEnergy += towCalibZNC[i];
    znaEnergy += towCalibZNA[i];
  }


  Double_t AvTowerGain[8] = {1., 1., 1., 1., 1., 1., 1., 1.};
  
  const Float_t x[4] = {-1.75,  1.75,-1.75, 1.75};
  const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};

  Float_t numXZNC=0., numYZNC=0., denZNC=0., wZNC;
  Float_t numXZNA=0., numYZNA=0., denZNA=0., wZNA;

  for(Int_t i=0; i<4; i++) {
    if(towCalibZNC[i+1]>0.) {
       wZNC = TMath::Power(towCalibZNC[i+1], fZDCGainAlpha)*AvTowerGain[i];
       numXZNC += x[i]*wZNC;
       numYZNC += y[i]*wZNC;
       denZNC  += wZNC;
    }

    if(towCalibZNA[i+1]>0.) {
       wZNA = TMath::Power(towCalibZNA[i+1], fZDCGainAlpha)*AvTowerGain[i+4];
       numXZNA += x[i]*wZNA;
       numYZNA += y[i]*wZNA;
       denZNA  += wZNA;
    }
  }

  if(denZNC!=0) {
    xyZNC[0] = numXZNC/denZNC;
    xyZNC[1] = numYZNC/denZNC;
  }
  else{
     xyZNC[0]  = 999.;
     xyZNC[1]  = 999.;
     zncEnergy =   0.;
  }
  if(denZNA!=0) {
     xyZNA[0] = numXZNA/denZNA;
     xyZNA[1] = numYZNA/denZNA;
  }
  else{
     xyZNA[0]  = 999.;
     xyZNA[1]  = 999.;
     znaEnergy =   0.;
  }


  
  xyZNA[0] = -1.*xyZNA[0]; //----- Important: zdcA_X = -zdcA_X ---------


  /*
  if(sqrt(xyZNC[0]*xyZNC[0] + xyZNC[1]*xyZNC[1])>1.5)  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(sqrt(xyZNA[0]*xyZNA[0] + xyZNA[1]*xyZNA[1])>1.5)  return;
  */

  fHist_Event_count->Fill(stepCount);
  stepCount++;
 
  Double_t psi1C = TMath::ATan2(xyZNC[1],xyZNC[0]);
  if(psi1C<0){
     psi1C += 2.*TMath::Pi();
  }
  Double_t psi1A = TMath::ATan2(xyZNA[1],xyZNA[0]);
  if(psi1A<0){
     psi1A += 2.*TMath::Pi();
  } 

  /*if(EvtCent>=5 && EvtCent<=45) {
    fHist_Psi1_ZDCC_wGainCorr->Fill(psi1C);
    fHist_Psi1_ZDCA_wGainCorr->Fill(psi1A);
  }*/

  if(fAnalysisSet=="DoGainEq") {

    tVertexBin1 = (Double_t) (indexVy-1)*vxBin + (Double_t)indexVx - 0.5 ; 

    if(!bApplyRecent) {
     if(bFillCosSin) {
       if(!bRunAveragedQn){
         fHist_znCx_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Cos(psi1C)); 
         fHist_znCy_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Sin(psi1C));
         fHist_znAx_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Cos(psi1A));
         fHist_znAy_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Sin(psi1A));
       }
       else{
         fHist_znCx_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Cos(psi1C)); 
         fHist_znCy_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Sin(psi1C));
         fHist_znAx_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Cos(psi1A));
         fHist_znAy_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Sin(psi1A));
       }
     }
     else{
       if(!bRunAveragedQn){
         fHist_znCx_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNC[0]); 
         fHist_znCy_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNC[1]);
         fHist_znAx_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNA[0]);
         fHist_znAy_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNA[1]);
       }
       else{
         fHist_znCx_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNC[0]); 
         fHist_znCy_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNC[1]);
         fHist_znAx_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNA[0]);
         fHist_znAy_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNA[1]);
       }
     }
    }
  }


  Double_t fRefMult = (Double_t) fEvent->GetReferenceMultiplicity();;

  if(fAnalysisSet!="FillGainEq" && bFillCosSin) {
    xyZNC[0] = TMath::Cos(psi1C);
    xyZNC[1] = TMath::Sin(psi1C);
    xyZNA[0] = TMath::Cos(psi1A);
    xyZNA[1] = TMath::Sin(psi1A);
  }



  /*  if(bFillZDCQAon){

  //Double_t  FillVsWith[5]  = {EvtCent,static_cast<Double_t>(nRefMult), Vxyz[0], Vxyz[1], Vxyz[2]};
    Double_t  FillVsWith[5]  = {EvtCent, fRefMult, Vxyz[0], Vxyz[1], Vxyz[2]};
    Double_t  FillValueQx[4] = {xyZNA[0],xyZNC[0],xyZNA[1],xyZNC[1]};
    Double_t  FillValueXX[4] = {xyZNA[0]*xyZNC[0],xyZNA[1]*xyZNC[1],xyZNC[0]*xyZNA[1],xyZNC[1]*xyZNA[0]}; //XaXc,YaYc,XcYa,YcXa

    for(int i=0;i<4;i++){
     for(int j=0;j<5;j++){
        if(j>0 && (EvtCent<5 || EvtCent>45)) continue;
        fHist_Qx_vs_Obs_woCorr[i][j]->Fill(FillVsWith[j],FillValueQx[i]);
        fHist_XX_vs_Obs_woCorr[i][j]->Fill(FillVsWith[j],FillValueXX[i]);
     }
    }
    } */


  if(fAnalysisSet=="FillGainEq" && bFillCosSin) {
    xyZNC[0] = TMath::Cos(psi1C);
    xyZNC[1] = TMath::Sin(psi1C);
    xyZNA[0] = TMath::Cos(psi1A);
    xyZNA[1] = TMath::Sin(psi1A);
  }


  Double_t meanCx = 0.;
  Double_t meanCy = 0.;
  Double_t meanAx = 0.;
  Double_t meanAy = 0.;

  
  //Apply the recentering:
  /*
  if(bApplyRecent) {
   Int_t tVertexBin2 = (indexVy-1)*vxBin + indexVx; 
    if(fListZDCQxy) {
      meanCx = fHist_Recenter_ZDCCx[indexVz-1]->GetBinContent(tVertexBin2,iCentBin);
      meanCy = fHist_Recenter_ZDCCy[indexVz-1]->GetBinContent(tVertexBin2,iCentBin);
      meanAx = fHist_Recenter_ZDCAx[indexVz-1]->GetBinContent(tVertexBin2,iCentBin);
      meanAy = fHist_Recenter_ZDCAy[indexVz-1]->GetBinContent(tVertexBin2,iCentBin);
    }
  }  

  xyZNC[0] = xyZNC[0] - meanCx;
  xyZNC[1] = xyZNC[1] - meanCy;
  xyZNA[0] = xyZNA[0] - meanAx;
  xyZNA[1] = xyZNA[1] - meanAy;
  */


  const AliAODVZERO *fvzero = aod->GetVZEROData();

  Float_t fV0Mult[64] = {0.,}; //0-31 v0c and 32-63 is v0A

  for(int i=0;i<64;i++){
    fV0Mult[i] = fvzero->GetMultiplicity(i); 
  }

  for(int i=0;i<64;i++){
   fHist_VZERO_Mult_vsRun[runindex]->Fill(EvtCent,i+0.5,fV0Mult[i]);
  }


  //use Jacopo's flowEvent ZDC Q-vectors:

  xyZNC[0] = vQarray[0].Px();
  xyZNC[1] = vQarray[0].Py();
  xyZNA[0] = vQarray[1].Px();
  xyZNA[1] = vQarray[1].Py();
  


  Double_t Qvect_ModC = TMath::Sqrt(xyZNC[0]*xyZNC[0] + xyZNC[1]*xyZNC[1]); 
  Double_t Qvect_ModA = TMath::Sqrt(xyZNA[0]*xyZNA[0] + xyZNA[1]*xyZNA[1]); 

  Double_t dTermXXYY = xyZNA[0]*xyZNC[0]  - xyZNA[1]*xyZNC[1];
  Double_t dTermXYXY = xyZNA[0]*xyZNC[1]  + xyZNA[1]*xyZNC[0];

  Double_t Qvect_ModSumAC = TMath::Sqrt(dTermXXYY*dTermXXYY + dTermXYXY*dTermXYXY); 



  fHist_Event_count->Fill(stepCount);
  stepCount++;

  Double_t Psi1C = TMath::ATan2(xyZNC[1],xyZNC[0]);
  if(Psi1C<0){
   Psi1C += 2.*TMath::Pi();
  }

  Double_t Psi1A = TMath::ATan2(xyZNA[1],xyZNA[0]);
  if(Psi1A<0){
   Psi1A += 2.*TMath::Pi();
  }
 
  Double_t PsiSumAC = Psi1C+Psi1A;

  if(PsiSumAC>2.*TMath::Pi()){
   PsiSumAC -= 2.*TMath::Pi();
  }


  fHist_PsiZDCA_ZeroQn->Fill(Psi1A);
  fHist_PsiZDCC_ZeroQn->Fill(Psi1C);
  fHist_PsiSumAC_ZeroQn->Fill(PsiSumAC);


  if(Qvect_ModC < 1.e-6 || Qvect_ModA < 1.e-6)   return;
  //std::cout<<"Cx = "<<xyZNC[0]<<"\tCy = "<<xyZNC[1]<<"\tAx = "<<xyZNA[0]<<"\tAy = "<<xyZNA[1]<<std::endl;


  fHist_Event_count->Fill(stepCount);
  stepCount++;

  //Shift Histograms in Vxy binning:
  //if(!bCentCutShift){  }
  //else if(bCentCutShift && EvtCent>=0 && EvtCent< 60){

  if(bCentCutShift && EvtCent>=0 && EvtCent< 60){
    for(int i=0;i<4;i++){
      fHist_ZDCC_AvgCos_VsRun[i]->Fill(EvtCent,runindex,TMath::Cos((i+1)*Psi1C));
      fHist_ZDCA_AvgCos_VsRun[i]->Fill(EvtCent,runindex,TMath::Cos((i+1)*Psi1A));
      fHist_ZDCC_AvgSin_VsRun[i]->Fill(EvtCent,runindex,TMath::Sin((i+1)*Psi1C));
      fHist_ZDCA_AvgSin_VsRun[i]->Fill(EvtCent,runindex,TMath::Sin((i+1)*Psi1A));

      fHist_ZDC_AvgCosPsiSum_VsRun[i]->Fill(EvtCent,runindex,TMath::Cos((i+1)*(Psi1C+Psi1A)));
      fHist_ZDC_AvgSinPsiSum_VsRun[i]->Fill(EvtCent,runindex,TMath::Sin((i+1)*(Psi1C+Psi1A)));

      fHist_ZDC_AvgCosPsiDif_VsRun[i]->Fill(EvtCent,runindex,TMath::Cos((i+1)*(Psi1C-Psi1A)));
      fHist_ZDC_AvgSinPsiDif_VsRun[i]->Fill(EvtCent,runindex,TMath::Sin((i+1)*(Psi1C-Psi1A)));
    }    
  }


  if(bFillZDCQAon){

    fHist_XXYY_vs_Cent_woCorr[0]->Fill(EvtCent, (xyZNA[0]*xyZNC[0] - xyZNA[1]*xyZNC[1]));
    fHist_XXYY_vs_Cent_woCorr[1]->Fill(EvtCent, (xyZNA[0]*xyZNC[1] + xyZNA[1]*xyZNC[0]));

  //Double_t  FillVsWith[5]  = {EvtCent,static_cast<Double_t>(nRefMult), Vxyz[0], Vxyz[1], Vxyz[2]};
    Double_t  FillVsWith[5]  = {EvtCent, fRefMult, Vxyz[0], Vxyz[1], Vxyz[2]};
    Double_t  FillValueQx[4] = {xyZNA[0],xyZNC[0],xyZNA[1],xyZNC[1]};
    Double_t  FillValueXX[4] = {xyZNA[0]*xyZNC[0],xyZNA[1]*xyZNC[1],xyZNC[0]*xyZNA[1],xyZNC[1]*xyZNA[0]}; //XaXc,YaYc,XcYa,YcXa

    Double_t  FillValueQxTrig[4] = {TMath::Cos(Psi1A),TMath::Cos(Psi1C),TMath::Sin(Psi1A),TMath::Sin(Psi1C)};
    Double_t  FillValueXXTrig[4] = {TMath::Cos(Psi1A)*TMath::Cos(Psi1C),TMath::Sin(Psi1A)*TMath::Sin(Psi1C),TMath::Cos(Psi1C)*TMath::Sin(Psi1A),TMath::Cos(Psi1A)*TMath::Sin(Psi1C)};

    Int_t centIndex = fHist_NormalCentralityBins->FindBin(EvtCent) - 1;
    if(centIndex<6) {
      for(int i=0;i<4;i++){
	if(fMagField>0){
	 fHist_XX_vs_QnC_2DwoCorr_PosMag[i][centIndex]->Fill(xyZNC[0],xyZNC[1],FillValueXX[i]);
	 fHist_XX_vs_QnA_2DwoCorr_PosMag[i][centIndex]->Fill(xyZNA[0],xyZNA[1],FillValueXX[i]);
	}
	else if(fMagField<0){
	 fHist_XX_vs_QnC_2DwoCorr_NegMag[i][centIndex]->Fill(xyZNC[0],xyZNC[1],FillValueXX[i]);
	 fHist_XX_vs_QnA_2DwoCorr_NegMag[i][centIndex]->Fill(xyZNA[0],xyZNA[1],FillValueXX[i]);
	}
      }
    }

    for(int i=0;i<4;i++){

      fHist_Qx_Trig_woCorr[i]->Fill(EvtCent,FillValueQxTrig[i]);
      fHist_XX_Trig_woCorr[i]->Fill(EvtCent,FillValueXXTrig[i]);

     for(int j=0;j<5;j++){
        if(j>0 && (EvtCent<5 || EvtCent>45)) continue;
        fHist_Qx_vs_Obs_woCorr[i][j]->Fill(FillVsWith[j],FillValueQx[i]);
        fHist_XX_vs_Obs_woCorr[i][j]->Fill(FillVsWith[j],FillValueXX[i]);
     }
    }
  }


  if(!bCentCutShift){
    fHist_Psi1_ZDCC_wGainCorr->Fill(psi1C);
    fHist_Psi1_ZDCA_wGainCorr->Fill(psi1A);
    fHist_PsiSumAC_woCorr->Fill(PsiSumAC);
  }
  else if(bCentCutShift && EvtCent>=5 && EvtCent<=45){
    fHist_Psi1_ZDCC_wGainCorr->Fill(psi1C);
    fHist_Psi1_ZDCA_wGainCorr->Fill(psi1A);
    fHist_PsiSumAC_woCorr->Fill(PsiSumAC);
  }

  //----------- Apply Shift correction: ----------
  Double_t ShiftCosC[4] = {0.,};
  Double_t ShiftSinC[4] = {0.,};

  Double_t ShiftCosA[4] = {0.,};
  Double_t ShiftSinA[4] = {0.,};

  if(bApplyShiftCorr) {
   if(!bShiftCombinedEP) {
     //Run dependent Shift Correction:
     for(int i=0;i<4;i++) {
       ShiftCosC[i] = fHist_ZDCC_AvgCosNPsi[i]->GetBinContent(iCentBin,runindex+1);
       ShiftSinC[i] = fHist_ZDCC_AvgSinNPsi[i]->GetBinContent(iCentBin,runindex+1); 
       ShiftCosA[i] = fHist_ZDCA_AvgCosNPsi[i]->GetBinContent(iCentBin,runindex+1);
       ShiftSinA[i] = fHist_ZDCA_AvgSinNPsi[i]->GetBinContent(iCentBin,runindex+1); 
     }
   }
   else if(bShiftCombinedEP) {
     for(int i=0;i<4;i++) {
        ShiftCosC[i] = fHist_ZDC_AvgCosNPsiAC[i]->GetBinContent(iCentBin,runindex+1);
        ShiftSinC[i] = fHist_ZDC_AvgSinNPsiAC[i]->GetBinContent(iCentBin,runindex+1);
     }
   }
  }

  Double_t shiftPsiC = 0;
  Double_t shiftPsiA = 0;

  PsiSumAC = Psi1C+Psi1A;

  if(bApplyShiftCorr) {
   if(!bShiftCombinedEP) {
     for(int i=0;i<4;i++) {
       shiftPsiC += 2.0/(i+1)*(-1.0*ShiftSinC[i]*TMath::Cos((i+1)*Psi1C) + ShiftCosC[i]*TMath::Sin((i+1)*Psi1C));
       shiftPsiA += 2.0/(i+1)*(-1.0*ShiftSinA[i]*TMath::Cos((i+1)*Psi1A) + ShiftCosA[i]*TMath::Sin((i+1)*Psi1A));
     }
     Psi1C += shiftPsiC;
     Psi1A += shiftPsiA;
   }
   else if(bShiftCombinedEP) {
     shiftPsiC = 0;
     for(int i=0;i<4;i++) {
       shiftPsiC += 2.0/(i+1)*(-1.0*ShiftSinC[i]*TMath::Cos((i+1)*PsiSumAC) + ShiftCosC[i]*TMath::Sin((i+1)*PsiSumAC));
     }
     PsiSumAC += shiftPsiC;
   }
  }


  if(Psi1C<0){
   Psi1C += 2.*TMath::Pi();
  }
  if(Psi1A<0){
   Psi1A += 2.*TMath::Pi();
  }
  if(PsiSumAC<0){
   PsiSumAC += 2.*TMath::Pi();
  }
  else if(PsiSumAC>2.*TMath::Pi()){
   PsiSumAC -= 2.*TMath::Pi();
  }



  //reconstruct the new ZDC-Qvect (with original Modulas).
  if(bApplyShiftCorr) {
   if(!bUseTrigonQn){
    xyZNC[0] = Qvect_ModC*TMath::Cos(Psi1C);
    xyZNC[1] = Qvect_ModC*TMath::Sin(Psi1C);
    xyZNA[0] = Qvect_ModA*TMath::Cos(Psi1A);
    xyZNA[1] = Qvect_ModA*TMath::Sin(Psi1A);
    //fill Qn Recentering histograms:
    fHist_ZDCC_AvgQx_VsRun->Fill(EvtCent,runindex,xyZNC[0]);
    fHist_ZDCC_AvgQy_VsRun->Fill(EvtCent,runindex,xyZNC[1]);
    fHist_ZDCA_AvgQx_VsRun->Fill(EvtCent,runindex,xyZNA[0]);
    fHist_ZDCA_AvgQy_VsRun->Fill(EvtCent,runindex,xyZNA[1]);
    fHist_ZDC_AvgXXminusYY_VsRun->Fill(EvtCent,runindex,(xyZNC[0]*xyZNA[0] - xyZNC[1]*xyZNA[1]));
    fHist_ZDC_AvgXYplusXY_VsRun->Fill (EvtCent,runindex,(xyZNC[1]*xyZNA[0] + xyZNC[0]*xyZNA[1]));

    dTermXXYY = Qvect_ModSumAC*TMath::Cos(PsiSumAC);
    dTermXYXY = Qvect_ModSumAC*TMath::Sin(PsiSumAC);
    fHist_ZDC_dTermXXYY_VsRun->Fill(EvtCent,runindex,dTermXXYY);
    fHist_ZDC_dTermXYXY_VsRun->Fill(EvtCent,runindex,dTermXYXY);
   }
   else if(bUseTrigonQn){//reconstruct the new ZDC-Qvect with Qn:[-1,1] 
    xyZNC[0] =  TMath::Cos(Psi1C);
    xyZNC[1] =  TMath::Sin(Psi1C);
    xyZNA[0] =  TMath::Cos(Psi1A);
    xyZNA[1] =  TMath::Sin(Psi1A);

    fHist_ZDCC_AvgQx_VsRun->Fill(EvtCent,runindex,xyZNC[0]);
    fHist_ZDCC_AvgQy_VsRun->Fill(EvtCent,runindex,xyZNC[1]);
    fHist_ZDCA_AvgQx_VsRun->Fill(EvtCent,runindex,xyZNA[0]);
    fHist_ZDCA_AvgQy_VsRun->Fill(EvtCent,runindex,xyZNA[1]);
    fHist_ZDC_AvgXXminusYY_VsRun->Fill(EvtCent,runindex,(xyZNC[0]*xyZNA[0] - xyZNC[1]*xyZNA[1]));
    fHist_ZDC_AvgXYplusXY_VsRun->Fill (EvtCent,runindex,(xyZNC[1]*xyZNA[0] + xyZNC[0]*xyZNA[1]));

    dTermXXYY = TMath::Cos(PsiSumAC);
    dTermXYXY = TMath::Sin(PsiSumAC);
    fHist_ZDC_dTermXXYY_VsRun->Fill(EvtCent,runindex,dTermXXYY);
    fHist_ZDC_dTermXYXY_VsRun->Fill (EvtCent,runindex,dTermXYXY);
   }
  }

  //do recentering after shift:
  Double_t CorrC[2] = {0.,};
  Double_t CorrA[2] = {0.,};

 if(bApplyShiftCorr) {
   if(!bShiftCombinedEP) {
     xyZNC[0] -= fHist_ZDCC_AvgQx_VsCR->GetBinContent(iCentBin,runindex+1);
     xyZNC[1] -= fHist_ZDCC_AvgQy_VsCR->GetBinContent(iCentBin,runindex+1);
     xyZNA[0] -= fHist_ZDCA_AvgQx_VsCR->GetBinContent(iCentBin,runindex+1);
     xyZNA[1] -= fHist_ZDCA_AvgQy_VsCR->GetBinContent(iCentBin,runindex+1);
   }
   if(bShiftCombinedEP){
     dTermXXYY -= fHist_ZDC_dTermXXYY_VsCR->GetBinContent(iCentBin,runindex+1);
     dTermXYXY -= fHist_ZDC_dTermXYXY_VsCR->GetBinContent(iCentBin,runindex+1);
   }
 }




  if(bFillZDCQAon){

    //fHist_XXYY_vs_Cent_wiCorr[0]->Fill(EvtCent, xyZNC[0]*xyZNA[0] - xyZNC[1]*xyZNA[1]);
    //fHist_XXYY_vs_Cent_wiCorr[1]->Fill(EvtCent, xyZNC[1]*xyZNA[0] + xyZNC[0]*xyZNA[1]);
    fHist_XXYY_vs_Cent_wiCorr[0]->Fill(EvtCent, dTermXXYY);
    fHist_XXYY_vs_Cent_wiCorr[1]->Fill(EvtCent, dTermXYXY);

    Double_t  FillVsWithNew[5]  = {EvtCent, fRefMult, Vxyz[0], Vxyz[1], Vxyz[2]};
    Double_t  FillValueQxNew[4] = {xyZNA[0],xyZNC[0],xyZNA[1],xyZNC[1]};
    Double_t  FillValueXXNew[4] = {xyZNA[0]*xyZNC[0],xyZNA[1]*xyZNC[1],xyZNC[0]*xyZNA[1],xyZNC[1]*xyZNA[0]}; //XaXc,YaYc,XcYa,YcXa

    Double_t  FillValueQxTrigNew[4] = {TMath::Cos(Psi1A),TMath::Cos(Psi1C),TMath::Sin(Psi1A),TMath::Sin(Psi1C)};
    Double_t  FillValueXXTrigNew[4] = {TMath::Cos(Psi1A)*TMath::Cos(Psi1C),TMath::Sin(Psi1A)*TMath::Sin(Psi1C),TMath::Cos(Psi1C)*TMath::Sin(Psi1A),TMath::Cos(Psi1A)*TMath::Sin(Psi1C)};

    for(int i=0;i<4;i++){
      if(EvtCent>=5 && EvtCent<=45){
        fHist_Qx_wiCorr_RunByRun[i]->Fill(runindex,FillValueQxNew[i]);
      }

      fHist_Qx_Trig_wiCorr[i]->Fill(EvtCent,FillValueQxTrigNew[i]);
      fHist_XX_Trig_wiCorr[i]->Fill(EvtCent,FillValueXXTrigNew[i]);

      for(int j=0;j<5;j++){
        if(j>0 && (EvtCent<5 || EvtCent>45)) continue;
        fHist_Qx_vs_Obs_wiCorr[i][j]->Fill(FillVsWithNew[j],FillValueQxNew[i]);
        fHist_XX_vs_Obs_wiCorr[i][j]->Fill(FillVsWithNew[j],FillValueXXNew[i]);
     }
    }
  }




  Int_t iTracks = fEvent->NumberOfTracks();
  AliFlowTrackSimple*   pTrack = NULL;
  Double_t         Qnx_TPC[4]  = {0,};
  Double_t         Qny_TPC[4]  = {0,};
  Double_t         psi2,dPhi,dPt,dEta;
  Double_t               pTwgt  = 1.0;
  Double_t               npoiMult = 0;
  Double_t                    dUx = 0;
  Double_t                    dUy = 0;
  Int_t                    ipTBin = 1;
  Int_t                   cIndex = -1;
  Double_t               nRefMult = 0;

  if(EvtCent<5.0) { cIndex  =     0; }
  else if(EvtCent>=5.0 && EvtCent<10){
   cIndex  = 1;
  }
  else{
   cIndex  =  abs(EvtCent/10.0)  +  1;
  }

  //Use Jacopo's notation:
  Double_t    Qtot, QRe, QIm;
  Double_t ZARe =   xyZNA[0];
  Double_t ZAIm =   xyZNA[1];
  Double_t ZCRe =   xyZNC[0];
  Double_t ZCIm =   xyZNC[1];


  Double_t fullTerm =    0;
  Double_t fullReso =    0; 
  Double_t CentWgt  =  1.0;
  Int_t dChrg       =    0;
  Double_t Qxpos[2] = {0.,};  //0=1*phi,1=2*phi
  Double_t Qypos[2] = {0.,};
  Double_t Qxneg[2] = {0.,};  //0=1*phi,1=2*phi
  Double_t Qyneg[2] = {0.,};

  Double_t QxEtaNeg[2] = {0.,};  //0=1*phi,1=2*phi
  Double_t QyEtaNeg[2] = {0.,};
  Double_t QxEtaPos[2] = {0.,};  //0=1*phi,1=2*phi
  Double_t QyEtaPos[2] = {0.,};

  Double_t mPosPoi = 0;
  Double_t mNegPoi = 0;
  Double_t mSubEtaPos = 0;
  Double_t mSubEtaNeg = 0;


  if(bApplyRecent && bApplyShiftCorr) {
   for(int i=0; i<iTracks; i++) {
     pTrack    =    fEvent->GetTrack(i);
     if(!pTrack)               continue;
     dPhi      =          pTrack->Phi();
     dPt       =          pTrack-> Pt();
     dEta      =          pTrack->Eta();
     dChrg     =       pTrack->Charge();
     if(fabs(dEta) > 0.8)      continue;
     if(dPt<0.20 || dPt>50.0)  continue;
     if(!pTrack->IsPOItype(1)) continue;
     nRefMult++;

     ipTBin = fFB_Efficiency_Cent[cIndex]->FindBin(dPt);
     pTwgt  = 1.0/fFB_Efficiency_Cent[cIndex]->GetBinContent(ipTBin);

     if(dEta>=0) {
      if(dChrg<0){
        Qxneg[0] += pTwgt*TMath::Cos(1.*dPhi);  //no Eta/Phi weight at this moment, improve later.
        Qyneg[0] += pTwgt*TMath::Sin(1.*dPhi);
        Qxneg[1] += pTwgt*pTwgt*TMath::Cos(2.*dPhi); //for auto-correlation
        Qyneg[1] += pTwgt*pTwgt*TMath::Sin(2.*dPhi);
        mNegPoi  += pTwgt;
      }
      else if(dChrg>0) {
        Qxpos[0] += pTwgt*TMath::Cos(1.*dPhi);  //no  Eta/Phi weight at this moment, improve later.
        Qypos[0] += pTwgt*TMath::Sin(1.*dPhi);
        Qxpos[1] += pTwgt*pTwgt*TMath::Cos(2.*dPhi); //for auto-correlation
        Qypos[1] += pTwgt*pTwgt*TMath::Sin(2.*dPhi);
        mPosPoi  += pTwgt;
      }

      QxEtaPos[0] += pTwgt*TMath::Cos(2.*dPhi);
      QyEtaPos[0] += pTwgt*TMath::Sin(2.*dPhi);
      QxEtaPos[1] += pTwgt*TMath::Cos(4.*dPhi);
      QyEtaPos[1] += pTwgt*TMath::Sin(4.*dPhi);
      mSubEtaPos  += pTwgt;
     }

     if(dEta<0) {
      QxEtaNeg[0] += pTwgt*TMath::Cos(2.*dPhi);
      QyEtaNeg[0] += pTwgt*TMath::Sin(2.*dPhi);
      QxEtaNeg[1] += pTwgt*TMath::Cos(4.*dPhi);
      QyEtaNeg[1] += pTwgt*TMath::Sin(4.*dPhi);
      mSubEtaNeg  += pTwgt;
     }


     dUx   =     TMath::Cos(dPhi);
     dUy   =     TMath::Sin(dPhi);
     Qnx_TPC[0] += dUx*pTwgt; 
     Qny_TPC[0] += dUy*pTwgt;

    
     fHist_v1xV1_ZDN_EtaDiff[0][cIndex]->Fill(dEta, dUx*xyZNC[0], pTwgt); //uxCx
     fHist_v1xV1_ZDN_EtaDiff[1][cIndex]->Fill(dEta, dUy*xyZNC[1], pTwgt); //uyCy
     fHist_v1xV1_ZDN_EtaDiff[2][cIndex]->Fill(dEta, dUx*xyZNA[0], pTwgt); //uxAx
     fHist_v1xV1_ZDN_EtaDiff[3][cIndex]->Fill(dEta, dUy*xyZNA[1], pTwgt); //uyAy
     fHist_v1xV1_ZDN_pTDiff[0][cIndex]->Fill(dPt, dUx*xyZNC[0], pTwgt); //uxCx
     fHist_v1xV1_ZDN_pTDiff[1][cIndex]->Fill(dPt, dUy*xyZNC[1], pTwgt); //uyCy
     fHist_v1xV1_ZDN_pTDiff[2][cIndex]->Fill(dPt, dUx*xyZNA[0], pTwgt); //uxAx
     fHist_v1xV1_ZDN_pTDiff[3][cIndex]->Fill(dPt, dUy*xyZNA[1], pTwgt); //uyAy 


     dUx   =   TMath::Cos(2.*dPhi);
     dUy   =   TMath::Sin(2.*dPhi);
     Qnx_TPC[1] += dUx*pTwgt;
     Qny_TPC[1] += dUy*pTwgt;

     if(cIndex<6){
       //fullTerm = dUx*(xyZNA[0]*xyZNC[0] - xyZNA[1]*xyZNC[1]) + dUy*(xyZNA[0]*xyZNC[1] + xyZNA[1]*xyZNC[0]);
       fullTerm = dUx*dTermXXYY + dUy*dTermXYXY;

       fHist_v2xV1_ZDN_pTDiff_All[cIndex]->Fill(dPt,fullTerm,pTwgt);
     }


     QRe   =   TMath::Cos(3.*dPhi);
     QIm   =   TMath::Sin(3.*dPhi);
     Qnx_TPC[2] += QRe*pTwgt;
     Qny_TPC[2] += QIm*pTwgt;

     if(cIndex<6){
       Qtot = QRe*(ZARe*ZARe*ZCRe-2.*ZARe*ZAIm*ZCIm-ZCRe*ZAIm*ZAIm) + QIm*(-ZAIm*ZAIm*ZCIm+2.*ZAIm*ZARe*ZCRe+ZCIm*ZARe*ZARe);
       fHist_v3xV1_ZDN_EtaDiff_Comb1[cIndex]->Fill(dEta,Qtot,pTwgt);
       Qtot = QRe*(ZCRe*ZCRe*ZARe-2.*ZCRe*ZCIm*ZAIm-ZARe*ZCIm*ZCIm) + QIm*(-ZCIm*ZCIm*ZAIm+2.*ZCIm*ZCRe*ZARe+ZAIm*ZCRe*ZCRe);
       fHist_v3xV1_ZDN_EtaDiff_Comb2[cIndex]->Fill(dEta,Qtot,pTwgt);
     }


     QRe   =   TMath::Cos(4.*dPhi);
     QIm   =   TMath::Sin(4.*dPhi);
     Qnx_TPC[3] += QRe*pTwgt;
     Qny_TPC[3] += QIm*pTwgt;

     if(cIndex<6){
       Qtot = QRe*((ZARe*ZARe-ZAIm*ZAIm)*(ZCRe*ZCRe-ZCIm*ZCIm)-(ZAIm*ZARe+ZARe*ZAIm)*(ZCIm*ZCRe+ZCRe*ZCIm)) + QIm*((ZAIm*ZARe+ZARe*ZAIm)*(ZCRe*ZCRe-ZCIm*ZCIm)+(ZCIm*ZCRe+ZCRe*ZCIm)*(ZARe*ZARe-ZAIm*ZAIm));
       fHist_v4xV1_ZDN_pTDiff_All[cIndex]->Fill(dPt,Qtot,pTwgt);
     }

     //cout<<"i = "<<i<<" iCentBin = "<<cIndex<<" pTBin "<<ipTBin<<" pt = "<<dPt<<"\t wgt = "<<pTwgt<<endl;
     npoiMult   +=     pTwgt;
   }
  }



  if(npoiMult>0){
   for(int i=0;i<4;i++){
      Qnx_TPC[i] = Qnx_TPC[i]/npoiMult;
      Qny_TPC[i] = Qny_TPC[i]/npoiMult;
    }
  }
  else{
   for(int i=0;i<4;i++){
      Qnx_TPC[i] = 0;
      Qny_TPC[i] = 0;
    }
  }


  //reduced q,Q for 3p Correlator:
  Double_t uPRe=0.,uNRe=0.,uPIm=0.,uNIm=0.,uN2Re=0.,uN2Im=0.,uP2Re=0.,uP2Im=0.,VCRe=0.,VCIm=0., VPRe=0.,VPIm=0.;
  Double_t uPM=0.,uNM=0.;

  if(mSubEtaPos>0 && mSubEtaNeg){
    VCRe = QxEtaNeg[0] /mSubEtaNeg;  //-ve eta for Event plane calculation;
    VCIm = QyEtaNeg[0] /mSubEtaNeg;

    VPRe = QxEtaPos[0] /mSubEtaPos;  //This is for resolution
    VPIm = QyEtaPos[0] /mSubEtaPos;
  }

  if( mNegPoi>0 && mPosPoi>0 ) {
    uPRe = Qxpos[0];
    uPIm = Qypos[0];
    uNRe = Qxneg[0];
    uNIm = Qyneg[0];

    uP2Re = Qxpos[1];
    uP2Im = Qypos[1];
    uN2Re = Qxneg[1];
    uN2Im = Qyneg[1];

    uPM = mPosPoi;
    uNM = mNegPoi;
  }

  //else{uPRe=0.; uPIm=0.; uNRe=0.; uNIm=0.;}

  Double_t TwoQpQnV = 0.;
  Double_t TwoQpQpV = 0.;
  Double_t TwoQnQnV = 0.;


  fullTerm = 0;
  fullReso = 0;

  if(bApplyRecent) {

   if( mNegPoi>0 && mPosPoi>0 ) {
    TwoQpQnV = ((uPRe*uNRe-uPIm*uNIm)*VCRe + (uPRe*uNIm+uPIm*uNRe)*VCIm) / (uPM*uNM) ;
    TwoQpQpV = ((uPRe*uPRe-uPIm*uPIm-uP2Re)*VCRe + (2.*uPRe*uPIm-uP2Im)*VCIm) / (uPM*(uPM-1.)) ;
    TwoQnQnV = ((uNRe*uNRe-uNIm*uNIm-uN2Re)*VCRe + (2.*uNRe*uNIm-uN2Im)*VCIm) / (uNM*(uNM-1.)) ;

    fHist_Corr3p_ZDN_Norm_PN->Fill(EvtCent, TwoQpQnV);
    fHist_Corr3p_ZDN_Norm_PP->Fill(EvtCent, TwoQpQpV);
    fHist_Corr3p_ZDN_Norm_NN->Fill(EvtCent, TwoQnQnV);

    fHist_Corr3p_ZDN_Cent_PN->Fill(EvtCent, TwoQpQnV);
    fHist_Corr3p_ZDN_Cent_PP->Fill(EvtCent, TwoQpQpV);
    fHist_Corr3p_ZDN_Cent_NN->Fill(EvtCent, TwoQnQnV);

    fHist_Reso2EP_TPC_Norm->Fill(EvtCent, (VCRe*VPRe+VCIm*VPIm));
    fHist_Reso2EP_TPC_Cent->Fill(EvtCent, (VCRe*VPRe+VCIm*VPIm));
   }



    if(EvtCent>=5 && EvtCent<=45) {
      fHist_Psi1_ZDCC_wRectCorr->Fill(Psi1C);
      fHist_Psi1_ZDCA_wRectCorr->Fill(Psi1A);

      fHist_PsiSumAC_wiCorr->Fill(PsiSumAC);
      fHist_Psi1_ZDCC_RunByRun->Fill(runindex,Psi1C);
      fHist_Psi1_ZDCA_RunByRun->Fill(runindex,Psi1A);
    }

    fHist_Psi1_ZDCC_wCorrFull->Fill(Psi1C);
    fHist_Psi1_ZDCA_wCorrFull->Fill(Psi1A);

    CentWgt  = fWeight_Cent->GetBinContent(iCentBin);

    fHist_ZDN_resol_Norm_Sep[0]->Fill(EvtCent, xyZNA[0]*xyZNC[0]);
    fHist_ZDN_resol_Cent_Sep[0]->Fill(EvtCent, xyZNA[0]*xyZNC[0]);

    fHist_ZDN_resol_Norm_Sep[1]->Fill(EvtCent, xyZNA[1]*xyZNC[1]);
    fHist_ZDN_resol_Cent_Sep[1]->Fill(EvtCent, xyZNA[1]*xyZNC[1]);

    Double_t v1FillTerms[4] = {Qnx_TPC[0]*xyZNC[0], Qny_TPC[0]*xyZNC[1], Qnx_TPC[0]*xyZNA[0], Qny_TPC[0]*xyZNA[1]};


    for(int i=0; i<4; i++){
      fHist_v2xV1_ZDN_Norm_Sep[i]->Fill(EvtCent, v1FillTerms[i]);
      fHist_v2xV1_ZDN_Cent_Sep[i]->Fill(EvtCent, v1FillTerms[i]);
    }


    //fullTerm = Qnx_TPC[1]*(xyZNA[0]*xyZNC[0] - xyZNA[1]*xyZNC[1]) + Qny_TPC[1]*(xyZNA[0]*xyZNC[1] + xyZNA[1]*xyZNC[0]);
    fullTerm = Qnx_TPC[1]*dTermXXYY + Qny_TPC[1]*dTermXYXY;

    fullReso = xyZNA[0]*xyZNC[0] + xyZNA[1]*xyZNC[1]; //no change in Resolution?

    fHist_v2xV1_ZDN_Norm_All->Fill( EvtCent, fullTerm, CentWgt);
    fHist_v2xV1_ZDN_Cent_All->Fill( EvtCent, fullTerm, CentWgt);
    fHist_v2xV1_ZDN_Refm_All->Fill(nRefMult, fullTerm, CentWgt);

    fHist_ZDN_resol_Norm_All->Fill( EvtCent, fullReso, CentWgt);
    fHist_ZDN_resol_Cent_All->Fill( EvtCent, fullReso, CentWgt);    
    fHist_ZDN_resol_Refm_All->Fill(nRefMult, fullReso, CentWgt);  

    fHist_ZDN_resol_Norm_XX->Fill( EvtCent, xyZNA[0]*xyZNC[0], CentWgt);
    fHist_ZDN_resol_Cent_XX->Fill( EvtCent, xyZNA[0]*xyZNC[0], CentWgt);   
    fHist_ZDN_resol_Norm_YY->Fill( EvtCent, xyZNA[1]*xyZNC[1], CentWgt);
    fHist_ZDN_resol_Cent_YY->Fill( EvtCent, xyZNA[1]*xyZNC[1], CentWgt);   


    QRe  = Qnx_TPC[2];
    QIm  = Qny_TPC[2];
    //integrated v3:
    Qtot = QRe*(ZARe*ZARe*ZCRe-2.*ZARe*ZAIm*ZCIm-ZCRe*ZAIm*ZAIm) + QIm*(-ZAIm*ZAIm*ZCIm+2.*ZAIm*ZARe*ZCRe+ZCIm*ZARe*ZARe);
    fHist_v3xV1_ZDN_Norm_Comb1->Fill( EvtCent, Qtot, CentWgt);
    fHist_v3xV1_ZDN_Cent_Comb1->Fill( EvtCent, Qtot, CentWgt);

    Qtot = QRe*(ZCRe*ZCRe*ZARe-2.*ZCRe*ZCIm*ZAIm-ZARe*ZCIm*ZCIm) + QIm*(-ZCIm*ZCIm*ZAIm+2.*ZCIm*ZCRe*ZARe+ZAIm*ZCRe*ZCRe);
    fHist_v3xV1_ZDN_Norm_Comb2->Fill( EvtCent, Qtot, CentWgt);
    fHist_v3xV1_ZDN_Cent_Comb2->Fill( EvtCent, Qtot, CentWgt);


    QRe  = Qnx_TPC[3];
    QIm  = Qny_TPC[3];
    //integrated v4:
    Qtot = QRe*((ZARe*ZARe-ZAIm*ZAIm)*(ZCRe*ZCRe-ZCIm*ZCIm)-(ZAIm*ZARe+ZARe*ZAIm)*(ZCIm*ZCRe+ZCRe*ZCIm)) 
         + QIm*((ZAIm*ZARe+ZARe*ZAIm)*(ZCRe*ZCRe-ZCIm*ZCIm)+(ZCIm*ZCRe+ZCRe*ZCIm)*(ZARe*ZARe-ZAIm*ZAIm));
    fHist_v4xV1_ZDN_Norm_Comb1->Fill( EvtCent, Qtot, CentWgt);
    fHist_v4xV1_ZDN_Cent_Comb1->Fill( EvtCent, Qtot, CentWgt);

  }




  fHist_Event_count->Fill(stepCount);
  stepCount++;


  fHist_Cent_wiZDCcut->Fill(EvtCent);



  //if(fievent%20==0) {
    //std::cout<<" cent= "<<EvtCent<<"\tV0C[1] = "<<fV0Mult[1]<<"\t V0A[63] = "<<fV0Mult[63]<<"\tRefMult = "<<nRefMult<<std::endl;
    //std::cout<<" <XXYY> = "<<CorrC[0]<<"\t <XYXY> = "<<CorrC[1]<<"\t MagField = "<<fMagField<<std::endl;
  //}


  PostData(1,fListHistos);

  if(!bApplyRecent && fAnalysisSet=="DoGainEq") {
    PostData(2,fListZDCQxy); 
  }
  else{
    PostData(2,fListDummy1);  //default if no 'fAnalysisSet' found.
  }


  fievent++;

  //printf("\n ... ::UserExec() is being called. Step last %d... Event %d  \n",stepCount,fievent);
  
}// UserExec ends




void AliAnalysisTaskZDCGainEq::Terminate(Option_t *)
{
      AliDebug(2,"\n ... AliAnalysisTaskZDCGainEq::Terminate() is being called ...  \n");
}


double AliAnalysisTaskZDCGainEq::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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

 Bool_t AliAnalysisTaskZDCGainEq::plpMV(const AliAODEvent* aod)
 {  // check for multi-vertexer pile-up
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if(!(nPlp=aod->GetNumberOfPileupVerticesTracks()))
  return kFALSE;

  vtPrm = aod->GetPrimaryVertex();
  if(vtPrm == aod->GetPrimaryVertexSPD())
  return kTRUE;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for(int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)aod->GetPileupVertexTracks(ipl);
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




