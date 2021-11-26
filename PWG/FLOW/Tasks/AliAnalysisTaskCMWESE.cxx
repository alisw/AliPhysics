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

//--------------------------------------------------------------------------------
// CME ESE analysis
// Contributor: Wenya Wu, <wenya.wu@cern.ch>, shanghai
//--------------------------------------------------------------------------------

#include <iostream>
#include <cstdlib>
#include <sys/time.h>
// ROOT classes
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TExMap.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TComplex.h"
#include "TSpline.h"
#include "TGrid.h"
// Alice analysis base class
#include "AliAnalysisTaskSE.h"
// Alice analysis additional classes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
// Alice AOD classes
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
// Alice classes
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
// Alice MC classes
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
// Alice "V" classes
#include "AliVParticle.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVVZERO.h"
// Alice PID classes
#include "AliAODPid.h"
#include "AliAODpidUtil.h"
#include "AliPID.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliOADBContainer.h"
#include "AliAnalysisTaskCMWESE.h"
#include "AliCentrality.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskCMWESE);

//---------------------------------------------------

AliAnalysisTaskCMWESE::AliAnalysisTaskCMWESE() : 
  AliAnalysisTaskSE(), 
  fDebug(0),
  fHarmonic(2.),
  fTrigger("kINT7"),
  fFltbit(1),
  fNclsCut(70),
  fChi2Hg(4.0),
  fChi2Lo(0.1),
  fDcaCutz(3.2),
  fDcaCutxy(2.4),
  fPtMin(0.2),
  fPtMax(5.0),
  fCbinHg(8),
  fCbinLo(0),
  fPeriod("LHC10h"),
  fMultComp("pileupByEDSTPC128"),
  fEtaGap(0.3),  
  fV0CalibOn(true),
  fTPCCalibOn(false),
  fQAV0(true),
  fQATPC(false),
  fDoNUE(true),
  fDoNUA(true),
  fDoQnESE(true),
  fDiffV2QA(true),
  fCalcV24VsAch(true),
  fDoTrivalCor(true),
  fCentCut(7.5),
  fRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCentBin(-999),
  fCent(-999),
  fQnBin(-999),
  fAchBin(-999),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fZvtxCut(10.0),
  fListNUE(NULL),
  fListNUA1(NULL),
  fListNUA2(NULL),
  fListNUA3(NULL), 
  fListVZEROCALIB(NULL),
  fNegEtaQ(-999., -999.),
  fNegEtaQStar(-999., -999.),
  fPosEtaQ(-999., -999.),
  fPosEtaQStar(-999., -999.),
  fNegEtaMQ(-999.),
  fPosEtaMQ(-999.), 
  hNUEweightPlus(NULL),
  hNUEweightMinus(NULL), 
  hNUAweightPlus(NULL),
  hNUAweightMinus(NULL),
  hCorrectNUAPos(NULL),
  hCorrectNUANeg(NULL),
  hMultV0Read(NULL),
  hQnPercentile(NULL),
  hQnPercentile_centThisEvt(NULL),
  sp(NULL),
  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),
  fMultCutPU(NULL),
  fOutputList(NULL),
  hEvtCount(NULL),
  hRunNumBin(NULL),
  hPt(NULL),
  hPDedx(NULL),
  hReQn_thisEvt(NULL),
  hImQn_thisEvt(NULL),
  hReQ2n_thisEvt(NULL), 
  hImQ2n_thisEvt(NULL), 
  hReQ2nw2_thisEvt(NULL),
  hImQ2nw2_thisEvt(NULL), 
  hReQnw3_thisEvt(NULL),
  hImQnw3_thisEvt(NULL),
  hMQ_thisEvt(NULL),
  hMQ_w1_thisEvt(NULL),
  hMQ_w2_thisEvt(NULL),
  hMQ_w3_thisEvt(NULL),
  hMQ_w4_thisEvt(NULL),
  hReQPos_thisEvt(NULL),
  hImQPos_thisEvt(NULL),
  hMQPos_thisEvt(NULL),
  hMQPos_w1_thisEvt(NULL),
  hReQNeg_thisEvt(NULL),
  hImQNeg_thisEvt(NULL),
  hMQNeg_thisEvt(NULL),
  hMQNeg_w1_thisEvt(NULL),
  hRepnw1_thisEvt(NULL),
  hImpnw1_thisEvt(NULL),
  hRep2nw2_thisEvt(NULL),
  hImp2nw2_thisEvt(NULL), 
  hRepnw3_thisEvt(NULL),
  hImpnw3_thisEvt(NULL),
  hMp_thisEvt(NULL), 
  hMp_w1_thisEvt(NULL),  
  hMp_w2_thisEvt(NULL),  
  hMp_w3_thisEvt(NULL),  
  hMp_w4_thisEvt(NULL),   
  pRefFlow_thisEvt(NULL),
  pIntd2_thisEvt(NULL),
  pDiffFlowpQStarGap_thisEvt(NULL),
  pDiffFlowpQStarGapPos_thisEvt(NULL),
  pDiffFlowpQStarGapNeg_thisEvt(NULL),
  pDiffFlowpQStarDirect_thisEvt(NULL),
  pDiffFlow4Direct_thisEvt(NULL)
{
  for (int i = 0; i < NRUNNUM; ++i) fRunNumList[i]="0";

  for (int i = 0; i < 3; ++i)   pV0XMeanRead[i]=NULL; 
  for (int i = 0; i < 3; ++i)   pV0YMeanRead[i]=NULL;

  for (int i = 0; i < NRUNNUM; ++i){
    hMultV0[i]=NULL; //Dobrin
    for (int j = 0; j < 2; ++j)   hQxnmV0[i][j]=NULL;
    for (int j = 0; j < 2; ++j)   hQynmV0[i][j]=NULL;
  }
  for (int i = 0; i < 64; ++i) fMultV0Ch[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0XMean[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0YMean[i]=-999.;

  for (int i = 0; i < 80; ++i) splQ2c[i]=NULL; 

  for (int i = 0; i < 2; ++i) hCent[i]=NULL;
  for (int i = 0; i < 2; ++i) hVz[i]=NULL;

  for (int i = 0; i < 8; ++i) hCentQA[i]=NULL;
  for (int i = 0; i < 2; ++i) hMultCentQA[i]=NULL;
  for (int i = 0; i < 6; ++i) hMultMultQA[i]=NULL;

  for (int i = 0; i < 2; ++i) hEta[i]=NULL;
  for (int i = 0; i < 2; ++i) hPhi[i]=NULL;
  for (int i = 0; i < 2; ++i) hEtaPhi[i]=NULL;
  for (int i = 0; i < 2; ++i) hDcaXy[i]=NULL;
  for (int i = 0; i < 2; ++i) hDcaZ[i]=NULL;
  for (int i = 0; i < 2; ++i) hNhits[i]=NULL;

  for (int i = 0; i < 3; ++i) hQxCentRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQxVtxRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQyCentRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQyVtxRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQnCentRecenter[i]=NULL;
  
  for (int i = 0; i < NCENTBINS; ++i){
    for (int j = 0; j < 3; ++j) hPsiV0Recenter[i][j]=NULL;
  } 

  for (int i = 0; i < NCENTBINS; ++i) pRefFlow[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pIntd2[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pIntd2Ach[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pAch[i]=NULL;

  for (int i = 0; i < NCENTBINS+1; ++i) {
    for (int j = 0; j < NQNBINS; ++j) hMult[i][j]=NULL;
  }
  for (int i = 0; i < NCENTBINS+1; ++i) {
    for (int j = 0; j < NQNBINS; ++j) hAch[i][j]=NULL;
  }

  for (int i = 0; i < NCENTBINS; ++i) pDiffFlowpQStarGap[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) hDiffFlowYield[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pDiffFlowpQStarGapPos[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pDiffFlowpQStarGapNeg[i]=NULL;
  for (int i = 0; i < NCENTBINS+1; ++i){
    for (int j = 0; j < NACHBINS; ++j) pRefFlowAch[i][j]=NULL;
    for (int j = 0; j < NACHBINS; ++j) pDiffFlowpQStarAch[i][j]=NULL;
    for (int j = 0; j < NACHBINS; ++j) pDiffFlow4Ach[i][j]=NULL;
    for (int j = 0; j < NACHBINS; ++j) hDiffFlowYieldAch[i][j]=NULL;
  }

}

//---------------------------------------------------

AliAnalysisTaskCMWESE::AliAnalysisTaskCMWESE(const char *name, TString _PR, bool _NUE, bool _NUA, bool _V0Calib) : 
  AliAnalysisTaskSE(name), 
  fDebug(0),
  fHarmonic(2.),
  fTrigger("kINT7"),
  fFltbit(1),
  fNclsCut(70),
  fChi2Hg(4.0),
  fChi2Lo(0.1),
  fDcaCutz(3.2),
  fDcaCutxy(2.4),
  fPtMin(0.2),
  fPtMax(5.0),
  fCbinHg(8),
  fCbinLo(0),
  fPeriod(_PR),
  fMultComp("pileupByEDSTPC128"),
  fEtaGap(0.3),  
  fV0CalibOn(_V0Calib),
  fTPCCalibOn(false),
  fQAV0(true),
  fQATPC(false),
  fDoNUE(_NUE),
  fDoNUA(_NUA),
  fDoQnESE(true),
  fDiffV2QA(true),
  fCalcV24VsAch(true),
  fDoTrivalCor(true),
  fCentCut(7.5),
  fRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCentBin(-999),
  fCent(-999),
  fQnBin(-999),
  fAchBin(-999),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fZvtxCut(10.0),
  fListNUE(NULL),
  fListNUA1(NULL),
  fListNUA2(NULL),
  fListNUA3(NULL), 
  fListVZEROCALIB(NULL),
  fNegEtaQ(-999., -999.),
  fNegEtaQStar(-999., -999.),
  fPosEtaQ(-999., -999.),
  fPosEtaQStar(-999., -999.),
  fNegEtaMQ(-999.),
  fPosEtaMQ(-999.), 
  hNUEweightPlus(NULL),
  hNUEweightMinus(NULL), 
  hNUAweightPlus(NULL),
  hNUAweightMinus(NULL),
  hCorrectNUAPos(NULL),
  hCorrectNUANeg(NULL),
  hMultV0Read(NULL),
  hQnPercentile(NULL),
  hQnPercentile_centThisEvt(NULL),
  sp(NULL),
  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),
  fMultCutPU(NULL),
  fOutputList(NULL),
  hEvtCount(NULL),
  hRunNumBin(NULL),
  hPt(NULL),
  hPDedx(NULL),
  hReQn_thisEvt(NULL),
  hImQn_thisEvt(NULL),
  hReQ2n_thisEvt(NULL), 
  hImQ2n_thisEvt(NULL), 
  hReQ2nw2_thisEvt(NULL),
  hImQ2nw2_thisEvt(NULL), 
  hReQnw3_thisEvt(NULL),
  hImQnw3_thisEvt(NULL),
  hMQ_thisEvt(NULL),
  hMQ_w1_thisEvt(NULL),
  hMQ_w2_thisEvt(NULL),
  hMQ_w3_thisEvt(NULL),
  hMQ_w4_thisEvt(NULL),
  hReQPos_thisEvt(NULL),
  hImQPos_thisEvt(NULL),
  hMQPos_thisEvt(NULL),
  hMQPos_w1_thisEvt(NULL),
  hReQNeg_thisEvt(NULL),
  hImQNeg_thisEvt(NULL),
  hMQNeg_thisEvt(NULL),
  hMQNeg_w1_thisEvt(NULL),
  hRepnw1_thisEvt(NULL),
  hImpnw1_thisEvt(NULL),
  hRep2nw2_thisEvt(NULL),
  hImp2nw2_thisEvt(NULL), 
  hRepnw3_thisEvt(NULL),
  hImpnw3_thisEvt(NULL),
  hMp_thisEvt(NULL), 
  hMp_w1_thisEvt(NULL),  
  hMp_w2_thisEvt(NULL),  
  hMp_w3_thisEvt(NULL),  
  hMp_w4_thisEvt(NULL),   
  pRefFlow_thisEvt(NULL),
  pIntd2_thisEvt(NULL),
  pDiffFlowpQStarGap_thisEvt(NULL),
  pDiffFlowpQStarGapPos_thisEvt(NULL),
  pDiffFlowpQStarGapNeg_thisEvt(NULL),
  pDiffFlowpQStarDirect_thisEvt(NULL),
  pDiffFlow4Direct_thisEvt(NULL)
{
  for (int i = 0; i < NRUNNUM; ++i) fRunNumList[i]="0";

  for (int i = 0; i < 3; ++i)   pV0XMeanRead[i]=NULL; 
  for (int i = 0; i < 3; ++i)   pV0YMeanRead[i]=NULL;

  for (int i = 0; i < NRUNNUM; ++i){
    hMultV0[i]=NULL; //Dobrin
    for (int j = 0; j < 2; ++j)   hQxnmV0[i][j]=NULL;
    for (int j = 0; j < 2; ++j)   hQynmV0[i][j]=NULL;
  }
  for (int i = 0; i < 64; ++i) fMultV0Ch[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0XMean[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0YMean[i]=-999.;

  for (int i = 0; i < 80; ++i) splQ2c[i]=NULL; 

  for (int i = 0; i < 2; ++i) hCent[i]=NULL;
  for (int i = 0; i < 2; ++i) hVz[i]=NULL;

  for (int i = 0; i < 8; ++i) hCentQA[i]=NULL;
  for (int i = 0; i < 2; ++i) hMultCentQA[i]=NULL;
  for (int i = 0; i < 6; ++i) hMultMultQA[i]=NULL;

  for (int i = 0; i < 2; ++i) hEta[i]=NULL;
  for (int i = 0; i < 2; ++i) hPhi[i]=NULL;
  for (int i = 0; i < 2; ++i) hEtaPhi[i]=NULL;
  for (int i = 0; i < 2; ++i) hDcaXy[i]=NULL;
  for (int i = 0; i < 2; ++i) hDcaZ[i]=NULL;
  for (int i = 0; i < 2; ++i) hNhits[i]=NULL;

  for (int i = 0; i < 3; ++i) hQxCentRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQxVtxRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQyCentRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQyVtxRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQnCentRecenter[i]=NULL;
  
  for (int i = 0; i < NCENTBINS; ++i){
    for (int j = 0; j < 3; ++j) hPsiV0Recenter[i][j]=NULL;
  } 

  for (int i = 0; i < NCENTBINS; ++i) pRefFlow[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pIntd2[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pIntd2Ach[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pAch[i]=NULL;

  for (int i = 0; i < NCENTBINS+1; ++i) {
    for (int j = 0; j < NQNBINS; ++j) hMult[i][j]=NULL;
  }
  for (int i = 0; i < NCENTBINS+1; ++i) {
    for (int j = 0; j < NQNBINS; ++j) hAch[i][j]=NULL;
  }

  for (int i = 0; i < NCENTBINS; ++i) pDiffFlowpQStarGap[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) hDiffFlowYield[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pDiffFlowpQStarGapPos[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pDiffFlowpQStarGapNeg[i]=NULL;
  for (int i = 0; i < NCENTBINS+1; ++i){
    for (int j = 0; j < NACHBINS; ++j) pRefFlowAch[i][j]=NULL;
    for (int j = 0; j < NACHBINS; ++j) pDiffFlowpQStarAch[i][j]=NULL;
    for (int j = 0; j < NACHBINS; ++j) pDiffFlow4Ach[i][j]=NULL;
    for (int j = 0; j < NACHBINS; ++j) hDiffFlowYieldAch[i][j]=NULL;
  }

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());

  int inputslot = 1;
  if (fPeriod.EqualTo("LHC10h")){
    if (fDoNUE) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    }
    if (fDoNUA) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    }
    if (fV0CalibOn) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    } 
  }; 
  if (fPeriod.EqualTo("LHC11h")){
    if (fDoNUE) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    }
    if (fDoNUA) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    }
    if (fV0CalibOn) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    } 
  };
  if (fPeriod.EqualTo("LHC15o")){
    if (fDoNUE) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    }
    if (fDoNUA) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    }
    if (fV0CalibOn) {
      DefineInput(inputslot, TList::Class());
      inputslot++;
    } 
  };

}

//------------------------------------------------

AliAnalysisTaskCMWESE::AliAnalysisTaskCMWESE(const char *name) : 
  AliAnalysisTaskSE(name), 
  fDebug(0),
  fHarmonic(2.),
  fTrigger("kINT7"),
  fFltbit(1),
  fNclsCut(70),
  fChi2Hg(4.0),
  fChi2Lo(0.1),
  fDcaCutz(3.2),
  fDcaCutxy(2.4),
  fPtMin(0.2),
  fPtMax(5.0),
  fCbinHg(8),
  fCbinLo(0),
  fPeriod("LHC10h"),
  fMultComp("pileupByEDSTPC128"),
  fEtaGap(0.3),  
  fV0CalibOn(true),
  fTPCCalibOn(false),
  fQAV0(true),
  fQATPC(false),
  fDoNUE(true),
  fDoNUA(true),
  fDoQnESE(true),
  fDiffV2QA(true),
  fCalcV24VsAch(true),
  fDoTrivalCor(true),
  fCentCut(7.5),
  fRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCentBin(-999),
  fCent(-999),
  fQnBin(-999),
  fAchBin(-999),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fZvtxCut(10.0),
  fListNUE(NULL),
  fListNUA1(NULL),
  fListNUA2(NULL),
  fListNUA3(NULL), 
  fListVZEROCALIB(NULL),
  fNegEtaQ(-999., -999.),
  fNegEtaQStar(-999., -999.),
  fPosEtaQ(-999., -999.),
  fPosEtaQStar(-999., -999.),
  fNegEtaMQ(-999.),
  fPosEtaMQ(-999.), 
  hNUEweightPlus(NULL),
  hNUEweightMinus(NULL), 
  hNUAweightPlus(NULL),
  hNUAweightMinus(NULL),
  hCorrectNUAPos(NULL),
  hCorrectNUANeg(NULL),
  hMultV0Read(NULL),
  hQnPercentile(NULL),
  hQnPercentile_centThisEvt(NULL),
  sp(NULL),
  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),
  fMultCutPU(NULL),
  fOutputList(NULL),
  hEvtCount(NULL),
  hRunNumBin(NULL),
  hPt(NULL),
  hPDedx(NULL),
  hReQn_thisEvt(NULL),
  hImQn_thisEvt(NULL),
  hReQ2n_thisEvt(NULL), 
  hImQ2n_thisEvt(NULL), 
  hReQ2nw2_thisEvt(NULL),
  hImQ2nw2_thisEvt(NULL), 
  hReQnw3_thisEvt(NULL),
  hImQnw3_thisEvt(NULL),
  hMQ_thisEvt(NULL),
  hMQ_w1_thisEvt(NULL),
  hMQ_w2_thisEvt(NULL),
  hMQ_w3_thisEvt(NULL),
  hMQ_w4_thisEvt(NULL),
  hReQPos_thisEvt(NULL),
  hImQPos_thisEvt(NULL),
  hMQPos_thisEvt(NULL),
  hMQPos_w1_thisEvt(NULL),
  hReQNeg_thisEvt(NULL),
  hImQNeg_thisEvt(NULL),
  hMQNeg_thisEvt(NULL),
  hMQNeg_w1_thisEvt(NULL),
  hRepnw1_thisEvt(NULL),
  hImpnw1_thisEvt(NULL),
  hRep2nw2_thisEvt(NULL),
  hImp2nw2_thisEvt(NULL), 
  hRepnw3_thisEvt(NULL),
  hImpnw3_thisEvt(NULL),
  hMp_thisEvt(NULL), 
  hMp_w1_thisEvt(NULL),  
  hMp_w2_thisEvt(NULL),  
  hMp_w3_thisEvt(NULL),  
  hMp_w4_thisEvt(NULL),   
  pRefFlow_thisEvt(NULL),
  pIntd2_thisEvt(NULL),
  pDiffFlowpQStarGap_thisEvt(NULL),
  pDiffFlowpQStarGapPos_thisEvt(NULL),
  pDiffFlowpQStarGapNeg_thisEvt(NULL),
  pDiffFlowpQStarDirect_thisEvt(NULL),
  pDiffFlow4Direct_thisEvt(NULL)
{
  for (int i = 0; i < NRUNNUM; ++i) fRunNumList[i]="0";

  for (int i = 0; i < 3; ++i)   pV0XMeanRead[i]=NULL; 
  for (int i = 0; i < 3; ++i)   pV0YMeanRead[i]=NULL;

  for (int i = 0; i < NRUNNUM; ++i){
    hMultV0[i]=NULL; //Dobrin
    for (int j = 0; j < 2; ++j)   hQxnmV0[i][j]=NULL;
    for (int j = 0; j < 2; ++j)   hQynmV0[i][j]=NULL;
  }
  for (int i = 0; i < 64; ++i) fMultV0Ch[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0XMean[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0YMean[i]=-999.;

  for (int i = 0; i < 80; ++i) splQ2c[i]=NULL; 

  for (int i = 0; i < 2; ++i) hCent[i]=NULL;
  for (int i = 0; i < 2; ++i) hVz[i]=NULL;

  for (int i = 0; i < 8; ++i) hCentQA[i]=NULL;
  for (int i = 0; i < 2; ++i) hMultCentQA[i]=NULL;
  for (int i = 0; i < 6; ++i) hMultMultQA[i]=NULL;

  for (int i = 0; i < 2; ++i) hEta[i]=NULL;
  for (int i = 0; i < 2; ++i) hPhi[i]=NULL;
  for (int i = 0; i < 2; ++i) hEtaPhi[i]=NULL;
  for (int i = 0; i < 2; ++i) hDcaXy[i]=NULL;
  for (int i = 0; i < 2; ++i) hDcaZ[i]=NULL;
  for (int i = 0; i < 2; ++i) hNhits[i]=NULL;

  for (int i = 0; i < 3; ++i) hQxCentRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQxVtxRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQyCentRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQyVtxRecenter[i]=NULL;
  for (int i = 0; i < 3; ++i) hQnCentRecenter[i]=NULL;
  
  for (int i = 0; i < NCENTBINS; ++i){
    for (int j = 0; j < 3; ++j) hPsiV0Recenter[i][j]=NULL;
  } 

  for (int i = 0; i < NCENTBINS; ++i) pRefFlow[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pIntd2[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pIntd2Ach[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pAch[i]=NULL;

  for (int i = 0; i < NCENTBINS+1; ++i) {
    for (int j = 0; j < NQNBINS; ++j) hMult[i][j]=NULL;
  }
  for (int i = 0; i < NCENTBINS+1; ++i) {
    for (int j = 0; j < NQNBINS; ++j) hAch[i][j]=NULL;
  }

  for (int i = 0; i < NCENTBINS; ++i) pDiffFlowpQStarGap[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) hDiffFlowYield[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pDiffFlowpQStarGapPos[i]=NULL;
  for (int i = 0; i < NCENTBINS; ++i) pDiffFlowpQStarGapNeg[i]=NULL;
  for (int i = 0; i < NCENTBINS+1; ++i){
    for (int j = 0; j < NACHBINS; ++j) pRefFlowAch[i][j]=NULL;
    for (int j = 0; j < NACHBINS; ++j) pDiffFlowpQStarAch[i][j]=NULL;
    for (int j = 0; j < NACHBINS; ++j) pDiffFlow4Ach[i][j]=NULL;
    for (int j = 0; j < NACHBINS; ++j) hDiffFlowYieldAch[i][j]=NULL;
  }

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());

}
  //------------------------------------------------

AliAnalysisTaskCMWESE::~AliAnalysisTaskCMWESE()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  if (fOutputList) delete fOutputList;
}

//_____________________________________________________________________________

void AliAnalysisTaskCMWESE::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate()");
}

//---------------------------------------------------

void AliAnalysisTaskCMWESE::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetName(GetName());
  fOutputList->SetOwner(kTRUE);
  if (fDebug) Printf("Start Create Outputs");
  //------------------
  // Info
  // ------------------
  if (fPeriod.EqualTo("LHC10h") ){
      TString RunList[90]={
        "139510","139507","139505","139503","139465","139438","139437","139360","139329","139328","139314","139310",
        "139309","139173","139107","139105","139038","139037","139036","139029","139028","138872","138871","138870",
        "138837","138732","138730","138666","138662","138653","138652","138638","138624","138621","138583","138582",
        "138579","138578","138534","138469","138442","138439","138438","138396","138364","138275","138225","138201",
        "138197","138192","138190","137848","137844","137752","137751","137724","137722","137718","137704","137693",
        "137692","137691","137686","137685","137639","137638","137608","137595","137549","137546","137544","137541",
        "137539","137531","137530","137443","137441","137440","137439","137434","137432","137431","137430","137243",
        "137236","137235","137232","137231","137162","137161"};
      hRunNumBin = new TH1I("runNumBin","",100,0,100);
      for (int i=0; i<90; ++i){    
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,RunList[i].Data());
        fRunNumList[i] = RunList[i];
      }
      fOutputList->Add(hRunNumBin);
  } else if (fPeriod.EqualTo("LHC11h") ){
    TString RunList[68]={"167915", "168115", "168460", "169035", "169238", "169859", "170228 ", "167920", "168310", "168464", "169091", "169411", 
      "169923", "170230", "167985", "168311", "168467", "169094", "169415", "170027", "170268", "167987", "168322", "168511", "169138", "169417", 
      "170081", "170269", "167988", "168325", "168512", "169144", "169835", "170155", "170270", "168069", "168341", "168514", "169145", "169837", 
      "170159", "170306", "168076", "168342", "168777", "169148", "169838", "170163", "170308", "168105", "168361", "168826", "169156", "169846", 
      "170193", "170309", "168107", "168362", "168988", "169160", "169855", "170203", "168108 ", "168458", "168992", "169167", "169858", "170204"};
      hRunNumBin = new TH1I("runNumBin","",100,0,100);
      for (int i=0; i<68; ++i){    
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,RunList[i].Data());
        fRunNumList[i] = RunList[i];
      }
      fOutputList->Add(hRunNumBin);
  } else if (fPeriod.EqualTo("LHC15o") ){
      TString RunList[138]={
       "246994","246991","246989","246984","246982","246948","246945","246928","246871","246870","246867","246865", 
       "246864","246859","246858","246851","246847","246846","246845","246844","246810","246809","246808","246807", 
       "246805","246804","246766","246765","246763","246760","246759","246758","246757","246751","246750","246434", 
       "246431","246424","246392","246391","246276","246275","246272","246271","246225","246222","246217","246185", 
       "246182","246181","246180","246178","246153","246152","246151","246148","246115","246113","246089","246087", 
       "246053","246052","246049","246048","246042","246037","246036","246012","246003","246001","245963","245954", 
       "245952","245949","245923","245833","245831","245829","245793","245785","245775","245766","245759","245752", 
       "245731","245729","245705","245702","245692","245683","245554","245545","245544","245543","245542","245540", 
       "245535","245507","245505","245504","245501","245497","245496","245454","245453","245450","245446","245441", 
       "245411","245410","245409","245407","245401","245397","245396","245353","245349","245347","245346","245345", 
       "245343","245259","245233","245232","245231","245152","245151","245146","245145","245068","245066","245064", 
       "244983","244982","244980","244975","244918","244917"};
      hRunNumBin = new TH1I("runNumBin","",150,0,150);
      for (int i=0; i<138; ++i){    
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,RunList[i].Data());
        fRunNumList[i] = RunList[i];
      }
      fOutputList->Add(hRunNumBin);
  }
  if (fDebug) Printf("runNum created!");

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC11h")){

      if (!gGrid) {
         TGrid::Connect("alien://");
      }

      TFile* fileV0Calib = NULL;
      if (fPeriod.EqualTo("LHC15o")) fileV0Calib = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc15o/calibV015oP2.root");
      if (fPeriod.EqualTo("LHC11h")) fileV0Calib = TFile::Open("alien:///alice/cern.ch/user/w/wenya/refData/reflhc11h/calibV011h_2P2.root");

      if(!fileV0Calib){
          printf("OADB V0 calibration file cannot be opened\n");
           return;
      }

      // Mult
      AliOADBContainer* contMult = (AliOADBContainer*) fileV0Calib->Get("hMultV0BefCorPfpx");
      if(!contMult){
          printf("OADB object hMultV0BefCorr is not available in the file\n");
          return;
      }
  
      // V0C Qx Mean 
      AliOADBContainer* contQxncm = (AliOADBContainer*) fileV0Calib->Get(Form("fqxc%im",(int)fHarmonic));
      if(!contQxncm){
          printf("OADB object fqxc2m is not available in the file\n");
          return;
      }
      
      // V0C Qy Mean 
      AliOADBContainer* contQyncm = (AliOADBContainer*) fileV0Calib->Get(Form("fqyc%im",(int)fHarmonic));
      if(!contQyncm){
          printf("OADB object fqyc2m is not available in the file\n");
          return;
      }
      
      // V0A Qx Mean 
      AliOADBContainer* contQxnam = (AliOADBContainer*) fileV0Calib->Get(Form("fqxa%im",(int)fHarmonic));
      if(!contQxnam){
          printf("OADB object fqxa2m is not available in the file\n");
          return;
      }
      
      // V0A Qy Mean 
      AliOADBContainer* contQynam = (AliOADBContainer*) fileV0Calib->Get(Form("fqya%im",(int)fHarmonic));
      if(!contQynam){
          printf("OADB object fqya2m is not available in the file\n");
          return;
      }

      int nRunnumThisPeriod = -1;
      if (fPeriod.EqualTo("LHC10h")) nRunnumThisPeriod = 90;
      if (fPeriod.EqualTo("LHC11h")) nRunnumThisPeriod = 68;
      if (fPeriod.EqualTo("LHC15o")) nRunnumThisPeriod = 138;
      for (int iRun = 0; iRun < nRunnumThisPeriod; ++iRun){

        if(!(contMult->GetObject(fRunNumList[iRun].Atoi() ) ) ){
            printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", fRunNumList[iRun].Atoi() );
            return;
        }
        hMultV0[iRun] = ((TH1D*) contMult->GetObject(fRunNumList[iRun].Atoi()));

         // V0C Qx Mean 
         if(!(contQxncm->GetObject(fRunNumList[iRun].Atoi() ) ) ){
             printf("OADB objectForm fqxcnm ,fHarmonic)is not available for run %i\n", fRunNumList[iRun].Atoi() );
             return;
         }
         hQxnmV0[iRun][1] = ((TH1D*) contQxncm->GetObject(fRunNumList[iRun].Atoi() ) );
         
         // V0C Qy Mean 
         if(!(contQyncm->GetObject(fRunNumList[iRun].Atoi() ) ) ){
             printf("OADB object fqycnm is not available for run %i\n", fRunNumList[iRun].Atoi() );
             return;
         }
         hQynmV0[iRun][1] = ((TH1D*) contQyncm->GetObject(fRunNumList[iRun].Atoi() ) );
         
         // V0A Qx Mean
         if(!(contQxnam->GetObject(fRunNumList[iRun].Atoi() ) ) ){
             printf("OADB object fqxanm is not available for run %i\n", fRunNumList[iRun].Atoi() );
             return;
         }
         hQxnmV0[iRun][2] = ((TH1D*) contQxnam->GetObject(fRunNumList[iRun].Atoi() ) );
         
         // V0A Qy Mean
         if(!(contQynam->GetObject(fRunNumList[iRun].Atoi() ) ) ){
             printf("OADB object fqyanm is not available for run %i\n", fRunNumList[iRun].Atoi() );
             return;
         }
         hQynmV0[iRun][2] = ((TH1D*) contQynam->GetObject(fRunNumList[iRun].Atoi() ) );
      }
  };

  // Copy TList
  Int_t inSlotCounter=1;
  if(fDoNUE) {
    if (fPeriod.EqualTo("LHC11h") || fPeriod.EqualTo("LHC10h") ) {
      fListNUE = (TList*) GetInputData(inSlotCounter);
      inSlotCounter++;
      if (!fListNUE) {
        AliError(Form("%s: Wrong NUE file Run1.", GetName()));
        return;
      };
    }
    else if (fPeriod.EqualTo("LHC15o")){
      fListNUE = (TList*) GetInputData(inSlotCounter);
      inSlotCounter++;
      if (!fListNUE) {
        AliError(Form("%s: Wrong NUE file 15o.", GetName()));
        return;
      };     
    }
  }

  if(fDoNUA) {
    if (fPeriod.EqualTo("LHC10h") ) {
      fListNUA1 = (TList*) GetInputData(inSlotCounter);
      inSlotCounter++;
      if (!fListNUA1) {
        AliError(Form("%s: Wrong NUA file 10h.", GetName()));
        return;
      };
    } else if (fPeriod.EqualTo("LHC11h") ){
      fListNUA1 = (TList*) GetInputData(inSlotCounter);
      inSlotCounter++;
      if (!fListNUA1) {        
        AliError(Form("%s: Wrong NUA file 11h.", GetName()));
        return;
      }
    } else if (fPeriod.EqualTo("LHC15o") ){ //Rihan's NUA
      fListNUA1 = (TList*) GetInputData(inSlotCounter);
      inSlotCounter++;
      if (!fListNUA1) {        
        AliError(Form("%s: Wrong NUA file 15o2.", GetName()));
        return;
      }
    } 
  }

  // Qn Calib
  if (fV0CalibOn){
      if (fPeriod.EqualTo("LHC10h") ) {
            fListVZEROCALIB = (TList*) GetInputData(inSlotCounter);
            inSlotCounter++;
            if (!fListVZEROCALIB) {
             AliError(Form("%s: Wrong V0 Calib file 10h.", GetName()));
              return;
            };      
            // Read GE Histo (x_y) : (iCh_runnumBin)
            hMultV0Read = (TH2D*)fListVZEROCALIB->FindObject("hMultV0");

            // Read Qx/y Mean Histo (x_y_z) : (runnumBin_centBin_vzBin)
            pV0XMeanRead[1] = (TProfile3D*)fListVZEROCALIB->FindObject("pV0CCosMean");
            pV0YMeanRead[1] = (TProfile3D*)fListVZEROCALIB->FindObject("pV0CSinMean"); 
            pV0XMeanRead[2] = (TProfile3D*)fListVZEROCALIB->FindObject("pV0ACosMean");
            pV0YMeanRead[2] = (TProfile3D*)fListVZEROCALIB->FindObject("pV0ASinMean");

            // Read Percentail histo (x_y) : (intCent_percentile)
            hQnPercentile = (TH2D*)fListVZEROCALIB->FindObject("h_qnPercentile");

      } else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC11h")) {
            fListVZEROCALIB = (TList*) GetInputData(inSlotCounter);
            inSlotCounter++;
            if (!fListVZEROCALIB) {
             AliError(Form("%s: Spline list for 15o cannot be opened.", GetName()));
              return;
            };
            //TSplines in 1% centralicities
            for (Int_t isp = 0; isp < 80; isp++) splQ2c[isp] = (TSpline3*)fListVZEROCALIB->FindObject(Form("sp_q2V0C_%d", isp));
      }
  }

  //------------------
  // QA
  //------------------
  // event-wise
  hEvtCount = new TH1D("hEvtCount", "Event Count", 25, 1, 26);
  hEvtCount->GetXaxis()->SetBinLabel(1,"all");
  hEvtCount->GetXaxis()->SetBinLabel(2,"readin");
  hEvtCount->GetXaxis()->SetBinLabel(3,"evt");
  hEvtCount->GetXaxis()->SetBinLabel(4,"runNum");
  hEvtCount->GetXaxis()->SetBinLabel(5,"vertex");
  hEvtCount->GetXaxis()->SetBinLabel(6,"centrality");
  hEvtCount->GetXaxis()->SetBinLabel(7,"pileup");
  hEvtCount->GetXaxis()->SetBinLabel(8,"analysis");
  hEvtCount->GetXaxis()->SetBinLabel(9,"loop1");
  hEvtCount->GetXaxis()->SetBinLabel(10,"cumulant");
  hEvtCount->GetXaxis()->SetBinLabel(11,"Ach");
  hEvtCount->GetXaxis()->SetBinLabel(12,"q-vector");
  hEvtCount->GetXaxis()->SetBinLabel(13,"loop2");
  hEvtCount->GetXaxis()->SetBinLabel(14,"fill");
  hEvtCount->GetXaxis()->SetBinLabel(20,"manager");
  hEvtCount->GetXaxis()->SetBinLabel(21,"handler");
  hEvtCount->GetXaxis()->SetBinLabel(22,"fAOD");
  hEvtCount->GetXaxis()->SetBinLabel(23,"fPID");
  hEvtCount->GetXaxis()->SetBinLabel(24,"fUtils");
  hEvtCount->GetXaxis()->SetBinLabel(25,"fMultSel");
  fOutputList->Add(hEvtCount);

  hCent[0] = new TH1D("hCentBefCut", "", 100, 0., 100.);
  hCent[1] = new TH1D("hCentAftCut", "", 100, 0., 100.);
  fOutputList->Add(hCent[0]);
  fOutputList->Add(hCent[1]);

  hVz[0] = new TH1D("hVzBeforeCut", "", 200, -50., 50.);
  hVz[1] = new TH1D("hVzAfterCut",  "", 200, -50., 50.);
  fOutputList->Add(hVz[0]);
  fOutputList->Add(hVz[1]);

  hCentQA[0] = new TH2D("hCentQAV0MSPD1BefCut", ";centV0M;centSPD1", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[0]);
  hCentQA[1] = new TH2D("hCentQAV0MSPD1AftCut", ";centV0M;centSPD1", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[1]);

  hCentQA[2] = new TH2D("hCentQAV0MTRKBefCut", ";centV0M;centTRK", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[2]);
  hCentQA[3] = new TH2D("hCentQAV0MTRKAftCut", ";centV0M;centTRK", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[3]);

  hCentQA[4] = new TH2D("hCentQAV0MSPD0BefCut", ";centV0M;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[4]);
  hCentQA[5] = new TH2D("hCentQAV0MSPD0AftCut", ";centV0M;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[5]);

  hCentQA[6] = new TH2D("hCentQASPD1SPD0BefCut", ";centSPD1;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[6]);
  hCentQA[7] = new TH2D("hCentQASPD1SPD0AftCut", ";centSPD1;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(hCentQA[7]);

  hMultCentQA[0] = new TH2D("hMultCentQABefCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
  hMultCentQA[1] = new TH2D("hMultCentQAAftCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
  fOutputList->Add(hMultCentQA[0]);
  fOutputList->Add(hMultCentQA[1]);

  if (fMultComp.EqualTo("pileupByEDSTPC128") || fMultComp.EqualTo("pileupByGlobalTPC1")){
    hMultMultQA[0] = new TH2D("hMultMultQAmTPCmESDPileupBefCut", "befCut;multTPC;multESD", 50, 0, 5000, 160, 0, 16000);
    hMultMultQA[1] = new TH2D("hMultMultQAmClobalmTPCEFPileupBefCut", "befCut;multGlobal;multTPCFE", 50, 0, 5000, 50, 0, 5000);
    hMultMultQA[2] = new TH2D("hMultMultQAmFB32mTrkTOFBefCut", "befCut;multTrkTOF;nTrk", 201, 0, 20000, 201, 0, 20000);
    hMultMultQA[3] = new TH2D("hMultMultQAmTPCmESDPileupAftCut", "aftCut;multTPC;multESD", 50, 0, 5000, 160, 0, 16000);
    hMultMultQA[4] = new TH2D("hMultMultQAmClobalmTPCEFPileupAftCut", "aftCut;multGlobal;multTPCFE", 50, 0, 5000, 50, 0, 5000);
    hMultMultQA[5] = new TH2D("hMultMultQAmFB32mTrkTOFAftCut", "aftCut;multTrkTOF;nTrk", 201, 0, 20000, 201, 0, 20000);
    fOutputList->Add(hMultMultQA[0]);
    fOutputList->Add(hMultMultQA[1]);
    fOutputList->Add(hMultMultQA[2]);
    fOutputList->Add(hMultMultQA[3]);
    fOutputList->Add(hMultMultQA[4]);
    fOutputList->Add(hMultMultQA[5]);
  }

  if (fDoQnESE){
    for (int i = 0; i < NCENTBINS+1; ++i) {
      for (int j = 0; j < NQNBINS; ++j){
        if (i==0) {
          hMult[i][j] = new TH1D(Form("hMultMbQn%i", j ), "", 100, 0., 5000.);
          hAch[i][j] = new TH1D(Form("hAchMbQn%i", j ), "", 200, -1., 1.);             
        }
        else {
          hMult[i][j] = new TH1D(Form("hMultCent%iQn%i", i-1, j ), "", 100, 0., 5000.);
          hAch[i][j] = new TH1D(Form("hAchCent%iQn%i", i-1, j ), "", 200, -1., 1.);      
        }
      }
    }
    for (int i = fCbinLo; i < fCbinHg+1; ++i) {
      for (int j = 0; j < NQNBINS; ++j) fOutputList->Add(hMult[i][j]); 
      for (int j = 0; j < NQNBINS; ++j) fOutputList->Add(hAch[i][j]); 
    }
  }

  // track-wise
  hPt = new TH1D("hPt", "", 200, 0., 20.);
  fOutputList->Add(hPt);
  hEta[0] = new TH1D("hEtaBeforeCut", "", 200, -10., 10.);
  hEta[1] = new TH1D("hEtaAfterCut",  "", 200, -10., 10.);
  fOutputList->Add(hEta[0]);
  fOutputList->Add(hEta[1]);
  hPhi[0] = new TH1D("hPhiBeforeCor", "", 50, 0, 6.283185);
  hPhi[1] = new TH1D("hPhiAfterCor", "", 50, 0, 6.283185);
  fOutputList->Add(hPhi[0]);
  fOutputList->Add(hPhi[1]);
  hEtaPhi[0] = new TH2D("hEtaPhiBeforeCor", "", 50, 0, 6.283185, 16,-0.8,0.8);
  hEtaPhi[1] = new TH2D("hEtaPhiAfterCor", "", 50, 0, 6.283185, 16,-0.8,0.8);
  fOutputList->Add(hEtaPhi[0]);
  fOutputList->Add(hEtaPhi[1]);
  hNhits[0] = new TH1D("hNhitsBeforeCut", "", 200, 0., 200.);
  hNhits[1] = new TH1D("hNhitsAfterCut",  "", 200, 0., 200.);
  fOutputList->Add(hNhits[0]);
  fOutputList->Add(hNhits[1]);
  hPDedx = new TH2D("hPDedx", "", 400, -10., 10., 400, 0, 1000);
  fOutputList->Add(hPDedx);
  if (fPeriod.EqualTo("LHC11h") || fPeriod.EqualTo("LHC10h")){
    hDcaXy[0] = new TH1D("hDcaXyBeforeCut", "", 100, 0., 10.);
    hDcaXy[1] = new TH1D("hDcaXyAfterCut",  "", 100, 0., 10.);
    fOutputList->Add(hDcaXy[0]);
    fOutputList->Add(hDcaXy[1]);
    hDcaZ[0] = new TH1D("hDcaZBeforeCut", "", 100, 0., 10.);
    hDcaZ[1] = new TH1D("hDcaZAfterCut",  "", 100, 0., 10.);
    fOutputList->Add(hDcaZ[0]);
    fOutputList->Add(hDcaZ[1]);
  }

  if (fDoQnESE || fQAV0){
    hQnCentRecenter[0] = new TH2D("hqnV0MRecenter","",100, 0, 100, 72, 0, 12); 
    hQnCentRecenter[1] = new TH2D("hqnV0CRecenter","",100, 0, 100, 72, 0, 12); 
    hQnCentRecenter[2] = new TH2D("hqnV0ARecenter","",100, 0, 100, 72, 0, 12); 
    fOutputList->Add(hQnCentRecenter[0]);
    fOutputList->Add(hQnCentRecenter[1]);
    fOutputList->Add(hQnCentRecenter[2]); 
    hQxCentRecenter[0] = new TH2D("hqxRecenterV0MvsCentSPD","",100, 0, 100, 32, -600, 600);
    hQxCentRecenter[1] = new TH2D("hqxRecenterV0CvsCentSPD","",100, 0, 100, 32, -600, 600);
    hQxCentRecenter[2] = new TH2D("hqxRecenterV0AvsCentSPD","",100, 0, 100, 32, -600, 600);
    hQyCentRecenter[0] = new TH2D("hqyRecenterV0MvsCentSPD","",100, 0, 100, 32, -600, 600);
    hQyCentRecenter[1] = new TH2D("hqyRecenterV0CvsCentSPD","",100, 0, 100, 32, -600, 600);
    hQyCentRecenter[2] = new TH2D("hqyRecenterV0AvsCentSPD","",100, 0, 100, 32, -600, 600);
    for (int i = 0; i < 3; ++i) fOutputList->Add(hQxCentRecenter[i]);
    for (int i = 0; i < 3; ++i) fOutputList->Add(hQyCentRecenter[i]); 
    hQxVtxRecenter[0] = new TH2D("hqxRecenterV0MvsVz","",20, -10, 10, 32, -600, 600);
    hQxVtxRecenter[1] = new TH2D("hqxRecenterV0CvsVz","",20, -10, 10, 32, -600, 600);
    hQxVtxRecenter[2] = new TH2D("hqxRecenterV0AvsVz","",20, -10, 10, 32, -600, 600);
    hQyVtxRecenter[0] = new TH2D("hqyRecenterV0MvsVz","",20, -10, 10, 32, -600, 600);
    hQyVtxRecenter[1] = new TH2D("hqyRecenterV0CvsVz","",20, -10, 10, 32, -600, 600);
    hQyVtxRecenter[2] = new TH2D("hqyRecenterV0AvsVz","",20, -10, 10, 32, -600, 600);
    for (int i = 0; i < 3; ++i) fOutputList->Add(hQxVtxRecenter[i]); 
    for (int i = 0; i < 3; ++i) fOutputList->Add(hQyVtxRecenter[i]);
    for (int iCent = 0; iCent < NCENTBINS; ++iCent){
      hPsiV0Recenter[iCent][0] = new TH1D(Form("hpsiV0MRecenterCent%i", iCent),"", 360, 0., 2*TMath::Pi());
      hPsiV0Recenter[iCent][1] = new TH1D(Form("hpsiV0CRecenterCent%i", iCent),"", 360, 0., 2*TMath::Pi());
      hPsiV0Recenter[iCent][2] = new TH1D(Form("hpsiV0ARecenterCent%i", iCent),"", 360, 0., 2*TMath::Pi());
    } 
    for (int iCent = fCbinLo; iCent < fCbinHg; ++iCent){
      for (int i = 0; i < 3; ++i) fOutputList->Add(hPsiV0Recenter[iCent][i]);    
    } 
  }
  
  //------------------
  // physics
  //------------------
  TH1::SetDefaultSumw2(kTRUE);

  if (fDoQnESE){
    hReQn_thisEvt            = new TH2D("ReQnthisEvt", "", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hImQn_thisEvt             = new TH2D("ImQnthisEvt", "", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMQ_thisEvt                = new TH2D("MQthisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMQ_w1_thisEvt         = new TH2D("MQWeight1thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    pRefFlow_thisEvt         = new TProfile("pRefFlowthisEvt", "", 20, 0, 20, "s");
    pIntd2_thisEvt              = new TProfile("pIntd2thisEvt", "", 30, 0, 30, "s");
    for (int i = 0; i < NCENTBINS; ++i) {
      pRefFlow[i]               = new TProfile(Form("pRefFlowCent%i",i), "", 20, 0, 20, "s");
      pIntd2[i]                   = new TProfile(Form("pIntd2Cent%i",i), "", 30, 0, 30, "s");
      pIntd2Ach[i]             = new TProfile(Form("pIntd2AchCent%i", i), "", 30, 0, 30, "s");
      pAch[i]                     = new TProfile(Form("pAchCent%i", i), "", 30, 0, 30, "s");
    }
    for (int i = fCbinLo; i < fCbinHg; ++i) {
      fOutputList->Add(pRefFlow[i]);
      fOutputList->Add(pIntd2[i]);
      fOutputList->Add(pIntd2Ach[i]);
      fOutputList->Add(pAch[i]);
    }
  };

  if (fDiffV2QA){
    hReQn_thisEvt            = new TH2D("ReQnthisEvt", "", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hImQn_thisEvt             = new TH2D("ImQnthisEvt", "", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMQ_thisEvt                = new TH2D("MQthisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMQ_w1_thisEvt         = new TH2D("MQWeight1thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    pRefFlow_thisEvt         = new TProfile("pRefFlowthisEvt", "", 20, 0, 20, "s");
    hRepnw1_thisEvt        =new TH3D("RepnWeight1thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hImpnw1_thisEvt         =new TH3D("ImpnWeight1thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMp_thisEvt                  =new TH3D("MpthisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8); 
    hMp_w1_thisEvt           =new TH3D("MpWeight1thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);  
    pDiffFlowpQStarGap_thisEvt  = new TProfile3D("pDiffFlowpQStarGapthisEvt","",NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8, "s");
    for (int i = fCbinLo; i < fCbinHg; ++i) {
      pDiffFlowpQStarGap[i] = new TProfile3D(Form("pDiffFlowpQStarGapCent%i",i), "", NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8, "s");
      hDiffFlowYield[i]           = new TH3D(Form("hDiffFlowYieldCent%i",i), "", NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    }
    for (int i = fCbinLo; i < fCbinHg; ++i) {
      fOutputList->Add(pDiffFlowpQStarGap[i]);
      fOutputList->Add(hDiffFlowYield[i]);
    }

    if (fDoTrivalCor){
      hReQPos_thisEvt        = new TH2D("ReQPosthisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
      hImQPos_thisEvt         = new TH2D("ImQPosthisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
      hMQPos_thisEvt          = new TH2D("MQPosthisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
      hMQPos_w1_thisEvt   = new TH2D("MQPosWeight1thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
      hReQNeg_thisEvt       = new TH2D("ReQNegthisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
      hImQNeg_thisEvt        = new TH2D("ImQNegthisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
      hMQNeg_thisEvt         = new TH2D("MQNegthisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
      hMQNeg_w1_thisEvt  = new TH2D("MQNegWeight1thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
      pDiffFlowpQStarGapPos_thisEvt  = new TProfile3D("pDiffFlowpQStarGapPosthisEvt","",NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8, "s");
      pDiffFlowpQStarGapNeg_thisEvt  = new TProfile3D("pDiffFlowpQStarGapNegthisEvt","",NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8, "s");
      for (int i = fCbinLo; i < fCbinHg; ++i) {
        pDiffFlowpQStarGapPos[i]  = new TProfile3D(Form("pDiffFlowpQStarGapPosCent%i", i), "", NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8, "s");
        pDiffFlowpQStarGapNeg[i]  = new TProfile3D(Form("pDiffFlowpQStarGapNegCent%i", i), "", NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8, "s");
      }
      for (int i = fCbinLo; i < fCbinHg; ++i) {
        fOutputList->Add(pDiffFlowpQStarGapPos[i]);
        fOutputList->Add(pDiffFlowpQStarGapNeg[i]);
      }
    } 
  };

  if (fCalcV24VsAch){
    hReQn_thisEvt            = new TH2D("ReQnthisEvt", "", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hImQn_thisEvt             = new TH2D("ImQnthisEvt", "", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hReQ2n_thisEvt          = new TH2D("ReQ2nthisEvt", "", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hImQ2n_thisEvt           = new TH2D("ImQ2nthisEvt", "", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hReQ2nw2_thisEvt      = new TH2D("ReQ2nWeight2thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hImQ2nw2_thisEvt      = new TH2D("ImQ2nWeight2thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8); 
    hReQnw3_thisEvt       = new TH2D("ReQ3nWeight2thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hImQnw3_thisEvt        = new TH2D("ImQ3nWeight2thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMQ_thisEvt                = new TH2D("MQthisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMQ_w1_thisEvt         = new TH2D("MQWeight1thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMQ_w2_thisEvt         = new TH2D("MQWeight2thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMQ_w3_thisEvt         = new TH2D("MQWeight3thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMQ_w4_thisEvt         = new TH2D("MQWeight4thisEvt","", 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    pRefFlow_thisEvt         = new TProfile("pRefFlowthisEvt", "", 20, 0, 20, "s");
    hRepnw1_thisEvt        =new TH3D("RepnWeight1thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hImpnw1_thisEvt         =new TH3D("ImpnWeight1thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMp_thisEvt                 =new TH3D("MpthisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8); 
    hMp_w1_thisEvt          =new TH3D("MpWeight1thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hRep2nw2_thisEvt       =new TH3D("Rep2nWeight2thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hImp2nw2_thisEvt        =new TH3D("Imp2nWeight2thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8); 
    hRepnw3_thisEvt         =new TH3D("RepnWeight3thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hImpnw3_thisEvt          =new TH3D("RepnWeight3thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    hMp_w2_thisEvt           =new TH3D("MpWeight2thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);  
    hMp_w3_thisEvt           =new TH3D("MpWeight3thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);  
    hMp_w4_thisEvt           =new TH3D("MpWeight4thisEvt", "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);  
    pDiffFlowpQStarDirect_thisEvt = new TProfile3D("pDiffFlowpQStarDirectthisEvt","",NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8, "s");
    pDiffFlow4Direct_thisEvt          = new TProfile3D("pDiffFlow4DirectthisEvt","",NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8, "s");
    for (int i = fCbinLo; i < fCbinHg; ++i) {
      for (int j = 0; j < NACHBINS; ++j) pRefFlowAch[i][j]            = new TProfile(Form("pRefFlowCent%iAch%i",i,j), "", 20, 0, 20, "s");
      for (int j = 0; j < NACHBINS; ++j) pDiffFlowpQStarAch[i][j] = new TProfile3D(Form("pDiffFlowpQStarCent%iAch%i",i,j), "",NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8, "s");
      for (int j = 0; j < NACHBINS; ++j) pDiffFlow4Ach[i][j]          = new TProfile3D(Form("pDiffFlow4Cent%iAch%i",i,j), "",NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8, "s"); 
      for (int j = 0; j < NACHBINS; ++j) hDiffFlowYieldAch[i][j]    = new TH3D(Form("hDiffFlowYieldCent%iAch%i",i,j), "" ,NPOIBINS, 0, NPOIBINS, 50, fPtMin, fPtMax, 32, -0.8, 0.8);
    }
    for (int i = fCbinLo; i < fCbinHg; ++i) {
      for (int j = 0; j < NACHBINS; ++j){
        fOutputList->Add(pRefFlowAch[i][j]);          
        fOutputList->Add(pDiffFlowpQStarAch[i][j]);
        fOutputList->Add(pDiffFlow4Ach[i][j]);        
        fOutputList->Add(hDiffFlowYieldAch[i][j]);
      }
    }
  };

  // Dobrin 15o pass2 V0 Calib
  fSPDCutPU = new TF1("fSPDCutPU", "450. + 3.9*x", 0, 50000);
    
  Double_t parV0[8] = {33.4237, 0.953516, 0.0712137, 227.923, 8.9239, -0.00319679, 0.000306314, -7.6627e-07};
  fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);
    
  Double_t parV0CL0[6] = {0.0193587, 0.975914, 0.675714, 0.0292263, -0.000549509, 5.86421e-06};
  fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);
    
  Double_t parFB32[9] = {-812.822, 6.41796, 5421.83, -0.382601, 0.0299686, -26.6249, 321.388, -0.82615, 0.0167828};
  fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
  fMultCutPU->SetParameters(parFB32);


  // Rihan Pile-up function 18?
  // fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

  // Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  // fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  // fV0CutPU->SetParameters(parV0);
  
  // Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  // fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  // fCenCutLowPU->SetParameters(parV0CL0);
  // fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  // fCenCutHighPU->SetParameters(parV0CL0);

  // Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  // fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  // fMultCutPU->SetParameters(parFB32);



  PostData(1,fOutputList); 
  if (fDebug) Printf("Post Data Success!");
}
//------------------------------------------------

void AliAnalysisTaskCMWESE::UserExec(Option_t *)
{
  ResetHists();
  hEvtCount->Fill(1);
  //----------------------------
  // Info
  //----------------------------
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else hEvtCount->Fill(20);
  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else hEvtCount->Fill(21);
  AliAODEvent* fAOD = (AliAODEvent*)InputEvent();
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  } else hEvtCount->Fill(22);
  AliPIDResponse* fPID = handler->GetPIDResponse();
  if (!fPID) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else hEvtCount->Fill(23);
  AliAnalysisUtils* fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else hEvtCount->Fill(24);
  if (fPeriod.EqualTo("LHC15o")){
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    if (!fMultSel) {
      AliError(Form("%s: Could not get AliMultSelection", GetName()));
    } else hEvtCount->Fill(25);
    if (!manager || !handler || !fAOD || !fPID || !fUtils || !fMultSel) return;
  }
  else if (!manager || !handler || !fAOD || !fPID || !fUtils) return;
  hEvtCount->Fill(2);
  if (fDebug) Printf("Info done!");

  //----------------------------
  // Trigger
  //----------------------------
  ULong64_t mask = handler->IsEventSelected();
  Bool_t isTrigselected = false;
        if (fTrigger.EqualTo("kINT7")) { // Run2
      isTrigselected = mask&AliVEvent::kINT7;
        } else if (fTrigger.EqualTo("kMB")) {
      isTrigselected = mask&AliVEvent::kMB; //10h
        } else if (fTrigger.EqualTo("kMB+kCentral+kSemiCentral")) {
      isTrigselected = mask&AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral; //11h
        }
  if(isTrigselected == false) return;
  hEvtCount->Fill(3);
  if (fDebug) Printf("trigger done!");

  //----------------------------
  // event-wise selection
  //----------------------------
  // run number
  fRunNum = fAOD->GetRunNumber();
  fRunNumBin = GetRunNumBin(fRunNum);
  if (fRunNumBin<0) return;
  hRunNumBin->Fill(fRunNumBin);
  hEvtCount->Fill(4);
  if (fDebug) Printf("run nummbr done!");

  // vertex
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  double vx    = fVtx->GetX();
  double vy    = fVtx->GetY();
  double vz    = fVtx->GetZ();
  double dz    = vz - fAOD->GetPrimaryVertexSPD()->GetZ();
  if (fabs(vz)>(double)fZvtxCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors()<1) return;
  hVz[0]->Fill(vz);
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
  // fEventCuts.SetCentralityEstimators("V0M","CL1");
  // if (!fEventCuts->AcceptEvent(fAOD) ) return;
  if (fPeriod.EqualTo("LHC10h")  || fPeriod.EqualTo("LHC11h") ) {if (fabs(dz)>0.5) return;}
  else if (fPeriod.EqualTo("LHC15o") ){
      double covTrc[6],covSPD[6];
      fVtx->GetCovarianceMatrix(covTrc);
      fAOD->GetPrimaryVertexSPD()->GetCovarianceMatrix(covSPD);
      double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
      double errTrc = TMath::Sqrt(covTrc[5]);
      double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
      if (fabs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return;
  }
  hVz[1]->Fill(vz); 
  for (int i = 0; i < 20; ++i) {
      if (vz > -10+i*1 && vz < -10+(i+1)*1) {fVzBin = i; break;}
  }
  if (fVzBin<-990) return;
  hEvtCount->Fill(5);
  if (fDebug) Printf("vertex done!");

  // centrality
  double centV0M = -1, centTRK = -1, centSPD0 = -1, centSPD1 = -1, centV0A = -1;
  if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC11h")){
    centV0M    =  fAOD->GetCentrality()->GetCentralityPercentile("V0M");
    centTRK    =  fAOD->GetCentrality()->GetCentralityPercentile("TRK");
    centSPD0  =  fAOD->GetCentrality()->GetCentralityPercentile("CL0");
    centSPD1  =  fAOD->GetCentrality()->GetCentralityPercentile("CL1");
    centV0A     =  fAOD->GetCentrality()->GetCentralityPercentile("V0A");
  }
  else if (fPeriod.EqualTo("LHC15o") ){
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    centV0M    =  fMultSel->GetMultiplicityPercentile("V0M");
    centTRK    =  fMultSel->GetMultiplicityPercentile("TRK");
    centSPD0 =  fMultSel->GetMultiplicityPercentile("CL0");
    centSPD1 =  fMultSel->GetMultiplicityPercentile("CL1");
  }
  fCent  = centV0M;
  hCentQA[0]->Fill(centV0M,centSPD1);
  hCentQA[2]->Fill(centV0M,centTRK);
  hCentQA[4]->Fill(centV0M,centSPD0);
  hCentQA[6]->Fill(centSPD1,centSPD0);
  if (fabs(fCent-centSPD1)>fCentCut) return;
  hCentQA[1]->Fill(centV0M,centSPD1);
  hCentQA[3]->Fill(centV0M,centTRK);
  hCentQA[5]->Fill(centV0M,centSPD0);
  hCentQA[7]->Fill(centSPD1,centSPD0);
  if (fCent<0) return;
  // cent bin
  fCentBin = (int)fCent/10;
  if (fCentBin<fCbinLo || fCentBin>=fCbinHg) return;
  hCent[0]->Fill(fCent);
  hEvtCount->Fill(6);
  if (fDebug) Printf("centrality done!");

  // pile up
  if (fPeriod.EqualTo("LHC10h")) {if(!RemovalForRun1(fAOD, fUtils) ) return;}
  if (fPeriod.EqualTo("LHC11h")) {if(!RemovalForRun1(fAOD, fUtils) ) return;}
  if (fPeriod.EqualTo("LHC15o")) {
    // hMultCentQA[0]->Fill(fCent, fAOD->GetNumberOfTracks()); // raw Trk Multi Vs Cent(V0M)
    // if (PileUpMultiVertex(fAOD)) return;
    // if (!RejectEvtMultComp(fAOD)) return;
    // hMultCentQA[1]->Fill(fCent, fAOD->GetNumberOfTracks()); // Mult_Cent QA
    // if (!AODPileupCheck (fAOD)) return;
    if (!RejectEvtTFFit (fAOD)) return; // 15o_pass2
  }
  hCent[1]->Fill(fCent);
  hEvtCount->Fill(7);
  if (fDebug) Printf("pile-up done!");

  // analysis 
  if (!AnalyzeAOD(fAOD, fVtx)) return;
  hEvtCount->Fill(8);

  //------------------
  // done
  //------------------
  if (fDebug) Printf("analysis done!");
  PostData(1,fOutputList);
}

//---------------------------------------------------

void AliAnalysisTaskCMWESE::ResetHists()
{
  if (fDoQnESE){
    hReQn_thisEvt->Reset();        
    hImQn_thisEvt->Reset();        
    hMQ_thisEvt->Reset();           
    hMQ_w1_thisEvt->Reset();     
    pRefFlow_thisEvt->Reset();    
    pIntd2_thisEvt->Reset();         
  };

  if (fDiffV2QA){
    hReQn_thisEvt->Reset();       
    hImQn_thisEvt->Reset();       
    hMQ_thisEvt->Reset();          
    hMQ_w1_thisEvt->Reset();    
    hRepnw1_thisEvt->Reset();   
    hImpnw1_thisEvt->Reset();    
    hMp_thisEvt->Reset();           
    hMp_w1_thisEvt->Reset(); 
    pDiffFlowpQStarGap_thisEvt->Reset();
    pRefFlow_thisEvt->Reset();     
    if (fDoTrivalCor){
      hReQPos_thisEvt->Reset();      
      hImQPos_thisEvt->Reset();       
      hMQPos_thisEvt->Reset();        
      hMQPos_w1_thisEvt->Reset(); 
      hReQNeg_thisEvt->Reset();      
      hImQNeg_thisEvt->Reset();      
      hMQNeg_thisEvt->Reset();       
      hMQNeg_w1_thisEvt->Reset(); 
      pDiffFlowpQStarGapPos_thisEvt->Reset();
      pDiffFlowpQStarGapNeg_thisEvt->Reset();
    }
  };

  if (fCalcV24VsAch){
    hReQn_thisEvt->Reset();        
    hImQn_thisEvt->Reset();        
    hReQ2n_thisEvt->Reset();      
    hImQ2n_thisEvt->Reset();      
    hReQ2nw2_thisEvt->Reset(); 
    hImQ2nw2_thisEvt->Reset();  
    hReQnw3_thisEvt->Reset();   
    hImQnw3_thisEvt->Reset();    
    hMQ_thisEvt->Reset();           
    hMQ_w1_thisEvt->Reset();     
    hMQ_w2_thisEvt->Reset();     
    hMQ_w3_thisEvt->Reset();     
    hMQ_w4_thisEvt->Reset();
    pRefFlow_thisEvt->Reset();   
    hRepnw1_thisEvt->Reset();    
    hImpnw1_thisEvt->Reset();     
    hMp_thisEvt->Reset();            
    hMp_w1_thisEvt->Reset();      
    hRep2nw2_thisEvt->Reset();  
    hImp2nw2_thisEvt->Reset();   
    hRepnw3_thisEvt->Reset();    
    hImpnw3_thisEvt->Reset();     
    hMp_w2_thisEvt->Reset();      
    hMp_w3_thisEvt->Reset();      
    hMp_w4_thisEvt->Reset();      
    pDiffFlowpQStarDirect_thisEvt->Reset();
    pDiffFlow4Direct_thisEvt->Reset();         
  };
  fRunNum=-999;
  fRunNumBin=-999;
  fVzBin=-999;
  fCentBin=-999;
  fCent=-999;
  fQnBin=-999;
  fNegEtaQ=(-999., -999.);
  fNegEtaQStar=(-999., -999.);
  fPosEtaQ=(-999., -999.);
  fPosEtaQStar=(-999., -999.);
  fNegEtaMQ=-999.;
  fPosEtaMQ=-999.; 
  return;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::AnalyzeAOD(AliAODEvent* fAOD, AliAODVertex* fVtx)
{
  const int MAXNTRK = 1e5;
  int nTrk = fAOD->GetNumberOfTracks();
  if (nTrk<10 || nTrk>MAXNTRK) return false; 
  double nPosTrk = 0., nNegTrk = 0.;
  double isSltTrk[MAXNTRK] = {0};
  double weightTrk[MAXNTRK] = {0};
  double Mult = 0.;
  for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrk);
    if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }
    if (!track->TestFilterBit(fFltbit)) continue; 
    if (!AcceptAODTrack(fAOD, track, fVtx)) continue;
    int    charge = track->Charge();
    //------------------
    // NUE & NUA
    //------------------
    double phi    = track->Phi();
    double weight=1;
    double vz    = fVtx->GetZ();
    double pt     = track->Pt();
    double eta    = track->Eta();
    hPhi[0]->Fill(phi);
    hEtaPhi[0]->Fill(phi,eta);
    if (fDoNUE) {
      double wEffi = GetNUECor(charge, pt);
      if (wEffi<0) continue;
      else weight *= wEffi;
    }
    if (fDoNUA){
      double wAcc = GetNUACor(charge, phi, eta, vz);
      if (wAcc<0) continue;
      else weight *= wAcc;
      hPhi[1]->Fill(phi, wAcc);
      hEtaPhi[1]->Fill(phi,eta,wAcc);
    } 
    weightTrk[iTrk] = weight;
    //------------------
    // for ach
    //------------------
    if (charge>0) nPosTrk++; 
    if (charge<0) nNegTrk++; 
    //------------------
    // fill QC hists
    //------------------
    double cosnphi = weight * TMath::Cos(fHarmonic*phi);
    double sinnphi = weight * TMath::Sin(fHarmonic*phi);
    double cos2nphi = weight * TMath::Cos(2*fHarmonic*phi); 
    double sin2nphi = weight * TMath::Sin(2*fHarmonic*phi); 

    // (Qx, Qy)
    hReQn_thisEvt            ->Fill(pt, eta, cosnphi); // w^1
    hImQn_thisEvt             ->Fill(pt, eta, sinnphi);
    if (fCalcV24VsAch){
      hReQ2n_thisEvt          ->Fill(pt,eta,cos2nphi); 
      hImQ2n_thisEvt           ->Fill(pt, eta, sin2nphi); 
      hReQ2nw2_thisEvt        ->Fill(pt,eta,cos2nphi *weight ); // w^2
      hImQ2nw2_thisEvt         ->Fill(pt, eta, sin2nphi *weight ); 
      hReQnw3_thisEvt          ->Fill(pt, eta, cosnphi *pow(weight,2.) ); // w^3
      hImQnw3_thisEvt           ->Fill(pt, eta, sinnphi *pow(weight,2.) );
    }
    //MQ
    hMQ_thisEvt               ->Fill(pt, eta);
    hMQ_w1_thisEvt        ->Fill(pt, eta, weight);               // w^1 
    if (fCalcV24VsAch){
      hMQ_w2_thisEvt        ->Fill(pt, eta, pow(weight,2.) ); // w^2
      hMQ_w3_thisEvt        ->Fill(pt, eta, pow(weight,3.) ); // w^3
      hMQ_w4_thisEvt        ->Fill(pt, eta, pow(weight,4.) ); // w^4
    }
    if (fDiffV2QA && fDoTrivalCor){    
      if (charge>0) {
        hReQPos_thisEvt              ->Fill(pt, eta, cosnphi);
        hImQPos_thisEvt               ->Fill(pt, eta, sinnphi);
        hMQPos_thisEvt                ->Fill(pt, eta);
        hMQPos_w1_thisEvt         ->Fill(pt, eta, weight);
      } else if (charge<0){
        hReQNeg_thisEvt              ->Fill(pt, eta, cosnphi);
        hImQNeg_thisEvt               ->Fill(pt, eta, sinnphi);
        hMQNeg_thisEvt                ->Fill(pt, eta);
        hMQNeg_w1_thisEvt         ->Fill(pt, eta, weight);
      }
    };

    if (fDiffV2QA || fCalcV24VsAch){
      for (int iPOIBin = 0; iPOIBin < NPOIBINS; ++iPOIBin) {
        bool isRightPID = false;
        if (iPOIBin==0 && charge>0) isRightPID = true;
        if (iPOIBin==1 && charge<0) isRightPID = true;
        if (!isRightPID) continue;
        // (px, py)           
        hRepnw1_thisEvt                  ->Fill(iPOIBin, pt, eta, cosnphi); // w^1
        hImpnw1_thisEvt                   ->Fill(iPOIBin, pt, eta, sinnphi);
        hRep2nw2_thisEvt                ->Fill(iPOIBin, pt, eta, cos2nphi * weight); // w^2
        hImp2nw2_thisEvt                 ->Fill(iPOIBin, pt, eta, sin2nphi * weight); 
        hRepnw3_thisEvt                  ->Fill(iPOIBin, pt, eta, cosnphi *pow(weight,2.) ); // w^3
        hImpnw3_thisEvt                   ->Fill(iPOIBin, pt, eta, sinnphi *pow(weight,2.) );
        // Mp
        hMp_thisEvt                     ->Fill(iPOIBin, pt, eta); 
        hMp_w1_thisEvt              ->Fill(iPOIBin, pt, eta, weight);  
        hMp_w2_thisEvt              ->Fill(iPOIBin, pt, eta, pow(weight,2.) );  
        hMp_w3_thisEvt              ->Fill(iPOIBin, pt, eta, pow(weight,3.) );  
        hMp_w4_thisEvt              ->Fill(iPOIBin, pt, eta, pow(weight,4.) );      
      }      
    };

    Mult++;
    isSltTrk[iTrk] = 1;
  }  
  hEvtCount->Fill(9);

  //------------------
  // Ref Cumulant
  //------------------
  if (fDoQnESE){
    if (!DoCumulantswtGap()) return false;
    hEvtCount->Fill(10);
  };

  if (fDiffV2QA){
    if (!DoCumulantswtGap()) return false;
    hEvtCount->Fill(10); 
    pDiffFlowpQStarGap[fCentBin]                                     ->Add(pDiffFlowpQStarGap_thisEvt); 
    hDiffFlowYield[fCentBin]                                               ->Add(hMp_w1_thisEvt); 
    if (fDoTrivalCor) pDiffFlowpQStarGapPos[fCentBin]    ->Add(pDiffFlowpQStarGapPos_thisEvt); 
    if (fDoTrivalCor) pDiffFlowpQStarGapNeg[fCentBin]    ->Add(pDiffFlowpQStarGapNeg_thisEvt); 
    hEvtCount->Fill(14);
    return true;
  };

  //------------------
  // Ach
  //------------------
  if (nPosTrk<1 || nNegTrk<1) return false;
  double mAch = (double)(nPosTrk-nNegTrk)/(nPosTrk+nNegTrk);
  if (fabs(mAch)>1) return false;
  fAchBin = GetAchBin(mAch);
  if (fAchBin < 0 ) return false;
  hEvtCount->Fill(11);

  //------------------------
  // Direct v2{2}/v2{4}
  //------------------------
  if (fCalcV24VsAch){
    if (!DoCumulantsDirect()) return false;
    hEvtCount->Fill(10);
    pRefFlowAch[fCentBin][fAchBin]                          ->Add(pRefFlow_thisEvt);
    pDiffFlowpQStarAch[fCentBin][fAchBin]               ->Add(pDiffFlowpQStarDirect_thisEvt); 
    pDiffFlow4Ach[fCentBin][fAchBin]                        ->Add(pDiffFlow4Direct_thisEvt); 
    hDiffFlowYieldAch[fCentBin][fAchBin]                  ->Add(hMp_w1_thisEvt);
    hEvtCount->Fill(14);
    return true;
  };

  //------------------
  // Qn vector
  //------------------
  if (fDoQnESE){
    if (!CalcQnVectorV0(fAOD, fVtx, mAch, Mult)) return false;
    hEvtCount->Fill(12);  
  
    //------------------
    // loop track again
    //------------------
    for (int iTrk = 0; iTrk < nTrk; ++iTrk) {
      AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrk);
      if (!track) {
        AliError(Form("%s: Could not get Track", GetName()));
        continue;
      }
      if (!track->TestFilterBit(fFltbit)) continue; 
      if (isSltTrk[iTrk] == 0) continue;
  
      //------------------
      // PID
      //------------------
      // double nSigmaTPCPi = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kPion);
      // double nSigmaTPCK  = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kKaon);
      // double nSigmaTPCP  = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kProton);
      // double nSigmaTPCDe = fPID->NumberOfSigmasTPC(track, (AliPID::EParticleType)AliPID::kDeuteron); 
      // double nSigmaTOFPi = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kPion);
      // double nSigmaTOFK  = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kKaon);
      // double nSigmaTOFP  = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kProton);
      // double nSigmaTOFDe = fPID->NumberOfSigmasTOF(track, (AliPID::EParticleType)AliPID::kDeuteron); 
      // bool isPi = false, isK = false, isP = false, isDe = false;
      // // if (fabs(nSigmaTPCPi) <= 2. && fabs(nSigmaTOFPi) <= 2.) isPi = true;
      // // if (fabs(nSigmaTPCK)  <= 2. && fabs(nSigmaTOFK)  <= 2.) isK  = true;
      // // if (fabs(nSigmaTPCP)  <= 2. && fabs(nSigmaTOFP)  <= 2.) isP  = true;
      // if (pt<0.6 && fabs(nSigmaTPCPi)<=2.) isPi = true;
      // if (pt>0.6 && pt<2. && fabs(nSigmaTPCPi)<=2.5 && fabs(nSigmaTOFPi)<=2.) isPi = true;
      // if (pt<0.6 && fabs(nSigmaTPCK)<=2.) isK = true;
      // if (pt>0.6 && pt<2 && fabs(nSigmaTPCK)<=2.5 && fabs(nSigmaTOFK)<=2.) isK = true;
      // if (pt<0.7 && fabs(nSigmaTPCP)<=2.) isP = true;
      // if (pt>0.7 && pt<3 && fabs(nSigmaTPCP)<=2.5 && fabs(nSigmaTOFP)<=2.) isP = true;
      // if (pt<1 && fabs(nSigmaTPCDe)<=2.) isDe = true; 
      // if (pt>1 && pt<5 &&fabs(nSigmaTPCDe)<=2.5 && fabs(nSigmaTOFDe)<=2. ) isDe = true; 
      double phi    = track->Phi();
      double eta    = track->Eta();
      int charge     = track->Charge(); 
      double weight1 = weightTrk[iTrk];
      double cosnphi = weight1 * TMath::Cos(fHarmonic*phi);
      double sinnphi = weight1 * TMath::Sin(fHarmonic*phi);
      double cos2nphi = weight1*weight1 * TMath::Cos(2*fHarmonic*phi); 
      double sin2nphi = weight1*weight1 * TMath::Sin(2*fHarmonic*phi); 
  
      // p pStar track-wise
      TComplex p_iTrk(cosnphi, sinnphi);
      TComplex pStar_iTrk = TComplex::Conjugate(p_iTrk);
      double d2_iTrk = 0., mM = 0.;
      if (eta < 0.) {
        mM = weight1* fPosEtaMQ;
        d2_iTrk = ((p_iTrk*fPosEtaQStar).Re()) / mM;
      }
      else{
        mM = weight1* fNegEtaMQ;
        d2_iTrk = ((p_iTrk*fNegEtaQStar).Re()) / mM;
      } 
  
      if (charge>0){
        pIntd2_thisEvt->Fill((double)(2*fQnBin+0.5), d2_iTrk, weight1); 
        pAch[fCentBin]->Fill((double)(2*fQnBin+0.5), mAch, weight1);
        pIntd2_thisEvt->Fill(28.5, d2_iTrk, weight1); //unbiased
        pAch[fCentBin]->Fill(28.5, mAch, weight1);
        // keep Ach entries consistent with d2Ach & d2 etc.
      } 
      else{
        pIntd2_thisEvt->Fill((double)(2*fQnBin+1.5), d2_iTrk, weight1); 
        pAch[fCentBin]->Fill((double)(2*fQnBin+1.5), mAch, weight1);
        pIntd2_thisEvt->Fill(29.5, d2_iTrk, weight1); //unbiased
        pAch[fCentBin]->Fill(29.5, mAch, weight1);
      } 
  
    } // loop iTrk end
    hEvtCount->Fill(13);
    //------------------
    // covariance
    //------------------
    pRefFlow[fCentBin]                                                  ->Add(pRefFlow_thisEvt);
    pIntd2[fCentBin]                                                       ->Add(pIntd2_thisEvt);
    pIntd2_thisEvt                                                          ->Scale(mAch); 
    pIntd2Ach[fCentBin]                                                 ->Add(pIntd2_thisEvt);
    hEvtCount->Fill(14);
    return true;
  }
  return true;
}

//---------------------------------------------------

double AliAnalysisTaskCMWESE::GetNUECor(int charge, double pt)
{
  double weightNUE = 1;
  if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC11h")){
    // hNUEweightPlus = (TH1D*)fListNUE->FindObject("effVsPtPlus");
    // hNUEweightMinus = (TH1D*)fListNUE->FindObject("effVsfPtMinus");
    if (charge>0){
      hNUEweightPlus = (TH1D*)fListNUE->FindObject(Form("effVsPt_cent%iPlus",fCentBin));
      if (!hNUEweightPlus) return -1;
      int ptBin = hNUEweightPlus->GetXaxis()->FindBin(pt);
      if (hNUEweightPlus->GetBinContent(ptBin)>0){
        weightNUE = 1./ hNUEweightPlus->GetBinContent(ptBin); 
      } 
      else return -1;
    }
    if (charge<0){
      hNUEweightMinus = (TH1D*)fListNUE->FindObject(Form("effVsPt_cent%iMinus",fCentBin));
      if (!hNUEweightMinus) return -1;
      int ptBin = hNUEweightMinus->GetXaxis()->FindBin(pt);
      if (hNUEweightMinus->GetBinContent(ptBin)>0){
        weightNUE = 1./ hNUEweightMinus->GetBinContent(ptBin); 
      }  
      else return -1;  
    }
  } 
  else if (fPeriod.EqualTo("LHC15o")){
    if (charge>0){
      hNUEweightPlus = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgPos");
      if (!hNUEweightPlus) return -1;
      int ptBin = hNUEweightPlus->GetXaxis()->FindBin(pt);
      if (hNUEweightPlus->GetBinContent(ptBin)>0){
        weightNUE = 1./ hNUEweightPlus->GetBinContent(ptBin); 
      } 
      else return -1;
    }
    if (charge<0){
      hNUEweightMinus = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgNeg");
      if (!hNUEweightMinus) return -1;
      int ptBin = hNUEweightMinus->GetXaxis()->FindBin(pt);
      if (hNUEweightMinus->GetBinContent(ptBin)>0){
        weightNUE = 1./ hNUEweightMinus->GetBinContent(ptBin); 
      } 
      else return -1;
    }
  }
  return weightNUE;
}

//---------------------------------------------------

double AliAnalysisTaskCMWESE::GetNUACor(int charge, double phi, double eta, double vz)
{
  double weightNUA = 1;
  if (fVzBin<0 || fCentBin<0 || fRunNum<0) return -1;
  if (fPeriod.EqualTo("LHC10h")){
    hNUAweightPlus = (TH2D*)fListNUA1->FindObject(Form("weightdPhidEta_run%i_cent0_vz%i_plus",fRunNum, fVzBin));
    hNUAweightMinus = (TH2D*)fListNUA1->FindObject(Form("weightdPhidEta_run%i_cent0_vz%i_minus",fRunNum, fVzBin));
    if(!hNUAweightPlus || !hNUAweightMinus) return -1;
    if (charge>0){ 
      int phiBin = hNUAweightPlus->GetXaxis()->FindBin(phi);
      int etaBin = hNUAweightPlus->GetYaxis()->FindBin(eta);  
      if (hNUAweightPlus->GetBinContent(phiBin, etaBin)>0) weightNUA = hNUAweightPlus->GetBinContent(phiBin, etaBin);
      return weightNUA;
    } else if (charge<0){
      int phiBin = hNUAweightMinus->GetXaxis()->FindBin(phi);
      int etaBin = hNUAweightMinus->GetYaxis()->FindBin(eta);  
      if (hNUAweightMinus->GetBinContent(phiBin, etaBin)>0) weightNUA = hNUAweightMinus->GetBinContent(phiBin, etaBin);
      return weightNUA;
    } 
  } else if (fPeriod.EqualTo("LHC11h")){
    if (charge>0){ 
      hNUAweightPlus = (TH2D*) fListNUA1->FindObject(Form("weightdPhidEta_run170163_cent%i_vz%i_plus",fCentBin,fVzBin));
      if (!hNUAweightPlus) return -1;
      int iBinNUA = hNUAweightPlus->FindBin(phi,eta); 
      if (hNUAweightPlus->GetBinContent(iBinNUA)>0) weightNUA = (double)hNUAweightPlus->GetBinContent(iBinNUA);
      return  weightNUA;
    } else if (charge<0){
      hNUAweightMinus = (TH2D*) fListNUA1->FindObject(Form("weightdPhidEta_run170163_cent%i_vz%i_minus",fCentBin,fVzBin));
      if (!hNUAweightMinus) return -1;
      int iBinNUA = hNUAweightMinus->FindBin(phi,eta); 
      if (hNUAweightMinus->GetBinContent(iBinNUA)>0) weightNUA = (double)hNUAweightMinus->GetBinContent(iBinNUA);
      return weightNUA;
    } 
  } else if (fPeriod.EqualTo("LHC15o")){ // Rihan and Protty 's NUA Results
    if (charge>0){ 
      hCorrectNUAPos = (TH3F*) fListNUA1->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Pos_Cent0_Run%i",fRunNum));
      if (!hCorrectNUAPos) return -1;
      int iBinNUA = hCorrectNUAPos->FindBin(vz,phi,eta); 
      if (hCorrectNUAPos->GetBinContent(iBinNUA)>0) weightNUA = (double)hCorrectNUAPos->GetBinContent(iBinNUA);
      return  weightNUA;
    } else if (charge<0){
      hCorrectNUANeg = (TH3F*) fListNUA1->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Neg_Cent0_Run%i",fRunNum));
      if (!hCorrectNUANeg) return -1;
      int iBinNUA = hCorrectNUANeg->FindBin(vz,phi,eta); 
      if (hCorrectNUANeg->GetBinContent(iBinNUA)>0) weightNUA = (double)hCorrectNUANeg->GetBinContent(iBinNUA);
      return weightNUA;
    } 
    // In Rihan and Protty 's NUA results, the phi distribution is independent on centrality and particle charge
  }
  return weightNUA;
}

//---------------------------------------------------

int AliAnalysisTaskCMWESE::GetRunNumBin(int runNum)
{
  int runNBin=-1;
  if (fPeriod.EqualTo("LHC10h") ){
      int IRunNumList[90]={
        139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310,
        139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870,
        138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582,
        138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201,
        138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693,
        137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541,
        137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243,
        137236, 137235, 137232, 137231, 137162, 137161};
      for (int i = 0; i < 90; ++i) {
        if (runNum==IRunNumList[i]) {runNBin=i; break;}
        else continue;
      }
  } else if (fPeriod.EqualTo("LHC11h") ){
      int IRunNumList[68]={
        167915, 168115, 168460, 169035, 169238, 169859, 170228 , 167920, 168310, 168464, 169091, 169411, 
        169923, 170230, 167985, 168311, 168467, 169094, 169415, 170027, 170268, 167987, 168322, 168511, 
        169138, 169417, 170081, 170269, 167988, 168325, 168512, 169144, 169835, 170155, 170270, 168069, 
        168341, 168514, 169145, 169837, 170159, 170306, 168076, 168342, 168777, 169148, 169838, 170163, 
        170308, 168105, 168361, 168826, 169156, 169846, 170193, 170309, 168107, 168362, 168988, 169160, 
        169855, 170203, 168108 , 168458, 168992, 169167, 169858, 170204};
        for (int i = 0; i < 68; ++i) {
          if (runNum==IRunNumList[i]) {runNBin=i; break;}
          else continue;
        }
  } else if (fPeriod.EqualTo("LHC15o") ){
      int IRunNumList[138]={
        246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 
        246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 
        246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246434, 
        246431, 246424, 246392, 246391, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 
        246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 
        246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 
        245952, 245949, 245923, 245833, 245831, 245829, 245793, 245785, 245775, 245766, 245759, 245752, 
        245731, 245729, 245705, 245702, 245692, 245683, 245554, 245545, 245544, 245543, 245542, 245540, 
        245535, 245507, 245505, 245504, 245501, 245497, 245496, 245454, 245453, 245450, 245446, 245441, 
        245411, 245410, 245409, 245407, 245401, 245397, 245396, 245353, 245349, 245347, 245346, 245345, 
        245343, 245259, 245233, 245232, 245231, 245152, 245151, 245146, 245145, 245068, 245066, 245064, 
        244983, 244982, 244980, 244975, 244918, 244917};
      for (int i = 0; i < 138; ++i) {
        if (runNum==IRunNumList[i]) {runNBin=i; break;}
        else continue;
      }
  }
  return runNBin;
}

//---------------------------------------------------
int AliAnalysisTaskCMWESE::GetAchBin(double mAch)
{
  // split equally
  int AchBin = -1;
  for (int i = 0; i < NACHBINS; ++i) {
    if (mAch>-0.12+0.02*i && mAch<-0.12+0.02*(i+1)) {
      AchBin=i; break;
    }
  }
  return AchBin;
}

//---------------------------------------------------

double AliAnalysisTaskCMWESE::GetEventPlane(double qx, double qy) 
{
  double psi = TMath::ATan2(qy, qx)/fHarmonic;

  if (psi < 0) return psi+TMath::Pi();
  else return psi;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::RemovalForRun1(AliAODEvent* fAOD, AliAnalysisUtils* fUtils)
{
  // pileup
  fUtils->SetUseOutOfBunchPileUp(true);
  fUtils->SetUseMVPlpSelection(true);
  // fUtils->SetMinPlpContribMV(5);
  bool isPileup = true;
  if (fPeriod.EqualTo("LHC10h")) isPileup = fUtils->IsPileUpEvent(fAOD);
  if (fPeriod.EqualTo("LHC11h")) isPileup = fAOD->IsPileupFromSPD();
  // bool isPileup = fUtils->IsPileUpMV(fAOD); // pp, p-Pb
  if (isPileup) return false;
  return true;   
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::RejectEvtMultComp(AliAODEvent* fAOD) // 15o_pass1, old pile-up
{
   // TPC cluster cut
    Int_t multEsd = ((AliAODHeader*)fAOD->GetHeader())->GetNumberOfESDTracks(); // multESD
    const Int_t nTrk = fAOD->GetNumberOfTracks();
    Int_t multTPC=0; // FB128 + Common Track Cuts
    Int_t multTPCFE=0; // FB1 + Common Track Cuts + chi2 > 0.2
    Int_t multGlobal=0; // FB16 + Common Track Cuts + Dca Cuts
    for (Int_t it1 = 0; it1 < nTrk; it1++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it1);
      if (!aodTrk) continue;
      if (aodTrk->TestFilterBit(128)) multTPC++;
      double Eta  = aodTrk->Eta();
      double Pt    = aodTrk->Pt();
      double Phi  = aodTrk->Phi();
      double charge = aodTrk->Charge();
      if(Pt<0.2 || Pt>5.0 || TMath::Abs(Eta)>0.8 || aodTrk->GetTPCNcls()<fNclsCut || aodTrk->GetTPCsignal()<10.0) continue;
      if(aodTrk->TestFilterBit(1) && aodTrk->Chi2perNDF()>0.2)  multTPCFE++;
      if(!aodTrk->TestFilterBit(16) || aodTrk->Chi2perNDF()<0.1)   continue;
      Double_t dca[2] = {-99., -99.};
      Double_t cov[3] = {-99., -99., -99.};
      Double_t magField = fAOD->GetMagneticField();
      if(magField!=0){     
        if(aodTrk->PropagateToDCA(fAOD->GetPrimaryVertex(), magField, 100., dca, cov) && TMath::Abs(dca[0]) < 0.3 && TMath::Abs(dca[1]) < 0.3) multGlobal++;    
      }
    }

    hMultMultQA[0]->Fill(multTPC,multEsd );
    hMultMultQA[1]->Fill(multGlobal,multTPCFE );

    if(fMultComp.EqualTo("pileupByEDSTPC128") ){ // Rihan
      if ((Double_t)(multEsd*(1/3.45) - 90.) < (Double_t)multTPC) 
      {
        hMultMultQA[3]->Fill(multTPC,multEsd );
        return true;
      }
      else return false;
    }

    if(fMultComp.EqualTo("pileupByGlobalTPC1") ) { // A.Dobrin
      if (multTPCFE-1.78*multGlobal<62.87 && multTPCFE-1.48*multGlobal>-36.73){
        hMultMultQA[4]->Fill(multGlobal,multTPCFE );
        return true;  
      } 
      else return false;
    } 

    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::RejectEvtTFFit(AliAODEvent* fAOD)
{
  Float_t centV0M=-1.;
  Float_t centCL1=-1.;
  Float_t centCL0=-1.;

  AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
 
  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
  }
  centV0M = (Float_t) fCent;
  centCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

  Int_t nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();
    
  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTrk = 0;
  for (Int_t it = 0; it < nTracks; it++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
      if (!aodTrk){
          delete aodTrk;
          continue;
      }
        
      if (aodTrk->TestFilterBit(32)){
        if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
      }
  }

  hMultCentQA[0]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)

  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;

  // pile-up cuts
  if (centCL0 < fCenCutLowPU->Eval(centV0M)) return false;        
  if (centCL0 > fCenCutHighPU->Eval(centV0M)) return false;
           
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) return false;
           
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;

  if (Float_t(multTrk) < fMultCutPU->Eval(centV0M)) return false;

  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;

  if (fAOD->IsIncompleteDAQ()) return false;

  hMultCentQA[1]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::RejectEvtTPCITSfb32TOF (AliAODEvent* fAOD)
{
    //TOD+FB32 pile-up removal
    // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
    Int_t multTrk=0;
    Int_t multTrkTOF=0;
    int nTrk = fAOD->GetNumberOfTracks();
    for (Int_t it2 = 0; it2 < nTrk; it2++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it2);
      if (!aodTrk) continue;
      if (aodTrk->TestFilterBit(32)){
        multTrk++;
        if ( TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10 && aodTrk->GetTOFsignal() >= 12000 && aodTrk->GetTOFsignal() <= 25000) multTrkTOF++;
        else return false; 
      }
    }
    hMultMultQA[2]->Fill(multTrkTOF, nTrk); 
    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::AODPileupCheck (AliAODEvent* fAOD)
{
  Int_t isPileup = fAOD->IsPileupFromSPD(3);
  if (isPileup !=0 && fPeriod.EqualTo("LHC16t")) return false; // LHC16t : pPb
  if (fAOD->IsIncompleteDAQ()) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fPeriod.EqualTo("LHC15o")){
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    if(!fMultSel->GetThisEventIsNotPileup()) return false;
    if(!fMultSel->GetThisEventIsNotPileupMV()) return false;
    if(!fMultSel->GetThisEventIsNotPileupInMultBins()) return false;
    if(!fMultSel->GetThisEventHasNoInconsistentVertices()) return false;
    if(!fMultSel->GetThisEventPassesTrackletVsCluster()) return false;
    if(!fMultSel->GetThisEventIsNotIncompleteDAQ()) return false;
    if(!fMultSel->GetThisEventHasGoodVertex2016()) return false;
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::PileUpMultiVertex(AliAODEvent* fAOD)
{
  // check for multi-vertexer pile-up
  const int        kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if(!(nPlp=fAOD->GetNumberOfPileupVerticesTracks()))
  return false;

  vtPrm = fAOD->GetPrimaryVertex();
  if(vtPrm == fAOD->GetPrimaryVertexSPD())
  return true;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for(int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)fAOD->GetPileupVertexTracks(ipl);
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF()    > kMaxPlpChi2)    continue;
    //int bcPlp = vtPlp->GetBC();
    //if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2)
    // return kTRUE; // pile-up from other BC

    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist)        continue;

    return true; // pile-up: well separated vertices
  }
  return false;
}

//---------------------------------------------------

double AliAnalysisTaskCMWESE::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
    // calculate sqrt of weighted distance to other vertex
    if (!v0 || !v1) {
        printf("One of vertices is not valid\n");
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
    if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
    dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
        +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
    return dist>0 ? TMath::Sqrt(dist) : -1;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::RemovalForLHC18 (AliAODEvent* fAOD)
{}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::RemovalForpPb (AliAODEvent* fAOD)
{}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::AcceptAODTrack(AliAODEvent* fAOD, AliAODTrack *track, AliAODVertex* fVtx)
{
    //------------------
    // track cut  
    //------------------
    double pt     = track->Pt();
    double eta    = track->Eta();
    int    nhits  = track->GetTPCNcls();
    double dedx   = track->GetTPCsignal();
    double chi2   = track->Chi2perNDF();
    int    charge = track->Charge();
    hEta[0]->Fill(eta);
    hNhits[0]->Fill(nhits);
    if(pt < fPtMin || pt > fPtMax || fabs(eta)>fEtaCut || fabs(nhits)<fNclsCut || chi2<fChi2Lo || chi2>fChi2Hg || dedx<fDedxCut) return false;
    hPt->Fill(pt);
    hEta[1]->Fill(eta);
    hNhits[1]->Fill(nhits);
    hPDedx->Fill(track->P()*charge, dedx);
    if (fPeriod.EqualTo("LHC11h") || fPeriod.EqualTo("LHC10h")){
      //------------------
      // dca cut
      //------------------
      double mag = fAOD->GetMagneticField(); 
      double dcaxy  = 999.;
      double dcaz   = 999.;
      double r[3];
      double dca[2];
      double cov[3];
      double vx    = fVtx->GetX();
      double vy    = fVtx->GetY();
      double vz    = fVtx->GetZ();
      bool proptodca = track->PropagateToDCA(fVtx, mag, 100., dca, cov);
      if (track->GetXYZ(r)) {
        dcaxy = r[0];
        dcaz  = r[1];
      } else {
        double dcax = r[0] - vx;
        double dcay = r[1] - vy;
        dcaz  = r[2] - vz;
        dcaxy = sqrt(dcax*dcax + dcay*dcay);
        // dcaxy = dca[0];
      }
      hDcaXy[0]->Fill(dcaxy);
      if (fabs(dcaxy)>fDcaCutxy) return false;
      hDcaXy[1]->Fill(dcaxy);
      hDcaZ[0]->Fill(dcaz);
      if (fabs(dcaz)>fDcaCutz) return false;
      hDcaZ[1]->Fill(dcaz);     
    }

    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::DoCumulantswtGap()
{
  if (hMQ_thisEvt->Integral(1, 100, 1, 32) < 2) return false;
  double allM  =  hMQ_w1_thisEvt->Integral(1, 100, 1, 32);

  // Q
  double allReQ  = hReQn_thisEvt ->Integral(1, 100, 1, 32);
  double allImQ  = hImQn_thisEvt ->Integral(1, 100, 1, 32);
  TComplex Q(allReQ, allImQ);
  TComplex QStar = TComplex::Conjugate(Q);

  double posEtaRe = hReQn_thisEvt->Integral(1, 100, hReQn_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  double posEtaIm = hImQn_thisEvt->Integral(1, 100, hImQn_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  if (hMQ_thisEvt->Integral(1, 100,  hMQ_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32) < 2) return false;
  fPosEtaMQ =  hMQ_w1_thisEvt->Integral(1, 100,  hMQ_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  fPosEtaQ = TComplex(posEtaRe, posEtaIm);
  fPosEtaQStar = TComplex::Conjugate(fPosEtaQ);

  double negEtaRe = hReQn_thisEvt->Integral(1, 100, 1, hReQn_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  double negEtaIm = hImQn_thisEvt->Integral(1, 100, 1, hImQn_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  if (hMQ_thisEvt->Integral(1, 100, 1,  hMQ_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6)) < 2) return false;
  fNegEtaMQ =  hMQ_w1_thisEvt->Integral(1, 100, 1,  hMQ_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  fNegEtaQ = TComplex(negEtaRe, negEtaIm);
  fNegEtaQStar = TComplex::Conjugate(fNegEtaQ);

  // c{2} wt gap
  pRefFlow_thisEvt->Fill(0.,  (fNegEtaQ*fPosEtaQStar).Re()       / (fNegEtaMQ*fPosEtaMQ),       (fNegEtaMQ*fPosEtaMQ));
  // d{2} wt gap
  if (fDiffV2QA){
    for (Int_t iPtBin = 1; iPtBin <= 100; ++iPtBin) {
      double pt_thisBin  = hRepnw1_thisEvt->GetYaxis()->GetBinCenter(iPtBin);
      for (Int_t iEtaBin = 1; iEtaBin <= 32; ++iEtaBin) {
        double eta_thisBin = hRepnw1_thisEvt->GetZaxis()->GetBinCenter(iEtaBin);
        for (Int_t iPOIBin = 1; iPOIBin <= NPOIBINS; ++iPOIBin) {
          double poi_thisBin = hRepnw1_thisEvt->GetXaxis()->GetBinCenter(iPOIBin);

          if (hMp_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin) < 1) continue;
          TComplex p(hRepnw1_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin), hImpnw1_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin));
          double mp = hMp_w1_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin);
  
          if (eta_thisBin > 0.) {
            pDiffFlowpQStarGap_thisEvt->Fill(poi_thisBin, pt_thisBin, eta_thisBin, (p*fNegEtaQStar).Re()/(mp*fNegEtaMQ),       (mp*fNegEtaMQ));
          } else {
            pDiffFlowpQStarGap_thisEvt->Fill(poi_thisBin, pt_thisBin, eta_thisBin, (p*fPosEtaQStar).Re()/(mp*fPosEtaMQ),       (mp*fPosEtaMQ));
          }
        }
      }
    }
  };

  if (fDoTrivalCor){ 

    // Q+
    double negEtaRePos = hReQPos_thisEvt->Integral(1, 100, 1, hReQPos_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
    double negEtaImPos = hImQPos_thisEvt->Integral(1, 100, 1, hImQPos_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
    if (hMQPos_thisEvt->Integral(1, 100, 1,  hMQPos_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6)) < 2) return false;
    double negEtaMQPos =  hMQPos_w1_thisEvt->Integral(1, 100, 1,  hMQPos_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
    TComplex negEtaQPos(negEtaRePos, negEtaImPos);
    TComplex negEtaQStarPos = TComplex::Conjugate(negEtaQPos);
  
    double posEtaRePos = hReQPos_thisEvt->Integral(1, 100, hReQPos_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
    double posEtaImPos = hImQPos_thisEvt->Integral(1, 100, hImQPos_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
    if (hMQPos_thisEvt->Integral(1, 100,  hMQPos_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32) < 2) return false;
    double posEtaMQPos =  hMQPos_w1_thisEvt->Integral(1, 100,  hMQPos_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
    TComplex posEtaQPos(posEtaRePos, posEtaImPos);
    TComplex posEtaQStarPos = TComplex::Conjugate(posEtaQPos);
  
    // Q-
    double negEtaReNeg = hReQNeg_thisEvt->Integral(1, 100, 1, hReQNeg_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
    double negEtaImNeg = hImQNeg_thisEvt->Integral(1, 100, 1, hImQNeg_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
    if (hMQNeg_thisEvt->Integral(1, 100, 1,  hMQNeg_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6)) < 2) return false;
    double negEtaMQNeg =  hMQNeg_w1_thisEvt->Integral(1, 100, 1,  hMQNeg_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
    TComplex negEtaQNeg(negEtaReNeg, negEtaImNeg);
    TComplex negEtaQStarNeg = TComplex::Conjugate(negEtaQNeg);
  
    double posEtaReNeg = hReQNeg_thisEvt->Integral(1, 100, hReQNeg_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
    double posEtaImNeg = hImQNeg_thisEvt->Integral(1, 100, hImQNeg_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
    if (hMQNeg_thisEvt->Integral(1, 100,  hMQNeg_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32) < 2) return false;
    double posEtaMQNeg =  hMQNeg_w1_thisEvt->Integral(1, 100,  hMQNeg_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
    TComplex posEtaQNeg(posEtaReNeg, posEtaImNeg);
    TComplex posEtaQStarNeg = TComplex::Conjugate(posEtaQNeg);
  
    pRefFlow_thisEvt->Fill(1.,  (fNegEtaQ*posEtaQStarPos).Re()    / (fNegEtaMQ*posEtaMQPos),    (fNegEtaMQ*posEtaMQPos));
    pRefFlow_thisEvt->Fill(2.,  (negEtaQPos*fPosEtaQStar).Re()    / (negEtaMQPos*fPosEtaMQ),    (negEtaMQPos*fPosEtaMQ));
    pRefFlow_thisEvt->Fill(3.,  (fNegEtaQ*posEtaQStarNeg).Re()    / (fNegEtaMQ*posEtaMQNeg),    (fNegEtaMQ*posEtaMQNeg));
    pRefFlow_thisEvt->Fill(4.,  (negEtaQNeg*fPosEtaQStar).Re()    / (negEtaMQNeg*fPosEtaMQ),    (negEtaMQNeg*fPosEtaMQ));
    pRefFlow_thisEvt->Fill(5.,  (negEtaQPos*posEtaQStarPos).Re() / (negEtaMQPos*posEtaMQPos), (negEtaMQPos*posEtaMQPos));
    pRefFlow_thisEvt->Fill(6.,  (negEtaQNeg*posEtaQStarNeg).Re() / (negEtaMQNeg*posEtaMQNeg), (negEtaMQNeg*posEtaMQNeg));
    pRefFlow_thisEvt->Fill(7.,  (negEtaQPos*posEtaQStarNeg).Re() / (negEtaMQPos*posEtaMQNeg), (negEtaMQPos*posEtaMQNeg));
    pRefFlow_thisEvt->Fill(8.,  (negEtaQNeg*posEtaQStarPos).Re() / (negEtaMQNeg*posEtaMQPos), (negEtaMQNeg*posEtaMQPos));

    // d{2} wt gap
    for (Int_t iPtBin = 1; iPtBin <= 100; ++iPtBin) {
      double pt_thisBin  = hRepnw1_thisEvt->GetYaxis()->GetBinCenter(iPtBin);
      for (Int_t iEtaBin = 1; iEtaBin <= 32; ++iEtaBin) {
        double eta_thisBin = hRepnw1_thisEvt->GetZaxis()->GetBinCenter(iEtaBin);
        for (Int_t iPOIBin = 1; iPOIBin <= NPOIBINS; ++iPOIBin) {
          double poi_thisBin = hRepnw1_thisEvt->GetXaxis()->GetBinCenter(iPOIBin);

          if (hMp_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin) < 1) continue;
          TComplex p(hRepnw1_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin), hImpnw1_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin));
          double mp = hMp_w1_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin);
  
          if (eta_thisBin > 0.) {
            pDiffFlowpQStarGapPos_thisEvt->Fill(poi_thisBin, pt_thisBin, eta_thisBin, (p*negEtaQStarPos).Re()/(mp*negEtaMQPos), (mp*negEtaMQPos));
            pDiffFlowpQStarGapNeg_thisEvt->Fill(poi_thisBin, pt_thisBin, eta_thisBin, (p*negEtaQStarNeg).Re()/(mp*negEtaMQNeg), (mp*negEtaMQNeg));
          } else {
            pDiffFlowpQStarGapPos_thisEvt->Fill(poi_thisBin, pt_thisBin, eta_thisBin, (p*posEtaQStarPos).Re()/(mp*posEtaMQPos), (mp*posEtaMQPos));
            pDiffFlowpQStarGapNeg_thisEvt->Fill(poi_thisBin, pt_thisBin, eta_thisBin, (p*posEtaQStarNeg).Re()/(mp*posEtaMQNeg), (mp*posEtaMQNeg));
          }
        }
      }
    }
  };

 return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::DoCumulantsDirect() 
{
  if (hMQ_thisEvt->Integral(1, 100, 1, 32) < 4) return false;

  double allReQnw1      = hReQn_thisEvt ->Integral(1, 100, 1, 32);
  double allImQnw1       = hImQn_thisEvt ->Integral(1, 100, 1, 32);
  double allReQ2nw2    = hReQ2nw2_thisEvt->Integral(1, 100, 1, 32);
  double allImQ2nw2     = hImQ2nw2_thisEvt->Integral(1, 100, 1, 32);
  double allReQnw3      = hReQnw3_thisEvt->Integral(1, 100, 1, 32);
  double allImQnw3      = hImQnw3_thisEvt->Integral(1, 100, 1, 32);
  TComplex Qnw1(allReQnw1, allImQnw1);
  TComplex Q2nw2(allReQ2nw2, allImQ2nw2);
  TComplex Qnw3(allReQnw3, allImQnw3);
  TComplex Qnw1Star = TComplex::Conjugate(Qnw1);
  TComplex Q2nw2Star = TComplex::Conjugate(Q2nw2);
  TComplex Qnw3Star = TComplex::Conjugate(Qnw3);

  double S11 = hMQ_w1_thisEvt->Integral(1,100, 1, 32);
  double S12 = hMQ_w2_thisEvt->Integral(1,100, 1, 32);
  double S13 = hMQ_w3_thisEvt->Integral(1,100, 1, 32);
  double S14 = hMQ_w4_thisEvt->Integral(1,100, 1, 32);
  double S21 = pow(hMQ_w1_thisEvt->Integral(1,100, 1, 32) ,2. );
  double S31 = pow(hMQ_w1_thisEvt->Integral(1,100, 1, 32) ,3. );
  double S41 = pow(hMQ_w1_thisEvt->Integral(1,100, 1, 32) ,4. );
  double S22 = pow(hMQ_w2_thisEvt->Integral(1,100, 1, 32) ,2. );

  // c{2} no gap
  pRefFlow_thisEvt->Fill(9.,  ((Qnw1*Qnw1Star).Re() - S12 ) / (S21 - S12), (S21 - S12)); 

  // c{4} 
  double M1111 = S41 - 6*S12*S21 + 8*S13*S11 + 3*S22 - 6*S14;
  double c4 = ((Qnw1*Qnw1*Qnw1Star*Qnw1Star).Re() + (Q2nw2*Q2nw2Star).Re() - 2*(Q2nw2*Qnw1Star*Qnw1Star).Re()
                        + 8* (Qnw3* Qnw1Star).Re() - 4* S12*(Qnw1*Qnw1Star).Re() - 6*S14 +2*S22 ) / M1111; //QC paper Eq.(18)
  pRefFlow_thisEvt->Fill(10.,  c4,  M1111);

  // d{4} 
  for (Int_t iPtBin = 1; iPtBin <= 100; ++iPtBin) {
    double pt_thisBin  = hRepnw1_thisEvt->GetYaxis()->GetBinCenter(iPtBin);
    for (Int_t iEtaBin = 1; iEtaBin <= 32; ++iEtaBin) {
      double eta_thisBin = hRepnw1_thisEvt->GetZaxis()->GetBinCenter(iEtaBin);
      for (Int_t iPOIBin = 1; iPOIBin <= NPOIBINS; ++iPOIBin) {
        double poi_thisBin = hRepnw1_thisEvt->GetXaxis()->GetBinCenter(iPOIBin);
        if (hMp_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin) < 1) continue;

        double s11 = hMp_w1_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin);
        double s12 = hMp_w2_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin);
        double s13 = hMp_w3_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin);
        double s14 = hMp_w4_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin);
        TComplex pnw1(hRepnw1_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin), hImpnw1_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin));
        TComplex pnw3(hRepnw3_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin), hImpnw3_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin));
        TComplex p2nw2(hRep2nw2_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin), hImp2nw2_thisEvt->GetBinContent(iPOIBin, iPtBin, iEtaBin));

        TComplex pnw1Star = TComplex::Conjugate(pnw1);
        TComplex pnw3Star = TComplex::Conjugate(pnw3);

        pDiffFlowpQStarDirect_thisEvt->Fill(poi_thisBin, pt_thisBin, eta_thisBin, ((pnw1*Qnw1Star).Re() - s12)/ (s11*S11 - s12),     (s11*S11 - s12));        

        double Mp1111 = s11* (S31 - 3*S11*S12 +2*S13) -  3*(s12*(S21 - S12) + 2*(s14 - s13*S11));

        double d4 = ((pnw1*Qnw1*Qnw1Star*Qnw1Star).Re() - (p2nw2*Qnw1Star*Qnw1Star).Re() - (pnw1*Qnw1*Q2nw2Star).Re()
         - 2*S12*(pnw1*Qnw1Star).Re() - 2*s12*(Qnw1*Qnw1Star).Re() + 7*(pnw3*Qnw1Star).Re()
         - (Qnw1*pnw3Star).Re() + (p2nw2*Q2nw2Star).Re() + 2*(pnw1*Qnw3Star).Re() + 2*s12*S12 - 6*s14) / Mp1111; 

        pDiffFlow4Direct_thisEvt->Fill(poi_thisBin, pt_thisBin, eta_thisBin, d4, Mp1111); 
      }
    }
  };
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::GetV0CalibHisto(AliAODEvent* fAOD, AliAODVertex* fVtx)
{
  if (fPeriod.EqualTo("LHC10h") ){
      for(int iCh = 0; iCh < 64; ++iCh) {
        fMultV0Ch[iCh] = hMultV0Read->GetBinContent(iCh+1, fRunNumBin+1);
      } 
      for(int i = 1; i < 3; ++i) {   // [0]: M; [1]: C; [2]: A;
        fV0XMean[i] = pV0XMeanRead[i]->GetBinContent(fRunNumBin+1, fCentBin+1, fVzBin+1);
        fV0YMean[i] = pV0YMeanRead[i]->GetBinContent(fRunNumBin+1, fCentBin+1, fVzBin+1);
      }    
      return true;
  } else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC11h")){ // A.Dobrin's V0Calib; "calibV0HIR.root"
      for(int iCh = 0; iCh < 64; ++iCh) {
        fMultV0Ch[iCh] = hMultV0[fRunNumBin]->GetBinContent(iCh+1);
      } 
      // AliCentrality* centrality = ((AliAODHeader*)fAOD->GetHeader())->GetCentralityP();
      // Double_t centSPD = centrality->GetCentralityPercentile("CL1");
      
      Int_t iCentSPD = -1;
      if (fPeriod.EqualTo("LHC15o")){
        AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
        if(!fMultSelection) {
          printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
          exit(1);
        }
        Double_t centCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
        iCentSPD = (Int_t)centCL1;
      }
      else if (fPeriod.EqualTo("LHC11h")){
        Double_t centCL1  =  fAOD->GetCentrality()->GetCentralityPercentile("CL1");
        iCentSPD = (Int_t)centCL1;
      }

      if (iCentSPD >= 90) return false;
      fV0XMean[0] = -999.; fV0YMean[0] = -999.; 
      for(int i = 0; i < 2; ++i) {   // [1]: C; [2]: A;
        fV0XMean[i+1] = hQxnmV0[fRunNumBin][i]->GetBinContent(iCentSPD+1);
        fV0YMean[i+1] = hQynmV0[fRunNumBin][i]->GetBinContent(iCentSPD+1);
      }   
  }
  else return false;
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::CalcQnVectorV0(AliAODEvent* fAOD, AliAODVertex* fVtx, double mAch, double Mult)
{
  if (!GetV0CalibHisto(fAOD, fVtx)) return false;
  double qxGE[3] = {0}, qyGE[3] = {0}, qxRecenter[3] = {0}, qyRecenter[3] = {0}; 
  double multRingGE[3] = {0};
  // [0]: M; [1]: C; [2]: A;
  for(int iCh = 0; iCh < 64; ++iCh) {
    if (TMath::IsNaN(fMultV0Ch[iCh]) || fMultV0Ch[iCh]<=0) continue;

    double phi = TMath::Pi()/8. + TMath::Pi()/4.*(iCh%8);
    double multCh = 0.;
    // double multCh = fAOD->GetVZEROEqMultiplicity(iCh);
    if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC11h")) multCh= fAOD->GetVZEROEqMultiplicity(iCh);
    else if (fPeriod.EqualTo("LHC15o")) {
      AliAODVZERO* aodV0 = fAOD->GetVZEROData();
      multCh = aodV0->GetMultiplicity(iCh);
    }
    if (iCh<32) { // C
      double multChGEC = -1;
      if (iCh < 8)
          multChGEC = multCh/fMultV0Ch[iCh] * fMultV0Ch[0];
      else if (iCh >= 8 && iCh < 16)
          multChGEC = multCh/fMultV0Ch[iCh] * fMultV0Ch[8];
      else if (iCh >= 16 && iCh < 24)
          multChGEC = multCh/fMultV0Ch[iCh] * fMultV0Ch[16];
      else if (iCh >= 24 && iCh < 32)
          multChGEC = multCh/fMultV0Ch[iCh] * fMultV0Ch[24];
      if (multChGEC<0 || TMath::IsNaN(multChGEC)) continue;

      qxGE[1] += multChGEC*TMath::Cos(fHarmonic*phi);
      qyGE[1] += multChGEC*TMath::Sin(fHarmonic*phi);
      multRingGE[1] += multChGEC;
      qxGE[0] += multChGEC*TMath::Cos(fHarmonic*phi);
      qyGE[0] += multChGEC*TMath::Sin(fHarmonic*phi); 
      multRingGE[0] += multChGEC;        
    } else if (iCh>=32 && iCh<64) { // A
      double multChGEA = -1;
      if (iCh >= 32 && iCh < 40)
          multChGEA = multCh/fMultV0Ch[iCh] * fMultV0Ch[32];
      else if (iCh >= 40 && iCh < 48)
          multChGEA = multCh/fMultV0Ch[iCh] * fMultV0Ch[40];
      else if (iCh >= 48 && iCh < 56)
          multChGEA = multCh/fMultV0Ch[iCh] * fMultV0Ch[48];
      else if (iCh >= 56 && iCh < 64)
          multChGEA = multCh/fMultV0Ch[iCh] * fMultV0Ch[56];
      if (multChGEA<0 || TMath::IsNaN(multChGEA)) continue;

      qxGE[2] += multChGEA*TMath::Cos(fHarmonic*phi);
      qyGE[2] += multChGEA*TMath::Sin(fHarmonic*phi);
      multRingGE[2] += multChGEA;  
      qxGE[0] += multChGEA*TMath::Cos(fHarmonic*phi);
      qyGE[0] += multChGEA*TMath::Sin(fHarmonic*phi); 
      multRingGE[0] += multChGEA;
    }
  }
  if (multRingGE[0] < 1e-6 || multRingGE[1] < 1e-6) return false;
  double Qn_thisEvt[3];
  for (int i = 1; i < 3; ++i) {
      Qn_thisEvt[i] = -1;
      double qxMean = fV0XMean[i];
      double qyMean = fV0YMean[i];  
      if (TMath::IsNaN(qxMean) || TMath::IsNaN(qyMean) ) continue;
      if (qyMean < -900 || qxMean < -900) continue; 
      // For 10 h, we've stored the qx/y of V0M, and they cannot been found in A.Dorbin's calib file for 15o period!
      qxRecenter[i] = qxGE[i] - qxMean;
      qyRecenter[i] = qyGE[i] - qyMean;
      Qn_thisEvt[i] = sqrt(qxRecenter[i]*qxRecenter[i] + qyRecenter[i]*qyRecenter[i])/ sqrt(multRingGE[i]);
      hQnCentRecenter[i]->Fill(fCent, Qn_thisEvt[i]);
      // psiRecenter
      double psiRecenter = GetEventPlane(qxRecenter[i], qyRecenter[i]);
      if (psiRecenter>TMath::Pi())  psiRecenter -=TMath::Pi();
      hPsiV0Recenter[fCentBin][i]->Fill(psiRecenter); 
      if (fQAV0){
        double vz    = fVtx->GetZ();  
        double centSPD  =0.;
        if (fPeriod.EqualTo("LHC15o")){
          AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
          centSPD = fMultSelection->GetMultiplicityPercentile("CL1");
        }
        else if (fPeriod.EqualTo("LHC11h") || fPeriod.EqualTo("LHC10h") ){
          centSPD =  fAOD->GetCentrality()->GetCentralityPercentile("CL1");
        }
        hQxCentRecenter[i]->Fill(centSPD, qxRecenter[i]);
        hQyCentRecenter[i]->Fill(centSPD, qyRecenter[i]);
        hQxVtxRecenter[i]->Fill(vz, qxRecenter[i]);
        hQyVtxRecenter[i]->Fill(vz, qyRecenter[i]);
      }   
  }

  if (Qn_thisEvt[1]<0) return false; // V0C
  fQnBin = GetQnPercV0(Qn_thisEvt[1], mAch, Mult); 

  if(fQnBin >=0) return true;
  else if (fQnBin < 0)  {Printf("Wrong QnBin separation!"); return false;}

  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCMWESE::CalcQnVectorTPC(AliAODEvent* fAOD)
{ // return true;
}

//---------------------------------------------------

int AliAnalysisTaskCMWESE::GetQnPercV0(double qn_thisEvt, double mAch, double Mult)
{
  if (fPeriod.EqualTo("LHC10h")){
    int NbinX = hQnPercentile->GetNbinsX();
    hQnPercentile->GetXaxis()->SetRange((int)fCent+1,(int)fCent+1);  
    hQnPercentile_centThisEvt = (TH1D*)hQnPercentile->ProjectionY("qnPercentileCentiEvt", 0, -1, "e");
    sp = new TSpline3(hQnPercentile_centThisEvt);
    double pec = 100* sp->Eval(qn_thisEvt);
    int QnBin = GetPercCode(pec);
    hAch[0][QnBin]->Fill(mAch);
    hMult[0][QnBin]->Fill(Mult); //mb

    hAch[fCentBin+1][QnBin]->Fill(mAch);
    hMult[fCentBin+1][QnBin]->Fill(Mult);  
    
    hQnPercentile->GetXaxis()->SetRange(1,NbinX); 
    return QnBin;   
  } 

  else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC11h")){
    double pec = 100.* splQ2c[(int)fCent]->Eval(qn_thisEvt);
    int QnBin = GetPercCode(pec);
    hAch[0][QnBin]->Fill(mAch);
    hMult[0][QnBin]->Fill(Mult); //mb

    hAch[fCentBin+1][QnBin]->Fill(mAch);
    hMult[fCentBin+1][QnBin]->Fill(Mult);
    
    return QnBin;
  }
  return -1;
}

//---------------------------------------------------

int AliAnalysisTaskCMWESE::GetPercCode(double perc)
{
    int percCode = -1;

    if ((perc >= 0) && (perc <= 10.0))
        percCode = 0;
    else if ((perc > 10.0) && (perc <= 20.0))
        percCode = 1;
    else if ((perc > 20.0) && (perc <= 30.0))
        percCode = 2;
    else if ((perc > 30.0) && (perc <= 40.0))
        percCode = 3;
    else if ((perc > 40.0) && (perc <= 50.0))
        percCode = 4;
    else if ((perc > 50.0) && (perc <= 60.0))
        percCode = 5;
    else if ((perc > 60.0) && (perc <= 70.0))
        percCode = 6;
    else if ((perc > 70.0) && (perc <= 80.0))
        percCode = 7;
    else if ((perc > 80.0) && (perc <= 90.0))
        percCode = 8;
    else if (perc > 90.0)
        percCode = 9;  

    return percCode;

}
