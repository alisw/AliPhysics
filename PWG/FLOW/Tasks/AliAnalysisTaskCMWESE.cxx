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
  fTrigger("kMB"),
  fCentEst("V0M"),
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
  fMultComp("none"),
  fEtaGap(0.3),  
  fV0CalibOn(true),
  fTPCCalibOn(false),
  fQAV0(true),
  fQATPC(false),
  fDoNUE(true),
  fDoNUA(true),
  fCentCut(7.5),
  fRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCentBin(-999),
  fCent(-999),
  fQnBin(-999),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fZvtxCut(10.0),
  fPUSyst(0),
  fListNUE(NULL),
  fListNUA(NULL), 
  fListVZEROCALIB(NULL),
  fListwCent(NULL),
  hwCent(NULL),
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
  fHCorrectV0ChWeghts(NULL),
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
  hReQ_thisEvt(NULL),
  hImQ_thisEvt(NULL),
  hReQ2_thisEvt(NULL), 
  hImQ2_thisEvt(NULL), 
  hMQ_thisEvt(NULL),
  hMQ_weight_thisEvt(NULL),
  hReQPos_thisEvt(NULL),
  hImQPos_thisEvt(NULL),
  hMQPos_thisEvt(NULL),
  hMQPos_weight_thisEvt(NULL),
  hReQNeg_thisEvt(NULL),
  hImQNeg_thisEvt(NULL),
  hMQNeg_thisEvt(NULL),
  hMQNeg_weight_thisEvt(NULL),
  pRefFlow_thisEvt(NULL), 
  pIntd2_thisEvt(NULL),
  pbufferPt_thisEvt(NULL),
  pbufferEta_thisEvt(NULL)
{

  for (int i = 0; i < 3; ++i)   pV0XMeanRead[i]=NULL; 
  for (int i = 0; i < 3; ++i)   pV0YMeanRead[i]=NULL;

  hMultV0=NULL; //Dobrin
  for (int i = 0; i < 2; ++i)   hQxnmV0[i]=NULL;
  for (int i = 0; i < 2; ++i)   hQynmV0[i]=NULL;

  for (int i = 0; i < 64; ++i) fMultV0Ch[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0XMean[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0YMean[i]=-999.;

  for (int i = 0; i < 90; ++i) splQ2c[i]=NULL; 

  for (int i = 0; i < 2; ++i) hCent[i]=NULL;
  for (int i = 0; i < 138; ++i) hCentRBR[i]=NULL; 
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

  for (int i = 0; i < NQNBINS; ++i) hMultMb[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) hAchMb[i]=NULL;

  for (int i = 0; i < NQNBINS; ++i) pPosPtAch[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) pNegPtAch[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) pPosEtaAch[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) pNegEtaAch[i]=NULL;
}

//---------------------------------------------------

AliAnalysisTaskCMWESE::AliAnalysisTaskCMWESE(const char *name, TString _PR, bool _NUE, bool _NUA, bool _V0Calib) : 
  AliAnalysisTaskSE(name), 
  fDebug(0),
  fHarmonic(2.),
  fTrigger("kMB"),
  fCentEst("V0M"),
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
  fMultComp("none"),
  fEtaGap(0.3),  
  fV0CalibOn(_V0Calib),
  fTPCCalibOn(false),
  fQAV0(true),
  fQATPC(false),
  fDoNUE(_NUE),
  fDoNUA(_NUA),
  fCentCut(7.5),
  fRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCentBin(-999),
  fCent(-999),
  fQnBin(-999),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fZvtxCut(10.0),
  fPUSyst(0),
  fListNUE(NULL),
  fListNUA(NULL), 
  fListVZEROCALIB(NULL),
  fListwCent(NULL),
  hwCent(NULL),
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
  fHCorrectV0ChWeghts(NULL),
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
  hReQ_thisEvt(NULL),
  hImQ_thisEvt(NULL),
  hReQ2_thisEvt(NULL), 
  hImQ2_thisEvt(NULL), 
  hMQ_thisEvt(NULL),
  hMQ_weight_thisEvt(NULL),
  hReQPos_thisEvt(NULL),
  hImQPos_thisEvt(NULL),
  hMQPos_thisEvt(NULL),
  hMQPos_weight_thisEvt(NULL),
  hReQNeg_thisEvt(NULL),
  hImQNeg_thisEvt(NULL),
  hMQNeg_thisEvt(NULL),
  hMQNeg_weight_thisEvt(NULL),
  pRefFlow_thisEvt(NULL), 
  pIntd2_thisEvt(NULL),
  pbufferPt_thisEvt(NULL),
  pbufferEta_thisEvt(NULL)
{

  for (int i = 0; i < 3; ++i)   pV0XMeanRead[i]=NULL; 
  for (int i = 0; i < 3; ++i)   pV0YMeanRead[i]=NULL;

  hMultV0=NULL; //Dobrin
  for (int i = 0; i < 2; ++i)   hQxnmV0[i]=NULL;
  for (int i = 0; i < 2; ++i)   hQynmV0[i]=NULL;

  for (int i = 0; i < 64; ++i) fMultV0Ch[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0XMean[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0YMean[i]=-999.;

  for (int i = 0; i < 90; ++i) splQ2c[i]=NULL; 

  for (int i = 0; i < 2; ++i) hCent[i]=NULL;
  for (int i = 0; i < 138; ++i) hCentRBR[i]=NULL; 
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

  for (int i = 0; i < NQNBINS; ++i) hMultMb[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) hAchMb[i]=NULL;

  for (int i = 0; i < NQNBINS; ++i) pPosPtAch[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) pNegPtAch[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) pPosEtaAch[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) pNegEtaAch[i]=NULL;

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
  } 
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
      DefineInput(inputslot, TList::Class());
      inputslot++;   
    } 
  }
  if (fPeriod.EqualTo("LHC18q")){
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
  }
  if (fPeriod.EqualTo("LHC18r")){
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
  }

}

//------------------------------------------------

AliAnalysisTaskCMWESE::AliAnalysisTaskCMWESE(const char *name) : 
  AliAnalysisTaskSE(name), 
  fDebug(0),
  fHarmonic(2.),
  fTrigger("kMB"),
  fCentEst("V0M"),
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
  fMultComp("none"),
  fEtaGap(0.3),  
  fV0CalibOn(true),
  fTPCCalibOn(false),
  fQAV0(true),
  fQATPC(false),
  fDoNUE(true),
  fDoNUA(true),
  fCentCut(7.5),
  fRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCentBin(-999),
  fCent(-999),
  fQnBin(-999),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fZvtxCut(10.0),
  fPUSyst(0),
  fListNUE(NULL),
  fListNUA(NULL), 
  fListVZEROCALIB(NULL),
  fListwCent(NULL),
  hwCent(NULL),
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
  fHCorrectV0ChWeghts(NULL),
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
  hReQ_thisEvt(NULL),
  hImQ_thisEvt(NULL),
  hReQ2_thisEvt(NULL), 
  hImQ2_thisEvt(NULL), 
  hMQ_thisEvt(NULL),
  hMQ_weight_thisEvt(NULL),
  hReQPos_thisEvt(NULL),
  hImQPos_thisEvt(NULL),
  hMQPos_thisEvt(NULL),
  hMQPos_weight_thisEvt(NULL),
  hReQNeg_thisEvt(NULL),
  hImQNeg_thisEvt(NULL),
  hMQNeg_thisEvt(NULL),
  hMQNeg_weight_thisEvt(NULL),
  pRefFlow_thisEvt(NULL), 
  pIntd2_thisEvt(NULL),
  pbufferPt_thisEvt(NULL),
  pbufferEta_thisEvt(NULL)
{

  for (int i = 0; i < 3; ++i)   pV0XMeanRead[i]=NULL; 
  for (int i = 0; i < 3; ++i)   pV0YMeanRead[i]=NULL;

  hMultV0=NULL; //Dobrin
  for (int i = 0; i < 2; ++i)   hQxnmV0[i]=NULL;
  for (int i = 0; i < 2; ++i)   hQynmV0[i]=NULL;

  for (int i = 0; i < 64; ++i) fMultV0Ch[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0XMean[i]=-999.;
  for (int i = 0; i < 3; ++i)   fV0YMean[i]=-999.;

  for (int i = 0; i < 90; ++i) splQ2c[i]=NULL; 

  for (int i = 0; i < 2; ++i) hCent[i]=NULL;
  for (int i = 0; i < 138; ++i) hCentRBR[i]=NULL; 
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

  for (int i = 0; i < NQNBINS; ++i) hMultMb[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) hAchMb[i]=NULL;

  for (int i = 0; i < NQNBINS; ++i) pPosPtAch[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) pNegPtAch[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) pPosEtaAch[i]=NULL;
  for (int i = 0; i < NQNBINS; ++i) pNegEtaAch[i]=NULL;

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

  //------------------
  // Info
  // ------------------
  if (fPeriod.EqualTo("LHC10h") ){
      TString runNumList[90]={
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
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,runNumList[i].Data());
      }
      fOutputList->Add(hRunNumBin);
  } else if (fPeriod.EqualTo("LHC11h") ){
      TString runNumList[39]={
        "170387","170040","170268","170228","170207","169838","170159","170204","170311","170084",
        "169835","170088","170593","170203","170270","169846","170163","170388","170155","170083",
        "170572","169837","169855","170306","170269","170089","170309","170091","170081","170230",
        "170085","170315","170027","170193","170312","170313","170308","169858","169859"};
      hRunNumBin = new TH1I("runNumBin","",100,0,100);
      for (int i=0; i<39; ++i){    
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,runNumList[i].Data());
      }
      fOutputList->Add(hRunNumBin);
  } else if (fPeriod.EqualTo("LHC15o") ){
      TString runNumList[138]={
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
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,runNumList[i].Data());
        if (fPUSyst==3){
          hCentRBR[i] = new TH1D(Form("hCentDist%i", (int)runNumList[i].Atoi()), "", 100, 0., 100.);
          fOutputList->Add(hCentRBR[i]);
        }
      }
      fOutputList->Add(hRunNumBin);
  } else if (fPeriod.EqualTo("LHC18q") ){
      TString runNumList[125]={
        "296623", "296622", "296621", "296619", "296618", "296616", "296615", "296594", "296553", "296552", "296551", "296550",
        "296548", "296547", "296516", "296512", "296511", "296510", "296509", "296472", "296433", "296424", "296423", "296420",
        "296419", "296415", "296414", "296383", "296381", "296380", "296379", "296378", "296377", "296376", "296375", "296312",
        "296309", "296304", "296303", "296280", "296279", "296273", "296270", "296269", "296247", "296246", "296244", "296243",
        "296242", "296241", "296240", "296198", "296197", "296196", "296195", "296194", "296192", "296191", "296143", "296142",
        "296135", "296134", "296133", "296132", "296123", "296074", "296066", "296065", "296063", "296062", "296060", "296016",
        "295942", "295941", "295937", "295936", "295913", "295910", "295909", "295861", "295860", "295859", "295856", "295855",
        "295854", "295853", "295831", "295829", "295826", "295825", "295822", "295819", "295818", "295816", "295791", "295788",
        "295786", "295763", "295762", "295759", "295758", "295755", "295754", "295725", "295723", "295721", "295719", "295718",
        "295717", "295714", "295712", "295676", "295675", "295673", "295668", "295667", "295666", "295615", "295612", "295611",
        "295610", "295589", "295588", "295586", "295585"
      };
      hRunNumBin = new TH1I("runNumBin","",150,0,150);
      for (int i=0; i<125; ++i){    
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,runNumList[i].Data());
      }
      fOutputList->Add(hRunNumBin);
  }  else if (fPeriod.EqualTo("LHC18r") ){
      TString runNumList[89]={
        "297595", "297590", "297588", "297558", "297544", "297542", "297541", "297540", "297537", "297512", "297483", "297479",
        "297452", "297451", "297450", "297446", "297442", "297441", "297415", "297414", "297413", "297406", "297405", "297380",
        "297379", "297372", "297367", "297366", "297363", "297336", "297335", "297333", "297332", "297317", "297311", "297310",
        "297278", "297222", "297221", "297218", "297196", "297195", "297193", "297133", "297132", "297129", "297128", "297124",
        "297123", "297119", "297118", "297117", "297085", "297035", "297031", "296966", "296941", "296938", "296935", "296934",
        "296932", "296931", "296930", "296903", "296900", "296899", "296894", "296852", "296851", "296850", "296848", "296839",
        "296838", "296836", "296835", "296799", "296794", "296793", "296790", "296787", "296786", "296785", "296784", "296781",
        "296752", "296694", "296693", "296691", "296690"
      };
      hRunNumBin = new TH1I("runNumBin","",150,0,150);
      for (int i=0; i<89; ++i){    
        hRunNumBin->GetXaxis()->SetBinLabel(i+1,runNumList[i].Data());
      }
      fOutputList->Add(hRunNumBin);
  }

  // Copy TList
  Int_t inSlotCounter=1;
  if(fDoNUE) {
    if (fPeriod.EqualTo("LHC11h") || fPeriod.EqualTo("LHC10h") ) {
      fListNUE = (TList*) GetInputData(inSlotCounter);
      if (fDebug) fListNUE->ls();
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
    else if (fPeriod.EqualTo("LHC18q")){
      fListNUE = (TList*) GetInputData(inSlotCounter);
      inSlotCounter++;
      if (!fListNUE) {
        AliError(Form("%s: Wrong NUE file 15o.", GetName()));
        return;
      };     
    }
    else if (fPeriod.EqualTo("LHC18r")){
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
      fListNUA = (TList*) GetInputData(inSlotCounter);
      if (fDebug) fListNUA->ls();
      inSlotCounter++;
      if (!fListNUA) {
        AliError(Form("%s: Wrong NUA file 10h.", GetName()));
        return;
      };
    } else if (fPeriod.EqualTo("LHC11h") ){
      fListNUA = (TList*) GetInputData(inSlotCounter);
      inSlotCounter++;
      if (!fListNUA) {        
        AliError(Form("%s: Wrong NUA file 11h.", GetName()));
        return;
      }
    } else if (fPeriod.EqualTo("LHC15o") ){ //Rihan's NUA
      fListNUA = (TList*) GetInputData(inSlotCounter);
      if (fDebug) fListNUA->ls();
      inSlotCounter++;
      if (!fListNUA) {        
        AliError(Form("%s: Wrong NUA file 15op2.", GetName()));
        return;
      }
    }  else if (fPeriod.EqualTo("LHC18q") ){ //Rihan's NUA
      fListNUA = (TList*) GetInputData(inSlotCounter);
      if (fDebug) fListNUA->ls();
      inSlotCounter++;
      if (!fListNUA) {        
        AliError(Form("%s: Wrong NUA file 18qp3.", GetName()));
        return;
      }
    }  else if (fPeriod.EqualTo("LHC18r") ){ //Rihan's NUA
      fListNUA = (TList*) GetInputData(inSlotCounter);
      if (fDebug) fListNUA->ls();
      inSlotCounter++;
      if (!fListNUA) {        
        AliError(Form("%s: Wrong NUA file 18rp3.", GetName()));
        return;
      }
    } 
  }

  // Qn Calib
  if (fV0CalibOn){
      if (fPeriod.EqualTo("LHC10h") ) {
            fListVZEROCALIB = (TList*) GetInputData(inSlotCounter);
            if (fDebug) fListVZEROCALIB->ls();
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
            hQnPercentile = (TH2D*)fListVZEROCALIB->FindObject("h_qncPercentile");

      } else if (fPeriod.EqualTo("LHC15o") ) {
            fListVZEROCALIB = (TList*) GetInputData(inSlotCounter);
            if (fDebug) fListVZEROCALIB->ls();
            inSlotCounter++;
            if (!fListVZEROCALIB) {
             AliError(Form("%s: Spline list for 15o cannot be opened.", GetName()));
              return;
            };
            //TSplines in 1% centralicities
            for (Int_t isp = 0; isp < 90; isp++) splQ2c[isp] = (TSpline3*)fListVZEROCALIB->FindObject(Form("sp_q2V0C_%d", isp));

            fListwCent = (TList*) GetInputData(inSlotCounter);   
            inSlotCounter++;   
            if (!fListwCent) {
             AliError(Form("%s: weight of centrality for 15o cannot be opened.", GetName()));
              return;
            };  
            hwCent = (TH1D*)fListwCent->FindObject("weightCent");
      } else if (fPeriod.EqualTo("LHC18q") ) {
            fListVZEROCALIB = (TList*) GetInputData(inSlotCounter);
            if (fDebug) fListVZEROCALIB->ls();
            inSlotCounter++;
            if (!fListVZEROCALIB) {
             AliError(Form("%s: Spline list for 18q cannot be opened.", GetName()));
              return;
            };
            //TSplines in 1% centralicities
            hQnPercentile = (TH2D*)fListVZEROCALIB->FindObject("h_qncPercentile");

      } else if (fPeriod.EqualTo("LHC18r") ) {
            fListVZEROCALIB = (TList*) GetInputData(inSlotCounter);
            if (fDebug) fListVZEROCALIB->ls();
            inSlotCounter++;
            if (!fListVZEROCALIB) {
             AliError(Form("%s: Spline list for 18r cannot be opened.", GetName()));
              return;
            };
            //TSplines in 1% centralicities
            hQnPercentile = (TH2D*)fListVZEROCALIB->FindObject("h_qncPercentile");
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
  hEvtCount->GetXaxis()->SetBinLabel(15,"ptetaAch");
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

  hMultCentQA[0] = new TH2D("hMultCentQABefCut", ";centV0M;multFB32", 100, 0, 100, 200, 0, 5000);
  hMultCentQA[1] = new TH2D("hMultCentQAAftCut", ";centV0M;multFB32", 100, 0, 100, 200, 0, 5000);
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

  for (int i = 0; i < NQNBINS; ++i){
    hMultMb[i] = new TH1D(Form("hMultMbQn%i", i ), "", 100, 0., 5000.);
    hAchMb[i] = new TH1D(Form("hAchMbQn%i", i ), "", 200, -1., 1.);                
  }
  for (int i = 0; i < NQNBINS; ++i) fOutputList->Add(hMultMb[i]); 
  for (int i = 0; i < NQNBINS; ++i) fOutputList->Add(hAchMb[i]); 

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
  hDcaXy[0] = new TH1D("hDcaXyBeforeCut", "", 100, 0., 10.);
  hDcaXy[1] = new TH1D("hDcaXyAfterCut",  "", 100, 0., 10.);
  fOutputList->Add(hDcaXy[0]);
  fOutputList->Add(hDcaXy[1]);
  hDcaZ[0] = new TH1D("hDcaZBeforeCut", "", 100, 0., 10.);
  hDcaZ[1] = new TH1D("hDcaZAfterCut",  "", 100, 0., 10.);
  fOutputList->Add(hDcaZ[0]);
  fOutputList->Add(hDcaZ[1]);

  if (fQAV0){
    hQnCentRecenter[0] = new TH2D("hqnV0MRecenter","",100, 0, 100, 570, 0, 12); 
    hQnCentRecenter[1] = new TH2D("hqnV0CRecenter","",100, 0, 100, 570, 0, 12); 
    hQnCentRecenter[2] = new TH2D("hqnV0ARecenter","",100, 0, 100, 570, 0, 12); 
    // fOutputList->Add(hQnCentRecenter[0]);
    fOutputList->Add(hQnCentRecenter[1]);
    fOutputList->Add(hQnCentRecenter[2]); 
    hQxCentRecenter[0] = new TH2D("hqxRecenterV0MvsCentSPD","",100, 0, 100, 32, -600, 600);
    hQxCentRecenter[1] = new TH2D("hqxRecenterV0CvsCentSPD","",100, 0, 100, 32, -600, 600);
    hQxCentRecenter[2] = new TH2D("hqxRecenterV0AvsCentSPD","",100, 0, 100, 32, -600, 600);
    hQyCentRecenter[0] = new TH2D("hqyRecenterV0MvsCentSPD","",100, 0, 100, 32, -600, 600);
    hQyCentRecenter[1] = new TH2D("hqyRecenterV0CvsCentSPD","",100, 0, 100, 32, -600, 600);
    hQyCentRecenter[2] = new TH2D("hqyRecenterV0AvsCentSPD","",100, 0, 100, 32, -600, 600);
    for (int i = 1; i < 3; ++i) fOutputList->Add(hQxCentRecenter[i]);
    for (int i = 1; i < 3; ++i) fOutputList->Add(hQyCentRecenter[i]); 
    hQxVtxRecenter[0] = new TH2D("hqxRecenterV0MvsVz","",20, -10, 10, 32, -600, 600);
    hQxVtxRecenter[1] = new TH2D("hqxRecenterV0CvsVz","",20, -10, 10, 32, -600, 600);
    hQxVtxRecenter[2] = new TH2D("hqxRecenterV0AvsVz","",20, -10, 10, 32, -600, 600);
    hQyVtxRecenter[0] = new TH2D("hqyRecenterV0MvsVz","",20, -10, 10, 32, -600, 600);
    hQyVtxRecenter[1] = new TH2D("hqyRecenterV0CvsVz","",20, -10, 10, 32, -600, 600);
    hQyVtxRecenter[2] = new TH2D("hqyRecenterV0AvsVz","",20, -10, 10, 32, -600, 600);
    for (int i = 1; i < 3; ++i) fOutputList->Add(hQxVtxRecenter[i]); 
    for (int i = 1; i < 3; ++i) fOutputList->Add(hQyVtxRecenter[i]);
    for (int iCent = 0; iCent < NCENTBINS; ++iCent){
      hPsiV0Recenter[iCent][0] = new TH1D(Form("hpsiV0MRecenterCent%i", iCent),"", 360, 0., 2*TMath::Pi());
      hPsiV0Recenter[iCent][1] = new TH1D(Form("hpsiV0CRecenterCent%i", iCent),"", 360, 0., 2*TMath::Pi());
      hPsiV0Recenter[iCent][2] = new TH1D(Form("hpsiV0ARecenterCent%i", iCent),"", 360, 0., 2*TMath::Pi());
    } 
    for (int iCent = fCbinLo; iCent < fCbinHg; ++iCent){
      for (int i = 1; i < 3; ++i) fOutputList->Add(hPsiV0Recenter[iCent][i]);    
    } 
  }
  
  //------------------
  // physics
  //------------------
  TH1::SetDefaultSumw2(kTRUE);

  hReQ_thisEvt             = new TH2D("ReQthisEvt", "", 100, 0., 5., 32, -0.8, 0.8);
  hImQ_thisEvt              = new TH2D("ImQthisEvt", "", 100, 0., 5., 32, -0.8, 0.8);
  hReQ2_thisEvt           = new TH2D("ReQ2thisEvt", "", 100, 0., 5., 32, -0.8, 0.8);
  hImQ2_thisEvt            = new TH2D("ImQ2thisEvt", "", 100, 0., 5., 32, -0.8, 0.8);
  hMQ_thisEvt               = new TH2D("mQthisEvt",  "", 100, 0., 5., 32, -0.8, 0.8);  
  hMQ_weight_thisEvt  = new TH2D("mQ_weightthisEvt",  "", 100, 0., 5., 32, -0.8, 0.8);  
  hReQPos_thisEvt       = new TH2D("ReQPosthisEvt", "", 100, 0., 5., 32, -0.8, 0.8);
  hImQPos_thisEvt        = new TH2D("ImQPosthisEvt", "", 100, 0., 5., 32, -0.8, 0.8);
  hMQPos_thisEvt        = new TH2D("mQPosthisEvt",  "", 100, 0., 5., 32, -0.8, 0.8);
  hMQPos_weight_thisEvt = new TH2D("mQPosWeightthisEvt",  "", 100, 0., 5., 32, -0.8, 0.8);
  hReQNeg_thisEvt     = new TH2D("ReQNegthisEvt", "", 100, 0., 5., 32, -0.8, 0.8);
  hImQNeg_thisEvt      = new TH2D("ImQNegthisEvt", "", 100, 0., 5., 32, -0.8, 0.8);
  hMQNeg_thisEvt      = new TH2D("mQNegthisEvt",  "", 100, 0., 5., 32, -0.8, 0.8);  
  hMQNeg_weight_thisEvt  = new TH2D("mQNegWeightthisEvt",  "", 100, 0., 5., 32, -0.8, 0.8);  
  pRefFlow_thisEvt      = new TProfile("pRefFlowthisEvt", "", 10, 0, 10, "s");  // as function of 
  pIntd2_thisEvt            = new TProfile("pIntd2thisEvt", "", 30, 0, 30, "s");

  for (int i = 0; i < NCENTBINS; ++i) {
    pRefFlow[i]               = new TProfile(Form("pRefFlowCent%i",i), "", 10, 0, 10, "s");  // as function of delta eta
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

  // pt,eta Vs. Ach
  pbufferPt_thisEvt               = new TProfile("pbufferPt_thisEvt", "", 10, 0, 10, "s");
  pbufferEta_thisEvt             = new TProfile("pbufferEta_thisEvt", "", 10, 0, 10, "s");
  for (int i = 0; i < NCENTBINS; ++i) {
    pPosPtAch[i]       = new TProfile(Form("pPosPtAch_cent%i", i), "", 31, -1, 1, "s");
    pNegPtAch[i]      = new TProfile(Form("pNegPtAch_cent%i", i), "", 31, -1, 1, "s");    
    pPosEtaAch[i]     = new TProfile(Form("pPosEtaAch_cent%i", i), "", 31, -1, 1, "s");
    pNegEtaAch[i]     = new TProfile(Form("pNegEtaAch_cent%i", i), "", 31, -1, 1, "s");  
  }
  for (int i = fCbinLo; i < fCbinHg; ++i) {
    fOutputList->Add(pPosPtAch[i]);
    fOutputList->Add(pNegPtAch[i]);  
    fOutputList->Add(pPosEtaAch[i]);
    fOutputList->Add(pNegEtaAch[i]);  
  }

  // Dobrin 15o pass2 V0 Calib
  if (fPeriod.EqualTo("LHC15o") && (fPUSyst!=4) ){
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
  }

  // Dobrin 15o pass2 V0 Calib; new Para
  if (fPeriod.EqualTo("LHC15o")&& fPUSyst==4){
   fSPDCutPU = new TF1("fSPDCutPU", "450. + 3.9*x", 0, 50000);
   Double_t parV0[8] = {33.4127, 0.95351, 0.0716651, 227.19, 8.78312, -0.000549686, 0.000275397, -6.68358e-07};

    fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);

    Double_t parV0CL0[6] = {0.0224258, 0.975612, 0.676298, 0.0292538, -0.000550616, 5.87181e-06};
    fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 6.*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);

    Double_t parFB32[9] = {-814.463, 6.43163, 5423.19, -0.382271, 0.0299549, -26.5402, 321.52, -0.8256, 0.0168231};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCutPU->SetParameters(parFB32);
  }

  // if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){
  //   // Rihan Pile-up function 18
  //   fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

  //   Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  //   fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  //   fV0CutPU->SetParameters(parV0);
    
  //   Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  //   fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  //   fCenCutLowPU->SetParameters(parV0CL0);
  //   fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  //   fCenCutHighPU->SetParameters(parV0CL0);

  //   Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  //   fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  //   fMultCutPU->SetParameters(parFB32);
  // }

  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){
    // ADobirn PU for 18
    fSPDCutPU = new TF1("fSPDCutPU", "480. + 3.95*x", 0, 50000);
    
    Double_t parV0[8] = {42.4921, 0.823255, 0.0824939, 139.826, 7.27032, 0.0488425, -0.00045769, 1.40891e-06};
    fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);
    
    Double_t parV0CL0[6] = {0.317973, 0.961823, 1.02383, 0.0330231, -0.000721551, 6.92564e-06};
    fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);
    
    Double_t parFB32[9] = {-817.169, 6.40836, 5380.3, -0.394358, 0.0295209, -25.9573, 316.586, -0.843951, 0.0165442};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCutPU->SetParameters(parFB32);
  }

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
        if (fTrigger.EqualTo("kINT7")) { // 15o
      isTrigselected = mask&AliVEvent::kINT7;
        } else if (fTrigger.EqualTo("kMB")) {
      isTrigselected = mask&AliVEvent::kMB; //10h
        } else if (fTrigger.EqualTo("kMB+kCentral+kSemiCentral")) {
      isTrigselected = mask&AliVEvent::kMB+AliVEvent::kCentral+AliVEvent::kSemiCentral; //11h
        } else if (fTrigger.EqualTo("kINT7+kCentral+kSemiCentral")) {
      isTrigselected = mask&AliVEvent::kINT7+AliVEvent::kCentral+AliVEvent::kSemiCentral; //18q/r
        }
  if(isTrigselected == false) return;
  hEvtCount->Fill(3);
  if (fDebug) Printf("trigger done!");

  //----------------------------
  // event-wise selection
  //----------------------------
  // run number
  Int_t run = fAOD->GetRunNumber();
  if (run < 200000 && fPeriod.EqualTo("LHC15o")) return;
  if (run > 200000 && fPeriod.EqualTo("LHC11h")) return;
  if (run < 290000 && fPeriod.EqualTo("LHC18q")) return;
  if (run < 290000 && fPeriod.EqualTo("LHC18r")) return;
  if(run != fRunNum){
      // Load the calibrations run dependent
      if (fPeriod.EqualTo("LHC15o") ) OpenInfoCalbration(run);
      if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r") ) OpenInfoCalbration18(run);
      fRunNum = run;
      fRunNumBin = GetRunNumBin(fRunNum);
      if (fRunNumBin<0) return;
  } 
  hRunNumBin->Fill(fRunNumBin);
  hEvtCount->Fill(4);
  if (fDebug) Printf("run nummbr done!");

  // vertex
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  double vx    = fVtx->GetX();
  double vy    = fVtx->GetY();
  double vz    = fVtx->GetZ();
  hVz[0]->Fill(vz);
  double dz    = vz - fAOD->GetPrimaryVertexSPD()->GetZ();
  if (fabs(vz)>(double)fZvtxCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors()<1) return;
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
  else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r") ){
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    centV0M    =  fMultSel->GetMultiplicityPercentile("V0M");
    centTRK    =  fMultSel->GetMultiplicityPercentile("TRK");
    centSPD0 =  fMultSel->GetMultiplicityPercentile("CL0");
    centSPD1 =  fMultSel->GetMultiplicityPercentile("CL1");
  }
  fCent  = centV0M; // default

  // Trigger check 
  // if (mask&AliVEvent::kCentral && fCent>10.0) return;
  // if (mask&AliVEvent::kSemiCentral && (fCent < 30.0 || fCent>50.0) ) return;

  if (fCentEst.EqualTo("V0M")) fCent  = centV0M;
  else if (fCentEst.EqualTo("SPD0")) fCent  = centSPD0;
  else if (fCentEst.EqualTo("SPD1")) fCent  = centSPD1;
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
  if (fPUSyst==0 && fPeriod.EqualTo("LHC10h")) {if(!RemovalForRun1(fAOD, fUtils) ) return;}
  if (fPeriod.EqualTo("LHC11h")) {if(!RemovalForRun1(fAOD, fUtils) ) return;}
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    // hMultCentQA[0]->Fill(fCent, fAOD->GetNumberOfTracks()); // raw Trk Multi Vs Cent(V0M)
    // if (PileUpMultiVertex(fAOD)) return;
    // if (!RejectEvtMultComp(fAOD)) return;
    // hMultCentQA[1]->Fill(fCent, fAOD->GetNumberOfTracks()); // Mult_Cent QA
    // if (!AODPileupCheck (fAOD)) return;
    if (!RejectEvtTFFit (fAOD)) return; // 15o_pass2
  }
  hCent[1]->Fill(fCent);
  if (fPeriod.EqualTo("LHC15o") && fPUSyst==3) hCentRBR[fRunNumBin]->Fill(fCent);
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
  hReQ_thisEvt->Reset();
  hImQ_thisEvt->Reset();
  hReQ2_thisEvt->Reset();
  hImQ2_thisEvt->Reset();
  hMQ_thisEvt->Reset();  
  hMQ_weight_thisEvt->Reset(); 
  hReQPos_thisEvt->Reset();
  hImQPos_thisEvt->Reset();
  hMQPos_thisEvt ->Reset();
  hMQPos_weight_thisEvt ->Reset();
  hReQNeg_thisEvt->Reset();
  hImQNeg_thisEvt->Reset();
  hMQNeg_thisEvt ->Reset();  
  hMQNeg_weight_thisEvt ->Reset();  
  pRefFlow_thisEvt ->Reset();  
  pIntd2_thisEvt->Reset();
  pbufferPt_thisEvt->Reset();
  pbufferEta_thisEvt->Reset();
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

    if (fPUSyst==2 && fPeriod.EqualTo("LHC15o")){
        if(TMath::Abs(track->GetTOFsignalDz())>10.)
            continue;
        if(track->GetTOFsignal() < 12000.)
            continue;
        if(track->GetTOFsignal() > 25000.)
            continue;
    }

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

    hReQ_thisEvt->Fill(pt, eta, cosnphi);
    hImQ_thisEvt->Fill(pt, eta, sinnphi);
    hReQ2_thisEvt->Fill(pt,eta,cos2nphi); 
    hImQ2_thisEvt->Fill(pt, eta, sin2nphi); 
    hMQ_thisEvt ->Fill(pt, eta);
    hMQ_weight_thisEvt ->Fill(pt, eta, weight);
    if (charge>0) {
      hReQPos_thisEvt->Fill(pt, eta, cosnphi);
      hImQPos_thisEvt->Fill(pt, eta, sinnphi);
      hMQPos_thisEvt ->Fill(pt, eta);
      hMQPos_weight_thisEvt->Fill(pt, eta, weight);
    } else if (charge<0){
      hReQNeg_thisEvt->Fill(pt, eta, cosnphi);
      hImQNeg_thisEvt->Fill(pt, eta, sinnphi);
      hMQNeg_thisEvt ->Fill(pt, eta);
      hMQNeg_weight_thisEvt ->Fill(pt, eta, weight);
    }
    //------------------
    // for pT,eta vs ach
    //------------------
    if (charge>0) pbufferPt_thisEvt->Fill(0.5, pt, weight);
    if (charge<0) pbufferPt_thisEvt->Fill(1.5, pt, weight);

    if (charge>0) pbufferEta_thisEvt->Fill(0.5, eta, weight);
    if (charge<0) pbufferEta_thisEvt->Fill(1.5, eta, weight);

    Mult++;
    isSltTrk[iTrk] = 1;
  }  

  hEvtCount->Fill(9);
  if (!DoCumulants()) return false;  // REF V2
  hEvtCount->Fill(10);
  //------------------
  // Ach
  //------------------
  if (nPosTrk<1 || nNegTrk<1) return false;
  double mAch = (double)(nPosTrk-nNegTrk)/(nPosTrk+nNegTrk);
  if (fabs(mAch)>1) return false;
  hEvtCount->Fill(11); 
  //------------------
  // Qn vector
  //------------------
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
  
  // if (fPeriod.EqualTo("LHC15o")&&fPUSyst==3) {
  //   double wCent = 1.0; 
  //   int centBin = hwCent->FindBin(fCent);
  //   wCent = hwCent->GetBinContent(centBin);
  //   pRefFlow[fCentBin]->Add(pRefFlow_thisEvt, wCent);
  //   pIntd2[fCentBin]->Add(pIntd2_thisEvt, wCent);
  //   pIntd2_thisEvt->Scale(mAch); 
  //   pIntd2Ach[fCentBin]->Add(pIntd2_thisEvt, wCent);
  //   hEvtCount->Fill(14);
  // }

  pRefFlow[fCentBin]->Add(pRefFlow_thisEvt);
  pIntd2[fCentBin]->Add(pIntd2_thisEvt);
  pIntd2_thisEvt->Scale(mAch); 
  pIntd2Ach[fCentBin]->Add(pIntd2_thisEvt);
  hEvtCount->Fill(14);

  pPosPtAch[fCentBin]->Fill(mAch, pbufferPt_thisEvt->GetBinContent(1) );    
  pNegPtAch[fCentBin]->Fill(mAch, pbufferPt_thisEvt->GetBinContent(2) );    
  pPosEtaAch[fCentBin]->Fill(mAch, pbufferEta_thisEvt->GetBinContent(1) );  
  pNegEtaAch[fCentBin]->Fill(mAch, pbufferEta_thisEvt->GetBinContent(2) ); 
  hEvtCount->Fill(15);

  return true;
}

//---------------------------------------------------

double AliAnalysisTaskCMWESE::GetNUECor(int charge, double pt)
{
  double weightNUE = 1;
  if (fPeriod.EqualTo("LHC10h")){
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
  else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){
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
    hNUAweightPlus = (TH2D*)fListNUA->FindObject(Form("weightdPhidEta_run%i_cent0_vz%i_plus",fRunNum, fVzBin));
    hNUAweightMinus = (TH2D*)fListNUA->FindObject(Form("weightdPhidEta_run%i_cent0_vz%i_minus",fRunNum, fVzBin));
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
  } else if (fPeriod.EqualTo("LHC15o")){ // Rihan and Protty 's NUA Results
    if (charge>0){ 
      hCorrectNUAPos = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Pos_Cent0_Run%i",fRunNum));
      if (!hCorrectNUAPos) return -1;
      int iBinNUA = hCorrectNUAPos->FindBin(vz,phi,eta); 
      if (hCorrectNUAPos->GetBinContent(iBinNUA)>0) weightNUA = (double)hCorrectNUAPos->GetBinContent(iBinNUA);
      return  weightNUA;
    } else if (charge<0){
      hCorrectNUANeg = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Neg_Cent0_Run%i",fRunNum));
      if (!hCorrectNUANeg) return -1;
      int iBinNUA = hCorrectNUANeg->FindBin(vz,phi,eta); 
      if (hCorrectNUANeg->GetBinContent(iBinNUA)>0) weightNUA = (double)hCorrectNUANeg->GetBinContent(iBinNUA);
      return weightNUA;
    } 
    // In Rihan and Protty 's NUA results, the phi distribution is independent on centrality and particle charge
  } else if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){ // Rihan and Protty 's NUA Results
    if (charge>0){ 
      hCorrectNUAPos = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d",0,fRunNum));
      if (!hCorrectNUAPos) return -1;
      int iBinNUA = hCorrectNUAPos->FindBin(vz,phi,eta); 
      if (hCorrectNUAPos->GetBinContent(iBinNUA)>0) weightNUA = (double)hCorrectNUAPos->GetBinContent(iBinNUA);
      return  weightNUA;
    } else if (charge<0){
      hCorrectNUANeg = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d",0,fRunNum));
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
      int runNumList[90]={
        139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310,
        139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870,
        138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582,
        138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201,
        138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693,
        137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541,
        137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243,
        137236, 137235, 137232, 137231, 137162, 137161};
      for (int i = 0; i < 90; ++i) {
        if (runNum==runNumList[i]) {runNBin=i; break;}
        else continue;
      }
  } else if (fPeriod.EqualTo("LHC11h") ){
      int runNumList[39]={
        170387,170040,170268,170228,170207,169838,170159,170204,170311,170084,
        169835,170088,170593,170203,170270,169846,170163,170388,170155,170083,
        170572,169837,169855,170306,170269,170089,170309,170091,170081,170230,
        170085,170315,170027,170193,170312,170313,170308,169858,169859};
        for (int i = 0; i < 39; ++i) {
          if (runNum==runNumList[i]) {runNBin=i; break;}
          else continue;
        }
  } else if (fPeriod.EqualTo("LHC15o") ){
      int runNumList[138]={
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
        if (runNum==runNumList[i]) {runNBin=i; break;}
        else continue;
      }
  } else if (fPeriod.EqualTo("LHC15o") ){
      int runNumList[138]={
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
        if (runNum==runNumList[i]) {runNBin=i; break;}
        else continue;
      }
  } else if (fPeriod.EqualTo("LHC18q") ){
      int runNumList[125]={
        296623, 296622, 296621, 296619, 296618, 296616, 296615, 296594, 296553, 296552, 296551, 296550,
        296548, 296547, 296516, 296512, 296511, 296510, 296509, 296472, 296433, 296424, 296423, 296420,
        296419, 296415, 296414, 296383, 296381, 296380, 296379, 296378, 296377, 296376, 296375, 296312,
        296309, 296304, 296303, 296280, 296279, 296273, 296270, 296269, 296247, 296246, 296244, 296243,
        296242, 296241, 296240, 296198, 296197, 296196, 296195, 296194, 296192, 296191, 296143, 296142,
        296135, 296134, 296133, 296132, 296123, 296074, 296066, 296065, 296063, 296062, 296060, 296016,
        295942, 295941, 295937, 295936, 295913, 295910, 295909, 295861, 295860, 295859, 295856, 295855,
        295854, 295853, 295831, 295829, 295826, 295825, 295822, 295819, 295818, 295816, 295791, 295788,
        295786, 295763, 295762, 295759, 295758, 295755, 295754, 295725, 295723, 295721, 295719, 295718,
        295717, 295714, 295712, 295676, 295675, 295673, 295668, 295667, 295666, 295615, 295612, 295611,
        295610, 295589, 295588, 295586, 295585
      };
      for (int i = 0; i < 125; ++i) {
        if (runNum==runNumList[i]) {runNBin=i; break;}
        else continue;
      }
  } else if (fPeriod.EqualTo("LHC18r") ){
      int runNumList[89]={
        297595, 297590, 297588, 297558, 297544, 297542, 297541, 297540, 297537, 297512, 297483, 297479,
        297452, 297451, 297450, 297446, 297442, 297441, 297415, 297414, 297413, 297406, 297405, 297380,
        297379, 297372, 297367, 297366, 297363, 297336, 297335, 297333, 297332, 297317, 297311, 297310,
        297278, 297222, 297221, 297218, 297196, 297195, 297193, 297133, 297132, 297129, 297128, 297124,
        297123, 297119, 297118, 297117, 297085, 297035, 297031, 296966, 296941, 296938, 296935, 296934,
        296932, 296931, 296930, 296903, 296900, 296899, 296894, 296852, 296851, 296850, 296848, 296839,
        296838, 296836, 296835, 296799, 296794, 296793, 296790, 296787, 296786, 296785, 296784, 296781,
        296752, 296694, 296693, 296691, 296690
      };
      for (int i = 0; i < 89; ++i) {
        if (runNum==runNumList[i]) {runNBin=i; break;}
        else continue;
      }
  }
  return runNBin;
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
  bool isPileup = fUtils->IsPileUpEvent(fAOD);
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
        // if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
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

  Int_t tpcClsTot = fAOD->GetNumberOfTPCClusters();
  // Float_t nclsDif = Float_t(tpcClsTot) - (53182.6 + 113.326*multV0Tot - 0.000831275*multV0Tot*multV0Tot); // 15o_p2
  Float_t nclsDif = Float_t(tpcClsTot) - (50580.2 + 67.6*multV0Tot - 0.00021*multV0Tot*multV0Tot); //18

  // pile-up cuts
  if (centCL0 < fCenCutLowPU->Eval(centV0M)) return false;        
  if (centCL0 > fCenCutHighPU->Eval(centV0M)) return false;
           
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) return false;
           
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;

  if (Float_t(multTrk) < fMultCutPU->Eval(centV0M)) return false;

  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;

  if (fAOD->IsIncompleteDAQ()) return false;

  if (fPUSyst==1){
     if (nclsDif > 200000)//can be varied to 150000, 200000
      return false;
  }

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
    if ((fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) && fDcaCutz==2.0){
      fDcaCutxy = 7*(0.0026 + 0.005/pow(pt, 1.01));
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

bool AliAnalysisTaskCMWESE::DoCumulants()
{
  if (hMQ_thisEvt->Integral(1, 100, 1, 32) < 2) return false;
  double allM  =  hMQ_weight_thisEvt->Integral(1, 100, 1, 32);

  //------------------
  // Q
  //------------------
  double allReQ  = hReQ_thisEvt ->Integral(1, 100, 1, 32);
  double allImQ  = hImQ_thisEvt ->Integral(1, 100, 1, 32);
  TComplex Q(allReQ, allImQ);
  TComplex QStar = TComplex::Conjugate(Q);

  double posEtaRe = hReQ_thisEvt->Integral(1, 100, hReQ_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  double posEtaIm = hImQ_thisEvt->Integral(1, 100, hImQ_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  if (hMQ_thisEvt->Integral(1, 100,  hMQ_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32) < 2) return false;
  fPosEtaMQ =  hMQ_weight_thisEvt->Integral(1, 100,  hMQ_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  fPosEtaQ = TComplex(posEtaRe, posEtaIm);
  fPosEtaQStar = TComplex::Conjugate(fPosEtaQ);

  double negEtaRe = hReQ_thisEvt->Integral(1, 100, 1, hReQ_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  double negEtaIm = hImQ_thisEvt->Integral(1, 100, 1, hImQ_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  if (hMQ_thisEvt->Integral(1, 100, 1,  hMQ_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6)) < 2) return false;
  fNegEtaMQ =  hMQ_weight_thisEvt->Integral(1, 100, 1,  hMQ_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  fNegEtaQ = TComplex(negEtaRe, negEtaIm);
  fNegEtaQStar = TComplex::Conjugate(fNegEtaQ);

  // //------------------
  // // Q+
  // //------------------
  // double negEtaRePos = hReQPos_thisEvt->Integral(1, 100, 1, hReQPos_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  // double negEtaImPos = hImQPos_thisEvt->Integral(1, 100, 1, hImQPos_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  // if (hMQPos_thisEvt->Integral(1, 100, 1,  hMQPos_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6)) < 2) return false;
  // double negEtaMQPos =  hMQPos_weight_thisEvt->Integral(1, 100, 1,  hMQPos_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  // TComplex negEtaQPos(negEtaRePos, negEtaImPos);
  // TComplex negEtaQStarPos = TComplex::Conjugate(negEtaQPos);

  // double posEtaRePos = hReQPos_thisEvt->Integral(1, 100, hReQPos_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  // double posEtaImPos = hImQPos_thisEvt->Integral(1, 100, hImQPos_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  // if (hMQPos_thisEvt->Integral(1, 100,  hMQPos_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32) < 2) return false;
  // double posEtaMQPos =  hMQPos_weight_thisEvt->Integral(1, 100,  hMQPos_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  // TComplex posEtaQPos(posEtaRePos, posEtaImPos);
  // TComplex posEtaQStarPos = TComplex::Conjugate(posEtaQPos);

  // //------------------
  // // Q-
  // //------------------
  // double negEtaReNeg = hReQNeg_thisEvt->Integral(1, 100, 1, hReQNeg_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  // double negEtaImNeg = hImQNeg_thisEvt->Integral(1, 100, 1, hImQNeg_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  // if (hMQNeg_thisEvt->Integral(1, 100, 1,  hMQNeg_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6)) < 2) return false;
  // double negEtaMQNeg =  hMQNeg_weight_thisEvt->Integral(1, 100, 1,  hMQNeg_thisEvt->GetYaxis()->FindBin(-fEtaGap-1e-6));
  // TComplex negEtaQNeg(negEtaReNeg, negEtaImNeg);
  // TComplex negEtaQStarNeg = TComplex::Conjugate(negEtaQNeg);

  // double posEtaReNeg = hReQNeg_thisEvt->Integral(1, 100, hReQNeg_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  // double posEtaImNeg = hImQNeg_thisEvt->Integral(1, 100, hImQNeg_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  // if (hMQNeg_thisEvt->Integral(1, 100,  hMQNeg_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32) < 2) return false;
  // double posEtaMQNeg =  hMQNeg_weight_thisEvt->Integral(1, 100,  hMQNeg_thisEvt->GetYaxis()->FindBin(fEtaGap+1e-6), 32);
  // TComplex posEtaQNeg(posEtaReNeg, posEtaImNeg);
  // TComplex posEtaQStarNeg = TComplex::Conjugate(posEtaQNeg);

  //------------------
  // c{2}
  //------------------
  pRefFlow_thisEvt->Fill(0.,  (fNegEtaQ*fPosEtaQStar).Re()       / (fNegEtaMQ*fPosEtaMQ),       (fNegEtaMQ*fPosEtaMQ));
  // pRefFlow_thisEvt->Fill(1.,  (fNegEtaQ*posEtaQStarPos).Re()    / (fNegEtaMQ*posEtaMQPos),    (fNegEtaMQ*posEtaMQPos));
  // pRefFlow_thisEvt->Fill(2.,  (negEtaQPos*fPosEtaQStar).Re()    / (negEtaMQPos*fPosEtaMQ),    (negEtaMQPos*fPosEtaMQ));
  // pRefFlow_thisEvt->Fill(3.,  (fNegEtaQ*posEtaQStarNeg).Re()    / (fNegEtaMQ*posEtaMQNeg),    (fNegEtaMQ*posEtaMQNeg));
  // pRefFlow_thisEvt->Fill(4.,  (negEtaQNeg*fPosEtaQStar).Re()    / (negEtaMQNeg*fPosEtaMQ),    (negEtaMQNeg*fPosEtaMQ));
  // pRefFlow_thisEvt->Fill(5.,  (negEtaQPos*posEtaQStarPos).Re() / (negEtaMQPos*posEtaMQPos), (negEtaMQPos*posEtaMQPos));
  // pRefFlow_thisEvt->Fill(6.,  (negEtaQNeg*posEtaQStarNeg).Re() / (negEtaMQNeg*posEtaMQNeg), (negEtaMQNeg*posEtaMQNeg));
  // pRefFlow_thisEvt->Fill(7.,  (negEtaQPos*posEtaQStarNeg).Re() / (negEtaMQPos*posEtaMQNeg), (negEtaMQPos*posEtaMQNeg));
  // pRefFlow_thisEvt->Fill(8.,  (negEtaQNeg*posEtaQStarPos).Re() / (negEtaMQNeg*posEtaMQPos), (negEtaMQNeg*posEtaMQPos));
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
  } else if (fPeriod.EqualTo("LHC15o") ){ // A.Dobrin's V0Calib; "calibV0HIR.root"
      for(int iCh = 0; iCh < 64; ++iCh) {
        fMultV0Ch[iCh] = hMultV0->GetBinContent(iCh+1);
      } 
      // AliCentrality* centrality = ((AliAODHeader*)fAOD->GetHeader())->GetCentralityP();
      // Double_t centSPD = centrality->GetCentralityPercentile("CL1");
      AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
      if(!fMultSelection) {
        printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
        exit(1);
      }
      Double_t centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
      Int_t iCentSPD = (Int_t)centrCL1;
      if (iCentSPD >= 90) return false;
      fV0XMean[0] = -999.; fV0YMean[0] = -999.; 
      for(int i = 0; i < 2; ++i) {   // [1]: C; [2]: A;
        fV0XMean[i+1] = hQxnmV0[i]->GetBinContent(iCentSPD+1);
        fV0YMean[i+1] = hQynmV0[i]->GetBinContent(iCentSPD+1);
      }   
  } else if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r") ){
      AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
      if(!fMultSelection) {
        printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
        exit(1);
      }
      Double_t centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
      Int_t iCentSPD = (Int_t)centrCL1;
      if (iCentSPD >= 90) return false;
      fV0XMean[0] = -999.; fV0YMean[0] = -999.; 
      for(int i = 0; i < 2; ++i) {   // [1]: C; [2]: A;
        fV0XMean[i+1] = hQxnmV0[i]->GetBinContent(iCentSPD+1);
        fV0YMean[i+1] = hQynmV0[i]->GetBinContent(iCentSPD+1);
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
    // if (TMath::IsNaN(fMultV0Ch[iCh]) || fMultV0Ch[iCh]<=0) continue;

    double phi = TMath::Pi()/8. + TMath::Pi()/4.*(iCh%8);
    double multCh = 0.;
    // double multCh = fAOD->GetVZEROEqMultiplicity(iCh);
    if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC11h")) multCh= fAOD->GetVZEROEqMultiplicity(iCh);
    else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      AliAODVZERO* aodV0 = fAOD->GetVZEROData();
      multCh = aodV0->GetMultiplicity(iCh);
    }
    if (iCh<32) { // C
      double multChGEC = -1;
      if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC15o")){
        if (iCh < 8)
            multChGEC = multCh/fMultV0Ch[iCh] * fMultV0Ch[0];
        else if (iCh >= 8 && iCh < 16)
            multChGEC = multCh/fMultV0Ch[iCh] * fMultV0Ch[8];
        else if (iCh >= 16 && iCh < 24)
            multChGEC = multCh/fMultV0Ch[iCh] * fMultV0Ch[16];
        else if (iCh >= 24 && iCh < 32)
            multChGEC = multCh/fMultV0Ch[iCh] * fMultV0Ch[24];
        if (multChGEC<0 || TMath::IsNaN(multChGEC)) continue;
      }
      if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){
        double vz = fVtx->GetZ();
        int ibinV0 = fHCorrectV0ChWeghts->FindBin(vz,iCh);
        double V0chGE = (double)fHCorrectV0ChWeghts->GetBinContent(ibinV0);
        multChGEC = multCh*V0chGE;
        if (multChGEC<0) continue;
      }
      qxGE[1] += multChGEC*TMath::Cos(fHarmonic*phi);
      qyGE[1] += multChGEC*TMath::Sin(fHarmonic*phi);
      multRingGE[1] += multChGEC;
      qxGE[0] += multChGEC*TMath::Cos(fHarmonic*phi);
      qyGE[0] += multChGEC*TMath::Sin(fHarmonic*phi); 
      multRingGE[0] += multChGEC; 

    } else if (iCh>=32 && iCh<64) { // A
      double multChGEA = -1;
      if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC15o")){
        if (iCh >= 32 && iCh < 40)
            multChGEA = multCh/fMultV0Ch[iCh] * fMultV0Ch[32];
        else if (iCh >= 40 && iCh < 48)
            multChGEA = multCh/fMultV0Ch[iCh] * fMultV0Ch[40];
        else if (iCh >= 48 && iCh < 56)
            multChGEA = multCh/fMultV0Ch[iCh] * fMultV0Ch[48];
        else if (iCh >= 56 && iCh < 64)
            multChGEA = multCh/fMultV0Ch[iCh] * fMultV0Ch[56];
        if (multChGEA<0 || TMath::IsNaN(multChGEA)) continue;
      }
      if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){
        double vz = fVtx->GetZ();  
        int ibinV0 = fHCorrectV0ChWeghts->FindBin(vz,iCh);
        double V0chGE = (double)fHCorrectV0ChWeghts->GetBinContent(ibinV0);
        multChGEA = multCh*V0chGE;
        if (multChGEA<0) continue;
      }
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
      if (TMath::IsNaN(qxMean) || TMath::IsNaN(qyMean)) continue;
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
        if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){
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
  if (fPeriod.EqualTo("LHC10h")  || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r") ){
    hQnPercentile->GetXaxis()->SetRange((int)fCent+1,(int)fCent+1);  
    hQnPercentile_centThisEvt = (TH1D*)hQnPercentile->ProjectionY("qnPercentileCentiEvt", 0, -1, "e");
    sp = new TSpline3(hQnPercentile_centThisEvt);
    double pec = 100* sp->Eval(qn_thisEvt);
    int QnBin = GetPercCode(pec);
    hAchMb[QnBin]->Fill(mAch);
    hMultMb[QnBin]->Fill(Mult); //mb

    hQnPercentile->GetXaxis()->SetRange(1,100); 
    return QnBin;   
  } 

  else if (fPeriod.EqualTo("LHC15o")){
    double pec = 100.* splQ2c[(int)fCent]->Eval(qn_thisEvt);
    int QnBin = GetPercCode(pec);
    hAchMb[QnBin]->Fill(mAch);
    hMultMb[QnBin]->Fill(Mult); //mb
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

//---------------------------------------------------

void AliAnalysisTaskCMWESE::OpenInfoCalbration(Int_t run)
{
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
      AliOADBContainer* contMult = (AliOADBContainer*) fileV0Calib->FindObjectAny("hMultV0BefCorPfpx");
      if(!contMult){
          printf("OADB object hMultV0BefCorr is not available in the file\n");
          return;
      }
  
      // V0C Qx Mean 
      AliOADBContainer* contQxncm = (AliOADBContainer*) fileV0Calib->FindObjectAny(Form("fqxc%im",(int)fHarmonic));
      if(!contQxncm){
          printf("OADB object fqxcnm is not available in the file\n");
          cout<<fHarmonic<<endl;
          return;
      }
  
      // V0C Qy Mean 
      AliOADBContainer* contQyncm = (AliOADBContainer*) fileV0Calib->Get(Form("fqyc%im",(int)fHarmonic));
      if(!contQyncm){
          printf("OADB object fqycnm is not available in the file\n");
          return;
      }
  
      // V0A Qx Mean 
      AliOADBContainer* contQxnam = (AliOADBContainer*) fileV0Calib->Get(Form("fqxa%im",(int)fHarmonic));
      if(!contQxnam){
          printf("OADB object fqxanm is not available in the file\n");
          return;
      }
  
      // V0A Qy Mean 
      AliOADBContainer* contQynam = (AliOADBContainer*) fileV0Calib->Get(Form("fqya%im",(int)fHarmonic));
      if(!contQynam){
          printf("OADB object fqyanm is not available in the file\n");
          return;
      }

        // Mult   
        if(!(contMult->GetObject(run) ) ){
            printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", run);
            return;
        }
        hMultV0 = ((TH1D*) contMult->GetObject(run ) );
  
        // V0C Qx Mean 
        if(!(contQxncm->GetObject(run) ) ){
            printf("OADB object fqxcnm is not available for run %i\n", run);
            return;
        }
        hQxnmV0[0] = ((TH1D*) contQxncm->GetObject(run) );
  
        // V0C Qy Mean 
        if(!(contQyncm->GetObject(run) ) ){
            printf("OADB object fqycnm is not available for run %i\n", run);
            return;
        }
        hQynmV0[0] = ((TH1D*) contQyncm->GetObject(run) );
  
        // V0A Qx Mean
        if(!(contQxnam->GetObject(run) ) ){
            printf("OADB object fqxanm is not available for run %i\n", run);
            return;
        }
        hQxnmV0[1] = ((TH1D*) contQxnam->GetObject(run) );
  
        // V0A Qy Mean
        if(!(contQynam->GetObject(run) ) ){
            printf("OADB object fqyanm is not available for run %i\n", run);
            return;
        }
        hQynmV0[1] = ((TH1D*) contQynam->GetObject(run) );    
}

//---------------------------------------------------

void AliAnalysisTaskCMWESE::OpenInfoCalbration18(Int_t run){ 
    
    fHCorrectV0ChWeghts = (TH2F *) fListVZEROCALIB->FindObject(Form("hWgtV0ChannelsvsVzRun%d",run));

    // V0C Qx Mean 
    AliOADBContainer* contQxncm = (AliOADBContainer*) fListVZEROCALIB->FindObject(Form("fqxc%im",(int)fHarmonic));
    if(!contQxncm){
        printf("OADB object fqxc2m is not available in the file\n");
        return;
    }
       
    // V0C Qy Mean 
    AliOADBContainer* contQyncm = (AliOADBContainer*) fListVZEROCALIB->FindObject(Form("fqyc%im",(int)fHarmonic));
    if(!contQyncm){
        printf("OADB object fqyc2m is not available in the file\n");
        return;
    }
       
    // V0A Qx Mean 
    AliOADBContainer* contQxnam = (AliOADBContainer*) fListVZEROCALIB->FindObject(Form("fqxa%im",(int)fHarmonic));
    if(!contQxnam){
        printf("OADB object fqxa2m is not available in the file\n");
        return;
    }
    
    // V0A Qy Mean 
    AliOADBContainer* contQynam = (AliOADBContainer*) fListVZEROCALIB->FindObject(Form("fqya%im",(int)fHarmonic));
    if(!contQynam){
        printf("OADB object fqya2m is not available in the file\n");
        return;
    }
  
    // V0C Qx Mean 
    if(!(contQxncm->GetObject(run) ) ){
        printf("OADB objectForm fqxcnm ,fHarmonic)is not available for run %i\n", run);
        return;
    }
    hQxnmV0[0] = ((TH1D*) contQxncm->GetObject(run) );
    
    // V0C Qy Mean 
    if(!(contQyncm->GetObject(run) ) ){
        printf("OADB object fqycnm is not available for run %i\n", run);
        return;
    }
    hQynmV0[0] = ((TH1D*) contQyncm->GetObject(run) );
    
    // V0A Qx Mean
    if(!(contQxnam->GetObject(run) ) ){
        printf("OADB object fqxanm is not available for run %i\n", run);
        return;
    }
    hQxnmV0[1] = ((TH1D*) contQxnam->GetObject(run) );
         
    // V0A Qy Mean
    if(!(contQynam->GetObject(run) ) ){
        printf("OADB object fqyanm is not available for run %i\n", run);
        return;
    }
    hQynmV0[1] = ((TH1D*) contQynam->GetObject(run) );

}
