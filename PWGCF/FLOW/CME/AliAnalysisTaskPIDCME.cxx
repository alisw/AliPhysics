/**************************************************************************
 * Copyright(c) 1998-2023, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------------------
// Task for PID CME Analysis
// Contributor: Chunzheng Wang, <chunzheng.wang@cern.ch> & Zhengqing Wang, <zhengqing.wang@cern.ch>, Fudan University, Shanghai
//-----------------------------------------------------------------------------------

#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
// ROOT classes
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TSpline.h"
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
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliOADBContainer.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskPIDCME.h"

ClassImp(AliAnalysisTaskPIDCME);

//---------------------------------------------------
AliAnalysisTaskPIDCME::AliAnalysisTaskPIDCME() :
  AliAnalysisTaskSE(),
  isUseTPCPlane(true),
  isUseVZEROPlane(true),
  isUseZDCPlane(false),
  isDoNUE(true),
  isDoNUA(false),
  isQATPC(true),
  isQAVZERO(true),
  isQAZDC(true),
  isNarrowDcaCuts768(false),
  isProtonCustomizedDCACut(false),
  isUsePionRejection(false),
  isCalculatePIDFlow(true),
  isCalculateDiffResult(true),
  isCalculateDeltaPhiSumPhi(false),
  isCalculatePionKaon(true),
  isCalculatePionProton(true),
  isCalculateKaonProton(true),
  isCalculatePionPion(true),
  isCalculateKaonKaon(false),
  isCalculateProtonProton(false),
  isCalculateHadronHadron(true),
  isUseOneSideTPCPlane(false),
  isTightPileUp(false),
  isOpenPIDSingletrk(true),
  isOpenSsandOsSelfCheck(false),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fPlanePtMin(0.2),
  fPlanePtMax(2.0),
  fEtaGapPos(0.1),
  fEtaGapNeg(-0.1),
  fFilterBit(768),
  fNclsCut(50),
  fChi2Max(4.0),
  fChi2Min(0.1),
  fDcaCutXY(2.4),
  fDcaCutZ(3.2),
  fPtMin(0.2),
  fPtMax(10.0),
  fEtaCut(0.8),
  fDedxCut(10.),
  fPionPtMin(0.2),  
  fPionPtMax(10.0),
  fKaonPtMin(0.2),
  fKaonPtMax(10.0),
  fProtonPtMin(0.2),
  fProtonPtMax(10.0),
  fAntiPionPtMin(0.2),  
  fAntiPionPtMax(10.0),
  fAntiKaonPtMin(0.2),
  fAntiKaonPtMax(10.0),
  fAntiProtonPtMin(0.2),
  fAntiProtonPtMax(10.0),
  fPtMinforRMSPion(0.5),
  fPtMinforRMSKaon(0.45),
  fPtMinforRMSProton(0.2),
  fNSigmaTPCCutPion(3),
  fNSigmaRMSCutPion(3),
  fNSigmaTPCCutKaon(3),
  fNSigmaRMSCutKaon(3),
  fNSigmaTPCCutProton(3),
  fNSigmaRMSCutProton(3),
  fNQ2Bins(6),
  fAOD(nullptr),
  fPIDResponse(nullptr),
  fUtils(nullptr),
  fMultSel(nullptr),
  fRunNum(-999),
  fOldRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCent(-999),
  fCentBin(-999),
  fCentV0M(-999),
  fCentTRK(-999),
  fCentSPD0(-999),
  fCentSPD1(-999),
  fSumQ2xTPCPos(0.),
  fSumQ2yTPCPos(0.),
  fWgtMultTPCPos(0.),
  fSumQ2xTPCNeg(0.),
  fSumQ2yTPCNeg(0.),
  fWgtMultTPCNeg(0.),
  fPsi2TPCPos(-999),
  fPsi2TPCNeg(-999),
  fPsi2V0C(-999),
  fPsi2V0A(-999),
  fPsi1ZNC(-999),
  fPsi1ZNA(-999),
  fPsi2TPC(-999),
  fQ2Bin(-999),
  isRightTPCPlane(false),
  isRightVZEROPlane(false),
  isRightZDCPlane(false),
  mapTPCPosTrksIDPhiWgt(0),
  mapTPCNegTrksIDPhiWgt(0),
  mapTPCTrksIDPhiWgt(0),
  vecParticle(0),
  fSPDCutPU(nullptr),
  fV0CutPU(nullptr),
  fCenCutLowPU(nullptr),
  fCenCutHighPU(nullptr),
  fMultCutPU(nullptr),
  fListNUE(nullptr),
  hNUEweightPlus(nullptr),
  hNUEweightMinus(nullptr),
  fListNUA(nullptr),
  hNUAweightPlus(nullptr),
  hNUAweightMinus(nullptr),
  hCorrectNUAPos(nullptr),
  hCorrectNUANeg(nullptr),
  fListVZEROCalib(nullptr),
  hMultV0Read(nullptr),
  hMultV0(nullptr),
  contMult(nullptr),
  contQxncm(nullptr),
  contQyncm(nullptr),
  contQxnam(nullptr),
  contQynam(nullptr),
  fHCorrectV0ChWeghts(nullptr),
  hQnPercentile(nullptr),
  fListZDCCalib(nullptr),
  tree(nullptr),
  fProfileForZNCGE(nullptr),
  fProfileForZNAGE(nullptr),
  fHn4DForZNCQxRC(nullptr),
  fHn4DForZNCQyRC(nullptr),
  fHn4DForZNCMtRC(nullptr),
  fHn4DForZNAQxRC(nullptr),
  fHn4DForZNAQyRC(nullptr),
  fHn4DForZNAMtRC(nullptr),
  fHn4DForZNCCountsRC(nullptr),
  fHn4DForZNACountsRC(nullptr),
  fProfile2DForCosC(nullptr),
  fProfile2DForSinC(nullptr),
  fProfile2DForCosA(nullptr),
  fProfile2DForSinA(nullptr),
  fHZDCCparameters(nullptr),
  fHZDCAparameters(nullptr),
  fQAList(nullptr),
  fEvtCount(nullptr),
  runNumList(nullptr),
  fHistRunNumBin(nullptr),
  fHistPt(nullptr),
  fHistEta(nullptr),
  fHistNhits(nullptr),
  fHist2PDedx(nullptr),
  fHistDcaXY(nullptr),
  fHistDcaZ(nullptr),
  fHistPionPt(nullptr),
  fHistPionEta(nullptr),
  fHistPionPhi(nullptr),
  fHistPionDcaXY(nullptr),
  fHistPionDcaZ(nullptr),
  fHist2PionSigTPC(nullptr),
  fHist2PionSigTOF(nullptr),
  fHist2PionSigRMS(nullptr),
  fHist3PionSigTPCTOFPt(nullptr),
  fHistAntiPionPt(nullptr),
  fHistAntiPionEta(nullptr),
  fHistAntiPionPhi(nullptr),
  fHistAntiPionDcaXY(nullptr),
  fHistAntiPionDcaZ(nullptr),
  fHist2AntiPionSigTPC(nullptr),
  fHist2AntiPionSigTOF(nullptr),
  fHist2AntiPionSigRMS(nullptr),
  fHist3AntiPionSigTPCTOFPt(nullptr),
  fHistKaonPt(nullptr),
  fHistKaonEta(nullptr),
  fHistKaonPhi(nullptr),
  fHistKaonDcaXY(nullptr),
  fHistKaonDcaZ(nullptr),
  fHist2KaonSigTPC(nullptr),
  fHist2KaonSigTOF(nullptr),
  fHist2KaonSigRMS(nullptr),
  fHist3KaonSigTPCTOFPt(nullptr),
  fHistAntiKaonPt(nullptr),
  fHistAntiKaonEta(nullptr),
  fHistAntiKaonPhi(nullptr),
  fHistAntiKaonDcaXY(nullptr),
  fHistAntiKaonDcaZ(nullptr),
  fHist2AntiKaonSigTPC(nullptr),
  fHist2AntiKaonSigTOF(nullptr),
  fHist2AntiKaonSigRMS(nullptr),
  fHist3AntiKaonSigTPCTOFPt(nullptr),
  fHistProtonPt(nullptr),
  fHistProtonEta(nullptr),
  fHistProtonPhi(nullptr),
  fHistProtonDcaXY(nullptr),
  fHistProtonDcaZ(nullptr),
  fHist2ProtonSigTPC(nullptr),
  fHist2ProtonSigTOF(nullptr),
  fHist2ProtonSigRMS(nullptr),
  fHist3ProtonSigTPCTOFPt(nullptr),
  fHistAntiProtonPt(nullptr),
  fHistAntiProtonEta(nullptr),
  fHistAntiProtonPhi(nullptr),
  fHistAntiProtonDcaXY(nullptr),
  fHistAntiProtonDcaZ(nullptr),
  fHist2AntiProtonSigTPC(nullptr),
  fHist2AntiProtonSigTOF(nullptr),
  fHist2AntiProtonSigRMS(nullptr),
  fHist3AntiProtonSigTPCTOFPt(nullptr),
  fResultsList(nullptr),
  fHist2Psi2TPCPosCent(nullptr),
  fHist2Psi2TPCNegCent(nullptr),
  fHist2Psi2V0CCent(nullptr),
  fHist2Psi2V0ACent(nullptr),
  fHist2Psi1ZNCCent(nullptr),
  fHist2Psi1ZNACent(nullptr),
  fHist2Psi2TPCCent(nullptr),
  fProfileTPCPsi2Correlation(nullptr),
  fProfileV0MPsi2Correlation(nullptr),
  fProfileZDCPsi1Correlation(nullptr),
  fProfileZDCPsi2Correlation(nullptr),
  fProfileV0CTPCPosPsi2Correlation(nullptr),
  fProfileV0ATPCPosPsi2Correlation(nullptr),
  fProfileV0CTPCNegPsi2Correlation(nullptr),
  fProfileV0ATPCNegPsi2Correlation(nullptr),
  fProfileV0CTPCPsi2Correlation(nullptr),
  fProfileV0ATPCPsi2Correlation(nullptr),
  fProfile2RawFlowHadronQ2(nullptr),
  fProfile2RawFlowPionQ2(nullptr),
  fProfile2RawFlowKaonQ2(nullptr),
  fProfile2RawFlowProtonQ2(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;
  for (int i = 0; i < 3; i++) pV0XMeanRead[i] = nullptr;
  for (int i = 0; i < 3; i++) pV0YMeanRead[i] = nullptr;
  for (int i = 0; i < 2; i++) hQx2mV0[i] = nullptr;
  for (int i = 0; i < 2; i++) hQy2mV0[i] = nullptr;
  for (int i = 0; i < 90; i++) splQ2c[i] = nullptr;
  for (int i = 0; i < 3; i++) vtxQuant1[i] = -999;
  for (int i = 0; i < 3; i++) vtxQuant2[i] = -999;
  for (int i = 0; i < 2; i++) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistVz[i] = nullptr;
  for (int i = 0; i < 8; i++) fHist2CentQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2MultCentQA[i] = nullptr;
  for (int i = 0; i < 6; i++) fHist2MultMultQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQxCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQyCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQxVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQyVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fHist2CalibPsi2V0CCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist3CalibQxQyCentV0C[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQxCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQyCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQxVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQyVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fHist2CalibPsi2V0ACent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist3CalibQxQyCentV0A[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileZNCTowerMeanEnegry[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileZNCQxCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileZNCQyCent[i] = nullptr;
  for (int i = 0; i < 3; i++) fHist2CalibPsi1ZNCCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileZNATowerMeanEnegry[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileZNAQxCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileZNAQyCent[i]= nullptr;
  for (int i = 0; i < 3; i++) fHist2CalibPsi1ZNACent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileZDCQxAQxCCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileZDCQxAQyCCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileZDCQyAQxCCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileZDCQyAQyCCent[i]= nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentHadron[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentPion[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentKaon[i][j] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaPionKaon[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaPionProton[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaKaonProton[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaPionPion[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaKaonKaon[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaProtonProton[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaHadronHadron[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaPionKaonSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaPionProtonSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaKaonProtonSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaPionPionSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaKaonKaonSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaProtonProtonSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaHadronHadronSplit[i] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaPionKaonSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaPionProtonSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaKaonProtonSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaPionPionSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaKaonKaonSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaProtonProtonSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaHadronHadronSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaPionKaon[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaPionProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaKaonProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaPionPion[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaKaonKaon[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaProtonProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaHadronHadron[i][j] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionKaonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionProtonDPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonProtonDPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionPionDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonKaonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaProtonProtonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaHadronHadronDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionKaonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionProtonSPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonProtonSPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionPionSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonKaonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaProtonProtonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaHadronHadronSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionKaonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionProtonDEta[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonProtonDEta[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionPionDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonKaonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaProtonProtonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaHadronHadronDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionKaonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionProtonQ2[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonProtonQ2[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionPionQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonKaonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaProtonProtonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaHadronHadronQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionKaonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionProtonDPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonProtonDPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionPionDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonKaonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaProtonProtonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaHadronHadronDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionKaonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionProtonSPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonProtonSPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionPionSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonKaonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaProtonProtonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaHadronHadronSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionKaonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionProtonDEta[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonProtonDEta[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionPionDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonKaonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaProtonProtonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaHadronHadronDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionKaonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionProtonQ2[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonProtonQ2[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionPionQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonKaonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaProtonProtonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaHadronHadronQ2[i] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiPionKaon[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiPionProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiKaonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiPionPion[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiKaonKaon[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiProtonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiHadronHadron[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiPionKaon[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiPionProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiKaonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiPionPion[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiKaonKaon[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiProtonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiHadronHadron[i][j] = nullptr;
}

//---------------------------------------------------
AliAnalysisTaskPIDCME::AliAnalysisTaskPIDCME(const char *name) :
  AliAnalysisTaskSE(name),
  isUseTPCPlane(true),
  isUseVZEROPlane(true),
  isUseZDCPlane(false),
  isDoNUE(true),
  isDoNUA(false),
  isQATPC(true),
  isQAVZERO(true),
  isQAZDC(true),
  isNarrowDcaCuts768(false),
  isProtonCustomizedDCACut(false),
  isUsePionRejection(false),
  isCalculatePIDFlow(true),
  isCalculateDiffResult(true),
  isCalculateDeltaPhiSumPhi(false),
  isCalculatePionKaon(true),
  isCalculatePionProton(true),
  isCalculateKaonProton(true),
  isCalculatePionPion(true),
  isCalculateKaonKaon(false),
  isCalculateProtonProton(false),
  isCalculateHadronHadron(true),
  isUseOneSideTPCPlane(false),
  isTightPileUp(false),
  isOpenPIDSingletrk(true),
  isOpenSsandOsSelfCheck(false),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fPlanePtMin(0.2),
  fPlanePtMax(2.0),
  fEtaGapPos(0.1),
  fEtaGapNeg(-0.1),
  fFilterBit(768),
  fNclsCut(50),
  fChi2Max(4.0),
  fChi2Min(0.1),
  fDcaCutXY(2.4),
  fDcaCutZ(3.2),
  fPtMin(0.2),
  fPtMax(10.0),
  fEtaCut(0.8),
  fDedxCut(10.),
  fPionPtMin(0.2),  
  fPionPtMax(10.0),
  fKaonPtMin(0.2),
  fKaonPtMax(10.0),
  fProtonPtMin(0.2),
  fProtonPtMax(10.0),
  fAntiPionPtMin(0.2),  
  fAntiPionPtMax(10.0),
  fAntiKaonPtMin(0.2),
  fAntiKaonPtMax(10.0),
  fAntiProtonPtMin(0.2),
  fAntiProtonPtMax(10.0),
  fPtMinforRMSPion(0.5),
  fPtMinforRMSKaon(0.45),
  fPtMinforRMSProton(0.2),
  fNSigmaTPCCutPion(3),
  fNSigmaRMSCutPion(3),
  fNSigmaTPCCutKaon(3),
  fNSigmaRMSCutKaon(3),
  fNSigmaTPCCutProton(3),
  fNSigmaRMSCutProton(3),
  fNQ2Bins(6),
  fAOD(nullptr),
  fPIDResponse(nullptr),
  fUtils(nullptr),
  fMultSel(nullptr),
  fRunNum(-999),
  fOldRunNum(-999),
  fRunNumBin(-999),
  fVzBin(-999),
  fCent(-999),
  fCentBin(-999),
  fCentV0M(-999),
  fCentTRK(-999),
  fCentSPD0(-999),
  fCentSPD1(-999),
  fSumQ2xTPCPos(0.),
  fSumQ2yTPCPos(0.),
  fWgtMultTPCPos(0.),
  fSumQ2xTPCNeg(0.),
  fSumQ2yTPCNeg(0.),
  fWgtMultTPCNeg(0.),
  fPsi2TPCPos(-999),
  fPsi2TPCNeg(-999),
  fPsi2V0C(-999),
  fPsi2V0A(-999),
  fPsi1ZNC(-999),
  fPsi1ZNA(-999),
  fPsi2TPC(-999),
  fQ2Bin(-999),
  isRightTPCPlane(false),
  isRightVZEROPlane(false),
  isRightZDCPlane(false),
  mapTPCPosTrksIDPhiWgt(0),
  mapTPCNegTrksIDPhiWgt(0),
  mapTPCTrksIDPhiWgt(0),
  vecParticle(0),
  fSPDCutPU(nullptr),
  fV0CutPU(nullptr),
  fCenCutLowPU(nullptr),
  fCenCutHighPU(nullptr),
  fMultCutPU(nullptr),
  fListNUE(nullptr),
  hNUEweightPlus(nullptr),
  hNUEweightMinus(nullptr),
  fListNUA(nullptr),
  hNUAweightPlus(nullptr),
  hNUAweightMinus(nullptr),
  hCorrectNUAPos(nullptr),
  hCorrectNUANeg(nullptr),
  fListVZEROCalib(nullptr),
  hMultV0Read(nullptr),
  hMultV0(nullptr),
  contMult(nullptr),
  contQxncm(nullptr),
  contQyncm(nullptr),
  contQxnam(nullptr),
  contQynam(nullptr),
  fHCorrectV0ChWeghts(nullptr),
  hQnPercentile(nullptr),
  fListZDCCalib(nullptr),
  tree(nullptr),
  fProfileForZNCGE(nullptr),
  fProfileForZNAGE(nullptr),
  fHn4DForZNCQxRC(nullptr),
  fHn4DForZNCQyRC(nullptr),
  fHn4DForZNCMtRC(nullptr),
  fHn4DForZNAQxRC(nullptr),
  fHn4DForZNAQyRC(nullptr),
  fHn4DForZNAMtRC(nullptr),
  fHn4DForZNCCountsRC(nullptr),
  fHn4DForZNACountsRC(nullptr),
  fProfile2DForCosC(nullptr),
  fProfile2DForSinC(nullptr),
  fProfile2DForCosA(nullptr),
  fProfile2DForSinA(nullptr),
  fHZDCCparameters(nullptr),
  fHZDCAparameters(nullptr),
  fQAList(nullptr),
  fEvtCount(nullptr),
  runNumList(nullptr),
  fHistRunNumBin(nullptr),
  fHistPt(nullptr),
  fHistEta(nullptr),
  fHistNhits(nullptr),
  fHist2PDedx(nullptr),
  fHistDcaXY(nullptr),
  fHistDcaZ(nullptr),
  fHistPionPt(nullptr),
  fHistPionEta(nullptr),
  fHistPionPhi(nullptr),
  fHistPionDcaXY(nullptr),
  fHistPionDcaZ(nullptr),
  fHist2PionSigTPC(nullptr),
  fHist2PionSigTOF(nullptr),
  fHist2PionSigRMS(nullptr),
  fHist3PionSigTPCTOFPt(nullptr),
  fHistAntiPionPt(nullptr),
  fHistAntiPionEta(nullptr),
  fHistAntiPionPhi(nullptr),
  fHistAntiPionDcaXY(nullptr),
  fHistAntiPionDcaZ(nullptr),
  fHist2AntiPionSigTPC(nullptr),
  fHist2AntiPionSigTOF(nullptr),
  fHist2AntiPionSigRMS(nullptr),
  fHist3AntiPionSigTPCTOFPt(nullptr),
  fHistKaonPt(nullptr),
  fHistKaonEta(nullptr),
  fHistKaonPhi(nullptr),
  fHistKaonDcaXY(nullptr),
  fHistKaonDcaZ(nullptr),
  fHist2KaonSigTPC(nullptr),
  fHist2KaonSigTOF(nullptr),
  fHist2KaonSigRMS(nullptr),
  fHist3KaonSigTPCTOFPt(nullptr),
  fHistAntiKaonPt(nullptr),
  fHistAntiKaonEta(nullptr),
  fHistAntiKaonPhi(nullptr),
  fHistAntiKaonDcaXY(nullptr),
  fHistAntiKaonDcaZ(nullptr),
  fHist2AntiKaonSigTPC(nullptr),
  fHist2AntiKaonSigTOF(nullptr),
  fHist2AntiKaonSigRMS(nullptr),
  fHist3AntiKaonSigTPCTOFPt(nullptr),
  fHistProtonPt(nullptr),
  fHistProtonEta(nullptr),
  fHistProtonPhi(nullptr),
  fHistProtonDcaXY(nullptr),
  fHistProtonDcaZ(nullptr),
  fHist2ProtonSigTPC(nullptr),
  fHist2ProtonSigTOF(nullptr),
  fHist2ProtonSigRMS(nullptr),
  fHist3ProtonSigTPCTOFPt(nullptr),
  fHistAntiProtonPt(nullptr),
  fHistAntiProtonEta(nullptr),
  fHistAntiProtonPhi(nullptr),
  fHistAntiProtonDcaXY(nullptr),
  fHistAntiProtonDcaZ(nullptr),
  fHist2AntiProtonSigTPC(nullptr),
  fHist2AntiProtonSigTOF(nullptr),
  fHist2AntiProtonSigRMS(nullptr),
  fHist3AntiProtonSigTPCTOFPt(nullptr),
  fResultsList(nullptr),
  fHist2Psi2TPCPosCent(nullptr),
  fHist2Psi2TPCNegCent(nullptr),
  fHist2Psi2V0CCent(nullptr),
  fHist2Psi2V0ACent(nullptr),
  fHist2Psi1ZNCCent(nullptr),
  fHist2Psi1ZNACent(nullptr),
  fHist2Psi2TPCCent(nullptr),
  fProfileTPCPsi2Correlation(nullptr),
  fProfileV0MPsi2Correlation(nullptr),
  fProfileZDCPsi1Correlation(nullptr),
  fProfileZDCPsi2Correlation(nullptr),
  fProfileV0CTPCPosPsi2Correlation(nullptr),
  fProfileV0ATPCPosPsi2Correlation(nullptr),
  fProfileV0CTPCNegPsi2Correlation(nullptr),
  fProfileV0ATPCNegPsi2Correlation(nullptr),
  fProfileV0CTPCPsi2Correlation(nullptr),
  fProfileV0ATPCPsi2Correlation(nullptr),
  fProfile2RawFlowHadronQ2(nullptr),
  fProfile2RawFlowPionQ2(nullptr),
  fProfile2RawFlowKaonQ2(nullptr),
  fProfile2RawFlowProtonQ2(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;
  for (int i = 0; i < 3; i++) pV0XMeanRead[i] = nullptr;
  for (int i = 0; i < 3; i++) pV0YMeanRead[i] = nullptr;
  for (int i = 0; i < 2; i++) hQx2mV0[i] = nullptr;
  for (int i = 0; i < 2; i++) hQy2mV0[i] = nullptr;
  for (int i = 0; i < 90; i++) splQ2c[i] = nullptr;
  for (int i = 0; i < 3; i++) vtxQuant1[i] = -999;
  for (int i = 0; i < 3; i++) vtxQuant2[i] = -999;
  for (int i = 0; i < 2; i++) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistVz[i] = nullptr;
  for (int i = 0; i < 8; i++) fHist2CentQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2MultCentQA[i] = nullptr;
  for (int i = 0; i < 6; i++) fHist2MultMultQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQxCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQyCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQxVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQyVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fHist2CalibPsi2V0CCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist3CalibQxQyCentV0C[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQxCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQyCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQxVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQyVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fHist2CalibPsi2V0ACent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist3CalibQxQyCentV0A[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileZNCTowerMeanEnegry[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileZNCQxCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileZNCQyCent[i] = nullptr;
  for (int i = 0; i < 3; i++) fHist2CalibPsi1ZNCCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileZNATowerMeanEnegry[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileZNAQxCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileZNAQyCent[i]= nullptr;
  for (int i = 0; i < 3; i++) fHist2CalibPsi1ZNACent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileZDCQxAQxCCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileZDCQxAQyCCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileZDCQyAQxCCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileZDCQyAQyCCent[i]= nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentHadron[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentPion[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentKaon[i][j] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaPionKaon[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaPionProton[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaKaonProton[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaPionPion[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaKaonKaon[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaProtonProton[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileDeltaHadronHadron[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaPionKaonSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaPionProtonSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaKaonProtonSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaPionPionSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaKaonKaonSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaProtonProtonSplit[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaHadronHadronSplit[i] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaPionKaonSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaPionProtonSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaKaonProtonSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaPionPionSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaKaonKaonSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaProtonProtonSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaHadronHadronSplit[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaPionKaon[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaPionProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaKaonProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaPionPion[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaKaonKaon[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaProtonProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfileGammaHadronHadron[i][j] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionKaonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionProtonDPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonProtonDPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionPionDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonKaonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaProtonProtonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaHadronHadronDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionKaonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionProtonSPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonProtonSPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionPionSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonKaonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaProtonProtonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaHadronHadronSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionKaonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionProtonDEta[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonProtonDEta[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionPionDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonKaonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaProtonProtonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaHadronHadronDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionKaonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionProtonQ2[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonProtonQ2[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaPionPionQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaKaonKaonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaProtonProtonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffDeltaHadronHadronQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionKaonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionProtonDPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonProtonDPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionPionDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonKaonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaProtonProtonDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaHadronHadronDPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionKaonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionProtonSPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonProtonSPt[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionPionSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonKaonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaProtonProtonSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaHadronHadronSPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionKaonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionProtonDEta[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonProtonDEta[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionPionDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonKaonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaProtonProtonDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaHadronHadronDEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionKaonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionProtonQ2[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonProtonQ2[i] = nullptr; 
  for (int i = 0; i < 2; i++) fProfile2DiffGammaPionPionQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaKaonKaonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaProtonProtonQ2[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfile2DiffGammaHadronHadronQ2[i] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiPionKaon[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiPionProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiKaonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiPionPion[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiKaonKaon[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiProtonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiHadronHadron[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiPionKaon[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiPionProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiKaonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiPionPion[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiKaonKaon[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiProtonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiHadronHadron[i][j] = nullptr;
  
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
}

//------------------------------------------------

AliAnalysisTaskPIDCME::~AliAnalysisTaskPIDCME()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  if (fQAList) delete fQAList;
  if (fResultsList) delete fResultsList;
}

//---------------------------------------------------

void AliAnalysisTaskPIDCME::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate");
}

//---------------------------------------------------

void AliAnalysisTaskPIDCME::UserCreateOutputObjects()
{
  ////////////////////////
  // Pile up Function
  ////////////////////////
  // Dobrin 15o pass2 Pile-up function
  if (fPeriod.EqualTo("LHC15o")) {
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

  // Rihan 18q/r Pile-up function
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

    Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
    fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);

    Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
    fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);

    Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
    fMultCutPU->SetParameters(parFB32);
  }

  ////////////////////////
  // NUE
  ////////////////////////
  // Load Calibration Files
  // The global read-in lists and hists are loaded here.
  // They do not need to be loaded run by run.
  if (isDoNUE) {
    if (!fListNUE) {
      std::cout<<("NUE list not found")<<std::endl;
      return;
    }
  }
  ////////////////////////
  // NUA
  ////////////////////////
  if (isDoNUA) {
    if (!fListNUA) {
      std::cout<<("NUA list not found")<<std::endl;
      return;
    }
    if (fPeriod.EqualTo("LHC10h")) {
      hNUAweightPlus  = (TH2D*)fListNUA->FindObject("hNUAweightPlus");
      hNUAweightMinus = (TH2D*)fListNUA->FindObject("hNUAweightMinus");
    }
    if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      hCorrectNUAPos = new TH3F();
      hCorrectNUANeg = new TH3F();
    }
  }
  ////////////////////////
  // VZERO
  ////////////////////////
  if (isUseVZEROPlane) {
    if (!fListVZEROCalib) {
     std::cout<<("VZERO calibration list not found")<<std::endl;
     return;
    }
    if (fPeriod.EqualTo("LHC10h")) {
      // Read GE Hists (x_y) : (iCh_runnumBin)
      hMultV0Read = (TH2D*)fListVZEROCalib->FindObject("hMultV0");
      // Read Qx/y Mean Hists (x_y_z) : (runnumBin_centBin_vzBin)
      pV0XMeanRead[1] = (TProfile3D*)fListVZEROCalib->FindObject("pV0CCosMean");
      pV0YMeanRead[1] = (TProfile3D*)fListVZEROCalib->FindObject("pV0CSinMean");
      pV0XMeanRead[2] = (TProfile3D*)fListVZEROCalib->FindObject("pV0ACosMean");
      pV0YMeanRead[2] = (TProfile3D*)fListVZEROCalib->FindObject("pV0ASinMean");
    }
    if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      contQxncm = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqxc%im",2)); // V0C Qx Mean
      contQyncm = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqyc%im",2)); // V0C Qy Mean
      contQxnam = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqxa%im",2)); // V0A Qx Mean
      contQynam = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqya%im",2)); // V0A Qy Mean
      for (int i = 0; i < 2; i++) {
        hQx2mV0[i] = new TH1D();
        hQy2mV0[i] = new TH1D();
      }
    }
    //15 V0 Mult
    if (fPeriod.EqualTo("LHC15o")) {
      contMult = (AliOADBContainer*)fListVZEROCalib->FindObject("hMultV0BefCorPfpx");
      hMultV0  = new TH1D();
      for (int isp = 0; isp < 90; isp++) splQ2c[isp] = (TSpline3*)fListVZEROCalib->FindObject(Form("sp_q2V0C_%d", isp));
    }
    if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      fHCorrectV0ChWeghts = new TH2F();
      hQnPercentile = (TH2D*)fListVZEROCalib->FindObject("h_qncPercentile");
    }
  }
  ////////////////////////
  // ZDC
  ////////////////////////
  if (isUseZDCPlane) {
    if ((fPeriod.EqualTo("LHC10h")||fPeriod.EqualTo("LHC18q")||fPeriod.EqualTo("LHC18r")) && !fListZDCCalib) {
      std::cout<<("ZDC calibration list not found")<<std::endl;
      return;
    }
    if((fPeriod.EqualTo("LHC15o"))) return;
    tree = new TTree();
    fProfileForZNCGE = new TProfile();
    fProfileForZNAGE = new TProfile();
    fHn4DForZNCQxRC = new THnSparseF();
    fHn4DForZNCQyRC = new THnSparseF();
    fHn4DForZNCMtRC = new THnSparseF();
    fHn4DForZNAQxRC = new THnSparseF();
    fHn4DForZNAQyRC = new THnSparseF();
    fHn4DForZNAMtRC = new THnSparseF();
    fHn4DForZNCCountsRC = new THnSparseI();
    fHn4DForZNACountsRC = new THnSparseI();
    fProfile2DForCosC = new TProfile2D();
    fProfile2DForSinC = new TProfile2D();
    fProfile2DForCosA = new TProfile2D();
    fProfile2DForSinA = new TProfile2D();
    fHZDCCparameters = new TH1D();
    fHZDCAparameters = new TH1D();
  }

  //------------------
  // QA
  //------------------
  fQAList = new TList();
  fQAList -> SetName("fQAList");
  fQAList -> SetOwner(kTRUE);
  // event-wise
  fEvtCount = new TH1D("EvtCount", "Event Count", 20, 1, 21);
  fEvtCount->GetXaxis()->SetBinLabel(1,"All");
  fEvtCount->GetXaxis()->SetBinLabel(2,"Manager");
  fEvtCount->GetXaxis()->SetBinLabel(3,"Handler");
  fEvtCount->GetXaxis()->SetBinLabel(4,"fAOD");
  fEvtCount->GetXaxis()->SetBinLabel(5,"fPID");
  fEvtCount->GetXaxis()->SetBinLabel(6,"fUtils");
  fEvtCount->GetXaxis()->SetBinLabel(7,"fMultSel");
  fEvtCount->GetXaxis()->SetBinLabel(8,"Trigger");
  fEvtCount->GetXaxis()->SetBinLabel(9,"Run Number");
  fEvtCount->GetXaxis()->SetBinLabel(10,"Read in");
  fEvtCount->GetXaxis()->SetBinLabel(11,"Vertex");
  fEvtCount->GetXaxis()->SetBinLabel(12,"Centrality");
  fEvtCount->GetXaxis()->SetBinLabel(13,"Pile up");
  fEvtCount->GetXaxis()->SetBinLabel(14,"Get VZERO Plane");
  fEvtCount->GetXaxis()->SetBinLabel(15,"Get ZDC Plane");
  fEvtCount->GetXaxis()->SetBinLabel(16,"Reset Vector");
  fEvtCount->GetXaxis()->SetBinLabel(17,"Loop Track");
  fEvtCount->GetXaxis()->SetBinLabel(18,"Get TPC Plane");
  fEvtCount->GetXaxis()->SetBinLabel(19,"Resolution");
  fEvtCount->GetXaxis()->SetBinLabel(20,"Pair TrkTrk");

  fQAList->Add(fEvtCount);

  ////////////////////////
  // Run Number Info
  ////////////////////////
  TString runNumList10h[90] = {
    "139510", "139507", "139505", "139503", "139465", "139438", "139437", "139360", "139329", "139328",
    "139314", "139310", "139309", "139173", "139107", "139105", "139038", "139037", "139036", "139029",
    "139028", "138872", "138871", "138870", "138837", "138732", "138730", "138666", "138662", "138653",
    "138652", "138638", "138624", "138621", "138583", "138582", "138578", "138534", "138469", "138442",
    "138439", "138438", "138396", "138364", "138275", "138225", "138201", "138197", "138192", "138190",
    "137848", "137844", "137752", "137751", "137724", "137722", "137718", "137704", "137693", "137692",
    "137691", "137686", "137685", "137639", "137638", "137608", "137595", "137549", "137546", "137544",
    "137541", "137539", "137531", "137530", "137443", "137441", "137440", "137439", "137434", "137432",
    "137431", "137430", "137243", "137236", "137235", "137232", "137231", "137230", "137162", "137161"};
  TString runNumList15o[138] = {
    "246994", "246991", "246989", "246984", "246982", "246948", "246945", "246928", "246871", "246870",
    "246867", "246865", "246864", "246859", "246858", "246851", "246847", "246846", "246845", "246844",
    "246810", "246809", "246808", "246807", "246805", "246804", "246766", "246765", "246763", "246760",
    "246759", "246758", "246757", "246751", "246750", "246434", "246431", "246424", "246392", "246391",
    "246276", "246275", "246272", "246271", "246225", "246222", "246217", "246185", "246182", "246181",
    "246180", "246178", "246153", "246152", "246151", "246148", "246115", "246113", "246089", "246087",
    "246053", "246052", "246049", "246048", "246042", "246037", "246036", "246012", "246003", "246001",
    "245963", "245954", "245952", "245949", "245923", "245833", "245831", "245829", "245793", "245785",
    "245775", "245766", "245759", "245752", "245731", "245729", "245705", "245702", "245692", "245683",
    "245554", "245545", "245544", "245543", "245542", "245540", "245535", "245507", "245505", "245504",
    "245501", "245497", "245496", "245454", "245453", "245450", "245446", "245441", "245411", "245410",
    "245409", "245407", "245401", "245397", "245396", "245353", "245349", "245347", "245346", "245345",
    "245343", "245259", "245233", "245232", "245231", "245152", "245151", "245146", "245145", "245068",
    "245066", "245064", "244983", "244982", "244980", "244975", "244918", "244917"};
  TString runNumList18q[125] = {
    "296623","296622","296621","296619","296618","296616","296615","296594","296553","296552",
    "296551","296550","296548","296547","296516","296512","296511","296510","296509","296472",
    "296433","296424","296423","296420","296419","296415","296414","296383","296381","296380",
    "296379","296378","296377","296376","296375","296312","296309","296304","296303","296280",
    "296279","296273","296270","296269","296247","296246","296244","296243","296242","296241",
    "296240","296198","296197","296196","296195","296194","296192","296191","296143","296142",
    "296135","296134","296133","296132","296123","296074","296066","296065","296063","296062",
    "296060","296016","295942","295941","295937","295936","295913","295910","295909","295861",
    "295860","295859","295856","295855","295854","295853","295831","295829","295826","295825",
    "295822","295819","295818","295816","295791","295788","295786","295763","295762","295759",
    "295758","295755","295754","295725","295723","295721","295719","295718","295717","295714",
    "295712","295676","295675","295673","295668","295667","295666","295615","295612","295611",
    "295610","295589","295588","295586","295585"};
  TString runNumList18r[89] = {
    "297595","297590","297588","297558","297544","297542","297541","297540","297537","297512",
    "297483","297479","297452","297451","297450","297446","297442","297441","297415","297414",
    "297413","297406","297405","297380","297379","297372","297367","297366","297363","297336",
    "297335","297333","297332","297317","297311","297310","297278","297222","297221","297218",
    "297196","297195","297193","297133","297132","297129","297128","297124","297123","297119",
    "297118","297117","297085","297035","297031","296966","296941","296938","296935","296934",
    "296932","296931","296930","296903","296900","296899","296894","296852","296851","296850",
    "296848","296839","296838","296836","296835","296799","296794","296793","296790","296787",
    "296786","296785","296784","296781","296752","296694","296693","296691","296690"};
  runNumList = new std::map<int,int>;
  if      (fPeriod.EqualTo("LHC10h")) for (int i = 0; i < 90; i++) runNumList->insert(std::pair<int,int>(runNumList10h[i].Atoi(),i+1));
  else if (fPeriod.EqualTo("LHC15o")) for (int i = 0; i <138; i++) runNumList->insert(std::pair<int,int>(runNumList15o[i].Atoi(),i+1));
  else if (fPeriod.EqualTo("LHC18q")) for (int i = 0; i <125; i++) runNumList->insert(std::pair<int,int>(runNumList18q[i].Atoi(),i+1));
  else if (fPeriod.EqualTo("LHC18r")) for (int i = 0; i < 89; i++) runNumList->insert(std::pair<int,int>(runNumList18r[i].Atoi(),i+1));
  else return;
  fHistRunNumBin = new TH1I("runNumBin","",(int)runNumList->size(),1,(int)runNumList->size()+1);
  std::map<int,int>::iterator iter;
  for (auto runNum : *runNumList) fHistRunNumBin->GetXaxis()->SetBinLabel(runNum.second, Form("%i",runNum.first));
  fQAList->Add(fHistRunNumBin);

  // Event-wise QA
  fHistCent[0] = new TH1D("fHistCentBfCut", "Dist. of Centrality Before Cut", 80, 0., 80.);
  fHistCent[1] = new TH1D("fHistCentAfCut", "Dist. of Centrality After Cut", 80, 0., 80.);
  fHistVz[0] = new TH1D("fHistVzBfCut", "Dist of Vz Before Cut", 200, -20., 20.);
  fHistVz[1] = new TH1D("fHistVzAfCut", "Dist of Vz After Cut", 200, -20., 20.);
  fHist2CentQA[0] = new TH2D("fHist2CentQA_V0M_SPD1_BfCut", ";centV0M;centSPD1", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[1] = new TH2D("fHist2CentQA_V0M_SPD1_AfCut", ";centV0M;centSPD1", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[2] = new TH2D("fHist2CentQA_V0M_TRK_BfCut", ";centV0M;centTRK", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[3] = new TH2D("fHist2CentQA_V0M_TRK_AfCut", ";centV0M;centTRK", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[4] = new TH2D("fHist2CentQA_V0M_SPD0_BfCut", ";centV0M;centSPD0", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[5] = new TH2D("fHist2CentQA_V0M_SPD0_AfCut", ";centV0M;centSPD0", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[6] = new TH2D("fHist2CentQA_SPD1_SPD0_BfCut", ";centSPD1;centSPD0", 80, 0, 80., 80, 0, 80.);
  fHist2CentQA[7] = new TH2D("fHist2CentQA_SPD1_SPD0_AfCut", ";centSPD1;centSPD0", 80, 0, 80., 80, 0, 80.);
  fHist2MultCentQA[0] = new TH2D("fHist2MultCentQA_BfCut", ";centV0M;multFB32", 80, 0, 80., 20, 0, 5000.);
  fHist2MultCentQA[1] = new TH2D("fHist2MultCentQA_AfCut", ";centV0M;multFB32", 80, 0, 80., 20, 0, 5000.);
  fQAList->Add(fHistCent[0]);
  fQAList->Add(fHistCent[1]);
  fQAList->Add(fHistVz[0]);
  fQAList->Add(fHistVz[1]);
  for (int i = 0; i < 8; i++) fQAList->Add(fHist2CentQA[i]);
  fQAList->Add(fHist2MultCentQA[0]);
  fQAList->Add(fHist2MultCentQA[1]);

  // if (fMultComp.EqualTo("pileupByEDSTPC128") || fMultComp.EqualTo("pileupByGlobalTPC1")) {
  //   hMultMultQA[0] = new TH2D("hMultMultQAmTPCmESDPileupBefCut", "befCut;multTPC;multESD", 50, 0, 5000, 160, 0, 16000);
  //   hMultMultQA[1] = new TH2D("hMultMultQAmClobalmTPCEFPileupBefCut", "befCut;multGlobal;multTPCFE", 50, 0, 5000, 50, 0, 5000);
  //   hMultMultQA[2] = new TH2D("hMultMultQAmFB32mTrkTOFBefCut", "befCut;multTrkTOF;nTrk", 201, 0, 20000, 201, 0, 20000);
  //   hMultMultQA[3] = new TH2D("hMultMultQAmTPCmESDPileupAftCut", "aftCut;multTPC;multESD", 50, 0, 5000, 160, 0, 16000);
  //   hMultMultQA[4] = new TH2D("hMultMultQAmClobalmTPCEFPileupAftCut", "aftCut;multGlobal;multTPCFE", 50, 0, 5000, 50, 0, 5000);
  //   hMultMultQA[5] = new TH2D("hMultMultQAmFB32mTrkTOFAftCut", "aftCut;multTrkTOF;nTrk", 201, 0, 20000, 201, 0, 20000);
  //   fOutputList->Add(hMultMultQA[0]);
  //   fOutputList->Add(hMultMultQA[1]);
  //   fOutputList->Add(hMultMultQA[2]);
  //   fOutputList->Add(hMultMultQA[3]);
  //   fOutputList->Add(hMultMultQA[4]);
  //   fOutputList->Add(hMultMultQA[5]);
  // }

  // track-wise QA
  fHistPt  = new TH1D("fHistPt", ";p_{T}", 200, 0., 20.);
  fHistEta = new TH1D("fHistEta", ";#eta", 100, -2., 2.);
  fHistNhits = new TH1D("fHistNhits", ";nhits", 200, 0., 200.);
  fHist2PDedx = new TH2D("fHist2PDedx", ";pdedx", 400, -10., 10., 400, 0, 1000);
  fHistDcaXY = new TH1D("fHistDcaXY", ";DcaXY", 500, 0., 5);
  fHistDcaZ  = new TH1D("fHistDcaZ",  ";DcaZ", 500, 0., 5);
  fHistPhi[0] = new TH1D("fHistPhi", ";#phi", 100, 0, TMath::TwoPi());
  fHistPhi[1] = new TH1D("fHistPhi_afterNUA", ";#phi", 100, 0, TMath::TwoPi());
  fQAList->Add(fHistPt);
  fQAList->Add(fHistEta);
  fQAList->Add(fHistNhits);
  fQAList->Add(fHist2PDedx);
  fQAList->Add(fHistDcaXY);
  fQAList->Add(fHistDcaZ);
  for (int i = 0; i < 2; i++) fQAList->Add(fHistPhi[i]);

  //Plane QA
  std::string charCalibStep;
  for (int i = 0; i < 2; i++) {
    if (i==0) charCalibStep = "GE";
    if (i==1) charCalibStep = "RC";
    if (isQAVZERO) {
      fProfileV0CQxCent[i] = new TProfile(Form("fProfileV0CQxCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileV0CQyCent[i] = new TProfile(Form("fProfileV0CQyCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileV0CQxVtx[i]  = new TProfile(Form("fProfileV0CQxVz%s",charCalibStep.data()), "", 20, -10, 10);
      fProfileV0CQyVtx[i]  = new TProfile(Form("fProfileV0CQyVz%s",charCalibStep.data()), "", 20, -10, 10);
      fHist2CalibPsi2V0CCent[i] = new TH2D(Form("fHist2CalibPsi2V0CCent%s",charCalibStep.data()), "", 80, 0, 80, 100, 0, TMath::Pi());
      fHist3CalibQxQyCentV0C[i] = new TH3D(Form("fHist3CalibQxQyCentV0C%s",charCalibStep.data()), "fHist3CalibQxQyCentV0C;Q_{x};Q_{y};centrality", 100, -10, 100, 10, -10, 10, 80, 0, 80);
      fQAList->Add(fProfileV0CQxCent[i]);
      fQAList->Add(fProfileV0CQyCent[i]);
      fQAList->Add(fProfileV0CQxVtx[i]);
      fQAList->Add(fProfileV0CQyVtx[i]);
      fQAList->Add(fHist2CalibPsi2V0CCent[i]);
      fQAList->Add(fHist3CalibQxQyCentV0C[i]);

      fProfileV0AQxCent[i] = new TProfile(Form("fProfileV0AQxCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileV0AQyCent[i] = new TProfile(Form("fProfileV0AQyCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileV0AQxVtx[i]  = new TProfile(Form("fProfileV0AQxVz%s",charCalibStep.data()), "", 20, -10, 10);
      fProfileV0AQyVtx[i]  = new TProfile(Form("fProfileV0AQyVz%s",charCalibStep.data()), "", 20, -10, 10);
      fHist2CalibPsi2V0ACent[i] = new TH2D(Form("fHist2CalibPsi2V0ACent%s",charCalibStep.data()), "", 80, 0, 80, 100, 0, TMath::Pi());
      fHist3CalibQxQyCentV0A[i] = new TH3D(Form("fHist3CalibQxQyCentV0A%s",charCalibStep.data()), "fHist3CalibQxQyCentV0A;Q_{x};Q_{y};centrality", 100, -10, 10, 100, -10, 10, 80, 0, 80);
      fQAList->Add(fProfileV0AQxCent[i]);
      fQAList->Add(fProfileV0AQyCent[i]);
      fQAList->Add(fProfileV0AQxVtx[i]);
      fQAList->Add(fProfileV0AQyVtx[i]);
      fQAList->Add(fHist2CalibPsi2V0ACent[i]);
      fQAList->Add(fHist3CalibQxQyCentV0A[i]);
    }

    if (isQAZDC) {
      fProfileZNCQxCent[i] = new TProfile(Form("fProfileZNCQxCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileZNCQyCent[i] = new TProfile(Form("fProfileZNCQyCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fHist2CalibPsi1ZNCCent[i] = new TH2D(Form("fHist2CalibPsi1ZNCCent%s",charCalibStep.data()), "", 80, 0, 80, 50, 0, TMath::Pi());
      fQAList->Add(fProfileZNCQxCent[i]);
      fQAList->Add(fProfileZNCQyCent[i]);
      fQAList->Add(fHist2CalibPsi1ZNCCent[i]);

      fProfileZNAQxCent[i] = new TProfile(Form("fProfileZNAQxCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileZNAQyCent[i] = new TProfile(Form("fProfileZNAQyCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fHist2CalibPsi1ZNACent[i] = new TH2D(Form("fHist2CalibPsi1ZNACent%s",charCalibStep.data()), "", 80, 0, 80, 50, 0, TMath::Pi());
      fQAList->Add(fProfileZNAQxCent[i]);
      fQAList->Add(fProfileZNAQyCent[i]);
      fQAList->Add(fHist2CalibPsi1ZNACent[i]);

      fProfileZDCQxAQxCCent[i] = new TProfile(Form("fProfileZDCQxAQxCCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileZDCQxAQyCCent[i] = new TProfile(Form("fProfileZDCQxAQyCCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileZDCQyAQxCCent[i] = new TProfile(Form("fProfileZDCQyAQxCCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileZDCQyAQyCCent[i] = new TProfile(Form("fProfileZDCQyAQyCCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fQAList->Add(fProfileZDCQxAQxCCent[i]);
      fQAList->Add(fProfileZDCQxAQyCCent[i]);
      fQAList->Add(fProfileZDCQyAQxCCent[i]);
      fQAList->Add(fProfileZDCQyAQyCCent[i]);
    }
  }

  if (isQAZDC) {
    fHist2CalibPsi1ZNCCent[2] = new TH2D("fHist2CalibPsi1ZNCCentSF", "", 16, 0, 80, 50, 0, TMath::Pi());
    fHist2CalibPsi1ZNACent[2] = new TH2D("fHist2CalibPsi1ZNACentSF", "", 16, 0, 80, 50, 0, TMath::Pi());
    fQAList->Add(fHist2CalibPsi1ZNCCent[2]);
    fQAList->Add(fHist2CalibPsi1ZNACent[2]);

    fProfileZNCTowerMeanEnegry[0] = new TProfile("fProfileZNCTowerMeanEnegryRW","",5,0,5);
    fProfileZNCTowerMeanEnegry[1] = new TProfile("fProfileZNCTowerMeanEnegryGE","",5,0,5);
    fProfileZNATowerMeanEnegry[0] = new TProfile("fProfileZNATowerMeanEnegryRW","",5,0,5);
    fProfileZNATowerMeanEnegry[1] = new TProfile("fProfileZNATowerMeanEnegryGE","",5,0,5);
    fQAList->Add(fProfileZNCTowerMeanEnegry[0]);
    fQAList->Add(fProfileZNCTowerMeanEnegry[1]);
    fQAList->Add(fProfileZNATowerMeanEnegry[0]);
    fQAList->Add(fProfileZNATowerMeanEnegry[1]);
  }
  //Pion QA
  fHistPionPt      = new TH1D("fHistPionPt"     , "fHistPionPt;p_{T}"      , 200, 0., 20.);
  fHistPionEta     = new TH1D("fHistPionEta"    , "fHistPionEta;#eta"      , 100, -2., 2.);
  fHistPionPhi     = new TH1D("fHistPionPhi"    , "fHistPionPhi;#phi"      , 360, 0., TMath::TwoPi());
  fHistPionDcaXY   = new TH1D("fHistPionDcaXY"  , "fHistPionDcaXY;DcaXY"   , 500, 0., 5);
  fHistPionDcaZ    = new TH1D("fHistPionDcaZ "  , "fHistPionDcaZ;DcaZ"     , 500, 0., 5);
  fHist2PionSigTPC = new TH2D("fHist2PionSigTPC", "fHist2PionSigTPC;p_{T};SigTPC", 25, 0., 5., 100, -5., 5.);
  fHist2PionSigTOF = new TH2D("fHist2PionSigTOF", "fHist2PionSigTOF;p_{T};SigTOF", 25, 0., 5., 100, -5., 5.);
  fHist2PionSigRMS = new TH2D("fHist2PionSigRMS", "fHist2PionSigRMS;p_{T};SigRMS", 25, 0., 5., 100, 0., 10.);
  fHist3PionSigTPCTOFPt = new TH3D("fHist3PionSigTPCTOFPt","fHist3PionSigTPCTOFPt;SigTPC;SigTOF;p_{T}", 100, -5., 5., 100, -5., 5., 20, 0., 10.);
  fQAList->Add(fHistPionPt);
  fQAList->Add(fHistPionEta);
  fQAList->Add(fHistPionPhi);
  fQAList->Add(fHistPionDcaXY);
  fQAList->Add(fHistPionDcaZ);
  fQAList->Add(fHist2PionSigTPC);
  fQAList->Add(fHist2PionSigTOF);
  fQAList->Add(fHist2PionSigRMS);
  fQAList->Add(fHist3PionSigTPCTOFPt);
  fHistAntiPionPt      = new TH1D("fHistAntiPionPt"     , "fHistAntiPionPt;p_{T}"      , 200, 0., 20.);
  fHistAntiPionEta     = new TH1D("fHistAntiPionEta"    , "fHistAntiPionEta;#eta"      , 100, -2., 2.);
  fHistAntiPionPhi     = new TH1D("fHistAntiPionPhi"    , "fHistAntiPionPhi;#phi"      , 360, 0., TMath::TwoPi());
  fHistAntiPionDcaXY   = new TH1D("fHistAntiPionDcaXY"  , "fHistAntiPionDcaXY;DcaXY"   , 500, 0., 5);
  fHistAntiPionDcaZ    = new TH1D("fHistAntiPionDcaZ "  , "fHistAntiPionDcaZ;DcaZ"     , 500, 0., 5);
  fHist2AntiPionSigTPC = new TH2D("fHist2AntiPionSigTPC", "fHist2AntiPionSigTPC;p_{T};SigTPC", 25, 0., 5., 100, -5., 5.);
  fHist2AntiPionSigTOF = new TH2D("fHist2AntiPionSigTOF", "fHist2AntiPionSigTOF;p_{T};SigTOF", 25, 0., 5., 100, -5., 5.);
  fHist2AntiPionSigRMS = new TH2D("fHist2AntiPionSigRMS", "fHist2AntiPionSigRMS;p_{T};SigRMS", 25, 0., 5., 100, 0., 10.);
  fHist3AntiPionSigTPCTOFPt = new TH3D("fHist3AntiPionSigTPCTOFPt","fHist3AntiPionSigTPCTOFPt;SigTPC;SigTOF;p_{T}", 100, -5., 5., 100, -5., 5., 20, 0., 10.);
  fQAList->Add(fHistAntiPionPt);
  fQAList->Add(fHistAntiPionEta);
  fQAList->Add(fHistAntiPionPhi);
  fQAList->Add(fHistAntiPionDcaXY);
  fQAList->Add(fHistAntiPionDcaZ);
  fQAList->Add(fHist2AntiPionSigTPC);
  fQAList->Add(fHist2AntiPionSigTOF);
  fQAList->Add(fHist2AntiPionSigRMS);
  fQAList->Add(fHist3AntiPionSigTPCTOFPt);

  //Kaon QA
  fHistKaonPt      = new TH1D("fHistKaonPt"     , "fHistKaonPt;p_{T}"      , 200, 0., 20.);
  fHistKaonEta     = new TH1D("fHistKaonEta"    , "fHistKaonEta;#eta"      , 100, -2., 2.);
  fHistKaonPhi     = new TH1D("fHistKaonPhi"    , "fHistKaonPhi;#phi"      , 360, 0., TMath::TwoPi());
  fHistKaonDcaXY   = new TH1D("fHistKaonDcaXY"  , "fHistKaonDcaXY;DcaXY"   , 500, 0., 5);
  fHistKaonDcaZ    = new TH1D("fHistKaonDcaZ "  , "fHistKaonDcaZ;DcaZ"     , 500, 0., 5);
  fHist2KaonSigTPC = new TH2D("fHist2KaonSigTPC", "fHist2KaonSigTPC;p_{T};SigTPC", 25, 0., 5., 100, -5., 5.);
  fHist2KaonSigTOF = new TH2D("fHist2KaonSigTOF", "fHist2KaonSigTOF;p_{T};SigTOF", 25, 0., 5., 100, -5., 5.);
  fHist2KaonSigRMS = new TH2D("fHist2KaonSigRMS", "fHist2KaonSigRMS;p_{T};SigRMS", 25, 0., 5., 100, 0., 10.);
  fHist3KaonSigTPCTOFPt = new TH3D("fHist3KaonSigTPCTOFPt","fHist3KaonSigTPCTOFPt;SigTPC;SigTOF;p_{T}", 100, -5., 5., 100, -5., 5., 20, 0., 10.);
  fQAList->Add(fHistKaonPt);
  fQAList->Add(fHistKaonEta);
  fQAList->Add(fHistKaonPhi);
  fQAList->Add(fHistKaonDcaXY);
  fQAList->Add(fHistKaonDcaZ);
  fQAList->Add(fHist2KaonSigTPC);
  fQAList->Add(fHist2KaonSigTOF);
  fQAList->Add(fHist2KaonSigRMS);
  fQAList->Add(fHist3KaonSigTPCTOFPt);
  fHistAntiKaonPt      = new TH1D("fHistAntiKaonPt"     , "fHistAntiKaonPt;p_{T}"      , 200, 0., 20.);
  fHistAntiKaonEta     = new TH1D("fHistAntiKaonEta"    , "fHistAntiKaonEta;#eta"      , 100, -2., 2.);
  fHistAntiKaonPhi     = new TH1D("fHistAntiKaonPhi"    , "fHistAntiKaonPhi;#phi"      , 360, 0., TMath::TwoPi());
  fHistAntiKaonDcaXY   = new TH1D("fHistAntiKaonDcaXY"  , "fHistAntiKaonDcaXY;DcaXY"   , 500, 0., 5);
  fHistAntiKaonDcaZ    = new TH1D("fHistAntiKaonDcaZ "  , "fHistAntiKaonDcaZ;DcaZ"     , 500, 0., 5);
  fHist2AntiKaonSigTPC = new TH2D("fHist2AntiKaonSigTPC", "fHist2AntiKaonSigTPC;p_{T};SigTPC", 25, 0., 5., 100, -5., 5.);
  fHist2AntiKaonSigTOF = new TH2D("fHist2AntiKaonSigTOF", "fHist2AntiKaonSigTOF;p_{T};SigTOF", 25, 0., 5., 100, -5., 5.);
  fHist2AntiKaonSigRMS = new TH2D("fHist2AntiKaonSigRMS", "fHist2AntiKaonSigRMS;p_{T};SigRMS", 25, 0., 5., 100, 0., 10.);
  fHist3AntiKaonSigTPCTOFPt = new TH3D("fHist3AntiKaonSigTPCTOFPt","fHist3AntiKaonSigTPCTOFPt;SigTPC;SigTOF;p_{T}", 100, -5., 5., 100, -5., 5., 20, 0., 10.);
  fQAList->Add(fHistAntiKaonPt);
  fQAList->Add(fHistAntiKaonEta);
  fQAList->Add(fHistAntiKaonPhi);
  fQAList->Add(fHistAntiKaonDcaXY);
  fQAList->Add(fHistAntiKaonDcaZ);
  fQAList->Add(fHist2AntiKaonSigTPC);
  fQAList->Add(fHist2AntiKaonSigTOF);
  fQAList->Add(fHist2AntiKaonSigRMS);
  fQAList->Add(fHist3AntiKaonSigTPCTOFPt);

  //Proton QA
  fHistProtonPt      = new TH1D("fHistProtonPt"     , "fHistProtonPt;p_{T}"      , 200, 0., 20.);
  fHistProtonEta     = new TH1D("fHistProtonEta"    , "fHistProtonEta;#eta"      , 100, -2., 2.);
  fHistProtonPhi     = new TH1D("fHistProtonPhi"    , "fHistProtonPhi;#phi"      , 360, 0., TMath::TwoPi());
  fHistProtonDcaXY   = new TH1D("fHistProtonDcaXY"  , "fHistProtonDcaXY;DcaXY"   , 500, 0., 5);
  fHistProtonDcaZ    = new TH1D("fHistProtonDcaZ "  , "fHistProtonDcaZ;DcaZ"     , 500, 0., 5);
  fHist2ProtonSigTPC = new TH2D("fHist2ProtonSigTPC", "fHist2ProtonSigTPC;p_{T};SigTPC", 25, 0., 5., 100, -5., 5.);
  fHist2ProtonSigTOF = new TH2D("fHist2ProtonSigTOF", "fHist2ProtonSigTOF;p_{T};SigTOF", 25, 0., 5., 100, -5., 5.);
  fHist2ProtonSigRMS = new TH2D("fHist2ProtonSigRMS", "fHist2ProtonSigRMS;p_{T};SigRMS", 25, 0., 5., 100, 0., 10.);
  fHist3ProtonSigTPCTOFPt = new TH3D("fHist3ProtonSigTPCTOFPt","fHist3ProtonSigTPCTOFPt;SigTPC;SigTOF;p_{T}", 100, -5., 5., 100, -5., 5., 20, 0., 10.);
  fQAList->Add(fHistProtonPt);
  fQAList->Add(fHistProtonEta);
  fQAList->Add(fHistProtonPhi);
  fQAList->Add(fHistProtonDcaXY);
  fQAList->Add(fHistProtonDcaZ);
  fQAList->Add(fHist2ProtonSigTPC);
  fQAList->Add(fHist2ProtonSigTOF);
  fQAList->Add(fHist2ProtonSigRMS);
  fQAList->Add(fHist3ProtonSigTPCTOFPt);
  fHistAntiProtonPt      = new TH1D("fHistAntiProtonPt"     , "fHistAntiProtonPt;p_{T}"      , 200, 0., 20.); 
  fHistAntiProtonEta     = new TH1D("fHistAntiProtonEta"    , "fHistAntiProtonEta;#eta"      , 100, -2., 2.); 
  fHistAntiProtonPhi     = new TH1D("fHistAntiProtonPhi"    , "fHistAntiProtonPhi;#phi"      , 360, 0., TMath::TwoPi()); 
  fHistAntiProtonDcaXY   = new TH1D("fHistAntiProtonDcaXY"  , "fHistAntiProtonDcaXY;DcaXY"   , 500, 0., 5); 
  fHistAntiProtonDcaZ    = new TH1D("fHistAntiProtonDcaZ "  , "fHistAntiProtonDcaZ;DcaZ"     , 500, 0., 5); 
  fHist2AntiProtonSigTPC = new TH2D("fHist2AntiProtonSigTPC", "fHist2AntiProtonSigTPC;p_{T};SigTPC", 25, 0., 5., 100, -5., 5.); 
  fHist2AntiProtonSigTOF = new TH2D("fHist2AntiProtonSigTOF", "fHist2AntiProtonSigTOF;p_{T};SigTOF", 25, 0., 5., 100, -5., 5.); 
  fHist2AntiProtonSigRMS = new TH2D("fHist2AntiProtonSigRMS", "fHist2AntiProtonSigRMS;p_{T};SigRMS", 25, 0., 5., 100, 0., 10.);
  fHist3AntiProtonSigTPCTOFPt = new TH3D("fHist3AntiProtonSigTPCTOFPt","fHist3AntiProtonSigTPCTOFPt;SigTPC;SigTOF;p_{T}", 100, -5., 5., 100, -5., 5., 20, 0., 10.);
  fQAList->Add(fHistAntiProtonPt);
  fQAList->Add(fHistAntiProtonEta);
  fQAList->Add(fHistAntiProtonPhi);
  fQAList->Add(fHistAntiProtonDcaXY);
  fQAList->Add(fHistAntiProtonDcaZ);
  fQAList->Add(fHist2AntiProtonSigTPC);
  fQAList->Add(fHist2AntiProtonSigTOF);
  fQAList->Add(fHist2AntiProtonSigRMS);
  fQAList->Add(fHist3AntiProtonSigTPCTOFPt);

  PostData(1,fQAList);
  if (fDebug) Printf("Post fQAList Data Success!");

  ////////////////////////
  // Results
  ////////////////////////
  fResultsList = new TList();
  fResultsList -> SetName("fResultsList");
  fResultsList -> SetOwner(kTRUE);

  // Plane
  fHist2Psi2TPCPosCent = new TH2D("fHist2Psi2TPCPosCent","fHist2Psi2TPCPosCent;centrality;#Psi(TPCPos)",80,0,80.,100,0.,TMath::Pi());
  fHist2Psi2TPCNegCent = new TH2D("fHist2Psi2TPCNegCent","fHist2Psi2TPCNegCent;centrality;#Psi(TPCNeg)",80,0,80.,100,0.,TMath::Pi());
  fHist2Psi2V0CCent    = new TH2D("fHist2Psi2V0CCent","fHist2Psi2V0CCent;centrality;#Psi(V0C)",80,0.,80.,100,0.,TMath::Pi());
  fHist2Psi2V0ACent    = new TH2D("fHist2Psi2V0ACent","fHist2Psi2V0ACent;centrality;#Psi(V0A)",80,0.,80.,100,0.,TMath::Pi());
  fHist2Psi1ZNCCent    = new TH2D("fHist2Psi1ZNCCent","fHist2Psi1ZNCCent;centrality;#Psi(ZNC)",80,0.,80.,100,0.,TMath::TwoPi());
  fHist2Psi1ZNACent    = new TH2D("fHist2Psi1ZNACent","fHist2Psi1ZNACent;centrality;#Psi(ZNA)",80,0.,80.,100,0.,TMath::TwoPi());
  fHist2Psi2TPCCent    = new TH2D("fHist2Psi2TPCCent","fHist2Psi2TPCCent;centrality;#Psi(TPC)",80,0.,80.,100,0.,TMath::Pi());
  
  fResultsList->Add(fHist2Psi2TPCPosCent);
  fResultsList->Add(fHist2Psi2TPCNegCent);
  fResultsList->Add(fHist2Psi2V0CCent);
  fResultsList->Add(fHist2Psi2V0ACent);
  fResultsList->Add(fHist2Psi1ZNCCent);
  fResultsList->Add(fHist2Psi1ZNACent);
  fResultsList->Add(fHist2Psi2TPCCent);

  // Res
  fProfileTPCPsi2Correlation       = new TProfile("fProfileTPCPsi2Correlation","fProfileTPCPsi2Correlation;centrality;Res",80,0.,80.);
  fProfileV0MPsi2Correlation       = new TProfile("fProfileV0MPsi2Correlation","fProfileV0MPsi2Correlation;centrality;Res",80,0.,80.);
  fProfileZDCPsi1Correlation       = new TProfile("fProfileZDCPsi1Correlation","fProfileZDCPsi1Correlation;centrality;Res",80,0.,80.);
  fProfileZDCPsi2Correlation       = new TProfile("fProfileZDCPsi2Correlation","fProfileZDCPsi2Correlation;centrality;Res",80,0.,80.);
  fProfileV0CTPCPosPsi2Correlation = new TProfile("fProfileV0CTPCPosPsi2Correlation","fProfileV0CTPCPosPsi2Correlation;centrality;Res",80,0.,80.);
  fProfileV0ATPCPosPsi2Correlation = new TProfile("fProfileV0ATPCPosPsi2Correlation","fProfileV0ATPCPosPsi2Correlation;centrality;Res",80,0.,80.);
  fProfileV0CTPCNegPsi2Correlation = new TProfile("fProfileV0CTPCNegPsi2Correlation","fProfileV0CTPCNegPsi2Correlation;centrality;Res",80,0.,80.);
  fProfileV0ATPCNegPsi2Correlation = new TProfile("fProfileV0ATPCNegPsi2Correlation","fProfileV0ATPCNegPsi2Correlation;centrality;Res",80,0.,80.);
  fProfileV0CTPCPsi2Correlation    = new TProfile("fProfileV0CTPCPsi2Correlation","fProfileV0CTPCPsi2Correlation;centrality;Res",80,0.,80.);
  fProfileV0ATPCPsi2Correlation    = new TProfile("fProfileV0ATPCPsi2Correlation","fProfileV0ATPCPsi2Correlation;centrality;Res",80,0.,80.);
  fResultsList->Add(fProfileTPCPsi2Correlation);
  fResultsList->Add(fProfileV0MPsi2Correlation);
  fResultsList->Add(fProfileZDCPsi1Correlation);
  fResultsList->Add(fProfileZDCPsi2Correlation);
  fResultsList->Add(fProfileV0CTPCPosPsi2Correlation);
  fResultsList->Add(fProfileV0ATPCPosPsi2Correlation);
  fResultsList->Add(fProfileV0CTPCNegPsi2Correlation);
  fResultsList->Add(fProfileV0ATPCNegPsi2Correlation);
  fResultsList->Add(fProfileV0CTPCPsi2Correlation);
  fResultsList->Add(fProfileV0ATPCPsi2Correlation);

  // Flow
  //[0]TPC [1]V0C [2]V0A [3]ZNC [4]ZNA
  //[0]pos [1]neg(anti)

  std::string charPlane;
  std::string charCharge;
  if (isCalculatePIDFlow) {
    for (int iPlane = 0; iPlane < 5; iPlane++) {
      if (iPlane == 0) charPlane = "_TPC";
      if (iPlane == 1) charPlane = "_V0C";
      if (iPlane == 2) charPlane = "_V0A";
      if (iPlane == 3) charPlane = "_ZNC";
      if (iPlane == 4) charPlane = "_ZNA";
      if (!isUseVZEROPlane && (iPlane == 1 || iPlane == 2)) continue;
      if (!isUseZDCPlane   && (iPlane == 3 || iPlane == 4)) continue;
      for (int jCharge = 0; jCharge < 2; jCharge++) {
        if (jCharge == 0) charCharge = "_Pos";
        if (jCharge == 1) charCharge = "_Neg";
          fProfile2RawFlowPtCentHadron[iPlane][jCharge] = new TProfile2D(Form("fProfile2RawFlowPtCentHadron%s%s",charPlane.data(),charCharge.data()),";centrality;pT;v2Obs",16,0.,80.,100,0.,10.);
          fProfile2RawFlowPtCentProton[iPlane][jCharge] = new TProfile2D(Form("fProfile2RawFlowPtCentProton%s%s",charPlane.data(),charCharge.data()),";centrality;pT;v2Obs",16,0.,80.,100,0.,10.);
          fProfile2RawFlowPtCentPion[iPlane][jCharge] = new TProfile2D(Form("fProfile2RawFlowPtCentPion%s%s",charPlane.data(),charCharge.data()),";centrality;pT;v2Obs",16,0.,80.,100,0.,10.);
          fProfile2RawFlowPtCentKaon[iPlane][jCharge] = new TProfile2D(Form("fProfile2RawFlowPtCentKaon%s%s",charPlane.data(),charCharge.data()),";centrality;pT;v2Obs",16,0.,80.,100,0.,10.);
          fResultsList->Add(fProfile2RawFlowPtCentHadron[iPlane][jCharge]);
          fResultsList->Add(fProfile2RawFlowPtCentProton[iPlane][jCharge]);
          fResultsList->Add(fProfile2RawFlowPtCentPion[iPlane][jCharge]);
          fResultsList->Add(fProfile2RawFlowPtCentKaon[iPlane][jCharge]);
      }
    }
  }
  fProfile2RawFlowHadronQ2 = new TProfile2D("fProfile2RawFlowHadronQ2",";centrality;q^{2}_{V0C};v2_{Hadron}", 16,0.,80, fNQ2Bins,0.,(double)fNQ2Bins);
  fProfile2RawFlowProtonQ2 = new TProfile2D("fProfile2RawFlowProtonQ2",";centrality;q^{2}_{V0C};v2_{Proton}", 16,0.,80, fNQ2Bins,0.,(double)fNQ2Bins);
  fProfile2RawFlowPionQ2 = new TProfile2D("fProfile2RawFlowPionQ2",";centrality;q^{2}_{V0C};v2_{Pion}", 16,0.,80, fNQ2Bins,0.,(double)fNQ2Bins);
  fProfile2RawFlowKaonQ2 = new TProfile2D("fProfile2RawFlowKaonQ2",";centrality;q^{2}_{V0C};v2_{Kaon}", 16,0.,80, fNQ2Bins,0.,(double)fNQ2Bins);
  fResultsList->Add(fProfile2RawFlowHadronQ2);
  fResultsList->Add(fProfile2RawFlowProtonQ2);
  fResultsList->Add(fProfile2RawFlowPionQ2);
  fResultsList->Add(fProfile2RawFlowKaonQ2);

  std::string Chargetype;
  for (int i = 0; i < 2; i++) {
    if(i == 0) Chargetype = "SS";
    if(i == 1) Chargetype = "OS";
    if (isCalculatePionKaon) {
      fProfileDeltaPionKaon[i] = new TProfile(Form("fProfileDeltaPionKaon_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
      fResultsList->Add(fProfileDeltaPionKaon[i]);
    }
    if (isCalculatePionProton) {
      fProfileDeltaPionProton[i] = new TProfile(Form("fProfileDeltaPionProton_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
      fResultsList->Add(fProfileDeltaPionProton[i]);
    }
    if (isCalculateKaonProton) {
      fProfileDeltaKaonProton[i] = new TProfile(Form("fProfileDeltaKaonProton_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
      fResultsList->Add(fProfileDeltaKaonProton[i]);
    }
    if (isCalculatePionPion) {
      fProfileDeltaPionPion[i] = new TProfile(Form("fProfileDeltaPionPion_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
      fResultsList->Add(fProfileDeltaPionPion[i]);
    }
    if (isCalculateKaonKaon) {
      fProfileDeltaKaonKaon[i] = new TProfile(Form("fProfileDeltaKaonKaon_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
      fResultsList->Add(fProfileDeltaKaonKaon[i]);
    }
    if (isCalculateProtonProton) {
      fProfileDeltaProtonProton[i] = new TProfile(Form("fProfileDeltaProtonProton_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
      fResultsList->Add(fProfileDeltaProtonProton[i]);
    }
    if (isCalculateHadronHadron) {
      fProfileDeltaHadronHadron[i] = new TProfile(Form("fProfileDeltaHadronHadron_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
      fResultsList->Add(fProfileDeltaHadronHadron[i]);
    }
  }

  if(isOpenSsandOsSelfCheck) {
      for (int i = 0; i < 4; i++) {
      if(i == 0) Chargetype = "PP";
      if(i == 1) Chargetype = "NN";
      if(i == 2) Chargetype = "PN";
      if(i == 3) Chargetype = "NP";
      if (isCalculatePionKaon) {
        fProfileDeltaPionKaonSplit[i] = new TProfile(Form("fProfileDeltaPionKaonSplit_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
        fResultsList->Add(fProfileDeltaPionKaonSplit[i]);
      }
      if (isCalculatePionProton) {
        fProfileDeltaPionProtonSplit[i] = new TProfile(Form("fProfileDeltaPionProtonSplit_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
        fResultsList->Add(fProfileDeltaPionProtonSplit[i]);
      }
      if (isCalculateKaonProton) {
        fProfileDeltaKaonProtonSplit[i] = new TProfile(Form("fProfileDeltaKaonProtonSplit_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
        fResultsList->Add(fProfileDeltaKaonProtonSplit[i]);
      }
      if (isCalculatePionPion) {
        fProfileDeltaPionPionSplit[i] = new TProfile(Form("fProfileDeltaPionPionSplit_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
        fResultsList->Add(fProfileDeltaPionPionSplit[i]);
      }
      if (isCalculateKaonKaon) {
        fProfileDeltaKaonKaonSplit[i] = new TProfile(Form("fProfileDeltaKaonKaonSplit_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
        fResultsList->Add(fProfileDeltaKaonKaonSplit[i]);
      }
      if (isCalculateProtonProton) {
        fProfileDeltaProtonProtonSplit[i] = new TProfile(Form("fProfileDeltaProtonProtonSplit_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
        fResultsList->Add(fProfileDeltaProtonProtonSplit[i]);
      }
      if (isCalculateHadronHadron) {
        fProfileDeltaHadronHadronSplit[i] = new TProfile(Form("fProfileDeltaHadronHadronSplit_%s",Chargetype.data()),";centrality;#delta",16,0.,80.);
        fResultsList->Add(fProfileDeltaHadronHadronSplit[i]);
      }
    }
  }

  for (int iPlane = 0; iPlane < 5; iPlane++) {
    if (iPlane == 0) charPlane = "_TPC";
    if (iPlane == 1) charPlane = "_V0C";
    if (iPlane == 2) charPlane = "_V0A";
    if (iPlane == 3) charPlane = "_ZNC";
    if (iPlane == 4) charPlane = "_ZNA";
    if (!isUseVZEROPlane && (iPlane == 1 || iPlane == 2)) continue;
    if (!isUseZDCPlane   && (iPlane == 3 || iPlane == 4)) continue;
    for (int iType = 0; iType < 2; iType++) {
      if(iType == 0) Chargetype = "SS";
      if(iType == 1) Chargetype = "OS";
      if (isCalculatePionKaon) {
        fProfileGammaPionKaon[iPlane][iType] = new TProfile(Form("fProfileGammaPionKaon%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
        fResultsList->Add(fProfileGammaPionKaon[iPlane][iType]);
      }
      if (isCalculatePionProton) {
        fProfileGammaPionProton[iPlane][iType] = new TProfile(Form("fProfileGammaPionProton%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
        fResultsList->Add(fProfileGammaPionProton[iPlane][iType]);
      }
      if (isCalculateKaonProton) {
        fProfileGammaKaonProton[iPlane][iType] = new TProfile(Form("fProfileGammaKaonProton%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
        fResultsList->Add(fProfileGammaKaonProton[iPlane][iType]);
      }
      if (isCalculatePionPion) {
        fProfileGammaPionPion[iPlane][iType] = new TProfile(Form("fProfileGammaPionPion%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
        fResultsList->Add(fProfileGammaPionPion[iPlane][iType]);
      }
      if (isCalculateKaonKaon) {
        fProfileGammaKaonKaon[iPlane][iType] = new TProfile(Form("fProfileGammaKaonKaon%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
        fResultsList->Add(fProfileGammaKaonKaon[iPlane][iType]);
      }
      if (isCalculateProtonProton) {
        fProfileGammaProtonProton[iPlane][iType] = new TProfile(Form("fProfileGammaProtonProton%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
        fResultsList->Add(fProfileGammaProtonProton[iPlane][iType]);
      }
      if (isCalculateHadronHadron) {
        fProfileGammaHadronHadron[iPlane][iType] = new TProfile(Form("fProfileGammaHadronHadron%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
        fResultsList->Add(fProfileGammaHadronHadron[iPlane][iType]);
      }
    }
    if(isOpenSsandOsSelfCheck) {
        for (int iType = 0; iType < 4; iType++) {
        if(iType == 0) Chargetype = "PP";
        if(iType == 1) Chargetype = "NN";
        if(iType == 2) Chargetype = "PN";
        if(iType == 3) Chargetype = "NP";
        if (isCalculatePionKaon) {
          fProfileGammaPionKaonSplit[iPlane][iType] = new TProfile(Form("fProfileGammaPionKaonSplit%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
          fResultsList->Add(fProfileGammaPionKaonSplit[iPlane][iType]);
        }
        if (isCalculatePionProton) {
          fProfileGammaPionProtonSplit[iPlane][iType] = new TProfile(Form("fProfileGammaPionProtonSplit%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
          fResultsList->Add(fProfileGammaPionProtonSplit[iPlane][iType]);
        }
        if (isCalculateKaonProton) {
          fProfileGammaKaonProtonSplit[iPlane][iType] = new TProfile(Form("fProfileGammaKaonProtonSplit%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
          fResultsList->Add(fProfileGammaKaonProtonSplit[iPlane][iType]);
        }
        if (isCalculatePionPion) {
          fProfileGammaPionPionSplit[iPlane][iType] = new TProfile(Form("fProfileGammaPionPionSplit%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
          fResultsList->Add(fProfileGammaPionPionSplit[iPlane][iType]);
        }
        if (isCalculateKaonKaon) {
          fProfileGammaKaonKaonSplit[iPlane][iType] = new TProfile(Form("fProfileGammaKaonKaonSplit%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
          fResultsList->Add(fProfileGammaKaonKaonSplit[iPlane][iType]);
        }
        if (isCalculateProtonProton) {
          fProfileGammaProtonProtonSplit[iPlane][iType] = new TProfile(Form("fProfileGammaProtonProtonSplit%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
          fResultsList->Add(fProfileGammaProtonProtonSplit[iPlane][iType]);
        }
        if (isCalculateHadronHadron) {
          fProfileGammaHadronHadronSplit[iPlane][iType] = new TProfile(Form("fProfileGammaHadronHadronSplit%s_%s",charPlane.data(),Chargetype.data()),";centrality;#gamma",16,0.,80.);
          fResultsList->Add(fProfileGammaHadronHadronSplit[iPlane][iType]);
        }
      }
    }
  }

  if (isCalculateDiffResult) {
    // δ
    for (int iType = 0; iType < 2; iType++) {
      if(iType == 0) Chargetype = "SS";
      if(iType == 1) Chargetype = "OS";
      if (isCalculatePionKaon) {
        fProfile2DiffDeltaPionKaonDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaPionKaonDPt_%s",Chargetype.data()),";centrality;#Delta pT;#delta",16,0.,80.,16,-8,8);
        fProfile2DiffDeltaPionKaonSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaPionKaonSPt_%s",Chargetype.data()),";centrality;#Sum pT;#delta",16,0.,80.,16,0,16);
        fProfile2DiffDeltaPionKaonDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaPionKaonDEta_%s",Chargetype.data()),";centrality;#Delta Eta;#delta",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffDeltaPionKaonQ2[iType] = new TProfile2D(Form("fProfile2DiffDeltaPionKaonQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffDeltaPionKaonDPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaPionKaonSPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaPionKaonDEta[iType]);
        fResultsList->Add(fProfile2DiffDeltaPionKaonQ2[iType]);
      }
      if (isCalculatePionProton) {
        fProfile2DiffDeltaPionProtonDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaPionProtonDPt_%s",Chargetype.data()),";centrality;#Delta pT;#delta",16,0.,80.,16,-8,8);
        fProfile2DiffDeltaPionProtonSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaPionProtonSPt_%s",Chargetype.data()),";centrality;#Sum pT;#delta",16,0.,80.,16,0,16);
        fProfile2DiffDeltaPionProtonDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaPionProtonDEta_%s",Chargetype.data()),";centrality;#Delta Eta;#delta",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffDeltaPionProtonQ2[iType]   = new TProfile2D(Form("fProfile2DiffDeltaPionProtonQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffDeltaPionProtonDPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaPionProtonSPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaPionProtonDEta[iType]);
        fResultsList->Add(fProfile2DiffDeltaPionProtonQ2[iType]);
      }
      if (isCalculateKaonProton) {
        fProfile2DiffDeltaKaonProtonDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaKaonProtonDPt_%s",Chargetype.data()),";centrality;#Delta pT;#delta",16,0.,80.,16,-8,8);
        fProfile2DiffDeltaKaonProtonSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaKaonProtonSPt_%s",Chargetype.data()),";centrality;#Sum pT;#delta",16,0.,80.,16,0,16);
        fProfile2DiffDeltaKaonProtonDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaKaonProtonDEta_%s",Chargetype.data()),";centrality;#Delta Eta;#delta",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffDeltaKaonProtonQ2[iType]   = new TProfile2D(Form("fProfile2DiffDeltaKaonProtonQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffDeltaKaonProtonDPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaKaonProtonSPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaKaonProtonDEta[iType]);
        fResultsList->Add(fProfile2DiffDeltaKaonProtonQ2[iType]);
      }
      if (isCalculatePionPion) {
        fProfile2DiffDeltaPionPionDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaPionPionDPt_%s",Chargetype.data()),";centrality;#Delta pT;#delta",16,0.,80.,16,-8,8);
        fProfile2DiffDeltaPionPionSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaPionPionSPt_%s",Chargetype.data()),";centrality;#Sum pT;#delta",16,0.,80.,16,0,16);
        fProfile2DiffDeltaPionPionDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaPionPionDEta_%s",Chargetype.data()),";centrality;#Delta Eta;#delta",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffDeltaPionPionQ2[iType] = new TProfile2D(Form("fProfile2DiffDeltaPionPionQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffDeltaPionPionDPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaPionPionSPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaPionPionDEta[iType]);
        fResultsList->Add(fProfile2DiffDeltaPionPionQ2[iType]);
      }
      if (isCalculateKaonKaon) {
        fProfile2DiffDeltaKaonKaonDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaKaonKaonDPt_%s",Chargetype.data()),";centrality;#Delta pT;#delta",16,0.,80.,16,-8,8);
        fProfile2DiffDeltaKaonKaonSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaKaonKaonSPt_%s",Chargetype.data()),";centrality;#Sum pT;#delta",16,0.,80.,16,0,16);
        fProfile2DiffDeltaKaonKaonDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaKaonKaonDEta_%s",Chargetype.data()),";centrality;#Delta Eta;#delta",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffDeltaKaonKaonQ2[iType]   = new TProfile2D(Form("fProfile2DiffDeltaKaonKaonQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffDeltaKaonKaonDPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaKaonKaonSPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaKaonKaonDEta[iType]);
        fResultsList->Add(fProfile2DiffDeltaKaonKaonQ2[iType]);
      }
      if (isCalculateProtonProton) {
        fProfile2DiffDeltaProtonProtonDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaProtonProtonDPt_%s",Chargetype.data()),";centrality;#Delta pT;#delta",16,0.,80.,16,-8,8);
        fProfile2DiffDeltaProtonProtonSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaProtonProtonSPt_%s",Chargetype.data()),";centrality;#Sum pT;#delta",16,0.,80.,16,0,16);
        fProfile2DiffDeltaProtonProtonDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaProtonProtonDEta_%s",Chargetype.data()),";centrality;#Delta Eta;#delta",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffDeltaProtonProtonQ2[iType]   = new TProfile2D(Form("fProfile2DiffDeltaProtonProtonQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffDeltaProtonProtonDPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaProtonProtonSPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaProtonProtonDEta[iType]);
        fResultsList->Add(fProfile2DiffDeltaProtonProtonQ2[iType]);
      }
      if (isCalculateHadronHadron) {
        fProfile2DiffDeltaHadronHadronDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaHadronHadronDPt_%s",Chargetype.data()),";centrality;#Delta pT;#delta",16,0.,80.,16,-8,8);
        fProfile2DiffDeltaHadronHadronSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaHadronHadronSPt_%s",Chargetype.data()),";centrality;#Sum pT;#delta",16,0.,80.,16,0,16);
        fProfile2DiffDeltaHadronHadronDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaHadronHadronDEta_%s",Chargetype.data()),";centrality;#Delta Eta;#delta",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffDeltaHadronHadronQ2[iType]   = new TProfile2D(Form("fProfile2DiffDeltaHadronHadronQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffDeltaHadronHadronDPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaHadronHadronSPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaHadronHadronDEta[iType]);
        fResultsList->Add(fProfile2DiffDeltaHadronHadronQ2[iType]);
      }

      // γ
      if (isCalculatePionKaon) {
        fProfile2DiffGammaPionKaonDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaPionKaonDPt_%s",Chargetype.data()),";centrality;#Delta pT;#gamma",16,0.,80.,16,-8,8);
        fProfile2DiffGammaPionKaonSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaPionKaonSPt_%s",Chargetype.data()),";centrality;#Sum pT;#gamma",16,0.,80.,16,0,16);
        fProfile2DiffGammaPionKaonDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaPionKaonDEta_%s",Chargetype.data()),";centrality;#Delta #Eta;#gamma",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffGammaPionKaonQ2[iType] = new TProfile2D(Form("fProfile2DiffGammaPionKaonQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffGammaPionKaonDPt[iType]);
        fResultsList->Add(fProfile2DiffGammaPionKaonSPt[iType]);
        fResultsList->Add(fProfile2DiffGammaPionKaonDEta[iType]);
        fResultsList->Add(fProfile2DiffGammaPionKaonQ2[iType]);
      }
      if (isCalculatePionProton) {
        fProfile2DiffGammaPionProtonDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaPionProtonDPt_%s",Chargetype.data()),";centrality;#Delta pT;#gamma",16,0.,80.,16,-8,8);
        fProfile2DiffGammaPionProtonSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaPionProtonSPt_%s",Chargetype.data()),";centrality;#Sum pT;#gamma",16,0.,80.,16,0,16);
        fProfile2DiffGammaPionProtonDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaPionProtonDEta_%s",Chargetype.data()),";centrality;#Delta #Eta;#gamma",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffGammaPionProtonQ2[iType] = new TProfile2D(Form("fProfile2DiffGammaPionProtonQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffGammaPionProtonDPt[iType]);
        fResultsList->Add(fProfile2DiffGammaPionProtonSPt[iType]);
        fResultsList->Add(fProfile2DiffGammaPionProtonDEta[iType]);
        fResultsList->Add(fProfile2DiffGammaPionProtonQ2[iType]);
      }
      if (isCalculateKaonProton) {
        fProfile2DiffGammaKaonProtonDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaKaonProtonDPt_%s",Chargetype.data()),";centrality;#Delta pT;#gamma",16,0.,80.,16,-8,8);
        fProfile2DiffGammaKaonProtonSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaKaonProtonSPt_%s",Chargetype.data()),";centrality;#Sum pT;#gamma",16,0.,80.,16,0,16);
        fProfile2DiffGammaKaonProtonDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaKaonProtonDEta_%s",Chargetype.data()),";centrality;#Delta #Eta;#gamma",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffGammaKaonProtonQ2[iType] = new TProfile2D(Form("fProfile2DiffGammaKaonProtonQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffGammaKaonProtonDPt[iType]);
        fResultsList->Add(fProfile2DiffGammaKaonProtonSPt[iType]);
        fResultsList->Add(fProfile2DiffGammaKaonProtonDEta[iType]);
        fResultsList->Add(fProfile2DiffGammaKaonProtonQ2[iType]);
      }
      if (isCalculatePionPion) {
        fProfile2DiffGammaPionPionDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaPionPionDPt_%s",Chargetype.data()),";centrality;#Delta pT;#gamma",16,0.,80.,16,-8,8);
        fProfile2DiffGammaPionPionSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaPionPionSPt_%s",Chargetype.data()),";centrality;#Sum pT;#gamma",16,0.,80.,16,0,16);
        fProfile2DiffGammaPionPionDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaPionPionDEta_%s",Chargetype.data()),";centrality;#Delta #Eta;#gamma",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffGammaPionPionQ2[iType] = new TProfile2D(Form("fProfile2DiffGammaPionPionQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffGammaPionPionDPt[iType]);
        fResultsList->Add(fProfile2DiffGammaPionPionSPt[iType]);
        fResultsList->Add(fProfile2DiffGammaPionPionDEta[iType]);
        fResultsList->Add(fProfile2DiffGammaPionPionQ2[iType]);
      }
      if (isCalculateKaonKaon) {
        fProfile2DiffGammaKaonKaonDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaKaonKaonDPt_%s",Chargetype.data()),";centrality;#Delta pT;#gamma",16,0.,80.,16,-8,8);
        fProfile2DiffGammaKaonKaonSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaKaonKaonSPt_%s",Chargetype.data()),";centrality;#Sum pT;#gamma",16,0.,80.,16,0,16);
        fProfile2DiffGammaKaonKaonDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaKaonKaonDEta_%s",Chargetype.data()),";centrality;#Delta #Eta;#gamma",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffGammaKaonKaonQ2[iType] = new TProfile2D(Form("fProfile2DiffGammaKaonKaonQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffGammaKaonKaonDPt[iType]);
        fResultsList->Add(fProfile2DiffGammaKaonKaonSPt[iType]);
        fResultsList->Add(fProfile2DiffGammaKaonKaonDEta[iType]);
        fResultsList->Add(fProfile2DiffGammaKaonKaonQ2[iType]);
      }
      if (isCalculateProtonProton) {
        fProfile2DiffGammaProtonProtonDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaProtonProtonDPt_%s",Chargetype.data()),";centrality;#Delta pT;#gamma",16,0.,80.,16,-8,8);
        fProfile2DiffGammaProtonProtonSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaProtonProtonSPt_%s",Chargetype.data()),";centrality;#Sum pT;#gamma",16,0.,80.,16,0,16);
        fProfile2DiffGammaProtonProtonDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaProtonProtonDEta_%s",Chargetype.data()),";centrality;#Delta #Eta;#gamma",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffGammaProtonProtonQ2[iType] = new TProfile2D(Form("fProfile2DiffGammaProtonProtonQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffGammaProtonProtonDPt[iType]);
        fResultsList->Add(fProfile2DiffGammaProtonProtonSPt[iType]);
        fResultsList->Add(fProfile2DiffGammaProtonProtonDEta[iType]);
        fResultsList->Add(fProfile2DiffGammaProtonProtonQ2[iType]);
      }
      if (isCalculateHadronHadron) {
        fProfile2DiffGammaHadronHadronDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaHadronHadronDPt_%s",Chargetype.data()),";centrality;#Delta pT;#gamma",16,0.,80.,16,-8,8);
        fProfile2DiffGammaHadronHadronSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaHadronHadronSPt_%s",Chargetype.data()),";centrality;#Sum pT;#gamma",16,0.,80.,16,0,16);
        fProfile2DiffGammaHadronHadronDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaHadronHadronDEta_%s",Chargetype.data()),";centrality;#Delta #Eta;#gamma",16,0.,80.,16,-1.6,1.6);
        fProfile2DiffGammaHadronHadronQ2[iType]   = new TProfile2D(Form("fProfile2DiffGammaHadronHadronQ2_%s",Chargetype.data()),";centrality;Q2Bin;#delta",16,0.,80.,fNQ2Bins,0,(double)fNQ2Bins);
        fResultsList->Add(fProfile2DiffGammaHadronHadronDPt[iType]);
        fResultsList->Add(fProfile2DiffGammaHadronHadronSPt[iType]);
        fResultsList->Add(fProfile2DiffGammaHadronHadronDEta[iType]);
        fResultsList->Add(fProfile2DiffGammaHadronHadronQ2[iType]);
      }
    }
  }
  
  if (isCalculateDeltaPhiSumPhi) {
    for (int iCent = 0; iCent < 8; iCent++) {
      for (int iType = 0; iType < 2; iType++) {
        // C(Δη,Δφ) [cent][pair type]
        if(iType == 0) Chargetype = "SS";
        if(iType == 1) Chargetype = "OS";
        if (isCalculatePionKaon) {
          fHist2DEtaDPhiPionKaon[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiPionKaon_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaDPhiPionKaon[iCent][iType]);
        }
        if (isCalculatePionProton) {
          fHist2DEtaDPhiPionProton[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiPionProton_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaDPhiPionProton[iCent][iType]);
        }
        if (isCalculateKaonProton) {
          fHist2DEtaDPhiKaonProton[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiKaonProton_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaDPhiKaonProton[iCent][iType]);
        }
        if (isCalculatePionPion) {
          fHist2DEtaDPhiPionPion[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiPionPion_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaDPhiPionPion[iCent][iType]);
        }
        if (isCalculateKaonKaon) {
          fHist2DEtaDPhiKaonKaon[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiKaonKaon_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaDPhiKaonKaon[iCent][iType]);
        }
        if (isCalculateProtonProton) {
          fHist2DEtaDPhiProtonProton[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiProtonProton_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaDPhiProtonProton[iCent][iType]);
        }
        if (isCalculateHadronHadron) {
          fHist2DEtaDPhiHadronHadron[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiHadronHadron_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaDPhiHadronHadron[iCent][iType]);
        }

        // C(Δη,sφ) [cent][pair type] only TPC Plane
        if (isCalculatePionKaon) {
          fHist2DEtaSPhiPionKaon[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiPionKaon_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaSPhiPionKaon[iCent][iType]);
        }
        if (isCalculatePionProton) {
          fHist2DEtaSPhiPionProton[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiPionProton_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaSPhiPionProton[iCent][iType]);
        }
        if (isCalculateKaonProton) {
          fHist2DEtaSPhiKaonProton[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiKaonProton_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaSPhiKaonProton[iCent][iType]);
        }
        if (isCalculatePionPion) {
          fHist2DEtaSPhiPionPion[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiPionPion_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaSPhiPionPion[iCent][iType]);
        }
        if (isCalculateKaonKaon) {
          fHist2DEtaSPhiKaonKaon[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiKaonKaon_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaSPhiKaonKaon[iCent][iType]);
        }
        if (isCalculateProtonProton) {
          fHist2DEtaSPhiProtonProton[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiProtonProton_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaSPhiProtonProton[iCent][iType]);
        }
        if (isCalculateHadronHadron) {
          fHist2DEtaSPhiHadronHadron[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiHadronHadron_cent%i_%s",iCent,Chargetype.data()),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaSPhiHadronHadron[iCent][iType]);
        }
      }
    }
  }
  PostData(2,fResultsList);
  if (fDebug) Printf("Post fResultsList Data Success!");
}

//------------------------------------------------

void AliAnalysisTaskPIDCME::UserExec(Option_t *)
{
  if (fDebug) Printf("===============================We are in UserExec!================================");
  fEvtCount->Fill(1);
  //----------------------------
  // Handle
  //----------------------------
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else fEvtCount->Fill(2);

  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else fEvtCount->Fill(3);

  fAOD = dynamic_cast <AliAODEvent*> (InputEvent());
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  } else fEvtCount->Fill(4);

  fPIDResponse = handler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else fEvtCount->Fill(5);

  fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else fEvtCount->Fill(6);

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    if (!fMultSel) {
      AliError(Form("%s: Could not get AliMultSelection", GetName()));
    } else fEvtCount->Fill(7);
    if (!manager || !handler || !fAOD || !fPIDResponse || !fUtils || !fMultSel) return;
  }

  if (!manager || !handler || !fAOD || !fPIDResponse || !fUtils) return;
  if (fDebug) Printf("Handles done!");

  //----------------------------
  // Trigger
  //----------------------------
  UInt_t mask = handler->IsEventSelected();
  bool isTrigselected = false;
  if (fTrigger.EqualTo("kMB"))
  isTrigselected = mask & AliVEvent::kMB;
  else if (fTrigger.EqualTo("kINT7"))
  isTrigselected = mask & AliVEvent::kINT7;
  else if (fTrigger.EqualTo("kINT7+kSemiCentral"))
  isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kSemiCentral);
  else if (fTrigger.EqualTo("kINT7+kCentral+kSemiCentral"))
  isTrigselected = mask & (AliVEvent::kINT7 + AliVEvent::kCentral + AliVEvent::kSemiCentral);
  if (isTrigselected == false) return;
  fEvtCount->Fill(8);
  if (fDebug) Printf("trigger done!");

  //----------------------------
  // Run Number
  //----------------------------
  fRunNum = fAOD->GetRunNumber();
  fEvtCount->Fill(9);
  if (fRunNum != fOldRunNum) {
     // Load the run dependent calibration hist
      if (!LoadCalibHistForThisRun()) return;
      fRunNumBin = runNumList->at(fRunNum);
      fOldRunNum = fRunNum;
      if (fRunNumBin < 0) return;
  }
  fHistRunNumBin->Fill(fRunNumBin);
  fEvtCount->Fill(10);
  if (fDebug) Printf("read in done!");

  //----------------------------
  // Vertex
  //----------------------------
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  fVtx -> GetXYZ(fVertex);
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  if (fabs(fVertex[0])<1e-6 || fabs(fVertex[1])<1e-6 || fabs(fVertex[2])<1e-6) return;
  double dz = fVertex[2] - fAOD->GetPrimaryVertexSPD()->GetZ();
  fHistVz[0]->Fill(fVertex[2]);
  if (fabs(fVertex[2]) > fVzCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors()<1) return;
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
  // fEventCuts.SetCentralityEstimators("V0M","CL1");
  // if (!fEventCuts->AcceptEvent(fAOD) ) return;
  if (fPeriod.EqualTo("LHC10h")) if (fabs(dz)>0.5) return;
  if (fPeriod.EqualTo("LHC15o")) {
      double covTrc[6],covSPD[6];
      fVtx->GetCovarianceMatrix(covTrc);
      fAOD->GetPrimaryVertexSPD()->GetCovarianceMatrix(covSPD);
      double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
      double errTrc = TMath::Sqrt(covTrc[5]);
      double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
      if (fabs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return;
  }
  fHistVz[1]->Fill(fVertex[2]);
  for (int i = 0; i < 20; ++i) {
      if (fVertex[2] > -10+i*1 && fVertex[2] < -10+(i+1)*1) {fVzBin = i; break;}
  }
  if (fVzBin<0) return;
  fEvtCount->Fill(11);
  if (fDebug) Printf("vertex done!");

  //----------------------------
  // Centrality
  //----------------------------
  if (fPeriod.EqualTo("LHC10h")) {
    fCentV0M  = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
    fCentTRK  = fAOD->GetCentrality()->GetCentralityPercentile("TRK");
    fCentSPD0 = fAOD->GetCentrality()->GetCentralityPercentile("CL0");
    fCentSPD1 = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
  } else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    fCentV0M  = fMultSel->GetMultiplicityPercentile("V0M");
    fCentTRK  = fMultSel->GetMultiplicityPercentile("TRK");
    fCentSPD0 = fMultSel->GetMultiplicityPercentile("CL0");
    fCentSPD1 = fMultSel->GetMultiplicityPercentile("CL1");
  } else return;
  //we use centV0M as the default centrality
  fCent = fCentV0M;
  fHist2CentQA[0]->Fill(fCentV0M,fCentSPD1);
  fHist2CentQA[2]->Fill(fCentV0M,fCentTRK);
  fHist2CentQA[4]->Fill(fCentV0M,fCentSPD0);
  fHist2CentQA[6]->Fill(fCentSPD1,fCentSPD0);
  if (fabs(fCent-fCentSPD1)>fCentDiffCut) return;
  fHist2CentQA[1]->Fill(fCentV0M,fCentSPD1);
  fHist2CentQA[3]->Fill(fCentV0M,fCentTRK);
  fHist2CentQA[5]->Fill(fCentV0M,fCentSPD0);
  fHist2CentQA[7]->Fill(fCentSPD1,fCentSPD0);
  if (fCent < 0 || fCent >= 80) return;
  // cent bin
  fCentBin = (int)fCent/10;
  fHistCent[0]->Fill(fCent);
  fEvtCount->Fill(12);
  if (fDebug) Printf("centrality done!");

  //----------------------------
  // Pile up
  //----------------------------
  if (fPeriod.EqualTo("LHC10h")) if (!RemovalForRun1()) return;
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    // hMultCentQA[0]->Fill(fCent, fAOD->GetNumberOfTracks()); // raw Trk Multi Vs Cent(V0M)
    // if (PileUpMultiVertex(fAOD)) return;
    // if (!RejectEvtMultComp(fAOD)) return;
    // hMultCentQA[1]->Fill(fCent, fAOD->GetNumberOfTracks()); // Mult_Cent QA
    // if (!AODPileupCheck (fAOD)) return;
    if (!RejectEvtTFFit()) return; // 15o_pass2
  }
  fHistCent[1]->Fill(fCent);
  fEvtCount->Fill(13);
  if (fDebug) Printf("pile-up done!");
  //----------------------------
  // VZERO Plane
  //----------------------------
  if (isUseVZEROPlane) {
    isRightVZEROPlane = GetVZEROPlane();
    if (isRightVZEROPlane) {
      fEvtCount->Fill(14);
      if (fDebug) Printf("Get VZERO Plane done!");
    }
  }
  //----------------------------
  // ZDC Plane
  //----------------------------
  if (isUseZDCPlane) {
    if (fPeriod.EqualTo("LHC10h")) 
      isRightZDCPlane = GetZDCPlane();
    if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r"))  
      isRightZDCPlane = GetZDCPlaneLsFit();
    if (fPeriod.EqualTo("LHC15o")) 
      isRightZDCPlane = false;

    if (isRightZDCPlane) {
      fEvtCount->Fill(15);
      if (fDebug) Printf("Get ZDC Plane done!");
    }
  }
  //----------------------------
  // Loop Tracks / Fill Vectors
  //----------------------------
  //Reset vectors
  ResetVectors();
  fEvtCount->Fill(16);
  // must loop tracks becasue we need the TPC plane
  if (!LoopTracks()) return;
  fEvtCount->Fill(17);
  if (fDebug) Printf("Loop Tracks done!");
  //----------------------------
  // TPC Plane
  //----------------------------
  if (!GetTPCPlane()) return;
  fEvtCount->Fill(18);
  if (fDebug) Printf("Get TPC Plane done!");
  //----------------------------
  // Fill Resolution
  //----------------------------
  if(isUseOneSideTPCPlane) {
    fHist2Psi2TPCPosCent -> Fill(fCent, fPsi2TPCPos);
    fHist2Psi2TPCNegCent -> Fill(fCent, fPsi2TPCNeg);
    fProfileTPCPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2TPCNeg - fPsi2TPCPos)));
  } else {
    fHist2Psi2TPCCent -> Fill(fCent, fPsi2TPC);
  }
  if (isUseVZEROPlane) {
    fHist2Psi2V0CCent -> Fill(fCentSPD1, fPsi2V0C);
    fHist2Psi2V0ACent -> Fill(fCentSPD1, fPsi2V0A);
    fProfileV0MPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0C - fPsi2V0A)));
    if(isUseOneSideTPCPlane) {
      fProfileV0CTPCPosPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0C - fPsi2TPCPos)));
      fProfileV0ATPCPosPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0A - fPsi2TPCPos)));
      fProfileV0CTPCNegPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0C - fPsi2TPCNeg)));
      fProfileV0ATPCNegPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0A - fPsi2TPCNeg)));
    } else {
      fProfileV0CTPCPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0C - fPsi2TPC)));
      fProfileV0ATPCPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0A - fPsi2TPC)));
    }
  }
  if (isUseZDCPlane) {
    fHist2Psi1ZNCCent    -> Fill(fCent, fPsi1ZNC);
    fHist2Psi1ZNACent    -> Fill(fCent, fPsi1ZNA);
    fProfileZDCPsi1Correlation -> Fill(fCent, TMath::Cos(1*(fPsi1ZNC - fPsi1ZNA)));
    fProfileZDCPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi1ZNC - fPsi1ZNA)));
  }
  fEvtCount->Fill(19);
  //----------------------------
  // Pair
  //----------------------------
  if (isCalculatePionKaon || isCalculatePionProton || isCalculateKaonProton || isCalculatePionPion || isCalculateKaonKaon|| isCalculateProtonProton || isCalculateHadronHadron) {
    if (!PairTrkTrk()) return;
    fEvtCount->Fill(20);
    if (fDebug) Printf("Pair Trk & Trk done!"); 
  }
  //------------------
  // Post output data.
  //------------------
  PostData(1,fQAList);
  PostData(2,fResultsList);
  if (fDebug) Printf("analysis done!");
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::GetVZEROPlane()
{
  double multV0Ch[64] = {0};
  double V0XMean[3] = {0};
  double V0YMean[3] = {0};

  // [0]: M; [1]: C; [2]: A;
  double qxGE[3] = {0}, qyGE[3] = {0};
  double qxRC[3] = {0}, qyRC[3] = {0};
  double multRingGE[3] = {0};
  double psi2GE[3] = {0};
  double psi2RC[3] = {0};

  //Load the GE and RC histograms
  if (fPeriod.EqualTo("LHC10h") ) {
    for (int iCh = 0; iCh < 64; ++iCh) multV0Ch[iCh] = hMultV0Read->GetBinContent(iCh+1, fRunNumBin+1);
    for (int i = 1; i < 3; ++i) {   // [0]: M; [1]: C; [2]: A;
      V0XMean[i] = pV0XMeanRead[i]->GetBinContent(fRunNumBin+1, fCentBin+1, fVzBin+1);
      V0YMean[i] = pV0YMeanRead[i]->GetBinContent(fRunNumBin+1, fCentBin+1, fVzBin+1);
    }
  }
  if (fPeriod.EqualTo("LHC15o"))  
  for (int iCh = 0; iCh < 64; ++iCh) multV0Ch[iCh] = hMultV0->GetBinContent(iCh+1);
  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    int iCentSPD = (int)fCentSPD1;
    if (iCentSPD >= 90) return false;
    V0XMean[0] = -999.;
    V0YMean[0] = -999.;
    for (int i = 0; i < 2; ++i) {   // [1]: C; [2]: A;
      V0XMean[i+1] = hQx2mV0[i]->GetBinContent(iCentSPD+1);
      V0YMean[i+1] = hQy2mV0[i]->GetBinContent(iCentSPD+1);
    }
  }

  //Loop Over VZERO Channels
  //Gain Equalization
  for (int iCh = 0; iCh < 64; ++iCh) {
    double phi = TMath::Pi()/8. + TMath::Pi()/4.*(iCh%8);
    double multCh = 0.;
    // double multCh = fAOD->GetVZEROEqMultiplicity(iCh);
    if (fPeriod.EqualTo("LHC10h")) multCh= fAOD->GetVZEROEqMultiplicity(iCh);
    else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      AliAODVZERO* aodV0 = fAOD->GetVZEROData();
      multCh = aodV0->GetMultiplicity(iCh);
    }
    if (iCh<32) { // C
      double multChGEC = -1;
      if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC15o")) {
        if      (iCh <  8)              multChGEC = multCh/multV0Ch[iCh] * multV0Ch[0];
        else if (iCh >= 8  && iCh < 16) multChGEC = multCh/multV0Ch[iCh] * multV0Ch[8];
        else if (iCh >= 16 && iCh < 24) multChGEC = multCh/multV0Ch[iCh] * multV0Ch[16];
        else if (iCh >= 24 && iCh < 32) multChGEC = multCh/multV0Ch[iCh] * multV0Ch[24];
      }
      if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){
        int ibinV0 = fHCorrectV0ChWeghts->FindBin(fVertex[2],iCh);
        double V0chGE = (double)fHCorrectV0ChWeghts->GetBinContent(ibinV0);
        multChGEC = multCh*V0chGE;
      }
      if (multChGEC<0) continue;
      //for V0C GE
      qxGE[1] += multChGEC*TMath::Cos(2*phi);
      qyGE[1] += multChGEC*TMath::Sin(2*phi);
      multRingGE[1] += multChGEC;
    } else if (iCh>=32 && iCh<64) { // A
      double multChGEA = -1;
      if (fPeriod.EqualTo("LHC10h") || fPeriod.EqualTo("LHC15o")) {
        if      (iCh >= 32 && iCh < 40) multChGEA = multCh/multV0Ch[iCh] * multV0Ch[32];
        else if (iCh >= 40 && iCh < 48) multChGEA = multCh/multV0Ch[iCh] * multV0Ch[40];
        else if (iCh >= 48 && iCh < 56) multChGEA = multCh/multV0Ch[iCh] * multV0Ch[48];
        else if (iCh >= 56 && iCh < 64) multChGEA = multCh/multV0Ch[iCh] * multV0Ch[56];
      }
      if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")){
        int ibinV0 = fHCorrectV0ChWeghts->FindBin(fVertex[2],iCh);
        double V0chGE = (double)fHCorrectV0ChWeghts->GetBinContent(ibinV0);
        multChGEA = multCh*V0chGE;
      }
        if (multChGEA<0) continue;
        //for V0A GE
        qxGE[2] += multChGEA*TMath::Cos(2*phi);
        qyGE[2] += multChGEA*TMath::Sin(2*phi);
        multRingGE[2] += multChGEA;
    }
  }
  if (multRingGE[1] < 1.e-6 || multRingGE[2] < 1.e-6) return false;

  //VZERO GE Plane
  for (int i = 1; i < 3; i++) {
    psi2GE[i] = GetEventPlane(qxGE[i], qyGE[i], 2.);
    if (TMath::IsNaN(psi2GE[i])) return false;
  }

  //VZERO Recenter
  for (int i = 1; i < 3; ++i) {
    double qxMean = V0XMean[i];
    double qyMean = V0YMean[i];
    if (TMath::IsNaN(qxMean) || TMath::IsNaN(qyMean)) continue;
    if (qyMean < -900 || qxMean < -900) continue;
    // For 10 h, we've stored the qx/y of V0M, and they cannot been found in A.Dorbin's calib file for 15o period!
    qxRC[i] = qxGE[i] - qxMean;
    qyRC[i] = qyGE[i] - qyMean;
    psi2RC[i] = GetEventPlane(qxRC[i], qyRC[i], 2.);
    if (TMath::IsNaN(psi2RC[i])) return false;
  }

  // VZERO QA
  if (isQAVZERO) {
    //V0C
    fProfileV0CQxCent[0]->Fill(fCentSPD1, qxGE[1]);
    fProfileV0CQyCent[0]->Fill(fCentSPD1, qyGE[1]);
    fProfileV0CQxVtx[0] ->Fill(fVertex[2], qxGE[1]);
    fProfileV0CQyVtx[0] ->Fill(fVertex[2], qyGE[1]);
    fHist2CalibPsi2V0CCent[0]->Fill(fCentSPD1, psi2GE[1]);
    fHist3CalibQxQyCentV0C[0]->Fill(qxGE[1], qyGE[1], fCentSPD1);

    fProfileV0CQxCent[1]->Fill(fCentSPD1, qxRC[1]);
    fProfileV0CQyCent[1]->Fill(fCentSPD1, qyRC[1]);
    fProfileV0CQxVtx[1]->Fill(fVertex[2], qxRC[1]);
    fProfileV0CQyVtx[1]->Fill(fVertex[2], qyRC[1]);
    fHist2CalibPsi2V0CCent[1]->Fill(fCentSPD1, psi2RC[1]);
    fHist3CalibQxQyCentV0C[0]->Fill(qxRC[1], qyRC[1], fCentSPD1);
    //V0A
    fProfileV0AQxCent[0]->Fill(fCentSPD1, qxGE[2]);
    fProfileV0AQyCent[0]->Fill(fCentSPD1, qyGE[2]);
    fProfileV0AQxVtx[0]->Fill(fVertex[2], qxGE[2]);
    fProfileV0AQyVtx[0]->Fill(fVertex[2], qyGE[2]);
    fHist2CalibPsi2V0ACent[0]->Fill(fCentSPD1, psi2GE[2]);
    fHist3CalibQxQyCentV0A[0]->Fill(qxGE[2], qyGE[2], fCentSPD1);

    fProfileV0AQxCent[1]->Fill(fCentSPD1, qxRC[2]);
    fProfileV0AQyCent[1]->Fill(fCentSPD1, qyRC[2]);
    fProfileV0AQxVtx[1]->Fill(fVertex[2], qxRC[2]);
    fProfileV0AQyVtx[1]->Fill(fVertex[2], qyRC[2]);
    fHist2CalibPsi2V0ACent[1]->Fill(fCentSPD1, psi2RC[2]);
    fHist3CalibQxQyCentV0A[1]->Fill(qxRC[2], qyRC[2], fCentSPD1);
  }
  fPsi2V0C = psi2RC[1];
  fPsi2V0A = psi2RC[2];

  double Q2V0C = TMath::Hypot(qxRC[1], qyRC[1])/sqrt(multRingGE[1]);
  fQ2Bin = GetQ2Bin(Q2V0C);

  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::GetZDCPlane()
{
  if(!fPeriod.EqualTo("LHC10h")) return true; //now we only have ZDC Calibration Files of LHC10h
  AliAODZDC* fZDC = fAOD -> GetZDCData();
  if (!fZDC) return false;
  const double x[4] = {-1.75, 1.75, -1.75, 1.75};
  const double y[4] = {-1.75, -1.75, 1.75, 1.75};

  const double* towerEnegryZNC = fZDC->GetZNCTowerEnergy();
  const double* towerEnegryZNA = fZDC->GetZNATowerEnergy();
  float towerEnegryZNCGE[5] = {0};
  float towerEnegryZNAGE[5] = {0};


  for (int iTower = 0; iTower < 5; ++iTower) {
    fProfileZNCTowerMeanEnegry[0]->Fill(iTower+0.5, towerEnegryZNC[iTower]);
    fProfileZNATowerMeanEnegry[0]->Fill(iTower+0.5, towerEnegryZNA[iTower]);
  }
  

  //Loop Over ZDC Towers
  //Gain Equalization
  for (int iTower = 0; iTower < 5; ++iTower) {
    towerEnegryZNCGE[iTower] = towerEnegryZNC[iTower] / (fProfileForZNCGE->GetBinContent(iTower + 1)) * (fProfileForZNCGE->GetBinContent(2));//here to select the ref Tower
    towerEnegryZNAGE[iTower] = towerEnegryZNA[iTower] / (fProfileForZNAGE->GetBinContent(iTower + 1)) * (fProfileForZNAGE->GetBinContent(2));
  }

  for (int iTower = 0; iTower < 5; ++iTower) {
    fProfileZNCTowerMeanEnegry[1]->Fill(iTower+0.5, towerEnegryZNCGE[iTower]);
    fProfileZNATowerMeanEnegry[1]->Fill(iTower+0.5, towerEnegryZNAGE[iTower]);
  }


  float QxC = 0, QyC = 0, MC  = 0;
  float QxA = 0, QyA = 0, MA  = 0;
  for (int iTower = 0; iTower < 4; ++iTower) {
    QxC += towerEnegryZNCGE[iTower + 1] * x[iTower];
    QyC += towerEnegryZNCGE[iTower + 1] * y[iTower];
    MC  += towerEnegryZNCGE[iTower + 1];
    QxA += towerEnegryZNAGE[iTower + 1] * x[iTower];
    QyA += towerEnegryZNAGE[iTower + 1] * y[iTower];
    MA  += towerEnegryZNAGE[iTower + 1];
  }
  if (fabs(MC) < 1.e-6 || fabs(MA) < 1.e-6) return false;
  QxC /= MC;
  QyC /= MC;
  QxA /= MA;
  QyA /= MA;

  double psiCGE = GetEventPlane(QxC,QyC,1);
  double psiAGE = GetEventPlane(-QxA,QyA,1);
  if (TMath::IsNaN(psiCGE) || TMath::IsNaN(psiAGE)) return false;

  if (isQAZDC) {
    //ZNC
    fProfileZNCQxCent[0] -> Fill(fCent, QxC);
    fProfileZNCQyCent[0] -> Fill(fCent, QyC);
    fHist2CalibPsi1ZNCCent[0] -> Fill(fCent, psiCGE);
    //ZNA
    fProfileZNAQxCent[0] -> Fill(fCent,-QxA);
    fProfileZNAQyCent[0] -> Fill(fCent, QyA);
    fHist2CalibPsi1ZNACent[0] -> Fill(fCent, psiAGE);
    //ZNC-ZNA
    fProfileZDCQxAQxCCent[0] -> Fill(fCent,-QxA*QxC);
    fProfileZDCQxAQyCCent[0] -> Fill(fCent,-QxA*QyC);
    fProfileZDCQyAQxCCent[0] -> Fill(fCent, QyA*QxC);
    fProfileZDCQyAQyCCent[0] -> Fill(fCent, QyA*QyC);
  }

  //ZDC Recenter
  double fillPosition[4] = {fCent,-999.,-999.,-999.};
  int vtxBin[3] = {-999,-999,-999};
  for (int i = 0; i < 3; i++) {
    if (fVertex[i] <= vtxQuant1[i])
    vtxBin[i]=1;
    else if (fVertex[i] > vtxQuant1[i] && fVertex[i] <= vtxQuant2[i])
    vtxBin[i]=2;
    else if (fVertex[i] > vtxQuant2[i])
    vtxBin[i]=3;
  }
  fillPosition[1] = vtxBin[0] - 0.5;
  fillPosition[2] = vtxBin[1] - 0.5;
  fillPosition[3] = vtxBin[2] - 0.5;

  float QxCMean = fHn4DForZNCQxRC -> GetBinContent(fHn4DForZNCQxRC->GetBin(fillPosition));
  float QyCMean = fHn4DForZNCQyRC -> GetBinContent(fHn4DForZNCQyRC->GetBin(fillPosition));
  float QxAMean = fHn4DForZNAQxRC -> GetBinContent(fHn4DForZNAQxRC->GetBin(fillPosition));
  float QyAMean = fHn4DForZNAQyRC -> GetBinContent(fHn4DForZNAQyRC->GetBin(fillPosition));
  int entriesC = fHn4DForZNCCountsRC -> GetBinContent(fHn4DForZNCCountsRC->GetBin(fillPosition));
  int entriesA = fHn4DForZNACountsRC -> GetBinContent(fHn4DForZNACountsRC->GetBin(fillPosition));

  if (fabs(entriesC) < 2 || fabs(entriesA) < 2) return false;
  QxCMean /= entriesC;
  QyCMean /= entriesC;
  QxAMean /= entriesA;
  QyAMean /= entriesA;

  QxC -= QxCMean;
  QyC -= QyCMean;
  QxA -= QxAMean;
  QyA -= QyAMean;

  double psiCRC = GetEventPlane(QxC,QyC,1);
  double psiARC = GetEventPlane(-QxA,QyA,1);
  if (TMath::IsNaN(psiCRC) || TMath::IsNaN(psiARC)) return false;

  if (isQAZDC) {
    //ZNC
    fProfileZNCQxCent[1] -> Fill(fCent, QxC);
    fProfileZNCQyCent[1] -> Fill(fCent, QyC);
    fHist2CalibPsi1ZNCCent[1] -> Fill(fCent, psiCRC);
    //ZNA
    fProfileZNAQxCent[1] -> Fill(fCent,-QxA);
    fProfileZNAQyCent[1] -> Fill(fCent, QyA);
    fHist2CalibPsi1ZNACent[1] -> Fill(fCent, psiARC);
    //ZNC-ZNA
    fProfileZDCQxAQxCCent[1] -> Fill(fCent,-QxA*QxC);
    fProfileZDCQxAQyCCent[1] -> Fill(fCent,-QxA*QyC);
    fProfileZDCQyAQxCCent[1] -> Fill(fCent, QyA*QxC);
    fProfileZDCQyAQyCCent[1] -> Fill(fCent, QyA*QyC);
  }

  //ZDC Shift
  double psiCSF = psiCRC;
  double psiASF = psiARC;
  for (int i = 1; i <= 20; i++) {
    double shiftCosC = fProfile2DForCosC->GetBinContent(fProfile2DForCosC->GetXaxis()->FindBin(fCent),i);
    double shiftSinC = fProfile2DForSinC->GetBinContent(fProfile2DForSinC->GetXaxis()->FindBin(fCent),i);
    double shiftCosA = fProfile2DForCosA->GetBinContent(fProfile2DForCosA->GetXaxis()->FindBin(fCent),i);
    double shiftSinA = fProfile2DForSinA->GetBinContent(fProfile2DForSinA->GetXaxis()->FindBin(fCent),i);
    psiCSF += (2./i) * (-shiftSinC * TMath::Cos(i*psiCSF) + shiftCosC * TMath::Sin(i*psiCSF));
    psiASF += (2./i) * (-shiftSinA * TMath::Cos(i*psiASF) + shiftCosA * TMath::Sin(i*psiASF));
  }

  if (TMath::IsNaN(psiCSF) || TMath::IsNaN(psiASF)) return false;
  if (isQAZDC) {
    fHist2CalibPsi1ZNCCent[2] -> Fill(fCent, psiCSF);
    fHist2CalibPsi1ZNACent[2] -> Fill(fCent, psiASF);
  }

  fPsi1ZNC = psiCSF;
  fPsi1ZNA = psiASF;
  return true;
}

//---------------------------------------------------
bool AliAnalysisTaskPIDCME::GetZDCPlaneLsFit()
{
  if(!(fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r"))) return true;
  AliAODZDC* fZDC = fAOD -> GetZDCData();
  if (!fZDC) return false;
  const double* fZNATowerRawAOD = fZDC->GetZNATowerEnergy();
  const double* fZNCTowerRawAOD = fZDC->GetZNCTowerEnergy();
  for (int iTower = 0; iTower < 5; iTower++){
    if(fZNATowerRawAOD[iTower] < 1.e-6) return false;
    if(fZNCTowerRawAOD[iTower] < 1.e-6) return false;
  }

  double towZNCraw1GainEq = 0, towZNCraw2GainEq = 0, towZNCraw3GainEq = 0, towZNCraw4GainEq = 0;
  towZNCraw1GainEq = fZNCTowerRawAOD[1] * fHZDCCparameters->GetBinContent(1);
  towZNCraw2GainEq = fZNCTowerRawAOD[2] * fHZDCCparameters->GetBinContent(2);
  towZNCraw3GainEq = fZNCTowerRawAOD[3] * fHZDCCparameters->GetBinContent(3);
  towZNCraw4GainEq = fZNCTowerRawAOD[4] * fHZDCCparameters->GetBinContent(4);

  double towZNAraw1GainEq = 0, towZNAraw2GainEq = 0, towZNAraw3GainEq = 0, towZNAraw4GainEq = 0;
  towZNAraw1GainEq = fZNATowerRawAOD[1] * fHZDCAparameters->GetBinContent(1);
  towZNAraw2GainEq = fZNATowerRawAOD[2] * fHZDCAparameters->GetBinContent(2);
  towZNAraw3GainEq = fZNATowerRawAOD[3] * fHZDCAparameters->GetBinContent(3);
  towZNAraw4GainEq = fZNATowerRawAOD[4] * fHZDCAparameters->GetBinContent(4);

  const double xZDCC[4] = {-1,  1, -1,  1}; // directional vector
  const double yZDCC[4] = {-1, -1,  1,  1};
  const double xZDCA[4] = { 1, -1,  1, -1};
  const double yZDCA[4] = {-1, -1,  1,  1};

  double towZNC[5] = {fZNCTowerRawAOD[0], towZNCraw1GainEq, towZNCraw2GainEq, towZNCraw3GainEq, towZNCraw4GainEq};
  double towZNA[5] = {fZNATowerRawAOD[0], towZNAraw1GainEq, towZNAraw2GainEq, towZNAraw3GainEq, towZNAraw4GainEq};
  
  //double EZNC = 0;
  double wZNC = 0, denZNC = 0, numXZNC = 0, numYZNC = 0;
  //double EZNA = 0;
  double wZNA = 0, denZNA = 0, numXZNA = 0, numYZNA = 0; 

  for(int i = 0; i < 4; i++) {
    // ZNC 
    // get energy
    // EZNC = towZNC[i+1];
    // build ZDCC centroid
    wZNC = TMath::Max(0., 4.0 + TMath::Log(towZNC[i+1]/fZNCTowerRawAOD[0]));
    numXZNC += xZDCC[i]*wZNC;
    numYZNC += yZDCC[i]*wZNC;
    denZNC += wZNC;
        
    // ZNA part
    // get energy
    // EZNA = towZNA[i+1];
    // build ZDCA centroid
    wZNA = TMath::Max(0., 4.0 + TMath::Log(towZNA[i+1]/fZNATowerRawAOD[0]));
    numXZNA += xZDCA[i]*wZNA;
    numYZNA += yZDCA[i]*wZNA;
    denZNA += wZNA;
  }
  if (fabs(denZNC) < 1.e-6) return false;
  if (fabs(denZNA) < 1.e-6) return false;

  double ZDCCxPosFromLogWeight = numXZNC/denZNC;
  double ZDCCyPosFromLogWeight = numYZNC/denZNC;
  double ZDCAxPosFromLogWeight = numXZNA/denZNA;
  double ZDCAyPosFromLogWeight = numYZNA/denZNA;

  //QA
  fProfileZNCQxCent[0] -> Fill(fCent,ZDCCxPosFromLogWeight);
  fProfileZNCQyCent[0] -> Fill(fCent,ZDCCyPosFromLogWeight);
  fProfileZNAQxCent[0] -> Fill(fCent,ZDCAxPosFromLogWeight);
  fProfileZNAQyCent[0] -> Fill(fCent,ZDCAyPosFromLogWeight);

  fProfileZDCQxAQxCCent[0] -> Fill(fCent,ZDCAxPosFromLogWeight*ZDCCxPosFromLogWeight);
  fProfileZDCQxAQyCCent[0] -> Fill(fCent,ZDCAxPosFromLogWeight*ZDCCyPosFromLogWeight);
  fProfileZDCQyAQxCCent[0] -> Fill(fCent,ZDCAyPosFromLogWeight*ZDCCxPosFromLogWeight);
  fProfileZDCQyAQyCCent[0] -> Fill(fCent,ZDCAyPosFromLogWeight*ZDCCyPosFromLogWeight);

  double psiZNCGE = GetEventPlane(ZDCCxPosFromLogWeight,ZDCCyPosFromLogWeight,1);
  double psiZNAGE = GetEventPlane(ZDCAxPosFromLogWeight,ZDCAyPosFromLogWeight,1);
  if (TMath::IsNaN(psiZNCGE) || TMath::IsNaN(psiZNAGE)) return false;

  fHist2CalibPsi1ZNCCent[0] ->Fill(fCent,psiZNCGE);
  fHist2CalibPsi1ZNACent[0] ->Fill(fCent,psiZNAGE);


  //Recenter
  double ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(5) *fCent + fHZDCCparameters->GetBinContent(6) *fVertex[0] + fHZDCCparameters->GetBinContent(7) *fVertex[1] + fHZDCCparameters->GetBinContent(8) *fVertex[2] + fHZDCCparameters->GetBinContent(9);
  double ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(10)*fCent + fHZDCCparameters->GetBinContent(11)*fVertex[0] + fHZDCCparameters->GetBinContent(12)*fVertex[1] + fHZDCCparameters->GetBinContent(13)*fVertex[2] + fHZDCCparameters->GetBinContent(14);
  double ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(5) *fCent + fHZDCAparameters->GetBinContent(6) *fVertex[0] + fHZDCAparameters->GetBinContent(7) *fVertex[1] + fHZDCAparameters->GetBinContent(8) *fVertex[2] + fHZDCAparameters->GetBinContent(9);
  double ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(10)*fCent + fHZDCAparameters->GetBinContent(11)*fVertex[0] + fHZDCAparameters->GetBinContent(12)*fVertex[1] + fHZDCAparameters->GetBinContent(13)*fVertex[2] + fHZDCAparameters->GetBinContent(14);
  
  double QxZNC = ZDCCxPosFromLogWeight - ZDCCAvgxPosFromVtxFit;
  double QyZNC = ZDCCyPosFromLogWeight - ZDCCAvgyPosFromVtxFit;
  
  double QxZNA = ZDCAxPosFromLogWeight - ZDCAAvgxPosFromVtxFit;
  double QyZNA = ZDCAyPosFromLogWeight - ZDCAAvgyPosFromVtxFit;

  //QA
  fProfileZNCQxCent[1] -> Fill(fCent, QxZNC);
  fProfileZNCQyCent[1] -> Fill(fCent, QyZNC);
  fProfileZNAQxCent[1] -> Fill(fCent, QxZNA);
  fProfileZNAQyCent[1] -> Fill(fCent, QyZNA);

  fProfileZDCQxAQxCCent[1] -> Fill(fCent,QxZNA*QxZNC);
  fProfileZDCQxAQyCCent[1] -> Fill(fCent,QxZNA*QyZNC);
  fProfileZDCQyAQxCCent[1] -> Fill(fCent,QyZNA*QxZNC);
  fProfileZDCQyAQyCCent[1] -> Fill(fCent,QyZNA*QyZNC);

  double psiZNCRC = GetEventPlane(QxZNC,QyZNC,1);
  double psiZNARC = GetEventPlane(QxZNA,QyZNA,1);
  if (TMath::IsNaN(psiZNCRC) || TMath::IsNaN(psiZNARC)) return false;

  fHist2CalibPsi1ZNCCent[1] ->Fill(fCent,psiZNCRC);
  fHist2CalibPsi1ZNACent[1] ->Fill(fCent,psiZNARC);

  fPsi1ZNC = psiZNCRC;
  fPsi1ZNA = psiZNARC;

  return true;
}
//---------------------------------------------------

bool AliAnalysisTaskPIDCME::LoopTracks()
{
  int nTrks = fAOD->GetNumberOfTracks();
  if (nTrks < 4) return false;
  for (int iTrk = 0; iTrk < nTrks; ++iTrk) {
    AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(iTrk);
    if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }
    if (!track->TestFilterBit(fFilterBit)) continue;
    if (!AcceptAODTrack(track)) continue;
    //------------------
    // NUE & NUA
    //------------------
    double phi  = track->Phi();
    double  pt  = track->Pt();
    double eta  = track->Eta();
    int charge  = track->Charge();
    int     id  = track->GetID();
    int  nhits  = track->GetTPCNcls();
    double dedx = track->GetTPCsignal();
    
    fHistPt->Fill(pt);
    fHistEta->Fill(eta);
    fHistNhits->Fill(nhits);
    fHist2PDedx->Fill(track->P()*charge, dedx);
    fHistPhi[0]->Fill(phi);

    double weight = 1.;
    if (isDoNUE) {
      double wEff = GetNUECor(charge, pt);
      if (wEff < 0) continue;
      else weight *= wEff;
    }
    if (isDoNUA) {
      double wAcc = GetNUACor(charge, phi, eta, fVertex[2]);
      if (wAcc < 0) continue;
      else weight *= wAcc;
      fHistPhi[1]->Fill(phi, wAcc);
    }

    //DCA Cut
    double dcaxy = -999, dcaz = -999;
    if(!GetDCA(dcaxy,dcaz,track)) continue;

    // if FB = 1, we need to cut dca for the plane
    if (fFilterBit == 1) {
      if (fabs(dcaz) > fDcaCutZ) continue;
      if (fabs(dcaxy) > fDcaCutXY) continue;
    }
    // if FB = 96 or 768, we don't need cut dca for the plane

    fHistDcaXY->Fill(fabs(dcaxy));
    fHistDcaZ->Fill(fabs(dcaz));

    if (pt > fPlanePtMin && pt < fPlanePtMax) {
      //Do we need to set pT as weight for Better resolution?
      if (isUseOneSideTPCPlane) {
        if (eta >= fEtaGapPos) {
          fSumQ2xTPCPos  += weight * TMath::Cos(2 * phi);
          fSumQ2yTPCPos  += weight * TMath::Sin(2 * phi);
          fWgtMultTPCPos += weight;
          std::vector<double> vec_phi_weight;
          vec_phi_weight.emplace_back(phi);
          vec_phi_weight.emplace_back(weight);
          mapTPCPosTrksIDPhiWgt[id] = vec_phi_weight;
        } else if (eta <= fEtaGapNeg) {
          fSumQ2xTPCNeg  += weight * TMath::Cos(2 * phi);
          fSumQ2yTPCNeg  += weight * TMath::Sin(2 * phi);
          fWgtMultTPCNeg += weight;
          std::vector<double> vec_phi_weight;
          vec_phi_weight.emplace_back(phi);
          vec_phi_weight.emplace_back(weight);
          mapTPCNegTrksIDPhiWgt[id] = vec_phi_weight;
        } else continue;
      } else {
        fSumQ2xTPC  += weight * TMath::Cos(2 * phi);
        fSumQ2yTPC  += weight * TMath::Sin(2 * phi);
        fWgtMultTPC += weight;
        std::vector<double> vec_phi_weight;
        vec_phi_weight.emplace_back(phi);
        vec_phi_weight.emplace_back(weight);
        mapTPCTrksIDPhiWgt[id] = vec_phi_weight;
      }
    }
    
    // but we need to set the dca cut for 768 when we start to choose the paiticle for pair(just for NarrowDCACut Set)
    if (fFilterBit == 768 && isNarrowDcaCuts768) { 
      if (fabs(dcaz) > 2.0) continue;
      if (fabs(dcaxy) > 7.0 * (0.0026 + 0.005/TMath::Power(pt, 1.01))) continue;
    }

    bool isItPiontrk = CheckPIDofParticle(track,1); // 1=pion
    bool isItKaontrk = CheckPIDofParticle(track,2); // 2=kaon
    bool isItProttrk = CheckPIDofParticle(track,3); // 3=proton
    // Proton Maybe need a customized DCA cut
    if(isProtonCustomizedDCACut) isItProttrk = isItProttrk && (fabs(dcaz) < 1. && fabs(dcaxy) < (0.0105 + 0.035/TMath::Power(pt,1.1))); 
    
    if(isOpenPIDSingletrk)  {
      int PIDPoints = 0;
      if(isItPiontrk && isItKaontrk) {PIDPoints = 1;}
      if(isItPiontrk && isItProttrk) {PIDPoints = 2;}
      if(isItKaontrk && isItProttrk) {PIDPoints = 3;}
      if(isItPiontrk && isItKaontrk && isItProttrk) {PIDPoints = 4;}
      int BestPIDPoints = 0;
      BestPIDPoints = DetermineTheBestTPCPID(track, PIDPoints);
      if(BestPIDPoints == 1) {isItPiontrk = true; isItKaontrk = false; isItProttrk = false;}
      if(BestPIDPoints == 2) {isItPiontrk = false; isItKaontrk = true; isItProttrk = false;}
      if(BestPIDPoints == 3) {isItPiontrk = false; isItKaontrk = false; isItProttrk = true;}
    }
    int code = 0;
    double pid_weight = 1.;
    if(charge > 0) {
      code = 999;
      isItProttrk = isItProttrk && (pt < fProtonPtMax && pt > fProtonPtMin);
      isItKaontrk = isItKaontrk && (pt < fKaonPtMax && pt > fKaonPtMin);
      isItPiontrk = isItPiontrk && (pt < fPionPtMax && pt > fPionPtMin);
      if(isItProttrk) {
        code =  2212;
        fHistProtonPt->Fill(pt);
        fHistProtonEta->Fill(eta);
        fHistProtonPhi->Fill(phi);
        fHistProtonDcaXY->Fill(fabs(dcaxy));
        fHistProtonDcaZ ->Fill(fabs(dcaz));
        float nSigTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        float nSigTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
        float nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
        fHist2ProtonSigTPC -> Fill(pt,nSigTPC);
        fHist2ProtonSigTOF -> Fill(pt,nSigTOF);
        fHist2ProtonSigRMS -> Fill(pt,nSigRMS);
        fHist3ProtonSigTPCTOFPt -> Fill(nSigTPC,nSigTOF,pt);
        if(isDoNUE) {
          double nue_pid_weight = GetPIDNUECor(2212,pt);
          if (nue_pid_weight > 0) pid_weight *= nue_pid_weight;
        }
        vecParticle.emplace_back(std::array<double,7>{pt,eta,phi,(double)id,(double)code,weight,pid_weight});
      }
      if(isItKaontrk) {
        code =  321;
        fHistKaonPt->Fill(pt);
        fHistKaonEta->Fill(eta);
        fHistKaonPhi->Fill(phi);
        fHistKaonDcaXY->Fill(fabs(dcaxy));
        fHistKaonDcaZ ->Fill(fabs(dcaz));
        float nSigTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        float nSigTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
        float nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
        fHist2KaonSigTPC -> Fill(pt,nSigTPC);
        fHist2KaonSigTOF -> Fill(pt,nSigTOF);
        fHist2KaonSigRMS -> Fill(pt,nSigRMS);
        fHist3KaonSigTPCTOFPt -> Fill(nSigTPC,nSigTOF,pt);
        if(isDoNUE) {
          double nue_pid_weight = GetPIDNUECor(321,pt);
          if (nue_pid_weight > 0) pid_weight *= nue_pid_weight;
        }
        vecParticle.emplace_back(std::array<double,7>{pt,eta,phi,(double)id,(double)code,weight,pid_weight});
      }
      if(isItPiontrk) {
        code =  211;
        fHistPionPt->Fill(pt);
        fHistPionEta->Fill(eta);
        fHistPionPhi->Fill(phi);
        fHistPionDcaXY->Fill(fabs(dcaxy));
        fHistPionDcaZ->Fill(fabs(dcaz));
        float nSigTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        float nSigTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
        float nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
        fHist2PionSigTPC -> Fill(pt,nSigTPC);
        fHist2PionSigTOF -> Fill(pt,nSigTOF);
        fHist2PionSigRMS -> Fill(pt,nSigRMS);
        fHist3PionSigTPCTOFPt -> Fill(nSigTPC,nSigTOF,pt);
        if(isDoNUE) {
          double nue_pid_weight = GetPIDNUECor(211,pt);
          if (nue_pid_weight > 0) pid_weight *= nue_pid_weight;
        }
        vecParticle.emplace_back(std::array<double,7>{pt,eta,phi,(double)id,(double)code,weight,pid_weight});
      }
      if(!(isItProttrk || isItKaontrk || isItPiontrk)) vecParticle.emplace_back(std::array<double,7>{pt,eta,phi,(double)id,(double)code,weight,pid_weight});
    } else {
      code =-999;
      isItProttrk = isItProttrk && (pt < fAntiProtonPtMax && pt > fAntiProtonPtMin);
      isItKaontrk = isItKaontrk && (pt < fAntiKaonPtMax && pt > fAntiKaonPtMin);
      isItPiontrk = isItPiontrk && (pt < fAntiPionPtMax && pt > fAntiPionPtMin);
      if(isItProttrk) {
        code = -2212;
        fHistAntiProtonPt->Fill(pt);
        fHistAntiProtonEta->Fill(eta);
        fHistAntiProtonPhi->Fill(phi);
        fHistAntiProtonDcaXY->Fill(fabs(dcaxy));
        fHistAntiProtonDcaZ ->Fill(fabs(dcaz));
        float nSigTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        float nSigTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
        float nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
        fHist2AntiProtonSigTPC -> Fill(pt, nSigTPC);
        fHist2AntiProtonSigTOF -> Fill(pt, nSigTOF);
        fHist2AntiProtonSigRMS -> Fill(pt, nSigRMS);
        fHist3AntiProtonSigTPCTOFPt -> Fill(nSigTPC,nSigTOF,pt);
        if(isDoNUE) {
          double nue_pid_weight = GetPIDNUECor(-2212,pt);
          if (nue_pid_weight > 0) pid_weight *= nue_pid_weight;
        }
        vecParticle.emplace_back(std::array<double,7>{pt,eta,phi,(double)id,(double)code,weight,pid_weight});
      }
      if(isItKaontrk) {
        code = -321;
        fHistAntiKaonPt->Fill(pt);
        fHistAntiKaonEta->Fill(eta);
        fHistAntiKaonPhi->Fill(phi);
        fHistAntiKaonDcaXY->Fill(fabs(dcaxy));
        fHistAntiKaonDcaZ->Fill(fabs(dcaz));
        float nSigTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        float nSigTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
        float nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
        fHist2AntiKaonSigTPC -> Fill(pt, nSigTPC);
        fHist2AntiKaonSigTOF -> Fill(pt, nSigTOF);
        fHist2AntiKaonSigRMS -> Fill(pt, nSigRMS);
        fHist3AntiKaonSigTPCTOFPt -> Fill(nSigTPC,nSigTOF,pt);
        if(isDoNUE) {
          double nue_pid_weight = GetPIDNUECor(-321,pt);
          if (nue_pid_weight > 0) pid_weight *= nue_pid_weight;
        }
        vecParticle.emplace_back(std::array<double,7>{pt,eta,phi,(double)id,(double)code,weight,pid_weight});
      }
      if(isItPiontrk) {
        code = -211;
        fHistAntiPionPt->Fill(pt);
        fHistAntiPionEta->Fill(eta);
        fHistAntiPionPhi->Fill(phi);
        fHistAntiPionDcaXY->Fill(fabs(dcaxy));
        fHistAntiPionDcaZ->Fill(fabs(dcaz));
        float nSigTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        float nSigTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
        float nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
        fHist2AntiPionSigTPC -> Fill(pt, nSigTPC);
        fHist2AntiPionSigTOF -> Fill(pt, nSigTOF);
        fHist2AntiPionSigRMS -> Fill(pt, nSigRMS);
        fHist3AntiPionSigTPCTOFPt -> Fill(nSigTPC,nSigTOF,pt);
        if(isDoNUE) {
          double nue_pid_weight = GetPIDNUECor(-211,pt);
          if (nue_pid_weight > 0) pid_weight *= nue_pid_weight;
        }
        vecParticle.emplace_back(std::array<double,7>{pt,eta,phi,(double)id,(double)code,weight,pid_weight});
      }
      if(!(isItProttrk || isItKaontrk || isItPiontrk)) vecParticle.emplace_back(std::array<double,7>{pt,eta,phi,(double)id,(double)code,weight,pid_weight});
    }
  }

  if(isUseOneSideTPCPlane) {
    if(fabs(fSumQ2xTPCNeg)<1.e-6 || fabs(fSumQ2yTPCNeg)<1.e-6 || fWgtMultTPCNeg < 1.e-6) return false;
    if(fabs(fSumQ2xTPCPos)<1.e-6 || fabs(fSumQ2yTPCPos)<1.e-6 || fWgtMultTPCPos < 1.e-6) return false;
  } else {
    if(fabs(fSumQ2xTPC)<1.e-6 || fabs(fSumQ2yTPC)<1.e-6 || fWgtMultTPC < 1.e-6) return false;
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::GetTPCPlane()
{
  double psi2TPCPos = -999;
  double psi2TPCNeg = -999;
  double psi2TPC    = -999;
  if(isUseOneSideTPCPlane) {
    psi2TPCPos = GetEventPlane(fSumQ2xTPCPos,fSumQ2yTPCPos,2);
    psi2TPCNeg = GetEventPlane(fSumQ2xTPCNeg,fSumQ2yTPCNeg,2);
  } else {
    psi2TPC = GetEventPlane(fSumQ2xTPC,fSumQ2yTPC,2);
  }
  fPsi2TPCPos = psi2TPCPos;
  fPsi2TPCNeg = psi2TPCNeg;
  fPsi2TPC    = psi2TPC;
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::PairTrkTrk()
{
  for (auto paticle_a : vecParticle) {
    double pt_a     = paticle_a[0];
    double eta_a    = paticle_a[1];
    double phi_a    = paticle_a[2];
    int    id_a     = (int)paticle_a[3];
    int    code_a   = (int)paticle_a[4];
    double weight_a = paticle_a[5];
    double pid_weight_a = paticle_a[6];

    // Get TPC Plane
    std::array<double, 2> arrOneSidePsi2TPCNoAuto = {nan(""), nan("")};
    if (isUseOneSideTPCPlane) {
      arrOneSidePsi2TPCNoAuto = GetOneSideTPCPlaneNoAutoCorr({id_a});
    }
    if (isCalculatePIDFlow) {
      double v2_tmp_plane[5] = {0.,0.,0.,0.,0.};
      if(isUseOneSideTPCPlane) {
        if(eta_a > 0.) v2_tmp_plane[0] = TMath::Cos(2 * (phi_a - fPsi2TPCNeg));
        else         v2_tmp_plane[0] = TMath::Cos(2 * (phi_a - fPsi2TPCPos));
      } else {
        double psi2TPC_forFlow = GetTPCPlaneNoAutoCorr({id_a});
        v2_tmp_plane[0] = TMath::Cos(2 * (phi_a - psi2TPC_forFlow));
      }
      v2_tmp_plane[1] = TMath::Cos(2 * (phi_a - fPsi2V0C));
      v2_tmp_plane[2] = TMath::Cos(2 * (phi_a - fPsi2V0A));
      v2_tmp_plane[3] = TMath::Cos(2 * (phi_a - fPsi1ZNC));
      v2_tmp_plane[4] = TMath::Cos(2 * (phi_a - fPsi1ZNA));
      fProfile2RawFlowHadronQ2 -> Fill(fCent, fQ2Bin, v2_tmp_plane[0]);
      if(code_a > 0) {
        fProfile2RawFlowPtCentHadron[0][0] -> Fill(fCent, pt_a, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentHadron[1][0] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentHadron[2][0] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentHadron[3][0] -> Fill(fCent, pt_a, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentHadron[4][0] -> Fill(fCent, pt_a, v2_tmp_plane[4]);
        }
      }
      if(code_a < 0) {
        fProfile2RawFlowPtCentHadron[0][1] -> Fill(fCent, pt_a, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentHadron[1][1] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentHadron[2][1] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentHadron[3][1] -> Fill(fCent, pt_a, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentHadron[4][1] -> Fill(fCent, pt_a, v2_tmp_plane[4]);
        }
      }
      if(code_a == 211) {
        fProfile2RawFlowPtCentPion[0][0] -> Fill(fCent, pt_a, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentPion[1][0] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentPion[2][0] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentPion[3][0] -> Fill(fCent, pt_a, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentPion[4][0] -> Fill(fCent, pt_a, v2_tmp_plane[4]);
        }
        fProfile2RawFlowPionQ2 -> Fill(fCent, fQ2Bin, v2_tmp_plane[0]);
      }
      if(code_a == -211) {
        fProfile2RawFlowPtCentPion[0][1] -> Fill(fCent, pt_a, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentPion[1][1] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentPion[2][1] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentPion[3][1] -> Fill(fCent, pt_a, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentPion[4][1] -> Fill(fCent, pt_a, v2_tmp_plane[4]);
        }
        fProfile2RawFlowPionQ2 -> Fill(fCent, fQ2Bin, v2_tmp_plane[0]);
      }
      if(code_a == 321) {
        fProfile2RawFlowPtCentKaon[0][0] -> Fill(fCent, pt_a, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentKaon[1][0] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentKaon[2][0] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentKaon[3][0] -> Fill(fCent, pt_a, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentKaon[4][0] -> Fill(fCent, pt_a, v2_tmp_plane[4]);
        }
        fProfile2RawFlowKaonQ2 -> Fill(fCent, fQ2Bin, v2_tmp_plane[0]);
      }
      if(code_a == -321) {
        fProfile2RawFlowPtCentKaon[0][1] -> Fill(fCent, pt_a, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentKaon[1][1] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentKaon[2][1] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentKaon[3][1] -> Fill(fCent, pt_a, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentKaon[4][1] -> Fill(fCent, pt_a, v2_tmp_plane[4]);
        }
        fProfile2RawFlowKaonQ2 -> Fill(fCent, fQ2Bin, v2_tmp_plane[0]);
      }
      if(code_a == 2212) {
        fProfile2RawFlowPtCentProton[0][0] -> Fill(fCent, pt_a, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentProton[1][0] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentProton[2][0] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentProton[3][0] -> Fill(fCent, pt_a, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentProton[4][0] -> Fill(fCent, pt_a, v2_tmp_plane[4]);
        }
        fProfile2RawFlowProtonQ2 -> Fill(fCent, fQ2Bin, v2_tmp_plane[0]);
      }
      if(code_a == -2212) {
        fProfile2RawFlowPtCentProton[0][1] -> Fill(fCent, pt_a, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentProton[1][1] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentProton[2][1] -> Fill(fCentSPD1, pt_a, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentProton[3][1] -> Fill(fCent, pt_a, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentProton[4][1] -> Fill(fCent, pt_a, v2_tmp_plane[4]);
        }
        fProfile2RawFlowProtonQ2 -> Fill(fCent, fQ2Bin, v2_tmp_plane[0]);
      }    
    } //PIDFlow over
    for (auto paticle_b : vecParticle) {
      double pt_b     = paticle_b[0];
      double eta_b    = paticle_b[1];
      double phi_b    = paticle_b[2];
      int    id_b     = (int)paticle_b[3];
      int    code_b   = (int)paticle_b[4];
      double weight_b = paticle_b[5];
      double pid_weight_b = paticle_b[6];
      if (id_a == id_b) continue;

      // Get TPC Plane For This Pair, Remove TPC Plane AutoCorrelation
      double psi2TPC_forThisPair = nan("");
      if (isUseOneSideTPCPlane) {
        if (eta_b > 0) psi2TPC_forThisPair = arrOneSidePsi2TPCNoAuto[1];
        else           psi2TPC_forThisPair = arrOneSidePsi2TPCNoAuto[0];
      } else {
        psi2TPC_forThisPair = GetTPCPlaneNoAutoCorr({id_a, id_b});
      }
      if (std::isnan(psi2TPC_forThisPair)) continue;

      double delta = TMath::Cos(phi_a - phi_b);
      double gamma[5] = {0.,0.,0.,0.,0.,};
      gamma[0] = TMath::Cos(phi_a + phi_b - 2 * psi2TPC_forThisPair);
      gamma[1] = TMath::Cos(phi_a + phi_b - 2 * fPsi2V0C);
      gamma[2] = TMath::Cos(phi_a + phi_b - 2 * fPsi2V0A);
      gamma[3] = TMath::Cos(phi_a + phi_b - 2 * fPsi1ZNC);
      gamma[4] = TMath::Cos(phi_a + phi_b - 2 * fPsi1ZNA);

      int nBits = 6;
      if (isCalculatePionKaon) {
        TBits bitsPionKaonPair(nBits);
        bitsPionKaonPair.SetBitNumber(0, (code_a ==  211 && code_b ==  321) || (code_a == -211 && code_b == -321));
        bitsPionKaonPair.SetBitNumber(1, (code_a ==  211 && code_b == -321) || (code_a == -211 && code_b ==  321));
        bitsPionKaonPair.SetBitNumber(2, (code_a ==  211 && code_b ==  321));
        bitsPionKaonPair.SetBitNumber(3, (code_a ==  -211 && code_b ==  -321));
        bitsPionKaonPair.SetBitNumber(4, (code_a ==  211 && code_b ==  -321));
        bitsPionKaonPair.SetBitNumber(5, (code_a ==  -211 && code_b ==  321));
        for (int iBits = 0; iBits < nBits; iBits++) {
         double weight_all = pid_weight_a * pid_weight_b;
          if ((iBits < 2) && bitsPionKaonPair.TestBitNumber(iBits)) {
           fProfileDeltaPionKaon[iBits] -> Fill(fCent, delta, weight_all);
           fProfileGammaPionKaon[0][iBits] -> Fill(fCent, gamma[0], weight_all);
           if(isUseVZEROPlane) {
             fProfileGammaPionKaon[1][iBits] -> Fill(fCentSPD1, gamma[1], weight_all);
             fProfileGammaPionKaon[2][iBits] -> Fill(fCentSPD1, gamma[2], weight_all);
           }
           if (isUseZDCPlane) {
             fProfileGammaPionKaon[3][iBits] -> Fill(fCent, gamma[3], weight_all);
             fProfileGammaPionKaon[4][iBits] -> Fill(fCent, gamma[4], weight_all);
           }
           if (isCalculateDiffResult) {
             fProfile2DiffDeltaPionKaonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, delta, weight_all);
             fProfile2DiffDeltaPionKaonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, delta, weight_all);
             fProfile2DiffDeltaPionKaonDEta[iBits] -> Fill(fCent, eta_a - eta_b, delta, weight_all);
             fProfile2DiffDeltaPionKaonQ2[iBits]   -> Fill(fCent,        fQ2Bin, delta, weight_all);
             fProfile2DiffGammaPionKaonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, gamma[1], weight_all);
             fProfile2DiffGammaPionKaonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, gamma[1], weight_all);
             fProfile2DiffGammaPionKaonDEta[iBits] -> Fill(fCent, eta_a - eta_b, gamma[1], weight_all);
             fProfile2DiffGammaPionKaonQ2[iBits]   -> Fill(fCent,        fQ2Bin, gamma[1], weight_all);
           }
           if (isCalculateDeltaPhiSumPhi) {
             double dphi_ranged = RangeDPhi(phi_a - phi_b);
             double sphi_ranged = RangeDPhi(phi_a + phi_b - 2 * psi2TPC_forThisPair);
             fHist2DEtaDPhiPionKaon[fCentBin][iBits] -> Fill(eta_a - eta_b, dphi_ranged);
             fHist2DEtaSPhiPionKaon[fCentBin][iBits] -> Fill(eta_a - eta_b, sphi_ranged);
           }
         }
         else if (isOpenSsandOsSelfCheck && (iBits >= 2) && bitsPionKaonPair.TestBitNumber(iBits)) {
            fProfileDeltaPionKaonSplit[iBits-2] -> Fill(fCent, delta, weight_all);
            fProfileGammaPionKaonSplit[0][iBits-2] -> Fill(fCent, gamma[0], weight_all);
            if(isUseVZEROPlane) {
              fProfileGammaPionKaonSplit[1][iBits-2] -> Fill(fCentSPD1, gamma[1], weight_all);
              fProfileGammaPionKaonSplit[2][iBits-2] -> Fill(fCentSPD1, gamma[2], weight_all);
            }
            if (isUseZDCPlane) {
              fProfileGammaPionKaonSplit[3][iBits-2] -> Fill(fCent, gamma[3], weight_all);
              fProfileGammaPionKaonSplit[4][iBits-2] -> Fill(fCent, gamma[4], weight_all);
            }
          }
        }
      }

      if (isCalculatePionProton) {
        TBits bitsPionProtonPair(nBits);
        bitsPionProtonPair.SetBitNumber(0, (code_a ==  211 && code_b ==  2212) || (code_a == -211 && code_b == -2212));
        bitsPionProtonPair.SetBitNumber(1, (code_a ==  211 && code_b == -2212) || (code_a == -211 && code_b ==  2212));
        bitsPionProtonPair.SetBitNumber(2, (code_a ==  211 && code_b ==  2212));
        bitsPionProtonPair.SetBitNumber(3, (code_a ==  -211 && code_b == -2212));
        bitsPionProtonPair.SetBitNumber(4, (code_a ==  211 && code_b ==  -2212));
        bitsPionProtonPair.SetBitNumber(5, (code_a ==  -211 && code_b == 2212));
        for (int iBits = 0; iBits < nBits; iBits++) {
         double weight_all = pid_weight_a * pid_weight_b;
          if ((iBits < 2) && bitsPionProtonPair.TestBitNumber(iBits)) {
           fProfileDeltaPionProton[iBits] -> Fill(fCent, delta, weight_all);
           fProfileGammaPionProton[0][iBits] -> Fill(fCent, gamma[0], weight_all);
           if(isUseVZEROPlane) {
             fProfileGammaPionProton[1][iBits] -> Fill(fCentSPD1, gamma[1], weight_all);
             fProfileGammaPionProton[2][iBits] -> Fill(fCentSPD1, gamma[2], weight_all);
           }
           if (isUseZDCPlane) {
             fProfileGammaPionProton[3][iBits] -> Fill(fCent, gamma[3], weight_all);
             fProfileGammaPionProton[4][iBits] -> Fill(fCent, gamma[4], weight_all);
           }
           if (isCalculateDiffResult) {
             fProfile2DiffDeltaPionProtonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, delta, weight_all);
             fProfile2DiffDeltaPionProtonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, delta, weight_all);
             fProfile2DiffDeltaPionProtonDEta[iBits] -> Fill(fCent, eta_a - eta_b, delta, weight_all);
             fProfile2DiffDeltaPionProtonQ2[iBits]   -> Fill(fCent,        fQ2Bin, delta, weight_all);
             fProfile2DiffGammaPionProtonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, gamma[1], weight_all);
             fProfile2DiffGammaPionProtonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, gamma[1], weight_all);
             fProfile2DiffGammaPionProtonDEta[iBits] -> Fill(fCent, eta_a - eta_b, gamma[1], weight_all);
             fProfile2DiffGammaPionProtonQ2[iBits]   -> Fill(fCent,        fQ2Bin, gamma[1], weight_all);
           }
           if (isCalculateDeltaPhiSumPhi) {
             double dphi_ranged = RangeDPhi(phi_a - phi_b);
             double sphi_ranged = RangeDPhi(phi_a + phi_b - 2 * psi2TPC_forThisPair);
             fHist2DEtaDPhiPionProton[fCentBin][iBits] -> Fill(eta_a - eta_b, dphi_ranged);
             fHist2DEtaSPhiPionProton[fCentBin][iBits] -> Fill(eta_a - eta_b, sphi_ranged);
           }
         }
         else if (isOpenSsandOsSelfCheck && (iBits >= 2) && bitsPionProtonPair.TestBitNumber(iBits)){
            fProfileDeltaPionProtonSplit[iBits-2] -> Fill(fCent, delta, weight_all);
            fProfileGammaPionProtonSplit[0][iBits-2] -> Fill(fCent, gamma[0], weight_all);
            if(isUseVZEROPlane) {
              fProfileGammaPionProtonSplit[1][iBits-2] -> Fill(fCentSPD1, gamma[1], weight_all);
              fProfileGammaPionProtonSplit[2][iBits-2] -> Fill(fCentSPD1, gamma[2], weight_all);
            }
            if (isUseZDCPlane) {
              fProfileGammaPionProtonSplit[3][iBits-2] -> Fill(fCent, gamma[3], weight_all);
              fProfileGammaPionProtonSplit[4][iBits-2] -> Fill(fCent, gamma[4], weight_all);
            }
          }
        }
      }

      if (isCalculateKaonProton) {
        TBits bitsKaonProtonPair(nBits);
        bitsKaonProtonPair.SetBitNumber(0, (code_a ==  321 && code_b ==  2212) || (code_a == -321 && code_b == -2212));
        bitsKaonProtonPair.SetBitNumber(1, (code_a ==  321 && code_b == -2212) || (code_a == -321 && code_b ==  2212));
        bitsKaonProtonPair.SetBitNumber(2, (code_a ==  321 && code_b ==  2212));
        bitsKaonProtonPair.SetBitNumber(3, (code_a ==  -321 && code_b == -2212));
        bitsKaonProtonPair.SetBitNumber(4, (code_a ==  321 && code_b ==  -2212));
        bitsKaonProtonPair.SetBitNumber(5, (code_a ==  -321 && code_b == 2212));
        for (int iBits = 0; iBits < nBits; iBits++) {
         double weight_all = pid_weight_a * pid_weight_b;
          if ((iBits < 2) && bitsKaonProtonPair.TestBitNumber(iBits)) {
           fProfileDeltaKaonProton[iBits] -> Fill(fCent, delta, weight_all);
           fProfileGammaKaonProton[0][iBits] -> Fill(fCent, gamma[0], weight_all);
           if(isUseVZEROPlane) {
             fProfileGammaKaonProton[1][iBits] -> Fill(fCentSPD1, gamma[1], weight_all);
             fProfileGammaKaonProton[2][iBits] -> Fill(fCentSPD1, gamma[2], weight_all);
           }
           if (isUseZDCPlane) {
             fProfileGammaKaonProton[3][iBits] -> Fill(fCent, gamma[3], weight_all);
             fProfileGammaKaonProton[4][iBits] -> Fill(fCent, gamma[4], weight_all);
           }
           if (isCalculateDiffResult) {
             fProfile2DiffDeltaKaonProtonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, delta, weight_all);
             fProfile2DiffDeltaKaonProtonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, delta, weight_all);
             fProfile2DiffDeltaKaonProtonDEta[iBits] -> Fill(fCent, eta_a - eta_b, delta, weight_all);
             fProfile2DiffDeltaKaonProtonQ2[iBits]   -> Fill(fCent,        fQ2Bin, delta, weight_all);
             fProfile2DiffGammaKaonProtonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, gamma[1], weight_all);
             fProfile2DiffGammaKaonProtonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, gamma[1], weight_all);
             fProfile2DiffGammaKaonProtonDEta[iBits] -> Fill(fCent, eta_a - eta_b, gamma[1], weight_all);
             fProfile2DiffGammaKaonProtonQ2[iBits]   -> Fill(fCent,        fQ2Bin, gamma[1], weight_all);
           }
           if (isCalculateDeltaPhiSumPhi) {
             double dphi_ranged = RangeDPhi(phi_a - phi_b);
             double sphi_ranged = RangeDPhi(phi_a + phi_b - 2 * psi2TPC_forThisPair);
             fHist2DEtaDPhiKaonProton[fCentBin][iBits] -> Fill(eta_a - eta_b, dphi_ranged);
             fHist2DEtaSPhiKaonProton[fCentBin][iBits] -> Fill(eta_a - eta_b, sphi_ranged);
           }
         }
         else if(isOpenSsandOsSelfCheck && (iBits >= 2) && bitsKaonProtonPair.TestBitNumber(iBits)){
            fProfileDeltaKaonProtonSplit[iBits-2] -> Fill(fCent, delta, weight_all);
            fProfileGammaKaonProtonSplit[0][iBits-2] -> Fill(fCent, gamma[0], weight_all);
            if(isUseVZEROPlane) {
              fProfileGammaKaonProtonSplit[1][iBits-2] -> Fill(fCentSPD1, gamma[1], weight_all);
              fProfileGammaKaonProtonSplit[2][iBits-2] -> Fill(fCentSPD1, gamma[2], weight_all);
            }
            if (isUseZDCPlane) {
              fProfileGammaKaonProtonSplit[3][iBits-2] -> Fill(fCent, gamma[3], weight_all);
              fProfileGammaKaonProtonSplit[4][iBits-2] -> Fill(fCent, gamma[4], weight_all);
            }
          }
        }
      }

      if (isCalculatePionPion) {
        TBits bitsPionPionPair(nBits);
        bitsPionPionPair.SetBitNumber(0, (code_a ==  211 && code_b ==  211) || (code_a == -211 && code_b == -211));
        bitsPionPionPair.SetBitNumber(1, (code_a ==  211 && code_b == -211) || (code_a == -211 && code_b ==  211));
        bitsPionPionPair.SetBitNumber(2, (code_a ==  211 && code_b ==  211));
        bitsPionPionPair.SetBitNumber(3, (code_a ==  -211 && code_b == -211));
        bitsPionPionPair.SetBitNumber(4, (code_a ==  211 && code_b == -211));
        bitsPionPionPair.SetBitNumber(5, (code_a ==  -211 && code_b == 211));
        for (int iBits = 0; iBits < nBits; iBits++) {
         double weight_all = pid_weight_a * pid_weight_b;
          if ((iBits < 2 ) && bitsPionPionPair.TestBitNumber(iBits)) {
           fProfileDeltaPionPion[iBits] -> Fill(fCent, delta, weight_all);
           fProfileGammaPionPion[0][iBits] -> Fill(fCent, gamma[0], weight_all);
           if(isUseVZEROPlane) {
             fProfileGammaPionPion[1][iBits] -> Fill(fCentSPD1, gamma[1], weight_all);
             fProfileGammaPionPion[2][iBits] -> Fill(fCentSPD1, gamma[2], weight_all);
           }
           if (isUseZDCPlane) {
             fProfileGammaPionPion[3][iBits] -> Fill(fCent, gamma[3], weight_all);
             fProfileGammaPionPion[4][iBits] -> Fill(fCent, gamma[4], weight_all);
           }
           if (isCalculateDiffResult) {
             fProfile2DiffDeltaPionPionDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, delta, weight_all);
             fProfile2DiffDeltaPionPionSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, delta, weight_all);
             fProfile2DiffDeltaPionPionDEta[iBits] -> Fill(fCent, eta_a - eta_b, delta, weight_all);
             fProfile2DiffDeltaPionPionQ2[iBits]   -> Fill(fCent,        fQ2Bin, delta, weight_all);
             fProfile2DiffGammaPionPionDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, gamma[1], weight_all);
             fProfile2DiffGammaPionPionSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, gamma[1], weight_all);
             fProfile2DiffGammaPionPionDEta[iBits] -> Fill(fCent, eta_a - eta_b, gamma[1], weight_all);
             fProfile2DiffGammaPionPionQ2[iBits]   -> Fill(fCent,        fQ2Bin, gamma[1], weight_all);
           }
           if (isCalculateDeltaPhiSumPhi) {
             double dphi_ranged = RangeDPhi(phi_a - phi_b);
             double sphi_ranged = RangeDPhi(phi_a + phi_b - 2 * psi2TPC_forThisPair);
             fHist2DEtaDPhiPionPion[fCentBin][iBits] -> Fill(eta_a - eta_b, dphi_ranged);
             fHist2DEtaSPhiPionPion[fCentBin][iBits] -> Fill(eta_a - eta_b, sphi_ranged);
           }
         }
          else if(isOpenSsandOsSelfCheck && (iBits >= 2) && bitsPionPionPair.TestBitNumber(iBits)) {
            fProfileDeltaPionPionSplit[iBits-2] -> Fill(fCent, delta, weight_all);
            fProfileGammaPionPionSplit[0][iBits-2] -> Fill(fCent, gamma[0], weight_all);
            if(isUseVZEROPlane) {
              fProfileGammaPionPionSplit[1][iBits-2] -> Fill(fCentSPD1, gamma[1], weight_all);
              fProfileGammaPionPionSplit[2][iBits-2] -> Fill(fCentSPD1, gamma[2], weight_all);
            }
            if (isUseZDCPlane) {
              fProfileGammaPionPionSplit[3][iBits-2] -> Fill(fCent, gamma[3], weight_all);
              fProfileGammaPionPionSplit[4][iBits-2] -> Fill(fCent, gamma[4], weight_all);
            }
          }
        }
      }

      if (isCalculateKaonKaon) {
        TBits bitsKaonKaonPair(nBits);
        bitsKaonKaonPair.SetBitNumber(0, (code_a ==  321 && code_b ==  321) || (code_a == -321 && code_b == -321));
        bitsKaonKaonPair.SetBitNumber(1, (code_a ==  321 && code_b == -321) || (code_a == -321 && code_b ==  321));
        bitsKaonKaonPair.SetBitNumber(2, (code_a ==  321 && code_b ==  321));
        bitsKaonKaonPair.SetBitNumber(3, (code_a ==  -321 && code_b == -321));
        bitsKaonKaonPair.SetBitNumber(4, (code_a ==  321 && code_b ==  -321));
        bitsKaonKaonPair.SetBitNumber(5, (code_a ==  -321 && code_b == 321));
        for (int iBits = 0; iBits < nBits; iBits++) {
         double weight_all = pid_weight_a * pid_weight_b;
          if ((iBits < 2) && bitsKaonKaonPair.TestBitNumber(iBits)) {
           fProfileDeltaKaonKaon[iBits] -> Fill(fCent, delta, weight_all);
           fProfileGammaKaonKaon[0][iBits] -> Fill(fCent, gamma[0], weight_all);
           if(isUseVZEROPlane) {
             fProfileGammaKaonKaon[1][iBits] -> Fill(fCentSPD1, gamma[1], weight_all);
             fProfileGammaKaonKaon[2][iBits] -> Fill(fCentSPD1, gamma[2], weight_all);
           }
           if (isUseZDCPlane) {
             fProfileGammaKaonKaon[3][iBits] -> Fill(fCent, gamma[3], weight_all);
             fProfileGammaKaonKaon[4][iBits] -> Fill(fCent, gamma[4], weight_all);
           }
           if (isCalculateDiffResult) {
             fProfile2DiffDeltaKaonKaonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, delta, weight_all);
             fProfile2DiffDeltaKaonKaonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, delta, weight_all);
             fProfile2DiffDeltaKaonKaonDEta[iBits] -> Fill(fCent, eta_a - eta_b, delta, weight_all);
             fProfile2DiffDeltaKaonKaonQ2[iBits]   -> Fill(fCent,        fQ2Bin, delta, weight_all);
             fProfile2DiffGammaKaonKaonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, gamma[1], weight_all);
             fProfile2DiffGammaKaonKaonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, gamma[1], weight_all);
             fProfile2DiffGammaKaonKaonDEta[iBits] -> Fill(fCent, eta_a - eta_b, gamma[1], weight_all);
             fProfile2DiffGammaKaonKaonQ2[iBits]   -> Fill(fCent,        fQ2Bin, gamma[1], weight_all);
           }
           if (isCalculateDeltaPhiSumPhi) {
             double dphi_ranged = RangeDPhi(phi_a - phi_b);
             double sphi_ranged = RangeDPhi(phi_a + phi_b - 2 * psi2TPC_forThisPair);
             fHist2DEtaDPhiKaonKaon[fCentBin][iBits] -> Fill(eta_a - eta_b, dphi_ranged);
             fHist2DEtaSPhiKaonKaon[fCentBin][iBits] -> Fill(eta_a - eta_b, sphi_ranged);
           }
         }
         else if(isOpenSsandOsSelfCheck && (iBits >= 2) && bitsKaonKaonPair.TestBitNumber(iBits)){
            fProfileDeltaKaonKaonSplit[iBits-2] -> Fill(fCent, delta, weight_all);
            fProfileGammaKaonKaonSplit[0][iBits-2] -> Fill(fCent, gamma[0], weight_all);
            if(isUseVZEROPlane) {
              fProfileGammaKaonKaonSplit[1][iBits-2] -> Fill(fCentSPD1, gamma[1], weight_all);
              fProfileGammaKaonKaonSplit[2][iBits-2] -> Fill(fCentSPD1, gamma[2], weight_all);
            }
            if (isUseZDCPlane) {
              fProfileGammaKaonKaonSplit[3][iBits-2] -> Fill(fCent, gamma[3], weight_all);
              fProfileGammaKaonKaonSplit[4][iBits-2] -> Fill(fCent, gamma[4], weight_all);
            }
          }
        }
      }

      if (isCalculateProtonProton) {
        TBits bitsProtonProtonPair(nBits);
        bitsProtonProtonPair.SetBitNumber(0, (code_a ==  2212 && code_b ==  2212) || (code_a == -2212 && code_b == -2212));
        bitsProtonProtonPair.SetBitNumber(1, (code_a ==  2212 && code_b == -2212) || (code_a == -2212 && code_b ==  2212));
        bitsProtonProtonPair.SetBitNumber(2, (code_a ==  2212 && code_b ==  2212));
        bitsProtonProtonPair.SetBitNumber(3, (code_a ==  -2212 && code_b == -2212));
        bitsProtonProtonPair.SetBitNumber(4, (code_a ==  2212 && code_b ==  -2212));
        bitsProtonProtonPair.SetBitNumber(5, (code_a ==  -2212 && code_b == 2212));
        for (int iBits = 0; iBits < nBits; iBits++) {
         double weight_all = pid_weight_a * pid_weight_b;
          if ((iBits < 2) && bitsProtonProtonPair.TestBitNumber(iBits)) {
           fProfileDeltaProtonProton[iBits] -> Fill(fCent, delta, weight_all);
           fProfileGammaProtonProton[0][iBits] -> Fill(fCent, gamma[0], weight_all);
           if(isUseVZEROPlane) {
             fProfileGammaProtonProton[1][iBits] -> Fill(fCentSPD1, gamma[1], weight_all);
             fProfileGammaProtonProton[2][iBits] -> Fill(fCentSPD1, gamma[2], weight_all);
           }
           if (isUseZDCPlane) {
             fProfileGammaProtonProton[3][iBits] -> Fill(fCent, gamma[3], weight_all);
             fProfileGammaProtonProton[4][iBits] -> Fill(fCent, gamma[4], weight_all);
           }
           if (isCalculateDiffResult) {
             fProfile2DiffDeltaProtonProtonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, delta, weight_all);
             fProfile2DiffDeltaProtonProtonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, delta, weight_all);
             fProfile2DiffDeltaProtonProtonDEta[iBits] -> Fill(fCent, eta_a - eta_b, delta, weight_all);
             fProfile2DiffDeltaProtonProtonQ2[iBits]   -> Fill(fCent,        fQ2Bin, delta, weight_all);
             fProfile2DiffGammaProtonProtonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, gamma[1], weight_all);
             fProfile2DiffGammaProtonProtonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, gamma[1], weight_all);
             fProfile2DiffGammaProtonProtonDEta[iBits] -> Fill(fCent, eta_a - eta_b, gamma[1], weight_all);
             fProfile2DiffGammaProtonProtonQ2[iBits]   -> Fill(fCent,        fQ2Bin, gamma[1], weight_all);
           }
           if (isCalculateDeltaPhiSumPhi) {
             double dphi_ranged = RangeDPhi(phi_a - phi_b);
             double sphi_ranged = RangeDPhi(phi_a + phi_b - 2 * psi2TPC_forThisPair);
             fHist2DEtaDPhiProtonProton[fCentBin][iBits] -> Fill(eta_a - eta_b, dphi_ranged);
             fHist2DEtaSPhiProtonProton[fCentBin][iBits] -> Fill(eta_a - eta_b, sphi_ranged);
           }
         }
         else if(isOpenSsandOsSelfCheck && (iBits >= 2) && bitsProtonProtonPair.TestBitNumber(iBits)){
           fProfileDeltaProtonProtonSplit[iBits-2] -> Fill(fCent, delta, weight_all);
           fProfileGammaProtonProtonSplit[0][iBits-2] -> Fill(fCent, gamma[0], weight_all);
           if(isUseVZEROPlane) {
             fProfileGammaProtonProtonSplit[1][iBits-2] -> Fill(fCentSPD1, gamma[1], weight_all);
             fProfileGammaProtonProtonSplit[2][iBits-2] -> Fill(fCentSPD1, gamma[2], weight_all);
           }
           if (isUseZDCPlane) {
             fProfileGammaProtonProtonSplit[3][iBits-2] -> Fill(fCent, gamma[3], weight_all);
             fProfileGammaProtonProtonSplit[4][iBits-2] -> Fill(fCent, gamma[4], weight_all);
           }
         }
        }
      }

      if (isCalculateHadronHadron) {
        TBits bitsHadronHadronPair(nBits);
        bitsHadronHadronPair.SetBitNumber(0, (code_a > 0 && code_b > 0) || (code_a < 0 && code_b < 0));
        bitsHadronHadronPair.SetBitNumber(1, (code_a > 0 && code_b < 0) || (code_a < 0 && code_b > 0));
        bitsHadronHadronPair.SetBitNumber(2, (code_a > 0 && code_b > 0));
        bitsHadronHadronPair.SetBitNumber(3, (code_a < 0 && code_b < 0));
        bitsHadronHadronPair.SetBitNumber(4, (code_a > 0 && code_b < 0));
        bitsHadronHadronPair.SetBitNumber(5, (code_a < 0 && code_b > 0));
        for (int iBits = 0; iBits < nBits; iBits++) {
          double weight_all = weight_a * weight_b;
          if ((iBits < 2) && bitsHadronHadronPair.TestBitNumber(iBits)) {
            fProfileDeltaHadronHadron[iBits] -> Fill(fCent, delta, weight_all);
            fProfileGammaHadronHadron[0][iBits] -> Fill(fCent, gamma[0], weight_all);
            if(isUseVZEROPlane) {
              fProfileGammaHadronHadron[1][iBits] -> Fill(fCentSPD1, gamma[1], weight_all);
              fProfileGammaHadronHadron[2][iBits] -> Fill(fCentSPD1, gamma[2], weight_all);
            }
            if (isUseZDCPlane) {
              fProfileGammaHadronHadron[3][iBits] -> Fill(fCent, gamma[3], weight_all);
              fProfileGammaHadronHadron[4][iBits] -> Fill(fCent, gamma[4], weight_all);
            }
            if (isCalculateDiffResult) {
              fProfile2DiffDeltaHadronHadronDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, delta, weight_all);
              fProfile2DiffDeltaHadronHadronSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, delta, weight_all);
              fProfile2DiffDeltaHadronHadronDEta[iBits] -> Fill(fCent, eta_a - eta_b, delta, weight_all);
              fProfile2DiffDeltaHadronHadronQ2[iBits]   -> Fill(fCent,        fQ2Bin, delta, weight_all);
              fProfile2DiffGammaHadronHadronDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, gamma[1], weight_all);
              fProfile2DiffGammaHadronHadronSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, gamma[1], weight_all);
              fProfile2DiffGammaHadronHadronDEta[iBits] -> Fill(fCent, eta_a - eta_b, gamma[1], weight_all);
              fProfile2DiffGammaHadronHadronQ2[iBits]   -> Fill(fCent,        fQ2Bin, gamma[1], weight_all);
            }
            if (isCalculateDeltaPhiSumPhi) {
              double dphi_ranged = RangeDPhi(phi_a - phi_b);
              double sphi_ranged = RangeDPhi(phi_a + phi_b - 2 * psi2TPC_forThisPair);
              fHist2DEtaDPhiHadronHadron[fCentBin][iBits] -> Fill(eta_a - eta_b, dphi_ranged);
              fHist2DEtaSPhiHadronHadron[fCentBin][iBits] -> Fill(eta_a - eta_b, sphi_ranged);
            }
          }
          else if(isOpenSsandOsSelfCheck && (iBits >= 2) && bitsHadronHadronPair.TestBitNumber(iBits)){
            fProfileDeltaHadronHadronSplit[iBits-2] -> Fill(fCent, delta, weight_all);
            fProfileGammaHadronHadronSplit[0][iBits-2] -> Fill(fCent, gamma[0], weight_all);
            if(isUseVZEROPlane) {
              fProfileGammaHadronHadronSplit[1][iBits-2] -> Fill(fCentSPD1, gamma[1], weight_all);
              fProfileGammaHadronHadronSplit[2][iBits-2] -> Fill(fCentSPD1, gamma[2], weight_all);
            }
            if (isUseZDCPlane) {
              fProfileGammaHadronHadronSplit[3][iBits-2] -> Fill(fCent, gamma[3], weight_all);
              fProfileGammaHadronHadronSplit[4][iBits-2] -> Fill(fCent, gamma[4], weight_all);
            }
          }
        }
      }
    }
  }
  return true;
}

//---------------------------------------------------

void AliAnalysisTaskPIDCME::ResetVectors()
{
  fSumQ2xTPCPos = 0.;
  fSumQ2yTPCPos = 0.;
  fWgtMultTPCPos = 0.;
  fSumQ2xTPCNeg = 0.;
  fSumQ2yTPCNeg = 0.;
  fWgtMultTPCNeg = 0.;
  fSumQ2xTPC = 0.;
  fSumQ2yTPC = 0.;
  fWgtMultTPC = 0.;
  std::unordered_map<int, std::vector<double>>().swap(mapTPCPosTrksIDPhiWgt);
  std::unordered_map<int, std::vector<double>>().swap(mapTPCNegTrksIDPhiWgt);
  std::unordered_map<int, std::vector<double>>().swap(mapTPCTrksIDPhiWgt);
  std::vector<std::array<double,7>>().swap(vecParticle);
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::LoadCalibHistForThisRun()
{
  if (fPeriod.EqualTo("LHC10h")) {
    // 10h VZERO Calibration Histograms is Global
    // 10h ZDC Calibration Histograms
    if (isUseZDCPlane) {
      tree -> Reset();
      fProfileForZNCGE    ->Reset();
      fProfileForZNAGE    ->Reset();
      fHn4DForZNCQxRC     ->Reset();
      fHn4DForZNCQyRC     ->Reset();
      fHn4DForZNCMtRC     ->Reset();
      fHn4DForZNAQxRC     ->Reset();
      fHn4DForZNAQyRC     ->Reset();
      fHn4DForZNAMtRC     ->Reset();
      fHn4DForZNCCountsRC ->Reset();
      fHn4DForZNACountsRC ->Reset();
      fProfile2DForCosC   ->Reset();
      fProfile2DForSinC   ->Reset();
      fProfile2DForCosA   ->Reset();
      fProfile2DForSinA   ->Reset();

      float vtxQuant1Tmp[3];
      float vtxQuant2Tmp[3];
      tree = (TTree*)fListZDCCalib->FindObject(Form("run%d",fRunNum));
      tree->SetBranchAddress("VtxQuant1",vtxQuant1Tmp);
      tree->SetBranchAddress("VtxQuant2",vtxQuant2Tmp);
      tree->GetEntry(0);
      for (int i = 0; i < 3; i++) {
        vtxQuant1[i] = vtxQuant1Tmp[i];
        vtxQuant2[i] = vtxQuant2Tmp[i];
      }
      fProfileForZNCGE    = (TProfile*)fListZDCCalib   -> FindObject(Form("ZNCTowerEnegryMean_run%d",fRunNum));
      fProfileForZNAGE    = (TProfile*)fListZDCCalib   -> FindObject(Form("ZNATowerEnegryMean_run%d",fRunNum));
      fHn4DForZNCQxRC     = (THnSparseF*)fListZDCCalib -> FindObject(Form("Hn4DQxZNCCentVxVyVz_run%d",fRunNum));
      fHn4DForZNCQyRC     = (THnSparseF*)fListZDCCalib -> FindObject(Form("Hn4DQyZNCCentVxVyVz_run%d",fRunNum));
      fHn4DForZNCMtRC     = (THnSparseF*)fListZDCCalib -> FindObject(Form("Hn4DMtZNCCentVxVyVz_run%d",fRunNum));
      fHn4DForZNAQxRC     = (THnSparseF*)fListZDCCalib -> FindObject(Form("Hn4DQxZNACentVxVyVz_run%d",fRunNum));
      fHn4DForZNAQyRC     = (THnSparseF*)fListZDCCalib -> FindObject(Form("Hn4DQyZNACentVxVyVz_run%d",fRunNum));
      fHn4DForZNAMtRC     = (THnSparseF*)fListZDCCalib -> FindObject(Form("Hn4DMtZNACentVxVyVz_run%d",fRunNum));
      fHn4DForZNCCountsRC = (THnSparseI*)fListZDCCalib -> FindObject(Form("Hn4DCountsZNCCentVxVyVz_run%d",fRunNum));
      fHn4DForZNACountsRC = (THnSparseI*)fListZDCCalib -> FindObject(Form("Hn4DCountsZNACentVxVyVz_run%d",fRunNum));
      fProfile2DForCosC   = (TProfile2D*)fListZDCCalib -> FindObject(Form("Profile2DShiftCosC_run%d",fRunNum));
      fProfile2DForSinC   = (TProfile2D*)fListZDCCalib -> FindObject(Form("Profile2DShiftSinC_run%d",fRunNum));
      fProfile2DForCosA   = (TProfile2D*)fListZDCCalib -> FindObject(Form("Profile2DShiftCosA_run%d",fRunNum));
      fProfile2DForSinA   = (TProfile2D*)fListZDCCalib -> FindObject(Form("Profile2DShiftSinA_run%d",fRunNum));
      if (!fProfileForZNCGE)    return false;
      if (!fProfileForZNAGE)    return false;
      if (!fHn4DForZNCQxRC)     return false;
      if (!fHn4DForZNCQyRC)     return false;
      if (!fHn4DForZNCMtRC)     return false;
      if (!fHn4DForZNAQxRC)     return false;
      if (!fHn4DForZNAQyRC)     return false;
      if (!fHn4DForZNAMtRC)     return false;
      if (!fHn4DForZNCCountsRC) return false;
      if (!fHn4DForZNACountsRC) return false;
      if (!fProfile2DForCosC)   return false;
      if (!fProfile2DForSinC)   return false;
      if (!fProfile2DForCosA)   return false;
      if (!fProfile2DForSinA)   return false;
    }
  }

  if (fPeriod.EqualTo("LHC15o")) {
    // 15o VZERO Calibration Histograms
    if (isUseVZEROPlane) {
      hMultV0 -> Reset();
      for (int i = 0; i < 2; i++) {
        hQx2mV0[i] -> Reset();
        hQy2mV0[i] -> Reset();
      }
      hMultV0    = ((TH1D*) contMult ->GetObject(fRunNum));
      hQx2mV0[0] = ((TH1D*) contQxncm->GetObject(fRunNum));
      hQy2mV0[0] = ((TH1D*) contQyncm->GetObject(fRunNum));
      hQx2mV0[1] = ((TH1D*) contQxnam->GetObject(fRunNum));
      hQy2mV0[1] = ((TH1D*) contQynam->GetObject(fRunNum));
      if (!hMultV0)    return false;
      if (!hQx2mV0[0]) return false;
      if (!hQy2mV0[0]) return false;
      if (!hQx2mV0[1]) return false;
      if (!hQy2mV0[1]) return false;
    }
    //15o NUA
    if (isDoNUA) {
      hCorrectNUAPos -> Reset();
      hCorrectNUANeg -> Reset();
      hCorrectNUAPos = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Pos_Cent0_Run%d",fRunNum));
      hCorrectNUANeg = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Neg_Cent0_Run%d",fRunNum));
      if (!hCorrectNUAPos) return false;
      if (!hCorrectNUANeg) return false;
    }
  }

  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    //18q/r VZERO
    if (isUseVZEROPlane) {
      for (int i = 0; i < 2; i++) {
        hQx2mV0[i] ->Reset();
        hQy2mV0[i] ->Reset();
      }
      hQx2mV0[0] = ((TH1D*) contQxncm->GetObject(fRunNum));
      hQy2mV0[0] = ((TH1D*) contQyncm->GetObject(fRunNum));
      hQx2mV0[1] = ((TH1D*) contQxnam->GetObject(fRunNum));
      hQy2mV0[1] = ((TH1D*) contQynam->GetObject(fRunNum));
      for (int i = 0; i < 2; i++) {
        if (!hQx2mV0[i]) return false;
        if (!hQy2mV0[i]) return false;
      }
      fHCorrectV0ChWeghts -> Reset();
      fHCorrectV0ChWeghts = (TH2F *) fListVZEROCalib->FindObject(Form("hWgtV0ChannelsvsVzRun%d",fRunNum));
      if (!fHCorrectV0ChWeghts) return false;
    }
    //18q/r NUA
    if (isDoNUA) {
      hCorrectNUAPos -> Reset();
      hCorrectNUANeg -> Reset();
      hCorrectNUAPos = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d",0,fRunNum));
      hCorrectNUANeg = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d",0,fRunNum));
      if (!hCorrectNUAPos) return false;
      if (!hCorrectNUANeg) return false;
    }
    //18q/r ZDC
    if (isUseZDCPlane) {
      fHZDCCparameters -> Reset();
      fHZDCAparameters -> Reset();
      fHZDCCparameters = (TH1D*)(fListZDCCalib->FindObject(Form("Run %d", fRunNum))->FindObject(Form("fZDCCparameters[%d]",fRunNum)));
      fHZDCAparameters = (TH1D*)(fListZDCCalib->FindObject(Form("Run %d", fRunNum))->FindObject(Form("fZDCAparameters[%d]",fRunNum)));
      if(fHZDCCparameters && fHZDCAparameters) std::cout<<"\n ===========> Info:: ZDC Channel Weights Found for Run "<<fRunNum<<std::endl;
      if (!fHZDCCparameters) return false;
      if (!fHZDCAparameters) return false;
    }
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::RemovalForRun1()
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

bool AliAnalysisTaskPIDCME::RejectEvtMultComp() // 15o_pass1, old pile-up
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
      // double Phi  = aodTrk->Phi();
      if (Pt<0.2 || Pt>5.0 || TMath::Abs(Eta)>0.8 || aodTrk->GetTPCNcls()<fNclsCut || aodTrk->GetTPCsignal()<10.0) continue;
      if (aodTrk->TestFilterBit(1) && aodTrk->Chi2perNDF()>0.2)  multTPCFE++;
      if (!aodTrk->TestFilterBit(16) || aodTrk->Chi2perNDF()<0.1)   continue;
      Double_t dca[2] = {-99., -99.};
      Double_t cov[3] = {-99., -99., -99.};
      Double_t magField = fAOD->GetMagneticField();
      if (magField!=0) {
        if (aodTrk->PropagateToDCA(fAOD->GetPrimaryVertex(), magField, 100., dca, cov) && TMath::Abs(dca[0]) < 0.3 && TMath::Abs(dca[1]) < 0.3) multGlobal++;
      }
    }

    fHist2MultMultQA[0]->Fill(multTPC,multEsd);
    fHist2MultMultQA[1]->Fill(multGlobal,multTPCFE);

    TString fMultComp = "pileupByGlobalTPC1";

    if (fMultComp.EqualTo("pileupByEDSTPC128")) { // Rihan
      if ((Double_t)(multEsd*(1/3.45) - 90.) < (Double_t)multTPC)
      {
        fHist2MultMultQA[3]->Fill(multTPC,multEsd);
        return true;
      }
      else return false;
    }

    if (fMultComp.EqualTo("pileupByGlobalTPC1")) { // A.Dobrin
      if (multTPCFE-1.78*multGlobal<62.87 && multTPCFE-1.48*multGlobal>-36.73) {
        fHist2MultMultQA[4]->Fill(multGlobal,multTPCFE);
        return true;
      }
      else return false;
    }

    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::RejectEvtTFFit()
{
  Int_t nITSClsLy0 = fAOD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fAOD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTrk = 0;
  for (Int_t it = 0; it < nTracks; it++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
      if (!aodTrk) {
          delete aodTrk;
          continue;
      }
      if (aodTrk->TestFilterBit(32)) {
        if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
      }
  }

  fHist2MultCentQA[0]->Fill(fCentV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)

  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  Float_t  multV0a = aodV0->GetMTotV0A();
  Float_t  multV0c = aodV0->GetMTotV0C();
  Float_t  multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;

  // pile-up cuts
  if (fCentSPD0 < fCenCutLowPU->Eval(fCentV0M)) return false;
  if (fCentSPD0 > fCenCutHighPU->Eval(fCentV0M)) return false;
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) return false;
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;
  if (Float_t(multTrk) < fMultCutPU->Eval(fCentV0M)) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fAOD->IsIncompleteDAQ()) return false;

  if (isTightPileUp == true){
    Int_t tpcClsTot = fAOD->GetNumberOfTPCClusters();
    Float_t nclsDif = Float_t(tpcClsTot) - (53182.6 + 113.326*multV0Tot - 0.000831275*multV0Tot*multV0Tot);
    if (nclsDif > 200000)//can be varied to 150000, 200000
    return false;
  }

  fHist2MultCentQA[1]->Fill(fCentV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::RejectEvtTPCITSfb32TOF ()
{
    //TOD+FB32 pile-up removal
    // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsEventProp
    Int_t multTrk=0;
    Int_t multTrkTOF=0;
    int nTrk = fAOD->GetNumberOfTracks();
    for (Int_t it2 = 0; it2 < nTrk; it2++) {
      AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it2);
      if (!aodTrk) continue;
      if (aodTrk->TestFilterBit(32)) {
        multTrk++;
        if ( TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10 && aodTrk->GetTOFsignal() >= 12000 && aodTrk->GetTOFsignal() <= 25000) multTrkTOF++;
        else return false;
      }
    }
    fHist2MultMultQA[2]->Fill(multTrkTOF, nTrk);
    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::AODPileupCheck()
{
  Int_t isPileup = fAOD->IsPileupFromSPD(3);
  if (isPileup !=0 && fPeriod.EqualTo("LHC16t")) return false; // LHC16t : pPb
  if (fAOD->IsIncompleteDAQ()) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fPeriod.EqualTo("LHC15o")) {
    if (!fMultSel->GetThisEventIsNotPileup()) return false;
    if (!fMultSel->GetThisEventIsNotPileupMV()) return false;
    if (!fMultSel->GetThisEventIsNotPileupInMultBins()) return false;
    if (!fMultSel->GetThisEventHasNoInconsistentVertices()) return false;
    if (!fMultSel->GetThisEventPassesTrackletVsCluster()) return false;
    if (!fMultSel->GetThisEventIsNotIncompleteDAQ()) return false;
    if (!fMultSel->GetThisEventHasGoodVertex2016()) return false;
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::PileUpMultiVertex()
{
  // check for multi-vertexer pile-up
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if (!(nPlp=fAOD->GetNumberOfPileupVerticesTracks()))
  return false;

  vtPrm = fAOD->GetPrimaryVertex();
  if (vtPrm == fAOD->GetPrimaryVertexSPD())
  return true;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for (int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)fAOD->GetPileupVertexTracks(ipl);
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF()    > kMaxPlpChi2)    continue;
    //int bcPlp = vtPlp->GetBC();
    //if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2)
    // return kTRUE; // pile-up from other BC

    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist) continue;

    return true; // pile-up: well separated vertices
  }
  return false;
}

//---------------------------------------------------

double AliAnalysisTaskPIDCME::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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

bool AliAnalysisTaskPIDCME::AcceptAODTrack(AliAODTrack *track)
{
    //------------------
    // track cut
    //------------------
    double pt = track->Pt();
    if(pt < fPtMin) return false;
    if(pt > fPtMax) return false;

    double eta = track->Eta();
    if(fabs(eta)  > fEtaCut) return false;

    int nhits  = track->GetTPCNcls();
    if(fabs(nhits) < fNclsCut) return false;

    double dedx = track->GetTPCsignal();
    if(dedx < fDedxCut) return false;

    double chi2   = track->Chi2perNDF();
    if(chi2 < fChi2Min) return false;
    if(chi2 > fChi2Max) return false;

    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck)
{
  if (pidToCheck==0) return kTRUE;    //// Charge Particles do not need PID check

  if (!fPIDResponse) {
    Printf("\n Could Not access PIDResponse Task, please Add the Task...\n return with kFALSE pid\n");
    return kFALSE;
  }

  /// Rihan todo: To set the low pT cuts for nSigmaTPC from AddTaskMacro!
  /// Although someone barely needs to change it given the purity..

  float nSigTPC = 0, nSigTOF = 0, nSigRMS = 0;
  double trkPtPID  = ftrack->Pt();

  ///Pion =>
  if (pidToCheck==1) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);//Some warning show here (***TDatabasePDG::AddParicle: particle with PDGcode = 3124 already defind),I don't understand what happended. --chunzheng
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if (trkPtPID<=fPtMinforRMSPion && TMath::Abs(nSigTPC)<=fNSigmaTPCCutPion) return true;
    if (trkPtPID> fPtMinforRMSPion && TMath::Abs(nSigRMS)<=fNSigmaRMSCutPion) return true;
    return false;
  }
  ///Kaon =>
  else if (pidToCheck==2) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kKaon);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kKaon);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if (trkPtPID<=fPtMinforRMSKaon && TMath::Abs(nSigTPC)<=fNSigmaTPCCutKaon) return true;
    if (trkPtPID> fPtMinforRMSKaon && TMath::Abs(nSigRMS)<=fNSigmaRMSCutKaon) return true;
    return false;
  }
  ///proton =>
  else if (pidToCheck==3) {///
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kProton);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kProton);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
    
    bool isProton = false;
    if (trkPtPID <= fPtMinforRMSProton && TMath::Abs(nSigTPC) <= fNSigmaTPCCutProton) isProton = true;
    if (trkPtPID >  fPtMinforRMSProton && TMath::Abs(nSigRMS) <= fNSigmaRMSCutProton) isProton = true;

    if(isUsePionRejection) {
      float nSigTPCPion = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);
      float nSigTOFPion = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
      float nSigRMSPion = TMath::Sqrt(nSigTPCPion*nSigTPCPion + nSigTOFPion*nSigTOFPion);
      double mom = ftrack->P();
      bool isPassPionRejection = (mom > 0.5 && TMath::Abs(nSigRMSPion) > 3.);
      isPassPionRejection = isPassPionRejection || (mom < 0.5 && TMath::Abs(nSigTPCPion) > 3.);
      isProton = isProton && isPassPionRejection;
    }
    return isProton;
  } else {
    Printf("\n -Ve number not allowed! Choose among: 0,1,2,3 (Charge Pion, Kaon, Proton)\n return with kFALSE \n");
    return false;
  }
}

//---------------------------------------------------

double AliAnalysisTaskPIDCME::GetNUECor(int charge, double pt)
{
  if(!fListNUE) return -1;
  if(charge == 0) return -1;
  TString histName;
  if (fPeriod.EqualTo("LHC10h")) {
    histName = (charge > 0) ? Form("effVsPt_cent%iPlus", fCentBin) : Form("effVsPt_cent%iMinus", fCentBin);
  }
  else if (fPeriod.EqualTo("LHC15o")) {
    histName = (charge > 0) ? "trkEfficiencyChrgPos" : "trkEfficiencyChrgNeg";
  } else if(fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    histName = (charge > 0) ? "trkEfficiencyChrgPos" : "trkEfficiencyChrgNeg";
  }
  else return -1;


  TH1D* efficiencyHist = (TH1D*)fListNUE->FindObject(histName);
  if (!efficiencyHist) return -1;
  int ptBin = efficiencyHist->GetXaxis()->FindBin(pt);
  double binContent = efficiencyHist->GetBinContent(ptBin);

  if (binContent > 1.e-5) {
    return 1.0 / binContent;
  } else return -1;
}


//---------------------------------------------------
double AliAnalysisTaskPIDCME::GetPIDNUECor(int pdgcode, double pt)
{
  if(!fListNUE) return -1;
  TString histName;
  if (fPeriod.EqualTo("LHC10h")) {
    std::cout<< "No PID NUE correction for LHC10h" << std::endl;
    return -1;
  }
  else if (fPeriod.EqualTo("LHC15o")) {
    if (pdgcode == 211) histName = "trkEfficiencyPionPos";
    else if (pdgcode == 321) histName = "trkEfficiencyKaonPos";
    else if (pdgcode == 2212) histName = "trkEfficiencyProtPos";
    else if (pdgcode == -221) histName = "trkEfficiencyPionNeg";
    else if (pdgcode == -321) histName = "trkEfficiencyKaonNeg";
    else if (pdgcode == -2212) histName = "trkEfficiencyProtNeg";
    else return -1;
  } 
  else if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    if (pdgcode == 2212)       histName = "h_eff_pos_proton";
    else if (pdgcode == -2212) histName = "h_eff_neg_proton";
    else return -1;
  }

  TH1D* efficiencyHist = (TH1D*)fListNUE->FindObject(histName);
  if (!efficiencyHist) return -1;
  int ptBin = efficiencyHist->GetXaxis()->FindBin(pt);
  double binContent = efficiencyHist->GetBinContent(ptBin);

  if (binContent > 1.e-5) {
    return 1.0 / binContent;
  } else return -1;
}


//---------------------------------------------------

double AliAnalysisTaskPIDCME::GetNUACor(int charge, double phi, double eta, double vz)
{
  double weightNUA = 1;
  if (fVzBin<0 || fCentBin<0 || fRunNum<0) return -1;
  if (fPeriod.EqualTo("LHC10h")) {
    hNUAweightPlus = (TH2D*)fListNUA->FindObject(Form("weightdPhidEta_run%i_cent0_vz%i_plus",fRunNum, fVzBin));
    hNUAweightMinus = (TH2D*)fListNUA->FindObject(Form("weightdPhidEta_run%i_cent0_vz%i_minus",fRunNum, fVzBin));
    if (!hNUAweightPlus || !hNUAweightMinus) return -1;
    if (charge>0) {
      int phiBin = hNUAweightPlus->GetXaxis()->FindBin(phi);
      int etaBin = hNUAweightPlus->GetYaxis()->FindBin(eta);
      if (hNUAweightPlus->GetBinContent(phiBin, etaBin)>0) weightNUA = hNUAweightPlus->GetBinContent(phiBin, etaBin);
      return weightNUA;
    } else if (charge<0) {
      int phiBin = hNUAweightMinus->GetXaxis()->FindBin(phi);
      int etaBin = hNUAweightMinus->GetYaxis()->FindBin(eta);
      if (hNUAweightMinus->GetBinContent(phiBin, etaBin)>0) weightNUA = hNUAweightMinus->GetBinContent(phiBin, etaBin);
      return weightNUA;
    }
  }
  if (fPeriod.EqualTo("LHC15o")|| fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) { // Rihan and Protty 's NUA Results
    if (charge>0) {
      if (!hCorrectNUAPos) return -1;
      int iBinNUA = hCorrectNUAPos->FindBin(vz,phi,eta);
      if (hCorrectNUAPos->GetBinContent(iBinNUA)>0) weightNUA = (double)hCorrectNUAPos->GetBinContent(iBinNUA);
      return  weightNUA;
    } else if (charge<0) {
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

inline double AliAnalysisTaskPIDCME::GetEventPlane(double qx, double qy, double harmonic)
{
  double psi = (1./harmonic)*TMath::ATan2(qy,qx);
  if (psi < 0) return psi += TMath::TwoPi()/harmonic;
  else return psi;
}

//---------------------------------------------------

inline double AliAnalysisTaskPIDCME::RangeDPhi(double dphi)
{
  while (dphi >= 1.5*TMath::Pi()) dphi -= TMath::TwoPi();
  while (dphi < -0.5*TMath::Pi()) dphi += TMath::TwoPi();
  return dphi;
}

//---------------------------------------------------

//[0]: TPC Pos [1] TPC Neg [2] TPC
std::array<double,2> AliAnalysisTaskPIDCME::GetOneSideTPCPlaneNoAutoCorr(std::vector<int> vec_id)
{
  double tempSumQ2xPos  = fSumQ2xTPCPos;
  double tempSumQ2yPos  = fSumQ2yTPCPos;
  double tempWgtMultPos = fWgtMultTPCPos;

  double tempSumQ2xNeg  = fSumQ2xTPCNeg;
  double tempSumQ2yNeg  = fSumQ2yTPCNeg;
  double tempWgtMultNeg = fWgtMultTPCNeg;

  double repeQ2xPos = 0., repeQ2yPos = 0., repeWgtMultPos = 0.;
  double repeQ2xNeg = 0., repeQ2yNeg = 0., repeWgtMultNeg = 0.;

  std::vector<std::unordered_map<int, std::vector<double>> ::iterator> vec_it_pos;
  std::vector<std::unordered_map<int, std::vector<double>> ::iterator> vec_it_neg;

  for (auto id : vec_id) {
    vec_it_pos.push_back(mapTPCPosTrksIDPhiWgt.find(id));
    vec_it_neg.push_back(mapTPCNegTrksIDPhiWgt.find(id));
  }

  for (auto it : vec_it_pos) {
    if(it != mapTPCPosTrksIDPhiWgt.end()) {
      repeQ2xPos     += it->second[1] * TMath::Cos(2 * (it->second)[0]);
      repeQ2yPos     += it->second[1] * TMath::Sin(2 * (it->second)[0]);
      repeWgtMultPos += it->second[1];
    }
  }

  for (auto it : vec_it_neg) {
    if(it != mapTPCNegTrksIDPhiWgt.end()) {
      repeQ2xNeg     += it->second[1] * TMath::Cos(2 * (it->second)[0]);
      repeQ2yNeg     += it->second[1] * TMath::Sin(2 * (it->second)[0]);
      repeWgtMultNeg += it->second[1];
    }
  }

  tempSumQ2xPos  -= repeQ2xPos;
  tempSumQ2yPos  -= repeQ2yPos;
  tempWgtMultPos -= repeWgtMultPos;

  tempSumQ2xNeg  -= repeQ2xNeg;
  tempSumQ2yNeg  -= repeQ2yNeg;
  tempWgtMultNeg -= repeWgtMultNeg;

  double psiNoAutoPos = nan(""), psiNoAutoNeg = nan("");
  if (tempWgtMultPos > 1.e-6) {
    tempSumQ2xPos /= tempWgtMultPos;
    tempSumQ2yPos /= tempWgtMultPos;
    psiNoAutoPos = GetEventPlane(tempSumQ2xPos,tempSumQ2yPos,2.);
  } else psiNoAutoPos = nan("");

  if (tempWgtMultNeg > 1.e-6) {
    tempSumQ2xNeg /= tempWgtMultNeg;
    tempSumQ2yNeg /= tempWgtMultNeg;
    psiNoAutoNeg = GetEventPlane(tempSumQ2xNeg,tempSumQ2yNeg,2.);
  } else psiNoAutoNeg = nan("");

  return {psiNoAutoPos, psiNoAutoNeg};
}

//---------------------------------------------------
double AliAnalysisTaskPIDCME::GetTPCPlaneNoAutoCorr(std::vector<int> vec_id) {
  double psiNoAuto = nan("");

  double tempSumQ2x = fSumQ2xTPC;
  double tempSumQ2y = fSumQ2yTPC;
  double tempWgtMult = fWgtMultTPC;
  double repeQ2x = 0., repeQ2y = 0., repeWgtMult = 0.;

  std::vector<std::unordered_map<int, std::vector<double>> ::iterator> vec_it;
  for (int id : vec_id) {
    vec_it.push_back(mapTPCTrksIDPhiWgt.find(id));
  }

  for (auto it : vec_it) {
    if(it != mapTPCTrksIDPhiWgt.end()) {
      repeQ2x     += it->second[1] * TMath::Cos(2 * (it->second)[0]);
      repeQ2y     += it->second[1] * TMath::Sin(2 * (it->second)[0]);
      repeWgtMult += it->second[1];
    }
  }

  tempSumQ2x  -= repeQ2x;
  tempSumQ2y  -= repeQ2y;
  tempWgtMult -= repeWgtMult;

  if (tempWgtMult > 1.e-6) {
    tempSumQ2x /= tempWgtMult;
    tempSumQ2y /= tempWgtMult;
    psiNoAuto = GetEventPlane(tempSumQ2x,tempSumQ2y,2.);
  } else psiNoAuto = nan("");

  return psiNoAuto;
}

//---------------------------------------------------

bool AliAnalysisTaskPIDCME::GetDCA(double &dcaxy, double &dcaz, AliAODTrack* track) {
  if(!track) return false;
  double r[3];
  if (track->GetXYZ(r)) {
    dcaxy = r[0];
    dcaz  = r[1];
  } else {
    double dcax = r[0] - fVertex[0];
    double dcay = r[1] - fVertex[1];
    dcaz  = r[2] - fVertex[2];
    dcaxy = sqrt(dcax*dcax + dcay*dcay);
  }
  return true;
}
//---------------------------------------------------

double AliAnalysisTaskPIDCME::GetQ2Bin(double q2) {
  double q2Bin = -1;
  double quantile = -999;
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r") ){
    TH1D* hQnPercentile_centThisEvt = (TH1D*)hQnPercentile->ProjectionY("qnPercentileCentiEvt", (int)fCent+1, (int)fCent+1, "e");
    TSpline3* sp = new TSpline3(hQnPercentile_centThisEvt);
    quantile = sp->Eval(q2);
    delete sp;
    sp = nullptr;
    delete hQnPercentile_centThisEvt;
    hQnPercentile_centThisEvt = nullptr;
  } else if (fPeriod.EqualTo("LHC15o")){
    quantile = splQ2c[(int)fCent]->Eval(q2);
  }

  double bin_size = 1.0/fNQ2Bins;
  for(int i = 0; i < fNQ2Bins; i++) {
    if(quantile >= i*bin_size && quantile < (i+1)*bin_size) {
        q2Bin = i + 0.5;
    }
  }
  return q2Bin;
}
//---------------------------------------------------
int AliAnalysisTaskPIDCME::DetermineTheBestTPCPID(AliAODTrack* track, int PIDPoints)
{
  float nSigmaTPCPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  float nSigmaTPCKaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  float nSigmaTPCProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  float nSigmaTPCPIDarray[3] = {nSigmaTPCPion, nSigmaTPCKaon, nSigmaTPCProton};
  if(PIDPoints == 1){
    return (nSigmaTPCPIDarray[0] < nSigmaTPCPIDarray[1]) ? 1 : 2;
  }
  if(PIDPoints == 2){
    return (nSigmaTPCPIDarray[0] < nSigmaTPCPIDarray[2]) ? 1 : 3;
  }
  if(PIDPoints == 3){
    return (nSigmaTPCPIDarray[1] < nSigmaTPCPIDarray[2]) ? 2 : 3;
  }
  if(PIDPoints == 4){
      if(nSigmaTPCPIDarray[0] < nSigmaTPCPIDarray[1] && nSigmaTPCPIDarray[0] < nSigmaTPCPIDarray[2]) return 1;
      else return (nSigmaTPCPIDarray[1] < nSigmaTPCPIDarray[2]) ? 2 : 3;
  }
  return -1;
}
