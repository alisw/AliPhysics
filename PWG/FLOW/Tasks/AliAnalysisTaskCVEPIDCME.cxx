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

//--------------------------------------------------------------------------------
// CVE and PID CME analysis
// Contributor: Chunzheng Wang, <chunzheng.wang@cern.ch>, Shanghai
//--------------------------------------------------------------------------------

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
#include "AliAnalysisTaskCVEPIDCME.h"

ClassImp(AliAnalysisTaskCVEPIDCME);

//---------------------------------------------------
AliAnalysisTaskCVEPIDCME::AliAnalysisTaskCVEPIDCME() :
  AliAnalysisTaskSE(),
  isUseTPCPlane(true),
  isUseVZEROPlane(true),
  isUseZDCPlane(false),
  isDoNUE(true),
  isDoNUA(false),
  isV0DaughterUseTOF(false),
  isQATPC(true),
  isQAVZERO(true),
  isQAZDC(true),
  isCalculatePIDFlow(true),
  isCalculateDeltaGamma(true),
  isCalculateDiffResult(true),
  isCalculateDeltaPhiSumPhi(false),
  isCalculateLambdaLambda(false),
  isCalculateProtonProton(false),
  isCalculateHadronHadron(false),
  isNarrowDcaCuts768(false),
  isStrictestProtonCut(false),
  isCheckDaughterProtonPassAllCuts(false),
  isUsePionRejection(true),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fPlanePtMin(0.2),
  fPlanePtMax(2.0),
  fEtaGapPos(0.1),
  fEtaGapNeg(-0.1),
  fFilterBit(768),
  fNclsCut(70),
  fChi2Max(4.0),
  fChi2Min(0.1),
  fDcaCutXY(2.4),
  fDcaCutZ(3.2),
  fPtMin(0.2),
  fPtMax(5.),
  fEtaCut(0.8),
  fDedxCut(10.),
  fProtonPtMin(0.2),
  fProtonPtMax(5.0),
  fAntiProtonPtMin(0.2),
  fAntiProtonPtMax(5.0),
  fNSigmaTPCCut(3.0),
  fNSigmaTOFCut(3.0),
  fV0PtMin(0.5),
  fV0CPAMin(0.997),
  fV0RapidityMax(0.5),
  fV0DecayLengthMin(3.0),
  fV0DecayLengthMax(100.),
  fV0DCAToPrimVtxMax(1.5),
  fV0DcaBetweenDaughtersMax(0.5),
  fDaughtersPtMax(20.0),
  fDaughtersEtaMax(0.8),
  fDaughtersTPCNclsMin(70),
  fDaughtersDCAToPrimVtxMin(0.05),
  fV0PosProtonTPCNsigma(3.0),
  fV0NegPionTPCNsigma(3.0),
  fV0NegProtonTPCNsigma(3.0),
  fV0PosPionTPCNsigma(3.0),
  fV0PosProtonTOFNsigma(4.0),
  fV0NegPionTOFNsigma(4.0),
  fV0NegProtonTOFNsigma(4.0),
  fV0PosPionTOFNsigma(4.0),
  fLambdaMassMean(1.115683),
  fLambdaMassRightCut(0.005),
  fLambdaMassLeftCut(0.005),
  fAntiLambdaMassRightCut(0.005),
  fAntiLambdaMassLeftCut(0.005),
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
  isRightTPCPlane(false),
  isRightVZEROPlane(false),
  isRightZDCPlane(false),
  mapTPCPosTrksIDPhiWgt(0),
  mapTPCNegTrksIDPhiWgt(0),
  vecParticle(0),
  vecParticleV0(0),
  vecParticleFromDecay(0),
  vecParticleFromDecayPassAllCuts(0),
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
  fHistProtonPt(nullptr),
  fHistProtonEta(nullptr),
  fHistProtonPhi(nullptr),
  fHistProtonDcaXY(nullptr),
  fHistProtonDcaZ(nullptr),
  fHist2ProtonSigTPC(nullptr),
  fHist2ProtonSigTOF(nullptr),
  fHistAntiProtonPt(nullptr),
  fHistAntiProtonEta(nullptr),
  fHistAntiProtonPhi(nullptr),
  fHistAntiProtonDcaXY(nullptr),
  fHistAntiProtonDcaZ(nullptr),
  fHist2AntiProtonSigTPC(nullptr),
  fHist2AntiProtonSigTOF(nullptr),
  fHistV0Pt(nullptr),
  fHistV0Eta(nullptr),
  fHistV0DcatoPrimVertex(nullptr),
  fHistV0CPA(nullptr),
  fHistV0DecayLength(nullptr),
  fResultsList(nullptr),
  fHist2Psi2TPCPosCent(nullptr),
  fHist2Psi2TPCNegCent(nullptr),
  fHist2Psi2V0CCent(nullptr),
  fHist2Psi2V0ACent(nullptr),
  fHist2Psi1ZNCCent(nullptr),
  fHist2Psi1ZNACent(nullptr),
  fProfileTPCPsi2Correlation(nullptr),
  fProfileV0MPsi2Correlation(nullptr),
  fProfileZDCPsi1Correlation(nullptr),
  fProfileZDCPsi2Correlation(nullptr),
  fProfileV0CTPCPosPsi2Correlation(nullptr),
  fProfileV0ATPCPosPsi2Correlation(nullptr),
  fProfileV0CTPCNegPsi2Correlation(nullptr),
  fProfileV0ATPCNegPsi2Correlation(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;
  for (int i = 0; i < 3; i++) pV0XMeanRead[i] = nullptr;
  for (int i = 0; i < 3; i++) pV0YMeanRead[i] = nullptr;
  for (int i = 0; i < 2; i++) hQx2mV0[i] = nullptr;
  for (int i = 0; i < 2; i++) hQy2mV0[i] = nullptr;
  for (int i = 0; i < 3; i++) vtxQuant1[i] = -999;
  for (int i = 0; i < 3; i++) vtxQuant2[i] = -999;
  for (int i = 0; i < 2; i++) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistVz[i] = nullptr;
  for (int i = 0; i < 8; i++) fHist2CentQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2MultCentQA[i] = nullptr;
  for (int i = 0; i < 6; i++) fHist2MultMultQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2EtaPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQxCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQyCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQxVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQyVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fHist2CalibPsi2V0CCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQxCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQyCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQxVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQyVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fHist2CalibPsi2V0ACent[i] = nullptr;
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
  for (int i = 0; i < 2; i++) fHistLambdaPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaCPA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaDecayLength[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaMass[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2LambdaMassPtY[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaCPA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaDecayLength[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaMass[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2AntiLambdaMassPtY[i] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentHadron[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentLambda[i][j] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaProton[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaHadron[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaLambda[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaProtonProton[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaHadronHadron[i] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaLambdaProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaLambdaHadron[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaLambdaLambda[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaProtonProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaHadronHadron[i][j] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaProtonDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaHadronDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaLambdaDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaProtonProtonDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaHadronHadronDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaProtonSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaHadronSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaLambdaSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaProtonProtonSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaHadronHadronSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaProtonDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaHadronDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaLambdaDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaProtonProtonDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaHadronHadronDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaProtonMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaHadronMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaLambdaMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaProtonDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaHadronDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaLambdaDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaProtonProtonDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaHadronHadronDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaProtonSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaHadronSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaLambdaSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaProtonProtonSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaHadronHadronSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaProtonDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaHadronDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaLambdaDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaProtonProtonDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaHadronHadronDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaProtonMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaHadronMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaLambdaMass[i] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiLambdaProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiLambdaHadron[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiLambdaLambda[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiProtonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiHadronHadron[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiLambdaProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiLambdaHadron[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiLambdaLambda[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiProtonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiHadronHadron[i][j] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaProtonDecay[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaPionDecay[i]   = nullptr;
  for (int i = 0; i < 4; i++) fProfileGammaLambdaProtonDecay[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileGammaLambdaPionDecay[i]   = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaProtonDecayPassAllCuts[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaPionDecayPassAllCuts[i]   = nullptr;
  for (int i = 0; i < 4; i++) fProfileGammaLambdaProtonDecayPassAllCuts[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileGammaLambdaPionDecayPassAllCuts[i]   = nullptr;
}

//---------------------------------------------------
AliAnalysisTaskCVEPIDCME::AliAnalysisTaskCVEPIDCME(const char *name) :
  AliAnalysisTaskSE(name),
  isUseTPCPlane(true),
  isUseVZEROPlane(true),
  isUseZDCPlane(false),
  isDoNUE(true),
  isDoNUA(false),
  isV0DaughterUseTOF(false),
  isQATPC(true),
  isQAVZERO(true),
  isQAZDC(true),
  isCalculatePIDFlow(true),
  isCalculateDeltaGamma(true),
  isCalculateDiffResult(true),
  isCalculateDeltaPhiSumPhi(false),
  isCalculateLambdaLambda(false),
  isCalculateProtonProton(false),
  isCalculateHadronHadron(false),
  isNarrowDcaCuts768(false),
  isStrictestProtonCut(false),
  isCheckDaughterProtonPassAllCuts(false),
  isUsePionRejection(true),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fPlanePtMin(0.2),
  fPlanePtMax(2.0),
  fEtaGapPos(0.1),
  fEtaGapNeg(-0.1),
  fFilterBit(768),
  fNclsCut(70),
  fChi2Max(4.0),
  fChi2Min(0.1),
  fDcaCutXY(2.4),
  fDcaCutZ(3.2),
  fPtMin(0.2),
  fPtMax(5.),
  fEtaCut(0.8),
  fDedxCut(10.),
  fProtonPtMin(0.2),
  fProtonPtMax(5.0),
  fAntiProtonPtMin(0.2),
  fAntiProtonPtMax(5.0),
  fNSigmaTPCCut(3.0),
  fNSigmaTOFCut(3.0),
  fV0PtMin(0.5),
  fV0CPAMin(0.997),
  fV0RapidityMax(0.5),
  fV0DecayLengthMin(3.0),
  fV0DecayLengthMax(100.),
  fV0DCAToPrimVtxMax(1.5),
  fV0DcaBetweenDaughtersMax(0.5),
  fDaughtersPtMax(20.0),
  fDaughtersEtaMax(0.8),
  fDaughtersTPCNclsMin(70),
  fDaughtersDCAToPrimVtxMin(0.05),
  fV0PosProtonTPCNsigma(3.0),
  fV0NegPionTPCNsigma(3.0),
  fV0NegProtonTPCNsigma(3.0),
  fV0PosPionTPCNsigma(3.0),
  fV0PosProtonTOFNsigma(4.0),
  fV0NegPionTOFNsigma(4.0),
  fV0NegProtonTOFNsigma(4.0),
  fV0PosPionTOFNsigma(4.0),
  fLambdaMassMean(1.115683),
  fLambdaMassRightCut(0.005),
  fLambdaMassLeftCut(0.005),
  fAntiLambdaMassRightCut(0.005),
  fAntiLambdaMassLeftCut(0.005),
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
  isRightTPCPlane(false),
  isRightVZEROPlane(false),
  isRightZDCPlane(false),
  mapTPCPosTrksIDPhiWgt(0),
  mapTPCNegTrksIDPhiWgt(0),
  vecParticle(0),
  vecParticleV0(0),
  vecParticleFromDecay(0),
  vecParticleFromDecayPassAllCuts(0),
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
  fHistProtonPt(nullptr),
  fHistProtonEta(nullptr),
  fHistProtonPhi(nullptr),
  fHistProtonDcaXY(nullptr),
  fHistProtonDcaZ(nullptr),
  fHist2ProtonSigTPC(nullptr),
  fHist2ProtonSigTOF(nullptr),
  fHistAntiProtonPt(nullptr),
  fHistAntiProtonEta(nullptr),
  fHistAntiProtonPhi(nullptr),
  fHistAntiProtonDcaXY(nullptr),
  fHistAntiProtonDcaZ(nullptr),
  fHist2AntiProtonSigTPC(nullptr),
  fHist2AntiProtonSigTOF(nullptr),
  fHistV0Pt(nullptr),
  fHistV0Eta(nullptr),
  fHistV0DcatoPrimVertex(nullptr),
  fHistV0CPA(nullptr),
  fHistV0DecayLength(nullptr),
  fResultsList(nullptr),
  fHist2Psi2TPCPosCent(nullptr),
  fHist2Psi2TPCNegCent(nullptr),
  fHist2Psi2V0CCent(nullptr),
  fHist2Psi2V0ACent(nullptr),
  fHist2Psi1ZNCCent(nullptr),
  fHist2Psi1ZNACent(nullptr),
  fProfileTPCPsi2Correlation(nullptr),
  fProfileV0MPsi2Correlation(nullptr),
  fProfileZDCPsi1Correlation(nullptr),
  fProfileZDCPsi2Correlation(nullptr),
  fProfileV0CTPCPosPsi2Correlation(nullptr),
  fProfileV0ATPCPosPsi2Correlation(nullptr),
  fProfileV0CTPCNegPsi2Correlation(nullptr),
  fProfileV0ATPCNegPsi2Correlation(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;
  for (int i = 0; i < 3; i++) pV0XMeanRead[i] = nullptr;
  for (int i = 0; i < 3; i++) pV0YMeanRead[i] = nullptr;
  for (int i = 0; i < 2; i++) hQx2mV0[i] = nullptr;
  for (int i = 0; i < 2; i++) hQy2mV0[i] = nullptr;
  for (int i = 0; i < 3; i++) vtxQuant1[i] = -999;
  for (int i = 0; i < 3; i++) vtxQuant2[i] = -999;
  for (int i = 0; i < 2; i++) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistVz[i] = nullptr;
  for (int i = 0; i < 8; i++) fHist2CentQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2MultCentQA[i] = nullptr;
  for (int i = 0; i < 6; i++) fHist2MultMultQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2EtaPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQxCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQyCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQxVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0CQyVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fHist2CalibPsi2V0CCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQxCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQyCent[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQxVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fProfileV0AQyVtx[i]= nullptr;
  for (int i = 0; i < 2; i++) fHist2CalibPsi2V0ACent[i] = nullptr;
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
  for (int i = 0; i < 2; i++) fHistLambdaPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaCPA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaDecayLength[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaMass[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2LambdaMassPtY[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaCPA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaDecayLength[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaMass[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2AntiLambdaMassPtY[i] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentHadron[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 2; j++) fProfile2RawFlowPtCentLambda[i][j] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaProton[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaHadron[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaLambda[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaProtonProton[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaHadronHadron[i] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaLambdaProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaLambdaHadron[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaLambdaLambda[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaProtonProton[i][j] = nullptr;
  for (int i = 0; i < 5; i++) for (int j = 0; j < 4; j++) fProfileGammaHadronHadron[i][j] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaProtonDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaHadronDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaLambdaDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaProtonProtonDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaHadronHadronDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaProtonSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaHadronSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaLambdaSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaProtonProtonSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaHadronHadronSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaProtonDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaHadronDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaLambdaDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaProtonProtonDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaHadronHadronDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaProtonMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaHadronMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffDeltaLambdaLambdaMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaProtonDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaHadronDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaLambdaDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaProtonProtonDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaHadronHadronDPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaProtonSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaHadronSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaLambdaSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaProtonProtonSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaHadronHadronSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaProtonDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaHadronDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaLambdaDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaProtonProtonDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaHadronHadronDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaProtonMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaHadronMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DiffGammaLambdaLambdaMass[i] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiLambdaProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiLambdaHadron[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiLambdaLambda[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiProtonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaDPhiHadronHadron[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiLambdaProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiLambdaHadron[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiLambdaLambda[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiProtonProton[i][j] = nullptr;
  for (int i = 0; i < 8; i++) for (int j = 0; j < 4; j++) fHist2DEtaSPhiHadronHadron[i][j] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaProtonDecay[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaPionDecay[i]   = nullptr;
  for (int i = 0; i < 4; i++) fProfileGammaLambdaProtonDecay[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileGammaLambdaPionDecay[i]   = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaProtonDecayPassAllCuts[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileDeltaLambdaPionDecayPassAllCuts[i]   = nullptr;
  for (int i = 0; i < 4; i++) fProfileGammaLambdaProtonDecayPassAllCuts[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfileGammaLambdaPionDecayPassAllCuts[i]   = nullptr;
  
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
}

//------------------------------------------------

AliAnalysisTaskCVEPIDCME::~AliAnalysisTaskCVEPIDCME()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  if (fQAList) delete fQAList;
  if (fResultsList) delete fResultsList;
}

//---------------------------------------------------

void AliAnalysisTaskCVEPIDCME::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate");
}

//---------------------------------------------------

void AliAnalysisTaskCVEPIDCME::UserCreateOutputObjects()
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
    if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      hNUEweightPlus  = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgPos");
      hNUEweightMinus = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgNeg");
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
    }
    if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      fHCorrectV0ChWeghts = new TH2F();
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
  fEvtCount = new TH1D("EvtCount", "Event Count", 23, 1, 24);
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
  fEvtCount->GetXaxis()->SetBinLabel(20,"Loop V0");
  fEvtCount->GetXaxis()->SetBinLabel(21,"Pair V0Trk");
  fEvtCount->GetXaxis()->SetBinLabel(22,"Pair V0V0");
  fEvtCount->GetXaxis()->SetBinLabel(23,"Pair TrkTrk");

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
  fHist2EtaPhi[0] = new TH2D("fHistEtaPhi", ";#eta;#phi", 16,-0.8,0.8, 100, 0, TMath::TwoPi());
  fHist2EtaPhi[1] = new TH2D("fHistEtaPhi_afterfNUA", ";#eta;#phi", 16,-0.8,0.8, 100, 0, TMath::TwoPi());
  fQAList->Add(fHistPt);
  fQAList->Add(fHistEta);
  fQAList->Add(fHistNhits);
  fQAList->Add(fHist2PDedx);
  fQAList->Add(fHistDcaXY);
  fQAList->Add(fHistDcaZ);
  for (int i = 0; i < 2; i++) fQAList->Add(fHistPhi[i]);
  for (int i = 0; i < 2; i++) fQAList->Add(fHist2EtaPhi[i]);

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
      fHist2CalibPsi2V0CCent[i] = new TH2D(Form("fHist2CalibPsi2V0CCent%s",charCalibStep.data()), "", 16, 0, 80, 50, 0, TMath::TwoPi());
      fQAList->Add(fProfileV0CQxCent[i]);
      fQAList->Add(fProfileV0CQyCent[i]);
      fQAList->Add(fProfileV0CQxVtx[i]);
      fQAList->Add(fProfileV0CQyVtx[i]);
      fQAList->Add(fHist2CalibPsi2V0CCent[i]);

      fProfileV0AQxCent[i] = new TProfile(Form("fProfileV0AQxCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileV0AQyCent[i] = new TProfile(Form("fProfileV0AQyCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileV0AQxVtx[i]  = new TProfile(Form("fProfileV0AQxVz%s",charCalibStep.data()), "", 20, -10, 10);
      fProfileV0AQyVtx[i]  = new TProfile(Form("fProfileV0AQyVz%s",charCalibStep.data()), "", 20, -10, 10);
      fHist2CalibPsi2V0ACent[i] = new TH2D(Form("fHist2CalibPsi2V0ACent%s",charCalibStep.data()), "", 16, 0, 80, 50, 0, TMath::Pi());
      fQAList->Add(fProfileV0AQxCent[i]);
      fQAList->Add(fProfileV0AQyCent[i]);
      fQAList->Add(fProfileV0AQxVtx[i]);
      fQAList->Add(fProfileV0AQyVtx[i]);
      fQAList->Add(fHist2CalibPsi2V0ACent[i]);
    }

    if (isQAZDC) {
      fProfileZNCQxCent[i] = new TProfile(Form("fProfileZNCQxCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileZNCQyCent[i] = new TProfile(Form("fProfileZNCQyCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fHist2CalibPsi1ZNCCent[i] = new TH2D(Form("fHist2CalibPsi1ZNCCent%s",charCalibStep.data()), "", 16, 0, 80, 50, 0, TMath::Pi());
      fQAList->Add(fProfileZNCQxCent[i]);
      fQAList->Add(fProfileZNCQyCent[i]);
      fQAList->Add(fHist2CalibPsi1ZNCCent[i]);

      fProfileZNAQxCent[i] = new TProfile(Form("fProfileZNAQxCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fProfileZNAQyCent[i] = new TProfile(Form("fProfileZNAQyCent%s",charCalibStep.data()), "", 80, 0, 80.);
      fHist2CalibPsi1ZNACent[i] = new TH2D(Form("fHist2CalibPsi1ZNACent%s",charCalibStep.data()), "", 16, 0, 80, 50, 0, TMath::Pi());
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

  //Proton QA
  fHistProtonPt      = new TH1D("fHistProtonPt"     , "fHistProtonPt;p_{T}"      , 200, 0., 20.);
  fHistProtonEta     = new TH1D("fHistProtonEta"    , "fHistProtonEta;#eta"      , 100, -2., 2.);
  fHistProtonPhi     = new TH1D("fHistProtonPhi"    , "fHistProtonPhi;#phi"      , 100, 0., TMath::TwoPi());
  fHistProtonDcaXY   = new TH1D("fHistProtonDcaXY"  , "fHistProtonDcaXY;DcaXY"   , 500, 0., 5);
  fHistProtonDcaZ    = new TH1D("fHistProtonDcaZ "  , "fHistProtonDcaZ;DcaZ"     , 500, 0., 5);
  fHist2ProtonSigTPC = new TH2D("fHist2ProtonSigTPC", "fHist2ProtonSigTPC;SigTPC", 25, 0., 5., 100, 0.,10.);
  fHist2ProtonSigTOF = new TH2D("fHist2ProtonSigTOF", "fHist2ProtonSigTOF;SigTOF", 25, 0., 5., 100, 0.,10.);
  fQAList->Add(fHistProtonPt);
  fQAList->Add(fHistProtonEta);
  fQAList->Add(fHistProtonPhi);
  fQAList->Add(fHistProtonDcaXY);
  fQAList->Add(fHistProtonDcaZ);
  fQAList->Add(fHist2ProtonSigTPC);
  fQAList->Add(fHist2ProtonSigTOF);
  fHistAntiProtonPt      = new TH1D("fHistAntiProtonPt"     , "fHistAntiProtonPt;p_{T}"      , 200, 0., 20.); 
  fHistAntiProtonEta     = new TH1D("fHistAntiProtonEta"    , "fHistAntiProtonEta;#eta"      , 100, -2., 2.); 
  fHistAntiProtonPhi     = new TH1D("fHistAntiProtonPhi"    , "fHistAntiProtonPhi;#phi"      , 100, 0., TMath::TwoPi()); 
  fHistAntiProtonDcaXY   = new TH1D("fHistAntiProtonDcaXY"  , "fHistAntiProtonDcaXY;DcaXY"   , 500, 0., 5); 
  fHistAntiProtonDcaZ    = new TH1D("fHistAntiProtonDcaZ "  , "fHistAntiProtonDcaZ;DcaZ"     , 500, 0., 5); 
  fHist2AntiProtonSigTPC = new TH2D("fHist2AntiProtonSigTPC", "fHist2AntiProtonSigTPC;p_{T};SigTPC", 25, 0., 5., 100, 0.,10.); 
  fHist2AntiProtonSigTOF = new TH2D("fHist2AntiProtonSigTOF", "fHist2AntiProtonSigTOF;p_{T};SigTOF", 25, 0., 5., 100, 0.,10.); 
  fQAList->Add(fHistAntiProtonPt);
  fQAList->Add(fHistAntiProtonEta);
  fQAList->Add(fHistAntiProtonPhi);
  fQAList->Add(fHistAntiProtonDcaXY);
  fQAList->Add(fHistAntiProtonDcaZ);
  fQAList->Add(fHist2AntiProtonSigTPC);
  fQAList->Add(fHist2AntiProtonSigTOF);

  //V0s QA
  fHistV0Pt = new TH1D("hV0Pt","", 200, 0., 20.);
  fHistV0Eta = new TH1D("hV0Eta","", 200, -3., 3.);
  fHistV0DcatoPrimVertex = new TH1D("hV0DcaToPrimVertex","",200, 0., 20.);
  fHistV0CPA = new TH1D("hV0CPA","", 1000, 0.9, 1.);
  fHistV0DecayLength = new TH1D("hV0DecayLength","",500,0,500.);
  fQAList->Add(fHistV0Pt);
  fQAList->Add(fHistV0Eta);
  fQAList->Add(fHistV0DcatoPrimVertex);
  fQAList->Add(fHistV0CPA);
  fQAList->Add(fHistV0DecayLength);

  //Lambda QA
  std::string charLambdaMassCut;
  for (int i=0; i<2; i++) {
    ///// Case 0 = before cut, case 1 = afterCut.
    if (i==0) charLambdaMassCut = "Bf";
    if (i==1) charLambdaMassCut = "Af";
    /// Lambda:
    fHistLambdaPt[i] = new TH1D(Form("hLambdaPt_%sMassCut",charLambdaMassCut.data()),";pT", 200, 0., 20.);
    fHistLambdaEta[i] = new TH1D(Form("hLambdaEta_%sMassCut",charLambdaMassCut.data()),";#eta",200, -3., 3.);
    fHistLambdaPhi[i] = new TH1D(Form("hLambdaPhi_%sMassCut",charLambdaMassCut.data()),";phi", 360, 0., TMath::TwoPi());
    fHistLambdaDcaToPrimVertex[i] = new TH1D(Form("hLambdaDcaToPrimVertex_%sMassCut",charLambdaMassCut.data()),"DcatoPV",200, 0., 20.);
    fHistLambdaCPA[i] = new TH1D(Form("hLambdaCPA_%sMassCut",charLambdaMassCut.data()),";CPA",200, 0.9, 1.);
    fHistLambdaDecayLength[i] = new TH1D(Form("hLambdaDecayLength_%sMassCut",charLambdaMassCut.data()),";DecayLength", 250, 0., 500.);
    fHistLambdaMass[i] = new TH1D(Form("hLambdaMass_%sMassCut",charLambdaMassCut.data()),";Mass",1000,1.,1.25); //  Current bin size = 0.00025
    fHist2LambdaMassPtY[i] = new TH2D(Form("fHist2LambdaMassPtY_%sMassCut",charLambdaMassCut.data()),";pT;yH",12,0.,6.,20,-1.,1.);
    fQAList->Add(fHistLambdaPt[i]);
    fQAList->Add(fHistLambdaEta[i]);
    fQAList->Add(fHistLambdaPhi[i]);
    fQAList->Add(fHistLambdaDcaToPrimVertex[i]);
    fQAList->Add(fHistLambdaCPA[i]);
    fQAList->Add(fHistLambdaDecayLength[i]);
    fQAList->Add(fHistLambdaMass[i]);
    fQAList->Add(fHist2LambdaMassPtY[i]);

    // AntiLambda
    fHistAntiLambdaPt[i] = new TH1D(Form("hAntiLambdaPt_%sMassCut",charLambdaMassCut.data()),";pT", 200, 0., 20.);
    fHistAntiLambdaEta[i] = new TH1D(Form("hAntiLambdaEta_%sMassCut",charLambdaMassCut.data()),";#eta",200, -3., 3.);
    fHistAntiLambdaPhi[i] = new TH1D(Form("hAntiLambdaPhi_%sMassCut",charLambdaMassCut.data()),";phi", 360, 0., TMath::TwoPi());
    fHistAntiLambdaDcaToPrimVertex[i] = new TH1D(Form("hAntiLambdaDcaToPrimVertex_%sMassCut",charLambdaMassCut.data()),"DcatoPV",200, 0., 20.);
    fHistAntiLambdaCPA[i] = new TH1D(Form("hAntiLambdaCPA_%sMassCut",charLambdaMassCut.data()),";CPA",200, 0.9, 1.);
    fHistAntiLambdaDecayLength[i] = new TH1D(Form("hAntiLambdaDecayLength_%sMassCut",charLambdaMassCut.data()),";DecayLength", 250, 0., 500.);
    fHistAntiLambdaMass[i] = new TH1D(Form("hAntiLambdaMass_%sMassCut",charLambdaMassCut.data()),";Mass",1000,1.,1.25); //  Current bin size = 0.00025
    fHist2AntiLambdaMassPtY[i] = new TH2D(Form("fHist2AntiLambdaMassPtY_%sMassCut",charLambdaMassCut.data()),";pT;yH",12,0.,6.,20,-1.,1.);
    fQAList->Add(fHistAntiLambdaPt[i]);
    fQAList->Add(fHistAntiLambdaEta[i]);
    fQAList->Add(fHistAntiLambdaPhi[i]);
    fQAList->Add(fHistAntiLambdaDcaToPrimVertex[i]);
    fQAList->Add(fHistAntiLambdaCPA[i]);
    fQAList->Add(fHistAntiLambdaDecayLength[i]);
    fQAList->Add(fHistAntiLambdaMass[i]);
    fQAList->Add(fHist2AntiLambdaMassPtY[i]);
  }
  PostData(1,fQAList);
  if (fDebug) Printf("Post fQAList Data Success!");

  ////////////////////////
  // Results
  ////////////////////////
  fResultsList = new TList();
  fResultsList -> SetName("fResultsList");
  fResultsList -> SetOwner(kTRUE);

  // Plane
  fHist2Psi2TPCPosCent = new TH2D("fHist2Psi2TPCPosCent","fHist2Psi2TPCPosCent;centrality;#Psi(TPCPos)",8,0,80.,100,0.,TMath::Pi());
  fHist2Psi2TPCNegCent = new TH2D("fHist2Psi2TPCNegCent","fHist2Psi2TPCNegCent;centrality;#Psi(TPCNeg)",8,0,80.,100,0.,TMath::Pi());
  fHist2Psi2V0CCent    = new TH2D("fHist2Psi2V0CCent","fHist2Psi2V0CCent;centrality;#Psi(V0C)",8,0.,80.,100,0.,TMath::Pi());
  fHist2Psi2V0ACent    = new TH2D("fHist2Psi2V0ACent","fHist2Psi2V0ACent;centrality;#Psi(V0A)",8,0.,80.,100,0.,TMath::Pi());
  fHist2Psi1ZNCCent    = new TH2D("fHist2Psi1ZNCCent","fHist2Psi1ZNCCent;centrality;#Psi(ZNC)",8,0.,80.,100,0.,TMath::TwoPi());
  fHist2Psi1ZNACent    = new TH2D("fHist2Psi1ZNACent","fHist2Psi1ZNACent;centrality;#Psi(ZNA)",8,0.,80.,100,0.,TMath::TwoPi());
  
  fResultsList->Add(fHist2Psi2TPCPosCent);
  fResultsList->Add(fHist2Psi2TPCNegCent);
  fResultsList->Add(fHist2Psi2V0CCent);
  fResultsList->Add(fHist2Psi2V0ACent);
  fResultsList->Add(fHist2Psi1ZNCCent);
  fResultsList->Add(fHist2Psi1ZNACent);

  // Res
  fProfileTPCPsi2Correlation       = new TProfile("fProfileTPCPsi2Correlation","fProfileTPCPsi2Correlation;centrality;Res",8,0.,80.);
  fProfileV0MPsi2Correlation       = new TProfile("fProfileV0MPsi2Correlation","fProfileV0MPsi2Correlation;centrality;Res",8,0.,80.);
  fProfileZDCPsi1Correlation       = new TProfile("fProfileZDCPsi1Correlation","fProfileZDCPsi1Correlation;centrality;Res",8,0.,80.);
  fProfileZDCPsi2Correlation       = new TProfile("fProfileZDCPsi2Correlation","fProfileZDCPsi2Correlation;centrality;Res",8,0.,80.);
  fProfileV0CTPCPosPsi2Correlation = new TProfile("fProfileV0CTPCPosPsi2Correlation","fProfileV0CTPCPosPsi2Correlation;centrality;Res",8,0.,80.);
  fProfileV0ATPCPosPsi2Correlation = new TProfile("fProfileV0ATPCPosPsi2Correlation","fProfileV0ATPCPosPsi2Correlation;centrality;Res",8,0.,80.);
  fProfileV0CTPCNegPsi2Correlation = new TProfile("fProfileV0CTPCNegPsi2Correlation","fProfileV0CTPCNegPsi2Correlation;centrality;Res",8,0.,80.);
  fProfileV0ATPCNegPsi2Correlation = new TProfile("fProfileV0ATPCNegPsi2Correlation","fProfileV0ATPCNegPsi2Correlation;centrality;Res",8,0.,80.);
  fResultsList->Add(fProfileTPCPsi2Correlation);
  fResultsList->Add(fProfileV0MPsi2Correlation);
  fResultsList->Add(fProfileZDCPsi1Correlation);
  fResultsList->Add(fProfileZDCPsi2Correlation);
  fResultsList->Add(fProfileV0CTPCPosPsi2Correlation);
  fResultsList->Add(fProfileV0ATPCPosPsi2Correlation);
  fResultsList->Add(fProfileV0CTPCNegPsi2Correlation);
  fResultsList->Add(fProfileV0ATPCNegPsi2Correlation);

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
          fProfile2RawFlowPtCentHadron[iPlane][jCharge] = new TProfile2D(Form("fProfile2RawFlowPtCentHadron%s%s",charPlane.data(),charCharge.data()),";centrality;pT;v2",8,0.,80.,25,0.,5.);
          fProfile2RawFlowPtCentProton[iPlane][jCharge] = new TProfile2D(Form("fProfile2RawFlowPtCentProton%s%s",charPlane.data(),charCharge.data()),";centrality;pT;v2",8,0.,80.,25,0.,5.);
          fProfile2RawFlowPtCentLambda[iPlane][jCharge] = new TProfile2D(Form("fProfile2RawFlowPtCentLambda%s%s",charPlane.data(),charCharge.data()),";centrality;pT;v2",8,0.,80.,25,0.,5.);
          fResultsList->Add(fProfile2RawFlowPtCentHadron[iPlane][jCharge]);
          fResultsList->Add(fProfile2RawFlowPtCentProton[iPlane][jCharge]);
          fResultsList->Add(fProfile2RawFlowPtCentLambda[iPlane][jCharge]);
      }
    }
  }

  if (isCalculateDeltaGamma) {
    for (int i = 0; i < 4; i++) {
      fProfileDeltaLambdaProton[i] = new TProfile(Form("fProfileDeltaLambdaProton_type%i",i),";centrality;#delta",8,0.,80.);
      fProfileDeltaLambdaHadron[i] = new TProfile(Form("fProfileDeltaLambdaHadron_type%i",i),";centrality;#delta",8,0.,80.);
      fResultsList->Add(fProfileDeltaLambdaProton[i]);
      fResultsList->Add(fProfileDeltaLambdaHadron[i]);
      if (isCalculateLambdaLambda) {
        fProfileDeltaLambdaLambda[i] = new TProfile(Form("fProfileDeltaLambdaLambda_type%i",i),";centrality;#delta",8,0.,80.);
        fResultsList->Add(fProfileDeltaLambdaLambda[i]);
      }
      if (isCalculateProtonProton) { 
        fProfileDeltaProtonProton[i] = new TProfile(Form("fProfileDeltaProtonProton_type%i",i),";centrality;#delta",8,0.,80.);
        fResultsList->Add(fProfileDeltaProtonProton[i]);
      }
      if (isCalculateHadronHadron) {
        fProfileDeltaHadronHadron[i] = new TProfile(Form("fProfileDeltaHadronHadron_type%i",i),";centrality;#delta",8,0.,80.);
        fResultsList->Add(fProfileDeltaHadronHadron[i]);
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
      for (int iType = 0; iType < 4; iType++) {
        fProfileGammaLambdaProton[iPlane][iType] = new TProfile(Form("fProfileGammaLambdaProton%s_%i",charPlane.data(),iType),";centrality;#gamma",8,0.,80.);
        fProfileGammaLambdaHadron[iPlane][iType] = new TProfile(Form("fProfileGammaLambdaHadron%s_%i",charPlane.data(),iType),";centrality;#gamma",8,0.,80.);
        fResultsList->Add(fProfileGammaLambdaProton[iPlane][iType]);
        fResultsList->Add(fProfileGammaLambdaHadron[iPlane][iType]);

        if (isCalculateLambdaLambda) {
          fProfileGammaLambdaLambda[iPlane][iType] = new TProfile(Form("fProfileGammaLambdaLambda%s_%i",charPlane.data(),iType),";centrality;#gamma",8,0.,80.);
          fResultsList->Add(fProfileGammaLambdaLambda[iPlane][iType]);
        }
        if (isCalculateProtonProton) {
          fProfileGammaProtonProton[iPlane][iType] = new TProfile(Form("fProfileGammaProtonProton%s_%i",charPlane.data(),iType),";centrality;#gamma",8,0.,80.);
          fResultsList->Add(fProfileGammaProtonProton[iPlane][iType]);
        }
        if (isCalculateHadronHadron) {
          fProfileGammaHadronHadron[iPlane][iType] = new TProfile(Form("fProfileGammaHadronHadron%s_%i",charPlane.data(),iType),";centrality;#gamma",8,0.,80.);
          fResultsList->Add(fProfileGammaHadronHadron[iPlane][iType]);
        }
      }
    }
  }

  if (isCalculateDiffResult) {
    // δ
    for (int iType = 0; iType < 4; iType++) {
      fProfile2DiffDeltaLambdaProtonDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaLambdaProtonDPt_%i",iType),";centrality;#Delta pT;#delta",8,0.,80.,16,-8,8);
      fProfile2DiffDeltaLambdaProtonSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaLambdaProtonSPt_%i",iType),";centrality;#sum pT;#delta",8,0.,80.,16,0,16);
      fProfile2DiffDeltaLambdaProtonDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaLambdaProtonDEta_%i",iType),";centrality;#Delta Eta;#delta",8,0.,80.,16,-1.6,1.6);
      fProfile2DiffDeltaLambdaProtonMass[iType] = new TProfile2D(Form("fProfile2DiffDeltaLambdaProtonMass_%i",iType),";centrality;#Lambda Mass;#delta",8,0.,80,10,fLambdaMassMean-fLambdaMassRightCut,fLambdaMassMean+fLambdaMassLeftCut);
      fResultsList->Add(fProfile2DiffDeltaLambdaProtonDPt[iType]);
      fResultsList->Add(fProfile2DiffDeltaLambdaProtonSPt[iType]);
      fResultsList->Add(fProfile2DiffDeltaLambdaProtonDEta[iType]);
      fResultsList->Add(fProfile2DiffDeltaLambdaProtonMass[iType]);
      fProfile2DiffDeltaLambdaHadronDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaLambdaHadronDPt_%i",iType),";centrality;#Delta pT;#delta",8,0.,80.,16,-8,8);
      fProfile2DiffDeltaLambdaHadronSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaLambdaHadronSPt_%i",iType),";centrality;#sum pT;#delta",8,0.,80.,16,0,16);
      fProfile2DiffDeltaLambdaHadronDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaLambdaHadronDEta_%i",iType),";centrality;#Delta Eta;#delta",8,0.,80.,16,-1.6,1.6);
      fProfile2DiffDeltaLambdaHadronMass[iType] = new TProfile2D(Form("fProfile2DiffDeltaLambdaHadronMass_%i",iType),";centrality;#Lambda Mass;#delta",8,0.,80,10,fLambdaMassMean-fLambdaMassRightCut,fLambdaMassMean+fLambdaMassLeftCut);
      fResultsList->Add(fProfile2DiffDeltaLambdaHadronDPt[iType]);
      fResultsList->Add(fProfile2DiffDeltaLambdaHadronSPt[iType]);
      fResultsList->Add(fProfile2DiffDeltaLambdaHadronDEta[iType]);
      fResultsList->Add(fProfile2DiffDeltaLambdaHadronMass[iType]);

      if (isCalculateLambdaLambda) {
        fProfile2DiffDeltaLambdaLambdaDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaLambdaLambdaDPt_%i",iType),";centrality;#Delta pT;#delta",8,0.,80.,16,-8,8);
        fProfile2DiffDeltaLambdaLambdaSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaLambdaLambdaSPt_%i",iType),";centrality;#sum pT;#delta",8,0.,80.,16,0,16);
        fProfile2DiffDeltaLambdaLambdaDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaLambdaLambdaDEta_%i",iType),";centrality;#Delta Eta;#delta",8,0.,80.,16,-1.6,1.6);
        fProfile2DiffDeltaLambdaLambdaMass[iType] = new TProfile2D(Form("fProfile2DiffDeltaLambdaLambdaMass_%i",iType),";centrality;#Lambda Mass;#delta",8,0.,80,10,fLambdaMassMean-fLambdaMassRightCut,fLambdaMassMean+fLambdaMassLeftCut);
        fResultsList->Add(fProfile2DiffDeltaLambdaLambdaDPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaLambdaLambdaSPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaLambdaLambdaDEta[iType]);
        fResultsList->Add(fProfile2DiffDeltaLambdaLambdaMass[iType]);
      }
      if (isCalculateProtonProton) {
        fProfile2DiffDeltaProtonProtonDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaProtonProtonDPt_%i",iType),";centrality;#Delta pT;#delta",8,0.,80.,16,-8,8);
        fProfile2DiffDeltaProtonProtonSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaProtonProtonSPt_%i",iType),";centrality;#sum pT;#delta",8,0.,80.,16,0,16);
        fProfile2DiffDeltaProtonProtonDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaProtonProtonDEta_%i",iType),";centrality;#Delta Eta;#delta",8,0.,80.,16,-1.6,1.6);
        fResultsList->Add(fProfile2DiffDeltaProtonProtonDPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaProtonProtonSPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaProtonProtonDEta[iType]);
      }
      if (isCalculateHadronHadron) {
        fProfile2DiffDeltaHadronHadronDPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaHadronHadronDPt_%i",iType),";centrality;#Delta pT;#delta",8,0.,80.,16,-8,8);
        fProfile2DiffDeltaHadronHadronSPt[iType]  = new TProfile2D(Form("fProfile2DiffDeltaHadronHadronSPt_%i",iType),";centrality;#sum pT;#delta",8,0.,80.,16,0,16);
        fProfile2DiffDeltaHadronHadronDEta[iType] = new TProfile2D(Form("fProfile2DiffDeltaHadronHadronDEta_%i",iType),";centrality;#Delta Eta;#delta",8,0.,80.,16,-1.6,1.6);
        fResultsList->Add(fProfile2DiffDeltaHadronHadronDPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaHadronHadronSPt[iType]);
        fResultsList->Add(fProfile2DiffDeltaHadronHadronDEta[iType]);
      }

      // γ
      fProfile2DiffGammaLambdaProtonDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaLambdaProtonDPt_%i",iType),";centrality;#Delta pT;#gamma",8,0.,80.,16,-8,8);
      fProfile2DiffGammaLambdaProtonSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaLambdaProtonSPt_%i",iType),";centrality;#sum pT;#gamma",8,0.,80.,16,0,16);
      fProfile2DiffGammaLambdaProtonDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaLambdaProtonDEta_%i",iType),";centrality;#Delta #Eta;#gamma",8,0.,80.,16,-1.6,1.6);
      fProfile2DiffGammaLambdaProtonMass[iType] = new TProfile2D(Form("fProfile2DiffGammaLambdaProtonMass_%i",iType),";centrality;#Lambda Mass;#gamma",8,0.,80,10,fLambdaMassMean-fLambdaMassRightCut,fLambdaMassMean+fLambdaMassLeftCut);
      fResultsList->Add(fProfile2DiffGammaLambdaProtonDPt[iType]);
      fResultsList->Add(fProfile2DiffGammaLambdaProtonSPt[iType]);
      fResultsList->Add(fProfile2DiffGammaLambdaProtonDEta[iType]);
      fResultsList->Add(fProfile2DiffGammaLambdaProtonMass[iType]);
      fProfile2DiffGammaLambdaHadronDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaLambdaHadronDPt_%i",iType),";centrality;#Delta pT;#gamma",8,0.,80.,16,-8,8);
      fProfile2DiffGammaLambdaHadronSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaLambdaHadronSPt_%i",iType),";centrality;#sum pT;#gamma",8,0.,80.,16,0,16);
      fProfile2DiffGammaLambdaHadronDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaLambdaHadronDEta_%i",iType),";centrality;#Delta #Eta;#gamma",8,0.,80.,16,-1.6,1.6);
      fProfile2DiffGammaLambdaHadronMass[iType] = new TProfile2D(Form("fProfile2DiffGammaLambdaHadronMass_%i",iType),";centrality;#Lambda Mass;#gamma",8,0.,80,10,fLambdaMassMean-fLambdaMassRightCut,fLambdaMassMean+fLambdaMassLeftCut);
      fResultsList->Add(fProfile2DiffGammaLambdaHadronDPt[iType]);
      fResultsList->Add(fProfile2DiffGammaLambdaHadronSPt[iType]);
      fResultsList->Add(fProfile2DiffGammaLambdaHadronDEta[iType]);
      fResultsList->Add(fProfile2DiffGammaLambdaHadronMass[iType]);
      
      if (isCalculateLambdaLambda) {
        fProfile2DiffGammaLambdaLambdaDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaLambdaLambdaDPt_%i",iType),";centrality;#Delta pT;#gamma",8,0.,80.,16,-8,8);
        fProfile2DiffGammaLambdaLambdaSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaLambdaLambdaSPt_%i",iType),";centrality;#sum pT;#gamma",8,0.,80.,16,0,16);
        fProfile2DiffGammaLambdaLambdaDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaLambdaLambdaDEta_%i",iType),";centrality;#Delta #Eta;#gamma",8,0.,80.,16,-1.6,1.6);
        fProfile2DiffGammaLambdaLambdaMass[iType] = new TProfile2D(Form("fProfile2DiffGammaLambdaLambdaMass_%i",iType),";centrality;#Lambda Mass;#gamma",8,0.,80,10,fLambdaMassMean-fLambdaMassRightCut,fLambdaMassMean+fLambdaMassLeftCut);
        fResultsList->Add(fProfile2DiffGammaLambdaLambdaDPt[iType]);
        fResultsList->Add(fProfile2DiffGammaLambdaLambdaSPt[iType]);
        fResultsList->Add(fProfile2DiffGammaLambdaLambdaDEta[iType]);
        fResultsList->Add(fProfile2DiffGammaLambdaLambdaMass[iType]);
      }
      if (isCalculateProtonProton) {
        fProfile2DiffGammaProtonProtonDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaProtonProtonDPt_%i",iType),";centrality;#Delta pT;#gamma",8,0.,80.,16,-8,8);
        fProfile2DiffGammaProtonProtonSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaProtonProtonSPt_%i",iType),";centrality;#sum pT;#gamma",8,0.,80.,16,0,16);
        fProfile2DiffGammaProtonProtonDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaProtonProtonDEta_%i",iType),";centrality;#Delta #Eta;#gamma",8,0.,80.,16,-1.6,1.6);
        fResultsList->Add(fProfile2DiffGammaProtonProtonDPt[iType]);
        fResultsList->Add(fProfile2DiffGammaProtonProtonSPt[iType]);
        fResultsList->Add(fProfile2DiffGammaProtonProtonDEta[iType]);
      }
      if (isCalculateHadronHadron) {
        fProfile2DiffGammaHadronHadronDPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaHadronHadronDPt_%i",iType),";centrality;#Delta pT;#gamma",8,0.,80.,16,-8,8);
        fProfile2DiffGammaHadronHadronSPt[iType]  = new TProfile2D(Form("fProfile2DiffGammaHadronHadronSPt_%i",iType),";centrality;#sum pT;#gamma",8,0.,80.,16,0,16);
        fProfile2DiffGammaHadronHadronDEta[iType] = new TProfile2D(Form("fProfile2DiffGammaHadronHadronSEta_%i",iType),";centrality;#Delta #Eta;#gamma",8,0.,80.,16,-1.6,1.6);
        fResultsList->Add(fProfile2DiffGammaHadronHadronDPt[iType]);
        fResultsList->Add(fProfile2DiffGammaHadronHadronSPt[iType]);
        fResultsList->Add(fProfile2DiffGammaHadronHadronDEta[iType]);
      }
    }
  }
  
  if (isCalculateDeltaPhiSumPhi) {
    for (int iCent = 0; iCent < 8; iCent++) {
      for (int iType = 0; iType < 4; iType++) {
        // C(Δη,Δφ) [cent][pair type]
        fHist2DEtaDPhiLambdaProton[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiLambdaProton_cent%i_%i",iCent,iType),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        fResultsList->Add(fHist2DEtaDPhiLambdaProton[iCent][iType]);
        fHist2DEtaDPhiLambdaHadron[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiLambdaHadron_cent%i_%i",iCent,iType),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        fResultsList->Add(fHist2DEtaDPhiLambdaHadron[iCent][iType]);

        if (isCalculateLambdaLambda) {
          fHist2DEtaDPhiLambdaLambda[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiLambdaLambda_cent%i_%i",iCent,iType),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaDPhiLambdaLambda[iCent][iType]);
        }
        if (isCalculateProtonProton) {
          fHist2DEtaDPhiProtonProton[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiProtonProton_cent%i_%i",iCent,iType),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaDPhiProtonProton[iCent][iType]);
        }
        if (isCalculateHadronHadron) {
          fHist2DEtaDPhiHadronHadron[iCent][iType] = new TH2D(Form("fHist2DEtaDPhiHadronHadron_cent%i_%i",iCent,iType),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaDPhiHadronHadron[iCent][iType]);
        }

     // C(Δη,sφ) [cent][pair type] only TPC Plane
        fHist2DEtaSPhiLambdaProton[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiLambdaProton_cent%i_%i",iCent,iType),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        fResultsList->Add(fHist2DEtaSPhiLambdaProton[iCent][iType]);
        fHist2DEtaSPhiLambdaHadron[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiLambdaHadron_cent%i_%i",iCent,iType),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        fResultsList->Add(fHist2DEtaSPhiLambdaHadron[iCent][iType]);
        if (isCalculateLambdaLambda) {
          fHist2DEtaSPhiLambdaLambda[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiLambdaLambda_cent%i_%i",iCent,iType),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaSPhiLambdaLambda[iCent][iType]);
        }
        if (isCalculateProtonProton) {
          fHist2DEtaSPhiProtonProton[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiProtonProton_cent%i_%i",iCent,iType),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaSPhiProtonProton[iCent][iType]);
        }
        if (isCalculateHadronHadron) {
          fHist2DEtaSPhiHadronHadron[iCent][iType] = new TH2D(Form("fHist2DEtaSPhiHadronHadron_cent%i_%i",iCent,iType),";deta;dphi", 20, -1.6, 1.6, 36, -0.5*TMath::Pi(), 1.5*TMath::Pi());
          fResultsList->Add(fHist2DEtaSPhiHadronHadron[iCent][iType]);
        }
      }
    }
  }

  if (isCalculateDeltaGamma) {
    for (int iType = 0; iType < 4; iType++) {
      fProfileDeltaLambdaProtonDecay[iType] = new TProfile(Form("fProfileDeltaLambdaProtonDecay_%i",iType),";centrality;#delta", 8,0.,80);
      fProfileDeltaLambdaPionDecay[iType]   = new TProfile(Form("fProfileDeltaLambdaPionDecay_%i",iType)  ,";centrality;#delta", 8,0.,80);
      fProfileGammaLambdaProtonDecay[iType] = new TProfile(Form("fProfileGammaLambdaProtonDecay_%i",iType),";centrality;#gamma", 8,0.,80);
      fProfileGammaLambdaPionDecay[iType]   = new TProfile(Form("fProfileGammaLambdaPionDecay_%i",iType)  ,";centrality;#gamma", 8,0.,80);
      fResultsList->Add(fProfileDeltaLambdaProtonDecay[iType]);
      fResultsList->Add(fProfileDeltaLambdaPionDecay[iType]);
      fResultsList->Add(fProfileGammaLambdaProtonDecay[iType]);
      fResultsList->Add(fProfileGammaLambdaPionDecay[iType]);
    }
    if (isCheckDaughterProtonPassAllCuts) {
      for (int iType = 0; iType < 4; iType++) {
        fProfileDeltaLambdaProtonDecayPassAllCuts[iType] = new TProfile(Form("fProfileDeltaLambdaProtonDecayPassAllCuts_%i",iType),";centrality;#delta", 8,0.,80);
        fProfileDeltaLambdaPionDecayPassAllCuts[iType]   = new TProfile(Form("fProfileDeltaLambdaPionDecayPassAllCuts_%i",iType)  ,";centrality;#delta", 8,0.,80);
        fProfileGammaLambdaProtonDecayPassAllCuts[iType] = new TProfile(Form("fProfileGammaLambdaProtonDecayPassAllCuts_%i",iType),";centrality;#gamma", 8,0.,80);
        fProfileGammaLambdaPionDecayPassAllCuts[iType]   = new TProfile(Form("fProfileGammaLambdaPionDecayPassAllCuts_%i",iType)  ,";centrality;#gamma", 8,0.,80);
        fResultsList->Add(fProfileDeltaLambdaProtonDecayPassAllCuts[iType]);
        fResultsList->Add(fProfileDeltaLambdaPionDecayPassAllCuts[iType]);
        fResultsList->Add(fProfileGammaLambdaProtonDecayPassAllCuts[iType]);
        fResultsList->Add(fProfileGammaLambdaPionDecayPassAllCuts[iType]);
      }
    }
  }

  PostData(2,fResultsList);
  if (fDebug) Printf("Post fResultsList Data Success!");
}

//------------------------------------------------

void AliAnalysisTaskCVEPIDCME::UserExec(Option_t *)
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
  fHist2Psi2TPCPosCent -> Fill(fCent, fPsi2TPCPos);
  fHist2Psi2TPCNegCent -> Fill(fCent, fPsi2TPCNeg);
  fProfileTPCPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2TPCNeg - fPsi2TPCPos)));
  if (isUseVZEROPlane) {
    fHist2Psi2V0CCent -> Fill(fCentSPD1, fPsi2V0C);
    fHist2Psi2V0ACent -> Fill(fCentSPD1, fPsi2V0A);
    fProfileV0MPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0C - fPsi2V0A)));
    fProfileV0CTPCPosPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0C - fPsi2TPCPos)));
    fProfileV0ATPCPosPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0A - fPsi2TPCPos)));
    fProfileV0CTPCNegPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0C - fPsi2TPCNeg)));
    fProfileV0ATPCNegPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0A - fPsi2TPCNeg)));
  }
  if (isUseZDCPlane) {
    fHist2Psi1ZNCCent    -> Fill(fCent, fPsi1ZNC);
    fHist2Psi1ZNACent    -> Fill(fCent, fPsi1ZNA);
    fProfileZDCPsi1Correlation -> Fill(fCent, TMath::Cos(1*(fPsi1ZNC - fPsi1ZNA)));
    fProfileZDCPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi1ZNC - fPsi1ZNA)));
  }
  fEvtCount->Fill(19);
  //----------------------------
  // Get Lambda Vector
  //----------------------------
  if (!LoopV0s()) return;
  fEvtCount->Fill(20);
  if (fDebug) Printf("Get Lambda Vector done!");
  //----------------------------
  // Pair
  //----------------------------
  if (!PairV0Trk()) return;
  fEvtCount->Fill(21);
  if (fDebug) Printf("Pair V0 & Trk done!");

  if (!PairV0V0()) return;
  if (fDebug) Printf("Pair V0 & V0 done!");
  fEvtCount->Fill(22);

  if (!PairTrkTrk()) return;
  if (fDebug) Printf("Pair Trk & Trk done!");
  fEvtCount->Fill(23);
  //------------------
  // Post output data.
  //------------------
  PostData(1,fQAList);
  PostData(2,fResultsList);
  if (fDebug) Printf("analysis done!");
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCME::GetVZEROPlane()
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

    fProfileV0CQxCent[1]->Fill(fCentSPD1, qxRC[1]);
    fProfileV0CQyCent[1]->Fill(fCentSPD1, qyRC[1]);
    fProfileV0CQxVtx[1]->Fill(fVertex[2], qxRC[1]);
    fProfileV0CQyVtx[1]->Fill(fVertex[2], qyRC[1]);
    fHist2CalibPsi2V0CCent[1]->Fill(fCentSPD1, psi2GE[1]);
    //V0A
    fProfileV0AQxCent[0]->Fill(fCentSPD1, qxGE[2]);
    fProfileV0AQyCent[0]->Fill(fCentSPD1, qyGE[2]);
    fProfileV0AQxVtx[0]->Fill(fVertex[2], qxGE[2]);
    fProfileV0AQyVtx[0]->Fill(fVertex[2], qyGE[2]);
    fHist2CalibPsi2V0ACent[0]->Fill(fCentSPD1, psi2GE[2]);

    fProfileV0AQxCent[1]->Fill(fCentSPD1, qxRC[2]);
    fProfileV0AQyCent[1]->Fill(fCentSPD1, qyRC[2]);
    fProfileV0AQxVtx[1]->Fill(fVertex[2], qxRC[2]);
    fProfileV0AQyVtx[1]->Fill(fVertex[2], qyRC[2]);
    fHist2CalibPsi2V0ACent[1]->Fill(fCentSPD1, psi2RC[2]);
  }
  fPsi2V0C = psi2RC[1];
  fPsi2V0A = psi2RC[2];
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCME::GetZDCPlane()
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
bool AliAnalysisTaskCVEPIDCME::GetZDCPlaneLsFit()
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

bool AliAnalysisTaskCVEPIDCME::LoopTracks()
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
    fHist2EtaPhi[0]->Fill(eta,phi);

    double weight=1;
    if (isDoNUE) {
      double wEffi = GetNUECor(charge, pt);
      if (wEffi<0) continue;
      else weight *= wEffi;
    }
    if (isDoNUA) {
      double wAcc = GetNUACor(charge, phi, eta, fVertex[2]);
      if (wAcc<0) continue;
      else weight *= wAcc;
      fHistPhi[1]->Fill(phi, wAcc);
      fHist2EtaPhi[1]->Fill(eta,phi,wAcc);
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
      if (eta >= fEtaGapPos) {
        fSumQ2xTPCPos  += weight * TMath::Cos(2 * phi);
        fSumQ2yTPCPos  += weight * TMath::Sin(2 * phi);
        fWgtMultTPCPos += weight;
        std::vector<double> vec_phi_weight(2);
        vec_phi_weight.emplace_back(phi);
        vec_phi_weight.emplace_back(weight);
        mapTPCPosTrksIDPhiWgt[id] = vec_phi_weight;
      } else if (eta <= fEtaGapNeg) {
        fSumQ2xTPCNeg  += weight * TMath::Cos(2 * phi);
        fSumQ2yTPCNeg  += weight * TMath::Sin(2 * phi);
        fWgtMultTPCNeg += weight;
        std::vector<double> vec_phi_weight(2);
        vec_phi_weight.emplace_back(phi);
        vec_phi_weight.emplace_back(weight);
        mapTPCNegTrksIDPhiWgt[id] = vec_phi_weight;
      }
    }
    
    // but we need to set the dca cut for 768 when we start to choose the paiticle for pair(just for NarrowDCACut)
    if (fFilterBit == 768 && isNarrowDcaCuts768) { 
      if (fabs(dcaz) > 2.0) continue;
      if (fabs(dcaxy) > 7.0 * (0.0026 + 0.005/TMath::Power(pt, 1.01))) continue;
    }

    bool isItProttrk = CheckPIDofParticle(track,3); // 3=proton
    isItProttrk = isItProttrk && (pt < fProtonPtMax && pt > fProtonPtMin);
    if(isStrictestProtonCut) {
      // Proton need a customized dca cut
      isItProttrk = isItProttrk && (fabs(dcaz) < 1. && fabs(dcaxy) < (0.0105 + 0.035/TMath::Power(pt,1.1))); 
    }

    int code = 0;
    if(charge > 0) {
      code = 999;
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
        fHist2ProtonSigTOF -> Fill(pt,nSigRMS);
      }
    } else {
      code =-999;
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
        fHist2AntiProtonSigTOF -> Fill(pt, nSigRMS);
      }
    }

    vecParticle.emplace_back(std::array<double,6>{pt,eta,phi,(double)id,(double)code,weight});

  }
  if(fabs(fSumQ2xTPCNeg)<1.e-6 || fabs(fSumQ2yTPCNeg)<1.e-6) return false;
  if(fabs(fSumQ2xTPCPos)<1.e-6 || fabs(fSumQ2yTPCPos)<1.e-6) return false;
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCME::GetTPCPlane()
{
  double psi2TPCPos = GetEventPlane(fSumQ2xTPCPos,fSumQ2yTPCPos,2);
  double psi2TPCNeg = GetEventPlane(fSumQ2xTPCNeg,fSumQ2yTPCNeg,2);
  if (TMath::IsNaN(psi2TPCPos) || TMath::IsNaN(psi2TPCNeg)) return false;
  fPsi2TPCPos = psi2TPCPos;
  fPsi2TPCNeg = psi2TPCNeg;
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCME::LoopV0s()
{
  int nV0s = fAOD->GetNumberOfV0s();
  for (int iV0 = 0; iV0 < nV0s; iV0++) {
    AliAODv0 *v0 = fAOD->GetV0(iV0);
    if (!v0) {
      AliError(Form("%s: Could not get v0s", GetName()));
      continue;
    }
    //Basic kinematic variable
    double pt      = v0->Pt();
    double eta     = v0->PseudoRapV0();
    double dcaToPV = v0->DcaV0ToPrimVertex();//DCA to Primary Vertex
    double CPA     = v0->CosPointingAngle(fVertex);//cosine pointing angle
    double dl      = v0->DecayLengthV0(fVertex);
    fHistV0Pt              -> Fill(pt);
    fHistV0Eta             -> Fill(eta);
    fHistV0DcatoPrimVertex -> Fill(dcaToPV);
    fHistV0CPA             -> Fill(CPA);
    fHistV0DecayLength     -> Fill(dl);
    //V0 cut
    if (!IsGoodV0(v0)) continue;
    //V0 daughters cut
    AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));
    AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));
    if (!(IsGoodDaughterTrack(nTrack)) || !(IsGoodDaughterTrack(pTrack))) continue;
    float nDcaPV = v0->DcaNegToPrimVertex();
    float pDcaPV = v0->DcaPosToPrimVertex();
    if ( nDcaPV<fDaughtersDCAToPrimVtxMin || pDcaPV<fDaughtersDCAToPrimVtxMin) continue;
    int code = GetLambdaCode(pTrack,nTrack);
    if (TMath::Abs(code) != 3122) continue;
    double phi = v0->Phi();
    int id_daughter_1 = v0->GetPosID();
    int id_daughter_2 = v0->GetNegID();
    double rapLambda = v0->RapLambda();

    if (code == 3122) {
      double mass  = v0->MassLambda();
      fHistLambdaPt[0]              -> Fill(pt);
      fHistLambdaEta[0]             -> Fill(eta);
      fHistLambdaPhi[0]             -> Fill(phi);
      fHistLambdaDcaToPrimVertex[0] -> Fill(dcaToPV);
      fHistLambdaCPA[0]             -> Fill(CPA);
      fHistLambdaDecayLength[0]     -> Fill(dl);
      fHistLambdaMass[0]            -> Fill(mass);
      fHist2LambdaMassPtY[0]        -> Fill(pt, rapLambda);

      bool isLambda = true;
      isLambda = isLambda && (mass > fLambdaMassMean - fLambdaMassLeftCut);
      isLambda = isLambda && (mass < fLambdaMassMean + fLambdaMassRightCut);

      if (isLambda) {
        // if a particle has been used as daughter particle before(It happends), we have to refuse a new one.
        // if (find(vecDaughterPosID.begin(), vecDaughterPosID.end(), id_daughter_1) != vecDaughterPosID.end()) continue;
        // if (find(vecDaughterNegID.begin(), vecDaughterNegID.end(), id_daughter_2) != vecDaughterNegID.end()) continue;
        fHistLambdaPt[1]              -> Fill(pt);
        fHistLambdaEta[1]             -> Fill(eta);
        fHistLambdaPhi[1]             -> Fill(phi);
        fHistLambdaDcaToPrimVertex[1] -> Fill(dcaToPV);
        fHistLambdaCPA[1]             -> Fill(CPA);
        fHistLambdaDecayLength[1]     -> Fill(dl);
        fHistLambdaMass[1]            -> Fill(mass);
        fHist2LambdaMassPtY[1]        -> Fill(pt, rapLambda);

        vecParticleV0.emplace_back(std::array<double,9>{pt,eta,phi,0.,(double)code,1,mass,(double)id_daughter_1,(double)id_daughter_2});
        vecParticleFromDecay.emplace_back(std::array<double,6>{pTrack->Pt(),pTrack->Eta(),pTrack->Phi(),(double)id_daughter_1,(double)2212,1});
        vecParticleFromDecay.emplace_back(std::array<double,6>{nTrack->Pt(),nTrack->Eta(),nTrack->Phi(),(double)id_daughter_2,(double)-211,1});

        if(isCheckDaughterProtonPassAllCuts) {
          int id_proton = id_daughter_1;
          bool found = std::any_of(vecParticle.begin(), vecParticle.end(), [id_proton](std::array<double, 6>& arr){
              return (int)arr[3] == id_proton && (int)arr[4] == 2212;
          });
          if(found) {
            vecParticleFromDecayPassAllCuts.emplace_back(std::array<double,6>{pTrack->Pt(),pTrack->Eta(),pTrack->Phi(),(double)id_daughter_1,(double)2212,1});
            vecParticleFromDecayPassAllCuts.emplace_back(std::array<double,6>{nTrack->Pt(),nTrack->Eta(),nTrack->Phi(),(double)id_daughter_2,(double)-211,1});
          }
        }

      }
    }

    if (code == -3122) {
      double mass  = v0->MassAntiLambda();
      fHistAntiLambdaPt[0]              -> Fill(pt);
      fHistAntiLambdaEta[0]             -> Fill(eta);
      fHistAntiLambdaPhi[0]             -> Fill(phi);
      fHistAntiLambdaDcaToPrimVertex[0] -> Fill(dcaToPV);
      fHistAntiLambdaCPA[0]             -> Fill(CPA);
      fHistAntiLambdaDecayLength[0]     -> Fill(dl);
      fHistAntiLambdaMass[0]            -> Fill(mass);
      fHist2AntiLambdaMassPtY[0]        -> Fill(pt, rapLambda);

      bool isAntiLambda = true;
      isAntiLambda = isAntiLambda && (mass > fLambdaMassMean - fAntiLambdaMassLeftCut);
      isAntiLambda = isAntiLambda && (mass < fLambdaMassMean + fAntiLambdaMassRightCut);

      if (isAntiLambda) {
        //if a particle has been used as daughter particle before(It happends), we have to refuse a new one.
        fHistAntiLambdaPt[1]              -> Fill(pt);
        fHistAntiLambdaEta[1]             -> Fill(eta);
        fHistAntiLambdaPhi[1]             -> Fill(phi);
        fHistAntiLambdaDcaToPrimVertex[1] -> Fill(dcaToPV);
        fHistAntiLambdaCPA[1]             -> Fill(CPA);
        fHistAntiLambdaDecayLength[1]     -> Fill(dl);
        fHistAntiLambdaMass[1]            -> Fill(mass);
        fHist2AntiLambdaMassPtY[1]        -> Fill(pt, rapLambda);

        vecParticleV0.emplace_back(std::array<double,9>{pt,eta,phi,0.,(double)code,1,mass,(double)id_daughter_1,(double)id_daughter_2});
        vecParticleFromDecay.emplace_back(std::array<double,6>{nTrack->Pt(),nTrack->Eta(),nTrack->Phi(),(double)id_daughter_2,(double)-2212,1});
        vecParticleFromDecay.emplace_back(std::array<double,6>{pTrack->Pt(),pTrack->Eta(),pTrack->Phi(),(double)id_daughter_1,(double)211,1});

        if(isCheckDaughterProtonPassAllCuts) {
          int id_proton = id_daughter_2;
          bool found = std::any_of(vecParticle.begin(), vecParticle.end(), [id_proton](const std::array<double, 6>& arr){
              return (int)arr[3] == id_proton && (int)arr[4] == -2212;
          });
          if(found) {
            vecParticleFromDecayPassAllCuts.emplace_back(std::array<double,6>{nTrack->Pt(),nTrack->Eta(),nTrack->Phi(),(double)id_daughter_2,(double)-2212,1});
            vecParticleFromDecayPassAllCuts.emplace_back(std::array<double,6>{pTrack->Pt(),pTrack->Eta(),pTrack->Phi(),(double)id_daughter_1,(double)211,1});
          }
        }
      }
    }
  }//loop V0 end
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCME::PairV0Trk()
{
  //Lambda - X
  for (auto lambda : vecParticleV0) {
    double pt_lambda     = lambda[0];
    double eta_lambda    = lambda[1];
    double phi_lambda    = lambda[2];
    int    id_lambda     = (int)lambda[3];
    int    code_lambda   = (int)lambda[4];
    double weight_lambda = lambda[5];
    double mass_lambda   = lambda[6];
    int    id_daughter_1 = (int)lambda[7];
    int    id_daughter_2 = (int)lambda[8];

    // NOTE: Remove AutoCorrelation:
    // Get the Total sum of Pos TPC Qx,Qy locally, then Remove AutoCorr if needed.
    double repeQ2xTPCPos = 0., repeQ2yTPCPos = 0. , repeWgtMultTPCPos = 0.;
    double repeQ2xTPCNeg = 0., repeQ2yTPCNeg = 0. , repeWgtMultTPCNeg = 0.;
    // use EP from opposite eta than the charged track! One way to remove AutoCorrelation.
    double tempSumQ2xTPCPos   = fSumQ2xTPCPos;
    double tempSumQ2yTPCPos   = fSumQ2yTPCPos;
    double tempWgtMultTPCPos = fWgtMultTPCPos;

    double tempSumQ2xTPCNeg   = fSumQ2xTPCNeg;
    double tempSumQ2yTPCNeg   = fSumQ2yTPCNeg;
    double tempWgtMultTPCNeg = fWgtMultTPCNeg;

    //for TPCPos
    std::unordered_map<int, std::vector<double>> ::iterator itDaughter1TPCPos = mapTPCPosTrksIDPhiWgt.find(id_daughter_1);
    std::unordered_map<int, std::vector<double>> ::iterator itDaughter2TPCPos = mapTPCPosTrksIDPhiWgt.find(id_daughter_2);

    //for TPCNeg
    std::unordered_map<int, std::vector<double>> ::iterator itDaughter1TPCNeg = mapTPCNegTrksIDPhiWgt.find(id_daughter_1);
    std::unordered_map<int, std::vector<double>> ::iterator itDaughter2TPCNeg = mapTPCNegTrksIDPhiWgt.find(id_daughter_2);

    if(itDaughter1TPCPos != mapTPCPosTrksIDPhiWgt.end()) {
      repeQ2xTPCPos     += itDaughter1TPCPos->second[1] * TMath::Cos(2* (itDaughter1TPCPos->second)[0]);
      repeQ2yTPCPos     += itDaughter1TPCPos->second[1] * TMath::Sin(2* (itDaughter1TPCPos->second)[0]);
      repeWgtMultTPCPos += itDaughter1TPCPos->second[1];
    } 
    if(itDaughter2TPCPos != mapTPCPosTrksIDPhiWgt.end()) {
      repeQ2xTPCPos     += itDaughter2TPCPos->second[1] * TMath::Cos(2* (itDaughter2TPCPos->second)[0]);
      repeQ2yTPCPos     += itDaughter2TPCPos->second[1] * TMath::Sin(2* (itDaughter2TPCPos->second)[0]);
      repeWgtMultTPCPos += itDaughter2TPCPos->second[1];
    }

    if(itDaughter1TPCNeg != mapTPCNegTrksIDPhiWgt.end()) {
      repeQ2xTPCNeg     += itDaughter1TPCNeg->second[1] * TMath::Cos(2* (itDaughter1TPCNeg->second)[0]);
      repeQ2yTPCNeg     += itDaughter1TPCNeg->second[1] * TMath::Sin(2* (itDaughter1TPCNeg->second)[0]);
      repeWgtMultTPCNeg += itDaughter1TPCNeg->second[1];
    }

    if(itDaughter2TPCNeg != mapTPCNegTrksIDPhiWgt.end()) {
      repeQ2xTPCNeg     += itDaughter2TPCNeg->second[1] * TMath::Cos(2* (itDaughter2TPCNeg->second)[0]);
      repeQ2yTPCNeg     += itDaughter2TPCNeg->second[1] * TMath::Sin(2* (itDaughter2TPCNeg->second)[0]);
      repeWgtMultTPCNeg += itDaughter2TPCNeg->second[1];
    } 

    tempSumQ2xTPCPos  -= repeQ2xTPCPos;
    tempSumQ2yTPCPos  -= repeQ2yTPCPos;
    tempWgtMultTPCPos -= repeWgtMultTPCPos;
    tempSumQ2xTPCNeg  -= repeQ2xTPCNeg;
    tempSumQ2yTPCNeg  -= repeQ2yTPCNeg;
    tempWgtMultTPCNeg -= repeWgtMultTPCNeg;

    if (tempWgtMultTPCPos > 1.e-6 && tempWgtMultTPCNeg > 1.e-6) {
      tempSumQ2xTPCPos /= tempWgtMultTPCPos;
      tempSumQ2yTPCPos /= tempWgtMultTPCPos;
      tempSumQ2xTPCNeg /= tempWgtMultTPCNeg;
      tempSumQ2yTPCNeg /= tempWgtMultTPCNeg;
    } else continue;

    double fPsi2PosTPCNoAuto = GetEventPlane(tempSumQ2xTPCPos,tempSumQ2yTPCPos,2.);   //AutoCorrelation Removed Pos EP.
    double fPsi2NegTPCNoAuto = GetEventPlane(tempSumQ2xTPCNeg,tempSumQ2yTPCNeg,2.);   //AutoCorrelation Removed Neg EP.

    //Check the PID Flow
    if (isCalculatePIDFlow) {
      if (code_lambda == 3122) {
        fProfile2RawFlowPtCentLambda[0][0] -> Fill(fCent, pt_lambda , TMath::Cos(2 * (phi_lambda - fPsi2PosTPCNoAuto)));
        if(isUseVZEROPlane) {
          fProfile2RawFlowPtCentLambda[1][0] -> Fill(fCentSPD1, pt_lambda , TMath::Cos(2 * (phi_lambda - fPsi2V0C)));
          fProfile2RawFlowPtCentLambda[2][0] -> Fill(fCentSPD1, pt_lambda , TMath::Cos(2 * (phi_lambda - fPsi2V0A)));
        }
        if(isUseZDCPlane) {
          fProfile2RawFlowPtCentLambda[3][0] -> Fill(fCent, pt_lambda , TMath::Cos(2 * (phi_lambda - fPsi1ZNC)));
          fProfile2RawFlowPtCentLambda[4][0] -> Fill(fCent, pt_lambda , TMath::Cos(2 * (phi_lambda - fPsi1ZNA)));
        }
      }

      if (code_lambda == -3122) {
        fProfile2RawFlowPtCentLambda[0][1] -> Fill(fCent, pt_lambda , TMath::Cos(2 * (phi_lambda - fPsi2PosTPCNoAuto)));
        if(isUseVZEROPlane) {
          fProfile2RawFlowPtCentLambda[1][1] -> Fill(fCentSPD1, pt_lambda , TMath::Cos(2 * (phi_lambda - fPsi2V0C)));
          fProfile2RawFlowPtCentLambda[2][1] -> Fill(fCentSPD1, pt_lambda , TMath::Cos(2 * (phi_lambda - fPsi2V0A)));
        }
        if(isUseZDCPlane) {
          fProfile2RawFlowPtCentLambda[3][1] -> Fill(fCent, pt_lambda , TMath::Cos(2 * (phi_lambda - fPsi1ZNC)));
          fProfile2RawFlowPtCentLambda[4][1] -> Fill(fCent, pt_lambda , TMath::Cos(2 * (phi_lambda - fPsi1ZNA)));
        }
      }
    }

    for (auto particle : vecParticle) {
      double pt     = particle[0];
      double eta    = particle[1];
      double phi    = particle[2];
      int    id     = (int)particle[3];
      int    code   = (int)particle[4];
      double weight = particle[5];
      if (id == id_daughter_1 || id == id_daughter_2) continue;

      double psi2TPCNoAuto = -999;
      if (eta > 0.) psi2TPCNoAuto = fPsi2NegTPCNoAuto;
      else          psi2TPCNoAuto = fPsi2PosTPCNoAuto;

      double delta = TMath::Cos(phi_lambda - phi);
      double gamma[5] = {0.,0.,0.,0.,0.,};
      gamma[0] = TMath::Cos(phi_lambda + phi - 2 * psi2TPCNoAuto);
      gamma[1] = TMath::Cos(phi_lambda + phi - 2 * fPsi2V0C);
      gamma[2] = TMath::Cos(phi_lambda + phi - 2 * fPsi2V0A);
      gamma[3] = TMath::Cos(phi_lambda + phi - 2 * fPsi1ZNC);
      gamma[4] = TMath::Cos(phi_lambda + phi - 2 * fPsi1ZNA);

      int nBits = 4;
      TBits bitsLambdaProtonPair(nBits);
      bitsLambdaProtonPair.SetBitNumber(0, code_lambda ==  3122 && code ==  2212);
      bitsLambdaProtonPair.SetBitNumber(1, code_lambda ==  3122 && code == -2212);
      bitsLambdaProtonPair.SetBitNumber(2, code_lambda == -3122 && code ==  2212);
      bitsLambdaProtonPair.SetBitNumber(3, code_lambda == -3122 && code == -2212);

      for (int iBits = 0; iBits < nBits; iBits++) {
         if (bitsLambdaProtonPair.TestBitNumber(iBits)) {
          fProfileDeltaLambdaProton[iBits] -> Fill(fCent, delta);
          fProfileGammaLambdaProton[0][iBits] -> Fill(fCent, gamma[0]);
          if(isUseVZEROPlane) {
            fProfileGammaLambdaProton[1][iBits] -> Fill(fCentSPD1, gamma[1]);
            fProfileGammaLambdaProton[2][iBits] -> Fill(fCentSPD1, gamma[2]);
          }
          if (isUseZDCPlane) {
            fProfileGammaLambdaProton[3][iBits] -> Fill(fCent, gamma[3]);
            fProfileGammaLambdaProton[4][iBits] -> Fill(fCent, gamma[4]);
          }
          if (isCalculateDiffResult) {
            fProfile2DiffDeltaLambdaProtonDPt[iBits]  -> Fill(fCent,   pt_lambda - pt, delta);
            fProfile2DiffDeltaLambdaProtonSPt[iBits]  -> Fill(fCent,   pt_lambda + pt, delta);
            fProfile2DiffDeltaLambdaProtonDEta[iBits] -> Fill(fCent, eta_lambda - eta, delta);
            fProfile2DiffDeltaLambdaProtonMass[iBits] -> Fill(fCent,      mass_lambda, delta);
            fProfile2DiffGammaLambdaProtonDPt[iBits]  -> Fill(fCent,   pt_lambda - pt, gamma[0]);
            fProfile2DiffGammaLambdaProtonSPt[iBits]  -> Fill(fCent,   pt_lambda + pt, gamma[0]);
            fProfile2DiffGammaLambdaProtonDEta[iBits] -> Fill(fCent, eta_lambda - eta, gamma[0]);
            fProfile2DiffGammaLambdaProtonMass[iBits] -> Fill(fCent,      mass_lambda, gamma[0]);
          }
          if (isCalculateDeltaPhiSumPhi) {
            double dphi_ranged = RangeDPhi(phi_lambda - phi);
            double sphi_ranged = RangeDPhi(phi_lambda + phi - 2 * psi2TPCNoAuto);
            fHist2DEtaDPhiLambdaProton[fCentBin][iBits] -> Fill(eta_lambda - eta, dphi_ranged);
            fHist2DEtaSPhiLambdaProton[fCentBin][iBits] -> Fill(eta_lambda - eta, sphi_ranged);
          }
        }
      }

      TBits bitsLambdaHadronPair(nBits);
      bitsLambdaHadronPair.SetBitNumber(0, code_lambda ==  3122 && code > 0);
      bitsLambdaHadronPair.SetBitNumber(1, code_lambda ==  3122 && code < 0);
      bitsLambdaHadronPair.SetBitNumber(2, code_lambda == -3122 && code > 0);
      bitsLambdaHadronPair.SetBitNumber(3, code_lambda == -3122 && code < 0);
      for (int iBits = 0; iBits < nBits; iBits++) {
         if (bitsLambdaHadronPair.TestBitNumber(iBits)) {
          fProfileDeltaLambdaHadron[iBits] -> Fill(fCent, delta);
          fProfileGammaLambdaHadron[0][iBits] -> Fill(fCent, gamma[0]);
          if(isUseVZEROPlane) {
            fProfileGammaLambdaHadron[1][iBits] -> Fill(fCentSPD1, gamma[1]);
            fProfileGammaLambdaHadron[2][iBits] -> Fill(fCentSPD1, gamma[2]);
          }
          if (isUseZDCPlane) {
            fProfileGammaLambdaHadron[3][iBits] -> Fill(fCent, gamma[3]);
            fProfileGammaLambdaHadron[4][iBits] -> Fill(fCent, gamma[4]);
          }
          if (isCalculateDiffResult) {
            fProfile2DiffDeltaLambdaHadronDPt[iBits]  -> Fill(fCent,   pt_lambda - pt, delta);
            fProfile2DiffDeltaLambdaHadronSPt[iBits]  -> Fill(fCent,   pt_lambda + pt, delta);
            fProfile2DiffDeltaLambdaHadronDEta[iBits] -> Fill(fCent, eta_lambda - eta, delta);
            fProfile2DiffDeltaLambdaHadronMass[iBits] -> Fill(fCent,      mass_lambda, delta);
            fProfile2DiffGammaLambdaHadronDPt[iBits]  -> Fill(fCent,   pt_lambda - pt, gamma[0]);
            fProfile2DiffGammaLambdaHadronSPt[iBits]  -> Fill(fCent,   pt_lambda + pt, gamma[0]);
            fProfile2DiffGammaLambdaHadronDEta[iBits] -> Fill(fCent, eta_lambda - eta, gamma[0]);
            fProfile2DiffGammaLambdaHadronMass[iBits] -> Fill(fCent,      mass_lambda, gamma[0]);
          }
          if (isCalculateDeltaPhiSumPhi) {
            double dphi_ranged = RangeDPhi(phi_lambda - phi);
            double sphi_ranged = RangeDPhi(phi_lambda + phi - 2 * psi2TPCNoAuto);
            fHist2DEtaDPhiLambdaHadron[fCentBin][iBits] -> Fill(eta_lambda - eta, dphi_ranged);
            fHist2DEtaSPhiLambdaHadron[fCentBin][iBits] -> Fill(eta_lambda - eta, sphi_ranged);
          }
        }
      } 
    }

    for (auto particle : vecParticleFromDecay) {
      double pt     = particle[0];
      double eta    = particle[1];
      double phi    = particle[2];
      int    id     = (int)particle[3];
      int    code   = (int)particle[4];
      double weight = particle[5];
      if (id == id_daughter_1 || id == id_daughter_2) continue;

      double psi2TPCNoAuto = -999;
      if (eta > 0.) psi2TPCNoAuto = fPsi2NegTPCNoAuto;
      else          psi2TPCNoAuto = fPsi2PosTPCNoAuto;

      double delta = TMath::Cos(phi_lambda - phi);
      double gamma = TMath::Cos(phi_lambda + phi - 2 * psi2TPCNoAuto);

      int nBits = 4;
      TBits bitsLambdaProtonDecayPair(nBits);
      bitsLambdaProtonDecayPair.SetBitNumber(0, code_lambda ==  3122 && code ==  2212);
      bitsLambdaProtonDecayPair.SetBitNumber(1, code_lambda ==  3122 && code == -2212);
      bitsLambdaProtonDecayPair.SetBitNumber(2, code_lambda == -3122 && code ==  2212);
      bitsLambdaProtonDecayPair.SetBitNumber(3, code_lambda == -3122 && code == -2212);

      for (int iBits = 0; iBits < nBits; iBits++) {
        if (bitsLambdaProtonDecayPair.TestBitNumber(iBits)) {
          fProfileDeltaLambdaProtonDecay[iBits]->Fill(fCent, delta);
          fProfileGammaLambdaProtonDecay[iBits]->Fill(fCent, gamma);
        }
      }

      TBits bitsLambdaPionDecayPair(nBits);
      bitsLambdaPionDecayPair.SetBitNumber(0, code_lambda ==  3122 && code ==  211);
      bitsLambdaPionDecayPair.SetBitNumber(1, code_lambda ==  3122 && code == -211);
      bitsLambdaPionDecayPair.SetBitNumber(2, code_lambda == -3122 && code ==  211);
      bitsLambdaPionDecayPair.SetBitNumber(3, code_lambda == -3122 && code == -211);

      for (int iBits = 0; iBits < nBits; iBits++) {
        if (bitsLambdaPionDecayPair.TestBitNumber(iBits)) {
          fProfileDeltaLambdaPionDecay[iBits]->Fill(fCent, delta);
          fProfileGammaLambdaPionDecay[iBits]->Fill(fCent, gamma);
        }
      }
    }


    if(isCheckDaughterProtonPassAllCuts) {
      for (auto particle : vecParticleFromDecayPassAllCuts) {
        double pt     = particle[0];
        double eta    = particle[1];
        double phi    = particle[2];
        int    id     = (int)particle[3];
        int    code   = (int)particle[4];
        double weight = particle[5];
        if (id == id_daughter_1 || id == id_daughter_2) continue;
  
        double psi2TPCNoAuto = -999;
        if (eta > 0.) psi2TPCNoAuto = fPsi2NegTPCNoAuto;
        else          psi2TPCNoAuto = fPsi2PosTPCNoAuto;
  
        double delta = TMath::Cos(phi_lambda - phi);
        double gamma = TMath::Cos(phi_lambda + phi - 2 * psi2TPCNoAuto);
  
        int nBits = 4;
        TBits bitsLambdaProtonDecayPair(nBits);
        bitsLambdaProtonDecayPair.SetBitNumber(0, code_lambda ==  3122 && code ==  2212);
        bitsLambdaProtonDecayPair.SetBitNumber(1, code_lambda ==  3122 && code == -2212);
        bitsLambdaProtonDecayPair.SetBitNumber(2, code_lambda == -3122 && code ==  2212);
        bitsLambdaProtonDecayPair.SetBitNumber(3, code_lambda == -3122 && code == -2212);
  
        for (int iBits = 0; iBits < nBits; iBits++) {
          if (bitsLambdaProtonDecayPair.TestBitNumber(iBits)) {
            fProfileDeltaLambdaProtonDecayPassAllCuts[iBits]->Fill(fCent, delta);
            fProfileGammaLambdaProtonDecayPassAllCuts[iBits]->Fill(fCent, gamma);
          }
        }
  
        TBits bitsLambdaPionDecayPair(nBits);
        bitsLambdaPionDecayPair.SetBitNumber(0, code_lambda ==  3122 && code ==  211);
        bitsLambdaPionDecayPair.SetBitNumber(1, code_lambda ==  3122 && code == -211);
        bitsLambdaPionDecayPair.SetBitNumber(2, code_lambda == -3122 && code ==  211);
        bitsLambdaPionDecayPair.SetBitNumber(3, code_lambda == -3122 && code == -211);
  
        for (int iBits = 0; iBits < nBits; iBits++) {
          if (bitsLambdaPionDecayPair.TestBitNumber(iBits)) {
            fProfileDeltaLambdaPionDecayPassAllCuts[iBits]->Fill(fCent, delta);
            fProfileGammaLambdaPionDecayPassAllCuts[iBits]->Fill(fCent, gamma);
          }
        }
      }
    }
  }

  if (isCalculatePIDFlow) {
    for (auto particle : vecParticle) {
      double pt     = particle[0];
      double eta    = particle[1];
      double phi    = particle[2];
      int    id     = (int)particle[3];
      int    code   = (int)particle[4];
      double weight = particle[5];

      double v2_tmp_plane[5] = {0.,0.,0.,0.,0.};
      if(eta > 0.) v2_tmp_plane[0] = TMath::Cos(2 * (phi - fPsi2TPCNeg));
      else         v2_tmp_plane[0] = TMath::Cos(2 * (phi - fPsi2TPCPos));
      v2_tmp_plane[1] = TMath::Cos(2 * (phi - fPsi2V0C));
      v2_tmp_plane[2] = TMath::Cos(2 * (phi - fPsi2V0A));
      v2_tmp_plane[3] = TMath::Cos(2 * (phi - fPsi2V0C));
      v2_tmp_plane[4] = TMath::Cos(2 * (phi - fPsi2V0A));

      if(code > 0) {
        fProfile2RawFlowPtCentHadron[0][0] -> Fill(fCent, pt, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentHadron[1][0] -> Fill(fCentSPD1, pt, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentHadron[2][0] -> Fill(fCentSPD1, pt, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentHadron[3][0] -> Fill(fCent, pt, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentHadron[4][0] -> Fill(fCent, pt, v2_tmp_plane[4]);
        }
      }
      if(code < 0) {
        fProfile2RawFlowPtCentHadron[0][1] -> Fill(fCent, pt, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentHadron[1][1] -> Fill(fCentSPD1, pt, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentHadron[2][1] -> Fill(fCentSPD1, pt, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentHadron[3][1] -> Fill(fCent, pt, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentHadron[4][1] -> Fill(fCent, pt, v2_tmp_plane[4]);
        }
      }
      if(code == 2212) {
        fProfile2RawFlowPtCentProton[0][0] -> Fill(fCent, pt, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentProton[1][0] -> Fill(fCentSPD1, pt, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentProton[2][0] -> Fill(fCentSPD1, pt, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentProton[3][0] -> Fill(fCent, pt, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentProton[4][0] -> Fill(fCent, pt, v2_tmp_plane[4]);
        }
      }
      if(code == -2212) {
        fProfile2RawFlowPtCentProton[0][1] -> Fill(fCent, pt, v2_tmp_plane[0]);
        if(isUseVZEROPlane){
          fProfile2RawFlowPtCentProton[1][1] -> Fill(fCentSPD1, pt, v2_tmp_plane[1]);
          fProfile2RawFlowPtCentProton[2][1] -> Fill(fCentSPD1, pt, v2_tmp_plane[2]);
        }
        if(isUseZDCPlane){
          fProfile2RawFlowPtCentProton[3][1] -> Fill(fCent, pt, v2_tmp_plane[3]);
          fProfile2RawFlowPtCentProton[4][1] -> Fill(fCent, pt, v2_tmp_plane[4]);
        }
      }    
    }
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCME::PairV0V0()
{
  if (isCalculateLambdaLambda) {
    for (auto lambda_a : vecParticleV0) {
      double pt_a     = lambda_a[0];
      double eta_a    = lambda_a[1];
      double phi_a    = lambda_a[2];
      int    id_a     = (int)lambda_a[3];
      int    code_a   = (int)lambda_a[4];
      double weight_a = lambda_a[5];
      double mass_a   = lambda_a[6];
      int    id_a_daughter_1 = (int)lambda_a[7];
      int    id_a_daughter_2 = (int)lambda_a[8];
  
      for (auto lambda_b : vecParticleV0) {
        double pt_b     = lambda_b[0];
        double eta_b    = lambda_b[1];
        double phi_b    = lambda_b[2];
        int    id_b     = (int)lambda_b[3];
        int    code_b   = (int)lambda_b[4];
        double weight_b = lambda_b[5];
        double mass_b   = lambda_b[6];
        int    id_b_daughter_1 = (int)lambda_b[7];
        int    id_b_daughter_2 = (int)lambda_b[8];
  
        bool b_nice_pair = (id_a_daughter_1 != id_b_daughter_1);
        b_nice_pair = b_nice_pair && (id_a_daughter_1 != id_b_daughter_2);
        b_nice_pair = b_nice_pair && (id_a_daughter_2 != id_b_daughter_1);
        b_nice_pair = b_nice_pair && (id_a_daughter_2 != id_b_daughter_2);
        if(!b_nice_pair) continue;
  
        int nBits = 4;
        TBits bitsLambdaLambdaPair(nBits);
        bitsLambdaLambdaPair.SetBitNumber(0, code_a ==  3122 && code_b ==  3122);
        bitsLambdaLambdaPair.SetBitNumber(1, code_a ==  3122 && code_b == -3122);
        bitsLambdaLambdaPair.SetBitNumber(2, code_a == -3122 && code_b ==  3122);
        bitsLambdaLambdaPair.SetBitNumber(3, code_a == -3122 && code_b == -3122);
  
        double psi2TPCNoAuto = GetTPCPlaneNoAutoCorr(id_a_daughter_1,id_a_daughter_2,id_b_daughter_1,id_b_daughter_2); 
        if(TMath::IsNaN(psi2TPCNoAuto)) continue;
  
        double delta = TMath::Cos(phi_a - phi_b);
        double gamma[5] = {0.,0.,0.,0.,0.,};
        gamma[0] = TMath::Cos(phi_a + phi_b - 2 * psi2TPCNoAuto);
        gamma[1] = TMath::Cos(phi_a + phi_b - 2 * fPsi2V0C);
        gamma[2] = TMath::Cos(phi_a + phi_b - 2 * fPsi2V0A);
        gamma[3] = TMath::Cos(phi_a + phi_b - 2 * fPsi1ZNC);
        gamma[4] = TMath::Cos(phi_a + phi_b - 2 * fPsi1ZNA);
  
        for (int iBits = 0; iBits < nBits; iBits++) {
          if (bitsLambdaLambdaPair.TestBitNumber(iBits)) {
            fProfileDeltaLambdaLambda[iBits] -> Fill(fCent, delta);
            fProfileGammaLambdaLambda[0][iBits] -> Fill(fCent, gamma[0]);
            if(isUseVZEROPlane) {
              fProfileGammaLambdaLambda[1][iBits] -> Fill(fCentSPD1, gamma[1]);
              fProfileGammaLambdaLambda[2][iBits] -> Fill(fCentSPD1, gamma[2]);
            }
            if (isUseZDCPlane) {
              fProfileGammaLambdaLambda[3][iBits] -> Fill(fCent, gamma[3]);
              fProfileGammaLambdaLambda[4][iBits] -> Fill(fCent, gamma[4]);
            }
            if (isCalculateDiffResult) {
              fProfile2DiffDeltaLambdaLambdaDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, delta);
              fProfile2DiffDeltaLambdaLambdaSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, delta);
              fProfile2DiffDeltaLambdaLambdaDEta[iBits] -> Fill(fCent, eta_a - eta_b, delta);
              fProfile2DiffDeltaLambdaLambdaMass[iBits] -> Fill(fCent,        mass_a, delta);
              fProfile2DiffGammaLambdaLambdaDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, gamma[0]);
              fProfile2DiffGammaLambdaLambdaSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, gamma[0]);
              fProfile2DiffGammaLambdaLambdaDEta[iBits] -> Fill(fCent, eta_a - eta_b, gamma[0]);
              fProfile2DiffGammaLambdaLambdaMass[iBits] -> Fill(fCent,        mass_a, gamma[0]);
            }
            if (isCalculateDeltaPhiSumPhi) {
              double dphi_ranged = RangeDPhi(phi_a - phi_b);
              double sphi_ranged = RangeDPhi(phi_a + phi_b - 2 * psi2TPCNoAuto);
              fHist2DEtaDPhiLambdaLambda[fCentBin][iBits] -> Fill(eta_a - eta_b, dphi_ranged);
              fHist2DEtaSPhiLambdaLambda[fCentBin][iBits] -> Fill(eta_a - eta_b, sphi_ranged);
            }
          }
        }
      }
    }
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCME::PairTrkTrk()
{
  if (isCalculateProtonProton || isCalculateHadronHadron) {
    for (auto paticle_a : vecParticle) {
      double pt_a     = paticle_a[0];
      double eta_a    = paticle_a[1];
      double phi_a    = paticle_a[2];
      int    id_a     = (int)paticle_a[3];
      int    code_a   = (int)paticle_a[4];
      double weight_a = paticle_a[5];
  
      for (auto paticle_b : vecParticle) {
        double pt_b     = paticle_b[0];
        double eta_b    = paticle_b[1];
        double phi_b    = paticle_b[2];
        int    id_b     = (int)paticle_b[3];
        int    code_b   = (int)paticle_b[4];
        double weight_b = paticle_b[5];

        if (id_a == id_b) continue;
      
        double psi2TPCNoAuto = -999;
        if      ( eta_a > 0. && eta_b > 0.) psi2TPCNoAuto = fPsi2TPCNeg; 
        else if ( eta_a < 0. && eta_b < 0.) psi2TPCNoAuto = fPsi2TPCNeg;
        else continue;

        double delta = TMath::Cos(phi_a - phi_b);
        double gamma[5] = {0.,0.,0.,0.,0.,};
        gamma[0] = TMath::Cos(phi_a + phi_b - 2 * psi2TPCNoAuto);
        gamma[1] = TMath::Cos(phi_a + phi_b - 2 * fPsi2V0C);
        gamma[2] = TMath::Cos(phi_a + phi_b - 2 * fPsi2V0A);
        gamma[3] = TMath::Cos(phi_a + phi_b - 2 * fPsi1ZNC);
        gamma[4] = TMath::Cos(phi_a + phi_b - 2 * fPsi1ZNA);

        int nBits = 4;
        if (isCalculateProtonProton) {
          TBits bitsProtonProtonPair(nBits);
          bitsProtonProtonPair.SetBitNumber(0, code_a ==  2212 && code_b ==  2212);
          bitsProtonProtonPair.SetBitNumber(1, code_a ==  2212 && code_b == -2212);
          bitsProtonProtonPair.SetBitNumber(2, code_a == -2212 && code_b ==  2212);
          bitsProtonProtonPair.SetBitNumber(3, code_a == -2212 && code_b == -2212);
          for (int iBits = 0; iBits < nBits; iBits++) {
           if (bitsProtonProtonPair.TestBitNumber(iBits)) {
             fProfileDeltaProtonProton[iBits] -> Fill(fCent, delta);
             fProfileGammaProtonProton[0][iBits] -> Fill(fCent, gamma[0]);
             if(isUseVZEROPlane) {
               fProfileGammaProtonProton[1][iBits] -> Fill(fCentSPD1, gamma[1]);
               fProfileGammaProtonProton[2][iBits] -> Fill(fCentSPD1, gamma[2]);
             }
             if (isUseZDCPlane) {
               fProfileGammaProtonProton[3][iBits] -> Fill(fCent, gamma[3]);
               fProfileGammaProtonProton[4][iBits] -> Fill(fCent, gamma[4]);
             }
             if (isCalculateDiffResult) {
               fProfile2DiffDeltaProtonProtonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, delta);
               fProfile2DiffDeltaProtonProtonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, delta);
               fProfile2DiffDeltaProtonProtonDEta[iBits] -> Fill(fCent, eta_a - eta_b, delta);
               fProfile2DiffGammaProtonProtonDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, gamma[0]);
               fProfile2DiffGammaProtonProtonSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, gamma[0]);
               fProfile2DiffGammaProtonProtonDEta[iBits] -> Fill(fCent, eta_a - eta_b, gamma[0]);
             }
             if (isCalculateDeltaPhiSumPhi) {
               double dphi_ranged = RangeDPhi(phi_a - phi_b);
               double sphi_ranged = RangeDPhi(phi_a + phi_b - 2 * psi2TPCNoAuto);
               fHist2DEtaDPhiProtonProton[fCentBin][iBits] -> Fill(eta_a - eta_b, dphi_ranged);
               fHist2DEtaSPhiProtonProton[fCentBin][iBits] -> Fill(eta_a - eta_b, sphi_ranged);
             }
           }
          }
        }

        if (isCalculateHadronHadron) {
          TBits bitsHadronHadronPair(nBits);
          bitsHadronHadronPair.SetBitNumber(0, code_a > 0 && code_b > 0);
          bitsHadronHadronPair.SetBitNumber(1, code_a > 0 && code_b < 0);
          bitsHadronHadronPair.SetBitNumber(2, code_a < 0 && code_b > 0);
          bitsHadronHadronPair.SetBitNumber(3, code_a < 0 && code_b < 0);
          for (int iBits = 0; iBits < nBits; iBits++) {
            if (bitsHadronHadronPair.TestBitNumber(iBits)) {
              fProfileDeltaHadronHadron[iBits] -> Fill(fCent, delta);
              fProfileGammaHadronHadron[0][iBits] -> Fill(fCent, gamma[0]);
              if(isUseVZEROPlane) {
                fProfileGammaHadronHadron[1][iBits] -> Fill(fCentSPD1, gamma[1]);
                fProfileGammaHadronHadron[2][iBits] -> Fill(fCentSPD1, gamma[2]);
              }
              if (isUseZDCPlane) {
                fProfileGammaHadronHadron[3][iBits] -> Fill(fCent, gamma[3]);
                fProfileGammaHadronHadron[4][iBits] -> Fill(fCent, gamma[4]);
              }
              if (isCalculateDiffResult) {
                fProfile2DiffDeltaHadronHadronDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, delta);
                fProfile2DiffDeltaHadronHadronSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, delta);
                fProfile2DiffDeltaHadronHadronDEta[iBits] -> Fill(fCent, eta_a - eta_b, delta);
                fProfile2DiffGammaHadronHadronDPt[iBits]  -> Fill(fCent,   pt_a - pt_b, gamma[0]);
                fProfile2DiffGammaHadronHadronSPt[iBits]  -> Fill(fCent,   pt_a + pt_b, gamma[0]);
                fProfile2DiffGammaHadronHadronDEta[iBits] -> Fill(fCent, eta_a - eta_b, gamma[0]);
              }
              if (isCalculateDeltaPhiSumPhi) {
                double dphi_ranged = RangeDPhi(phi_a - phi_b);
                double sphi_ranged = RangeDPhi(phi_a + phi_b - 2 * psi2TPCNoAuto);
                fHist2DEtaDPhiHadronHadron[fCentBin][iBits] -> Fill(eta_a - eta_b, dphi_ranged);
                fHist2DEtaSPhiHadronHadron[fCentBin][iBits] -> Fill(eta_a - eta_b, sphi_ranged);
              }
            }
          }
        }
      }
    }
  }
  return true;
}

//---------------------------------------------------

void AliAnalysisTaskCVEPIDCME::ResetVectors()
{
  fSumQ2xTPCPos = 0.;
  fSumQ2yTPCPos = 0.;
  fWgtMultTPCPos = 0.;
  fSumQ2xTPCNeg = 0.;
  fSumQ2yTPCNeg = 0.;
  fWgtMultTPCNeg = 0.;
  std::unordered_map<int, std::vector<double>>().swap(mapTPCPosTrksIDPhiWgt);
  std::unordered_map<int, std::vector<double>>().swap(mapTPCNegTrksIDPhiWgt);
  std::vector<std::array<double,6>>().swap(vecParticle);
  std::vector<std::array<double,9>>().swap(vecParticleV0);
  std::vector<std::array<double,6>>().swap(vecParticleFromDecay);
  std::vector<std::array<double,6>>().swap(vecParticleFromDecayPassAllCuts);
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCME::LoadCalibHistForThisRun()
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

bool AliAnalysisTaskCVEPIDCME::RemovalForRun1()
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

bool AliAnalysisTaskCVEPIDCME::RejectEvtMultComp() // 15o_pass1, old pile-up
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

bool AliAnalysisTaskCVEPIDCME::RejectEvtTFFit()
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
  
  // //
  // Int_t tpcClsTot = fAOD->GetNumberOfTPCClusters();
  // Float_t nclsDif = Float_t(tpcClsTot) - (53182.6 + 113.326*multV0Tot - 0.000831275*multV0Tot*multV0Tot);

  // pile-up cuts
  if (fCentSPD0 < fCenCutLowPU->Eval(fCentV0M)) return false;
  if (fCentSPD0 > fCenCutHighPU->Eval(fCentV0M)) return false;
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) return false;
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;
  if (Float_t(multTrk) < fMultCutPU->Eval(fCentV0M)) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fAOD->IsIncompleteDAQ()) return false;

  fHist2MultCentQA[1]->Fill(fCentV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCME::RejectEvtTPCITSfb32TOF ()
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

bool AliAnalysisTaskCVEPIDCME::AODPileupCheck()
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

bool AliAnalysisTaskCVEPIDCME::PileUpMultiVertex()
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

double AliAnalysisTaskCVEPIDCME::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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

bool AliAnalysisTaskCVEPIDCME::AcceptAODTrack(AliAODTrack *track)
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

bool AliAnalysisTaskCVEPIDCME::CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck)
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
  int  trkChargePID = ftrack->Charge();

  ///Pion =>
  if (pidToCheck==1) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);//Some warning show here (***TDatabasePDG::AddParicle: particle with PDGcode = 3124 already defind),I don't understand what happended. --chunzheng
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if (trkPtPID<=0.5 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut) return true;
    if (trkPtPID> 0.5 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut) return true;
    return false;
  }
  ///Kaon =>
  else if (pidToCheck==2) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kKaon);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kKaon);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if (trkPtPID<=0.45 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut) return true;
    if (trkPtPID> 0.45 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut) return true;
    return false;
  }
  ///proton =>
  else if (pidToCheck==3) {///
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kProton);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kProton);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if(isStrictestProtonCut) { //if Set the Strictest Cut to proton
      if(trkPtPID > 0.5 && trkPtPID < 3.) {
        if (isUsePionRejection) {
          float nSigTPCPion = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);
          float nSigTOFPion = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
          float nSigRMSPion = TMath::Sqrt(nSigTPCPion*nSigTPCPion + nSigTOFPion*nSigTOFPion);
          double mom = ftrack->P();
          if(mom < 1.e-6) return false;
          if(mom > 0.5 && TMath::Abs(nSigRMS) < 2. && TMath::Abs(nSigRMSPion) > 2.) return true;
          if(mom < 0.5 && TMath::Abs(nSigRMS) < 2. && TMath::Abs(nSigTPCPion) > 2.) return true;
          return false;
        } else {
          if(nSigRMS < 2.) return true;
          return false;
        }
      }
      return false;
    } else {
      if(trkChargePID>0 && trkPtPID<0.4) return false;
      if(trkPtPID < 0.6) {
        if(TMath::Abs(nSigTPC)<=fNSigmaTPCCut) return true;
      } else {
        if(TMath::Abs(nSigRMS)<=fNSigmaTOFCut) return true;
      }
      return false;
    }
  } else {
    Printf("\n -Ve number not allowed! Choose among: 0,1,2,3 (Charge Pion, Kaon, Proton)\n return with kFALSE \n");
    return false;
  }
}

//---------------------------------------------------

double AliAnalysisTaskCVEPIDCME::GetNUECor(int charge, double pt)
{
  double weightNUE = 1;
  if (fPeriod.EqualTo("LHC10h")) {
    if (charge>0) {
      hNUEweightPlus = (TH1D*)fListNUE->FindObject(Form("effVsPt_cent%iPlus",fCentBin));
      if (!hNUEweightPlus) return -1;
      int ptBin = hNUEweightPlus->GetXaxis()->FindBin(pt);
      if (hNUEweightPlus->GetBinContent(ptBin)>0) {
        weightNUE = 1./ hNUEweightPlus->GetBinContent(ptBin);
      }
      else return -1;
    }
    if (charge<0) {
      hNUEweightMinus = (TH1D*)fListNUE->FindObject(Form("effVsPt_cent%iMinus",fCentBin));
      if (!hNUEweightMinus) return -1;
      int ptBin = hNUEweightMinus->GetXaxis()->FindBin(pt);
      if (hNUEweightMinus->GetBinContent(ptBin)>0) {
        weightNUE = 1./ hNUEweightMinus->GetBinContent(ptBin);
      }
      else return -1;
    }
  }
  else if (fPeriod.EqualTo("LHC15o")) {
    if (charge>0) {
      //hNUEweightPlus = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgPos");
      if (!hNUEweightPlus) return -1;
      int ptBin = hNUEweightPlus->GetXaxis()->FindBin(pt);
      if (hNUEweightPlus->GetBinContent(ptBin)>0) {
        weightNUE = 1./ hNUEweightPlus->GetBinContent(ptBin);
      }
      else return -1;
    }
    if (charge<0) {
      //hNUEweightMinus = (TH1D*)fListNUE->FindObject("trkEfficiencyChrgNeg");
      if (!hNUEweightMinus) return -1;
      int ptBin = hNUEweightMinus->GetXaxis()->FindBin(pt);
      if (hNUEweightMinus->GetBinContent(ptBin)>0) {
        weightNUE = 1./ hNUEweightMinus->GetBinContent(ptBin);
      }
      else return -1;
    }
  }
  return weightNUE;
}

//---------------------------------------------------

double AliAnalysisTaskCVEPIDCME::GetNUACor(int charge, double phi, double eta, double vz)
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

bool AliAnalysisTaskCVEPIDCME::IsGoodV0(AliAODv0 *aodV0)
{
  // Offline reconstructed V0 only
  if (aodV0->GetOnFlyStatus()) return false;
  // Get daughters and check them
  AliAODTrack *myTrackPosTest = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
  AliAODTrack *myTrackNegTest = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));
  if (!myTrackPosTest || !myTrackNegTest) {
    Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
    return false;
  }
  // Unlike signs of daughters
  if (myTrackNegTest->Charge() * myTrackPosTest->Charge() > 0) return false;
  // Cosinus of pointing angle < 0.997
  double dCPA = aodV0->CosPointingAngle(fVertex);
  if (dCPA < fV0CPAMin) return false;
  // DCA of V0 < 1.5 cm
  double dV0Dca = aodV0->DcaV0ToPrimVertex();
  if (TMath::Abs(dV0Dca) > fV0DCAToPrimVtxMax) return false;
  // V0 path length before decay 3-100 cm
  double dDecayLength = aodV0->DecayLengthV0(fVertex);
  if (dDecayLength > fV0DecayLengthMax) return false;
  if (dDecayLength < fV0DecayLengthMin) return false;
  // DCA between daughters < 0.5cm
  double dDCA = aodV0->DcaV0Daughters();
  if (dDCA > fV0DcaBetweenDaughtersMax) return false;
  // Transverse momentum > 0.5 GeV/c
  double dPt = aodV0->Pt();
  if (dPt < fV0PtMin ) return false;
  // Pseudorapidity < 0.5
  double dRapidity = aodV0->RapLambda();
  if (TMath::Abs(dRapidity) > fV0RapidityMax) return false;
  return kTRUE;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCME::IsGoodDaughterTrack(const AliAODTrack *track)
{
  // TPC refit
  if (!track->IsOn(AliAODTrack::kTPCrefit)) return false;
  // No kinks
  if (int(track->GetProdVertex()->GetType()) == AliAODVertex::kKink) return false;
  // Maximum value of transverse momentum
  double dPt = track->Pt();
  if (dPt > fDaughtersPtMax) return false;
  // Maximum value of pseudorapidity
  double dEta = track->Eta();
  if (TMath::Abs(dEta) > fDaughtersEtaMax) return false;
  // Minimum number of clusters
  float nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < fDaughtersTPCNclsMin) return false;
  // Findable clusters > 0
  int findable = track->GetTPCNclsF();
  if (findable <= 0) return false;
  // [number of crossed rows]>0.8  [number of findable clusters].
  if (nCrossedRowsTPC/findable < 0.8) return false;
  return true;
}

//---------------------------------------------------

int AliAnalysisTaskCVEPIDCME::GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *nTrack)
{
  bool isLambda     = kFALSE;
  bool isAntiLambda = kFALSE;
  int  code = 0;

  //Λ-->(p+)+(π-)
  float nSigTPCPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton));//TPC p+
  float nSigTPCNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion));//TPC π-
  //(Λ-)-->(p-)+(π+)
  float nSigTPCPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion));//TPC π+
  float nSigTPCNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton));//TPC p-

  isLambda     = (nSigTPCPosProton < fV0PosProtonTPCNsigma) && (nSigTPCNegPion < fV0NegPionTPCNsigma);
  isAntiLambda = (nSigTPCNegProton < fV0NegProtonTPCNsigma) && (nSigTPCPosPion < fV0PosPionTPCNsigma);

  if (isV0DaughterUseTOF) {
    float nSigTOFPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kProton));//TOF p+
    float nSigTOFNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kPion));//TOF π-
    float nSigTOFPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kPion));//TOF π+
    float nSigTOFNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kProton));//TOF p-

    isLambda     = isLambda     && (nSigTOFPosProton < fV0PosProtonTOFNsigma) && (nSigTOFNegPion < fV0NegPionTOFNsigma);
    isAntiLambda = isAntiLambda && (nSigTOFNegProton < fV0NegProtonTOFNsigma) && (nSigTOFPosPion < fV0PosPionTOFNsigma);
  }

  if (isLambda)     code =  3122;
  if (isAntiLambda) code = -3122;
  if (isLambda && isAntiLambda) code = 0;
  return code;
}

//---------------------------------------------------

inline double AliAnalysisTaskCVEPIDCME::GetEventPlane(double qx, double qy, double harmonic)
{
  double psi = (1./harmonic)*TMath::ATan2(qy,qx);
  if (psi < 0) return psi += TMath::TwoPi()/harmonic;
  else return psi;
}

//---------------------------------------------------

inline double AliAnalysisTaskCVEPIDCME::RangeDPhi(double dphi)
{
  while (dphi >= 1.5*TMath::Pi()) dphi -= TMath::TwoPi();
  while (dphi < -0.5*TMath::Pi()) dphi += TMath::TwoPi();
  return dphi;
}

//---------------------------------------------------

double AliAnalysisTaskCVEPIDCME::GetTPCPlaneNoAutoCorr(int id_0, int id_1, int id_2, int id_3)
{
  double tempSumQ2x  = fSumQ2xTPCPos;
  double tempSumQ2y  = fSumQ2yTPCPos;
  double tempWgtMult = fWgtMultTPCPos;
  double repeQ2x = 0., repeQ2y = 0., repeWgtMult = 0.;

  std::vector<std::unordered_map<int, std::vector<double>> ::iterator> vec_it;  
  vec_it.push_back(mapTPCPosTrksIDPhiWgt.find(id_0));
  vec_it.push_back(mapTPCPosTrksIDPhiWgt.find(id_1));
  vec_it.push_back(mapTPCPosTrksIDPhiWgt.find(id_2));
  vec_it.push_back(mapTPCPosTrksIDPhiWgt.find(id_3));

  for (auto it : vec_it) {
    if(it != mapTPCPosTrksIDPhiWgt.end()) {
      repeQ2x     += it->second[1] * TMath::Cos(2 * (it->second)[0]);
      repeQ2y     += it->second[1] * TMath::Sin(2 * (it->second)[0]);
      repeWgtMult += it->second[1];
    }
  }
    
  tempSumQ2x  -= repeQ2x;
  tempSumQ2y  -= repeQ2y;
  tempWgtMult -= repeWgtMult;

  double psiNoAuto = -999;

  if (tempWgtMult > 1.e-6) {
    tempSumQ2x /= tempWgtMult;
    tempSumQ2y /= tempWgtMult;
    psiNoAuto = GetEventPlane(tempSumQ2x,tempSumQ2y,2.);
  } else psiNoAuto = nan("");
  return psiNoAuto;
}

bool AliAnalysisTaskCVEPIDCME::GetDCA(double &dcaxy, double &dcaz, AliAODTrack* track) {
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

