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
// CVE analysis
// Contributor: Chunzheng Wang, <chunzheng.wang@cern.ch>, Shanghai
//--------------------------------------------------------------------------------

#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <algorithm>
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
#include "AliAnalysisTaskLambdaProtonCVE.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskLambdaProtonCVE);

//---------------------------------------------------
AliAnalysisTaskLambdaProtonCVE::AliAnalysisTaskLambdaProtonCVE() :
  AliAnalysisTaskSE(),
  fDebug(0),
  fTrigger("kMB"),
  fPeriod("LHC10h"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  IsVZEROCalibOn(false),
  IsZDCCalibOn(false),
  IsTPCCalibOn(false),
  IsQAVZERO(false),
  IsQAZDC(false),
  IsQATPC(false),
  fPlanePtMin(0.2),
  fPlanePtMax(2.0),
  fEtaGapPos( 0.1),
  fEtaGapNeg(-0.1),
  fFilterBit(768),
  fNclsCut(70),
  fChi2Max(4.0),
  fChi2Min(0.1),
  fDcaCutz(3.2),
  fDcaCutxy(2.4),
  fPtMin(0.2),
  fPtMax(5.0),
  IsDoNUE(true),
  IsDoNUA(true),
  fProtonPtMax(3.0),
  fNSigmaTPCCut(4),
  fNSigmaTOFCut(4),
  fV0PtMin(0.5),
  fV0CPAMin(0.995),
  fV0RapidityMax(0.5),
  fV0DecayLengthMin(3.),
  fV0DecayLengthMax(100.),
  fV0DCAToPrimVtxMax(1.5),
  fV0DcaBetweenDaughtersMax(1.),
  fDaughtersPtMax(20.),
  fDaughtersEtaMax(0.8),
  fDaughtersTPCNclsMin(70),
  fDaughtersDCAToPrimVtxMin(0.02),
  fV0PosProtonTPCNsigma(4.),
  fV0NegPionTPCNsigma(4.),
  fV0NegProtonTPCNsigma(4.),
  fV0PosPionTPCNsigma(4.),
  IsV0DaughterUseTOF(false),
  fV0PosProtonTOFNsigma(4.),
  fV0NegPionTOFNsigma(4.),
  fV0NegProtonTOFNsigma(4.),
  fV0PosPionTOFNsigma(4.),
  fLambdaMassCut(0.005),
  IsCheckPIDFlow(false),
  fMassMean(1.115683),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fAOD(nullptr),         //! aod Event
  fPIDResponse(nullptr), //! PID Handler
  fUtils(nullptr),       //! Event Selection Options
  runNumList(0),
  fRunNum(-999), // runnumber
  fOldRunNum(-999), // old runnumber
  fRunNumBin(-999), // runnumer bin, 10:139510..., 11:170387..., 15HIR:246994...
  fVzBin(-999), // vertex z bin
  fCentBin(-999), // centrality bin: 0-10
  fCent(-999), // value of centrality
  fPsi1ZNC(-999),
  fPsi1ZNA(-999),
  fPsi2V0C(-999),
  fPsi2V0A(-999),
  fPsi2TPCPos(-999),
  fPsi2TPCNeg(-999),
  SumQ2xTPCPos(0.),
  SumQ2yTPCPos(0.),
  fWgtMultTPCPos(0.),
  SumQ2xTPCNeg(0.),
  SumQ2yTPCNeg(0.),
  fWgtMultTPCNeg(0.),
  vecPosEPTrkID(0),
  vecNegEPTrkID(0),
  vecPosEPTrkPhi(0),
  vecNegEPTrkPhi(0),
  vecPosEPTrkNUAWgt(0),
  vecNegEPTrkNUAWgt(0),
  vecPDGCode(0),
  vecID(0),
  vecPhi(0),
  vecEta(0),
  vecPt(0),
  vecNUAWeight(0),
  vecNUEWeight(0),
  vecNUAWeightPID(0),
  vecNUEWeightPID(0),
  vecLambdaCode(0),
  vecLambdaPhi(0),
  vecLambdaPt(0),
  vecDaughterPosID(0),
  vecDaughterNegID(0),
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
  fOutputList(nullptr),
  fEvtCount(nullptr),
  fHistRunNumBin(nullptr),
  fHistPt(nullptr),
  fHistEta(nullptr),
  fHistNhits(nullptr),
  fHist2DPDedx(nullptr),
  fHistDcaXY(nullptr),
  fHistDcaZ(nullptr),
  fHistV0Pt(nullptr),
  fHistV0Eta(nullptr),
  fHistV0DcatoPrimVertex(nullptr),
  fHistV0CPA(nullptr),
  fHistV0DecayLength(nullptr),
  fHist2DPsi2TPCPosCent(nullptr),
  fHist2DPsi2TPCNegCent(nullptr),
  fHist2DPsi2V0CCent(nullptr),
  fHist2DPsi2V0ACent(nullptr),
  fHist2DPsi1ZNCCent(nullptr),
  fHist2DPsi1ZNACent(nullptr),
  fProfileTPCPsi2Correlation(nullptr),
  fProfileV0MPsi2Correlation(nullptr),
  fProfileZDCPsi1Correlation(nullptr),
  fProfileZDCPsi2Correlation(nullptr),
  fProfileV0CTPCPosPsi2Correlation(nullptr),
  fProfileV0ATPCPosPsi2Correlation(nullptr),
  fProfileV0CTPCNegPsi2Correlation(nullptr),
  fProfileV0ATPCNegPsi2Correlation(nullptr),
  fProfileDelta_Lambda_hPos(nullptr),
  fProfileDelta_Lambda_hNeg(nullptr),
  fProfileDelta_Lambda_Proton(nullptr),
  fProfileDelta_Lambda_AntiProton(nullptr),
  fProfileDelta_AntiLambda_hPos(nullptr),
  fProfileDelta_AntiLambda_hNeg(nullptr),
  fProfileDelta_AntiLambda_Proton(nullptr),
  fProfileDelta_AntiLambda_AntiProton(nullptr),
  fProfileGammaTPC_Lambda_hPos(nullptr),
  fProfileGammaTPC_Lambda_hNeg(nullptr),
  fProfileGammaTPC_Lambda_Proton(nullptr),
  fProfileGammaTPC_Lambda_AntiProton(nullptr),
  fProfileGammaTPC_AntiLambda_hPos(nullptr),
  fProfileGammaTPC_AntiLambda_hNeg(nullptr),
  fProfileGammaTPC_AntiLambda_Proton(nullptr),
  fProfileGammaTPC_AntiLambda_AntiProton(nullptr),
  fProfileGammaV0C_Lambda_hPos(nullptr),
  fProfileGammaV0C_Lambda_hNeg(nullptr),
  fProfileGammaV0C_Lambda_Proton(nullptr),
  fProfileGammaV0C_Lambda_AntiProton(nullptr),
  fProfileGammaV0C_AntiLambda_hPos(nullptr),
  fProfileGammaV0C_AntiLambda_hNeg(nullptr),
  fProfileGammaV0C_AntiLambda_Proton(nullptr),
  fProfileGammaV0C_AntiLambda_AntiProton(nullptr),
  fProfileGammaV0A_Lambda_hPos(nullptr),
  fProfileGammaV0A_Lambda_hNeg(nullptr),
  fProfileGammaV0A_Lambda_Proton(nullptr),
  fProfileGammaV0A_Lambda_AntiProton(nullptr),
  fProfileGammaV0A_AntiLambda_hPos(nullptr),
  fProfileGammaV0A_AntiLambda_hNeg(nullptr),
  fProfileGammaV0A_AntiLambda_Proton(nullptr),
  fProfileGammaV0A_AntiLambda_AntiProton(nullptr),
  fProfileGammaZNC_Lambda_hPos(nullptr),
  fProfileGammaZNC_Lambda_hNeg(nullptr),
  fProfileGammaZNC_Lambda_Proton(nullptr),
  fProfileGammaZNC_Lambda_AntiProton(nullptr),
  fProfileGammaZNC_AntiLambda_hPos(nullptr),
  fProfileGammaZNC_AntiLambda_hNeg(nullptr),
  fProfileGammaZNC_AntiLambda_Proton(nullptr),
  fProfileGammaZNC_AntiLambda_AntiProton(nullptr),
  fProfileGammaZNA_Lambda_hPos(nullptr),
  fProfileGammaZNA_Lambda_hNeg(nullptr),
  fProfileGammaZNA_Lambda_Proton(nullptr),
  fProfileGammaZNA_Lambda_AntiProton(nullptr),
  fProfileGammaZNA_AntiLambda_hPos(nullptr),
  fProfileGammaZNA_AntiLambda_hNeg(nullptr),
  fProfileGammaZNA_AntiLambda_Proton(nullptr),
  fProfileGammaZNA_AntiLambda_AntiProton(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;

  for (int i = 0; i < 3; ++i) pV0XMeanRead[i] = nullptr;
  for (int i = 0; i < 3; ++i) pV0YMeanRead[i] = nullptr;

  for (int i = 0; i < 2; ++i) hQx2mV0[i] = nullptr;
  for (int i = 0; i < 2; ++i) hQy2mV0[i] = nullptr;

  for (int i = 0; i < 3; ++i) vtxQuant1[i] = -999;
  for (int i = 0; i < 3; ++i) vtxQuant2[i] = -999;

  for (int i = 0; i < 2; ++i) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fHistVz[i]   = nullptr;
  for (int i = 0; i < 8; ++i) fHist2DCentQA[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHist2DMultCentQA[i] = nullptr;
  for (int i = 0; i < 6; ++i) fHist2DMultMultQA[i] = nullptr;

  for (int i = 0; i < 2; ++i) fHistPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2DEtaPhi[i] = nullptr;

  for (int i = 0; i < 2; ++i) fProfileV0CQxCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0CQyCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0CQxVtx[i]  = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0CQyVtx[i]  = nullptr;
  for (int i = 0; i < 2; ++i) fHist2DCalibPsi2V0CCent[i] = nullptr;

  for (int i = 0; i < 2; ++i) fProfileV0AQxCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0AQyCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0AQxVtx[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0AQyVtx[i] = nullptr;
  for (int i = 0; i < 2; ++i) fHist2DCalibPsi2V0ACent[i] = nullptr;

  for (int i = 0; i < 2; ++i) fProfileZNCTowerMeanEnegry[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZNCQxCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZNCQyCent[i] = nullptr;
  for (int i = 0; i < 3; ++i) fHist2DCalibPsi1ZNCCent[i] = nullptr;

  for (int i = 0; i < 2; ++i) fProfileZNATowerMeanEnegry[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZNAQxCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZNAQyCent[i] = nullptr;
  for (int i = 0; i < 3; ++i) fHist2DCalibPsi1ZNACent[i] = nullptr;

  for (int i = 0; i < 2; ++i) fProfileZDCQxAQxCCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZDCQxAQyCCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZDCQyAQxCCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZDCQyAQyCCent[i] = nullptr;

  for (int i = 0; i < 2; ++i) fHistLambdaPt[i]                  = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaEta[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaPhi[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaDcaToPrimVertex[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaCPA[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaDecayLength[i]         = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaMass[i]                = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaPt[i]              = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaEta[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaPhi[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaCPA[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaDecayLength[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaMass[i]            = nullptr;
  for (int i = 0; i < 2; ++i) fProfileLambdaMassVsPt[i]         = nullptr;
  for (int i = 0; i < 2; ++i) fProfileAntiLambdaMassVsPt[i]     = nullptr;

  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPthPos[i]       = nullptr;
  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPthNeg[i]       = nullptr;
  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPtProton[i]     = nullptr;
  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPtAntiProton[i] = nullptr;
  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPtLambda[i]     = nullptr;
  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPtAntiLambda[i] = nullptr;
}

//---------------------------------------------------
AliAnalysisTaskLambdaProtonCVE::AliAnalysisTaskLambdaProtonCVE(const char *name) :
  AliAnalysisTaskSE(name),
  fDebug(0),
  fTrigger("kMB"),
  fPeriod("LHC10h"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  IsVZEROCalibOn(false),
  IsZDCCalibOn(false),
  IsTPCCalibOn(false),
  IsQAVZERO(false),
  IsQAZDC(false),
  IsQATPC(false),
  fPlanePtMin(0.2),
  fPlanePtMax(2.0),
  fEtaGapPos( 0.1),
  fEtaGapNeg(-0.1),
  fFilterBit(768),
  fNclsCut(70),
  fChi2Max(4.0),
  fChi2Min(0.1),
  fDcaCutz(3.2),
  fDcaCutxy(2.4),
  fPtMin(0.2),
  fPtMax(5.0),
  IsDoNUE(true),
  IsDoNUA(true),
  fProtonPtMax(3.0),
  fNSigmaTPCCut(4),
  fNSigmaTOFCut(4),
  fV0PtMin(0.5),
  fV0CPAMin(0.995),
  fV0RapidityMax(0.5),
  fV0DecayLengthMin(3.),
  fV0DecayLengthMax(100.),
  fV0DCAToPrimVtxMax(1.5),
  fV0DcaBetweenDaughtersMax(1.),
  fDaughtersPtMax(20.),
  fDaughtersEtaMax(0.8),
  fDaughtersTPCNclsMin(70),
  fDaughtersDCAToPrimVtxMin(0.02),
  fV0PosProtonTPCNsigma(4.),
  fV0NegPionTPCNsigma(4.),
  fV0NegProtonTPCNsigma(4.),
  fV0PosPionTPCNsigma(4.),
  IsV0DaughterUseTOF(false),
  fV0PosProtonTOFNsigma(4.),
  fV0NegPionTOFNsigma(4.),
  fV0NegProtonTOFNsigma(4.),
  fV0PosPionTOFNsigma(4.),
  fLambdaMassCut(0.005),
  IsCheckPIDFlow(false),
  fMassMean(1.115683),
  fEtaCut(0.8),
  fDedxCut(10.0),
  fAOD(nullptr),         //! aod Event
  fPIDResponse(nullptr), //! PID Handler
  fUtils(nullptr),       //! Event Selection Options
  runNumList(0),
  fRunNum(-999), // runnumber
  fOldRunNum(-999), // old runnumber
  fRunNumBin(-999), // runnumer bin, 10:139510..., 11:170387..., 15HIR:246994...
  fVzBin(-999), // vertex z bin
  fCentBin(-999), // centrality bin: 0-10
  fCent(-999), // value of centrality
  fPsi1ZNC(-999),
  fPsi1ZNA(-999),
  fPsi2V0C(-999),
  fPsi2V0A(-999),
  fPsi2TPCPos(-999),
  fPsi2TPCNeg(-999),
  SumQ2xTPCPos(0.),
  SumQ2yTPCPos(0.),
  fWgtMultTPCPos(0.),
  SumQ2xTPCNeg(0.),
  SumQ2yTPCNeg(0.),
  fWgtMultTPCNeg(0.),
  vecPosEPTrkID(0),
  vecNegEPTrkID(0),
  vecPosEPTrkPhi(0),
  vecNegEPTrkPhi(0),
  vecPosEPTrkNUAWgt(0),
  vecNegEPTrkNUAWgt(0),
  vecPDGCode(0),
  vecID(0),
  vecPhi(0),
  vecEta(0),
  vecPt(0),
  vecNUAWeight(0),
  vecNUEWeight(0),
  vecNUAWeightPID(0),
  vecNUEWeightPID(0),
  vecLambdaCode(0),
  vecLambdaPhi(0),
  vecLambdaPt(0),
  vecDaughterPosID(0),
  vecDaughterNegID(0),
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
  fOutputList(nullptr),
  fEvtCount(nullptr),
  fHistRunNumBin(nullptr),
  fHistPt(nullptr),
  fHistEta(nullptr),
  fHistNhits(nullptr),
  fHist2DPDedx(nullptr),
  fHistDcaXY(nullptr),
  fHistDcaZ(nullptr),
  fHistV0Pt(nullptr),
  fHistV0Eta(nullptr),
  fHistV0DcatoPrimVertex(nullptr),
  fHistV0CPA(nullptr),
  fHistV0DecayLength(nullptr),
  fHist2DPsi2TPCPosCent(nullptr),
  fHist2DPsi2TPCNegCent(nullptr),
  fHist2DPsi2V0CCent(nullptr),
  fHist2DPsi2V0ACent(nullptr),
  fHist2DPsi1ZNCCent(nullptr),
  fHist2DPsi1ZNACent(nullptr),
  fProfileTPCPsi2Correlation(nullptr),
  fProfileV0MPsi2Correlation(nullptr),
  fProfileZDCPsi1Correlation(nullptr),
  fProfileZDCPsi2Correlation(nullptr),
  fProfileV0CTPCPosPsi2Correlation(nullptr),
  fProfileV0ATPCPosPsi2Correlation(nullptr),
  fProfileV0CTPCNegPsi2Correlation(nullptr),
  fProfileV0ATPCNegPsi2Correlation(nullptr),
  fProfileDelta_Lambda_hPos(nullptr),
  fProfileDelta_Lambda_hNeg(nullptr),
  fProfileDelta_Lambda_Proton(nullptr),
  fProfileDelta_Lambda_AntiProton(nullptr),
  fProfileDelta_AntiLambda_hPos(nullptr),
  fProfileDelta_AntiLambda_hNeg(nullptr),
  fProfileDelta_AntiLambda_Proton(nullptr),
  fProfileDelta_AntiLambda_AntiProton(nullptr),
  fProfileGammaTPC_Lambda_hPos(nullptr),
  fProfileGammaTPC_Lambda_hNeg(nullptr),
  fProfileGammaTPC_Lambda_Proton(nullptr),
  fProfileGammaTPC_Lambda_AntiProton(nullptr),
  fProfileGammaTPC_AntiLambda_hPos(nullptr),
  fProfileGammaTPC_AntiLambda_hNeg(nullptr),
  fProfileGammaTPC_AntiLambda_Proton(nullptr),
  fProfileGammaTPC_AntiLambda_AntiProton(nullptr),
  fProfileGammaV0C_Lambda_hPos(nullptr),
  fProfileGammaV0C_Lambda_hNeg(nullptr),
  fProfileGammaV0C_Lambda_Proton(nullptr),
  fProfileGammaV0C_Lambda_AntiProton(nullptr),
  fProfileGammaV0C_AntiLambda_hPos(nullptr),
  fProfileGammaV0C_AntiLambda_hNeg(nullptr),
  fProfileGammaV0C_AntiLambda_Proton(nullptr),
  fProfileGammaV0C_AntiLambda_AntiProton(nullptr),
  fProfileGammaV0A_Lambda_hPos(nullptr),
  fProfileGammaV0A_Lambda_hNeg(nullptr),
  fProfileGammaV0A_Lambda_Proton(nullptr),
  fProfileGammaV0A_Lambda_AntiProton(nullptr),
  fProfileGammaV0A_AntiLambda_hPos(nullptr),
  fProfileGammaV0A_AntiLambda_hNeg(nullptr),
  fProfileGammaV0A_AntiLambda_Proton(nullptr),
  fProfileGammaV0A_AntiLambda_AntiProton(nullptr),
  fProfileGammaZNC_Lambda_hPos(nullptr),
  fProfileGammaZNC_Lambda_hNeg(nullptr),
  fProfileGammaZNC_Lambda_Proton(nullptr),
  fProfileGammaZNC_Lambda_AntiProton(nullptr),
  fProfileGammaZNC_AntiLambda_hPos(nullptr),
  fProfileGammaZNC_AntiLambda_hNeg(nullptr),
  fProfileGammaZNC_AntiLambda_Proton(nullptr),
  fProfileGammaZNC_AntiLambda_AntiProton(nullptr),
  fProfileGammaZNA_Lambda_hPos(nullptr),
  fProfileGammaZNA_Lambda_hNeg(nullptr),
  fProfileGammaZNA_Lambda_Proton(nullptr),
  fProfileGammaZNA_Lambda_AntiProton(nullptr),
  fProfileGammaZNA_AntiLambda_hPos(nullptr),
  fProfileGammaZNA_AntiLambda_hNeg(nullptr),
  fProfileGammaZNA_AntiLambda_Proton(nullptr),
  fProfileGammaZNA_AntiLambda_AntiProton(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;

  for (int i = 0; i < 3; ++i) pV0XMeanRead[i] = nullptr;
  for (int i = 0; i < 3; ++i) pV0YMeanRead[i] = nullptr;

  for (int i = 0; i < 2; ++i) hQx2mV0[i] = nullptr;
  for (int i = 0; i < 2; ++i) hQy2mV0[i] = nullptr;

  for (int i = 0; i < 3; ++i) vtxQuant1[i] = -999;
  for (int i = 0; i < 3; ++i) vtxQuant2[i] = -999;

  for (int i = 0; i < 2; ++i) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fHistVz[i]   = nullptr;
  for (int i = 0; i < 8; ++i) fHist2DCentQA[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHist2DMultCentQA[i] = nullptr;
  for (int i = 0; i < 6; ++i) fHist2DMultMultQA[i] = nullptr;

  for (int i = 0; i < 2; ++i) fHistPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2DEtaPhi[i] = nullptr;

  for (int i = 0; i < 2; ++i) fProfileV0CQxCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0CQyCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0CQxVtx[i]  = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0CQyVtx[i]  = nullptr;
  for (int i = 0; i < 2; ++i) fHist2DCalibPsi2V0CCent[i] = nullptr;

  for (int i = 0; i < 2; ++i) fProfileV0AQxCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0AQyCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0AQxVtx[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileV0AQyVtx[i] = nullptr;
  for (int i = 0; i < 2; ++i) fHist2DCalibPsi2V0ACent[i] = nullptr;

  for (int i = 0; i < 2; ++i) fProfileZNCTowerMeanEnegry[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZNCQxCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZNCQyCent[i] = nullptr;
  for (int i = 0; i < 3; ++i) fHist2DCalibPsi1ZNCCent[i] = nullptr;

  for (int i = 0; i < 2; ++i) fProfileZNATowerMeanEnegry[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZNAQxCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZNAQyCent[i] = nullptr;
  for (int i = 0; i < 3; ++i) fHist2DCalibPsi1ZNACent[i] = nullptr;

  for (int i = 0; i < 2; ++i) fProfileZDCQxAQxCCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZDCQxAQyCCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZDCQyAQxCCent[i] = nullptr;
  for (int i = 0; i < 2; ++i) fProfileZDCQyAQyCCent[i] = nullptr;

  for (int i = 0; i < 2; ++i) fHistLambdaPt[i]                  = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaEta[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaPhi[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaDcaToPrimVertex[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaCPA[i]                 = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaDecayLength[i]         = nullptr;
  for (int i = 0; i < 2; ++i) fHistLambdaMass[i]                = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaPt[i]              = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaEta[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaPhi[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaCPA[i]             = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaDecayLength[i]     = nullptr;
  for (int i = 0; i < 2; ++i) fHistAntiLambdaMass[i]            = nullptr;
  for (int i = 0; i < 2; ++i) fProfileLambdaMassVsPt[i]         = nullptr;
  for (int i = 0; i < 2; ++i) fProfileAntiLambdaMassVsPt[i]     = nullptr;

  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPthPos[i]       = nullptr;
  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPthNeg[i]       = nullptr;
  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPtProton[i]     = nullptr;
  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPtAntiProton[i] = nullptr;
  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPtLambda[i]     = nullptr;
  for (int i = 0; i < 5; ++i) fProfile2DRawFlowCentPtAntiLambda[i] = nullptr;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//------------------------------------------------

AliAnalysisTaskLambdaProtonCVE::~AliAnalysisTaskLambdaProtonCVE()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  if (fOutputList) delete fOutputList;
}

//---------------------------------------------------

void AliAnalysisTaskLambdaProtonCVE::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate");
}

//---------------------------------------------------

void AliAnalysisTaskLambdaProtonCVE::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList -> SetName(GetName());
  fOutputList -> SetOwner(kTRUE);
  ////////////////////////
  // Run Number Info
  ////////////////////////
  TString runNumList10h[90] = {
    "139510","139507","139505","139503","139465","139438","139437","139360","139329","139328","139314","139310",
    "139309","139173","139107","139105","139038","139037","139036","139029","139028","138872","138871","138870",
    "138837","138732","138730","138666","138662","138653","138652","138638","138624","138621","138583","138582",
    "138579","138578","138534","138469","138442","138439","138438","138396","138364","138275","138225","138201",
    "138197","138192","138190","137848","137844","137752","137751","137724","137722","137718","137704","137693",
    "137692","137691","137686","137685","137639","137638","137608","137595","137549","137546","137544","137541",
    "137539","137531","137530","137443","137441","137440","137439","137434","137432","137431","137430","137243",
    "137236","137235","137232","137231","137162","137161"};
  TString runNumList15o[138] = {
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
  TString runNumList18q[125] = {
    "296623","296622","296621","296619","296618","296616","296615","296594","296553","296552","296551","296550",
    "296548","296547","296516","296512","296511","296510","296509","296472","296433","296424","296423","296420",
    "296419","296415","296414","296383","296381","296380","296379","296378","296377","296376","296375","296312",
    "296309","296304","296303","296280","296279","296273","296270","296269","296247","296246","296244","296243",
    "296242","296241","296240","296198","296197","296196","296195","296194","296192","296191","296143","296142",
    "296135","296134","296133","296132","296123","296074","296066","296065","296063","296062","296060","296016",
    "295942","295941","295937","295936","295913","295910","295909","295861","295860","295859","295856","295855",
    "295854","295853","295831","295829","295826","295825","295822","295819","295818","295816","295791","295788",
    "295786","295763","295762","295759","295758","295755","295754","295725","295723","295721","295719","295718",
    "295717","295714","295712","295676","295675","295673","295668","295667","295666","295615","295612","295611",
    "295610","295589","295588","295586","295585"};
  TString runNumList18r[89] = {
    "297595","297590","297588","297558","297544","297542","297541","297540","297537","297512","297483","297479",
    "297452","297451","297450","297446","297442","297441","297415","297414","297413","297406","297405","297380",
    "297379","297372","297367","297366","297363","297336","297335","297333","297332","297317","297311","297310",
    "297278","297222","297221","297218","297196","297195","297193","297133","297132","297129","297128","297124",
    "297123","297119","297118","297117","297085","297035","297031","296966","296941","296938","296935","296934",
    "296932","296931","296930","296903","296900","296899","296894","296852","296851","296850","296848","296839",
    "296838","296836","296835","296799","296794","296793","296790","296787","296786","296785","296784","296781",
    "296752","296694","296693","296691","296690"};
  runNumList = new std::map<int,int>;
  if      (fPeriod.EqualTo("LHC10h")) for (int i = 0; i < 90; i++) runNumList->insert(std::pair<int,int>(runNumList10h[i].Atoi(),i));
  else if (fPeriod.EqualTo("LHC15o")) for (int i = 0; i <138; i++) runNumList->insert(std::pair<int,int>(runNumList15o[i].Atoi(),i));
  else if (fPeriod.EqualTo("LHC18q")) for (int i = 0; i <125; i++) runNumList->insert(std::pair<int,int>(runNumList18q[i].Atoi(),i));
  else if (fPeriod.EqualTo("LHC18r")) for (int i = 0; i < 89; i++) runNumList->insert(std::pair<int,int>(runNumList18r[i].Atoi(),i));
  else return;

  fHistRunNumBin = new TH1I("runNumBin","",(int)runNumList->size(),0,(int)runNumList->size());
  std::map<int,int>::iterator iter;
  for (auto runNum : *runNumList) fHistRunNumBin->GetXaxis()->SetBinLabel(runNum.second, Form("%i",runNum.first));
  fOutputList->Add(fHistRunNumBin);

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
  if (IsDoNUE) {
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
  if (IsDoNUA) {
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
  if (IsVZEROCalibOn) {
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
      // V0C Qx Mean
      // V0C Qy Mean
      // V0A Qx Mean
      // V0A Qy Mean
      contQxncm = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqxc%im",2));
      contQyncm = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqyc%im",2));
      contQxnam = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqxa%im",2));
      contQynam = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqya%im",2));
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
  if (IsZDCCalibOn) {
    if ((fPeriod.EqualTo("LHC10h")||fPeriod.EqualTo("LHC18q")||fPeriod.EqualTo("LHC18r")) && !fListZDCCalib) {
      std::cout<<("ZDC calibration list not found")<<std::endl;
      return;
    }
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
  // event-wise
  fEvtCount = new TH1D("EvtCount", "Event Count", 25, 1, 26);
  fEvtCount->GetXaxis()->SetBinLabel(1,"All");
  fEvtCount->GetXaxis()->SetBinLabel(2,"Read in");
  fEvtCount->GetXaxis()->SetBinLabel(3,"Event");
  fEvtCount->GetXaxis()->SetBinLabel(4,"Run Number");
  fEvtCount->GetXaxis()->SetBinLabel(5,"Vertex");
  fEvtCount->GetXaxis()->SetBinLabel(6,"Centrality");
  fEvtCount->GetXaxis()->SetBinLabel(7,"Pile up");
  fEvtCount->GetXaxis()->SetBinLabel(8,"Get VZERO Plane");
  fEvtCount->GetXaxis()->SetBinLabel(9,"Get ZDC Plane");
  fEvtCount->GetXaxis()->SetBinLabel(10,"Reset Vector");
  fEvtCount->GetXaxis()->SetBinLabel(11,"Loop Track");
  fEvtCount->GetXaxis()->SetBinLabel(12,"Get TPC Plane");
  fEvtCount->GetXaxis()->SetBinLabel(13,"Resolution");
  fEvtCount->GetXaxis()->SetBinLabel(14,"Loop V0");
  fEvtCount->GetXaxis()->SetBinLabel(15,"Pair Lambda");
  fEvtCount->GetXaxis()->SetBinLabel(20,"Manager");
  fEvtCount->GetXaxis()->SetBinLabel(21,"Handler");
  fEvtCount->GetXaxis()->SetBinLabel(22,"fAOD");
  fEvtCount->GetXaxis()->SetBinLabel(23,"fPID");
  fEvtCount->GetXaxis()->SetBinLabel(24,"fUtils");
  fEvtCount->GetXaxis()->SetBinLabel(25,"fMultSel");
  fOutputList->Add(fEvtCount);

  ////////////
  //QA Plots//
  ////////////
  // Event-wise QA
  fHistCent[0] = new TH1D("fHistCentBfCut", " Dist. of Centrality Before Cut ", 100, 0., 100.);
  fHistCent[1] = new TH1D("fHistCentAfCut", " Dist. of Centrality After Cut ", 100, 0., 100.);
  fOutputList->Add(fHistCent[0]);
  fOutputList->Add(fHistCent[1]);

  fHistVz[0] = new TH1D("fHistVzBfCut", "Dist of Centrality Before Cut", 200, -50., 50.);
  fHistVz[1] = new TH1D("fHistVzAfCut", "Dist of Centrality After Cut", 200, -50., 50.);
  fOutputList->Add(fHistVz[0]);
  fOutputList->Add(fHistVz[1]);

  fHist2DCentQA[0] = new TH2D("fHist2DCentQA_V0M_SPD1_BfCut", ";centV0M;centSPD1", 100, 0, 100, 100, 0, 100);
  fHist2DCentQA[1] = new TH2D("fHist2DCentQA_V0M_SPD1_AfCut", ";centV0M;centSPD1", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(fHist2DCentQA[0]);
  fOutputList->Add(fHist2DCentQA[1]);

  fHist2DCentQA[2] = new TH2D("fHist2DCentQA_V0M_TRK_BfCut", ";centV0M;centTRK", 100, 0, 100, 100, 0, 100);
  fHist2DCentQA[3] = new TH2D("fHist2DCentQA_V0M_TRK_AfCut", ";centV0M;centTRK", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(fHist2DCentQA[2]);
  fOutputList->Add(fHist2DCentQA[3]);

  fHist2DCentQA[4] = new TH2D("fHist2DCentQA_V0M_SPD0_BfCut", ";centV0M;centSPD0", 100, 0, 100, 100, 0, 100);
  fHist2DCentQA[5] = new TH2D("fHist2DCentQA_V0M_SPD0_AfCut", ";centV0M;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(fHist2DCentQA[4]);
  fOutputList->Add(fHist2DCentQA[5]);

  fHist2DCentQA[6] = new TH2D("fHist2DCentQA_SPD1_SPD0_BfCut", ";centSPD1;centSPD0", 100, 0, 100, 100, 0, 100);
  fHist2DCentQA[7] = new TH2D("fHist2DCentQA_SPD1_SPD0_AfCut", ";centSPD1;centSPD0", 100, 0, 100, 100, 0, 100);
  fOutputList->Add(fHist2DCentQA[6]);
  fOutputList->Add(fHist2DCentQA[7]);

  fHist2DMultCentQA[0] = new TH2D("fHist2DMultCentQA_BfCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
  fHist2DMultCentQA[1] = new TH2D("fHist2DMultCentQA_AfCut", ";centV0M;multFB32", 100, 0, 100, 20, 0, 5000);
  fOutputList->Add(fHist2DMultCentQA[0]);
  fOutputList->Add(fHist2DMultCentQA[1]);

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
  fHistPt  = new TH1D("fHistPt", "", 200, 0., 20.);
  fHistEta = new TH1D("fHistEta", "", 200, -10., 10.);
  fHistNhits = new TH1D("fHistNhits", "", 200, 0., 200.);
  fHist2DPDedx = new TH2D("fHist2DPDedx", "", 400, -10., 10., 400, 0, 1000);
  fHistDcaXY = new TH1D("fHistDcaXY", "", 100, 0., 10.);
  fHistDcaZ  = new TH1D("fHistDcaZ", "", 100, 0., 10.);
  fHistPhi[0] = new TH1D("fHistPhi", "", 100, 0, TMath::TwoPi());
  fHistPhi[1] = new TH1D("fHistPhi_afterNUA", "", 100, 0, TMath::TwoPi());
  fHist2DEtaPhi[0] = new TH2D("fHistEtaPhi", "", 16,-0.8,0.8, 100, 0, TMath::TwoPi());
  fHist2DEtaPhi[1] = new TH2D("fHistEtaPhi_afterfNUA", "", 16,-0.8,0.8, 100, 0, TMath::TwoPi());
  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistEta);
  fOutputList->Add(fHistNhits);
  fOutputList->Add(fHist2DPDedx);
  fOutputList->Add(fHistDcaXY);
  fOutputList->Add(fHistDcaZ);
  for (int i = 0; i < 2; i++) fOutputList->Add(fHistPhi[i]);
  for (int i = 0; i < 2; i++) fOutputList->Add(fHist2DEtaPhi[i]);

  //Plane QA
  char chCalibStep[5];
  for (int i = 0; i < 2; i++) {
    if (i==0) sprintf(chCalibStep,"GE");
    if (i==1) sprintf(chCalibStep,"RC");
    if (IsQAVZERO) {
      fProfileV0CQxCent[i] = new TProfile(Form("fProfileV0CQxCent%s",chCalibStep), "", 100, 0, 100);
      fProfileV0CQyCent[i] = new TProfile(Form("fProfileV0CQyCent%s",chCalibStep), "", 100, 0, 100);
      fProfileV0CQxVtx[i]  = new TProfile(Form("fProfileV0CQxVtx%s",chCalibStep), "", 20, -10, 10);
      fProfileV0CQyVtx[i]  = new TProfile(Form("fProfileV0CQyVtx%s",chCalibStep), "", 20, -10, 10);
      fHist2DCalibPsi2V0CCent[i] = new TH2D(Form("fHist2DCalibPsi2V0CCent%s",chCalibStep), "", 20, 0, 100, 100, 0, TMath::TwoPi());
      fOutputList->Add(fProfileV0CQxCent[i]);
      fOutputList->Add(fProfileV0CQyCent[i]);
      fOutputList->Add(fProfileV0CQxVtx[i]);
      fOutputList->Add(fProfileV0CQyVtx[i]);
      fOutputList->Add(fHist2DCalibPsi2V0CCent[i]);

      fProfileV0AQxCent[i] = new TProfile(Form("fProfileV0AQxCent%s",chCalibStep), "", 100, 0, 100);
      fProfileV0AQyCent[i] = new TProfile(Form("fProfileV0AQyCent%s",chCalibStep), "", 100, 0, 100);
      fProfileV0AQxVtx[i]  = new TProfile(Form("fProfileV0AQxVtx%s",chCalibStep), "", 20, -10, 10);
      fProfileV0AQyVtx[i]  = new TProfile(Form("fProfileV0AQyVtx%s",chCalibStep), "", 20, -10, 10);
      fHist2DCalibPsi2V0ACent[i] = new TH2D(Form("fHist2DCalibPsi2V0ACent%s",chCalibStep), "", 20, 0, 100, 100, 0, TMath::TwoPi());
      fOutputList->Add(fProfileV0AQxCent[i]);
      fOutputList->Add(fProfileV0AQyCent[i]);
      fOutputList->Add(fProfileV0AQxVtx[i]);
      fOutputList->Add(fProfileV0AQyVtx[i]);
      fOutputList->Add(fHist2DCalibPsi2V0ACent[i]);
    }

    if (IsQAZDC) {
      fProfileZNCQxCent[i] = new TProfile(Form("fProfileZNCQxCent%s",chCalibStep), "", 100, 0, 100);
      fProfileZNCQyCent[i] = new TProfile(Form("fProfileZNCQyCent%s",chCalibStep), "", 100, 0, 100);
      fHist2DCalibPsi1ZNCCent[i] = new TH2D(Form("fHist2DCalibPsi1ZNCCent%s",chCalibStep), "", 20, 0, 100, 100, 0, TMath::TwoPi());
      fOutputList->Add(fProfileZNCQxCent[i]);
      fOutputList->Add(fProfileZNCQyCent[i]);
      fOutputList->Add(fHist2DCalibPsi1ZNCCent[i]);

      fProfileZNAQxCent[i] = new TProfile(Form("fProfileZNAQxCent%s",chCalibStep), "", 100, 0, 100);
      fProfileZNAQyCent[i] = new TProfile(Form("fProfileZNAQyCent%s",chCalibStep), "", 100, 0, 100);
      fHist2DCalibPsi1ZNACent[i] = new TH2D(Form("fHist2DCalibPsi1ZNACent%s",chCalibStep), "", 20, 0, 100, 100, 0, TMath::TwoPi());
      fOutputList->Add(fProfileZNAQxCent[i]);
      fOutputList->Add(fProfileZNAQyCent[i]);
      fOutputList->Add(fHist2DCalibPsi1ZNACent[i]);

      fProfileZDCQxAQxCCent[i] = new TProfile(Form("fProfileZDCQxAQxCCent%s",chCalibStep), "", 100, 0, 100);
      fProfileZDCQxAQyCCent[i] = new TProfile(Form("fProfileZDCQxAQyCCent%s",chCalibStep), "", 100, 0, 100);
      fProfileZDCQyAQxCCent[i] = new TProfile(Form("fProfileZDCQyAQxCCent%s",chCalibStep), "", 100, 0, 100);
      fProfileZDCQyAQyCCent[i] = new TProfile(Form("fProfileZDCQyAQyCCent%s",chCalibStep), "", 100, 0, 100);
      fOutputList->Add(fProfileZDCQxAQxCCent[i]);
      fOutputList->Add(fProfileZDCQxAQyCCent[i]);
      fOutputList->Add(fProfileZDCQyAQxCCent[i]);
      fOutputList->Add(fProfileZDCQyAQyCCent[i]);
    }
  }

  if (IsQAZDC) {
    fHist2DCalibPsi1ZNCCent[2] = new TH2D("fHist2DCalibPsi1ZNCCentSF", "", 20, 0, 100, 100, 0, TMath::TwoPi());
    fHist2DCalibPsi1ZNACent[2] = new TH2D("fHist2DCalibPsi1ZNACentSF", "", 20, 0, 100, 100, 0, TMath::TwoPi());
    fOutputList->Add(fHist2DCalibPsi1ZNCCent[2]);
    fOutputList->Add(fHist2DCalibPsi1ZNACent[2]);

    fProfileZNCTowerMeanEnegry[0] = new TProfile("fProfileZNCTowerMeanEnegryRW","",5,0,5);
    fProfileZNCTowerMeanEnegry[1] = new TProfile("fProfileZNCTowerMeanEnegryGE","",5,0,5);
    fProfileZNATowerMeanEnegry[0] = new TProfile("fProfileZNATowerMeanEnegryRW","",5,0,5);
    fProfileZNATowerMeanEnegry[1] = new TProfile("fProfileZNATowerMeanEnegryGE","",5,0,5);
    fOutputList->Add(fProfileZNCTowerMeanEnegry[0]);
    fOutputList->Add(fProfileZNCTowerMeanEnegry[1]);
    fOutputList->Add(fProfileZNATowerMeanEnegry[0]);
    fOutputList->Add(fProfileZNATowerMeanEnegry[1]);
  }

  //V0s QA
  fHistV0Pt = new TH1D("hV0Pt","", 200, 0., 20.);
  fHistV0Eta = new TH1D("hV0Eta","", 200, -10., 10.);
  fHistV0DcatoPrimVertex = new TH1D("hV0DcaToPrimVertex","",200, 0., 20.);
  fHistV0CPA = new TH1D("hV0CPA","", 1000, 0.9, 1.);
  fHistV0DecayLength = new TH1D("hV0DecayLength","",500,0,500.);
  fOutputList->Add(fHistV0Pt);
  fOutputList->Add(fHistV0Eta);
  fOutputList->Add(fHistV0DcatoPrimVertex);
  fOutputList->Add(fHistV0CPA);
  fOutputList->Add(fHistV0DecayLength);

  //Lambda QA
  char chCut[5];
  for (int i=0; i<2; i++) {
    ///// Case 0 = before cut, case 1 = afterCut.
    if (i==0) sprintf(chCut,"Bf");
    if (i==1) sprintf(chCut,"Af");
    /// Lambdas:
    fHistLambdaPt[i] = new TH1D(Form("hLambdaPt_%sMassCut",chCut),"", 200, 0., 20.);
    fHistLambdaEta[i] = new TH1D(Form("hLambdaEta_%sMassCut",chCut),"",200, -10., 10.);
    fHistLambdaPhi[i] = new TH1D(Form("hLambdaPhi_%sMassCut",chCut),"", 360, 0., TMath::TwoPi());
    fHistLambdaDcaToPrimVertex[i] = new TH1D(Form("hLambdaDcaToPrimVertex_%sMassCut",chCut),"",200, 0., 20.);
    fHistLambdaCPA[i] = new TH1D(Form("hLambdaCPA_%sMassCut",chCut),"",200, 0.9, 1.);
    fHistLambdaDecayLength[i] = new TH1D(Form("hLambdaDecayLength_%sMassCut",chCut),"", 250, 0., 500.);
    fHistLambdaMass[i] = new TH1D(Form("hLambdaMass_%sMassCut",chCut),"",1000,1.,1.25); //  Current bin size = 0.00025
    fProfileLambdaMassVsPt[i] = new TProfile(Form("pLambdaMassVsPt_%sMassCut",chCut),"",200,0,20);
    fOutputList->Add(fHistLambdaPt[i]);
    fOutputList->Add(fHistLambdaEta[i]);
    fOutputList->Add(fHistLambdaPhi[i]);
    fOutputList->Add(fHistLambdaDcaToPrimVertex[i]);
    fOutputList->Add(fHistLambdaCPA[i]);
    fOutputList->Add(fHistLambdaDecayLength[i]);
    fOutputList->Add(fHistLambdaMass[i]);
    fOutputList->Add(fProfileLambdaMassVsPt[i]);

    // AntiLambdas
    fHistAntiLambdaPt[i] = new TH1D(Form("hAntiLambdaPt_%sMassCut",chCut),"", 200, 0., 20.);
    fHistAntiLambdaEta[i] = new TH1D(Form("hAntiLambdaEta_%sMassCut",chCut),"",200, -10., 10.);
    fHistAntiLambdaPhi[i] = new TH1D(Form("hAntiLambdaPhi_%sMassCut",chCut),"", 360, 0, TMath::TwoPi());
    fHistAntiLambdaDcaToPrimVertex[i] = new TH1D(Form("hAntiLambdaDcaToPrimVertex_%sMassCut",chCut),"",200, 0., 20.);
    fHistAntiLambdaCPA[i] = new TH1D(Form("hAntiLambdaCPA_%sMassCut",chCut),"",200, 0.9, 1.);
    fHistAntiLambdaDecayLength[i] = new TH1D(Form("hAntiLambdaDecayLength_%sMassCut",chCut),"", 250, 0., 500.);
    fHistAntiLambdaMass[i] = new TH1D(Form("hAntiLambdaMass_%sMassCut",chCut),"",1000,1.,1.25); // Current bin size = 0.00025
    fProfileAntiLambdaMassVsPt[i] = new TProfile(Form("pAntiLambdaMassVsPt_%sMassCut",chCut),"",200,0,20);
    fOutputList->Add(fHistAntiLambdaPt[i]);
    fOutputList->Add(fHistAntiLambdaEta[i]);
    fOutputList->Add(fHistAntiLambdaPhi[i]);
    fOutputList->Add(fHistAntiLambdaDcaToPrimVertex[i]);
    fOutputList->Add(fHistAntiLambdaCPA[i]);
    fOutputList->Add(fHistAntiLambdaDecayLength[i]);
    fOutputList->Add(fHistAntiLambdaMass[i]);
    fOutputList->Add(fProfileAntiLambdaMassVsPt[i]);
  }

  //Flow QA
  char chPlaneType[5];
  for (int i = 0; i < 5; i++) {
    if (i==0) sprintf(chPlaneType,"TPC");
    if (i==1) sprintf(chPlaneType,"V0C");
    if (i==2) sprintf(chPlaneType,"V0A");
    if (i==3) sprintf(chPlaneType,"ZNC");
    if (i==4) sprintf(chPlaneType,"ZNA");
    fProfile2DRawFlowCentPthPos[i] = new TProfile2D(Form("fProfile2DRawFlowCentPthPos_%s",chPlaneType),"",10,0,100,25,0,5);
    fProfile2DRawFlowCentPthNeg[i] = new TProfile2D(Form("fProfile2DRawFlowCentPthNeg_%s",chPlaneType),"",10,0,100,25,0,5);
    fProfile2DRawFlowCentPtProton[i] = new TProfile2D(Form("fProfile2DRawFlowCentPtProton_%s",chPlaneType),"",10,0,100,25,0,5);
    fProfile2DRawFlowCentPtAntiProton[i] = new TProfile2D(Form("fProfile2DRawFlowCentPtAntiProton_%s",chPlaneType),"",10,0,100,25,0,5);
    fProfile2DRawFlowCentPtLambda[i] = new TProfile2D(Form("fProfile2DRawFlowCentPtLambda_%s",chPlaneType),"",10,0,100,25,0,5);
    fProfile2DRawFlowCentPtAntiLambda[i] = new TProfile2D(Form("fProfile2DRawFlowCentPtAntiLambda_%s",chPlaneType),"",10,0,100,25,0,5);
    fOutputList->Add(fProfile2DRawFlowCentPthPos[i]);
    fOutputList->Add(fProfile2DRawFlowCentPthNeg[i]);
    fOutputList->Add(fProfile2DRawFlowCentPtProton[i]);
    fOutputList->Add(fProfile2DRawFlowCentPtAntiProton[i]);
    fOutputList->Add(fProfile2DRawFlowCentPtLambda[i]);
    fOutputList->Add(fProfile2DRawFlowCentPtAntiLambda[i]);
  }


  ////////////////////////
  // Results
  ////////////////////////
  // Plane we used
  fHist2DPsi2TPCPosCent = new TH2D("fHist2DPsi2TPCPosCent","",20,0,100,100,0,TMath::TwoPi());
  fHist2DPsi2TPCNegCent = new TH2D("fHist2DPsi2TPCNegCent","",20,0,100,100,0,TMath::TwoPi());
  fHist2DPsi2V0CCent    = new TH2D("fHist2DPsi2V0CCent","",20,0,100,100,0,TMath::TwoPi());
  fHist2DPsi2V0ACent    = new TH2D("fHist2DPsi2V0ACent","",20,0,100,100,0,TMath::TwoPi());
  fHist2DPsi1ZNCCent    = new TH2D("fHist2DPsi1ZNCCent","",20,0,100,100,0,TMath::TwoPi());
  fHist2DPsi1ZNACent    = new TH2D("fHist2DPsi1ZNACent","",20,0,100,100,0,TMath::TwoPi());
  fOutputList->Add(fHist2DPsi2TPCPosCent);
  fOutputList->Add(fHist2DPsi2TPCNegCent);
  fOutputList->Add(fHist2DPsi2V0CCent);
  fOutputList->Add(fHist2DPsi2V0ACent);
  fOutputList->Add(fHist2DPsi1ZNCCent);
  fOutputList->Add(fHist2DPsi1ZNACent);

  //Resolution
  ///SideC-SideA Event Plane Correlations for Resolution Estimation
  fProfileTPCPsi2Correlation = new TProfile("fProfileTPCPsi2Correlation","TPCPos-TPCNeg Psi2 Res. vs Cent; Cent; Resolution",100,0,100.);
  fProfileV0MPsi2Correlation = new TProfile("fProfileV0MPsi2Correlation","V0C-V0A Psi2 Res. vs Cent; Cent; Resolution",100,0,100.);
  fProfileZDCPsi1Correlation = new TProfile("fProfileZDCPsi1Correlation","ZNC-ZNA Psi1 Res. vs Cent; Cent; Resolution",100,0,100.);
  fProfileZDCPsi2Correlation = new TProfile("fProfileZDCPsi2Correlation","ZNC-ZNA Psi2 Res. vs Cent; Cent; Resolution",100,0,100.);
  fOutputList->Add(fProfileTPCPsi2Correlation);
  fOutputList->Add(fProfileV0MPsi2Correlation);
  fOutputList->Add(fProfileZDCPsi1Correlation);
  fOutputList->Add(fProfileZDCPsi2Correlation);
  ///V0X-TPC Event Plane Correlations for Resolution:
  fProfileV0CTPCPosPsi2Correlation = new TProfile("fProfileV0CTPCPosPsi2Correlation","V0C-TPCPos Psi2; Cent; Resolution",100,0,100.);
  fProfileV0ATPCPosPsi2Correlation = new TProfile("fProfileV0ATPCPosPsi2Correlation","V0A-TPCPos Psi2; Cent; Resolution",100,0,100.);
  fProfileV0CTPCNegPsi2Correlation = new TProfile("fProfileV0CTPCNegPsi2Correlation","V0C-TPCNeg Psi2; Cent; Resolution",100,0,100.);
  fProfileV0ATPCNegPsi2Correlation = new TProfile("fProfileV0ATPCNegPsi2Correlation","V0A-TPCNeg Psi2; Cent; Resolution",100,0,100.);
  fOutputList->Add(fProfileV0CTPCPosPsi2Correlation);
  fOutputList->Add(fProfileV0ATPCPosPsi2Correlation);
  fOutputList->Add(fProfileV0CTPCNegPsi2Correlation);
  fOutputList->Add(fProfileV0ATPCNegPsi2Correlation);

  ///Lambda-X correlators
  ///Delta Correlators:
  ///Lambda - X
  fProfileDelta_Lambda_hPos = new TProfile("fProfileDelta_Lambda_hPos", "", 20, 0., 100.);
  fProfileDelta_Lambda_hNeg = new TProfile("fProfileDelta_Lambda_hNeg", "", 20, 0., 100.);
  fProfileDelta_Lambda_Proton = new TProfile("fProfileDelta_Lambda_Proton", "", 20, 0., 100.);
  fProfileDelta_Lambda_AntiProton = new TProfile("fProfileDelta_Lambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileDelta_Lambda_hPos);
  fOutputList->Add(fProfileDelta_Lambda_hNeg);
  fOutputList->Add(fProfileDelta_Lambda_Proton);
  fOutputList->Add(fProfileDelta_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileDelta_AntiLambda_hPos = new TProfile("fProfileDelta_AntiLambda_hPos", "", 20, 0., 100.);
  fProfileDelta_AntiLambda_hNeg = new TProfile("fProfileDelta_AntiLambda_hNeg", "", 20, 0., 100.);
  fProfileDelta_AntiLambda_Proton = new TProfile("fProfileDelta_AntiLambda_Proton", "", 20, 0., 100.);
  fProfileDelta_AntiLambda_AntiProton = new TProfile("fProfileDelta_AntiLambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileDelta_AntiLambda_hPos);
  fOutputList->Add(fProfileDelta_AntiLambda_hNeg);
  fOutputList->Add(fProfileDelta_AntiLambda_Proton);
  fOutputList->Add(fProfileDelta_AntiLambda_AntiProton);

  ///Gamma Correlators:
  //TPC Plane
  ///Lambda - X
  fProfileGammaTPC_Lambda_hPos = new TProfile("fProfileGammaTPC_Lambda_hPos", "", 20, 0., 100.);
  fProfileGammaTPC_Lambda_hNeg = new TProfile("fProfileGammaTPC_Lambda_hNeg", "", 20, 0., 100.);
  fProfileGammaTPC_Lambda_Proton = new TProfile("fProfileGammaTPC_Lambda_Proton", "", 20, 0., 100.);
  fProfileGammaTPC_Lambda_AntiProton = new TProfile("fProfileGammaTPC_Lambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileGammaTPC_Lambda_hPos);
  fOutputList->Add(fProfileGammaTPC_Lambda_hNeg);
  fOutputList->Add(fProfileGammaTPC_Lambda_Proton);
  fOutputList->Add(fProfileGammaTPC_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaTPC_AntiLambda_hPos = new TProfile("fProfileGammaTPC_AntiLambda_hPos", "", 20, 0., 100.);
  fProfileGammaTPC_AntiLambda_hNeg = new TProfile("fProfileGammaTPC_AntiLambda_hNeg", "", 20, 0., 100.);
  fProfileGammaTPC_AntiLambda_Proton = new TProfile("fProfileGammaTPC_AntiLambda_Proton", "", 20, 0., 100.);
  fProfileGammaTPC_AntiLambda_AntiProton = new TProfile("fProfileGammaTPC_AntiLambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileGammaTPC_AntiLambda_hPos);
  fOutputList->Add(fProfileGammaTPC_AntiLambda_hNeg);
  fOutputList->Add(fProfileGammaTPC_AntiLambda_Proton);
  fOutputList->Add(fProfileGammaTPC_AntiLambda_AntiProton);

  //V0C Plane
  ///Lambda - X
  fProfileGammaV0C_Lambda_hPos = new TProfile("fProfileGammaV0C_Lambda_hPos", "", 20, 0., 100.);
  fProfileGammaV0C_Lambda_hNeg = new TProfile("fProfileGammaV0C_Lambda_hNeg", "", 20, 0., 100.);
  fProfileGammaV0C_Lambda_Proton = new TProfile("fProfileGammaV0C_Lambda_Proton", "", 20, 0., 100.);
  fProfileGammaV0C_Lambda_AntiProton = new TProfile("fProfileGammaV0C_Lambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileGammaV0C_Lambda_hPos);
  fOutputList->Add(fProfileGammaV0C_Lambda_hNeg);
  fOutputList->Add(fProfileGammaV0C_Lambda_Proton);
  fOutputList->Add(fProfileGammaV0C_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaV0C_AntiLambda_hPos = new TProfile("fProfileGammaV0C_AntiLambda_hPos", "", 20, 0., 100.);
  fProfileGammaV0C_AntiLambda_hNeg = new TProfile("fProfileGammaV0C_AntiLambda_hNeg", "", 20, 0., 100.);
  fProfileGammaV0C_AntiLambda_Proton = new TProfile("fProfileGammaV0C_AntiLambda_Proton", "", 20, 0., 100.);
  fProfileGammaV0C_AntiLambda_AntiProton = new TProfile("fProfileGammaV0C_AntiLambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileGammaV0C_AntiLambda_hPos);
  fOutputList->Add(fProfileGammaV0C_AntiLambda_hNeg);
  fOutputList->Add(fProfileGammaV0C_AntiLambda_Proton);
  fOutputList->Add(fProfileGammaV0C_AntiLambda_AntiProton);

  //V0A Plane
  ///Lambda - X
  fProfileGammaV0A_Lambda_hPos = new TProfile("fProfileGammaV0A_Lambda_hPos", "", 20, 0., 100.);
  fProfileGammaV0A_Lambda_hNeg = new TProfile("fProfileGammaV0A_Lambda_hNeg", "", 20, 0., 100.);
  fProfileGammaV0A_Lambda_Proton = new TProfile("fProfileGammaV0A_Lambda_Proton", "", 20, 0., 100.);
  fProfileGammaV0A_Lambda_AntiProton = new TProfile("fProfileGammaV0A_Lambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileGammaV0A_Lambda_hPos);
  fOutputList->Add(fProfileGammaV0A_Lambda_hNeg);
  fOutputList->Add(fProfileGammaV0A_Lambda_Proton);
  fOutputList->Add(fProfileGammaV0A_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaV0A_AntiLambda_hPos = new TProfile("fProfileGammaV0A_AntiLambda_hPos", "", 20, 0., 100.);
  fProfileGammaV0A_AntiLambda_hNeg = new TProfile("fProfileGammaV0A_AntiLambda_hNeg", "", 20, 0., 100.);
  fProfileGammaV0A_AntiLambda_Proton = new TProfile("fProfileGammaV0A_AntiLambda_Proton", "", 20, 0., 100.);
  fProfileGammaV0A_AntiLambda_AntiProton = new TProfile("fProfileGammaV0A_AntiLambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileGammaV0A_AntiLambda_hPos);
  fOutputList->Add(fProfileGammaV0A_AntiLambda_hNeg);
  fOutputList->Add(fProfileGammaV0A_AntiLambda_Proton);
  fOutputList->Add(fProfileGammaV0A_AntiLambda_AntiProton);

  //ZNC Plane
  ///Lambda - X
  fProfileGammaZNC_Lambda_hPos = new TProfile("fProfileGammaZNC_Lambda_hPos", "", 20, 0., 100.);
  fProfileGammaZNC_Lambda_hNeg = new TProfile("fProfileGammaZNC_Lambda_hNeg", "", 20, 0., 100.);
  fProfileGammaZNC_Lambda_Proton = new TProfile("fProfileGammaZNC_Lambda_Proton", "", 20, 0., 100.);
  fProfileGammaZNC_Lambda_AntiProton = new TProfile("fProfileGammaZNC_Lambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileGammaZNC_Lambda_hPos);
  fOutputList->Add(fProfileGammaZNC_Lambda_hNeg);
  fOutputList->Add(fProfileGammaZNC_Lambda_Proton);
  fOutputList->Add(fProfileGammaZNC_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaZNC_AntiLambda_hPos = new TProfile("fProfileGammaZNC_AntiLambda_hPos", "", 20, 0., 100.);
  fProfileGammaZNC_AntiLambda_hNeg = new TProfile("fProfileGammaZNC_AntiLambda_hNeg", "", 20, 0., 100.);
  fProfileGammaZNC_AntiLambda_Proton = new TProfile("fProfileGammaZNC_AntiLambda_Proton", "", 20, 0., 100.);
  fProfileGammaZNC_AntiLambda_AntiProton = new TProfile("fProfileGammaZNC_AntiLambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileGammaZNC_AntiLambda_hPos);
  fOutputList->Add(fProfileGammaZNC_AntiLambda_hNeg);
  fOutputList->Add(fProfileGammaZNC_AntiLambda_Proton);
  fOutputList->Add(fProfileGammaZNC_AntiLambda_AntiProton);

  //ZNA Plane
  ///Lambda - X
  fProfileGammaZNA_Lambda_hPos = new TProfile("fProfileGammaZNA_Lambda_hPos", "", 20, 0., 100.);
  fProfileGammaZNA_Lambda_hNeg = new TProfile("fProfileGammaZNA_Lambda_hNeg", "", 20, 0., 100.);
  fProfileGammaZNA_Lambda_Proton = new TProfile("fProfileGammaZNA_Lambda_Proton", "", 20, 0., 100.);
  fProfileGammaZNA_Lambda_AntiProton = new TProfile("fProfileGammaZNA_Lambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileGammaZNA_Lambda_hPos);
  fOutputList->Add(fProfileGammaZNA_Lambda_hNeg);
  fOutputList->Add(fProfileGammaZNA_Lambda_Proton);
  fOutputList->Add(fProfileGammaZNA_Lambda_AntiProton);

  ///AntiLambda - X
  fProfileGammaZNA_AntiLambda_hPos = new TProfile("fProfileGammaZNA_AntiLambda_hPos", "", 20, 0., 100.);
  fProfileGammaZNA_AntiLambda_hNeg = new TProfile("fProfileGammaZNA_AntiLambda_hNeg", "", 20, 0., 100.);
  fProfileGammaZNA_AntiLambda_Proton = new TProfile("fProfileGammaZNA_AntiLambda_Proton", "", 20, 0., 100.);
  fProfileGammaZNA_AntiLambda_AntiProton = new TProfile("fProfileGammaZNA_AntiLambda_AntiProton", "", 20, 0., 100.);
  fOutputList->Add(fProfileGammaZNA_AntiLambda_hPos);
  fOutputList->Add(fProfileGammaZNA_AntiLambda_hNeg);
  fOutputList->Add(fProfileGammaZNA_AntiLambda_Proton);
  fOutputList->Add(fProfileGammaZNA_AntiLambda_AntiProton);

  PostData(1,fOutputList);
  if (fDebug) Printf("UserCreateOutputObjects() Post Data Success!");
}

//------------------------------------------------

void AliAnalysisTaskLambdaProtonCVE::UserExec(Option_t *)
{
  if (fDebug) Printf("===============================We are in UserExec!================================");
  fEvtCount->Fill(1);
  //----------------------------
  // Handle
  //----------------------------
  AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
  if (!manager) {
    AliError(Form("%s: Could not get Analysis Manager", GetName()));
  } else fEvtCount->Fill(20);

  AliAODInputHandler* handler = (AliAODInputHandler*)manager->GetInputEventHandler();
  if (!handler) {
    AliError(Form("%s: Could not get Input Handler", GetName()));
  } else fEvtCount->Fill(21);

  fAOD = dynamic_cast <AliAODEvent*> (InputEvent());
  if (!fAOD) {
    AliError(Form("%s: Could not get AOD event", GetName()));
  } else fEvtCount->Fill(22);

  fPIDResponse = handler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError(Form("%s: Could not get PIDResponse", GetName()));
  } else fEvtCount->Fill(23);

  fUtils = new AliAnalysisUtils();
  if (!fUtils) {
    AliError(Form("%s: Could not get AliAnalysisUtils", GetName()));
  } else fEvtCount->Fill(24);

  if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    if (!fMultSel) {
      AliError(Form("%s: Could not get AliMultSelection", GetName()));
    } else fEvtCount->Fill(25);
    if (!manager || !handler || !fAOD || !fPIDResponse || !fUtils || !fMultSel) return;
  }

  if (!manager || !handler || !fAOD || !fPIDResponse || !fUtils) return;
  fEvtCount->Fill(2);
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
  fEvtCount->Fill(3);
  if (fDebug) Printf("trigger done!");

  //----------------------------
  // Run Number
  //----------------------------
  fRunNum = fAOD->GetRunNumber();
  if (fRunNum != fOldRunNum) {
     // Load the run dependent calibration hist
      if (!LoadCalibHistForThisRun()) return;
      fRunNumBin = runNumList->at(fRunNum);
      fOldRunNum = fRunNum;
      if (fRunNumBin < 0) return;
  }
  fHistRunNumBin->Fill(fRunNumBin);
  fEvtCount->Fill(4);
  if (fDebug) Printf("run nummbr done!");

  //----------------------------
  // Vertex
  //----------------------------
  AliAODVertex* fVtx = fAOD->GetPrimaryVertex();
  fVtx -> GetXYZ(fVertex);
  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  double vx = fVertex[0];
  double vy = fVertex[1];
  double vz = fVertex[2];
  if (fabs(fVertex[0])<1e-6 || fabs(fVertex[1])<1e-6 || fabs(fVertex[2])<1e-6) return;
  double dz = vz - fAOD->GetPrimaryVertexSPD()->GetZ();
  if (fabs(vz) > fVzCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors()<1) return;
  fHistVz[0]->Fill(vz);
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
  fHistVz[1]->Fill(vz);
  for (int i = 0; i < 20; ++i) {
      if (vz > -10+i*1 && vz < -10+(i+1)*1) {fVzBin = i; break;}
  }
  if (fVzBin<-990) return;
  fEvtCount->Fill(5);
  if (fDebug) Printf("vertex done!");

  //----------------------------
  // Centrality
  //----------------------------
  double centV0M = -1, centTRK = -1, centSPD0 = -1, centSPD1 = -1, centV0A = -1;
  if (fPeriod.EqualTo("LHC10h")) {
    centV0M  = fAOD->GetCentrality()->GetCentralityPercentile("V0M");
    centTRK  = fAOD->GetCentrality()->GetCentralityPercentile("TRK");
    centSPD0 = fAOD->GetCentrality()->GetCentralityPercentile("CL0");
    centSPD1 = fAOD->GetCentrality()->GetCentralityPercentile("CL1");
    centV0A  = fAOD->GetCentrality()->GetCentralityPercentile("V0A");
  } else if (fPeriod.EqualTo("LHC15o") || fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
    centV0M  = fMultSel->GetMultiplicityPercentile("V0M");
    centTRK  = fMultSel->GetMultiplicityPercentile("TRK");
    centSPD0 = fMultSel->GetMultiplicityPercentile("CL0");
    centSPD1 = fMultSel->GetMultiplicityPercentile("CL1");
  }
  //we use centV0M as the default centrality
  fCent = centV0M;
  fHist2DCentQA[0]->Fill(centV0M,centSPD1);
  fHist2DCentQA[2]->Fill(centV0M,centTRK);
  fHist2DCentQA[4]->Fill(centV0M,centSPD0);
  fHist2DCentQA[6]->Fill(centSPD1,centSPD0);
  if (fabs(fCent-centSPD1)>fCentDiffCut) return;
  fHist2DCentQA[1]->Fill(centV0M,centSPD1);
  fHist2DCentQA[3]->Fill(centV0M,centTRK);
  fHist2DCentQA[5]->Fill(centV0M,centSPD0);
  fHist2DCentQA[7]->Fill(centSPD1,centSPD0);
  if (fCent < 0 || fCent >= 80) return;
  // cent bin
  fCentBin = (int)fCent/10;
  fHistCent[0]->Fill(fCent);
  fEvtCount->Fill(6);
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
  fEvtCount->Fill(7);
  if (fDebug) Printf("pile-up done!");

  //----------------------------
  // VZERO Plane
  //----------------------------
  if (!GetVZEROPlane()) return;
  fEvtCount->Fill(8);
  if (fDebug) Printf("Get VZERO Plane done!");
  //----------------------------
  // ZDC Plane
  //----------------------------
  if (fPeriod.EqualTo("LHC10h")) if (!GetZDCPlane()) return;
  if (fPeriod.EqualTo("LHC18q")
  || fPeriod.EqualTo("LHC18r"))  if (!GetZDCPlaneLsFit()) return;

  fEvtCount->Fill(9);
  if (fDebug) Printf("Get ZDC Plane done!");
  //----------------------------
  // Loop Tracks / Fill Vectors
  //----------------------------
  //Reset vectors
  ResetVectors();
  fEvtCount->Fill(10);
  if (!LoopTracks()) return;
  fEvtCount->Fill(11);
  if (fDebug) Printf("Loop Tracks done!");
  //----------------------------
  // TPC Plane
  //----------------------------
  if (!GetTPCPlane()) return;
  fEvtCount->Fill(12);
  if (fDebug) Printf("Get TPC Plane done!");
  //----------------------------
  // Fill Resolution
  //----------------------------
  fHist2DPsi2TPCPosCent -> Fill(fCent,    fPsi2TPCPos);
  fHist2DPsi2TPCNegCent -> Fill(fCent,    fPsi2TPCNeg);
  fHist2DPsi2V0CCent    -> Fill(centSPD1, fPsi2V0C);
  fHist2DPsi2V0ACent    -> Fill(centSPD1, fPsi2V0A);
  fHist2DPsi1ZNCCent    -> Fill(fCent,    fPsi1ZNC);
  fHist2DPsi1ZNACent    -> Fill(fCent,    fPsi1ZNA);

  fProfileTPCPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2TPCNeg - fPsi2TPCPos)));
  fProfileV0MPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0C - fPsi2V0A)));
  fProfileZDCPsi1Correlation -> Fill(fCent, TMath::Cos(1*(fPsi1ZNC - fPsi1ZNA)));
  fProfileZDCPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi1ZNC - fPsi1ZNA)));

  fProfileV0CTPCPosPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0C - fPsi2TPCPos)));
  fProfileV0ATPCPosPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0A - fPsi2TPCPos)));
  fProfileV0CTPCNegPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0C - fPsi2TPCNeg)));
  fProfileV0ATPCNegPsi2Correlation -> Fill(fCent, TMath::Cos(2*(fPsi2V0A - fPsi2TPCNeg)));
  fEvtCount->Fill(13);
  //----------------------------
  // Get Lambda Vector
  //----------------------------
  if (!LoopV0s()) return;
  fEvtCount->Fill(14);
  if (fDebug) Printf("Get Lambda Vector done!");
  //----------------------------
  // Pair
  //----------------------------
  if (!PairLambda()) return;
  fEvtCount->Fill(15);
  if (fDebug) Printf("Pair done!");
  //------------------
  // Post output data.
  //------------------
  if (fDebug) Printf("analysis done!");
  PostData(1,fOutputList);
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::GetVZEROPlane()
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
    AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
    double centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
    int iCentSPD = (int)centrCL1;
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
  if (IsQAVZERO) {
    double centSPD  = -999.;
    if (fPeriod.EqualTo("LHC15o")||fPeriod.EqualTo("LHC18q")||fPeriod.EqualTo("LHC18r")) {
      AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
      centSPD = fMultSelection->GetMultiplicityPercentile("CL1");
    }
    else if (fPeriod.EqualTo("LHC10h") ) {
      centSPD =  fAOD->GetCentrality()->GetCentralityPercentile("CL1");
    }
    //V0C
    fProfileV0CQxCent[0]->Fill(centSPD, qxGE[1]);
    fProfileV0CQyCent[0]->Fill(centSPD, qyGE[1]);
    fProfileV0CQxVtx[0] ->Fill(fVertex[2], qxGE[1]);
    fProfileV0CQyVtx[0] ->Fill(fVertex[2], qyGE[1]);
    fHist2DCalibPsi2V0CCent[0]->Fill(centSPD, psi2GE[1]);

    fProfileV0CQxCent[1]->Fill(centSPD, qxRC[1]);
    fProfileV0CQyCent[1]->Fill(centSPD, qyRC[1]);
    fProfileV0CQxVtx[1]->Fill(fVertex[2], qxRC[1]);
    fProfileV0CQyVtx[1]->Fill(fVertex[2], qyRC[1]);
    fHist2DCalibPsi2V0CCent[1]->Fill(centSPD, psi2GE[1]);
    //V0A
    fProfileV0AQxCent[0]->Fill(centSPD, qxGE[2]);
    fProfileV0AQyCent[0]->Fill(centSPD, qyGE[2]);
    fProfileV0AQxVtx[0]->Fill(fVertex[2], qxGE[2]);
    fProfileV0AQyVtx[0]->Fill(fVertex[2], qyGE[2]);
    fHist2DCalibPsi2V0ACent[0]->Fill(centSPD, psi2GE[2]);

    fProfileV0AQxCent[1]->Fill(centSPD, qxRC[2]);
    fProfileV0AQyCent[1]->Fill(centSPD, qyRC[2]);
    fProfileV0AQxVtx[1]->Fill(fVertex[2], qxRC[2]);
    fProfileV0AQyVtx[1]->Fill(fVertex[2], qyRC[2]);
    fHist2DCalibPsi2V0ACent[1]->Fill(centSPD, psi2RC[2]);
  }
  fPsi2V0C = psi2RC[1];
  fPsi2V0A = psi2RC[2];
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::GetZDCPlane()
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

  if(IsZDCCalibOn) {
    for (int iTower = 0; iTower < 5; ++iTower) {
      fProfileZNCTowerMeanEnegry[0]->Fill(iTower+0.5, towerEnegryZNC[iTower]);
      fProfileZNATowerMeanEnegry[0]->Fill(iTower+0.5, towerEnegryZNA[iTower]);
    }
  }

  //Loop Over ZDC Towers
  //Gain Equalization
  for (int iTower = 0; iTower < 5; ++iTower) {
    towerEnegryZNCGE[iTower] = towerEnegryZNC[iTower] / (fProfileForZNCGE->GetBinContent(iTower + 1)) * (fProfileForZNCGE->GetBinContent(2));//here to select the ref Tower
    towerEnegryZNAGE[iTower] = towerEnegryZNA[iTower] / (fProfileForZNAGE->GetBinContent(iTower + 1)) * (fProfileForZNAGE->GetBinContent(2));
  }

  if(IsZDCCalibOn) {
    for (int iTower = 0; iTower < 5; ++iTower) {
      fProfileZNCTowerMeanEnegry[1]->Fill(iTower+0.5, towerEnegryZNCGE[iTower]);
      fProfileZNATowerMeanEnegry[1]->Fill(iTower+0.5, towerEnegryZNAGE[iTower]);
    }
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

  if (IsQAZDC) {
    //ZNC
    fProfileZNCQxCent[0] -> Fill(fCent, QxC);
    fProfileZNCQyCent[0] -> Fill(fCent, QyC);
    fHist2DCalibPsi1ZNCCent[0] -> Fill(fCent, psiCGE);
    //ZNA
    fProfileZNAQxCent[0] -> Fill(fCent,-QxA);
    fProfileZNAQyCent[0] -> Fill(fCent, QyA);
    fHist2DCalibPsi1ZNACent[0] -> Fill(fCent, psiAGE);
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

  if (IsQAZDC) {
    //ZNC
    fProfileZNCQxCent[1] -> Fill(fCent, QxC);
    fProfileZNCQyCent[1] -> Fill(fCent, QyC);
    fHist2DCalibPsi1ZNCCent[1] -> Fill(fCent, psiCRC);
    //ZNA
    fProfileZNAQxCent[1] -> Fill(fCent,-QxA);
    fProfileZNAQyCent[1] -> Fill(fCent, QyA);
    fHist2DCalibPsi1ZNACent[1] -> Fill(fCent, psiARC);
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
  if (IsQAZDC) {
    fHist2DCalibPsi1ZNCCent[2] -> Fill(fCent, psiCSF);
    fHist2DCalibPsi1ZNACent[2] -> Fill(fCent, psiASF);
  }

  fPsi1ZNC = psiCSF;
  fPsi1ZNA = psiASF;
  return true;
}

//---------------------------------------------------
bool AliAnalysisTaskLambdaProtonCVE::GetZDCPlaneLsFit()
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
  
  double EZNC = 0, wZNC = 0, denZNC = 0, numXZNC = 0, numYZNC = 0;
  double EZNA = 0, wZNA = 0, denZNA = 0, numXZNA = 0, numYZNA = 0; 

  for(int i = 0; i < 4; i++) {
    // ZNC 
    // get energy
    //EZNC = towZNC[i+1];
    // build ZDCC centroid
    wZNC = TMath::Max(0., 4.0 + TMath::Log(towZNC[i+1]/fZNCTowerRawAOD[0]));
    numXZNC += xZDCC[i]*wZNC;
    numYZNC += yZDCC[i]*wZNC;
    denZNC += wZNC;
        
    // ZNA part
    // get energy
    //EZNA = towZNA[i+1];
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

  fHist2DCalibPsi1ZNCCent[0] ->Fill(fCent,psiZNCGE);
  fHist2DCalibPsi1ZNACent[0] ->Fill(fCent,psiZNAGE);


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

  fHist2DCalibPsi1ZNCCent[1] ->Fill(fCent,psiZNCRC);
  fHist2DCalibPsi1ZNACent[1] ->Fill(fCent,psiZNARC);

  fPsi1ZNC = psiZNCRC;
  fPsi1ZNA = psiZNARC;

  return true;
}
//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::LoopTracks()
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
    double phi = track->Phi();
    double  pt = track->Pt();
    double eta = track->Eta();
    int charge = track->Charge();
    int     id = track->GetID();

    fHistPhi[0]->Fill(phi);
    fHist2DEtaPhi[0]->Fill(eta,phi);

    double weight=1;
    if (IsDoNUE) {
      double wEffi = GetNUECor(charge, pt);
      if (wEffi<0) continue;
      else weight *= wEffi;
    }
    if (IsDoNUA) {
      double wAcc = GetNUACor(charge, phi, eta, fVertex[2]);
      if (wAcc<0) continue;
      else weight *= wAcc;
      fHistPhi[1]->Fill(phi, wAcc);
      fHist2DEtaPhi[1]->Fill(eta,phi,wAcc);
    }

    if (pt > fPlanePtMin && pt < fPlanePtMax) {
      ///TODO Use Pt as weight for Better resolution?
      if (eta >= fEtaGapPos) {
        SumQ2xTPCPos   += weight * TMath::Cos(2 * phi);
        SumQ2yTPCPos   += weight * TMath::Sin(2 * phi);
        fWgtMultTPCPos += weight;
        vecPosEPTrkID.push_back(id);
        vecPosEPTrkPhi.push_back(phi);
        vecPosEPTrkNUAWgt.push_back(weight);
      } else if (eta <= fEtaGapNeg) {
        SumQ2xTPCNeg   += weight * TMath::Cos(2 * phi);
        SumQ2yTPCNeg   += weight * TMath::Sin(2 * phi);
        fWgtMultTPCNeg += weight;
        vecNegEPTrkID.push_back(id);
        vecNegEPTrkPhi.push_back(phi);
        vecNegEPTrkNUAWgt.push_back(weight);
      }
    }

    bool isItProttrk = CheckPIDofParticle(track,3); // 3=proton
    isItProttrk *= (pt < fProtonPtMax);

    int code = 0;
    if (charge > 0) {
      code = 999;
      if (isItProttrk) code = 2212;
    } else {  /// charge < 0
      code = -999;
      if (isItProttrk) code = -2212;
    }

    vecPDGCode.push_back(code);
    vecPhi.push_back(phi);
    vecEta.push_back(eta);
    vecPt.push_back(pt);
    vecID.push_back(id);
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::GetTPCPlane()
{
  double psi2TPCPos = GetEventPlane(SumQ2xTPCPos,SumQ2yTPCPos,2);
  double psi2TPCNeg = GetEventPlane(SumQ2xTPCNeg,SumQ2yTPCNeg,2);
  if (TMath::IsNaN(psi2TPCPos) || TMath::IsNaN(psi2TPCNeg)) return false;
  fPsi2TPCPos = psi2TPCPos;
  fPsi2TPCNeg = psi2TPCNeg;
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::LoopV0s()
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
    TVector2 Vt(v0->MomV0X(), v0->MomV0Y());
    double phi = Vt.Phi() > 0 ? Vt.Phi() : Vt.Phi() + TMath::TwoPi();
    int id_posDaughter = v0->GetPosID();
    int id_negDaughter = v0->GetNegID();

    if (code == 3122) {
      double massLambda  = v0->MassLambda();
      fHistLambdaPt[0]              -> Fill(pt);
      fHistLambdaEta[0]             -> Fill(eta);
      fHistLambdaPhi[0]             -> Fill(phi);
      fHistLambdaDcaToPrimVertex[0] -> Fill(dcaToPV);
      fHistLambdaCPA[0]             -> Fill(CPA);
      fHistLambdaDecayLength[0]     -> Fill(dl);
      fHistLambdaMass[0]            -> Fill(massLambda);
      fProfileLambdaMassVsPt[0]     -> Fill(pt, massLambda);

      if (TMath::Abs(massLambda - fMassMean) < fLambdaMassCut) {
        //if a particle has been used as daughter particle before(It happends), we have to refuse a new one.
        if (find(vecDaughterPosID.begin(), vecDaughterPosID.end(), id_posDaughter) != vecDaughterPosID.end()) continue;
        if (find(vecDaughterNegID.begin(), vecDaughterNegID.end(), id_negDaughter) != vecDaughterNegID.end()) continue;
        fHistLambdaPt[1]              -> Fill(pt);
        fHistLambdaEta[1]             -> Fill(eta);
        fHistLambdaPhi[1]             -> Fill(phi);
        fHistLambdaDcaToPrimVertex[1] -> Fill(dcaToPV);
        fHistLambdaCPA[1]             -> Fill(CPA);
        fHistLambdaDecayLength[1]     -> Fill(dl);
        fHistLambdaMass[1]            -> Fill(massLambda);
        fProfileLambdaMassVsPt[1]     -> Fill(pt, massLambda);
        vecLambdaCode.push_back(code);
        vecLambdaPhi.push_back(phi);
        vecLambdaPt.push_back(pt);
        vecDaughterPosID.push_back(id_posDaughter);
        vecDaughterNegID.push_back(id_negDaughter);
      }
    }

    if (code == -3122) {
      double massAntiLambda  = v0->MassAntiLambda();
      fHistAntiLambdaPt[0]              -> Fill(pt);
      fHistAntiLambdaEta[0]             -> Fill(eta);
      fHistAntiLambdaPhi[0]             -> Fill(phi);
      fHistAntiLambdaDcaToPrimVertex[0] -> Fill(dcaToPV);
      fHistAntiLambdaCPA[0]             -> Fill(CPA);
      fHistAntiLambdaDecayLength[0]     -> Fill(dl);
      fHistAntiLambdaMass[0]            -> Fill(massAntiLambda);
      fProfileAntiLambdaMassVsPt[0]     -> Fill(pt, massAntiLambda);

      if (TMath::Abs(massAntiLambda - fMassMean) < fLambdaMassCut) {
        if (find(vecDaughterPosID.begin(), vecDaughterPosID.end(), id_posDaughter) != vecDaughterPosID.end()) continue;
        if (find(vecDaughterNegID.begin(), vecDaughterNegID.end(), id_negDaughter) != vecDaughterNegID.end()) continue;
        fHistAntiLambdaPt[1]              -> Fill(pt);
        fHistAntiLambdaEta[1]             -> Fill(eta);
        fHistAntiLambdaPhi[1]             -> Fill(phi);
        fHistAntiLambdaDcaToPrimVertex[1] -> Fill(dcaToPV);
        fHistAntiLambdaCPA[1]             -> Fill(CPA);
        fHistAntiLambdaDecayLength[1]     -> Fill(dl);
        fHistAntiLambdaMass[1]            -> Fill(massAntiLambda);
        fProfileAntiLambdaMassVsPt[1]     -> Fill(pt, massAntiLambda);

        vecLambdaCode.push_back(code);
        vecLambdaPhi.push_back(phi);
        vecLambdaPt.push_back(pt);
        vecDaughterPosID.push_back(id_posDaughter);
        vecDaughterNegID.push_back(id_negDaughter);
      }
    }
  }//loop V0 end
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::PairLambda()
{
  //Lambda - X
  for (std::vector<double>::size_type jLambda = 0; jLambda < vecLambdaPhi.size(); jLambda++) {
    int    code_lambda    = vecLambdaCode[jLambda];
    double phi_lambda     = vecLambdaPhi[jLambda];
    double pt_lambda      = vecLambdaPt[jLambda];
    int    id_posDaughter = vecDaughterPosID[jLambda];
    int    id_negDaughter = vecDaughterNegID[jLambda];
    //Remove AutoCorrelation:
    ///Get the Total sum of Pos TPC Qx,Qy locally, then Remove AutoCorr if needed.
    double fPosTPCQxTemp = 0., fPosTPCQyTemp = 0.;
    double fNegTPCQxTemp = 0., fNegTPCQyTemp = 0.;
    double fPosTPCMult   = 0., fNegTPCMult   = 0.;
    double qxPosTPC = 0., qyPosTPC = 0.;
    double qxNegTPC = 0., qyNegTPC = 0.;
    // use EP from opposite eta than the charged track! One way to remove AutoCorrelation.
    fNegTPCQxTemp = SumQ2xTPCNeg;
    fNegTPCQyTemp = SumQ2yTPCNeg;
    fNegTPCMult   = fWgtMultTPCNeg;

    fPosTPCQxTemp = SumQ2xTPCPos;
    fPosTPCQyTemp = SumQ2yTPCPos;
    fPosTPCMult   = fWgtMultTPCPos;

    //for NegEP
    std::vector<int>::iterator iterPosDaughterNegTPC = find(vecNegEPTrkID.begin(), vecNegEPTrkID.end(), id_posDaughter);
    std::vector<int>::iterator iterNegDaughterNegTPC = find(vecNegEPTrkID.begin(), vecNegEPTrkID.end(), id_negDaughter);

    //for PosEP
    std::vector<int>::iterator iterPosDaughterPosTPC = find(vecPosEPTrkID.begin(), vecPosEPTrkID.end(), id_posDaughter);
    std::vector<int>::iterator iterNegDaughterPosTPC = find(vecPosEPTrkID.begin(), vecPosEPTrkID.end(), id_negDaughter);

    if (iterPosDaughterNegTPC != vecNegEPTrkID.end()) {
      int iDaughter = distance(vecNegEPTrkID.begin(), iterPosDaughterNegTPC);
      qxNegTPC +=  vecNegEPTrkNUAWgt[iDaughter] * TMath::Cos(2 * vecNegEPTrkPhi[iDaughter]);
      qyNegTPC +=  vecNegEPTrkNUAWgt[iDaughter] * TMath::Sin(2 * vecNegEPTrkPhi[iDaughter]);
      fNegTPCMult -= vecNegEPTrkNUAWgt[iDaughter];
    }

    if (iterNegDaughterNegTPC != vecNegEPTrkID.end()) {
      int iDaughter = distance(vecNegEPTrkID.begin(), iterNegDaughterNegTPC);
      qxNegTPC += vecNegEPTrkNUAWgt[iDaughter] * TMath::Cos(2 * vecNegEPTrkPhi[iDaughter]);
      qyNegTPC += vecNegEPTrkNUAWgt[iDaughter] * TMath::Sin(2 * vecNegEPTrkPhi[iDaughter]);
      fNegTPCMult -= vecNegEPTrkNUAWgt[iDaughter];
    }

    if (iterPosDaughterPosTPC != vecPosEPTrkID.end()) {
      int iDaughter = distance(vecPosEPTrkID.begin(), iterPosDaughterPosTPC);
      qxPosTPC += vecPosEPTrkNUAWgt[iDaughter] * TMath::Cos(2 * vecPosEPTrkPhi[iDaughter]);
      qyPosTPC += vecPosEPTrkNUAWgt[iDaughter] * TMath::Sin(2 * vecPosEPTrkPhi[iDaughter]);
      fPosTPCMult -= vecPosEPTrkNUAWgt[iDaughter];
    }

    if (iterNegDaughterPosTPC != vecPosEPTrkID.end()) {
      int iDaughter = distance(vecPosEPTrkID.begin(), iterNegDaughterPosTPC);
      qxPosTPC +=  vecPosEPTrkNUAWgt[iDaughter] * TMath::Cos(2 * vecPosEPTrkPhi[iDaughter]);
      qyPosTPC +=  vecPosEPTrkNUAWgt[iDaughter] * TMath::Sin(2 * vecPosEPTrkPhi[iDaughter]);
      fPosTPCMult -= vecPosEPTrkNUAWgt[iDaughter];
    }

    fPosTPCQxTemp -= qxPosTPC;   /// qx=0,qy=0 if Lambda daughters are on the opposite eta of the EP used..
    fPosTPCQyTemp -= qyPosTPC;

    fNegTPCQxTemp -= qxNegTPC;   /// qx=0,qy=0 if Lambda daughters are on the opposite eta of the EP used..
    fNegTPCQyTemp -= qyNegTPC;

    if (fPosTPCMult > 0. && fNegTPCMult > 0.) {
      fPosTPCQxTemp = fPosTPCQxTemp/fPosTPCMult;
      fPosTPCQyTemp = fPosTPCQyTemp/fPosTPCMult;
      fNegTPCQxTemp = fNegTPCQxTemp/fNegTPCMult;
      fNegTPCQyTemp = fNegTPCQyTemp/fNegTPCMult;
    } else continue;

    double fPsi2PosTPCNoAuto = GetEventPlane(fPosTPCQxTemp,fPosTPCQyTemp,2.);   //AutoCorrelation Removed Pos EP.
    double fPsi2NegTPCNoAuto = GetEventPlane(fNegTPCQxTemp,fNegTPCQyTemp,2.);   //AutoCorrelation Removed Neg EP.
    if (TMath::IsNaN(fPsi2PosTPCNoAuto) || TMath::IsNaN(fPsi2NegTPCNoAuto)) continue;

    //Check the PID Flow
    if (IsCheckPIDFlow) {
      if (code_lambda == 3122) {
        fProfile2DRawFlowCentPtLambda[0]->Fill(fCent, pt_lambda, TMath::Cos(2 * (phi_lambda - fPsi2PosTPCNoAuto)));
        fProfile2DRawFlowCentPtLambda[1]->Fill(fCent, pt_lambda, TMath::Cos(2 * (phi_lambda - fPsi2V0C)));
        fProfile2DRawFlowCentPtLambda[2]->Fill(fCent, pt_lambda, TMath::Cos(2 * (phi_lambda - fPsi2V0A)));
        fProfile2DRawFlowCentPtLambda[3]->Fill(fCent, pt_lambda, TMath::Cos(2 * (phi_lambda - fPsi1ZNC)));
        fProfile2DRawFlowCentPtLambda[4]->Fill(fCent, pt_lambda, TMath::Cos(2 * (phi_lambda - fPsi1ZNA)));
      }

      if (code_lambda == -3122) {
        fProfile2DRawFlowCentPtAntiLambda[0]->Fill(fCent, pt_lambda, TMath::Cos(2 * (phi_lambda - fPsi2PosTPCNoAuto)));
        fProfile2DRawFlowCentPtAntiLambda[1]->Fill(fCent, pt_lambda, TMath::Cos(2 * (phi_lambda - fPsi2V0C)));
        fProfile2DRawFlowCentPtAntiLambda[2]->Fill(fCent, pt_lambda, TMath::Cos(2 * (phi_lambda - fPsi2V0A)));
        fProfile2DRawFlowCentPtAntiLambda[3]->Fill(fCent, pt_lambda, TMath::Cos(2 * (phi_lambda - fPsi1ZNC)));
        fProfile2DRawFlowCentPtAntiLambda[4]->Fill(fCent, pt_lambda, TMath::Cos(2 * (phi_lambda - fPsi1ZNA)));
      }
    }// check PID flow over.

    for (std::vector<double>::size_type iTrk = 0; iTrk < vecPhi.size(); iTrk++) {
      int    code = vecPDGCode[iTrk];
      int    id   = vecID[iTrk];
      double pt   = vecPt[iTrk];
      double eta  = vecEta[iTrk];
      double phi  = vecPhi[iTrk];
      if (id == id_posDaughter || id == id_negDaughter) continue;

      double fPsiNTPCNoAuto = -999;
      if (eta > 0.) fPsiNTPCNoAuto = fPsi2NegTPCNoAuto;
      else          fPsiNTPCNoAuto = fPsi2PosTPCNoAuto;

      double delta     = TMath::Cos(phi_lambda - phi);
      double gammaTPC  = TMath::Cos(phi_lambda + phi - 2 *fPsiNTPCNoAuto);
      double gammaV0C  = TMath::Cos(phi_lambda + phi - 2 *fPsi2V0C);
      double gammaV0A  = TMath::Cos(phi_lambda + phi - 2 *fPsi2V0A);
      double gammaZNC  = TMath::Cos(phi_lambda + phi - 2 *fPsi1ZNC);
      double gammaZNA  = TMath::Cos(phi_lambda + phi - 2 *fPsi1ZNA);

      if (code > 0 && code_lambda ==  3122) {
        fProfileDelta_Lambda_hPos        -> Fill(fCent, delta);
        fProfileGammaTPC_Lambda_hPos     -> Fill(fCent, gammaTPC);
        fProfileGammaV0C_Lambda_hPos     -> Fill(fCent, gammaV0C);
        fProfileGammaV0A_Lambda_hPos     -> Fill(fCent, gammaV0A);
        fProfileGammaZNC_Lambda_hPos     -> Fill(fCent, gammaZNC);
        fProfileGammaZNA_Lambda_hPos     -> Fill(fCent, gammaZNA);
      }
      if (code < 0 && code_lambda ==  3122) {
        fProfileDelta_Lambda_hNeg        -> Fill(fCent, delta);
        fProfileGammaTPC_Lambda_hNeg     -> Fill(fCent, gammaTPC);
        fProfileGammaV0C_Lambda_hNeg     -> Fill(fCent, gammaV0C);
        fProfileGammaV0A_Lambda_hNeg     -> Fill(fCent, gammaV0A);
        fProfileGammaZNC_Lambda_hNeg     -> Fill(fCent, gammaZNC);
        fProfileGammaZNA_Lambda_hNeg     -> Fill(fCent, gammaZNA);
      }
      if (code > 0 && code_lambda == -3122) {
        fProfileDelta_AntiLambda_hPos    -> Fill(fCent, delta);
        fProfileGammaTPC_AntiLambda_hPos -> Fill(fCent, gammaTPC);
        fProfileGammaV0C_AntiLambda_hPos -> Fill(fCent, gammaV0C);
        fProfileGammaV0A_AntiLambda_hPos -> Fill(fCent, gammaV0A);
        fProfileGammaZNC_AntiLambda_hPos -> Fill(fCent, gammaZNC);
        fProfileGammaZNA_AntiLambda_hPos -> Fill(fCent, gammaZNA);
      }
      if (code < 0 && code_lambda == -3122) {
        fProfileDelta_AntiLambda_hNeg    -> Fill(fCent, delta);
        fProfileGammaTPC_AntiLambda_hNeg -> Fill(fCent, gammaTPC);
        fProfileGammaV0C_AntiLambda_hNeg -> Fill(fCent, gammaV0C);
        fProfileGammaV0A_AntiLambda_hNeg -> Fill(fCent, gammaV0A);
        fProfileGammaZNC_AntiLambda_hNeg -> Fill(fCent, gammaZNC);
        fProfileGammaZNA_AntiLambda_hNeg -> Fill(fCent, gammaZNA);
      }
      if (code ==  2212 && code_lambda ==  3122) {
        fProfileDelta_Lambda_Proton      -> Fill(fCent, delta);
        fProfileGammaTPC_Lambda_Proton   -> Fill(fCent, gammaTPC);
        fProfileGammaV0C_Lambda_Proton   -> Fill(fCent, gammaV0C);
        fProfileGammaV0A_Lambda_Proton   -> Fill(fCent, gammaV0A);
        fProfileGammaZNC_Lambda_Proton   -> Fill(fCent, gammaZNC);
        fProfileGammaZNA_Lambda_Proton   -> Fill(fCent, gammaZNA);
      }
      if (code == -2212 && code_lambda ==  3122) {
        fProfileDelta_Lambda_AntiProton    -> Fill(fCent, delta);
        fProfileGammaTPC_Lambda_AntiProton -> Fill(fCent, gammaTPC);
        fProfileGammaV0C_Lambda_AntiProton -> Fill(fCent, gammaV0C);
        fProfileGammaV0A_Lambda_AntiProton -> Fill(fCent, gammaV0A);
        fProfileGammaZNC_Lambda_AntiProton -> Fill(fCent, gammaZNC);
        fProfileGammaZNA_Lambda_AntiProton -> Fill(fCent, gammaZNA);
      }
      if (code ==  2212 && code_lambda == -3122) {
        fProfileDelta_AntiLambda_Proton    -> Fill(fCent, delta);
        fProfileGammaTPC_AntiLambda_Proton -> Fill(fCent, gammaTPC);
        fProfileGammaV0C_AntiLambda_Proton -> Fill(fCent, gammaV0C);
        fProfileGammaV0A_AntiLambda_Proton -> Fill(fCent, gammaV0A);
        fProfileGammaZNC_AntiLambda_Proton -> Fill(fCent, gammaZNC);
        fProfileGammaZNA_AntiLambda_Proton -> Fill(fCent, gammaZNA);
      }
      if (code == -2212 && code_lambda == -3122) {
        fProfileDelta_AntiLambda_AntiProton    -> Fill(fCent, delta);
        fProfileGammaTPC_AntiLambda_AntiProton -> Fill(fCent, gammaTPC);
        fProfileGammaV0C_AntiLambda_AntiProton -> Fill(fCent, gammaV0C);
        fProfileGammaV0A_AntiLambda_AntiProton -> Fill(fCent, gammaV0A);
        fProfileGammaZNC_AntiLambda_AntiProton -> Fill(fCent, gammaZNC);
        fProfileGammaZNA_AntiLambda_AntiProton -> Fill(fCent, gammaZNA);
      }
    }// (Anti)Lambda-X pair done
  }/// loop over charge particle array

  if (IsCheckPIDFlow) {
    for (std::vector<double>::size_type iTrk = 0; iTrk < vecPhi.size(); iTrk++) {
      int    code = vecPDGCode[iTrk];
      double pt   = vecPt[iTrk];
      double eta  = vecEta[iTrk];
      double phi  = vecPhi[iTrk];
      if (code > 0) {
        if (eta > 0) fProfile2DRawFlowCentPthPos[0]->Fill(fCent, pt, TMath::Cos(2 * (phi-fPsi2TPCNeg)));
        else         fProfile2DRawFlowCentPthPos[0]->Fill(fCent, pt, TMath::Cos(2 * (phi-fPsi2TPCPos)));
        fProfile2DRawFlowCentPthPos[1]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi2V0C)));
        fProfile2DRawFlowCentPthPos[2]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi2V0A)));
        fProfile2DRawFlowCentPthPos[3]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi1ZNC)));
        fProfile2DRawFlowCentPthPos[4]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi1ZNA)));
      }
      if (code < 0) {
        if (eta > 0) fProfile2DRawFlowCentPthNeg[0]->Fill(fCent, pt, TMath::Cos(2 * (phi-fPsi2TPCNeg)));
        else         fProfile2DRawFlowCentPthNeg[0]->Fill(fCent, pt, TMath::Cos(2 * (phi-fPsi2TPCPos)));
        fProfile2DRawFlowCentPthNeg[1]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi2V0C)));
        fProfile2DRawFlowCentPthNeg[2]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi2V0A)));
        fProfile2DRawFlowCentPthNeg[3]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi1ZNC)));
        fProfile2DRawFlowCentPthNeg[4]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi1ZNA)));
      }
      if (code == 2212) {
        if (eta > 0) fProfile2DRawFlowCentPtProton[0]->Fill(fCent, pt, TMath::Cos(2 * (phi-fPsi2TPCNeg)));
        else         fProfile2DRawFlowCentPtProton[0]->Fill(fCent, pt, TMath::Cos(2 * (phi-fPsi2TPCPos)));
        fProfile2DRawFlowCentPtProton[1]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi2V0C)));
        fProfile2DRawFlowCentPtProton[2]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi2V0A)));
        fProfile2DRawFlowCentPtProton[3]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi1ZNC)));
        fProfile2DRawFlowCentPtProton[4]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi1ZNA)));
      }
      if (code == -2212) {
        if (eta > 0) fProfile2DRawFlowCentPtAntiProton[0]->Fill(fCent, pt, TMath::Cos(2 * (phi-fPsi2TPCNeg)));
        else         fProfile2DRawFlowCentPtAntiProton[0]->Fill(fCent, pt, TMath::Cos(2 * (phi-fPsi2TPCPos)));
        fProfile2DRawFlowCentPtAntiProton[1]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi2V0C)));
        fProfile2DRawFlowCentPtAntiProton[2]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi2V0A)));
        fProfile2DRawFlowCentPtAntiProton[3]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi1ZNC)));
        fProfile2DRawFlowCentPtAntiProton[4]->Fill(fCent, pt, TMath::Cos(2 * (phi - fPsi1ZNA)));
      }
    }
  }////--------- PID Flow hist are Filled ----------
  return true;
}

//---------------------------------------------------

void AliAnalysisTaskLambdaProtonCVE::ResetVectors()
{
  SumQ2xTPCPos   = 0.;
  SumQ2yTPCPos   = 0.;
  fWgtMultTPCPos = 0.;

  SumQ2xTPCNeg   = 0.;
  SumQ2yTPCNeg   = 0.;
  fWgtMultTPCNeg = 0.;

  std::vector<int>().swap(vecPosEPTrkID);
  std::vector<int>().swap(vecNegEPTrkID);
  std::vector<double>().swap(vecPosEPTrkPhi);
  std::vector<double>().swap(vecNegEPTrkPhi);
  std::vector<double>().swap(vecPosEPTrkNUAWgt);
  std::vector<double>().swap(vecNegEPTrkNUAWgt);

  std::vector<int>().swap(vecPDGCode);
  std::vector<int>().swap(vecID);
  std::vector<double>().swap(vecPhi);
  std::vector<double>().swap(vecEta);
  std::vector<double>().swap(vecPt);
  std::vector<double>().swap(vecNUAWeight);
  std::vector<double>().swap(vecNUEWeight);
  std::vector<double>().swap(vecNUAWeightPID);
  std::vector<double>().swap(vecNUEWeightPID);

  std::vector<int>().swap(vecLambdaCode);
  std::vector<double>().swap(vecLambdaPhi);
  std::vector<double>().swap(vecLambdaPt);
  std::vector<int>().swap(vecDaughterPosID);
  std::vector<int>().swap(vecDaughterNegID);
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::LoadCalibHistForThisRun()
{
  if (fPeriod.EqualTo("LHC10h")) {
    // 10h VZERO Calibration Histograms is Global
    // 10h ZDC Calibration Histograms
    if (IsZDCCalibOn) {
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
    if (IsVZEROCalibOn) {
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
    if (IsDoNUA) {
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
    if (IsVZEROCalibOn) {
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
    if (IsDoNUA) {
      hCorrectNUAPos -> Reset();
      hCorrectNUANeg -> Reset();
      hCorrectNUAPos = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d",0,fRunNum));
      hCorrectNUANeg = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d",0,fRunNum));
      if (!hCorrectNUAPos) return false;
      if (!hCorrectNUANeg) return false;
    }
    //18q/r ZDC
    if (IsZDCCalibOn) {
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

bool AliAnalysisTaskLambdaProtonCVE::RemovalForRun1()
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

bool AliAnalysisTaskLambdaProtonCVE::RejectEvtMultComp() // 15o_pass1, old pile-up
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

    fHist2DMultMultQA[0]->Fill(multTPC,multEsd);
    fHist2DMultMultQA[1]->Fill(multGlobal,multTPCFE);

    //TODO
    TString fMultComp = "pileupByGlobalTPC1";

    if (fMultComp.EqualTo("pileupByEDSTPC128")) { // Rihan
      if ((Double_t)(multEsd*(1/3.45) - 90.) < (Double_t)multTPC)
      {
        fHist2DMultMultQA[3]->Fill(multTPC,multEsd);
        return true;
      }
      else return false;
    }

    if (fMultComp.EqualTo("pileupByGlobalTPC1")) { // A.Dobrin
      if (multTPCFE-1.78*multGlobal<62.87 && multTPCFE-1.48*multGlobal>-36.73) {
        fHist2DMultMultQA[4]->Fill(multGlobal,multTPCFE);
        return true;
      }
      else return false;
    }

    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::RejectEvtTFFit()
{
  Float_t centV0M = -999;
  Float_t centCL1 = -999;
  Float_t centCL0 = -999;

  AliMultSelection* fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
  if (!fMultSelection) {
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
      if (!aodTrk) {
          delete aodTrk;
          continue;
      }
      if (aodTrk->TestFilterBit(32)) {
        if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
      }
  }

  fHist2DMultCentQA[0]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)

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
  if (centCL0 < fCenCutLowPU->Eval(centV0M)) return false;
  if (centCL0 > fCenCutHighPU->Eval(centV0M)) return false;
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) return false;
  if (multV0On < fV0CutPU->Eval(multV0Tot)) return false;
  if (Float_t(multTrk) < fMultCutPU->Eval(centV0M)) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fAOD->IsIncompleteDAQ()) return false;

  fHist2DMultCentQA[1]->Fill(centV0M, multTrk); //  Mult(FB32) Vs Cent(V0M)
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::RejectEvtTPCITSfb32TOF ()
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
    fHist2DMultMultQA[2]->Fill(multTrkTOF, nTrk);
    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::AODPileupCheck()
{
  Int_t isPileup = fAOD->IsPileupFromSPD(3);
  if (isPileup !=0 && fPeriod.EqualTo("LHC16t")) return false; // LHC16t : pPb
  if (fAOD->IsIncompleteDAQ()) return false;
  if (((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) return false;
  if (fPeriod.EqualTo("LHC15o")) {
    AliMultSelection* fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
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

bool AliAnalysisTaskLambdaProtonCVE::PileUpMultiVertex()
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
    if (wDst<kMinWDist)        continue;

    return true; // pile-up: well separated vertices
  }
  return false;
}

//---------------------------------------------------

double AliAnalysisTaskLambdaProtonCVE::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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

bool AliAnalysisTaskLambdaProtonCVE::AcceptAODTrack(AliAODTrack *track)
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

    if ( pt < fPtMin
      || pt > fPtMax
      || fabs(eta) > fEtaCut
      || fabs(nhits) < fNclsCut
      || chi2 < fChi2Min 
      || chi2 > fChi2Max
      || dedx < fDedxCut) 
    return false;

    if (fFilterBit != 768){
      if (fPeriod.EqualTo("LHC10h")) {
        //------------------
        // dca cut
        //------------------
        double mag = fAOD->GetMagneticField();
        double dcaxy  = 999.;
        double dcaz   = 999.;
        double r[3];
        double dca[2];
        double cov[3];
        double vx = fVertex[0];
        double vy = fVertex[1];
        double vz = fVertex[2];
        bool proptodca = track->PropagateToDCA(fAOD->GetPrimaryVertex(), mag, 100., dca, cov);
        if (track->GetXYZ(r)) {
          dcaxy = r[0];
          dcaz  = r[1];
        } else {
          double dcax = r[0] - vx;
          double dcay = r[1] - vy;
          dcaz  = r[2] - vz;
          dcaxy = sqrt(dcax*dcax + dcay*dcay);
        }
        if (fabs(dcaxy)>fDcaCutxy) return false;
        if (fabs(dcaz)>fDcaCutz) return false;
        fHistDcaXY->Fill(dcaxy);
        fHistDcaZ->Fill(dcaz);
      }
    }

    fHistPt->Fill(pt);
    fHistEta->Fill(eta);
    fHistNhits->Fill(nhits);
    fHist2DPDedx->Fill(track->P()*charge, dedx);

    return true;
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck)
{
  if (pidToCheck==0) return kTRUE;    //// Charge Particles do not need PID check
  bool bPIDokay = kFALSE;

  if (!fPIDResponse) {
    Printf("\n Could Not access PIDResponse Task, please Add the Task...\n return with kFALSE pid\n");
    return kFALSE;
  }

  /// Rihan todo: To set the low pT cuts for nSigmaTPC from AddTaskMacro!
  /// Although someone barely needs to change it given the purity..

  float nSigTPC = 0, nSigTOF = 0, nSigRMS = 0;
  bool trkPtPID  = ftrack->Pt();
  int  trkChargePID = ftrack->Charge();

  ///Pion =>
  if (pidToCheck==1) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);//Some warning show here (***TDatabasePDG::AddParicle: particle with PDGcode = 3124 already defind),I don't understand what happended. --chunzheng
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if (trkPtPID<=0.5 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut)     bPIDokay = kTRUE;
    // Using TPCTOF RMS cut for higher pt:
    else if (trkPtPID>0.5 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut) bPIDokay = kTRUE;
    return bPIDokay;
  }
  ///Kaon =>
  else if (pidToCheck==2) {
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kKaon);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kKaon);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if (trkPtPID<=0.45 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut)     bPIDokay = kTRUE;
    else if (trkPtPID>0.45 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut) bPIDokay = kTRUE;
    return bPIDokay;
  }
  ///proton =>
  else if (pidToCheck==3) {///
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kProton);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kProton);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if (trkPtPID<=0.6 && TMath::Abs(nSigTPC)<=fNSigmaTPCCut) {
      bPIDokay = kTRUE;
      if (trkChargePID>0 && trkPtPID<0.4) bPIDokay = kFALSE;
    }
    else if (trkPtPID>0.6 && TMath::Abs(nSigRMS)<=fNSigmaTOFCut) {
      bPIDokay = kTRUE;
    }
    return bPIDokay;
  }
  else{
    Printf("\n -Ve number not allowed! Choose among: 0,1,2,3 (Charge Pion, Kaon, Proton)\n return with kFALSE \n");
    return kFALSE;
  }

  return kFALSE;
}

//---------------------------------------------------

double AliAnalysisTaskLambdaProtonCVE::GetNUECor(int charge, double pt)
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

double AliAnalysisTaskLambdaProtonCVE::GetNUACor(int charge, double phi, double eta, double vz)
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

bool AliAnalysisTaskLambdaProtonCVE::IsGoodV0(AliAODv0 *aodV0)
{
  // Offline reconstructed V0 only
  if (aodV0->GetOnFlyStatus()) return false;
  // Get daughters and check them
  AliAODTrack *myTrackNegTest = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(1));
  AliAODTrack *myTrackPosTest = dynamic_cast<AliAODTrack*>(aodV0->GetDaughter(0));
  if (!myTrackPosTest || !myTrackNegTest) {
    Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
    return false;
  }
  // Unlike signs of daughters
  if (myTrackNegTest->Charge() == myTrackPosTest->Charge()) return false;
  // Cosinus of pointing angle
  double dCPA = aodV0->CosPointingAngle(fVertex);
  // cut on Cosinus of pointing angle
  if (dCPA < fV0CPAMin) return false;
  // DCA of V0
  double dV0Dca = aodV0->DcaV0ToPrimVertex();
  if (TMath::Abs(dV0Dca) > fV0DCAToPrimVtxMax) return false;
  // V0 path length before decay
  double dDecayLength = aodV0->DecayLengthV0(fVertex);
  if (dDecayLength > fV0DecayLengthMax) return false;
  if (dDecayLength < fV0DecayLengthMin) return false;
  // DCA between daughters
  double dDCA = aodV0->DcaV0Daughters();
  if (dDCA > fV0DcaBetweenDaughtersMax) return false;
  double dPt = aodV0->Pt();
  if (dPt < fV0PtMin ) return false;
  double dRapidity = aodV0->RapLambda();
  if (TMath::Abs(dRapidity) > fV0RapidityMax) return false;
  return kTRUE;
}

//---------------------------------------------------

bool AliAnalysisTaskLambdaProtonCVE::IsGoodDaughterTrack(const AliAODTrack *track)
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

int AliAnalysisTaskLambdaProtonCVE::GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *nTrack)
{
  bool IsLambda     = kFALSE;
  bool IsAntiLambda = kFALSE;
  int  code = 0;

  //Λ-->(p+)+(π-)
  float nSigTPCPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton));//TPC p+
  float nSigTPCNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion));//TPC π-
  //(Λ-)-->(p-)+(π+)
  float nSigTPCPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion));//TPC π+
  float nSigTPCNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton));//TPC p-

  IsLambda     = (nSigTPCPosProton < fV0PosProtonTPCNsigma) && (nSigTPCNegPion < fV0NegPionTPCNsigma);
  IsAntiLambda = (nSigTPCNegProton < fV0NegProtonTPCNsigma) && (nSigTPCPosPion < fV0PosPionTPCNsigma);

  if (IsV0DaughterUseTOF) {
    float nSigTOFPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kProton));//TOF p+
    float nSigTOFNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kPion));//TOF π-
    float nSigTOFPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(pTrack, AliPID::kPion));//TOF π+
    float nSigTOFNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(nTrack, AliPID::kProton));//TOF p-

    IsLambda     *= (nSigTOFPosProton < fV0PosProtonTOFNsigma && nSigTOFNegPion < fV0NegPionTOFNsigma);
    IsAntiLambda *= (nSigTOFNegProton < fV0NegProtonTOFNsigma && nSigTOFPosPion < fV0PosPionTOFNsigma);
  }

  if (IsLambda)     code =  3122;
  if (IsAntiLambda) code = -3122;
  if (IsLambda && IsAntiLambda) code = 0;
  return code;
}

//---------------------------------------------------

double AliAnalysisTaskLambdaProtonCVE::GetEventPlane(double qx, double qy, double harmonic)
{
  double psi = (1./harmonic)*TMath::ATan2(qy,qx);
  if (psi < 0) return psi += TMath::TwoPi()/harmonic;
  else return psi;
}

//---------------------------------------------------

