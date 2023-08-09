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
#include "AliAnalysisTaskCVEPIDCMEDiff.h"

ClassImp(AliAnalysisTaskCVEPIDCMEDiff);

//---------------------------------------------------
AliAnalysisTaskCVEPIDCMEDiff::AliAnalysisTaskCVEPIDCMEDiff() :
  AliAnalysisTaskSE(),
  isCalculateLambdaHadron(false),
  isDoNUE(true),
  isDoLambdaNUE(false),
  isV0DaughterUseTOF(false),
  isNarrowDcaCuts768(false),
  isProtonCustomizedDCACut(false),
  isUsePionRejection(false),
  isTightPileUp(false),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fPlaneEstimator("TPC"),
  fPlanePtMin(0.2),
  fPlanePtMax(2.0),
  fFilterBit(768),
  fNclsCut(70),
  fChi2Max(2.5),
  fChi2Min(0.1),
  fDcaCutXY(2.4),
  fDcaCutZ(3.2),
  fPtMin(0.2),
  fPtMax(5.),
  fEtaCut(0.8),
  fDedxCut(10.),
  fProtonPtMin(0.4),
  fProtonPtMax(5.0),
  fAntiProtonPtMin(0.4),
  fAntiProtonPtMax(5.0),
  fNSigmaTPCCut(3.0),
  fNSigmaTOFCut(3.0),
  fV0CPAMin(0.997),
  fV0DecayLengthMin(3.),
  fV0DecayLengthMax(100.),
  fV0DCAToPrimVtxMax(1.5),
  fV0DcaBetweenDaughtersMax(0.5),
  fDaughtersPtMax(20.0),
  fDaughtersEtaMax(0.8),
  fDaughtersTPCNclsMin(70),
  fDaughtersDCAToPrimVtxMin(0.05),
  fRatioCrossedRowsFindable(0.8),
  fV0PosProtonTPCNsigma(3.0),
  fV0NegPionTPCNsigma(3.0),
  fV0NegProtonTPCNsigma(3.0),
  fV0PosPionTPCNsigma(3.0),
  fV0PosProtonTOFNsigma(4.0),
  fV0NegPionTOFNsigma(4.0),
  fV0NegProtonTOFNsigma(4.0),
  fV0PosPionTOFNsigma(4.0),
  fLambdaPtMin(0.5),
  fLambdaPtMax(10.0),
  fAntiLambdaPtMin(0.5),
  fAntiLambdaPtMax(10.0),
  fLambdaRapidityMax(0.5),
  fAntiLambdaRapidityMax(0.5),
  fLambdaMassRightCut(0.02),
  fLambdaMassLeftCut(0.02),
  fAntiLambdaMassRightCut(0.02),
  fAntiLambdaMassLeftCut(0.02),
  fLambdaMassMean(1.115683),
  fNMassBins(30),
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
  fSumQ2xTPC(0.),
  fSumQ2yTPC(0.),
  fWgtMultTPC(0.),
  fPsi2(-999),
  mapTPCTrksIDPhiWgt(0),
  vecParticle(0),
  vecParticleV0(0),
  fSPDCutPU(nullptr),
  fV0CutPU(nullptr),
  fCenCutLowPU(nullptr),
  fCenCutHighPU(nullptr),
  fMultCutPU(nullptr),
  fListNUE(nullptr),
  hNUEweightPlus(nullptr),
  hNUEweightMinus(nullptr),
  fListNUA(nullptr),
  hCorrectNUAPos(nullptr),
  hCorrectNUANeg(nullptr),
  fListVZEROCalib(nullptr),
  contQxncm(nullptr),
  contQyncm(nullptr),
  hQx2mV0C(nullptr),
  hQy2mV0C(nullptr),
  fHCorrectV0ChWeghts(nullptr),
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
  fHistV0NegDaughterDca(nullptr),
  fHistV0PosDaughterDca(nullptr),
  fResultsList(nullptr),
  fHist2Psi2(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;
  for (int i = 0; i < 2; i++) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistVz[i] = nullptr;
  for (int i = 0; i < 8; i++) fHist2CentQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2MultCentQA[i] = nullptr;
  for (int i = 0; i < 6; i++) fHist2MultMultQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaNegDaugtherDca[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaPosDaugtherDca[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaCPA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaDecayLength[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist3LambdaCentPtMass[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2LambdaMassPtY[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaNegDaugtherDca[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaPosDaugtherDca[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaCPA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaDecayLength[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist3AntiLambdaCentPtMass[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2AntiLambdaMassPtY[i] = nullptr;
  for (int i = 0; i < 4; i++) fHist2LambdaMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fHist3LambdaProtonMassSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fHist3LambdaProtonMassDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffDeltaLambdaProtonMassSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffDeltaLambdaProtonMassDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffGammaLambdaProtonMassSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffGammaLambdaProtonMassDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DDeltaLambdaHadronMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DDeltaLambdaHadronMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DGammaLambdaHadronMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DGammaLambdaHadronMass[i] = nullptr;

  for (int i = 0; i < 4; i++) fHist3LambdaProtonMassDCA[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffDeltaLambdaProtonMassDCA[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffGammaLambdaProtonMassDCA[i] = nullptr;
}

//---------------------------------------------------
AliAnalysisTaskCVEPIDCMEDiff::AliAnalysisTaskCVEPIDCMEDiff(const char *name) :
  AliAnalysisTaskSE(name),
  isCalculateLambdaHadron(false),
  isDoNUE(true),
  isDoLambdaNUE(false),
  isV0DaughterUseTOF(false),
  isNarrowDcaCuts768(false),
  isProtonCustomizedDCACut(false),
  isUsePionRejection(false),
  isTightPileUp(false),
  fTrigger("kINT7"),
  fPeriod("LHC18q"),
  fVzCut(10.0),
  fCentDiffCut(7.5),
  fPlaneEstimator("TPC"),
  fPlanePtMin(0.2),
  fPlanePtMax(2.0),
  fFilterBit(768),
  fNclsCut(70),
  fChi2Max(2.5),
  fChi2Min(0.1),
  fDcaCutXY(2.4),
  fDcaCutZ(3.2),
  fPtMin(0.2),
  fPtMax(5.),
  fEtaCut(0.8),
  fDedxCut(10.),
  fProtonPtMin(0.4),
  fProtonPtMax(5.0),
  fAntiProtonPtMin(0.4),
  fAntiProtonPtMax(5.0),
  fNSigmaTPCCut(3.0),
  fNSigmaTOFCut(3.0),
  fV0CPAMin(0.997),
  fV0DecayLengthMin(3.),
  fV0DecayLengthMax(100.),
  fV0DCAToPrimVtxMax(1.5),
  fV0DcaBetweenDaughtersMax(0.5),
  fDaughtersPtMax(20.0),
  fDaughtersEtaMax(0.8),
  fDaughtersTPCNclsMin(70),
  fDaughtersDCAToPrimVtxMin(0.05),
  fRatioCrossedRowsFindable(0.8),
  fV0PosProtonTPCNsigma(3.0),
  fV0NegPionTPCNsigma(3.0),
  fV0NegProtonTPCNsigma(3.0),
  fV0PosPionTPCNsigma(3.0),
  fV0PosProtonTOFNsigma(4.0),
  fV0NegPionTOFNsigma(4.0),
  fV0NegProtonTOFNsigma(4.0),
  fV0PosPionTOFNsigma(4.0),
  fLambdaPtMin(0.5),
  fLambdaPtMax(10.0),
  fAntiLambdaPtMin(0.5),
  fAntiLambdaPtMax(10.0),
  fLambdaRapidityMax(0.5),
  fAntiLambdaRapidityMax(0.5),
  fLambdaMassRightCut(0.02),
  fLambdaMassLeftCut(0.02),
  fAntiLambdaMassRightCut(0.02),
  fAntiLambdaMassLeftCut(0.02),
  fLambdaMassMean(1.115683),
  fNMassBins(30),
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
  fSumQ2xTPC(0.),
  fSumQ2yTPC(0.),
  fWgtMultTPC(0.),
  fPsi2(-999),
  mapTPCTrksIDPhiWgt(0),
  vecParticle(0),
  vecParticleV0(0),
  fSPDCutPU(nullptr),
  fV0CutPU(nullptr),
  fCenCutLowPU(nullptr),
  fCenCutHighPU(nullptr),
  fMultCutPU(nullptr),
  fListNUE(nullptr),
  hNUEweightPlus(nullptr),
  hNUEweightMinus(nullptr),
  fListNUA(nullptr),
  hCorrectNUAPos(nullptr),
  hCorrectNUANeg(nullptr),
  fListVZEROCalib(nullptr),
  contQxncm(nullptr),
  contQyncm(nullptr),
  hQx2mV0C(nullptr),
  hQy2mV0C(nullptr),
  fHCorrectV0ChWeghts(nullptr),
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
  fHistV0NegDaughterDca(nullptr),
  fHistV0PosDaughterDca(nullptr),
  fResultsList(nullptr),
  fHist2Psi2(nullptr)
{
  for (int i = 0; i < 3; i++) fVertex[i] = -999;
  for (int i = 0; i < 2; i++) fHistCent[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistVz[i] = nullptr;
  for (int i = 0; i < 8; i++) fHist2CentQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2MultCentQA[i] = nullptr;
  for (int i = 0; i < 6; i++) fHist2MultMultQA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaNegDaugtherDca[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaPosDaugtherDca[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaCPA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistLambdaDecayLength[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist3LambdaCentPtMass[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2LambdaMassPtY[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaPt[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaEta[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaPhi[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaDcaToPrimVertex[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaNegDaugtherDca[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaPosDaugtherDca[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaCPA[i] = nullptr;
  for (int i = 0; i < 2; i++) fHistAntiLambdaDecayLength[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist3AntiLambdaCentPtMass[i] = nullptr;
  for (int i = 0; i < 2; i++) fHist2AntiLambdaMassPtY[i] = nullptr;
  for (int i = 0; i < 4; i++) fHist2LambdaMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fHist3LambdaProtonMassSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fHist3LambdaProtonMassDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffDeltaLambdaProtonMassSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffDeltaLambdaProtonMassDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffGammaLambdaProtonMassSPt[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffGammaLambdaProtonMassDEta[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DDeltaLambdaHadronMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DDeltaLambdaHadronMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DGammaLambdaHadronMass[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile2DGammaLambdaHadronMass[i] = nullptr;

  for (int i = 0; i < 4; i++) fHist3LambdaProtonMassDCA[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffDeltaLambdaProtonMassDCA[i] = nullptr;
  for (int i = 0; i < 4; i++) fProfile3DDiffGammaLambdaProtonMassDCA[i] = nullptr;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
}

//------------------------------------------------

AliAnalysisTaskCVEPIDCMEDiff::~AliAnalysisTaskCVEPIDCMEDiff()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  if (fQAList) delete fQAList;
  if (fResultsList) delete fResultsList;
}

//---------------------------------------------------

void AliAnalysisTaskCVEPIDCMEDiff::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate");
}

//---------------------------------------------------

void AliAnalysisTaskCVEPIDCMEDiff::UserCreateOutputObjects()
{
  ////////////////////////
  // Pile up Function
  ////////////////////////
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
  if (fPlaneEstimator.EqualTo("TPC")) {
    if (!fListNUA) {
      std::cout<<("NUA list not found")<<std::endl;
      return;
    }
    if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
      hCorrectNUAPos = new TH3F();
      hCorrectNUANeg = new TH3F();
    }
  }
  ////////////////////////
  // V0C
  ////////////////////////
  if (fPlaneEstimator.EqualTo("V0C")) {
    if(!fListVZEROCalib) {
      std::cout<<("V0C list not found")<<std::endl;
      return;
    }
    contQxncm = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqxc%im",2)); // V0C Qx Mean
    contQyncm = (AliOADBContainer*)fListVZEROCalib->FindObject(Form("fqyc%im",2)); // V0C Qy Mean
    hQx2mV0C = new TH1D();
    hQy2mV0C = new TH1D();
    fHCorrectV0ChWeghts = new TH2F();
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
  if      (fPeriod.EqualTo("LHC18q")) for (int i = 0; i <125; i++) runNumList->insert(std::pair<int,int>(runNumList18q[i].Atoi(),i+1));
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
  fHist2MultCentQA[0] = new TH2D("fHist2MultCentQA_BfCut", ";centV0M;multFB32", 80, 0, 80., 50, 0, 5000.);
  fHist2MultCentQA[1] = new TH2D("fHist2MultCentQA_AfCut", ";centV0M;multFB32", 80, 0, 80., 50, 0, 5000.);
  fQAList->Add(fHistCent[0]);
  fQAList->Add(fHistCent[1]);
  fQAList->Add(fHistVz[0]);
  fQAList->Add(fHistVz[1]);
  for (int i = 0; i < 8; i++) fQAList->Add(fHist2CentQA[i]);
  fQAList->Add(fHist2MultCentQA[0]);
  fQAList->Add(fHist2MultCentQA[1]);

  // track-wise QA
  fHistPt  = new TH1D("fHistPt", ";p_{T}", 200, 0., 20.);
  fHistEta = new TH1D("fHistEta", ";#eta", 100, -2., 2.);
  fHistNhits = new TH1D("fHistNhits", ";nhits", 200, 0., 200.);
  fHist2PDedx = new TH2D("fHist2PDedx", ";pdedx", 400, -10., 10., 400, 0, 1000);
  fHistDcaXY = new TH1D("fHistDcaXY", ";DcaXY", 500, 0., 1);
  fHistDcaZ  = new TH1D("fHistDcaZ",  ";DcaZ", 500, 0., 1);
  fHistPhi[0] = new TH1D("fHistPhi", ";#phi", 100, 0, TMath::TwoPi());
  fHistPhi[1] = new TH1D("fHistPhi_afterNUA", ";#phi", 100, 0, TMath::TwoPi());
  fQAList->Add(fHistPt);
  fQAList->Add(fHistEta);
  fQAList->Add(fHistNhits);
  fQAList->Add(fHist2PDedx);
  fQAList->Add(fHistDcaXY);
  fQAList->Add(fHistDcaZ);
  for (int i = 0; i < 2; i++) fQAList->Add(fHistPhi[i]);

  //Proton QA
  fHistProtonPt      = new TH1D("fHistProtonPt"     , "fHistProtonPt;p_{T}"      , 200, 0., 20.);
  fHistProtonEta     = new TH1D("fHistProtonEta"    , "fHistProtonEta;#eta"      , 100, -2., 2.);
  fHistProtonPhi     = new TH1D("fHistProtonPhi"    , "fHistProtonPhi;#phi"      , 360, 0., TMath::TwoPi());
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
  fHistAntiProtonPhi     = new TH1D("fHistAntiProtonPhi"    , "fHistAntiProtonPhi;#phi"      , 360, 0., TMath::TwoPi()); 
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
  fHistV0NegDaughterDca = new TH1D("hV0NegDaughterDca","",200,0.,20.);
  fHistV0PosDaughterDca = new TH1D("hV0PosDaughterDca","",200,0.,20.);
  fQAList->Add(fHistV0Pt);
  fQAList->Add(fHistV0Eta);
  fQAList->Add(fHistV0DcatoPrimVertex);
  fQAList->Add(fHistV0CPA);
  fQAList->Add(fHistV0DecayLength);
  fQAList->Add(fHistV0NegDaughterDca);
  fQAList->Add(fHistV0PosDaughterDca);

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
    fHistLambdaNegDaugtherDca[i] = new TH1D(Form("hLambdaNegDaugtherDca_%sMassCut",charLambdaMassCut.data()),"NegDaughterDcatoPV",200, 0., 20.);
    fHistLambdaPosDaugtherDca[i] = new TH1D(Form("hLambdaPosDaugtherDca_%sMassCut",charLambdaMassCut.data()),"PosDaughterDcatoPV",200, 0., 20.);
    fHistLambdaCPA[i] = new TH1D(Form("hLambdaCPA_%sMassCut",charLambdaMassCut.data()),";CPA",200, 0.9, 1.);
    fHistLambdaDecayLength[i] = new TH1D(Form("hLambdaDecayLength_%sMassCut",charLambdaMassCut.data()),";DecayLength", 250, 0., 500.);
    fHist2LambdaMassPtY[i] = new TH2D(Form("fHist2LambdaMassPtY_%sMassCut",charLambdaMassCut.data()),";pT;yH",12,0.,6.,20,-1.,1.);
    fQAList->Add(fHistLambdaPt[i]);
    fQAList->Add(fHistLambdaEta[i]);
    fQAList->Add(fHistLambdaPhi[i]);
    fQAList->Add(fHistLambdaDcaToPrimVertex[i]);
    fQAList->Add(fHistLambdaNegDaugtherDca[i]);
    fQAList->Add(fHistLambdaPosDaugtherDca[i]);
    fQAList->Add(fHistLambdaCPA[i]);
    fQAList->Add(fHistLambdaDecayLength[i]);
    fQAList->Add(fHist2LambdaMassPtY[i]);

    // AntiLambda
    fHistAntiLambdaPt[i] = new TH1D(Form("hAntiLambdaPt_%sMassCut",charLambdaMassCut.data()),";pT", 200, 0., 20.);
    fHistAntiLambdaEta[i] = new TH1D(Form("hAntiLambdaEta_%sMassCut",charLambdaMassCut.data()),";#eta",200, -3., 3.);
    fHistAntiLambdaPhi[i] = new TH1D(Form("hAntiLambdaPhi_%sMassCut",charLambdaMassCut.data()),";phi", 360, 0., TMath::TwoPi());
    fHistAntiLambdaDcaToPrimVertex[i] = new TH1D(Form("hAntiLambdaDcaToPrimVertex_%sMassCut",charLambdaMassCut.data()),"DcatoPV",200, 0., 20.);
    fHistAntiLambdaNegDaugtherDca[i] = new TH1D(Form("hAntiLambdaNegDaugtherDca_%sMassCut",charLambdaMassCut.data()),"NegDaughterDcatoPV",200, 0., 20.);
    fHistAntiLambdaPosDaugtherDca[i] = new TH1D(Form("hAntiLambdaPosDaugtherDca_%sMassCut",charLambdaMassCut.data()),"PosDaughterDcatoPV",200, 0., 20.);
    fHistAntiLambdaCPA[i] = new TH1D(Form("hAntiLambdaCPA_%sMassCut",charLambdaMassCut.data()),";CPA",200, 0.9, 1.);
    fHistAntiLambdaDecayLength[i] = new TH1D(Form("hAntiLambdaDecayLength_%sMassCut",charLambdaMassCut.data()),";DecayLength", 250, 0., 500.);
    fHist2AntiLambdaMassPtY[i] = new TH2D(Form("fHist2AntiLambdaMassPtY_%sMassCut",charLambdaMassCut.data()),";pT;yH",12,0.,6.,20,-1.,1.);
    fQAList->Add(fHistAntiLambdaPt[i]);
    fQAList->Add(fHistAntiLambdaEta[i]);
    fQAList->Add(fHistAntiLambdaPhi[i]);
    fQAList->Add(fHistAntiLambdaDcaToPrimVertex[i]);
    fQAList->Add(fHistAntiLambdaNegDaugtherDca[i]);
    fQAList->Add(fHistAntiLambdaPosDaugtherDca[i]);
    fQAList->Add(fHistAntiLambdaCPA[i]);
    fQAList->Add(fHistAntiLambdaDecayLength[i]);
    fQAList->Add(fHist2AntiLambdaMassPtY[i]);
  }
  fHist3LambdaCentPtMass[0] = new TH3D("fHist3LambdaCentPtMass_BfMassCut",";Centrality;Pt;Mass",8,0,80,20,0,10,1000,1.,1.25); //  Current bin size = 0.00025
  fHist3LambdaCentPtMass[1] = new TH3D("fHist3LambdaCentPtMass_AfMassCut",";Centrality;Pt;Mass",8,0,80,20,0,10,fNMassBins,fLambdaMassMean-fLambdaMassRightCut,fLambdaMassMean+fLambdaMassLeftCut);
  fQAList->Add(fHist3LambdaCentPtMass[0]);
  fQAList->Add(fHist3LambdaCentPtMass[1]);

  fHist3AntiLambdaCentPtMass[0] = new TH3D("fHist3AntiLambdaCentPtMass_BfMassCut",";Centrality;Pt;Mass",8,0,80,20,0,10,1000,1.,1.25); //  Current bin size = 0.00025
  fHist3AntiLambdaCentPtMass[1] = new TH3D("fHist3AntiLambdaCentPtMass_AfMassCut",";Centrality;Pt;Mass",8,0,80,20,0,10,fNMassBins,fLambdaMassMean-fLambdaMassRightCut,fLambdaMassMean+fLambdaMassLeftCut);
  fQAList->Add(fHist3AntiLambdaCentPtMass[0]);
  fQAList->Add(fHist3AntiLambdaCentPtMass[1]);

  PostData(1,fQAList);
  if (fDebug) Printf("Post fQAList Data Success!");

  ////////////////////////
  // Results
  ////////////////////////
  fResultsList = new TList();
  fResultsList -> SetName("fResultsList");
  fResultsList -> SetOwner(kTRUE);
  fHist2Psi2 = new TH2D("fHist2Psi2","Psi2 Dist.;centrality;#Psi2",8,0.,80.,100,0.,TMath::Pi());
  fResultsList->Add(fHist2Psi2);

  for (int iType = 0; iType < 4; iType++) {
    // Inv Mass All
    fHist2LambdaMass[iType] = new TH2D(Form("fHist2LambdaMass_%i",iType),";Counts;Mass",7,0,70,fNMassBins,fLambdaMassMean-fLambdaMassRightCut,fLambdaMassMean+fLambdaMassLeftCut);
    fResultsList->Add(fHist2LambdaMass[iType]);
   
    if(!isCalculateLambdaHadron) {
      // Inv Mass different kinematic regions
      fHist3LambdaProtonMassSPt[iType]  = new TH3D(Form("fHist3LambdaProtonMassSPt_%i",iType),  Form("fHist3LambdaProtonMassSPt_%i",iType), 7, 0, 70, 3, 0, 3, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
      fHist3LambdaProtonMassDEta[iType] = new TH3D(Form("fHist3LambdaProtonMassDEta_%i",iType),Form("fHist3LambdaProtonMassDEta_%i",iType), 7, 0, 70, 2, 0, 2, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
      fResultsList->Add(fHist3LambdaProtonMassSPt[iType]);
      fResultsList->Add(fHist3LambdaProtonMassDEta[iType]);
      // Diff δ(ΔpT)
      fProfile3DDiffDeltaLambdaProtonMassSPt[iType] = new TProfile3D(Form("fProfile3DDiffDeltaLambdaProtonMassSPt_%i",iType), Form("fProfile3DDiffDeltaLambdaProtonMassSPt_%i",iType), 7, 0, 70, 3, 0, 3, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
      // Diff δ(Δη)
      fProfile3DDiffDeltaLambdaProtonMassDEta[iType] = new TProfile3D(Form("fProfile3DDiffDeltaLambdaProtonMassDEta_%i",iType), Form("fProfile3DDiffDeltaLambdaProtonMassDEta_%i",iType), 7, 0, 70, 2, 0, 2, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
      // Diff γ(SpT)(only TPC Plane)
      fProfile3DDiffGammaLambdaProtonMassSPt[iType] = new TProfile3D(Form("fProfile3DDiffGammaLambdaProtonMassSPt_%i",iType), Form("fProfile3DDiffGammaLambdaProtonMassSPt_%i",iType), 7, 0, 70, 3, 0, 3, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
      // Diff γ(Δη)(only TPC Plane)
      fProfile3DDiffGammaLambdaProtonMassDEta[iType] = new TProfile3D(Form("fProfile3DDiffGammaLambdaProtonMassDEta_%i",iType), Form("fProfile3DDiffGammaLambdaProtonMassDEta_%i",iType), 7, 0, 70, 2, 0, 2, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
      fResultsList->Add(fProfile3DDiffDeltaLambdaProtonMassSPt[iType]);
      fResultsList->Add(fProfile3DDiffDeltaLambdaProtonMassDEta[iType]);
      fResultsList->Add(fProfile3DDiffGammaLambdaProtonMassSPt[iType]);
      fResultsList->Add(fProfile3DDiffGammaLambdaProtonMassDEta[iType]);
    }

    if(isCalculateLambdaHadron) {
      fProfile2DDeltaLambdaHadronMass[iType] = new TProfile2D(Form("fProfile2DDeltaLambdaHadronMass_%i",iType), Form("fProfile2DDeltaLambdaHadronMass_%i",iType), 7, 0, 70, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
      fProfile2DGammaLambdaHadronMass[iType] = new TProfile2D(Form("fProfile2DGammaLambdaHadronMass_%i",iType), Form("fProfile2DGammaLambdaHadronMass_%i",iType), 7, 0, 70, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
      fResultsList->Add(fProfile2DDeltaLambdaHadronMass[iType]);
      fResultsList->Add(fProfile2DGammaLambdaHadronMass[iType]);
    }

    fHist3LambdaProtonMassDCA[iType] = new TH3D(Form("fHist3LambdaProtonMassDCA_%i",iType),Form("fHist3LambdaProtonMassDCA_%i",iType), 7, 0, 70, 4, 0, 4, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
    fProfile3DDiffDeltaLambdaProtonMassDCA[iType] = new TProfile3D(Form("fProfile3DDiffDeltaLambdaProtonMassDCA_%i",iType), Form("fProfile3DDiffDeltaLambdaProtonMassDCA_%i",iType), 7, 0, 70, 4, 0, 4, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
    fProfile3DDiffGammaLambdaProtonMassDCA[iType] = new TProfile3D(Form("fProfile3DDiffGammaLambdaProtonMassDCA_%i",iType), Form("fProfile3DDiffGammaLambdaProtonMassDCA_%i",iType), 7, 0, 70, 4, 0, 4, fNMassBins, fLambdaMassMean - fLambdaMassLeftCut, fLambdaMassMean + fLambdaMassRightCut);
    fResultsList->Add(fHist3LambdaProtonMassDCA[iType]);
    fResultsList->Add(fProfile3DDiffDeltaLambdaProtonMassDCA[iType]);
    fResultsList->Add(fProfile3DDiffGammaLambdaProtonMassDCA[iType]);
  }

  PostData(2,fResultsList);
  if (fDebug) Printf("Post fResultsList Data Success!");
}

//------------------------------------------------

void AliAnalysisTaskCVEPIDCMEDiff::UserExec(Option_t *)
{
  if (fDebug) Printf("===============================We are in UserExec!!!===================================");
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

  fMultSel = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
  if (!fMultSel) {
    AliError(Form("%s: Could not get AliMultSelection", GetName()));
  } else fEvtCount->Fill(7);

  if (!manager || !handler || !fAOD || !fPIDResponse || !fUtils || !fMultSel) return;
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
  fHistVz[0]->Fill(fVertex[2]);
  if (fabs(fVertex[2]) > fVzCut) return;
  if (!fVtx || fVtx->GetNContributors() < 2 || vtSPD->GetNContributors()<1) return;
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
  fCentV0M  = fMultSel->GetMultiplicityPercentile("V0M");
  fCentTRK  = fMultSel->GetMultiplicityPercentile("TRK");
  fCentSPD0 = fMultSel->GetMultiplicityPercentile("CL0");
  fCentSPD1 = fMultSel->GetMultiplicityPercentile("CL1");

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

  //PF-Preview comment
  if(fCent > 30 && fCent < 50) {
    if (!(mask & (AliVEvent::kINT7 + AliVEvent::kSemiCentral))) return;
  } else {
    if (!(mask & AliVEvent::kINT7)) return;
  }

  // cent bin
  fCentBin = (int)fCent/10;
  fHistCent[0]->Fill(fCent);
  fEvtCount->Fill(12);
  if (fDebug) Printf("centrality done!");

  //----------------------------
  // Pile up
  //----------------------------
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) {
    if (!RejectEvtTFFit()) return;
  }
  fHistCent[1]->Fill(fCent);
  fEvtCount->Fill(13);
  if (fDebug) Printf("pile-up done!");
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
  if (fPlaneEstimator.EqualTo("TPC")) {
    fPsi2 = GetTPCPlane();
  }
  else if (fPlaneEstimator.EqualTo("V0C")) {
    fPsi2 = GetV0CPlane();
  } else {
    std::cout<<"Wrong fPlaneEstimator!"<<std::endl;
    return;
  }
  if (TMath::IsNaN(fPsi2)) return;
  fEvtCount->Fill(18);
  if (fDebug) Printf("Get Plane done!");
  //----------------------------
  // Fill Resolution
  //----------------------------
  fHist2Psi2 -> Fill(fCent, fPsi2);
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
  //------------------
  // Post output data.
  //------------------
  PostData(1,fQAList);
  PostData(2,fResultsList);
  if (fDebug) Printf("analysis done!");
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::LoopTracks()
{
  int nTrks = fAOD->GetNumberOfTracks();
  if (nTrks < 10) return false;
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
    if (fPlaneEstimator.EqualTo("TPC")) {
      double wAcc = GetNUACor(charge, phi, eta, fVertex[2]);
      if (wAcc < 0) continue;
      else weight *= wAcc;
      fHistPhi[1]->Fill(phi, wAcc);
    }

    //DCA Cut
    double dcaxy = -999, dcaz = -999;
    if(!GetDCA(dcaxy,dcaz,track)) continue;
    // if FB = 96 or 768, we don't need special DCA cut

    if (pt > fPlanePtMin && pt < fPlanePtMax) {
      //Do we need to set pT as weight for Better resolution?
      fSumQ2xTPC  += weight * TMath::Cos(2 * phi);
      fSumQ2yTPC  += weight * TMath::Sin(2 * phi);
      fWgtMultTPC += weight;
      if (fPlaneEstimator.EqualTo("TPC")) {
        std::vector<double> vec_phi_weight;
        vec_phi_weight.emplace_back(phi);
        vec_phi_weight.emplace_back(weight);
        mapTPCTrksIDPhiWgt[id] = vec_phi_weight;
      }
    }

    // but we need to set the dca cut for 768 when we start to choose the paiticle for pair
    if (fabs(dcaxy) > fDcaCutXY) continue;
    if (fabs(dcaz) > fDcaCutZ) continue;
    
    fHistDcaXY->Fill(fabs(dcaxy));
    fHistDcaZ->Fill(fabs(dcaz));
    if (fFilterBit == 768 && isNarrowDcaCuts768) { 
      if (fabs(dcaz) > 2.0) continue;
      if (fabs(dcaxy) > 7.0 * (0.0026 + 0.005/TMath::Power(pt, 1.01))) continue;
    }

    bool isItProttrk = CheckPIDofParticle(track,3); // 3=proton
    
    // Proton Maybe need a customized DCA cut
    if(isProtonCustomizedDCACut) isItProttrk = isItProttrk && (fabs(dcaz) < 1. && fabs(dcaxy) < (0.0105 + 0.035/TMath::Power(pt,1.1))); 

    int code = 0;
    double pid_weight = 1.;
    if(charge > 0) {
      code = 999;
      isItProttrk = isItProttrk && (pt < fProtonPtMax && pt > fProtonPtMin);
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
        if(isDoNUE) {
          double nue_pid_weight = GetPIDNUECor(2212,pt);
          if (nue_pid_weight > 0) pid_weight *= nue_pid_weight;
        }
      }
    } else {
      code =-999;
      isItProttrk = isItProttrk && (pt < fAntiProtonPtMax && pt > fAntiProtonPtMin);
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
        if(isDoNUE) {
          double nue_pid_weight = GetPIDNUECor(-2212,pt);
          if (nue_pid_weight > 0) pid_weight *= nue_pid_weight;
        }
      }
    }

    vecParticle.emplace_back(std::array<double,8>{pt,eta,phi,(double)id,(double)code,weight,pid_weight,dcaz});
  }

  if(fabs(fSumQ2xTPC)<1.e-6 || fabs(fSumQ2yTPC)<1.e-6 || fWgtMultTPC < 1.e-5) return false;
  return true;
}

//---------------------------------------------------

double AliAnalysisTaskCVEPIDCMEDiff::GetTPCPlane()
{
  double psi2TPC = nan("");
  psi2TPC = GetEventPlane(fSumQ2xTPC,fSumQ2yTPC,2);
  return psi2TPC;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::LoopV0s()
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
    double dcaNeg  = v0->DcaNegToPrimVertex();
    double dcaPos  = v0->DcaPosToPrimVertex();
    fHistV0Pt              -> Fill(pt);
    fHistV0Eta             -> Fill(eta);
    fHistV0DcatoPrimVertex -> Fill(dcaToPV);
    fHistV0CPA             -> Fill(CPA);
    fHistV0DecayLength     -> Fill(dl);
    fHistV0NegDaughterDca  -> Fill(dcaNeg);
    fHistV0PosDaughterDca  -> Fill(dcaPos);
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
    double rap = v0->RapLambda();

    if (code == 3122) {
      double mass  = v0->MassLambda();
      fHistLambdaPt[0]              -> Fill(pt);
      fHistLambdaEta[0]             -> Fill(eta);
      fHistLambdaPhi[0]             -> Fill(phi);
      fHistLambdaDcaToPrimVertex[0] -> Fill(dcaToPV);
      fHistLambdaNegDaugtherDca[0]  -> Fill(nDcaPV);
      fHistLambdaPosDaugtherDca[0]  -> Fill(pDcaPV);
      fHistLambdaCPA[0]             -> Fill(CPA);
      fHistLambdaDecayLength[0]     -> Fill(dl);
      fHist3LambdaCentPtMass[0]     -> Fill(fCent, pt, mass);
      fHist2LambdaMassPtY[0]        -> Fill(pt, rap);

      bool isLambda = true;
      isLambda = isLambda && (pt > fLambdaPtMin && pt < fLambdaPtMax);
      isLambda = isLambda && (fabs(rap) < fLambdaRapidityMax);

      isLambda = isLambda && (mass > fLambdaMassMean - fLambdaMassLeftCut);
      isLambda = isLambda && (mass < fLambdaMassMean + fLambdaMassRightCut);

      if (isLambda) {
        // if a daughter has been used as a daughter from Lambda before(It happends), we have to refuse a new one.
        // [pt,eta,phi,id,pdgcode,weight,mass,id1,id2]
        bool found = std::any_of(vecParticleV0.begin(), vecParticleV0.end(), [id_daughter_1,id_daughter_2](std::array<double, 9>& arr){
            return ((int)arr[7] == id_daughter_1 || (int)arr[8] == id_daughter_2 || (int)arr[7] == id_daughter_2 || (int)arr[8] == id_daughter_1);
        });
        if (found) continue;
        fHistLambdaPt[1]              -> Fill(pt);
        fHistLambdaEta[1]             -> Fill(eta);
        fHistLambdaPhi[1]             -> Fill(phi);
        fHistLambdaDcaToPrimVertex[1] -> Fill(dcaToPV);
        fHistLambdaNegDaugtherDca[1]  -> Fill(nDcaPV);
        fHistLambdaPosDaugtherDca[1]  -> Fill(pDcaPV);
        fHistLambdaCPA[1]             -> Fill(CPA);
        fHistLambdaDecayLength[1]     -> Fill(dl);
        fHist3LambdaCentPtMass[1]     -> Fill(fCent, pt, mass);
        fHist2LambdaMassPtY[1]        -> Fill(pt, rap);

        double weight = 1.;
        if(isDoLambdaNUE) {
          double nue_eight = GetPIDNUECor(3122,pt);
          if (nue_eight > 0.) weight *= nue_eight;
        }

        vecParticleV0.emplace_back(std::array<double,9>{pt,eta,phi,0.,(double)code,weight,mass,(double)id_daughter_1,(double)id_daughter_2});
      }
    }

    if (code == -3122) {
      double mass  = v0->MassAntiLambda();
      fHistAntiLambdaPt[0]              -> Fill(pt);
      fHistAntiLambdaEta[0]             -> Fill(eta);
      fHistAntiLambdaPhi[0]             -> Fill(phi);
      fHistAntiLambdaDcaToPrimVertex[0] -> Fill(dcaToPV);
      fHistAntiLambdaNegDaugtherDca[0]  -> Fill(nDcaPV);
      fHistAntiLambdaPosDaugtherDca[0]  -> Fill(pDcaPV);
      fHistAntiLambdaCPA[0]             -> Fill(CPA);
      fHistAntiLambdaDecayLength[0]     -> Fill(dl);
      fHist3AntiLambdaCentPtMass[0]     -> Fill(fCent, pt, mass);
      fHist2AntiLambdaMassPtY[0]        -> Fill(pt, rap);

      bool isAntiLambda = true;
      isAntiLambda = isAntiLambda && (pt > fAntiLambdaPtMin && pt < fAntiLambdaPtMax);
      isAntiLambda = isAntiLambda && (fabs(rap) < fAntiLambdaRapidityMax);

      isAntiLambda = isAntiLambda && (mass > fLambdaMassMean - fAntiLambdaMassLeftCut);
      isAntiLambda = isAntiLambda && (mass < fLambdaMassMean + fAntiLambdaMassRightCut);

      if (isAntiLambda) {
        // if a daughter has been used as a daughter from Lambda before(It happends), we have to refuse a new one.
        // [pt,eta,phi,id,pdgcode,weight,mass,id1,id2]
        bool found = std::any_of(vecParticleV0.begin(), vecParticleV0.end(), [id_daughter_1,id_daughter_2](std::array<double, 9>& arr){
            return ((int)arr[7] == id_daughter_1 || (int)arr[8] == id_daughter_2 || (int)arr[7] == id_daughter_2 || (int)arr[8] == id_daughter_1);
        });
        if (found) continue;
        fHistAntiLambdaPt[1]              -> Fill(pt);
        fHistAntiLambdaEta[1]             -> Fill(eta);
        fHistAntiLambdaPhi[1]             -> Fill(phi);
        fHistAntiLambdaDcaToPrimVertex[1] -> Fill(dcaToPV);
        fHistAntiLambdaNegDaugtherDca[1]  -> Fill(nDcaPV);
        fHistAntiLambdaPosDaugtherDca[1]  -> Fill(pDcaPV);
        fHistAntiLambdaCPA[1]             -> Fill(CPA);
        fHistAntiLambdaDecayLength[1]     -> Fill(dl);
        fHist3AntiLambdaCentPtMass[1]     -> Fill(fCent, pt, mass);
        fHist2AntiLambdaMassPtY[1]        -> Fill(pt, rap);

        double weight = 1.;
        if(isDoLambdaNUE) {
          double nue_eight = GetPIDNUECor(-3122,pt);
          if (nue_eight > 0.) weight *= nue_eight;
        }

        vecParticleV0.emplace_back(std::array<double,9>{pt,eta,phi,0.,(double)code,weight,mass,(double)id_daughter_1,(double)id_daughter_2});
      }
    }
  }//loop V0 end
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::PairV0Trk()
{
  //Lambda - X
  for (auto lambda : vecParticleV0) {
    double pt_lambda     = lambda[0];
    double eta_lambda    = lambda[1];
    double phi_lambda    = lambda[2];
    // int    id_lambda     = (int)lambda[3];
    int    code_lambda   = (int)lambda[4];
    double weight_lambda = lambda[5];
    double mass_lambda   = lambda[6];
    int    id_daughter_1 = (int)lambda[7];
    int    id_daughter_2 = (int)lambda[8];

    if(code_lambda == 3122) {
      fHist2LambdaMass[0] -> Fill(fCent, mass_lambda);
      fHist2LambdaMass[1] -> Fill(fCent, mass_lambda);
    } else if (code_lambda == -3122) {
      fHist2LambdaMass[2] -> Fill(fCent, mass_lambda);
      fHist2LambdaMass[3] -> Fill(fCent, mass_lambda);
    } else {
      AliError("this is not Lambda or AntiLambda");
      return false;
    }

    for (auto particle : vecParticle) {
      double pt     = particle[0];
      double eta    = particle[1];
      double phi    = particle[2];
      int    id     = (int)particle[3];
      int    code   = (int)particle[4];
      double weight = particle[5];
      double pidweight = particle[6];
      double  dcaz  = particle[7];
      if (id == id_daughter_1 || id == id_daughter_2) continue;

      double psi2_forThisPair = nan("");
      if (fPlaneEstimator.EqualTo("TPC")) {
        psi2_forThisPair = GetTPCPlaneNoAutoCorr({id_daughter_1, id_daughter_2, id});
      } else if (fPlaneEstimator.EqualTo("V0C")) {
        psi2_forThisPair = fPsi2;
      } else {
        AliError("fPlaneEstimator is not set properly.");
        return false;
      }

      if (std::isnan(psi2_forThisPair)) continue;
      double delta = TMath::Cos(phi_lambda - phi);
      double gamma = TMath::Cos(phi_lambda + phi - 2 * psi2_forThisPair);

      int nBits = 4;
      if(!isCalculateLambdaHadron) {
        TBits bitsLambdaProtonPair(nBits);
        bitsLambdaProtonPair.SetBitNumber(0, code_lambda ==  3122 && code ==  2212);
        bitsLambdaProtonPair.SetBitNumber(1, code_lambda ==  3122 && code == -2212);
        bitsLambdaProtonPair.SetBitNumber(2, code_lambda == -3122 && code ==  2212);
        bitsLambdaProtonPair.SetBitNumber(3, code_lambda == -3122 && code == -2212);

        double sumPtBin = GetSumPtBin(pt_lambda + pt);
        double deltaEtaBin = GetDeltaEtaBin(TMath::Abs(eta_lambda - eta));
        double dcaBin = GetDCABin(dcaz);

        for (int iBits = 0; iBits < nBits; iBits++) {
          double weight_all = weight_lambda * pidweight;
           if (bitsLambdaProtonPair.TestBitNumber(iBits)) {
            fHist3LambdaProtonMassSPt[iBits]               -> Fill(fCent,    sumPtBin, mass_lambda);
            fProfile3DDiffDeltaLambdaProtonMassSPt[iBits]  -> Fill(fCent,    sumPtBin, mass_lambda, delta, weight_all);
            fProfile3DDiffGammaLambdaProtonMassSPt[iBits]  -> Fill(fCent,    sumPtBin, mass_lambda, gamma, weight_all);
  
            fHist3LambdaProtonMassDEta[iBits]              -> Fill(fCent, deltaEtaBin, mass_lambda);
            fProfile3DDiffDeltaLambdaProtonMassDEta[iBits] -> Fill(fCent, deltaEtaBin, mass_lambda, delta, weight_all);
            fProfile3DDiffGammaLambdaProtonMassDEta[iBits] -> Fill(fCent, deltaEtaBin, mass_lambda, gamma, weight_all);

            fHist3LambdaProtonMassDCA[iBits]               -> Fill(fCent,    dcaBin,   mass_lambda);
            fProfile3DDiffDeltaLambdaProtonMassDCA[iBits]  -> Fill(fCent,    dcaBin,   mass_lambda, delta, weight_all);
            fProfile3DDiffGammaLambdaProtonMassDCA[iBits]  -> Fill(fCent,    dcaBin,   mass_lambda, gamma, weight_all);
          }
        }
      }

      if(isCalculateLambdaHadron) {
        TBits bitsLambdaHadronPair(nBits);
        bitsLambdaHadronPair.SetBitNumber(0, code_lambda ==  3122 && code > 0);
        bitsLambdaHadronPair.SetBitNumber(1, code_lambda ==  3122 && code < 0);
        bitsLambdaHadronPair.SetBitNumber(2, code_lambda == -3122 && code > 0);
        bitsLambdaHadronPair.SetBitNumber(3, code_lambda == -3122 && code < 0);
        for (int iBits = 0; iBits < nBits; iBits++) {
          double weight_all = weight_lambda * weight;
          if (bitsLambdaHadronPair.TestBitNumber(iBits)) {
            fProfile2DDeltaLambdaHadronMass[iBits]  -> Fill(fCent, mass_lambda, delta, weight_all);
            fProfile2DGammaLambdaHadronMass[iBits]  -> Fill(fCent, mass_lambda, gamma, weight_all);
          }
        }
      }
    } 
  }
  return true;
}


//---------------------------------------------------

void AliAnalysisTaskCVEPIDCMEDiff::ResetVectors()
{
  fSumQ2xTPC = 0.;
  fSumQ2yTPC = 0.;
  fWgtMultTPC = 0.;
  std::unordered_map<int, std::vector<double>>().swap(mapTPCTrksIDPhiWgt);
  std::vector<std::array<double,8>>().swap(vecParticle);
  std::vector<std::array<double,9>>().swap(vecParticleV0);
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::LoadCalibHistForThisRun()
{
  //18q/r NUA
  if (fPlaneEstimator.EqualTo("TPC")) {
    hCorrectNUAPos -> Reset();
    hCorrectNUANeg -> Reset();
    hCorrectNUAPos = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dPos_Run%d",0,fRunNum));
    hCorrectNUANeg = (TH3F*) fListNUA->FindObject(Form("fHist_NUA_VzPhiEta_kPID%dNeg_Run%d",0,fRunNum));
    if (!hCorrectNUAPos) return false;
    if (!hCorrectNUANeg) return false;
  }
  //V0C
  if (fPlaneEstimator.EqualTo("V0C")) {
    // for recenter
    hQx2mV0C ->Reset();
    hQy2mV0C ->Reset();
    hQx2mV0C = ((TH1D*) contQxncm->GetObject(fRunNum));
    hQy2mV0C = ((TH1D*) contQyncm->GetObject(fRunNum));
    // for gain equalization
    fHCorrectV0ChWeghts -> Reset();
    fHCorrectV0ChWeghts = (TH2F *) fListVZEROCalib->FindObject(Form("hWgtV0ChannelsvsVzRun%d",fRunNum));
    if (!hQx2mV0C) return false;
    if (!hQy2mV0C) return false;
    if (!fHCorrectV0ChWeghts) return false;
  }
  return true;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::RejectEvtTFFit()
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

bool AliAnalysisTaskCVEPIDCMEDiff::AcceptAODTrack(AliAODTrack *track)
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

bool AliAnalysisTaskCVEPIDCMEDiff::CheckPIDofParticle(AliAODTrack* ftrack, int pidToCheck)
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

    bool isProton = TMath::Abs(nSigRMS) < fNSigmaTOFCut;

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

double AliAnalysisTaskCVEPIDCMEDiff::GetNUECor(int charge, double pt)
{
  if(!fListNUE) return -1;
  if(charge == 0) return -1;
  // TString histName = (charge > 0) ? "h_eff_pos_hadron" : "h_eff_neg_hadron";
  // TH1D* efficiencyHist = (TH1D*)fListNUE->FindObject(histName);
  // if (!efficiencyHist) return -1;
  // int ptBin = efficiencyHist->GetXaxis()->FindBin(pt);
  // double binContent = efficiencyHist->GetBinContent(ptBin);

  // if (binContent > 1.e-5) {
  //   return 1.0 / binContent;
  // } else return -1;
  double binContent = -1;
  TString histName = (charge > 0) ? "h2_eff_pos_hadron_dcazcut" : "h2_eff_neg_hadron_dcazcut";
  TH2D* h2_eff = (TH2D*)fListNUE->FindObject(histName);
  if (!h2_eff) return -1;
  double dcazcut = fDcaCutZ;
  if (dcazcut > 1.0) dcazcut = 1.0;
  dcazcut = dcazcut - 1.e-6;
  binContent = h2_eff->GetBinContent(h2_eff->FindBin(dcazcut, pt));

  if (binContent > 1.e-6) {
    return 1.0 / binContent;
  } else return -1;
}


//---------------------------------------------------
double AliAnalysisTaskCVEPIDCMEDiff::GetPIDNUECor(int pdgcode, double pt)
{
  if(!fListNUE) return -1;
  TString histName;

  // if (pdgcode == 2212)       histName = "h_eff_pos_proton";
  // else if (pdgcode == -2212) histName = "h_eff_neg_proton";
  // if (pdgcode == 3122)       histName = "h_eff_lambda";
  // else if (pdgcode == -3122) histName = "h_eff_antilambda";
  // else return -1;

  // TH1D* efficiencyHist = (TH1D*)fListNUE->FindObject(histName);
  // if (!efficiencyHist) return -1;
  // int ptBin = efficiencyHist->GetXaxis()->FindBin(pt);
  // double binContent = efficiencyHist->GetBinContent(ptBin);

  // if (binContent > 1.e-5) {
  //   return 1.0 / binContent;
  // } else return -1;

  double binContent = -1;
  if(abs(pdgcode) == 2212) {
    pdgcode == 2212 ? histName = "h2_eff_pos_proton_dcazcut" : "h2_eff_neg_proton_dcazcut";
    TH2D* h2_eff = (TH2D*)fListNUE->FindObject(histName);
    if (!h2_eff) return -1;
    double dcazcut = fDcaCutZ;
    if (dcazcut > 1.0) dcazcut = 1.0;
    dcazcut = dcazcut - 1.e-6;
    binContent = h2_eff->GetBinContent(h2_eff->FindBin(dcazcut, pt));
  } else if (abs(pdgcode == 3122)) {
    pdgcode == 3122 ? histName = "h_eff_lambda" : "h_eff_antilambda";
    TH1D* h_eff = (TH1D*)fListNUE->FindObject(histName);
    if (!h_eff) return -1;
    binContent = h_eff->GetBinContent(h_eff->FindBin(pt));
  } else return -1;

  if (binContent > 1.e-6) {
  return 1.0 / binContent;
  } else return -1;
}


//---------------------------------------------------

double AliAnalysisTaskCVEPIDCMEDiff::GetNUACor(int charge, double phi, double eta, double vz)
{
  double weightNUA = 1;
  if (fVzBin<0 || fCentBin<0 || fRunNum<0) return -1;
  if (fPeriod.EqualTo("LHC18q") || fPeriod.EqualTo("LHC18r")) { // Rihan and Protty 's NUA Results
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

bool AliAnalysisTaskCVEPIDCMEDiff::IsGoodV0(AliAODv0 *aodV0)
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
  return kTRUE;
}

//---------------------------------------------------

bool AliAnalysisTaskCVEPIDCMEDiff::IsGoodDaughterTrack(const AliAODTrack *track)
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
  if (nCrossedRowsTPC/findable < fRatioCrossedRowsFindable) return false;
  return true;
}

//---------------------------------------------------

int AliAnalysisTaskCVEPIDCMEDiff::GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *nTrack)
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

inline double AliAnalysisTaskCVEPIDCMEDiff::GetEventPlane(double qx, double qy, double harmonic)
{
  double psi = (1./harmonic)*TMath::ATan2(qy,qx);
  if (psi < 0) return psi += TMath::TwoPi()/harmonic;
  else return psi;
}

//---------------------------------------------------
double AliAnalysisTaskCVEPIDCMEDiff::GetTPCPlaneNoAutoCorr(std::vector<int> vec_id)
{
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

bool AliAnalysisTaskCVEPIDCMEDiff::GetDCA(double &dcaxy, double &dcaz, AliAODTrack* track)
{
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
inline double AliAnalysisTaskCVEPIDCMEDiff::GetSumPtBin(double sumPt)
{
  //1~3 GeV return 0.5
  //3~5 GeV return 1.5
  //5~8 GeV return 2.5
  //else return -1
  if (sumPt > 1. && sumPt < 3.) return 0.5;
  else if (sumPt > 3. && sumPt < 5.) return 1.5;
  else if (sumPt > 5. && sumPt < 8.) return 2.5;
  else return -1.;
}

//---------------------------------------------------
inline double AliAnalysisTaskCVEPIDCMEDiff::GetDeltaEtaBin(double deltaEta)
{
  //-0.6 ~ 0.6 return 0.5
  //-1.6 ~ -0.6 & 0.6 ~ 1.6 return 1.5
  //else return -1
  if (deltaEta > -0.6 && deltaEta < 0.6) return 0.5;
  else if ((deltaEta > -1.6 && deltaEta < -0.6) || (deltaEta > 0.6 && deltaEta < 1.6)) return 1.5;
  else return -1.;
}

//---------------------------------------------------
inline double AliAnalysisTaskCVEPIDCMEDiff::GetDCABin(double dca)
{
  dca = abs(dca);
  //0, 0.002 0.005 0.01
  if (dca > 0. && dca < 0.005) return 0.5;
  else if (dca > 0.005 && dca < 0.01) return 1.5;
  else if (dca > 0.01 && dca < 0.1)  return 2.5;
  else if (dca > 0.1) return 3.5;
  else return -1.;
}

//---------------------------------------------------

double AliAnalysisTaskCVEPIDCMEDiff::GetV0CPlane()
{

  double psi2V0C = nan("");
  double Qx2 = 0, Qy2 = 0, Mult = 0;

  //Loop Over VZERO Channels
  //Gain Equalization
  for (int iCh = 0; iCh < 32; ++iCh) { //just need C side
    //phi angle of the channel
    double phi = TMath::Pi()/8. + TMath::Pi()/4.*(iCh%8);
    //energy of the channel
    double energy_raw = 0.;
    AliAODVZERO* aodV0 = fAOD->GetVZEROData();
    energy_raw = aodV0->GetMultiplicity(iCh);

    //get the gain equlization factor
    int ibinV0 = fHCorrectV0ChWeghts->FindBin(fVertex[2],iCh);
    double factor_gain = (double)fHCorrectV0ChWeghts->GetBinContent(ibinV0);

    double energy_gain = energy_raw * factor_gain;
    if (energy_gain < 1.e-6) return nan("");

    Qx2 += energy_gain * TMath::Cos(2*phi);
    Qy2 += energy_gain * TMath::Sin(2*phi);
    Mult += energy_gain;
  }
  if (Mult < 1.e-6) return nan("");

  //VZERO Recenter
  // get the mean Q vector from input file
  int iCentSPD = (int)fCentSPD1;
  double Qx2Mean = hQx2mV0C->GetBinContent(iCentSPD+1);
  double Qy2Mean = hQy2mV0C->GetBinContent(iCentSPD+1);
  if (Qy2Mean < -900 || Qx2Mean < -900) return nan("");

  // recenter the Q vector
  Qx2 -= Qx2Mean;
  Qy2 -= Qy2Mean;

  psi2V0C = GetEventPlane(Qx2, Qy2, 2.);

  return psi2V0C;
}

//---------------------------------------------------
