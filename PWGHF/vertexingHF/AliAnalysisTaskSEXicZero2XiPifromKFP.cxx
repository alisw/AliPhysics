/**************************************************************************
 * Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////
// Author: Jianhui Zhu (1,2)
// (1) Central China Normal University
// (2) GSI Helmholtz Centre for Heavy Ion Research
// E-mail: zjh@mail.ccnu.edu.cn
/////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <TDatabasePDG.h>
#include <vector>
#include <TVector3.h>
#include "TChain.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLine.h"
#include "TList.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliNormalizationCounter.h"
#include "AliAnalysisTaskSEXicZero2XiPifromKFP.h"
#include "AliPIDResponse.h"
#include "AliVertexingHFUtils.h"

#include "AliAODMCParticle.h"

// includes added to play with KFParticle
#ifndef HomogeneousField
#define HomogeneousField 
#endif

using std::cout;
using std::endl;

class AliAnalysisTaskSEXicZero2XiPifromKFP;    // your analysis class

ClassImp(AliAnalysisTaskSEXicZero2XiPifromKFP) // classimp: necessary for root

std::ostream&  operator<<(std::ostream& os, const KFParticleBase& particle) {                            
  static const char *vn[14] = {"x","y","z","px","py","pz","E","S","M","t","p","Q","Chi2","NDF"};

  for (Int_t i = 0; i < 8; i++) {
    if (i == 6) continue;                                    // E
    if (i == 7 && particle.GetParameter(i) <= 0.0) continue; // S
    if (particle.GetParameter(i) == 0. && particle.GetCovariance(i,i) == 0) continue;
    if (particle.GetCovariance(i,i) > 0) 
      os << " " << vn[i]<<": "<< std::setw(8) << particle.GetParameter(i)<< " +/- " << std::setw(6) << sqrt(particle.GetCovariance(i,i));
    else 
      os << " " << vn[i] << ": " << std::setw(8) << particle.GetParameter(i);
  }
  float Mtp[3], MtpErr[3];
  particle.GetMass(Mtp[0], MtpErr[0]);     if (MtpErr[0] < 1e-7 || MtpErr[0] > 1e10) MtpErr[0] = -13;
  particle.GetLifeTime(Mtp[1], MtpErr[1]); if (MtpErr[1] <=   0 || MtpErr[1] > 1e10) MtpErr[1] = -13;
  particle.GetMomentum(Mtp[2], MtpErr[2]); if (MtpErr[2] <=   0 || MtpErr[2] > 1e10) MtpErr[2] = -13;
  for (Int_t i = 8; i < 11; i++) {
    if (i == 9 && Mtp[i-8] <= 0.0) continue; // t
    if (MtpErr[i-8] > 0 && MtpErr[i-8] < 1e10) os << " " << vn[i] << ": " << std::setw(8) << Mtp[i-8] << " +/-" << std::setw(7) << MtpErr[i-8];
    else                                       os << " " << vn[i] << ": " << std::setw(8) << Mtp[i-8];
  }
  os << " pdg:" << std::setw(5) << particle.GetPDG() << " Q: "<< std::setw(2) << int(particle.GetQ()) << " chi2/NDF: " << std::setw(8) << particle.GetChi2() << "/" << std::setw(2) << particle.GetNDF();
  return os;
}

AliAnalysisTaskSEXicZero2XiPifromKFP::AliAnalysisTaskSEXicZero2XiPifromKFP() :
  AliAnalysisTaskSE(),
  fIsMC(kFALSE),
  fPID(0),
  fAnaCuts(0),
  fpVtx(0),
  fMCEvent(0),
  fBzkG(0),
  fCentrality(0),
//  fRunNumber(0),
//  fEvNumberCounter(0),
  fAodTrackInd(0),
  fOutputList(0),
  fOutputWeight(0),
  fListCuts(0),
  fTree_Event(0),
  fVar_Event(0),
  fTree_Xic0(0),
  fVar_Xic0(0),
  fTree_Xic0MCGen(0),
  fVar_Xic0MCGen(0),
  fCounter(0),
  fHistMCGen_Lambda_Pt(0),
  fHistMCGen_AntiLambda_Pt(0),
  fHistMCGen_Lambda_Pt_wYcut(0),
  fHistMCGen_AntiLambda_Pt_wYcut(0),
  fCounterGen_Cuts_Lambda(0),
  fCounterGen_Cuts_AntiLambda(0),
  fCounterRecMC_Cuts_Lambda(0),
  fCounterRecMC_Cuts_AntiLambda(0),
  fCounterRec_Cuts_Lambda(0),
  fCounterRec_Cuts_AntiLambda(0),
  fHistMCGen_XiMinus_Pt(0),
  fHistMCGen_XiPlus_Pt(0),
  fHistMCGen_XiMinus_Pt_wYcut(0),
  fHistMCGen_XiPlus_Pt_wYcut(0),
  fCounterGen_Cuts_XiMinus(0),
  fCounterGen_Cuts_XiPlus(0),
  fCounterRecMC_Cuts_XiMinus(0),
  fCounterRecMC_Cuts_XiPlus(0),
  fCounterRec_Cuts_XiMinus(0),
  fCounterRec_Cuts_XiPlus(0),
  f2DCounterRecMC_CutsVsPt_Lambda(0),
  f2DCounterRecMC_CutsVsPt_AntiLambda(0),
  f2DCounterRecMC_CutsVsPt_XiMinus(0),
  f2DCounterRecMC_CutsVsPt_XiPlus(0),
  f2DHistMassPtLambda(0),
  f2DHistMassPtAntiLambda(0),
  f2DHistMassPtXiMinus(0),
  f2DHistMassPtXiPlus(0),
  f2DHistMassPtXicZero_woQtCut(0),
  f2DHistMassPtXicZero(0),
  f2DHistMassPtAntiXicZero(0),
  f2DHistMassPtXicZeroTot(0),
  f2DHistChi2vsNDF_Lambda(0),
  f2DHistChi2vsNDF_Lambda_Match(0),
  f2DHistChi2vsNDF_AntiLambda(0),
  f2DHistXKFMCvsPt_Proton(0),
  f2DHistYKFMCvsPt_Proton(0),
  f2DHistZKFMCvsPt_Proton(0),
  f2DHistXKFMCvsPt_Pion(0),
  f2DHistYKFMCvsPt_Pion(0),
  f2DHistZKFMCvsPt_Pion(0),
  f2DHistXPULLvsPt_Proton(0),
  f2DHistYPULLvsPt_Proton(0),
  f2DHistZPULLvsPt_Proton(0),
  f2DHistXPULLvsPt_Pion(0),
  f2DHistYPULLvsPt_Pion(0),
  f2DHistZPULLvsPt_Pion(0),
  f2DHistXRecMCvsPt_Lambda(0),
  f2DHistXRecMCvsPt_AntiLambda(0),
  f2DHistXV0MCvsPt_Lambda(0),
  f2DHistPtRecMCvsPt_Lambda(0),
  f2DHistPtV0MCvsPt_Lambda(0),
  f2DHistPtRecMCvsPt_AntiLambda(0),
  f2DHistMassRecMCvsPt_Lambda(0),
  f2DHistMassV0MCvsPt_Lambda(0),
  f2DHistMassRecMCvsPt_AntiLambda(0),
  f2DHistXPULLvsPt_Lambda(0),
  f2DHistXPULLvsPt_Lambda_V0(0),
  f2DHistXPULLvsPt_AntiLambda(0),
  f2DHistPtPULLvsPt_Lambda(0),
  f2DHistPtPULLvsPt_Lambda_V0(0),
  f2DHistPtPULLvsPt_AntiLambda(0),
  f2DHistMassPULLvsPt_Lambda(0),
  f2DHistMassPULLvsPt_Lambda_V0(0),
  f2DHistMassPULLvsPt_AntiLambda(0),
  f2DHistMassPULLvsRadius_Lambda(0),
  f2DHistMassPULLvsRadius_AntiLambda(0),
  f2DHistArmenterosPodolanski_FirstDaugPos(0),
  f2DHistArmenterosPodolanski_FirstDaugNeg(0),
  f2DHistArmenterosPodolanski_candidate(0),
  f2DHistArmenterosPodolanski_Lam(0),
  f2DHistArmenterosPodolanski_AntiLam(0),
  f2DHistChargeDaughters(0),
  f2DHistV0XY_OnFly(0),
  f2DHistV0XY_Offline(0),
  f2DHistV0XY_FirstDaugPos(0),
  f2DHistV0XY_FirstDaugNeg(0),
  f2DHistLambdaXY(0),
  f2DHistXiMinusXY_DV(0),
  f2DHistXiMinusXY_PV(0),
  f2DHistXiPlusXY_DV(0),
  f2DHistXiPlusXY_PV(0),
  fHistEvents(0),
  fHTrigger(0),
  fHistOnFlyStatus(0),
  fHistOnFlyStatus_FirstDaugPos(0),
  fHistOnFlyStatus_FirstDaugNeg(0),
  fHistChargeV0(0),
  fHistNDaughterV0(0),
  fHistNProngV0(0),
  fHistChargeFirstDaughter(0),
  fHistChargeSecondDaughter(0),
  fHistXtrkP(0),
  fHistYtrkP(0),
  fHistZtrkP(0),
  fHistXtrkP_XYZv(0),
  fHistYtrkP_XYZv(0),
  fHistZtrkP_XYZv(0),
  fHistXtrkP_Rec_MC(0),
  fHistYtrkP_Rec_MC(0),
  fHistZtrkP_Rec_MC(0),
  fHistXtrkP_Rec_MC_XYZv(0),
  fHistYtrkP_Rec_MC_XYZv(0),
  fHistZtrkP_Rec_MC_XYZv(0),
  fHistXtrkN(0),
  fHistYtrkN(0),
  fHistZtrkN(0),
  fHistXtrkN_XYZv(0),
  fHistYtrkN_XYZv(0),
  fHistZtrkN_XYZv(0),
  fHistXtrkN_Rec_MC(0),
  fHistYtrkN_Rec_MC(0),
  fHistZtrkN_Rec_MC(0),
  fHistXtrkN_Rec_MC_XYZv(0),
  fHistYtrkN_Rec_MC_XYZv(0),
  fHistZtrkN_Rec_MC_XYZv(0),
  fHistLDeltaLRec_Lambda(0),
  fHistLDeltaLRecMC_Lambda(0),
  fHistLDeltaLRecMC_LambdaFromXi(0),
  fHistLDeltaLRec_AntiLambda(0),
  fHistLDeltaLRecMC_AntiLambda(0),
  fHistLDeltaLRecMC_AntiLambdaFromXi(0),
  fHistXLambdaTot(0),
  fHistYLambdaTot(0),
  fHistZLambdaTot(0),
  fHistXLambda_KF_MC(0),
  fHistYLambda_KF_MC(0),
  fHistZLambda_KF_MC(0),
  fHistXProton_KF_MC(0),
  fHistYProton_KF_MC(0),
  fHistZProton_KF_MC(0),
  fHistXPion_KF_MC(0),
  fHistYPion_KF_MC(0),
  fHistZPion_KF_MC(0),
  fHistXLambda_V0_MC(0),
  fHistXAntiLambda_Rec_MC(0),
  fHistYAntiLambda_Rec_MC(0),
  fHistZAntiLambda_Rec_MC(0),
  fHistXLambda_PULL(0),
  fHistYLambda_PULL(0),
  fHistZLambda_PULL(0),
  fHistXProton_PULL(0),
  fHistYProton_PULL(0),
  fHistZProton_PULL(0),
  fHistXPion_PULL(0),
  fHistYPion_PULL(0),
  fHistZPion_PULL(0),
  fHistXAntiLambda_PULL(0),
  fHistYAntiLambda_PULL(0),
  fHistZAntiLambda_PULL(0),
  fHistXXiTot(0),
  fHistYXiTot(0),
  fHistZXiTot(0),
  fHistXXicZeroTot(0),
  fHistYXicZeroTot(0),
  fHistZXicZeroTot(0),
  fGenHistRapidity_Lambda(0),
  fGenHistRapidity_AntiLambda(0),
  fRecHistRapidity_Lambda_offline(0),
  fRecHistRapidity_Lambda_wSTD(0),
  fRecHistRapidity_AntiLambda_wSTD(0),
  fHistPtLambda(0),
  fHistPtAntiLambda(0),
  fHistPtLambdaTot(0),
  fHistPtXiMinus(0),
  fHistPtXiPlus(0),
  fHistPtXiTot(0),
  fHistPtXicZero(0),
  fHistPtAntiXicZero(0),
  fHistPtXicZeroTot(0),
  fHistMassK0S(0),
  fHistMassLambda_woCut(0),
  fHistMassAntiLambda_woCut(0),
  fHistMassLambdaCheck(0),
  fHistMassAntiLambdaCheck(0),
  fHistMassLambda_wSTDv0Cut(0),
  fHistMassAntiLambda_wSTDv0Cut(0),
  fHistMassLambda_BeforeSecSel(0),
  fHistMassAntiLambda_BeforeSecSel(0),
  fHistMassLambda(0),
  fHistMassLambda_woArmenterosPodolanskiCut(0),
  fHistMassLambda_woMassCut(0),
  fHistMassAntiLambda(0),
  fHistMassAntiLambda_woMassCut(0),
  fHistMassLambdaTot(0),
  fHistMassLambda_Match(0),
  fHistMassAntiLambda_Match(0),
  fHistMassLambdaTot_Match(0),
  fHistMassLambda_V0(0),
  fHistMassAntiLambda_V0(0),
  fHistMassLambdaTot_V0(0),
  fHistMassLambda_KF_V0(0),
  fHistMassLambda_KF_MC(0),
  fHistMassLambda_V0_MC(0),
  fHistMassAntiLambda_KF_V0(0),
  fHistMassAntiLambda_KF_MC(0),
  fHistMassLambda_PULL_KF(0),
  fHistMassAntiLambda_PULL_KF(0),
  fHistMassLambda_M(0),
  fHistMassAntiLambda_M(0),
  fHistMassLambdaTot_M(0),
  fHistMassLambda_MV(0),
  fHistMassAntiLambda_MV(0),
  fHistMassLambdaTot_MV(0),
  fHistMassXiMinus(0),
  fHistMassXiMinus_M(0),
  fHistMassXiMinus_Match(0),
  fHistMassXiPlus(0),
  fHistMassXiPlus_M(0),
  fHistMassXiPlus_Match(0),
  fHistMassXiTot(0),
  fHistMassXicZero_woQtCut(0),
  fHistMassXicZero(0),
  fHistMassAntiXicZero(0),
  fHistMassXicZeroTot(0),
  fHistQtDiffPionXiMinus(0),
  fHistQtDiffPionXiPlus(0),
  fHistMassXiMinus2(0),
  fHistChi2ndfProton(0),
  fHistChi2ndfPion(0),
  fHistChi2ndfLambda(0),
  fHistChi2ndfLambda_Match(0),
  fHistChi2ndfAntiLambda(0),
  fHistChi2ndfAntiLambda_Match(0),
  fHistChi2ndfLambdaTot(0),
  fHistProbProton(0),
  fHistProbPion(0),
  fHistProbLambda(0),
  fHistProbLambda_chi2cut(0),
  fHistProbLambda_Match(0),
  fHistProbAntiLambda(0),
  fHistProbAntiLambda_chi2cut(0),
  fHistProbAntiLambda_Match(0),
  fHistProbLambdaTot(0),
  fHistProbXiMinus(0),
  fHistProbXiMinus_chi2cut(0),
  fHistProbXiPlus(0),
  fHistProbXiPlus_chi2cut(0),
  fHistProbXicZero(0),
  fHistProbXicZero_chi2cut(0),
  fHistProbAntiXicZero(0),
  fHistProbAntiXicZero_chi2cut(0),
  fHistChi2ndfXiMinus(0),
  fHistChi2ndfXiPlus(0),
  fHistChi2ndfXiTot(0),
  fHistChi2ndfXicZero(0),
  fHistChi2ndfAntiXicZero(0),
  fHistChi2ndfXicZeroTot(0),
  fHistDecayLLambda(0),
  fHistDecayLAntiLambda(0),
  fHistDecayLLambdaTot(0),
  fHistDecayLXiMinus(0),
  fHistDecayLXiPlus(0),
  fHistDecayLXiTot(0),
  fHistDecayLXicZero(0),
  fHistDecayLAntiXicZero(0),
  fHistDecayLXicZeroTot(0),
  fHistCosPA_Lambda(0),
  fHistCosPA_AntiLambda(0),
  fHistCosPA_LambdaTot(0),
  fHistCosPA_XiMinus(0),
  fHistCosPA_XiPlus(0),
  fHistCosPA_XiTot(0),
  fHistPVx(0),
  fHistPVy(0),
  fHistPVz(0),
  fHCentrality(0),
  fHistMCGen_XicZeroTot(0),
  fHistMCGen_XicZero(0),
  fHistMCGen_AntiXicZero(0),
  fHistMCGen_PionTot(0),
  fHistMCGen_PionPlus(0),
  fHistMCGen_PionMinus(0),
  fHistMCGen_XiTot(0),
  fHistMCGen_XiMinus(0),
  fHistMCGen_XiPlus(0),
  fHistMCGen_Lambda(0),
  fHistMCGen_AntiLambda(0),
  fHistMCGen_PiXiInvMass(0),
  fHistMCGen_PiXiMassvsPiPt(0),
  fHistMCGen_PiXiMassvsPiPt_PionPlus(0),
  fHistMCGen_PiXiMassvsPiPt_PionMinus(0),
  fHistMCXicZeroDecayType(0),
  fHistMCXiDecayType(0),
  fHistMCpdg_All(0),
  fHistMCpdg_Dau_XicZero(0),
  fHistMCpdg_Dau_XicPM(0),
  fWriteXic0Tree(kFALSE),
  fWriteXic0MCGenTree(kFALSE),
  fWeight(0),
  fHistMCGen_Xic0Pt_weight(0),
  f2DHistMCRec_Xic0Pt_weight(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskSEXicZero2XiPifromKFP::AliAnalysisTaskSEXicZero2XiPifromKFP(const char* name, AliRDHFCutsKFP* cuts) :
  AliAnalysisTaskSE(name),
  fIsMC(kFALSE),
  fPID(0),
  fAnaCuts(cuts),
  fpVtx(0),
  fMCEvent(0),
  fBzkG(0),
  fCentrality(0),
//  fRunNumber(0),
//  fEvNumberCounter(0),
  fAodTrackInd(0),
  fOutputList(0),
  fOutputWeight(0),
  fListCuts(0),
  fTree_Event(0),
  fVar_Event(0),
  fTree_Xic0(0),
  fVar_Xic0(0),
  fTree_Xic0MCGen(0),
  fVar_Xic0MCGen(0),
  fCounter(0),
  fHistMCGen_Lambda_Pt(0),
  fHistMCGen_AntiLambda_Pt(0),
  fHistMCGen_Lambda_Pt_wYcut(0),
  fHistMCGen_AntiLambda_Pt_wYcut(0),
  fCounterGen_Cuts_Lambda(0),
  fCounterGen_Cuts_AntiLambda(0),
  fCounterRecMC_Cuts_Lambda(0),
  fCounterRecMC_Cuts_AntiLambda(0),
  fCounterRec_Cuts_Lambda(0),
  fCounterRec_Cuts_AntiLambda(0),
  fHistMCGen_XiMinus_Pt(0),
  fHistMCGen_XiPlus_Pt(0),
  fHistMCGen_XiMinus_Pt_wYcut(0),
  fHistMCGen_XiPlus_Pt_wYcut(0),
  fCounterGen_Cuts_XiMinus(0),
  fCounterGen_Cuts_XiPlus(0),
  fCounterRecMC_Cuts_XiMinus(0),
  fCounterRecMC_Cuts_XiPlus(0),
  fCounterRec_Cuts_XiMinus(0),
  fCounterRec_Cuts_XiPlus(0),
  f2DCounterRecMC_CutsVsPt_Lambda(0),
  f2DCounterRecMC_CutsVsPt_AntiLambda(0),
  f2DCounterRecMC_CutsVsPt_XiMinus(0),
  f2DCounterRecMC_CutsVsPt_XiPlus(0),
  f2DHistMassPtLambda(0),
  f2DHistMassPtAntiLambda(0),
  f2DHistMassPtXiMinus(0),
  f2DHistMassPtXiPlus(0),
  f2DHistMassPtXicZero_woQtCut(0),
  f2DHistMassPtXicZero(0),
  f2DHistMassPtAntiXicZero(0),
  f2DHistMassPtXicZeroTot(0),
  f2DHistChi2vsNDF_Lambda(0),
  f2DHistChi2vsNDF_Lambda_Match(0),
  f2DHistChi2vsNDF_AntiLambda(0),
  f2DHistXKFMCvsPt_Proton(0),
  f2DHistYKFMCvsPt_Proton(0),
  f2DHistZKFMCvsPt_Proton(0),
  f2DHistXKFMCvsPt_Pion(0),
  f2DHistYKFMCvsPt_Pion(0),
  f2DHistZKFMCvsPt_Pion(0),
  f2DHistXPULLvsPt_Proton(0),
  f2DHistYPULLvsPt_Proton(0),
  f2DHistZPULLvsPt_Proton(0),
  f2DHistXPULLvsPt_Pion(0),
  f2DHistYPULLvsPt_Pion(0),
  f2DHistZPULLvsPt_Pion(0),
  f2DHistXRecMCvsPt_Lambda(0),
  f2DHistXRecMCvsPt_AntiLambda(0),
  f2DHistXV0MCvsPt_Lambda(0),
  f2DHistPtRecMCvsPt_Lambda(0),
  f2DHistPtV0MCvsPt_Lambda(0),
  f2DHistPtRecMCvsPt_AntiLambda(0),
  f2DHistMassRecMCvsPt_Lambda(0),
  f2DHistMassV0MCvsPt_Lambda(0),
  f2DHistMassRecMCvsPt_AntiLambda(0),
  f2DHistXPULLvsPt_Lambda(0),
  f2DHistXPULLvsPt_Lambda_V0(0),
  f2DHistXPULLvsPt_AntiLambda(0),
  f2DHistPtPULLvsPt_Lambda(0),
  f2DHistPtPULLvsPt_Lambda_V0(0),
  f2DHistPtPULLvsPt_AntiLambda(0),
  f2DHistMassPULLvsPt_Lambda(0),
  f2DHistMassPULLvsPt_Lambda_V0(0),
  f2DHistMassPULLvsPt_AntiLambda(0),
  f2DHistMassPULLvsRadius_Lambda(0),
  f2DHistMassPULLvsRadius_AntiLambda(0),
  fHistEvents(0),
  fHTrigger(0),
  f2DHistArmenterosPodolanski_FirstDaugPos(0),
  f2DHistArmenterosPodolanski_FirstDaugNeg(0),
  f2DHistArmenterosPodolanski_candidate(0),
  f2DHistArmenterosPodolanski_Lam(0),
  f2DHistArmenterosPodolanski_AntiLam(0),
  f2DHistChargeDaughters(0),
  f2DHistV0XY_OnFly(0),
  f2DHistV0XY_Offline(0),
  f2DHistV0XY_FirstDaugPos(0),
  f2DHistV0XY_FirstDaugNeg(0),
  f2DHistLambdaXY(0),
  f2DHistXiMinusXY_DV(0),
  f2DHistXiMinusXY_PV(0),
  f2DHistXiPlusXY_DV(0),
  f2DHistXiPlusXY_PV(0),
  fHistOnFlyStatus(0),
  fHistOnFlyStatus_FirstDaugPos(0),
  fHistOnFlyStatus_FirstDaugNeg(0),
  fHistChargeV0(0),
  fHistNDaughterV0(0),
  fHistNProngV0(0),
  fHistChargeFirstDaughter(0),
  fHistChargeSecondDaughter(0),
  fHistXtrkP(0),
  fHistYtrkP(0),
  fHistZtrkP(0),
  fHistXtrkP_XYZv(0),
  fHistYtrkP_XYZv(0),
  fHistZtrkP_XYZv(0),
  fHistXtrkP_Rec_MC(0),
  fHistYtrkP_Rec_MC(0),
  fHistZtrkP_Rec_MC(0),
  fHistXtrkP_Rec_MC_XYZv(0),
  fHistYtrkP_Rec_MC_XYZv(0),
  fHistZtrkP_Rec_MC_XYZv(0),
  fHistXtrkN(0),
  fHistYtrkN(0),
  fHistZtrkN(0),
  fHistXtrkN_XYZv(0),
  fHistYtrkN_XYZv(0),
  fHistZtrkN_XYZv(0),
  fHistXtrkN_Rec_MC(0),
  fHistYtrkN_Rec_MC(0),
  fHistZtrkN_Rec_MC(0),
  fHistXtrkN_Rec_MC_XYZv(0),
  fHistYtrkN_Rec_MC_XYZv(0),
  fHistZtrkN_Rec_MC_XYZv(0),
  fHistLDeltaLRec_Lambda(0),
  fHistLDeltaLRecMC_Lambda(0),
  fHistLDeltaLRecMC_LambdaFromXi(0),
  fHistLDeltaLRec_AntiLambda(0),
  fHistLDeltaLRecMC_AntiLambda(0),
  fHistLDeltaLRecMC_AntiLambdaFromXi(0),
  fHistXLambdaTot(0),
  fHistYLambdaTot(0),
  fHistZLambdaTot(0),
  fHistXLambda_KF_MC(0),
  fHistYLambda_KF_MC(0),
  fHistZLambda_KF_MC(0),
  fHistXProton_KF_MC(0),
  fHistYProton_KF_MC(0),
  fHistZProton_KF_MC(0),
  fHistXPion_KF_MC(0),
  fHistYPion_KF_MC(0),
  fHistZPion_KF_MC(0),
  fHistXLambda_V0_MC(0),
  fHistXAntiLambda_Rec_MC(0),
  fHistYAntiLambda_Rec_MC(0),
  fHistZAntiLambda_Rec_MC(0),
  fHistXLambda_PULL(0),
  fHistYLambda_PULL(0),
  fHistZLambda_PULL(0),
  fHistXProton_PULL(0),
  fHistYProton_PULL(0),
  fHistZProton_PULL(0),
  fHistXPion_PULL(0),
  fHistYPion_PULL(0),
  fHistZPion_PULL(0),
  fHistXAntiLambda_PULL(0),
  fHistYAntiLambda_PULL(0),
  fHistZAntiLambda_PULL(0),
  fHistXXiTot(0),
  fHistYXiTot(0),
  fHistZXiTot(0),
  fHistXXicZeroTot(0),
  fHistYXicZeroTot(0),
  fHistZXicZeroTot(0),
  fGenHistRapidity_Lambda(0),
  fGenHistRapidity_AntiLambda(0),
  fRecHistRapidity_Lambda_offline(0),
  fRecHistRapidity_Lambda_wSTD(0),
  fRecHistRapidity_AntiLambda_wSTD(0),
  fHistPtLambda(0),
  fHistPtAntiLambda(0),
  fHistPtLambdaTot(0),
  fHistPtXiMinus(0),
  fHistPtXiPlus(0),
  fHistPtXiTot(0),
  fHistPtXicZero(0),
  fHistPtAntiXicZero(0),
  fHistPtXicZeroTot(0),
  fHistMassK0S(0),
  fHistMassLambda_woCut(0),
  fHistMassAntiLambda_woCut(0),
  fHistMassLambdaCheck(0),
  fHistMassAntiLambdaCheck(0),
  fHistMassLambda_wSTDv0Cut(0),
  fHistMassAntiLambda_wSTDv0Cut(0),
  fHistMassLambda_BeforeSecSel(0),
  fHistMassAntiLambda_BeforeSecSel(0),
  fHistMassLambda(0),
  fHistMassLambda_woArmenterosPodolanskiCut(0),
  fHistMassLambda_woMassCut(0),
  fHistMassAntiLambda(0),
  fHistMassAntiLambda_woMassCut(0),
  fHistMassLambdaTot(0),
  fHistMassLambda_Match(0),
  fHistMassAntiLambda_Match(0),
  fHistMassLambdaTot_Match(0),
  fHistMassLambda_V0(0),
  fHistMassAntiLambda_V0(0),
  fHistMassLambdaTot_V0(0),
  fHistMassLambda_KF_V0(0),
  fHistMassLambda_KF_MC(0),
  fHistMassLambda_V0_MC(0),
  fHistMassAntiLambda_KF_V0(0),
  fHistMassAntiLambda_KF_MC(0),
  fHistMassLambda_PULL_KF(0),
  fHistMassAntiLambda_PULL_KF(0),
  fHistMassLambda_M(0),
  fHistMassAntiLambda_M(0),
  fHistMassLambdaTot_M(0),
  fHistMassLambda_MV(0),
  fHistMassAntiLambda_MV(0),
  fHistMassLambdaTot_MV(0),
  fHistMassXiMinus(0),
  fHistMassXiMinus_M(0),
  fHistMassXiMinus_Match(0),
  fHistMassXiPlus(0),
  fHistMassXiPlus_M(0),
  fHistMassXiPlus_Match(0),
  fHistMassXiTot(0),
  fHistMassXicZero_woQtCut(0),
  fHistMassXicZero(0),
  fHistMassAntiXicZero(0),
  fHistMassXicZeroTot(0),
  fHistQtDiffPionXiMinus(0),
  fHistQtDiffPionXiPlus(0),
  fHistMassXiMinus2(0),
  fHistChi2ndfProton(0),
  fHistChi2ndfPion(0),
  fHistChi2ndfLambda(0),
  fHistChi2ndfLambda_Match(0),
  fHistChi2ndfAntiLambda(0),
  fHistChi2ndfAntiLambda_Match(0),
  fHistChi2ndfLambdaTot(0),
  fHistProbProton(0),
  fHistProbPion(0),
  fHistProbLambda(0),
  fHistProbLambda_chi2cut(0),
  fHistProbLambda_Match(0),
  fHistProbAntiLambda(0),
  fHistProbAntiLambda_chi2cut(0),
  fHistProbAntiLambda_Match(0),
  fHistProbLambdaTot(0),
  fHistProbXiMinus(0),
  fHistProbXiMinus_chi2cut(0),
  fHistProbXiPlus(0),
  fHistProbXiPlus_chi2cut(0),
  fHistProbXicZero(0),
  fHistProbXicZero_chi2cut(0),
  fHistProbAntiXicZero(0),
  fHistProbAntiXicZero_chi2cut(0),
  fHistChi2ndfXiMinus(0),
  fHistChi2ndfXiPlus(0),
  fHistChi2ndfXiTot(0),
  fHistChi2ndfXicZero(0),
  fHistChi2ndfAntiXicZero(0),
  fHistChi2ndfXicZeroTot(0),
  fHistDecayLLambda(0),
  fHistDecayLAntiLambda(0),
  fHistDecayLLambdaTot(0),
  fHistDecayLXiMinus(0),
  fHistDecayLXiPlus(0),
  fHistDecayLXiTot(0),
  fHistDecayLXicZero(0),
  fHistDecayLAntiXicZero(0),
  fHistDecayLXicZeroTot(0),
  fHistCosPA_Lambda(0),
  fHistCosPA_AntiLambda(0),
  fHistCosPA_LambdaTot(0),
  fHistCosPA_XiMinus(0),
  fHistCosPA_XiPlus(0),
  fHistCosPA_XiTot(0),
  fHistPVx(0),
  fHistPVy(0),
  fHistPVz(0),
  fHCentrality(0),
  fHistMCGen_XicZeroTot(0),
  fHistMCGen_XicZero(0),
  fHistMCGen_AntiXicZero(0),
  fHistMCGen_PionTot(0),
  fHistMCGen_PionPlus(0),
  fHistMCGen_PionMinus(0),
  fHistMCGen_XiTot(0),
  fHistMCGen_XiMinus(0),
  fHistMCGen_XiPlus(0),
  fHistMCGen_Lambda(0),
  fHistMCGen_AntiLambda(0),
  fHistMCGen_PiXiInvMass(0),
  fHistMCGen_PiXiMassvsPiPt(0),
  fHistMCGen_PiXiMassvsPiPt_PionPlus(0),
  fHistMCGen_PiXiMassvsPiPt_PionMinus(0),
  fHistMCXicZeroDecayType(0),
  fHistMCXiDecayType(0),
  fHistMCpdg_All(0),
  fHistMCpdg_Dau_XicZero(0),
  fHistMCpdg_Dau_XicPM(0),
  fWriteXic0Tree(kFALSE),
  fWriteXic0MCGenTree(kFALSE),
  fWeight(0),
  fHistMCGen_Xic0Pt_weight(0),
  f2DHistMCRec_Xic0Pt_weight(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
  DefineOutput(2, AliNormalizationCounter::Class());
  DefineOutput(3, TTree::Class()); // event
  DefineOutput(4, TTree::Class()); // Xic0
  DefineOutput(5, TTree::Class()); // Xic0 MCGen
  DefineOutput(6, TList::Class()); // Xic0 weight of MC pt shape

}
//_____________________________________________________________________________
AliAnalysisTaskSEXicZero2XiPifromKFP::~AliAnalysisTaskSEXicZero2XiPifromKFP()
{
    // destructor
    if (fOutputList) {
      delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
      fOutputList = 0;
    }

    if (fOutputWeight) {
      delete fOutputWeight;     // at the end of your task, it is deleted from memory by calling this function
      fOutputWeight = 0;
    }

    if (fListCuts) {
      delete fListCuts;
      fListCuts = 0;
    }

    if (fAnaCuts) {
      delete fAnaCuts;
      fAnaCuts = 0;
    }

    if (fTree_Event) {
      delete fTree_Event;
      fTree_Event = 0;
    }

    if (fVar_Event) {
      delete fVar_Event;
      fVar_Event = 0;
    }

    if (fTree_Xic0) {
      delete fTree_Xic0;
      fTree_Xic0 = 0;
    }

    if (fVar_Xic0) {
      delete fVar_Xic0;
      fVar_Xic0 = 0;
    }

    if (fTree_Xic0MCGen) {
      delete fTree_Xic0MCGen;
      fTree_Xic0MCGen = 0;
    }

    if (fVar_Xic0MCGen) {
      delete fVar_Xic0MCGen;
      fVar_Xic0MCGen = 0;
    }

    if (fCounter) {
      delete fCounter;
      fCounter = 0;
    }


}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::Init()
{
  // Initialization

  if (fDebug > 1) AliInfo("Init");

  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsKFP(*fAnaCuts));
  PostData(1, fListCuts);

  return;

}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
  fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)

  DefineAnaHist(); // define analysis histograms

    // example of a histogram
  fHistEvents = new TH1F("fHistEvents", "fHistEvents", 18, 0.5, 18.5);
  fHistEvents->GetXaxis()->SetBinLabel(1,"Analyzed events");
  fHistEvents->GetXaxis()->SetBinLabel(2,"AliAODVertex exists");
  fHistEvents->GetXaxis()->SetBinLabel(3,"TriggerOK");
  fHistEvents->GetXaxis()->SetBinLabel(4,"IsEventSelected");
  fHistEvents->GetXaxis()->SetBinLabel(5,"V0 exists");
  fHistEvents->GetXaxis()->SetBinLabel(6,"Cascade exists");

  fHistEvents->GetXaxis()->SetBinLabel(7,"MCarray exists");
  fHistEvents->GetXaxis()->SetBinLabel(8,"MCheader exists");
  fHistEvents->GetXaxis()->SetBinLabel(9,"triggerClass!=CINT1");
  fHistEvents->GetXaxis()->SetBinLabel(10,"triggerMask!=kAnyINT");
  fHistEvents->GetXaxis()->SetBinLabel(11,"triggerMask!=kAny");
  fHistEvents->GetXaxis()->SetBinLabel(12,"vtxTitle.Contains(Z)");
  fHistEvents->GetXaxis()->SetBinLabel(13,"vtxTitle.Contains(3D)");
  fHistEvents->GetXaxis()->SetBinLabel(14,"vtxTitle.Doesn'tContain(Z-3D)");
  fHistEvents->GetXaxis()->SetBinLabel(15,Form("zVtx<=%2.0fcm", fAnaCuts->GetMaxVtxZ()));
  fHistEvents->GetXaxis()->SetBinLabel(16,"!IsEventSelected");
  fHistEvents->GetXaxis()->SetBinLabel(17,"triggerMask!=kAnyINT || triggerClass!=CINT1");
  fHistEvents->GetXaxis()->SetBinLabel(18,Form("zVtxMC<=%2.0fcm",fAnaCuts->GetMaxVtxZ()));
  fHistEvents->GetYaxis()->SetTitle("counts");

  fHTrigger = new TH1F("fHTrigger", "counter", 18, -0.5, 17.5);                                      
  fHTrigger->SetStats(kTRUE);
  fHTrigger->GetXaxis()->SetBinLabel(1,"X1");
  fHTrigger->GetXaxis()->SetBinLabel(2,"kMB");
  fHTrigger->GetXaxis()->SetBinLabel(3,"kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(4,"kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(5,"kINT7");
  fHTrigger->GetXaxis()->SetBinLabel(6,"kEMC7");
  //fHTrigger->GetXaxis()->SetBinLabel(7,"Space");
  fHTrigger->GetXaxis()->SetBinLabel(8,"kMB|kSemiCentral|kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(9,"kINT7|kEMC7");
  fHTrigger->GetXaxis()->SetBinLabel(11,"kMB&kSemiCentral");
  fHTrigger->GetXaxis()->SetBinLabel(12,"kMB&kCentral");
  fHTrigger->GetXaxis()->SetBinLabel(13,"kINT7&kEMC7");

  fHistMCGen_Lambda_Pt = new TH1F("fHistMCGen_Lambda_Pt", "#it{p}_{T} distribution of #Lambda", 20, 0., 20.);
  fHistMCGen_AntiLambda_Pt = new TH1F("fHistMCGen_AntiLambda_Pt", "#it{p}_{T} distribution of #bar{#Lambda}", 20, 0., 20.);
  fHistMCGen_Lambda_Pt_wYcut = new TH1F("fHistMCGen_Lambda_Pt_wYcut", "#it{p}_{T} distribution of #Lambda", 20, 0., 20.);
  fHistMCGen_AntiLambda_Pt_wYcut = new TH1F("fHistMCGen_AntiLambda_Pt_wYcut", "#it{p}_{T} distribution of #bar{#Lambda}", 20, 0., 20.);
  fHistMCGen_XiMinus_Pt = new TH1F("fHistMCGen_XiMinus_Pt", "#it{p}_{T} distribution of #Xi^{-}", 20, 0., 20.);
  fHistMCGen_XiPlus_Pt = new TH1F("fHistMCGen_XiPlus_Pt", "#it{p}_{T} distribution of #Xi^{+}", 20, 0., 20.);
  fHistMCGen_XiMinus_Pt_wYcut = new TH1F("fHistMCGen_XiMinus_Pt_wYcut", "#it{p}_{T} distribution of #Xi^{-}", 20, 0., 20.);
  fHistMCGen_XiPlus_Pt_wYcut = new TH1F("fHistMCGen_XiPlus_Pt_wYcut", "#it{p}_{T} distribution of #Xi^{+}", 20, 0., 20.);


  fCounterGen_Cuts_Lambda = new TH1F("fCounterGen_Cuts_Lambda", "Counter of #Lambda (Gen)", 6, 0.5, 6.5);
  fCounterGen_Cuts_Lambda->GetXaxis()->SetBinLabel(1,"inclusive #Lambda");
  fCounterGen_Cuts_Lambda->GetXaxis()->SetBinLabel(2,"pri. #Lambda");
  fCounterGen_Cuts_Lambda->GetXaxis()->SetBinLabel(3,"physical pri. #Lambda");
  fCounterGen_Cuts_Lambda->GetXaxis()->SetBinLabel(4,"sec. from weak decay");
  fCounterGen_Cuts_Lambda->GetXaxis()->SetBinLabel(5,"sec. from material");
  fCounterGen_Cuts_Lambda->GetXaxis()->SetBinLabel(6,"from subsidiary");

  fCounterGen_Cuts_AntiLambda = new TH1F("fCounterGen_Cuts_AntiLambda", "Counter of #bar{#Lambda} (Gen)", 6, 0.5, 6.5);
  fCounterGen_Cuts_AntiLambda->GetXaxis()->SetBinLabel(1,"inclusive #bar{#Lambda}");
  fCounterGen_Cuts_AntiLambda->GetXaxis()->SetBinLabel(2,"pri. #bar{#Lambda}");
  fCounterGen_Cuts_AntiLambda->GetXaxis()->SetBinLabel(3,"physical pri. #bar{#Lambda}");
  fCounterGen_Cuts_AntiLambda->GetXaxis()->SetBinLabel(4,"sec. from weak decay");
  fCounterGen_Cuts_AntiLambda->GetXaxis()->SetBinLabel(5,"sec. from material");
  fCounterGen_Cuts_AntiLambda->GetXaxis()->SetBinLabel(6,"from subsidiary");

  fCounterGen_Cuts_XiMinus = new TH1F("fCounterGen_Cuts_XiMinus", "Counter of #Xi^{-} (Gen)", 6, 0.5, 6.5);
  fCounterGen_Cuts_XiMinus->GetXaxis()->SetBinLabel(1,"inclusive #Xi^{-}");
  fCounterGen_Cuts_XiMinus->GetXaxis()->SetBinLabel(2,"pri. #Xi^{-}");
  fCounterGen_Cuts_XiMinus->GetXaxis()->SetBinLabel(3,"physical pri. #Xi^{-}");
  fCounterGen_Cuts_XiMinus->GetXaxis()->SetBinLabel(4,"sec. from weak decay");
  fCounterGen_Cuts_XiMinus->GetXaxis()->SetBinLabel(5,"sec. from material");
  fCounterGen_Cuts_XiMinus->GetXaxis()->SetBinLabel(6,"from subsidiary");

  fCounterGen_Cuts_XiPlus = new TH1F("fCounterGen_Cuts_XiPlus", "Counter of #Xi^{+} (Gen)", 6, 0.5, 6.5);
  fCounterGen_Cuts_XiPlus->GetXaxis()->SetBinLabel(1,"inclusive #Xi^{+}");
  fCounterGen_Cuts_XiPlus->GetXaxis()->SetBinLabel(2,"pri. #Xi^{+}");
  fCounterGen_Cuts_XiPlus->GetXaxis()->SetBinLabel(3,"physical pri. #Xi^{+}");
  fCounterGen_Cuts_XiPlus->GetXaxis()->SetBinLabel(4,"sec. from weak decay");
  fCounterGen_Cuts_XiPlus->GetXaxis()->SetBinLabel(5,"sec. from material");
  fCounterGen_Cuts_XiPlus->GetXaxis()->SetBinLabel(6,"from subsidiary");

  fCounterRecMC_Cuts_Lambda = new TH1F("fCounterRecMC_Cuts_Lambda", "Counter of #Lambda (RecMC)", 8, 0.5, 8.5);
  fCounterRecMC_Cuts_Lambda->GetXaxis()->SetBinLabel(1,"offline v0s");
  fCounterRecMC_Cuts_Lambda->GetXaxis()->SetBinLabel(2,"N_daughter=2");
  fCounterRecMC_Cuts_Lambda->GetXaxis()->SetBinLabel(3,"daug track exist");
  fCounterRecMC_Cuts_Lambda->GetXaxis()->SetBinLabel(4,"daug cov. exist");
  fCounterRecMC_Cuts_Lambda->GetXaxis()->SetBinLabel(5,"daug cov. QA");
  fCounterRecMC_Cuts_Lambda->GetXaxis()->SetBinLabel(6,"STD v0 cuts");
  fCounterRecMC_Cuts_Lambda->GetXaxis()->SetBinLabel(7,"KFP QA");
  fCounterRecMC_Cuts_Lambda->GetXaxis()->SetBinLabel(8,"|y(#Lambda)|<0.5");
//  fCounterRecMC_Cuts_Lambda->GetXaxis()->SetBinLabel(12,"sec. #Lambda cut");

  fCounterRecMC_Cuts_AntiLambda = new TH1F("fCounterRecMC_Cuts_AntiLambda", "Counter of #bar{#Lambda} (RecMC)", 8, 0.5, 8.5);
  fCounterRecMC_Cuts_AntiLambda->GetXaxis()->SetBinLabel(1,"offline v0s");
  fCounterRecMC_Cuts_AntiLambda->GetXaxis()->SetBinLabel(2,"N_daughter=2");
  fCounterRecMC_Cuts_AntiLambda->GetXaxis()->SetBinLabel(3,"daug track exist");
  fCounterRecMC_Cuts_AntiLambda->GetXaxis()->SetBinLabel(4,"daug cov. exist");
  fCounterRecMC_Cuts_AntiLambda->GetXaxis()->SetBinLabel(5,"daug cov. QA");
  fCounterRecMC_Cuts_AntiLambda->GetXaxis()->SetBinLabel(6,"STD v0 cuts");
  fCounterRecMC_Cuts_AntiLambda->GetXaxis()->SetBinLabel(7,"KFP QA");
  fCounterRecMC_Cuts_AntiLambda->GetXaxis()->SetBinLabel(8,"|y(#bar{#Lambda})|<0.5");
//  fCounterRecMC_Cuts_AntiLambda->GetXaxis()->SetBinLabel(12,"sec. #bar{#Lambda} cut");

  fCounterRecMC_Cuts_XiMinus = new TH1F("fCounterRecMC_Cuts_XiMinus", "Counter of #Xi^{-} (RecMC)", 5, 0.5, 5.5);
  fCounterRecMC_Cuts_XiMinus->GetXaxis()->SetBinLabel(1,"#pi cov. QA");
  fCounterRecMC_Cuts_XiMinus->GetXaxis()->SetBinLabel(2,"STD sec. #pi cuts");
  fCounterRecMC_Cuts_XiMinus->GetXaxis()->SetBinLabel(3,"KFP QA");
  fCounterRecMC_Cuts_XiMinus->GetXaxis()->SetBinLabel(4,"PV constraint");
  fCounterRecMC_Cuts_XiMinus->GetXaxis()->SetBinLabel(5,"|y(#Xi^{-})|<0.5");

  fCounterRecMC_Cuts_XiPlus = new TH1F("fCounterRecMC_Cuts_XiPlus", "Counter of #Xi^{+} (RecMC)", 5, 0.5, 5.5);
  fCounterRecMC_Cuts_XiPlus->GetXaxis()->SetBinLabel(1,"#pi cov. QA");
  fCounterRecMC_Cuts_XiPlus->GetXaxis()->SetBinLabel(2,"STD sec. #pi cuts");
  fCounterRecMC_Cuts_XiPlus->GetXaxis()->SetBinLabel(3,"KFP QA");
  fCounterRecMC_Cuts_XiPlus->GetXaxis()->SetBinLabel(4,"PV constraint");
  fCounterRecMC_Cuts_XiPlus->GetXaxis()->SetBinLabel(5,"|y(#Xi^{+})|<0.5");

  fCounterRec_Cuts_Lambda = new TH1F("fCounterRec_Cuts_Lambda", "Counter of #Lambda candidates (Rec)", 8, 0.5, 8.5);
  fCounterRec_Cuts_Lambda->GetXaxis()->SetBinLabel(1,"offline v0s");
  fCounterRec_Cuts_Lambda->GetXaxis()->SetBinLabel(2,"N_daughter=2");
  fCounterRec_Cuts_Lambda->GetXaxis()->SetBinLabel(3,"daug track exist");
  fCounterRec_Cuts_Lambda->GetXaxis()->SetBinLabel(4,"daug cov. exist");
  fCounterRec_Cuts_Lambda->GetXaxis()->SetBinLabel(5,"daug cov. QA");
  fCounterRec_Cuts_Lambda->GetXaxis()->SetBinLabel(6,"STD v0 cuts");
  fCounterRec_Cuts_Lambda->GetXaxis()->SetBinLabel(7,"KFP QA");
  fCounterRec_Cuts_Lambda->GetXaxis()->SetBinLabel(8,"|y(#Lambda)|<0.5");
//  fCounterRec_Cuts_Lambda->GetXaxis()->SetBinLabel(12,"sec. #Lambda cut");

  fCounterRec_Cuts_AntiLambda = new TH1F("fCounterRec_Cuts_AntiLambda", "Counter of #bar{#Lambda} candidates (Rec)", 8, 0.5, 8.5);
  fCounterRec_Cuts_AntiLambda->GetXaxis()->SetBinLabel(1,"offline v0s");
  fCounterRec_Cuts_AntiLambda->GetXaxis()->SetBinLabel(2,"N_daughter=2");
  fCounterRec_Cuts_AntiLambda->GetXaxis()->SetBinLabel(3,"daug track exist");
  fCounterRec_Cuts_AntiLambda->GetXaxis()->SetBinLabel(4,"daug cov. exist");
  fCounterRec_Cuts_AntiLambda->GetXaxis()->SetBinLabel(5,"daug cov. QA");
  fCounterRec_Cuts_AntiLambda->GetXaxis()->SetBinLabel(6,"STD v0 cuts");
  fCounterRec_Cuts_AntiLambda->GetXaxis()->SetBinLabel(7,"KFP QA");
  fCounterRec_Cuts_AntiLambda->GetXaxis()->SetBinLabel(8,"|y(#bar{#Lambda})|<0.5");
//  fCounterRec_Cuts_AntiLambda->GetXaxis()->SetBinLabel(12,"sec. #bar{#Lambda} cut");

  fCounterRec_Cuts_XiMinus = new TH1F("fCounterRec_Cuts_XiMinus", "Counter of #Xi^{-} candidates (Rec)", 5, 0.5, 5.5);
  fCounterRec_Cuts_XiMinus->GetXaxis()->SetBinLabel(1,"#pi cov. QA");
  fCounterRec_Cuts_XiMinus->GetXaxis()->SetBinLabel(2,"STD sec. #pi");
  fCounterRec_Cuts_XiMinus->GetXaxis()->SetBinLabel(3,"KFP QA");
  fCounterRec_Cuts_XiMinus->GetXaxis()->SetBinLabel(4,"PV constraint");
  fCounterRec_Cuts_XiMinus->GetXaxis()->SetBinLabel(5,"|y(#Xi^{-})|<0.5");

  fCounterRec_Cuts_XiPlus = new TH1F("fCounterRec_Cuts_XiPlus", "Counter of #Xi^{+} candidates (Rec)", 5, 0.5, 5.5);
  fCounterRec_Cuts_XiPlus->GetXaxis()->SetBinLabel(1,"#pi cov. QA");
  fCounterRec_Cuts_XiPlus->GetXaxis()->SetBinLabel(2,"STD sec. #pi cuts");
  fCounterRec_Cuts_XiPlus->GetXaxis()->SetBinLabel(3,"KFP QA");
  fCounterRec_Cuts_XiPlus->GetXaxis()->SetBinLabel(4,"PV constraint");
  fCounterRec_Cuts_XiPlus->GetXaxis()->SetBinLabel(5,"|y(#Xi^{+})|<0.5");

  f2DCounterRecMC_CutsVsPt_Lambda = new TH2F("f2DCounterRecMC_CutsVsPt_Lambda", "", 20, 0.5, 20.5, 20, 0., 20.);
  f2DCounterRecMC_CutsVsPt_Lambda->GetXaxis()->SetBinLabel(1,"offline v0s");
  f2DCounterRecMC_CutsVsPt_Lambda->GetXaxis()->SetBinLabel(2,"N_daughter=2");
  f2DCounterRecMC_CutsVsPt_Lambda->GetXaxis()->SetBinLabel(3,"daug track exist");
  f2DCounterRecMC_CutsVsPt_Lambda->GetXaxis()->SetBinLabel(4,"daug cov. exist");
  f2DCounterRecMC_CutsVsPt_Lambda->GetXaxis()->SetBinLabel(5,"daug cov. QA");
  f2DCounterRecMC_CutsVsPt_Lambda->GetXaxis()->SetBinLabel(6,"STD v0 cuts");
  f2DCounterRecMC_CutsVsPt_Lambda->GetXaxis()->SetBinLabel(7,"KFP QA");
  f2DCounterRecMC_CutsVsPt_Lambda->GetXaxis()->SetBinLabel(8,"|y(#Lambda)|<0.5");

  f2DCounterRecMC_CutsVsPt_AntiLambda = new TH2F("f2DCounterRecMC_CutsVsPt_AntiLambda", "", 20, 0.5, 20.5, 20, 0., 20.);
  f2DCounterRecMC_CutsVsPt_AntiLambda->GetXaxis()->SetBinLabel(1,"offline v0s");
  f2DCounterRecMC_CutsVsPt_AntiLambda->GetXaxis()->SetBinLabel(2,"N_daughter=2");
  f2DCounterRecMC_CutsVsPt_AntiLambda->GetXaxis()->SetBinLabel(3,"daug track exist");
  f2DCounterRecMC_CutsVsPt_AntiLambda->GetXaxis()->SetBinLabel(4,"daug cov. exist");
  f2DCounterRecMC_CutsVsPt_AntiLambda->GetXaxis()->SetBinLabel(5,"daug cov. QA");
  f2DCounterRecMC_CutsVsPt_AntiLambda->GetXaxis()->SetBinLabel(6,"STD v0 cuts");
  f2DCounterRecMC_CutsVsPt_AntiLambda->GetXaxis()->SetBinLabel(7,"KFP QA");
  f2DCounterRecMC_CutsVsPt_AntiLambda->GetXaxis()->SetBinLabel(8,"|y(#bar{#Lambda})|<0.5");

  f2DCounterRecMC_CutsVsPt_XiMinus = new TH2F("f2DCounterRecMC_CutsVsPt_XiMinus", "", 20, 0.5, 20.5, 20, 0., 20.);
  f2DCounterRecMC_CutsVsPt_XiMinus->GetXaxis()->SetBinLabel(1,"#pi cov. QA");
  f2DCounterRecMC_CutsVsPt_XiMinus->GetXaxis()->SetBinLabel(2,"STD sec. #pi cuts");
  f2DCounterRecMC_CutsVsPt_XiMinus->GetXaxis()->SetBinLabel(3,"KFP QA");
  f2DCounterRecMC_CutsVsPt_XiMinus->GetXaxis()->SetBinLabel(4,"PV constraint");
  f2DCounterRecMC_CutsVsPt_XiMinus->GetXaxis()->SetBinLabel(5,"|y(#Xi^{-})|<0.5");

  f2DCounterRecMC_CutsVsPt_XiPlus = new TH2F("f2DCounterRecMC_CutsVsPt_XiPlus", "", 20, 0.5, 20.5, 20, 0., 20.);
  f2DCounterRecMC_CutsVsPt_XiPlus->GetXaxis()->SetBinLabel(1,"#pi cov. QA");
  f2DCounterRecMC_CutsVsPt_XiPlus->GetXaxis()->SetBinLabel(2,"STD sec. #pi cuts");
  f2DCounterRecMC_CutsVsPt_XiPlus->GetXaxis()->SetBinLabel(3,"KFP QA");
  f2DCounterRecMC_CutsVsPt_XiPlus->GetXaxis()->SetBinLabel(4,"PV constraint");
  f2DCounterRecMC_CutsVsPt_XiPlus->GetXaxis()->SetBinLabel(5,"|y(#Xi^{+})|<0.5");

  f2DHistMassPtLambda = new TH2F("f2DHistMassPtLambda", "#Lambda InvMass vs. P_{T}", 20, 0., 20., 1000, 1., 2.);
  f2DHistMassPtAntiLambda = new TH2F("f2DHistMassPtAntiLambda", "#bar{#Lambda} InvMass vs. P_{T}", 20, 0., 20., 1000, 1., 2.);
  f2DHistMassPtXiMinus = new TH2F("f2DHistMassPtXiMinus", "#Xi^{-} InvMass vs. P_{T}", 20, 0., 20., 1000, 1., 2.);
  f2DHistMassPtXiPlus = new TH2F("f2DHistMassPtXiPlus", "#Xi^{+} InvMass vs. P_{T}", 20, 0., 20., 1000, 1., 2.);
  f2DHistMassPtXicZero_woQtCut = new TH2F("f2DHistMassPtXicZero_woQtCut", "#Xi_{c}^{0} InvMass vs. P_{T}", 20, 0., 20., 1000, 2., 3.);
  f2DHistMassPtXicZero = new TH2F("f2DHistMassPtXicZero", "#Xi_{c}^{0} InvMass vs. P_{T}", 20, 0., 20., 1000, 2., 3.);
  f2DHistMassPtAntiXicZero = new TH2F("f2DHistMassPtAntiXicZero", "#bar{#Xi_{c}^{0}} InvMass vs. P_{T}", 20, 0., 20., 1000, 2., 3.);
  f2DHistMassPtXicZeroTot = new TH2F("f2DHistMassPtXicZeroTot", "#Xi_{c}^{0} + #bar{#Xi_{c}^{0}} InvMass vs. P_{T}", 20, 0., 20., 1000, 2., 3.);

  f2DHistChi2vsNDF_Lambda = new TH2F("f2DHistChi2vsNDF_Lambda", "#Lambda: #chi^{2} vs. NDF", 10, -5., 5., 2000, -100., 100.);
  f2DHistChi2vsNDF_Lambda_Match = new TH2F("f2DHistChi2vsNDF_Lambda_Match", "#Lambda: #chi^{2} vs. NDF", 10, -5., 5., 2000, -100., 100.);
  f2DHistChi2vsNDF_AntiLambda = new TH2F("f2DHistChi2vsNDF_AntiLambda", "#bar{#Lambda}: #chi^{2} vs. NDF", 10, -5., 5., 2000, -100., 100.);

  f2DHistXKFMCvsPt_Proton = new TH2F("f2DHistXKFMCvsPt_Proton", "Proton at decay vertex: (x_{Rec} - x_{MC}) vs. P_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistYKFMCvsPt_Proton = new TH2F("f2DHistYKFMCvsPt_Proton", "Proton at decay vertex: (y_{Rec} - y_{MC}) vs. P_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistZKFMCvsPt_Proton = new TH2F("f2DHistZKFMCvsPt_Proton", "Proton at decay vertex: (z_{Rec} - z_{MC}) vs. P_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistXKFMCvsPt_Pion = new TH2F("f2DHistXKFMCvsPt_Pion", "Pion at decay vertex: (x_{Rec} - x_{MC}) vs. P_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistYKFMCvsPt_Pion = new TH2F("f2DHistYKFMCvsPt_Pion", "Pion at decay vertex: (y_{Rec} - y_{MC}) vs. P_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistZKFMCvsPt_Pion = new TH2F("f2DHistZKFMCvsPt_Pion", "Pion at decay vertex: (z_{Rec} - z_{MC}) vs. P_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistXRecMCvsPt_Lambda = new TH2F("f2DHistXRecMCvsPt_Lambda", "#Lambda: (x_{Rec} - x_{MC}) vs. P_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistXRecMCvsPt_AntiLambda = new TH2F("f2DHistXRecMCvsPt_AntiLambda", "#bar{#Lambda}: (x_{Rec} - x_{MC}) vs. P_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistXV0MCvsPt_Lambda = new TH2F("f2DHistXV0MCvsPt_Lambda", "#Lambda: (x_{V0} - x_{MC}) vs. P_{T}", 20, 0., 20., 2000, -1., 1.);

  f2DHistPtRecMCvsPt_Lambda = new TH2F("f2DHistPtRecMCvsPt_Lambda", "#Lambda: ((p_{T}^{Rec} - p_{T}^{MC})/p_{T}^{Rec}) vs. p_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistPtV0MCvsPt_Lambda = new TH2F("f2DHistPtV0MCvsPt_Lambda", "#Lambda: ((p_{T}^{Rec} - p_{T}^{MC})/p_{T}^{Rec}) vs. p_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistPtRecMCvsPt_AntiLambda = new TH2F("f2DHistPtRecMCvsPt_AntiLambda", "#bar{#Lambda}: ((p_{T}^{Rec} - p_{T}^{MC})/p_{T}^{Rec}) vs. p_{T}", 20, 0., 20., 2000, -1., 1.);

  f2DHistMassRecMCvsPt_Lambda = new TH2F("f2DHistMassRecMCvsPt_Lambda", "#Lambda: (m_{Rec} - m_{MC}) vs. p_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistMassV0MCvsPt_Lambda = new TH2F("f2DHistMassV0MCvsPt_Lambda", "#Lambda: (m_{Rec} - m_{MC}) vs. p_{T}", 20, 0., 20., 2000, -1., 1.);
  f2DHistMassRecMCvsPt_AntiLambda = new TH2F("f2DHistMassRecMCvsPt_AntiLambda", "#bar{#Lambda}: (m_{Rec} - m_{MC}) vs. p_{T}", 20, 0., 20., 2000, -1., 1.);

  f2DHistXPULLvsPt_Proton = new TH2F("f2DHistXPULLvsPt_Proton", "Proton at decay vertex: x_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistYPULLvsPt_Proton = new TH2F("f2DHistYPULLvsPt_Proton", "Proton at decay vertex: y_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistZPULLvsPt_Proton = new TH2F("f2DHistZPULLvsPt_Proton", "Proton at decay vertex: z_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistXPULLvsPt_Pion = new TH2F("f2DHistXPULLvsPt_Pion", "Pion at decay vertex: x_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistYPULLvsPt_Pion = new TH2F("f2DHistYPULLvsPt_Pion", "Pion at decay vertex: y_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistZPULLvsPt_Pion = new TH2F("f2DHistZPULLvsPt_Pion", "Pion at decay vertex: z_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistXPULLvsPt_Lambda = new TH2F("f2DHistXPULLvsPt_Lambda", "#Lambda: x_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistXPULLvsPt_Lambda_V0 = new TH2F("f2DHistXPULLvsPt_Lambda_V0", "#Lambda: x_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistXPULLvsPt_AntiLambda = new TH2F("f2DHistXPULLvsPt_AntiLambda", "#bar{#Lambda}: x_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);

  f2DHistPtPULLvsPt_Lambda = new TH2F("f2DHistPtPULLvsPt_Lambda", "#Lambda: p_{T}^{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistPtPULLvsPt_Lambda_V0 = new TH2F("f2DHistPtPULLvsPt_Lambda_V0", "#Lambda: p_{T}^{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistPtPULLvsPt_AntiLambda = new TH2F("f2DHistPtPULLvsPt_AntiLambda", "#bar{#Lambda}: p_{T}^{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);

  f2DHistMassPULLvsPt_Lambda = new TH2F("f2DHistMassPULLvsPt_Lambda", "#Lambda: m_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistMassPULLvsPt_Lambda_V0 = new TH2F("f2DHistMassPULLvsPt_Lambda_V0", "#Lambda: m_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistMassPULLvsPt_AntiLambda = new TH2F("f2DHistMassPULLvsPt_AntiLambda", "#bar{#Lambda}: m_{PULL} vs. p_{T}", 20, 0., 20., 2000, -10., 10.);
  f2DHistMassPULLvsRadius_Lambda = new TH2F("f2DHistMassPULLvsRadius_Lambda", "#Lambda: m_{PULL} vs. Radius", 200, 0., 20., 2000, -10., 10.);
  f2DHistMassPULLvsRadius_AntiLambda = new TH2F("f2DHistMassPULLvsRadius_AntiLambda", "#bar{#Lambda}: m_{PULL} vs. Radius", 200, 0., 20., 2000, -10., 10.);

  f2DHistArmenterosPodolanski_FirstDaugPos = new TH2F("f2DHistArmenterosPodolanski_FirstDaugPos", "f2DHistArmenterosPodolanski_FirstDaugPos", 1000, -1., 1., 1000, 0., 1.);
  f2DHistArmenterosPodolanski_FirstDaugNeg = new TH2F("f2DHistArmenterosPodolanski_FirstDaugNeg", "f2DHistArmenterosPodolanski_FirstDaugNeg", 1000, -1., 1., 1000, 0., 1.);
  f2DHistArmenterosPodolanski_candidate = new TH2F("f2DHistArmenterosPodolanski_candidate", "f2DHistArmenterosPodolanski_candidate", 1000, -1., 1., 1000, 0., 1.);
  f2DHistArmenterosPodolanski_Lam = new TH2F("f2DHistArmenterosPodolanski_Lam", "f2DHistArmenterosPodolanski_Lam", 1000, -1., 1., 1000, 0., 1.);
  f2DHistArmenterosPodolanski_AntiLam = new TH2F("f2DHistArmenterosPodolanski_AntiLam", "f2DHistArmenterosPodolanski_AntiLam", 1000, -1., 1., 1000, 0., 1.);
  f2DHistChargeDaughters = new TH2F("f2DHistChargeDaughters", "f2DHistChargeDaughters", 11, -5.5, 5.5, 11., -5.5, 5.5);
  f2DHistV0XY_OnFly = new TH2F("f2DHistV0XY_OnFly", "f2DHistV0XY_OnFly", 1000, -10., 10., 1000, -10., 10.);
  f2DHistV0XY_Offline = new TH2F("f2DHistV0XY_Offline", "f2DHistV0XY_Offline", 1000, -10., 10., 1000, -10., 10.);
  f2DHistV0XY_FirstDaugPos = new TH2F("f2DHistV0XY_FirstDaugPos", "f2DHistV0XY_FirstDaugPos", 1000, -10., 10., 1000, -10., 10.);
  f2DHistV0XY_FirstDaugNeg = new TH2F("f2DHistV0XY_FirstDaugNeg", "f2DHistV0XY_FirstDaugNeg", 1000, -10., 10., 1000, -10., 10.);
  f2DHistLambdaXY = new TH2F("f2DHistLambdaXY", "f2DHistLambdaXY", 1000, -10., 10., 1000, -10., 10.);
  f2DHistXiMinusXY_DV = new TH2F("f2DHistXiMinusXY_DV", "f2DHistXiMinusXY_DV", 1000, -10., 10., 1000, -10., 10.);
  f2DHistXiMinusXY_PV = new TH2F("f2DHistXiMinusXY_PV", "f2DHistXiMinusXY_PV", 1000, -10., 10., 1000, -10., 10.);
  f2DHistXiPlusXY_DV = new TH2F("f2DHistXiPlusXY_DV", "f2DHistXiPlusXY_DV", 1000, -10., 10., 1000, -10., 10.);
  f2DHistXiPlusXY_PV = new TH2F("f2DHistXiPlusXY_PV", "f2DHistXiPlusXY_PV", 1000, -10., 10., 1000, -10., 10.);

  fHistChargeV0 = new TH1F("fHistChargeV0", "fHistChargeV0", 11, -5.5, 5.5);
  fHistOnFlyStatus = new TH1F("fHistOnFlyStatus"," fHistOnFlyStatus", 3, -1.5, 1.5);
  fHistOnFlyStatus_FirstDaugPos = new TH1F("fHistOnFlyStatus__FirstDaugPos"," fHistOnFlyStatus__FirstDaugPos", 3, -1.5, 1.5);
  fHistOnFlyStatus_FirstDaugNeg = new TH1F("fHistOnFlyStatus_FirstDaugNeg", "fHistOnFlyStatus_FirstDaugNeg", 3, -1.5, 1.5);
  fHistNDaughterV0 = new TH1F("fHistNDaughterV0", "fHistNDaughterV0", 11, -0.5, 10.5);
  fHistNProngV0 = new TH1F("fHistNProngV0", "fHistNProngV0", 11, -0.5, 10.5);
  fHistChargeFirstDaughter = new TH1F("fHistChargeFirstDaughter", "fHistChargeFirstDaughter", 11, -5.5, 5.5);
  fHistChargeSecondDaughter = new TH1F("fHistChargeSecondDaughter", "fHistChargeSecondDaughter", 11, -5.5, 5.5);
  fHistXtrkP       = new TH1F("fHistXtrkP", "fHistXtrkP", 400, -20., 20.);
  fHistYtrkP       = new TH1F("fHistYtrkP", "fHistYtrkP", 400, -20., 20.);
  fHistZtrkP       = new TH1F("fHistZtrkP", "fHistZtrkP", 400, -20., 20.);
  fHistXtrkP_XYZv  = new TH1F("fHistXtrkP_XYZv", "fHistXtrkP_XYZv", 400, -20., 20.);
  fHistYtrkP_XYZv  = new TH1F("fHistYtrkP_XYZv", "fHistYtrkP_XYZv", 400, -20., 20.);
  fHistZtrkP_XYZv  = new TH1F("fHistZtrkP_XYZv", "fHistZtrkP_XYZv", 400, -20., 20.);
  fHistXtrkP_Rec_MC = new TH1F("fHistXtrkP_Rec_MC", "x^{Rec}-x^{MC} (track_P)", 2000, -100., 100.);
  fHistYtrkP_Rec_MC = new TH1F("fHistYtrkP_Rec_MC", "y^{Rec}-y^{MC} (track_P)", 2000, -100., 100.);
  fHistZtrkP_Rec_MC = new TH1F("fHistZtrkP_Rec_MC", "z^{Rec}-z^{MC} (track_P)", 2000, -100., 100.);
  fHistXtrkP_Rec_MC_XYZv = new TH1F("fHistXtrkP_Rec_MC_XYZv", "x^{Rec}-x^{MC} (track_P)", 2000, -100., 100.);
  fHistYtrkP_Rec_MC_XYZv = new TH1F("fHistYtrkP_Rec_MC_XYZv", "y^{Rec}-y^{MC} (track_P)", 2000, -100., 100.);
  fHistZtrkP_Rec_MC_XYZv = new TH1F("fHistZtrkP_Rec_MC_XYZv", "z^{Rec}-z^{MC} (track_P)", 2000, -100., 100.);
  fHistXtrkN       = new TH1F("fHistXtrkN", "fHistXtrkN", 400, -20., 20.);
  fHistYtrkN       = new TH1F("fHistYtrkN", "fHistYtrkN", 400, -20., 20.);
  fHistZtrkN       = new TH1F("fHistZtrkN", "fHistZtrkN", 400, -20., 20.);
  fHistXtrkN_XYZv  = new TH1F("fHistXtrkN_XYZv", "fHistXtrkN_XYZv", 400, -20., 20.);
  fHistYtrkN_XYZv  = new TH1F("fHistYtrkN_XYZv", "fHistYtrkN_XYZv", 400, -20., 20.);
  fHistZtrkN_XYZv  = new TH1F("fHistZtrkN_XYZv", "fHistZtrkN_XYZv", 400, -20., 20.);
  fHistXtrkN_Rec_MC = new TH1F("fHistXtrkN_Rec_MC", "x^{Rec}-x^{MC} (track_N)", 2000, -100., 100.);
  fHistYtrkN_Rec_MC = new TH1F("fHistYtrkN_Rec_MC", "y^{Rec}-x^{MC} (track_N)", 2000, -100., 100.);
  fHistZtrkN_Rec_MC = new TH1F("fHistZtrkN_Rec_MC", "z^{Rec}-x^{MC} (track_N)", 2000, -100., 100.);
  fHistXtrkN_Rec_MC_XYZv = new TH1F("fHistXtrkN_Rec_MC_XYZv", "x^{Rec}-x^{MC} (track_N)", 2000, -100., 100.);
  fHistYtrkN_Rec_MC_XYZv = new TH1F("fHistYtrkN_Rec_MC_XYZv", "y^{Rec}-y^{MC} (track_N)", 2000, -100., 100.);
  fHistZtrkN_Rec_MC_XYZv = new TH1F("fHistZtrkN_Rec_MC_XYZv", "z^{Rec}-z^{MC} (track_N)", 2000, -100., 100.);
  fHistLDeltaLRec_Lambda  = new TH1F("fHistLDeltaLRec_Lambda", "primary #Lambda: l/#Deltal", 2000, -100., 100.);
  fHistLDeltaLRecMC_Lambda  = new TH1F("fHistLDeltaLRecMC_Lambda", "primary #Lambda: l/#Deltal", 2000, -100., 100.);
  fHistLDeltaLRecMC_LambdaFromXi  = new TH1F("fHistLDeltaLRecMC_LambdaFromXi", "#Lambda#leftarrow#Xi^{-}: l/#Deltal", 2000, -100., 100.);
  fHistLDeltaLRec_AntiLambda  = new TH1F("fHistLDeltaLRec_AntiLambda", "primary #bar{#Lambda}: l/#Deltal", 2000, -100., 100.);
  fHistLDeltaLRecMC_AntiLambda  = new TH1F("fHistLDeltaLRecMC_AntiLambda", "primary #bar{#Lambda}: l/#Deltal", 2000, -100., 100.);
  fHistLDeltaLRecMC_AntiLambdaFromXi  = new TH1F("fHistLDeltaLRecMC_AntiLambdaFromXi", "#bar{#Lambda}#leftarrow#Xi^{+}: l/#Deltal", 2000, -100., 100.);
  fHistXLambdaTot  = new TH1F("fHistXLambdaTot", "x_{Rec} at #Lambda decay vertex", 200, -100., 100.);
  fHistYLambdaTot  = new TH1F("fHistYLambdaTot", "y_{Rec} at #Lambda decay vertex", 200, -100., 100.);
  fHistZLambdaTot  = new TH1F("fHistZLambdaTot", "z_{Rec} at #Lambda decay vertex", 200, -100., 100.);
  fHistXLambda_KF_MC = new TH1F("fHistXLambda_KF_MC", "x^{Rec} - x^{MC} (#Lambda decay vertex)", 2000, -1., 1.);
  fHistYLambda_KF_MC = new TH1F("fHistYLambda_KF_MC", "y^{Rec} - y^{MC} (#Lambda decay vertex)", 2000, -1., 1.);
  fHistZLambda_KF_MC = new TH1F("fHistZLambda_KF_MC", "z^{Rec} - z^{MC} (#Lambda decay vertex)", 2000, -1., 1.);
  fHistXProton_KF_MC = new TH1F("fHistXProton_KF_MC", "x^{Rec} - x^{MC} (Proton at decay vertex)", 10000, -0.05, 0.05);
  fHistYProton_KF_MC = new TH1F("fHistYProton_KF_MC", "y^{Rec} - y^{MC} (Proton at decay vertex)", 10000, -0.05, 0.05);
  fHistZProton_KF_MC = new TH1F("fHistZProton_KF_MC", "z^{Rec} - z^{MC} (Proton at decay vertex)", 10000, -0.05, 0.05);
  fHistXPion_KF_MC = new TH1F("fHistXPion_KF_MC", "x^{Rec} - x^{MC} (Pion at decay vertex)", 10000, -0.05, 0.05);
  fHistYPion_KF_MC = new TH1F("fHistYPion_KF_MC", "y^{Rec} - y^{MC} (Pion at decay vertex)", 10000, -0.05, 0.05);
  fHistZPion_KF_MC = new TH1F("fHistZPion_KF_MC", "z^{Rec} - z^{MC} (Pion at decay vertex)", 10000, -0.05, 0.05);
  fHistXLambda_V0_MC = new TH1F("fHistXLambda_V0_MC", "x^{V0} - x^{MC} (#Lambda decay vertex)", 2000, -1., 1.);
  fHistXAntiLambda_Rec_MC = new TH1F("fHistXAntiLambda_Rec_MC", "x^{Rec} - x^{MC} (#bar{#Lambda} decay vertex)", 200, -1., 1.);
  fHistYAntiLambda_Rec_MC = new TH1F("fHistYAntiLambda_Rec_MC", "y^{Rec} - y^{MC} (#bar{#Lambda} decay vertex)", 200, -1., 1.);
  fHistZAntiLambda_Rec_MC = new TH1F("fHistZAntiLambda_Rec_MC", "z^{Rec} - z^{MC} (#bar{#Lambda} decay vertex)", 200, -1., 1.);
  fHistXProton_PULL = new TH1F("fHistXProton_PULL", "x^{PULL} (Proton at decay vertex)", 20000, -10., 10.);
  fHistYProton_PULL = new TH1F("fHistYProton_PULL", "y^{PULL} (Proton at decay vertex)", 20000, -10., 10.);
  fHistZProton_PULL = new TH1F("fHistZProton_PULL", "z^{PULL} (Proton at decay vertex)", 20000, -10., 10.);
  fHistXPion_PULL = new TH1F("fHistXPion_PULL", "x^{PULL} (Pion at decay vertex)", 20000, -10., 10.);
  fHistYPion_PULL = new TH1F("fHistYPion_PULL", "y^{PULL} (Pion at decay vertex)", 20000, -10., 10.);
  fHistZPion_PULL = new TH1F("fHistZPion_PULL", "z^{PULL} (Pion at decay vertex)", 20000, -10., 10.);
  fHistXLambda_PULL = new TH1F("fHistXLambda_PULL", "x^{PULL} (#Lambda decay vertex)", 2000, -10., 10.);
  fHistYLambda_PULL = new TH1F("fHistYLambda_PULL", "y^{PULL} (#Lambda decay vertex)", 2000, -10., 10.);
  fHistZLambda_PULL = new TH1F("fHistZLambda_PULL", "z^{PULL} (#Lambda decay vertex)", 2000, -10., 10.);
  fHistXAntiLambda_PULL = new TH1F("fHistXAntiLambda_PULL", "x^{PULL} (#bar{#Lambda} decay vertex)", 2000, -10., 10.);
  fHistYAntiLambda_PULL = new TH1F("fHistYAntiLambda_PULL", "y^{PULL} (#bar{#Lambda} decay vertex)", 2000, -10., 10.);
  fHistZAntiLambda_PULL = new TH1F("fHistZAntiLambda_PULL", "z^{PULL} (#bar{#Lambda} decay vertex)", 2000, -10., 10.);
  fHistXXiTot      = new TH1F("fHistXXiTot", "fHistXXiTot", 200, -100., 100.);
  fHistYXiTot      = new TH1F("fHistYXiTot", "fHistYXiTot", 200, -100., 100.);
  fHistZXiTot      = new TH1F("fHistZXiTot", "fHistZXiTot", 200, -100., 100.);
  fHistXXicZeroTot = new TH1F("fHistXXicZeroTot", "fHistXXicZeroTot", 200, -100., 100.);
  fHistYXicZeroTot = new TH1F("fHistYXicZeroTot", "fHistYXicZeroTot", 200, -100., 100.);
  fHistZXicZeroTot = new TH1F("fHistZXicZeroTot", "fHistZXicZeroTot", 200, -100., 100.);

  fGenHistRapidity_Lambda = new TH1F("fGenHistRapidity_Lambda", "", 200, -10., 10.);
  fGenHistRapidity_AntiLambda = new TH1F("fGenHistRapidity_AntiLambda", "", 200, -10., 10.);
  fRecHistRapidity_Lambda_offline = new TH1F("fRecHistRapidity_Lambda_offline", "", 200, -10., 10.);
  fRecHistRapidity_Lambda_wSTD = new TH1F("fRecHistRapidity_Lambda_wSTD", "", 200, -10., 10.);
  fRecHistRapidity_AntiLambda_wSTD = new TH1F("fRecHistRapidity_AntiLambda_wSTD", "", 200, -10., 10.);

  fHistPtLambda      = new TH1F("fHistPtLambda", "fHistPtLambda", 200, 0, 20);       // create your histogra
  fHistPtAntiLambda  = new TH1F("fHistPtAntiLambda", "fHistPtAntiLambda", 200, 0, 20);
  fHistPtLambdaTot  = new TH1F("fHistPtLambdaTot", "fHistPtLambdaTot", 200, 0, 20);
  fHistPtXiMinus     = new TH1F("fHistPtXiMinus", "fHistPtXiMinus", 200, 0, 20);
  fHistPtXiPlus      = new TH1F("fHistPtXiPlus", "fHistPtXiPlus", 200, 0, 20);
  fHistPtXiTot      = new TH1F("fHistPtXiTot", "fHistPtXiTot", 200, 0, 20);
  fHistPtXicZero     = new TH1F("fHistPtXicZero", "fHistPtXicZero", 200, 0, 20);
  fHistPtAntiXicZero = new TH1F("fHistPtAntiXicZero", "fHistPtAntiXicZero", 200, 0, 20);
  fHistPtXicZeroTot     = new TH1F("fHistPtXicZeroTot", "fHistPtXicZeroTot", 200, 0, 20);
  // your histogram in the output file, add it to the list!
  fHistMassK0S = new TH1F("fHistMassK0S", "HistMassK0S", 1000, 0., 1.);
  fHistMassLambda_woCut  = new TH1F("fHistMassLambda_woCut", "fHistMassLambda_woCut", 1000, 1., 2.);
  fHistMassAntiLambda_woCut  = new TH1F("fHistMassAntiLambda_woCut", "fHistMassAntiLambda_woCut", 1000, 1., 2.);
  fHistMassLambdaCheck  = new TH1F("fHistMassLambdaCheck", "fHistMassLambdaCheck", 20000, -10., 10.);
  fHistMassAntiLambdaCheck  = new TH1F("fHistMassAntiLambdaCheck", "fHistMassAntiLambdaCheck", 20000, -10., 10.);
  fHistMassLambda_wSTDv0Cut  = new TH1F("fHistMassLambda_wSTDv0Cut", "fHistMassLambda_wSTDv0Cut", 1000, 1., 2.);
  fHistMassAntiLambda_wSTDv0Cut  = new TH1F("fHistMassAntiLambda_wSTDv0Cut", "fHistMassAntiLambda_wSTDv0Cut", 1000, 1., 2.);
  fHistMassLambda_BeforeSecSel = new TH1F("fHistMassLambda_BeforeSecSel", "fHistMassLambda_BeforeSecSel", 1000, 1., 2.);
  fHistMassAntiLambda_BeforeSecSel = new TH1F("fHistMassAntiLambda_BeforeSecSel", "fHistMassAntiLambda_BeforeSecSel", 1000, 1., 2.);
  fHistMassLambda      = new TH1F("fHistMassLambda", "fHistMassLambda", 1000, 1., 2.);
  fHistMassLambda_woArmenterosPodolanskiCut = new TH1F("fHistMassLambda_woArmenterosPodolanskiCut", "fHistMassLambda_woArmenterosPodolanskiCut", 1000, 1., 2.);
  fHistMassLambda_woMassCut = new TH1F("fHistMassLambda_woMassCut", "fHistMassLambda_woMassCut", 1000, 1., 2.);
  fHistMassAntiLambda  = new TH1F("fHistMassAntiLambda", "fHistMassAntiLambda", 1000, 1., 2.);
  fHistMassAntiLambda_woMassCut = new TH1F("fHistMassAntiLambda_woMassCut", "fHistMassAntiLambda_woMassCut", 1000, 1., 2.);
  fHistMassLambdaTot   = new TH1F("fHistMassLambdaTot", "fHistMassLambdaTot", 1000, 1., 2.);
  fHistMassLambda_Match   = new TH1F("fHistMassLambda_Match", "fHistMassLambda_Match", 1000, 1., 2.);
  fHistMassAntiLambda_Match  = new TH1F("fHistMassAntiLambda_Match", "fHistMassAntiLambda", 1000, 1., 2.);
  fHistMassLambdaTot_Match   = new TH1F("fHistMassLambdaTot_Match", "fHistMassLambdaTot", 1000, 1., 2.);
  fHistMassLambda_V0      = new TH1F("fHistMassLambda_V0", "m (#Lambda) from V0", 1000, 1., 2.);
  fHistMassAntiLambda_V0  = new TH1F("fHistMassAntiLambda_V0", "m (#bar{#Lambda}) from V0", 1000, 1., 2.);
  fHistMassLambdaTot_V0   = new TH1F("fHistMassLambdaTot_V0", "fHistMassLambdaTot_V0", 1000, 1., 2.);
  fHistMassLambda_KF_V0      = new TH1F("fHistMassLambda_KF_V0", "m_{KF} - m_{V0} of #Lambda", 2000, -0.1, 0.1);
  fHistMassLambda_KF_MC      = new TH1F("fHistMassLambda_KF_MC", "m_{KF} - m_{MC} of #Lambda", 200, -0.1, 0.1);
  fHistMassLambda_V0_MC      = new TH1F("fHistMassLambda_V0_MC", "m_{V0} - m_{MC} of #Lambda", 200, -0.1, 0.1);
  fHistMassAntiLambda_KF_V0  = new TH1F("fHistMassAntiLambda_KF_V0", "m_{KF} - m_{V0} of #bar{#Lambda}", 200, -0.1, 0.1);
  fHistMassAntiLambda_KF_MC  = new TH1F("fHistMassAntiLambda_KF_MC", "m_{KF} - m_{MC} of #bar{#Lambda}", 200, -0.1, 0.1);
    fHistMassLambda_PULL_KF = new TH1F("fHistMassLambda_PULL_KF", "m_{PULL} (KF) of #Lambda", 2000, -10., 10.);
    fHistMassAntiLambda_PULL_KF = new TH1F("fHistMassAntiLambda_PULL_KF", "m_{PULL} (KF) of #bar{#Lambda}", 2000, -10., 10.);

  fHistMassLambda_M     = new TH1F("fHistMassLambda_M", "fHistMassLambda_M", 2000, 0., 2.);
  fHistMassAntiLambda_M = new TH1F("fHistMassAntiLambda_M", "fHistMassAntiLambda_M", 100, 1., 2.);
  fHistMassLambdaTot_M  = new TH1F("fHistMassLambdaTot_M", "fHistMassLambdaTot_MV", 100, 1., 2.);
  fHistMassLambda_MV     = new TH1F("fHistMassLambda_MV", "fHistMassLambda_MV", 100, 1., 2.);
  fHistMassAntiLambda_MV = new TH1F("fHistMassAntiLambda_MV", "fHistMassAntiLambda_MV", 100, 1., 2.);
  fHistMassLambdaTot_MV  = new TH1F("fHistMassLambdaTot_MV", "fHistMassLambdaTot_MV", 100, 1., 2.);
  fHistMassXiMinus     = new TH1F("fHistMassXiMinus", "fHistMassXiMinus", 1000, 1., 2.);
  fHistMassXiMinus_M     = new TH1F("fHistMassXiMinus_M", "fHistMassXiMinus_M", 1000, 1., 2.);
  fHistMassXiMinus_Match = new TH1F("fHistMassXiMinus_Match", "fHistMassXiMinus_Match", 1000, 1., 2.);
  fHistMassXiPlus      = new TH1F("fHistMassXiPlus", "fHistMassXiPlus", 1000, 1., 2.);
  fHistMassXiPlus_M     = new TH1F("fHistMassXiPlus_M", "fHistMassXiPlus_M", 2000, 0., 2.);
  fHistMassXiTot      = new TH1F("fHistMassXiTot", "fHistMassXiTot", 1000, 0., 10.);
  fHistMassXicZero_woQtCut  = new TH1F("fHistMassXicZero_woQtCut", "fHistMassXicZero_woQtCut", 10000, 0., 10.);
  fHistMassXicZero     = new TH1F("fHistMassXicZero", "fHistMassXicZero", 10000, 0., 10.);
  fHistMassAntiXicZero = new TH1F("fHistMassAntiXicZero", "fHistMassAntiXicZero", 10000, 0., 10.);
  fHistMassXicZeroTot   = new TH1F("fHistMassXicZeroTot", "fHistMassXicZeroTot", 10000, 0., 10.);
  fHistQtDiffPionXiMinus = new TH1F("fHistQtDiffPionXiMinus", "fHistQtDiffPionXiMinus", 20000, -0.01, 0.01);
  fHistQtDiffPionXiPlus = new TH1F("fHistQtDiffPionXiPlus", "fHistQtDiffPionXiPlus", 20000, -0.01, 0.01);
  fHistMassXiMinus2     = new TH1F("fHistMassXiMinus2", "fHistMassXiMinus2", 100, 1.2, 1.4);
  fHistChi2ndfProton   = new TH1F("fHistChi2ndfProton", "Proton at decay point: #chi^{2}/NDF", 2000, -100., 100.);
  fHistChi2ndfPion   = new TH1F("fHistChi2ndfPion", "Pion at decay point: #chi^{2}/NDF", 2000, -100., 100.);
  fHistChi2ndfLambda   = new TH1F("fHistChi2ndfLambda", "fHistChi2ndfLambda", 2000, -100., 100.);
  fHistChi2ndfLambda_Match   = new TH1F("fHistChi2ndfLambda_Match", "fHistChi2ndfLambda", 2000, -100., 100.);
  fHistChi2ndfAntiLambda   = new TH1F("fHistChi2ndfAntiLambda", "fHistChi2ndfAntiLambda", 2000, -100., 100.);
  fHistChi2ndfAntiLambda_Match   = new TH1F("fHistChi2ndfAntiLambda_Match", "fHistChi2ndfAntiLambda", 2000, -100., 100.);
  fHistChi2ndfLambdaTot   = new TH1F("fHistChi2ndfLambdaTot", "fHistChi2ndfLambdaTot", 2000, -100., 100.);
    fHistProbProton  = new TH1F("fHistProbProton", "Proton at decay point: prob.", 1000, 0., 1.);
    fHistProbPion  = new TH1F("fHistProbPion", "Pion at decay point: prob.", 1000, 0., 1.);
    fHistProbLambda  = new TH1F("fHistProbLambda", "Prob (#Lambda) KF", 1000, 0., 1.);
    fHistProbLambda_chi2cut  = new TH1F("fHistProbLambda_chi2cut", "Prob (#Lambda) KF", 1000, 0., 1.);
    fHistProbLambda_Match  = new TH1F("fHistProbLambda_Match", "Prob (#Lambda) KF", 1000, 0., 1.);
    fHistProbAntiLambda = new TH1F("fHistProbAntiLambda", "Prob (#bar{#Lambda}) KF", 1000, 0., 1.);
    fHistProbAntiLambda_chi2cut = new TH1F("fHistProbAntiLambda_chi2cut", "Prob (#bar{#Lambda}) KF", 1000, 0., 1.);
    fHistProbAntiLambda_Match = new TH1F("fHistProbAntiLambda_Match", "Prob (#bar{#Lambda}) KF", 1000, 0., 1.);
    fHistProbLambdaTot = new TH1F("fHistProbLambdaTot", "Prob (#Lambda+#bar{#Lambda}) KF", 1000, 0., 1.);
  fHistProbXiMinus = new TH1F("fHistProbXiMinus", "Prob (#Xi^{-}) KF", 1000, 0., 1.);
  fHistProbXiMinus_chi2cut = new TH1F("fHistProbXiMinus_chi2cut", "Prob (#Xi^{-}) KF", 1000, 0., 1.);
  fHistProbXiPlus  = new TH1F("fHistProbXiPlus", "Prob (#Xi^{+}) KF", 1000, 0., 1.);
  fHistProbXiPlus_chi2cut  = new TH1F("fHistProbXiPlus_chi2cut", "Prob (#Xi^{+}) KF", 1000, 0., 1.);
  fHistProbXicZero = new TH1F("fHistProbXicZero", "Prob (#Xi_{c}^{0}) KF", 1000, 0., 1.);
  fHistProbXicZero_chi2cut = new TH1F("fHistProbXicZero_chi2cut", "Prob (#Xi_{c}^{0}) KF", 1000, 0., 1.);
  fHistProbAntiXicZero = new TH1F("fHistProbAntiXicZero", "Prob (#bar{#Xi_{c}^{0}}) KF", 1000, 0., 1.);
  fHistProbAntiXicZero_chi2cut = new TH1F("fHistProbAntiXicZero_chi2cut", "Prob (#bar{#Xi_{c}^{0}}) KF", 1000, 0., 1.);
  fHistChi2ndfXiMinus   = new TH1F("fHistChi2ndfXiMinus", "fHistChi2ndfXiMinus", 2000, -100., 100.);
  fHistChi2ndfXiPlus   = new TH1F("fHistChi2ndfXiPlus", "fHistChi2ndfXiPlus", 2000, -100., 100.);
  fHistChi2ndfXiTot   = new TH1F("fHistChi2ndfXiTot", "fHistChi2ndfXiTot", 2000, -100., 100.);
  fHistChi2ndfXicZero   = new TH1F("fHistChi2ndfXicZero", "fHistChi2ndfXicZero", 2000, -100., 100.);
  fHistChi2ndfAntiXicZero   = new TH1F("fHistChi2ndfAntiXicZero", "fHistChi2ndfAntiXicZero", 2000, -100., 100.);
  fHistChi2ndfXicZeroTot   = new TH1F("fHistChi2ndfXicZeroTot", "fHistChi2ndfXicZeroTot", 2000, -100., 100.);
  fHistDecayLLambda = new TH1F("fHistDecayLLambda", "fHistDecayLLambda", 1000, -100., 100.);
  fHistDecayLAntiLambda = new TH1F("fHistDecayLAntiLambda", "fHistDecayLAntiLambda", 1000, -100., 100.);
  fHistDecayLLambdaTot = new TH1F("fHistDecayLLambdaTot", "fHistDecayLLambdaTot", 1000, -100., 100.);
  fHistDecayLXiMinus = new TH1F("fHistDecayLXiMinus", "fHistDecayLXiMinus", 1000, -100., 100.);
  fHistDecayLXiPlus = new TH1F("fHistDecayLXiPlus", "fHistDecayLXiPlus", 1000, -100., 100.);
  fHistDecayLXiTot = new TH1F("fHistDecayLXiTot", "fHistDecayLXiTot", 1000, -100., 100.);
  fHistDecayLXicZero = new TH1F("fHistDecayLXicZero", "fHistDecayLXicZero", 1000, -100., 100.);
  fHistDecayLAntiXicZero = new TH1F("fHistDecayLAntiXicZero", "fHistDecayLAntiXicZero", 1000, -100., 100.);
  fHistDecayLXicZeroTot = new TH1F("fHistDecayLXicZeroTot", "fHistDecayLXicZeroTot", 1000, -100., 100.);
  fHistCosPA_Lambda = new TH1F("fHistCosPA_Lambda", "CosPA_Lambda", 200, -1., 1.);
  fHistCosPA_AntiLambda = new TH1F("fHistCosPA_AntiLambda", "CosPA_AntiLambda", 200, -1., 1.);
  fHistCosPA_LambdaTot = new TH1F("fHistCosPA_LambdaTot", "CosPA_LambdaTot", 200, -1., 1.);
  fHistCosPA_XiMinus = new TH1F("fHistCosPA_XiMinus", "CosPA_XiMinus", 200, -1., 1.);
  fHistCosPA_XiPlus = new TH1F("fHistCosPA_XiPlus", "CosPA_XiPlus", 200, -1., 1.);
  fHistCosPA_XiTot = new TH1F("fHistCosPA_XiTot", "CosPA_XiTot", 200, -1., 1.);
  fHistPVx              = new TH1F("fHistPVx", "fHistPVx", 2000, -1., 1.);
  fHistPVy              = new TH1F("fHistPVy", "fHistPVy", 2000, -1., 1.);
  fHistPVz              = new TH1F("fHistPVz", "fHistPVz", 400, -20., 20.);
  fHCentrality          = new TH1F("fHCentrality", "counter", 100, 0., 100.);


  fOutputList->Add(fHistEvents); // don't forget to add it to the list! the list will be written to file, so if you want
  fOutputList->Add(fHTrigger);
  fOutputList->Add(fCounterGen_Cuts_Lambda);
  fOutputList->Add(fCounterGen_Cuts_AntiLambda);
  fOutputList->Add(fCounterRecMC_Cuts_Lambda);
  fOutputList->Add(fCounterRecMC_Cuts_AntiLambda);
  fOutputList->Add(fCounterRec_Cuts_Lambda);
  fOutputList->Add(fCounterRec_Cuts_AntiLambda);
  fOutputList->Add(fHistLDeltaLRec_Lambda);                                                         
  fOutputList->Add(fHistLDeltaLRecMC_Lambda);
  fOutputList->Add(fHistLDeltaLRecMC_LambdaFromXi);
  fOutputList->Add(fHistLDeltaLRec_AntiLambda);
  fOutputList->Add(fHistLDeltaLRecMC_AntiLambda);
  fOutputList->Add(fHistLDeltaLRecMC_AntiLambdaFromXi);
  fOutputList->Add(f2DHistArmenterosPodolanski_Lam);
  fOutputList->Add(f2DHistArmenterosPodolanski_AntiLam);
  fOutputList->Add(fHistMassLambda_M);
  fOutputList->Add(fHistMassXiPlus_M);
  fOutputList->Add(fHistMassK0S);
  /*
  fOutputList->Add(fGenHistRapidity_Lambda);
  fOutputList->Add(fRecHistRapidity_Lambda_offline);
  fOutputList->Add(fRecHistRapidity_Lambda_wSTD);
  fOutputList->Add(fGenHistRapidity_AntiLambda);
  fOutputList->Add(fRecHistRapidity_AntiLambda_wSTD);
  fOutputList->Add(fHistMCGen_Lambda_Pt);
  fOutputList->Add(fHistMCGen_AntiLambda_Pt);
  fOutputList->Add(fHistMCGen_Lambda_Pt_wYcut);
  fOutputList->Add(fHistMCGen_AntiLambda_Pt_wYcut);
  fOutputList->Add(fHistMCGen_XiMinus_Pt);
  fOutputList->Add(fHistMCGen_XiPlus_Pt);
  fOutputList->Add(fHistMCGen_XiMinus_Pt_wYcut);
  fOutputList->Add(fHistMCGen_XiPlus_Pt_wYcut);
  fOutputList->Add(fCounterGen_Cuts_XiMinus);
  fOutputList->Add(fCounterGen_Cuts_XiPlus);
  fOutputList->Add(fCounterRecMC_Cuts_XiMinus);
  fOutputList->Add(fCounterRecMC_Cuts_XiPlus);
  fOutputList->Add(fCounterRec_Cuts_XiMinus);
  fOutputList->Add(fCounterRec_Cuts_XiPlus);
  fOutputList->Add(f2DCounterRecMC_CutsVsPt_Lambda);
  fOutputList->Add(f2DCounterRecMC_CutsVsPt_AntiLambda);
  fOutputList->Add(f2DCounterRecMC_CutsVsPt_XiMinus);
  fOutputList->Add(f2DCounterRecMC_CutsVsPt_XiPlus);
  fOutputList->Add(f2DHistMassPtLambda);
  fOutputList->Add(f2DHistMassPtAntiLambda);
  fOutputList->Add(f2DHistMassPtXiMinus);
  fOutputList->Add(f2DHistMassPtXiPlus);
//  fOutputList->Add(f2DHistMassPtXicZero_woQtCut);
  fOutputList->Add(f2DHistMassPtXicZero);
  fOutputList->Add(f2DHistMassPtAntiXicZero);
  fOutputList->Add(f2DHistMassPtXicZeroTot);
  */
  /*
  fOutputList->Add(f2DHistChi2vsNDF_Lambda);
  fOutputList->Add(f2DHistChi2vsNDF_Lambda_Match);
  fOutputList->Add(f2DHistChi2vsNDF_AntiLambda);
  fOutputList->Add(f2DHistXKFMCvsPt_Proton);
  fOutputList->Add(f2DHistYKFMCvsPt_Proton);
  fOutputList->Add(f2DHistZKFMCvsPt_Proton);
  fOutputList->Add(f2DHistXPULLvsPt_Proton);
  fOutputList->Add(f2DHistYPULLvsPt_Proton);
  fOutputList->Add(f2DHistZPULLvsPt_Proton);
  fOutputList->Add(f2DHistXKFMCvsPt_Pion);
  fOutputList->Add(f2DHistYKFMCvsPt_Pion);
  fOutputList->Add(f2DHistZKFMCvsPt_Pion);
  fOutputList->Add(f2DHistXPULLvsPt_Pion);
  fOutputList->Add(f2DHistYPULLvsPt_Pion);
  fOutputList->Add(f2DHistZPULLvsPt_Pion);
  fOutputList->Add(f2DHistXRecMCvsPt_Lambda);
  fOutputList->Add(f2DHistXRecMCvsPt_AntiLambda);
  fOutputList->Add(f2DHistXV0MCvsPt_Lambda);
  fOutputList->Add(f2DHistPtRecMCvsPt_Lambda);
  fOutputList->Add(f2DHistPtV0MCvsPt_Lambda);
  fOutputList->Add(f2DHistPtRecMCvsPt_AntiLambda);
  fOutputList->Add(f2DHistMassRecMCvsPt_Lambda);
  fOutputList->Add(f2DHistMassV0MCvsPt_Lambda);
  fOutputList->Add(f2DHistMassRecMCvsPt_AntiLambda);
  fOutputList->Add(f2DHistXPULLvsPt_Lambda);
  fOutputList->Add(f2DHistXPULLvsPt_Lambda_V0);
  fOutputList->Add(f2DHistXPULLvsPt_AntiLambda);
  fOutputList->Add(f2DHistPtPULLvsPt_Lambda);
  fOutputList->Add(f2DHistPtPULLvsPt_Lambda_V0);
  fOutputList->Add(f2DHistPtPULLvsPt_AntiLambda);
  fOutputList->Add(f2DHistMassPULLvsPt_Lambda);
  fOutputList->Add(f2DHistMassPULLvsPt_Lambda_V0);
  fOutputList->Add(f2DHistMassPULLvsPt_AntiLambda);
  fOutputList->Add(f2DHistMassPULLvsRadius_Lambda);
  fOutputList->Add(f2DHistMassPULLvsRadius_AntiLambda);
  */
//  fOutputList->Add(f2DHistArmenterosPodolanski_FirstDaugPos);
//  fOutputList->Add(f2DHistArmenterosPodolanski_FirstDaugNeg);
//  fOutputList->Add(f2DHistArmenterosPodolanski_candidate);
  /*
  fOutputList->Add(fHistOnFlyStatus);
  fOutputList->Add(fHistOnFlyStatus_FirstDaugPos);
  fOutputList->Add(fHistOnFlyStatus_FirstDaugNeg);
  fOutputList->Add(fHistChargeV0);
  fOutputList->Add(fHistNDaughterV0);
  fOutputList->Add(fHistNProngV0);
  fOutputList->Add(fHistChargeFirstDaughter);
  fOutputList->Add(fHistChargeSecondDaughter);
  fOutputList->Add(f2DHistChargeDaughters);
  fOutputList->Add(f2DHistV0XY_OnFly);
  fOutputList->Add(f2DHistV0XY_Offline);
  */
  /*
  fOutputList->Add(f2DHistV0XY_FirstDaugPos);
  fOutputList->Add(f2DHistLambdaXY);
  fOutputList->Add(f2DHistXiMinusXY_DV);
  fOutputList->Add(f2DHistXiMinusXY_PV);
  fOutputList->Add(f2DHistXiPlusXY_DV);
  fOutputList->Add(f2DHistXiPlusXY_PV);
  */
  /*
  fOutputList->Add(f2DHistV0XY_FirstDaugNeg);
  fOutputList->Add(fHistXtrkP);
  fOutputList->Add(fHistYtrkP);
  fOutputList->Add(fHistZtrkP);
  fOutputList->Add(fHistXtrkP_XYZv);
  fOutputList->Add(fHistYtrkP_XYZv);
  fOutputList->Add(fHistZtrkP_XYZv);
  fOutputList->Add(fHistXtrkP_Rec_MC);
  fOutputList->Add(fHistYtrkP_Rec_MC);
  fOutputList->Add(fHistZtrkP_Rec_MC);
  fOutputList->Add(fHistXtrkP_Rec_MC_XYZv);
  fOutputList->Add(fHistYtrkP_Rec_MC_XYZv);
  fOutputList->Add(fHistZtrkP_Rec_MC_XYZv);
  fOutputList->Add(fHistXtrkN);
  fOutputList->Add(fHistYtrkN);
  fOutputList->Add(fHistZtrkN);
  fOutputList->Add(fHistXtrkN_XYZv);
  fOutputList->Add(fHistYtrkN_XYZv);
  fOutputList->Add(fHistZtrkN_XYZv);
  fOutputList->Add(fHistXtrkN_Rec_MC);
  fOutputList->Add(fHistYtrkN_Rec_MC);
  fOutputList->Add(fHistZtrkN_Rec_MC);
  fOutputList->Add(fHistXtrkN_Rec_MC_XYZv);
  fOutputList->Add(fHistYtrkN_Rec_MC_XYZv);
  fOutputList->Add(fHistZtrkN_Rec_MC_XYZv);
  fOutputList->Add(fHistXLambdaTot);
  fOutputList->Add(fHistYLambdaTot);
  fOutputList->Add(fHistZLambdaTot);
  fOutputList->Add(fHistXLambda_KF_MC);
  fOutputList->Add(fHistYLambda_KF_MC);
  fOutputList->Add(fHistZLambda_KF_MC);
  fOutputList->Add(fHistXProton_KF_MC);
  fOutputList->Add(fHistYProton_KF_MC);
  fOutputList->Add(fHistZProton_KF_MC);
  fOutputList->Add(fHistXPion_KF_MC);
  fOutputList->Add(fHistYPion_KF_MC);
  fOutputList->Add(fHistZPion_KF_MC);
  fOutputList->Add(fHistXLambda_V0_MC);
  fOutputList->Add(fHistXAntiLambda_Rec_MC);
  fOutputList->Add(fHistYAntiLambda_Rec_MC);
  fOutputList->Add(fHistZAntiLambda_Rec_MC);
  fOutputList->Add(fHistXLambda_PULL);
  fOutputList->Add(fHistYLambda_PULL);
  fOutputList->Add(fHistZLambda_PULL);
  fOutputList->Add(fHistXProton_PULL);
  fOutputList->Add(fHistYProton_PULL);
  fOutputList->Add(fHistZProton_PULL);
  fOutputList->Add(fHistXPion_PULL);
  fOutputList->Add(fHistYPion_PULL);
  fOutputList->Add(fHistZPion_PULL);
  fOutputList->Add(fHistXAntiLambda_PULL);
  fOutputList->Add(fHistYAntiLambda_PULL);
  fOutputList->Add(fHistZAntiLambda_PULL);
  fOutputList->Add(fHistXXiTot);
  fOutputList->Add(fHistYXiTot);
  fOutputList->Add(fHistZXiTot);
  fOutputList->Add(fHistXXicZeroTot);
  fOutputList->Add(fHistYXicZeroTot);
  fOutputList->Add(fHistZXicZeroTot);
  fOutputList->Add(fHistPtLambda);
  fOutputList->Add(fHistPtAntiLambda);
  fOutputList->Add(fHistPtLambdaTot);
  fOutputList->Add(fHistPtXiMinus);
  fOutputList->Add(fHistPtXiPlus);
  fOutputList->Add(fHistPtXiTot);
  fOutputList->Add(fHistPtXicZero);
  fOutputList->Add(fHistPtAntiXicZero);
  fOutputList->Add(fHistPtXicZeroTot);
  */
  
  /*
  fOutputList->Add(fHistMassLambda_woCut);
  fOutputList->Add(fHistMassAntiLambda_woCut);
//  fOutputList->Add(fHistMassLambdaCheck);
//  fOutputList->Add(fHistMassAntiLambdaCheck);
  fOutputList->Add(fHistMassLambda_wSTDv0Cut);
  fOutputList->Add(fHistMassAntiLambda_wSTDv0Cut);
  fOutputList->Add(fHistMassLambda_BeforeSecSel);
  fOutputList->Add(fHistMassAntiLambda_BeforeSecSel);
  fOutputList->Add(fHistMassLambda);
  fOutputList->Add(fHistMassLambda_woArmenterosPodolanskiCut);
  fOutputList->Add(fHistMassLambda_woMassCut);
  fOutputList->Add(fHistMassAntiLambda);
  fOutputList->Add(fHistMassAntiLambda_woMassCut);
  fOutputList->Add(fHistMassLambdaTot);
//  fOutputList->Add(fHistMassLambda_Match);
//  fOutputList->Add(fHistMassAntiLambda_Match);
//  fOutputList->Add(fHistMassLambdaTot_Match);
//  fOutputList->Add(fHistMassLambda_V0);
//  fOutputList->Add(fHistMassAntiLambda_V0);
//  fOutputList->Add(fHistMassLambdaTot_V0);
//  fOutputList->Add(fHistMassLambda_KF_V0);
//  fOutputList->Add(fHistMassLambda_KF_MC);
//  fOutputList->Add(fHistMassLambda_V0_MC);
//  fOutputList->Add(fHistMassAntiLambda_KF_V0);
//  fOutputList->Add(fHistMassAntiLambda_KF_MC);
//  fOutputList->Add(fHistMassLambda_PULL_KF);
//  fOutputList->Add(fHistMassAntiLambda_PULL_KF);
  fOutputList->Add(fHistMassLambda_M);
  fOutputList->Add(fHistMassAntiLambda_M);
//  fOutputList->Add(fHistMassLambdaTot_M);
//  fOutputList->Add(fHistMassLambda_MV);
//  fOutputList->Add(fHistMassAntiLambda_MV);
//  fOutputList->Add(fHistMassLambdaTot_MV);
  fOutputList->Add(fHistMassXiMinus);
  fOutputList->Add(fHistMassXiMinus_M);
//  fOutputList->Add(fHistMassXiMinus_Match);
  fOutputList->Add(fHistMassXiPlus);
  fOutputList->Add(fHistMassXiTot);
  fOutputList->Add(fHistMassXicZero_woQtCut);
  fOutputList->Add(fHistMassXicZero);
  fOutputList->Add(fHistMassAntiXicZero);
  fOutputList->Add(fHistMassXicZeroTot);

  fOutputList->Add(fHistQtDiffPionXiMinus);
  fOutputList->Add(fHistQtDiffPionXiPlus);
  
//  fOutputList->Add(fHistMassXiMinus2);
//  fOutputList->Add(fHistChi2ndfProton);
//  fOutputList->Add(fHistChi2ndfPion);
  fOutputList->Add(fHistChi2ndfLambda);
//  fOutputList->Add(fHistChi2ndfLambda_Match);
//  fOutputList->Add(fHistChi2ndfAntiLambda_Match);
//  fOutputList->Add(fHistChi2ndfLambdaTot);
//  fOutputList->Add(fHistProbProton);
//  fOutputList->Add(fHistProbPion);
  fOutputList->Add(fHistProbLambda);
  fOutputList->Add(fHistProbLambda_chi2cut);
//  fOutputList->Add(fHistProbLambda_Match);
  fOutputList->Add(fHistChi2ndfAntiLambda);
  fOutputList->Add(fHistProbAntiLambda);
  fOutputList->Add(fHistProbAntiLambda_chi2cut);
//  fOutputList->Add(fHistProbAntiLambda_Match);
  fOutputList->Add(fHistProbLambdaTot);
  fOutputList->Add(fHistChi2ndfXiMinus);
  fOutputList->Add(fHistProbXiMinus);
  fOutputList->Add(fHistProbXiMinus_chi2cut);
  fOutputList->Add(fHistChi2ndfXiPlus);
  fOutputList->Add(fHistProbXiPlus);
  fOutputList->Add(fHistProbXiPlus_chi2cut);
  fOutputList->Add(fHistChi2ndfXicZero);
  fOutputList->Add(fHistProbXicZero);
  fOutputList->Add(fHistProbXicZero_chi2cut);
  fOutputList->Add(fHistChi2ndfAntiXicZero);
  fOutputList->Add(fHistProbAntiXicZero);
  fOutputList->Add(fHistProbAntiXicZero_chi2cut);
//  fOutputList->Add(fHistChi2ndfXiPlus);
//  fOutputList->Add(fHistChi2ndfXiTot);
//  fOutputList->Add(fHistChi2ndfXicZero);
//  fOutputList->Add(fHistChi2ndfAntiXicZero);
//  fOutputList->Add(fHistChi2ndfXicZeroTot);
  fOutputList->Add(fHistDecayLLambda);
  fOutputList->Add(fHistDecayLAntiLambda);
  fOutputList->Add(fHistDecayLLambdaTot);
  fOutputList->Add(fHistDecayLXiMinus);
  fOutputList->Add(fHistDecayLXiPlus);
  fOutputList->Add(fHistDecayLXiTot);
  fOutputList->Add(fHistDecayLXicZero);
  fOutputList->Add(fHistDecayLAntiXicZero);
  fOutputList->Add(fHistDecayLXicZeroTot);
//  fOutputList->Add(fHistCosPA_Lambda);
//  fOutputList->Add(fHistCosPA_AntiLambda);
//  fOutputList->Add(fHistCosPA_LambdaTot);
//  fOutputList->Add(fHistCosPA_XiMinus);
//  fOutputList->Add(fHistCosPA_XiPlus);
//  fOutputList->Add(fHistCosPA_XiTot);
  fOutputList->Add(fHistPVx);
  fOutputList->Add(fHistPVy);
  fOutputList->Add(fHistPVz);
//  fOutputList->Add(fHCentrality);
  */
    
/*
  fListCuts = new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("ListCuts");
  fListCuts->Add(new AliRDHFCutsKFP(*fAnaCuts));
  PostData(3, fListCuts);
*/

  // Counter for Normalization
  TString normName="NormalizationCounter";
  AliAnalysisDataContainer *cont = GetOutputSlot(2)->GetContainer();
  if(cont) normName = (TString)cont->GetName();
  fCounter = new AliNormalizationCounter(normName.Data());
  fCounter->Init();
  PostData(2, fCounter);
  DefineEvent();
  PostData(3, fTree_Event);  // postdata will notify the analysis manager of changes / updates to the 

  DefineTreeRecXic0();
  PostData(4, fTree_Xic0);

  DefineTreeGenXic0();
  PostData(5, fTree_Xic0MCGen);

  fOutputWeight = new TList();
  fOutputWeight->SetOwner(kTRUE);
  fHistMCGen_Xic0Pt_weight = new TH1D("fHistMCGen_Xic0Pt_weight", "", 11, 1., 12.);
  f2DHistMCRec_Xic0Pt_weight = new TH2D("f2DHistMCRec_Xic0Pt_weight", "", 11, 1., 12., 495, 0.5, 50);
  fOutputWeight->Add(fHistMCGen_Xic0Pt_weight);
  fOutputWeight->Add(f2DHistMCRec_Xic0Pt_weight);
  PostData(6, fOutputWeight);

  return;
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::UserExec(Option_t *)
{
  // user exec
  // this function is called once for each event
  // the manager will take care of reading the events from file, and with the static function InputEvent() you 
  // have access to the current event. 
  // once you return from the UserExec function, the manager will retrieve the next event from the chain

  if (!fInputEvent) { // if the event is empty (getting it failed) skip this event
    AliError("NO EVENT FOUND!");
    return;
  }
  AliAODEvent* AODEvent = dynamic_cast<AliAODEvent*>(fInputEvent);    // get an event (called AODEvent) from the input file
                                                        // there's another event format (ESD) which works in a similar way
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with aod's

  fHistEvents->Fill(1);

  //--------------------------------------------------------------
  // First check if the event has magnetic field and proper vertex
  //--------------------------------------------------------------

  fBzkG = (Double_t)AODEvent->GetMagneticField();
  if (TMath::Abs(fBzkG)<0.001) return;
  KFParticle::SetField(fBzkG);

  fpVtx = (AliAODVertex*)AODEvent->GetPrimaryVertex();
  if (!fpVtx) return;
  fHistEvents->Fill(2);

  fCounter->StoreEvent(AODEvent,fAnaCuts,fIsMC);

  //------------------------------------------------
  // MC analysis setting                                                                    
  //------------------------------------------------

  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if (fIsMC) {
    fMCEvent = MCEvent(); // get the corresponding MC event fMCEvent
    if (!fMCEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }

    // MC array need for maching
    mcArray = dynamic_cast<TClonesArray*>(AODEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if ( !mcArray ) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    fHistEvents->Fill(7); // number of MC array exist

    // load MC header
    mcHeader = (AliAODMCHeader*)AODEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if ( !mcHeader ) {
      AliError("AliAnalysisTaskSEXicZero2XiPifromKFP::UserExec: MC header branch not found!\n");
      return;
    }
    fHistEvents->Fill(8); // number of MC header exist

    Double_t zMCvtx = mcHeader->GetVtxZ();
    if ( TMath::Abs(zMCvtx) > fAnaCuts->GetMaxVtxZ() ) {
      AliDebug(2,Form("Event rejected: fabs(zVtxMC)=%f > fAnaCuts->GetMaxVtxZ()=%f", zMCvtx, fAnaCuts->GetMaxVtxZ()));
      return;
    } else {
      fHistEvents->Fill(18);
    }
    if ((TMath::Abs(zMCvtx) < fAnaCuts->GetMaxVtxZ()) && (!fAnaCuts->IsEventRejectedDuePhysicsSelection()) && (!fAnaCuts->IsEventRejectedDueToTrigger())) {
      Bool_t selevt = MakeMCAnalysis(mcArray);
      if(!selevt) return;
    }
  }


  //------------------------------------------------
  // Event selection
  //------------------------------------------------
  Bool_t IsTriggerNotOK = fAnaCuts->IsEventRejectedDueToTrigger();
  Bool_t IsPhysSelNotOK = fAnaCuts->IsEventRejectedDuePhysicsSelection();
  Bool_t IsNoVertex = fAnaCuts->IsEventRejectedDueToNotRecoVertex();
  if( !IsTriggerNotOK && !IsPhysSelNotOK && !IsNoVertex && fabs(fpVtx->GetZ())<fAnaCuts->GetMaxVtxZ() ) fHistEvents->Fill(3);

  Bool_t IsEventSelected = fAnaCuts->IsEventSelected(AODEvent);
  if(!IsEventSelected) {
//    cout<<"Why: "<<fAnaCuts->GetWhyRejection()<<endl;
    return;
  }
  fHistEvents->Fill(4);


  Bool_t IsMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB)==(AliVEvent::kMB);
  Bool_t IsSemi = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kSemiCentral)==(AliVEvent::kSemiCentral);
  Bool_t IsCent = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kCentral)==(AliVEvent::kCentral);
  Bool_t IsINT7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7)==(AliVEvent::kINT7);
  Bool_t IsEMC7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMC7)==(AliVEvent::kEMC7);
  if(IsMB) fHTrigger->Fill(1);
  if(IsSemi) fHTrigger->Fill(2);
  if(IsCent) fHTrigger->Fill(3);
  if(IsINT7) fHTrigger->Fill(4);
  if(IsEMC7) fHTrigger->Fill(5);
  if(IsMB||IsSemi||IsCent) fHTrigger->Fill(7);
  if(IsINT7||IsEMC7) fHTrigger->Fill(8);
  if(IsMB&&IsSemi) fHTrigger->Fill(10);
  if(IsMB&&IsCent) fHTrigger->Fill(11);
  if(IsINT7&&IsEMC7) fHTrigger->Fill(12);

//  AliCentrality *cent = AODEvent->GetCentrality();
//  Float_t Centrality = cent->GetCentralityPercentile("V0M");
//  fHCentrality->Fill(Centrality);

  //------------------------------------------------
  // Check if the event has v0 candidate
  //------------------------------------------------
  Int_t num_v0 = AODEvent->GetNumberOfV0s();
  if (num_v0>0) fHistEvents->Fill(5);

  //------------------------------------------------
  // Check if the event has cascade candidate
  //------------------------------------------------
  Int_t num_casc = AODEvent->GetNumberOfCascades();
  if (num_casc<=0) return;
  fHistEvents->Fill(6);

  // set primary vertex
  KFPVertex pVertex;
  Double_t pos[3],cov[6];
  fpVtx->GetXYZ(pos);
  if ( fabs(pos[2])>10. ) return; // vertex cut on z-axis direction
  fpVtx->GetCovarianceMatrix(cov);
//  if ( !CheckVertexCov(fpVtx) ) cout << "Vertex Cov. is wrong!!!" << endl;
  pVertex.SetXYZ((Float_t)pos[0], (Float_t)pos[1], (Float_t)pos[2]);
  Float_t covF[6];
  for (Int_t i=0; i<6; i++) { covF[i] = (Float_t)cov[i]; }
  pVertex.SetCovarianceMatrix(covF);
  pVertex.SetChi2(fpVtx->GetChi2());
  pVertex.SetNDF(fpVtx->GetNDF());
  pVertex.SetNContributors(fpVtx->GetNContributors());

  KFParticle PV(pVertex);

  fHistPVx->Fill(pos[0]);
  fHistPVy->Fill(pos[1]);
  fHistPVz->Fill(pos[2]);

  if(!fAnaCuts) return;

  FillEventROOTObjects();

/*
    //--------------------------------------------------------------------------------------------------------
    // create a translation table: fAodTrackInd(mcTrackIndex) = aodTrackIndex, or = -1 if there is no aodTrack
    //--------------------------------------------------------------------------------------------------------

    fAodTrackInd.resize(0);
    fAodTrackInd.resize(fMCEvent->GetNumberOfTracks(), -1);

    Int_t nTracks = AODEvent->GetNumberOfTracks();
//  cout << "number of AOD tracks == " << nTracks << endl;
    
    // loop over all aod tracks and fill fAodTrackInd
    for(Int_t i=0; i < nTracks; i++) {                 // loop ove rall these tracks
      AliAODTrack* track = static_cast<AliAODTrack*>(AODEvent->GetTrack(i));         // get a track (type AliAODTrack) from the event
      if( !track || !track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA) || track->GetID()<0 ) continue;                            // if we failed, skip this track
//      cout << "Track Label:" << track->GetLabel() << endl;
//      cout << "i:" << i << endl;
//      cout << "Track ID:" << track->GetID() << endl;
      if (track->GetLabel()<0) continue;
      fAodTrackInd[track->GetLabel()]=i;
    } // continue until all the tracks are processed

    //--------------------------------------------------------------
    // we loop over all mc particles
    // check whether it is a Xi with two daughters
    // if the corresponding aod tracks exist
    // find the mother via kfparticle
    //--------------------------------------------------------------

    // loop over all primary MC particle
//    cout << "number of MC tracks == " << fMCEvent->GetNumberOfTracks() << endl;
//    cout << "number of MC particles == " << mcArray->GetEntriesFast() << endl;
    for (Long_t iPart=0; iPart < mcArray->GetEntriesFast(); iPart++) {
      AliAODMCParticle* mcPart = static_cast<AliAODMCParticle*>(mcArray->At(iPart));
      if (!mcPart) continue;
//  cout << "PDG CODE = " << mcPart->GetPdgCode() << endl;

      Int_t pdgDaughter[2] = {3312, 211};
      Int_t daughterIndex[2] = {0, 1};
//    MatchToMCXic0(mcArray, 4132, 2, daughterIndex, pdgDaughter);

    }
  }
*/


//------------------------------------------------
// Main analysis done in this function
//------------------------------------------------
  
  fPID = fInputHandler->GetPIDResponse();
//  MakeAnaXicZeroFromV0(AODEvent, mcArray, PV);
  MakeAnaXicZeroFromCasc(AODEvent, mcArray, PV);

  PostData(2, fCounter);
  PostData(3, fTree_Event);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
  PostData(4, fTree_Xic0);
  PostData(5, fTree_Xic0MCGen);

  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
/*
    TCanvas *c1 = new TCanvas();
    fHistDecayLLambda->Draw();
    TLine *lLPDG = new TLine(7.89, 0, 7.89, 1e10);
    lLPDG->SetLineColor(2);
    lLPDG->Draw();

    TCanvas *c2 = new TCanvas();
    fHistDecayLXiMinus->Draw();
    TLine *lXiPDG = new TLine(4.91, 0, 4.91, 1e10);
    lXiPDG->SetLineColor(2);
    lXiPDG->Draw();
*/
    return;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEXicZero2XiPifromKFP::MakeMCAnalysis(TClonesArray *mcArray)
{
  // Analyse AliAODMCParticle
  
  Double_t Num_Xic0  = 0;
  Int_t nmcpart = mcArray->GetEntriesFast();
  Int_t NDaughters = 2;

  for(Int_t i=0;i<nmcpart;i++) {
    AliAODMCParticle *mcpart = NULL;
    mcpart = (AliAODMCParticle*) mcArray->At(i);

    fHistMCpdg_All->Fill(mcpart->GetPdgCode());



    /*
    // =============== check generated particles in MC =====================
    
    if ( TMath::Abs(mcpart->GetPdgCode())==4232 ) { // Xic+
//      cout << "Particle: " << mcpart->GetPdgCode() << endl;
//      cout << "Daughters: " << endl;
      for(Int_t idau=mcpart->GetDaughterFirst();idau<=mcpart->GetDaughterLast();idau++) {
        if(idau<0) break;
        AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
        fHistMCpdg_Dau_XicPM->Fill(mcdau->GetPdgCode());
//        cout << mcdau->GetPdgCode() << endl;
      }
    }

    if ( TMath::Abs(mcpart->GetPdgCode())==4132 ) { // Xic0
      cout << "==Xic0 Flag: " << mcpart->GetFlag() << "==" << endl;
      cout << "Daughters: " << endl;
      for(Int_t idau=mcpart->GetDaughterFirst();idau<=mcpart->GetDaughterLast();idau++) {
        if(idau<0) break;
        AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
        fHistMCpdg_Dau_XicZero->Fill(mcdau->GetPdgCode());
//        cout << mcdau->GetPdgCode() << endl;
//        cout << "Flag:" << mcdau->GetFlag() << endl;
        if ( TMath::Abs(mcdau->GetPdgCode())==211 ) {
          cout << "  Pi from Xic0 Flag:" << mcdau->GetFlag() << endl;
        }

        if ( TMath::Abs(mcdau->GetPdgCode())==3312 ) {
          cout << "  Xi from Xic0 Flag:" << mcdau->GetFlag() << endl;
          for(Int_t idau=mcdau->GetDaughterFirst();idau<=mcdau->GetDaughterLast();idau++) {
            if(idau<0) break;
            AliAODMCParticle *mcdau1 = (AliAODMCParticle*) mcArray->At(idau);
            if ( TMath::Abs(mcdau1->GetPdgCode())==211 ) {
              cout << "    Pi from Xi Flag:" << mcdau1->GetFlag() << endl;
            }
            if ( TMath::Abs(mcdau1->GetPdgCode())==3122 ) {
              cout << "    Lambda from Xi Flag:" << mcdau1->GetFlag() << endl;
              for(Int_t idau=mcdau1->GetDaughterFirst();idau<=mcdau1->GetDaughterLast();idau++) {
                if(idau<0) break;
                AliAODMCParticle *mcdau2 = (AliAODMCParticle*) mcArray->At(idau);
                if ( TMath::Abs(mcdau2->GetPdgCode())==211 ) {
                  cout << "      Pi from Lambda Flag:" << mcdau2->GetFlag() << endl;
                }
                if ( TMath::Abs(mcdau2->GetPdgCode())==2212 ) {
                  cout << "      Proton from Lambda Flag:" << mcdau2->GetFlag() << endl;
                }
              }
            }
          }
        }
      }
    }
    // ======================================================================
    */


//    if ( mcpart->GetNDaughters() != NDaughters ) continue;


    // ======================================= Xic0 ====================================================
    if ( TMath::Abs(mcpart->GetPdgCode())==4132 && mcpart->GetNDaughters()==NDaughters ) { // 4132: Xic0
      Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcpart,kTRUE);
//      if (CheckOrigin==0) continue;
      Bool_t pifromXic_flag = kFALSE;
      Bool_t xi_flag = kFALSE;
      AliAODMCParticle *mcpipart = NULL;
      AliAODMCParticle *mccascpart = NULL;
      AliAODMCParticle *mcv0part = NULL;
      for(Int_t idau=mcpart->GetDaughterFirst();idau<=mcpart->GetDaughterLast();idau++) {
        if(idau<0) break;
        AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
        if(TMath::Abs(mcdau->GetPdgCode())==211) { // 211: pion
          pifromXic_flag = kTRUE;
          mcpipart = mcdau;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==3312) { // 3312: Xi
          xi_flag = kTRUE;
          mccascpart = mcdau;
        }
      }

      Int_t decaytype = -9999;
      if ( pifromXic_flag &&  xi_flag ) decaytype = 0;
      if ( pifromXic_flag &&  xi_flag) fHistMCXicZeroDecayType->Fill(1);
      if (!pifromXic_flag &&  xi_flag) fHistMCXicZeroDecayType->Fill(2);
      if ( pifromXic_flag && !xi_flag) fHistMCXicZeroDecayType->Fill(3);
      if (!pifromXic_flag && !xi_flag) fHistMCXicZeroDecayType->Fill(4);

      if (decaytype==0) {
        FillTreeGenXic0(mcpart, CheckOrigin);
        if (fabs(mcpart->Y())<0.8) fHistMCGen_Xic0Pt_weight->Fill(mcpart->Pt(), fWeight->Eval(mcpart->Pt()));
      }
    } // for Xic0
    // ======================================= Xi ====================================================
    if ( TMath::Abs(mcpart->GetPdgCode())==3312 && mcpart->GetNDaughters()==NDaughters ) { // 3312: Xi
      Bool_t pifromXi_flag = kFALSE;
      Bool_t v0_flag = kFALSE;
      Bool_t pifromLam_flag = kFALSE;
      Bool_t prfromLam_flag = kFALSE;
      AliAODMCParticle *mcv0part = NULL;
      for(Int_t idau=mcpart->GetDaughterFirst();idau<=mcpart->GetDaughterLast();idau++) {
        if(idau<0) break;
        AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
        if(TMath::Abs(mcdau->GetPdgCode())==211){ // 211: pion
          pifromXi_flag = kTRUE;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==3122 && mcdau->GetNDaughters()==NDaughters) { // 3122: lambda
          v0_flag = kTRUE;
          mcv0part = mcdau;
          for(Int_t jdau=mcv0part->GetDaughterFirst();jdau<=mcv0part->GetDaughterLast();jdau++) {
            if (jdau<0) break;
            AliAODMCParticle *mcdau_Lam = (AliAODMCParticle*) mcArray->At(jdau);
            if(TMath::Abs(mcdau_Lam->GetPdgCode())==211) pifromLam_flag = kTRUE;
            if(TMath::Abs(mcdau_Lam->GetPdgCode())==2212) prfromLam_flag = kTRUE;
          }
        }
      }

      Int_t decaytype = -9999;
      if ( pifromXi_flag &&  v0_flag && pifromLam_flag && prfromLam_flag) decaytype = 0;
      if ( pifromXi_flag &&  v0_flag) fHistMCXiDecayType->Fill(1);
      if (!pifromXi_flag &&  v0_flag) fHistMCXiDecayType->Fill(2);
      if ( pifromXi_flag && !v0_flag) fHistMCXiDecayType->Fill(3);
      if (!pifromXi_flag && !v0_flag) fHistMCXiDecayType->Fill(4);

      if (decaytype==0) {
//        FillMCCascROOTObjects(mcpart, mcArray);
        if (mcpart->GetPdgCode()>0) {
          fCounterGen_Cuts_XiMinus->Fill(1);
          if ( mcpart->IsPrimary() ) fCounterGen_Cuts_XiMinus->Fill(2);
          if ( mcpart->IsPhysicalPrimary() ) {
            fCounterGen_Cuts_XiMinus->Fill(3);
            fHistMCGen_XiMinus_Pt->Fill(mcpart->Pt());
            if (TMath::Abs(mcpart->Y())<0.5) {
              fHistMCGen_XiMinus_Pt_wYcut->Fill(mcpart->Pt());
            }
          }
          if ( mcpart->IsSecondaryFromWeakDecay() ) fCounterGen_Cuts_XiMinus->Fill(4);
          if ( mcpart->IsSecondaryFromMaterial() ) fCounterGen_Cuts_XiMinus->Fill(5);
          if ( mcpart->IsFromSubsidiaryEvent() ) fCounterGen_Cuts_XiMinus->Fill(6);
        }
        if (mcpart->GetPdgCode()<0) {
          fCounterGen_Cuts_XiPlus->Fill(1);
          if ( mcpart->IsPrimary() ) fCounterGen_Cuts_XiPlus->Fill(2);
          if ( mcpart->IsPhysicalPrimary() ) {
            fCounterGen_Cuts_XiPlus->Fill(3);
            fHistMCGen_XiPlus_Pt->Fill(mcpart->Pt());
            if (TMath::Abs(mcpart->Y())<0.5) {
              fHistMCGen_XiPlus_Pt_wYcut->Fill(mcpart->Pt());
            }
          }
          if ( mcpart->IsSecondaryFromWeakDecay() ) fCounterGen_Cuts_XiPlus->Fill(4);
          if ( mcpart->IsSecondaryFromMaterial() ) fCounterGen_Cuts_XiPlus->Fill(5);
          if ( mcpart->IsFromSubsidiaryEvent() ) fCounterGen_Cuts_XiPlus->Fill(6);
        }
      }
    } // for Xi
    // ======================================= Lambda ====================================================
    if ( TMath::Abs(mcpart->GetPdgCode())==3122 && mcpart->GetNDaughters()==NDaughters ) { // 3122: Lambda
      Bool_t pi_flag = kFALSE;
      Bool_t pr_flag = kFALSE;
      AliAODMCParticle *mcpi = NULL;
      AliAODMCParticle *mcpr = NULL;
      for(Int_t idau=mcpart->GetDaughterFirst();idau<=mcpart->GetDaughterLast();idau++) {
        if (idau<0) break;
        AliAODMCParticle *mcdau = (AliAODMCParticle*) mcArray->At(idau);
        if(TMath::Abs(mcdau->GetPdgCode())==211){
          pi_flag = kTRUE;
          mcpi = mcdau;
        }
        if(TMath::Abs(mcdau->GetPdgCode())==2212){
          pr_flag = kTRUE;
          mcpr = mcdau;
        }
      }
      Int_t decaytype = -9999;
      if ( pi_flag && pr_flag) decaytype = 0;
      if (decaytype==0) {
//        FillMCV0ROOTObjects(mcpart, mcArray);
        if (mcpart->GetPdgCode()>0) { // Lambda
          if (fabs(mcpart->Y())<0.8) fCounterGen_Cuts_Lambda->Fill(1);
//================================== pri. and sec. lambda =====================================================
          if ( mcpart->IsPrimary() && fabs(mcpart->Y())<0.8 ) fCounterGen_Cuts_Lambda->Fill(2);
          if ( mcpart->IsPhysicalPrimary() ) {
            fGenHistRapidity_Lambda->Fill(mcpart->Y());
            if (fabs(mcpart->Y())<0.8) {
              fCounterGen_Cuts_Lambda->Fill(3);
              fHistMCGen_Lambda_Pt->Fill(mcpart->Pt());
            }
            if (TMath::Abs(mcpart->Y())<0.5) fHistMCGen_Lambda_Pt_wYcut->Fill(mcpart->Pt());
          }
          if ( mcpart->IsSecondaryFromWeakDecay() && fabs(mcpart->Y())<0.8 ) fCounterGen_Cuts_Lambda->Fill(4);
          if ( mcpart->IsSecondaryFromMaterial() && fabs(mcpart->Y())<0.8 ) fCounterGen_Cuts_Lambda->Fill(5);
          if ( mcpart->IsFromSubsidiaryEvent() && fabs(mcpart->Y())<0.8 ) fCounterGen_Cuts_Lambda->Fill(6);

          /*
          AliAODMCParticle *mother = mcpart;
          Int_t IndexMother = -1;
          while ( mother->GetMother()>=0 ) {
             IndexMother = mother->GetMother();
             mother = (AliAODMCParticle*)mcArray->At(IndexMother);
             if (!mother) {
               printf("no MC mother particle\n");
               break;
             }
             Int_t pdgMother = mother->GetPdgCode();
             cout << "Index Mother==" << IndexMother << ", " << "PDG Mother==" << pdgMother << endl;
             if ( pdgMother== 3122 ) {
             }
          }
          cout << "***************************" << endl;
          */
        }

        if (mcpart->GetPdgCode()<0) { // Anti-Lambda
          if (fabs(mcpart->Y())<0.8) fCounterGen_Cuts_AntiLambda->Fill(1);
//================================== pri. and sec. anti-lambda =====================================================
          if ( mcpart->IsPrimary() && fabs(mcpart->Y())<0.8 ) fCounterGen_Cuts_AntiLambda->Fill(2);
          if ( mcpart->IsPhysicalPrimary() ) {
            fGenHistRapidity_AntiLambda->Fill(mcpart->Y());
            if (fabs(mcpart->Y())<0.8) {
              fCounterGen_Cuts_AntiLambda->Fill(3);
              fHistMCGen_AntiLambda_Pt->Fill(mcpart->Pt());
            }
            if (TMath::Abs(mcpart->Y())<0.5) fHistMCGen_AntiLambda_Pt_wYcut->Fill(mcpart->Pt());
          }
          if ( mcpart->IsSecondaryFromWeakDecay() && fabs(mcpart->Y())<0.8 ) fCounterGen_Cuts_AntiLambda->Fill(4);
          if ( mcpart->IsSecondaryFromMaterial() && fabs(mcpart->Y())<0.8 ) fCounterGen_Cuts_AntiLambda->Fill(5);
          if ( mcpart->IsFromSubsidiaryEvent() && fabs(mcpart->Y())<0.8 ) fCounterGen_Cuts_AntiLambda->Fill(6);
        }

      } // decaytype=0
    } // for lambda
  } // all loop of MC particles

  return kTRUE;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::FillTreeGenXic0(AliAODMCParticle *mcpart, Int_t CheckOrigin)
{
  // Fill histograms or tree depending

  for(Int_t i=0;i<3;i++){
    fVar_Xic0MCGen[i] = -9999.;
  }

  fVar_Xic0MCGen[0] = mcpart->Y();
  fVar_Xic0MCGen[1] = mcpart->Pt();
  fVar_Xic0MCGen[2] = CheckOrigin;

  if (fWriteXic0MCGenTree && fVar_Xic0MCGen[1]>0.9999) fTree_Xic0MCGen->Fill();

//  fVar_Xic0MCGen[ 0] = fCentrality;
//  fVar_Xic0MCGen[ 1] = decaytype;
//  if (mcpart->IsPrimary() && (!mcpart->IsPhysicalPrimary())) fVar_Xic0MCGen[2] = 1;
//  if (mcpart->IsPhysicalPrimary()) fVar_Xic0MCGen[2] = 2;
//  if (mcpart->IsSecondaryFromWeakDecay()) fVar_Xic0MCGen[2] = 3;
//  if (mcpart->IsSecondaryFromMaterial()) fVar_Xic0MCGen[2] = 4;
//  if (mcpart->IsFromSubsidiaryEvent()) fVar_Xic0MCGen[2] = 5;
//  fVar_Xic0MCGen[ 3] = mcpart->Eta();
//  fVar_Xic0MCGen[ 4] = mcpart->Y();
//  fVar_Xic0MCGen[ 5] = mcpart->Px();
//  fVar_Xic0MCGen[ 6] = mcpart->Py();
//  fVar_Xic0MCGen[ 7] = mcpart->Pz();
//  fVar_Xic0MCGen[ 8] = mcpipart->Px();
//  fVar_Xic0MCGen[ 9] = mcpipart->Py();
//  fVar_Xic0MCGen[10] = mcpipart->Pz();
//  fVar_Xic0MCGen[11] = mccascpart->Px();
//  fVar_Xic0MCGen[12] = mccascpart->Py();
//  fVar_Xic0MCGen[13] = mccascpart->Pz();
//  fVar_Xic0MCGen[14] = mcpart->GetPdgCode();
//  fVar_Xic0MCGen[15] = mcpipart->GetPdgCode();
//  fVar_Xic0MCGen[16] = mccascpart->GetPdgCode();
//  fVar_Xic0MCGen[17] = fRunNumber;
//  fVar_Xic0MCGen[18] = fEvNumberCounter;
//  fVar_Xic0MCGen[19] = CheckOrigin;

  /*
  const Double_t massPion = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t massXi   = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  Double_t pipx = mcpipart->Px();
  Double_t pipy = mcpipart->Py();
  Double_t pipz = mcpipart->Pz();
//  Double_t piE  = sqrt(pipx*pipx+pipy*pipy+pipz*pipz+0.000511*0.000511);
  Double_t piE  = sqrt(pipx*pipx+pipy*pipy+pipz*pipz+massPion*massPion);
  Double_t cascpx = mccascpart->Px();
  Double_t cascpy = mccascpart->Py();
  Double_t cascpz = mccascpart->Pz();
  Double_t cascE  = sqrt(cascpx*cascpx+cascpy*cascpy+cascpz*cascpz+massXi*massXi);

  Double_t InvMassPiXi = sqrt(pow(piE+cascE,2)-pow(pipx+cascpx,2)-pow(pipy+cascpy,2)-pow(pipz+cascpz,2));

  Double_t contXicZeroMC[3];
  contXicZeroMC[0] = mcpart->Pt();
  contXicZeroMC[1] = mcpart->Y();
  contXicZeroMC[2] = fCentrality;

  Double_t contPionMC[3];
  contPionMC[0] = mcpipart->Pt();
  contPionMC[1] = mcpipart->Eta();
  contPionMC[2] = fCentrality;

  Double_t contXiMC[3];
  contXiMC[0] = mccascpart->Pt();
  contXiMC[1] = mccascpart->Y();
  contXiMC[2] = fCentrality;

  Double_t contPiXiMassMCGen[3];
  contPiXiMassMCGen[0] = InvMassPiXi;
  contPiXiMassMCGen[1] = mcpart->Pt();
  contPiXiMassMCGen[2] = fCentrality;

  Double_t contPiXiMassvsPiPtMCGen[3];
  contPiXiMassvsPiPtMCGen[0] = InvMassPiXi;
  contPiXiMassvsPiPtMCGen[1] = mcpipart->Pt();
  contPiXiMassvsPiPtMCGen[2] = fCentrality;

  if (decaytype==0) {
    fHistMCGen_XicZeroTot->Fill(contXicZeroMC);
    if (mcpart->GetPdgCode()>0) fHistMCGen_XicZero->Fill(contXicZeroMC);
    if (mcpart->GetPdgCode()<0) fHistMCGen_AntiXicZero->Fill(contXicZeroMC);
    fHistMCGen_PionTot->Fill(contPionMC);
    if (mcpipart->GetPdgCode()>0) fHistMCGen_PionPlus->Fill(contPionMC);
    if (mcpipart->GetPdgCode()<0) fHistMCGen_PionMinus->Fill(contPionMC);
//    fHistMCGen_XiTot->Fill(contXiMC);
//    if (mccascpart->GetPdgCode()>0) fHistMCGen_XiMinus->Fill(contXiMC);
//    if (mccascpart->GetPdgCode()<0) fHistMCGen_XiPlus->Fill(contXiMC);
    fHistMCGen_PiXiInvMass->Fill(contPiXiMassMCGen);
    if (fabs(mcpipart->Eta())<fAnaCuts->GetProdTrackEtaRange()) {
      fHistMCGen_PiXiMassvsPiPt->Fill(contPiXiMassvsPiPtMCGen);
      if (mcpipart->GetPdgCode()>0) fHistMCGen_PiXiMassvsPiPt_PionPlus->Fill(contPiXiMassvsPiPtMCGen);
      if (mcpipart->GetPdgCode()<0) fHistMCGen_PiXiMassvsPiPt_PionMinus->Fill(contPiXiMassvsPiPtMCGen);
    }
  }
  */
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::MakeAnaXicZeroFromV0(AliAODEvent *AODEvent, TClonesArray *mcArray, KFParticle PV)
{
  // Main analysis called from "UserExec"

//  std::cout.setf(std::ios::fixed);
//  std::cout.setf(std::ios::showpoint);
//  std::cout.precision(3);

  // set the magnetic field
  KFParticle::SetField(fBzkG);

  const UInt_t nV0s = AODEvent->GetNumberOfV0s();
//  if (nV0s==0) return;

  Double_t xyzP[3], xyzN[3];
  Double_t xvyvzvP[3], xvyvzvN[3];
  Double_t covP[21], covN[21];
  const Int_t NDaughters = 2;
  const Float_t massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  const Float_t massXi = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  const Float_t massK0S = TDatabasePDG::Instance()->GetParticle(310)->Mass();

  const UInt_t nTracks = AODEvent->GetNumberOfTracks();
//  if (nTracks<4) return;

  // select good candidates for pion
  AliAODTrack *trackP[nTracks], *trackN[nTracks];
  Int_t flag_trkP = 0, flag_trkN = 0;

  for (UInt_t itrk=0; itrk<nTracks; itrk++) {
    AliAODTrack *trk = static_cast<AliAODTrack*>(AODEvent->GetTrack(itrk));
    Double_t covtest[21];
//    if ( !trk || trk->GetID()<0 || !trk->GetCovarianceXYZPxPyPz(covtest) || !CheckTrackCov(trk) || !fAnaCuts->SinglePionPoolCuts(trk) ) continue;
    if ( !trk || trk->GetID()<0 || !trk->GetCovarianceXYZPxPyPz(covtest) || !CheckTrackCov(trk) ) continue;

    if  (trk->Charge() > 0 ) {
      trackP[flag_trkP] = trk;
      flag_trkP++;
    }
    if  (trk->Charge() < 0 ) {
      trackN[flag_trkN] = trk;
      flag_trkN++;
    }
  }

  // at least one pion+ and one pion-
//  if ( flag_trkP < 1 || flag_trkN<1 ) return;

  // loop for v0
  for (UInt_t itrk=0; itrk<nV0s; itrk++) {
    // check v0
    AliAODv0 *v0 = AODEvent->GetV0(itrk);

//    if ( !v0 ) continue;
    if ( v0->GetOnFlyStatus() ) continue; // to select offline V0s

    fRecHistRapidity_Lambda_offline->Fill(v0->RapLambda());

    if ( fabs(v0->RapLambda())<0.8 ) {
      fCounterRec_Cuts_Lambda->Fill(1);
      fCounterRec_Cuts_AntiLambda->Fill(1);

      if ( fIsMC ) {
        AliAODTrack *trackP = (AliAODTrack*) (v0->GetDaughter(0));
        AliAODTrack *trackN = (AliAODTrack*) (v0->GetDaughter(1));
        if (trackP->Charge()<0) {
          trackP = (AliAODTrack*) (v0->GetDaughter(1));
          trackN = (AliAODTrack*) (v0->GetDaughter(0));
        }
        Int_t lab_Lam     = MatchToMCLambda(trackP, trackN, mcArray);
        Int_t lab_AntiLam = MatchToMCAntiLambda(trackN, trackP, mcArray);
        if (lab_Lam>-1) {
          fCounterRecMC_Cuts_Lambda->Fill(1);
          f2DCounterRecMC_CutsVsPt_Lambda->Fill(1, sqrt(v0->Pt2V0()));
        }
        if (lab_AntiLam>-1) {
          fCounterRecMC_Cuts_AntiLambda->Fill(1);
          f2DCounterRecMC_CutsVsPt_AntiLambda->Fill(1, sqrt(v0->Pt2V0()));
        }
      }
    }

    if ( v0->GetNDaughters()!=2 ) continue;

    if ( fabs(v0->RapLambda())<0.8 ) {
      fCounterRec_Cuts_Lambda->Fill(2);
      fCounterRec_Cuts_AntiLambda->Fill(2);

      if ( fIsMC ) {
        AliAODTrack *trackP = (AliAODTrack*) (v0->GetDaughter(0));
        AliAODTrack *trackN = (AliAODTrack*) (v0->GetDaughter(1));
        if (trackP->Charge()<0) {
          trackP = (AliAODTrack*) (v0->GetDaughter(1));
          trackN = (AliAODTrack*) (v0->GetDaughter(0));
        }
        Int_t lab_Lam     = MatchToMCLambda(trackP, trackN, mcArray);
        Int_t lab_AntiLam = MatchToMCAntiLambda(trackN, trackP, mcArray);
        if (lab_Lam>-1) {
          fCounterRecMC_Cuts_Lambda->Fill(2);
          f2DCounterRecMC_CutsVsPt_Lambda->Fill(2, sqrt(v0->Pt2V0()));
        }
        if (lab_AntiLam>-1) {
          fCounterRecMC_Cuts_AntiLambda->Fill(2);
          f2DCounterRecMC_CutsVsPt_AntiLambda->Fill(2, sqrt(v0->Pt2V0()));
        }
      }
    }

    // check daughters
    AliAODTrack *trk0 = (AliAODTrack*) (v0->GetDaughter(0));
    AliAODTrack *trk1 = (AliAODTrack*) (v0->GetDaughter(1));
    AliAODTrack *trkP=NULL, *trkN=NULL;
    if (!trk0 || !trk1 ) continue;
    if (trk0->Charge()>0) {
      trkP = (AliAODTrack*)v0->GetDaughter(0);
      trkN = (AliAODTrack*)v0->GetDaughter(1);
    }
    if (trk0->Charge()<0) {
      trkP = (AliAODTrack*)v0->GetDaughter(1);
      trkN = (AliAODTrack*)v0->GetDaughter(0);
    }

    if ( fabs(v0->RapLambda())<0.8 ) {
      fCounterRec_Cuts_Lambda->Fill(3);
      fCounterRec_Cuts_AntiLambda->Fill(3);

      if ( fIsMC ) {
        Int_t lab_Lam     = MatchToMCLambda(trkP, trkN, mcArray);
        Int_t lab_AntiLam = MatchToMCAntiLambda(trkN, trkP, mcArray);
        if (lab_Lam>-1) {
          fCounterRecMC_Cuts_Lambda->Fill(3);
          f2DCounterRecMC_CutsVsPt_Lambda->Fill(3, sqrt(v0->Pt2V0()));
        }
        if (lab_AntiLam>-1) {
          fCounterRecMC_Cuts_AntiLambda->Fill(3);
          f2DCounterRecMC_CutsVsPt_AntiLambda->Fill(3, sqrt(v0->Pt2V0()));
        }
      }
    }

    if ( !trkP->GetCovarianceXYZPxPyPz(covP) || !trkN->GetCovarianceXYZPxPyPz(covN) ) continue;

    if ( fabs(v0->RapLambda())<0.8 ) {
      fCounterRec_Cuts_Lambda->Fill(4);
      fCounterRec_Cuts_AntiLambda->Fill(4);

      if ( fIsMC ) {
        Int_t lab_Lam     = MatchToMCLambda(trkP, trkN, mcArray);
        Int_t lab_AntiLam = MatchToMCAntiLambda(trkN, trkP, mcArray);
        if (lab_Lam>-1) {
          fCounterRecMC_Cuts_Lambda->Fill(4);
          f2DCounterRecMC_CutsVsPt_Lambda->Fill(4, sqrt(v0->Pt2V0()));
        }
        if (lab_AntiLam>-1) {
          fCounterRecMC_Cuts_AntiLambda->Fill(4);
          f2DCounterRecMC_CutsVsPt_AntiLambda->Fill(4, sqrt(v0->Pt2V0()));
        }
      }
    }

    if ( !CheckTrackCov(trkP) || !CheckTrackCov(trkN) ) continue;

    if ( fabs(v0->RapLambda())<0.8 ) {
      fCounterRec_Cuts_Lambda->Fill(5);
      fCounterRec_Cuts_AntiLambda->Fill(5);

      if ( fIsMC ) {
        Int_t lab_Lam     = MatchToMCLambda(trkP, trkN, mcArray);
        Int_t lab_AntiLam = MatchToMCAntiLambda(trkN, trkP, mcArray);
        if (lab_Lam>-1) {
          fCounterRecMC_Cuts_Lambda->Fill(5);
          f2DCounterRecMC_CutsVsPt_Lambda->Fill(5, sqrt(v0->Pt2V0()));
        }
        if (lab_AntiLam>-1) {
          fCounterRecMC_Cuts_AntiLambda->Fill(5);
          f2DCounterRecMC_CutsVsPt_AntiLambda->Fill(5, sqrt(v0->Pt2V0()));
        }
      }
    }

    fHistMassLambda_woCut->Fill(v0->MassLambda());
    fHistMassAntiLambda_woCut->Fill(v0->MassAntiLambda());

    // v0 cut
    if ( !fAnaCuts->SingleV0LambdaTotCuts(v0) ) continue;

//    Float_t TPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
//    Float_t TOFnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTOF(track, AliPID::kElectron) : 1000;





    if (fAnaCuts->LambdaPIDCuts(v0) && fabs(v0->RapLambda())<0.8 ) {
      fRecHistRapidity_Lambda_wSTD->Fill(v0->RapLambda());
      fCounterRec_Cuts_Lambda->Fill(6);
      if ( fIsMC ) {
        Int_t lab_Lam     = MatchToMCLambda(trkP, trkN, mcArray);
        if (lab_Lam>-1) {
          fCounterRecMC_Cuts_Lambda->Fill(6);
          f2DCounterRecMC_CutsVsPt_Lambda->Fill(6, sqrt(v0->Pt2V0()));
        }
      }
    }
    if (fAnaCuts->AntiLambdaPIDCuts(v0) && fabs(v0->RapLambda())<0.8 ) {
      fRecHistRapidity_AntiLambda_wSTD->Fill(v0->RapLambda());
      fCounterRec_Cuts_AntiLambda->Fill(6);
      if ( fIsMC ) {
        Int_t lab_AntiLam = MatchToMCAntiLambda(trkN, trkP, mcArray);
        if (lab_AntiLam>-1) {
          fCounterRecMC_Cuts_AntiLambda->Fill(6);
          f2DCounterRecMC_CutsVsPt_AntiLambda->Fill(6, sqrt(v0->Pt2V0()));
        }
      }
    }

//    f2DHistChargeDaughters->Fill(trkP->Charge(), trkN->Charge());
//    fHistChargeV0->Fill(v0->Charge());
//    fHistNDaughterV0->Fill(v0->GetNDaughters());
//    fHistNProngV0->Fill(v0->GetNProngs());
//    fHistChargeFirstDaughter->Fill(trkP->Charge());
//    fHistChargeSecondDaughter->Fill(trkN->Charge());

    Double_t alpha_FirstDaugPos = (v0->MomPosAlongV0() - v0->MomNegAlongV0())/(v0->MomPosAlongV0() + v0->MomNegAlongV0());

    if (trk0->Charge()>0) {
      fHistMassLambda_wSTDv0Cut->Fill(v0->MassLambda());
      fHistMassAntiLambda_wSTDv0Cut->Fill(v0->MassAntiLambda());
    }

    if (trk0->Charge()<0) {
      alpha_FirstDaugPos = (v0->MomNegAlongV0() - v0->MomPosAlongV0())/(v0->MomPosAlongV0() + v0->MomNegAlongV0());
      fHistMassLambda_wSTDv0Cut->Fill(v0->MassAntiLambda());
      fHistMassAntiLambda_wSTDv0Cut->Fill(v0->MassLambda());
    }

    f2DHistArmenterosPodolanski_FirstDaugPos->Fill(alpha_FirstDaugPos, v0->PtArmV0());
    f2DHistV0XY_FirstDaugPos->Fill(v0->DecayVertexV0X(), v0->DecayVertexV0Y());

    trkP->GetXYZ(xyzP);
    trkP->XvYvZv(xvyvzvP);
    fHistXtrkP->Fill(xyzP[0]);
    fHistYtrkP->Fill(xyzP[1]);
    fHistZtrkP->Fill(xyzP[2]);
    fHistXtrkP_XYZv->Fill(xvyvzvP[0]);
    fHistYtrkP_XYZv->Fill(xvyvzvP[1]);
    fHistZtrkP_XYZv->Fill(xvyvzvP[2]);
    trkN->GetXYZ(xyzN);
    trkN->XvYvZv(xvyvzvN);
    fHistXtrkN->Fill(xyzN[0]);
    fHistYtrkN->Fill(xyzN[1]);
    fHistZtrkN->Fill(xyzN[2]);
    fHistXtrkN_XYZv->Fill(xvyvzvN[0]);
    fHistYtrkN_XYZv->Fill(xvyvzvN[1]);
    fHistZtrkN_XYZv->Fill(xvyvzvN[2]);

    KFParticle kfpProton     = CreateKFParticleFromAODtrack(trkP, 2212);
    KFParticle kfpAntiProton = CreateKFParticleFromAODtrack(trkN, -2212);
    KFParticle kfpPionMinus  = CreateKFParticleFromAODtrack(trkN, -211);
    KFParticle kfpPionPlus   = CreateKFParticleFromAODtrack(trkP, 211);
    KFParticle kfpLambda;
    KFParticle kfpAntiLambda;
    KFParticle kfpK0Short;
//    kfpLambda.SetPDG(3122);
//    kfpAntiLambda.SetPDG(-3122);
    const KFParticle *vDaughters[2]     = {&kfpProton, &kfpPionMinus};
    const KFParticle *vAntiDaughters[2] = {&kfpPionPlus, &kfpAntiProton};
    const KFParticle *vk0sDaughters[2]  = {&kfpPionPlus, &kfpPionMinus};

    kfpLambda.Construct(vDaughters, NDaughters);
    kfpAntiLambda.Construct(vAntiDaughters, NDaughters);
    kfpK0Short.Construct(vk0sDaughters, NDaughters);

    Float_t massK0S_Rec, err_massK0S;
    kfpK0Short.GetMass(massK0S_Rec, err_massK0S);
    fHistMassK0S->Fill(massK0S_Rec);

//    if ( fabs(massK0S_Rec-massK0S)<=fAnaCuts->GetProdMassTolKs0() ) continue;

    Float_t massLambda_Rec, err_massLambda;
    kfpLambda.GetMass(massLambda_Rec, err_massLambda);
    Float_t massAntiLambda_Rec, err_massAntiLambda;
    kfpAntiLambda.GetMass(massAntiLambda_Rec, err_massAntiLambda);

    // check rapidity of lambda && anti-lambda
    if ( TMath::Abs(kfpLambda.GetE())<=TMath::Abs(kfpLambda.GetPz()) && TMath::Abs(kfpAntiLambda.GetE())<=TMath::Abs(kfpAntiLambda.GetPz()) ) continue;

    // chi2>0 && NDF>0 for selecting Lambda and Anti-Lambda candidates
    if ( (kfpLambda.GetNDF()<=0 || kfpLambda.GetChi2()<=0) && (kfpAntiLambda.GetNDF()<=0 || kfpAntiLambda.GetChi2()<=0) ) continue;

    // check cov. of Lambda and Anti-Lambda candidates
    if ( !CheckKFParticleCov(kfpLambda) && !CheckKFParticleCov(kfpAntiLambda) ) continue;

    // err_mass>0 of Lambda and Anti-Lambda candidates
    if ( err_massLambda<=0 && err_massAntiLambda<=0 ) continue;

    if (fAnaCuts->LambdaPIDCuts(v0) && TMath::Abs(kfpLambda.GetE())>TMath::Abs(kfpLambda.GetPz()) && kfpLambda.GetNDF()>0 && kfpLambda.GetChi2()>0 && CheckKFParticleCov(kfpLambda) && err_massLambda>0) {
      if (fabs(kfpLambda.GetRapidity())<0.8) {
        fCounterRec_Cuts_Lambda->Fill(7);
        fHistChi2ndfLambda->Fill(kfpLambda.GetChi2()/kfpLambda.GetNDF());
        fHistProbLambda->Fill(TMath::Prob(kfpLambda.GetChi2(), kfpLambda.GetNDF()));
        if (TMath::Abs(kfpLambda.GetRapidity())<0.5) fCounterRec_Cuts_Lambda->Fill(8);
        if ( fIsMC ) {
          Int_t lab_Lam     = MatchToMCLambda(trkP, trkN, mcArray);
          if (lab_Lam>-1) {
            fCounterRecMC_Cuts_Lambda->Fill(7);
            f2DCounterRecMC_CutsVsPt_Lambda->Fill(7, kfpLambda.GetPt());
            if (TMath::Abs(kfpLambda.GetRapidity())<0.5) {
              fCounterRecMC_Cuts_Lambda->Fill(8);
              f2DCounterRecMC_CutsVsPt_Lambda->Fill(8, kfpLambda.GetPt());
            }
          }
        }
      }
    }
    if (fAnaCuts->AntiLambdaPIDCuts(v0) && TMath::Abs(kfpAntiLambda.GetE())>TMath::Abs(kfpAntiLambda.GetPz()) && kfpAntiLambda.GetNDF()>0 && kfpAntiLambda.GetChi2()>0 && CheckKFParticleCov(kfpAntiLambda) && err_massAntiLambda>0) {
      if (fabs(kfpAntiLambda.GetRapidity())<0.8) {
        fCounterRec_Cuts_AntiLambda->Fill(7);
        fHistChi2ndfAntiLambda->Fill(kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF());
        fHistProbAntiLambda->Fill(TMath::Prob(kfpAntiLambda.GetChi2(), kfpAntiLambda.GetNDF()));
        if (TMath::Abs(kfpAntiLambda.GetRapidity())<0.5) fCounterRec_Cuts_AntiLambda->Fill(8);
        if ( fIsMC ) {
          Int_t lab_AntiLam = MatchToMCAntiLambda(trkN, trkP, mcArray);
          if (lab_AntiLam>-1) {
            fCounterRecMC_Cuts_AntiLambda->Fill(7);
            f2DCounterRecMC_CutsVsPt_AntiLambda->Fill(7, kfpAntiLambda.GetPt());
            if (TMath::Abs(kfpAntiLambda.GetRapidity())<0.5) {
              fCounterRecMC_Cuts_AntiLambda->Fill(8);
              f2DCounterRecMC_CutsVsPt_AntiLambda->Fill(8, kfpAntiLambda.GetPt());
            }
          }
        }
      }
    }


//    fHistMassLambda_LDeltaL->Fill(kfpLambda.GetMass());
    f2DHistArmenterosPodolanski_candidate->Fill(alpha_FirstDaugPos, v0->PtArmV0());

    //************************** calculate l/Δl for Lambda *************************************
    Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
    Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
    Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
    Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
    Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
    if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
    dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
    Double_t nErr_l_Lambda = l_Lambda/dl_Lambda;
    //***************************************************************************************
    //************************** calculate l/Δl for Anti-Lambda *************************************
    Double_t dx_AntiLambda = PV.GetX()-kfpAntiLambda.GetX();
    Double_t dy_AntiLambda = PV.GetY()-kfpAntiLambda.GetY();
    Double_t dz_AntiLambda = PV.GetZ()-kfpAntiLambda.GetZ();
    Double_t l_AntiLambda = TMath::Sqrt(dx_AntiLambda*dx_AntiLambda + dy_AntiLambda*dy_AntiLambda + dz_AntiLambda*dz_AntiLambda);
    Double_t dl_AntiLambda = (PV.GetCovariance(0)+kfpAntiLambda.GetCovariance(0))*dx_AntiLambda*dx_AntiLambda + (PV.GetCovariance(2)+kfpAntiLambda.GetCovariance(2))*dy_AntiLambda*dy_AntiLambda + (PV.GetCovariance(5)+kfpAntiLambda.GetCovariance(5))*dz_AntiLambda*dz_AntiLambda + 2*( (PV.GetCovariance(1)+kfpAntiLambda.GetCovariance(1))*dx_AntiLambda*dy_AntiLambda + (PV.GetCovariance(3)+kfpAntiLambda.GetCovariance(3))*dx_AntiLambda*dz_AntiLambda + (PV.GetCovariance(4)+kfpAntiLambda.GetCovariance(4))*dy_AntiLambda*dz_AntiLambda );
    if ( fabs(l_AntiLambda)<1.e-8f ) l_AntiLambda = 1.e-8f;
    dl_AntiLambda = dl_AntiLambda<0. ? 1.e8f : sqrt(dl_AntiLambda)/l_AntiLambda;
    Double_t nErr_l_AntiLambda = l_AntiLambda/dl_AntiLambda;
    //***************************************************************************************


    //======================================== start to reconstruct XicZero ============================
    if ( fAnaCuts->LambdaPIDCuts(v0) && TMath::Abs(kfpLambda.GetE())>TMath::Abs(kfpLambda.GetPz()) && kfpLambda.GetNDF()>0 && kfpLambda.GetChi2()>0 && CheckKFParticleCov(kfpLambda) && err_massLambda>0 && (kfpLambda.GetChi2()/kfpLambda.GetNDF())<fAnaCuts->GetKFPLam_Chi2geoMax() && nErr_l_Lambda>fAnaCuts->GetKFPLam_lDeltalMin() ) {

      fHistProbLambda_chi2cut->Fill(TMath::Prob(kfpLambda.GetChi2(), kfpLambda.GetNDF()));
      fHistMassLambda_woArmenterosPodolanskiCut->Fill(kfpLambda.GetMass());
      fHistMassLambda_woMassCut->Fill(kfpLambda.GetMass());

      if (TMath::Abs(kfpLambda.GetRapidity())<0.5) {
        fHistLDeltaLRec_Lambda->Fill(nErr_l_Lambda);
        if (fIsMC) {
          Int_t lab_Lam = MatchToMCLambda(trkP, trkN, mcArray);
          if (lab_Lam>-1) fHistLDeltaLRecMC_Lambda->Fill(nErr_l_Lambda);
          Int_t lab_LamFromXi = MatchToMCLambdaFromXi(trkP, trkN, mcArray);

          if (lab_LamFromXi>-1) fHistLDeltaLRecMC_LambdaFromXi->Fill(nErr_l_Lambda);
        }
      }

//      Double_t y = v0->PtArmV0();
//      Double_t x = alpha_FirstDaugPos;
//      Double_t ycut_L = 1.3 * (-1.*x*x + 1.38*x - 0.420);
//      Double_t ycut_H = 1.3 * (-1.*x*x + 1.38*x - 0.380);
//      if ( y<ycut_L || y>ycut_H ) continue;


      fHistMassLambda_BeforeSecSel->Fill(kfpLambda.GetMass());

          // select secondary Lambda
          /*
          KFParticle kfpLambda_sec = kfpLambda;
//        cout << "PV:" << PV << endl;
//        cout << "Lambda:" << kfpLambda << endl;
          kfpLambda_sec.SetProductionVertex(PV);
          */

//      cout << "DM2=" << 4.* (kfpLambda.GetCovariance(9)*kfpLambda.GetPx()*kfpLambda.GetPx() + kfpLambda.GetCovariance(14)*kfpLambda.GetPy()*kfpLambda.GetPy() + kfpLambda.GetCovariance(20)*kfpLambda.GetPz()*kfpLambda.GetPz() + kfpLambda.GetCovariance(27)*kfpLambda.GetE()*kfpLambda.GetE() + float(2.f) * ( (kfpLambda.GetCovariance(13)*kfpLambda.GetPy() + kfpLambda.GetCovariance(18)*kfpLambda.GetPz() - kfpLambda.GetCovariance(24)*kfpLambda.GetE() )*kfpLambda.GetPx() + ( kfpLambda.GetCovariance(19)*kfpLambda.GetPz() - kfpLambda.GetCovariance(25)*kfpLambda.GetE() )*kfpLambda.GetPy() - kfpLambda.GetCovariance(26)*kfpLambda.GetPz()*kfpLambda.GetE() ) ) << endl;

//          if ( kfpLambda_sec.GetChi2()/kfpLambda_sec.GetNDF()>=3 ) {

            /*
            fCounterRec_Cuts_Lambda->Fill(13);
            f2DCounterRecMC_CutsVsPt_Lambda->Fill(kfpLambda.GetPt(), 13);
            if ( fIsMC ) {
              Int_t lab_Lam     = MatchToMCLambda(trkP, trkN, mcArray);
              if (lab_Lam>-1) fCounterRecMC_Cuts_Lambda->Fill(13);
            }
            */

        f2DHistLambdaXY->Fill(kfpLambda.GetX(), kfpLambda.GetY());
        fHistMassLambda->Fill(kfpLambda.GetMass());
        fHistMassLambdaTot->Fill(kfpLambda.GetMass());
        f2DHistMassPtLambda->Fill(kfpLambda.GetPt(), kfpLambda.GetMass());
        f2DHistChi2vsNDF_Lambda->Fill(kfpLambda.GetNDF(), kfpLambda.GetChi2());
//        fHistProbLambdaTot->Fill(TMath::Prob(kfpLambda.GetChi2(), kfpLambda.GetNDF()));
        Float_t ptLambda_Rec, err_ptLambda;
        kfpLambda.GetPt(ptLambda_Rec, err_ptLambda);
        fHistMassLambda_V0->Fill(v0->MassLambda());
        fHistMassLambdaTot_V0->Fill(v0->MassLambda());
        Double_t Lambda_KF_V0 = kfpLambda.GetMass() - v0->MassLambda();
        fHistMassLambda_KF_V0->Fill(Lambda_KF_V0);
        fHistPtLambda->Fill(kfpLambda.GetPt());
        fHistPtLambdaTot->Fill(kfpLambda.GetPt());
        
//        cout << "Proton" << kfpProton << endl;
//        cout << "Pion" << kfpPionMinus << endl;
//        cout << "Lambda:" << kfpLambda << endl;

/*
          if ( fIsMC && trkP->GetLabel()>0 && trkN->GetLabel()>0 ) {
            AliAODMCParticle* mcTrkP = static_cast<AliAODMCParticle*>(mcArray->At(trkP->GetLabel()));  
            AliAODMCParticle* mcTrkN = static_cast<AliAODMCParticle*>(mcArray->At(trkN->GetLabel()));  
            if ( mcTrkP->GetPdgCode()==2212 && mcTrkN->GetPdgCode()==-211 ) {
              Int_t pdgDgv0[2]={2212, 211};
              Int_t labelLambda = v0->MatchToMC(3122, mcArray, 2, pdgDgv0);
              if ( labelLambda>-1 ) fCounterRecMC_Cuts_Lambda->Fill(13);
            }
          }
*/

        //======================================== reconstruct Xi- ============================
        if ( TMath::Abs(massLambda_Rec-massLambda)<=(fAnaCuts->GetProdMassTolLambda()) ) {
          f2DHistArmenterosPodolanski_Lam->Fill(alpha_FirstDaugPos, v0->PtArmV0());
          KFParticle kfpLambda_m = kfpLambda;
          kfpLambda_m.SetNonlinearMassConstraint(massLambda);
//          if (kfpLambda_m.GetE()<kfpLambda_m.GetP()) cout << "E<P!!!!!!!!!!!!!!!!" << endl;
//          kfpLambda_m.SetMassConstraint(massLambda);
//          Float_t massLam, errormassLam;
//          kfpLambda.GetMass(massLam, errormassLam);
//          if (kfpLambda.GetMass()<0.1) {cout << "massLam==" << massLam << endl; cout << "errormassLam==" << errormassLam << endl;}
          fHistMassLambda_M->Fill(kfpLambda_m.GetMass());
          if ( CheckKFParticleCov(kfpLambda_m) && TMath::Abs(kfpLambda.GetE())>TMath::Abs(kfpLambda.GetPz()) ) {

            for (Int_t itrkPion2=0; itrkPion2<flag_trkN; itrkPion2++) { // Loop for bachelor in cascade pion-
              if ( trackN[itrkPion2]->GetID()==trkN->GetID() ) continue;
              fCounterRec_Cuts_XiMinus->Fill(1);
              if ( fIsMC ) {
                Int_t lab_XiMinus = MatchToMCXiMinus(trkP, trkN, trackN[itrkPion2], mcArray);
                if (lab_XiMinus>-1) {
                  fCounterRecMC_Cuts_XiMinus->Fill(1);
                }
              }

              if ( !fAnaCuts->PassedTrackQualityCuts_SecondaryPion(trackN[itrkPion2]) ) continue;

              KFParticle kfpPion2 = CreateKFParticleFromAODtrack(trackN[itrkPion2], -211);

              // select secondary pion tracks
              /*
              KFParticle kfpPion2_sec = kfpPion2;
              kfpPion2_sec.SetProductionVertex(PV);
              if ( kfpPion2_sec.GetChi2()<=18.4 ) continue;
              fCounterRec_Cuts_XiMinus->Fill(2);
              if ( fIsMC ) {
                Int_t lab_XiMinus = MatchToMCXiMinus(trkP, trkN, trackN[itrkPion2], mcArray);
                if (lab_XiMinus>-1) {
                  fCounterRecMC_Cuts_XiMinus->Fill(2);
                }
              }
              */

              KFParticle kfpXiMinus;
//              kfpXiMinus.SetPDG(3312);
              const KFParticle *vXiDs[2] = {&kfpPion2, &kfpLambda_m};
              kfpXiMinus.Construct(vXiDs, NDaughters);

              // check rapidity of Xi-
              if ( TMath::Abs(kfpXiMinus.GetE())<=TMath::Abs(kfpXiMinus.GetPz()) ) continue;

              // err_massXi > 0
              Float_t massXiMinus_Rec, err_massXiMinus;
              kfpXiMinus.GetMass(massXiMinus_Rec, err_massXiMinus);
              if ( err_massXiMinus<=0 ) continue;

              if ( fabs(kfpXiMinus.GetRapidity())<0.8 ) {
                fCounterRec_Cuts_XiMinus->Fill(2);
                if ( fIsMC ) {
                  Int_t lab_XiMinus = MatchToMCXiMinus(trkP, trkN, trackN[itrkPion2], mcArray);
                  if (lab_XiMinus>-1) {
                    fCounterRecMC_Cuts_XiMinus->Fill(2);
                    f2DCounterRecMC_CutsVsPt_XiMinus->Fill(2, kfpXiMinus.GetPt());
                  }
                }
              }

              // chi2>0 && NDF>0
              if ( kfpXiMinus.GetNDF()<=0 || kfpXiMinus.GetChi2()<=0 ) continue;

              // Prefilter
              if ( kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF() >= fAnaCuts->GetKFPXi_Chi2geoMax() ) continue;

              // check covariance matrix
              if ( !CheckKFParticleCov(kfpXiMinus) ) continue;

              if ( fabs(kfpXiMinus.GetRapidity())<0.8 ) {
                fCounterRec_Cuts_XiMinus->Fill(3);
                if ( fIsMC ) {
                  Int_t lab_XiMinus = MatchToMCXiMinus(trkP, trkN, trackN[itrkPion2], mcArray);
                  if (lab_XiMinus>-1) {
                    fCounterRecMC_Cuts_XiMinus->Fill(3);
                    f2DCounterRecMC_CutsVsPt_XiMinus->Fill(3, kfpXiMinus.GetPt());
                  }
                }
              }
              fHistChi2ndfXiMinus->Fill(kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF());
              fHistProbXiMinus->Fill(TMath::Prob(kfpXiMinus.GetChi2(), kfpXiMinus.GetNDF()));
              /*
              // 0<chi2/NDF<cut
              if ( kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF()<=0 || (kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF())>=(fAnaCuts->GetKFPXiMinusChi2Cut()) ) continue;
              fCounterRec_Cuts_XiMinus->Fill(5);
              if ( fIsMC ) {
                Int_t lab_XiMinus = MatchToMCXiMinus(trkP, trkN, trackN[itrkPion2], mcArray);
                if (lab_XiMinus>-1) {
                  fCounterRecMC_Cuts_XiMinus->Fill(5);
                  f2DCounterRecMC_CutsVsPt_XiMinus->Fill(5, kfpXiMinus.GetPt());
                }
              }
              fHistProbXiMinus_chi2cut->Fill(TMath::Prob(kfpXiMinus.GetChi2(), kfpXiMinus.GetNDF()));
              */
              // set production vertext to PV
//              KFParticle kfpXiMinus_pv = kfpXiMinus;
//              kfpXiMinus_pv.SetProductionVertex(PV);
              /*
              if ( kfpXiMinus_pv.GetChi2()/kfpXiMinus_pv.GetNDF()>=3 ) continue;
              fCounterRec_Cuts_XiMinus->Fill(6);
              if ( fIsMC ) {
                Int_t lab_XiMinus = MatchToMCXiMinus(trkP, trkN, trackN[itrkPion2], mcArray);
                if (lab_XiMinus>-1) {
                  fCounterRecMC_Cuts_XiMinus->Fill(6);
                  f2DCounterRecMC_CutsVsPt_XiMinus->Fill(6, kfpXiMinus_pv.GetPt());
                }
              }
              */
              // store mass of Xi-
              fHistMassXiMinus->Fill(kfpXiMinus.GetMass());
              fHistMassXiTot->Fill(kfpXiMinus.GetMass());
              f2DHistMassPtXiMinus->Fill(kfpXiMinus.GetPt(), kfpXiMinus.GetMass());
              
              if ( fabs(kfpXiMinus.GetRapidity())<0.8 ) {
                fCounterRec_Cuts_XiMinus->Fill(4);
                if (TMath::Abs(kfpXiMinus.GetRapidity())<0.5) fCounterRec_Cuts_XiMinus->Fill(5);
                if ( fIsMC ) {
                  Int_t lab_XiMinus = MatchToMCXiMinus(trkP, trkN, trackN[itrkPion2], mcArray);
                  if (lab_XiMinus>-1) {
                    fCounterRecMC_Cuts_XiMinus->Fill(4);
                    f2DCounterRecMC_CutsVsPt_XiMinus->Fill(4, kfpXiMinus.GetPt());
                    if (TMath::Abs(kfpXiMinus.GetRapidity())<0.5) {
                      fCounterRecMC_Cuts_XiMinus->Fill(5);
                      f2DCounterRecMC_CutsVsPt_XiMinus->Fill(5, kfpXiMinus.GetPt());
                    }
                  }
                }
              }
              
              /*
              // calculate proper decay length of Xi- and lambda
              if ( kfpLambda.GetP()>0 && kfpXiMinus_pv.GetP()>0 ) {
                Double_t properDecayL_XiMinus = kfpXiMinus_pv.GetDecayLength() * kfpXiMinus_pv.GetMass() / kfpXiMinus_pv.GetP();
                fHistDecayLXiMinus->Fill(properDecayL_XiMinus);
                f2DHistXiMinusXY_DV->Fill(kfpXiMinus.GetX(), kfpXiMinus.GetY());
                f2DHistXiMinusXY_PV->Fill(kfpXiMinus_pv.GetX(), kfpXiMinus_pv.GetY());
                // calculate proper decay length of lambda
                kfpLambda.SetProductionVertex(kfpXiMinus);
                Double_t properDecayL_Lambda = kfpLambda.GetDecayLength() * kfpLambda.GetMass() / kfpLambda.GetP();
                fHistDecayLLambda->Fill(properDecayL_Lambda);
              }
              */



              //======================================== reconstruct XicZero ============================
              if ( TMath::Abs(massXiMinus_Rec-massXi)<=(fAnaCuts->GetProdMassTolXi()) ) {
              // set nonlinear mass constraint
              KFParticle kfpXiMinus_m = kfpXiMinus;
              kfpXiMinus_m.SetNonlinearMassConstraint(massXi);
//              kfpXiMinus_m.SetMassConstraint(massXi);
              fHistMassXiMinus_M->Fill(kfpXiMinus_m.GetMass());
              if ( CheckKFParticleCov(kfpXiMinus_m) && TMath::Abs(kfpXiMinus.GetE())>TMath::Abs(kfpXiMinus.GetPz()) ) {

              for (Int_t itrkBP=0; itrkBP<flag_trkP; itrkBP++) { // Loop for first bachelor pion+
                if ( trackP[itrkBP]->GetID()==trkP->GetID() ) continue;

//            if ( !fAnaCuts->SingleTrkCuts(trackP[itrkBP]) ) continue;
                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP]) ) continue;

            // DCA of Cascade-bachelor to PV (cm)
//            Double_t d0z0bach[2],covd0z0bach[3];
//            if (!trackN[itrkPion2]->PropagateToDCA(fpVtx,fBzkG,kVeryBig,d0z0bach,covd0z0bach)) continue;
//            if ( d0z0bach[0] <= fAnaCuts->GetProdDcaBachToPrimVertexMin() ) continue;

                KFParticle kfpBP    = CreateKFParticleFromAODtrack(trackP[itrkBP], 211);

                // select primary pion tracks
                /*
                KFParticle kfpBP_prm = kfpBP;
                kfpBP_prm.SetProductionVertex(PV);
                if ( kfpBP_prm.GetChi2()>18.4 ) continue;
                */

//            if (fIsMC) {
//              Int_t labelXiMinus = MatchToMCXiMinus(trkP, trkN, trackN[itrkPion2], mcArray);
//              Int_t labelXiPlus = MatchToMCXiPlus(trkN, trkP, trackP[itrkBP], mcArray);
//              if (labelXiMinus>-1) {
//                fHistMassXiMinus_Match->Fill(kfpXiMinus.GetMass());
//              }
//            }

                // reconstruct XicZero
                KFParticle kfpXic0;
                const KFParticle *vXicZeroDs[2] = {&kfpBP, &kfpXiMinus_m};
                kfpXic0.Construct(vXicZeroDs, NDaughters);
                fHistProbXicZero->Fill(TMath::Prob(kfpXic0.GetChi2(), kfpXic0.GetNDF()));



                // chi2>0 && NDF>0
                if ( kfpXic0.GetNDF()<=0 || kfpXic0.GetChi2()<=0 ) continue;

                // Prefilter
                if ( kfpXic0.GetChi2()/kfpXic0.GetNDF() >= fAnaCuts->GetKFPXic0_Chi2geoMax() ) continue;
                if ( kfpXic0.GetPt() < fAnaCuts->GetPtMinXic0() ) continue;

                // check rapidity of XicZero
                if ( TMath::Abs(kfpXic0.GetE())<=TMath::Abs(kfpXic0.GetPz()) ) continue;

                // check covariance matrix
                if ( !CheckKFParticleCov(kfpXic0) ) continue;

                // err_massXicZero > 0
                Float_t massXicZero_Rec, err_massXicZero;
                kfpXic0.GetMass(massXicZero_Rec, err_massXicZero);
                if ( err_massXicZero<=0 ) continue;

                /*
                // 0<chi2/NDF<cut
                if ( kfpXic0.GetChi2()/kfpXic0.GetNDF()<=0 || (kfpXic0.GetChi2()/kfpXic0.GetNDF())>=(fAnaCuts->GetKFPXicZeroChi2Cut()) ) continue;
                fHistProbXicZero_chi2cut->Fill(TMath::Prob(kfpXic0.GetChi2(), kfpXic0.GetNDF()));
                */

                // set production vertext to PV
//                KFParticle kfpXic0_pv = kfpXic0;
//                kfpXic0_pv.SetProductionVertex(PV);
                /*
                if ( kfpXic0_pv.GetChi2()/kfpXic0_pv.GetNDF()>=3 ) continue;
                */
                /*
                // store mass of XicZero
                fHistMassXicZero->Fill(kfpXic0_pv.GetMass());
                fHistMassXicZeroTot->Fill(kfpXic0_pv.GetMass());
                f2DHistMassPtXicZero->Fill(kfpXic0_pv.GetPt(), kfpXic0_pv.GetMass());
                f2DHistMassPtXicZeroTot->Fill(kfpXic0_pv.GetPt(), kfpXic0_pv.GetMass());


                TVector3 momXiMinus(kfpXiMinus.GetPx(), kfpXiMinus.GetPy(), kfpXiMinus.GetPz());
//                TVector3 momBP(kfpBP_prm.GetPx(), kfpBP_prm.GetPy(), kfpBP_prm.GetPz());
                TVector3 momBP(kfpBP.GetPx(), kfpBP.GetPy(), kfpBP.GetPz());
                TVector3 momXicZero(kfpXic0_pv.GetPx(), kfpXic0_pv.GetPy(), kfpXic0_pv.GetPz());
                Float_t QtDiffPionXiMinus = momBP.Perp(momXicZero) - momXiMinus.Perp(momXicZero);
                fHistQtDiffPionXiMinus->Fill(QtDiffPionXiMinus);
                */

//            fHistMassXicZero_woQtCut->Fill(kfpXic0_pv.GetMass());
//            f2DHistMassPtXicZero_woQtCut->Fill(kfpXic0_pv.GetPt(), kfpXic0_pv.GetMass());
//            if (TMath::Abs(QtDiffPionXiMinus)>1.e-5) continue;

                if (fWriteXic0Tree) {
                  Int_t lab_Xic0 = -9999;
                  if (fIsMC) {
                    lab_Xic0 = MatchToMCXic0(trkP, trkN, trackN[itrkPion2], trackP[itrkBP], mcArray);
                  }
                  FillTreeRecXic0FromV0(kfpXic0, trackP[itrkBP], kfpBP, kfpXiMinus, kfpXiMinus_m, trackN[itrkPion2], v0, kfpK0Short, kfpLambda, kfpLambda_m, trkP, trkN, PV, mcArray, lab_Xic0);
                }

                kfpXic0.Clear();
                kfpBP.Clear();
              } // loop for XicZero
              }
              kfpXiMinus_m.Clear();
              } // Xi- mass cut
            kfpXiMinus.Clear();
            kfpPion2.Clear();
            } // loop for Xi-
          } // check cov. of lambda && |E|>|Pz| after nonlinear mass constrain
        kfpLambda_m.Clear();
        } // lambda mass cut
//          } // select secondary lambda
    } // separate lambda & anti-lambda and check cov. of lambda


    //======================================== start to reconstruct Anti-XicZero ============================
    if ( fAnaCuts->AntiLambdaPIDCuts(v0) && TMath::Abs(kfpAntiLambda.GetE())>TMath::Abs(kfpAntiLambda.GetPz()) && kfpAntiLambda.GetNDF()>0 && kfpAntiLambda.GetChi2()>0 && CheckKFParticleCov(kfpAntiLambda) && err_massAntiLambda>0 && (kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF())<fAnaCuts->GetKFPLam_Chi2geoMax() && nErr_l_AntiLambda>fAnaCuts->GetKFPLam_lDeltalMin() ) {
      fHistProbAntiLambda_chi2cut->Fill(TMath::Prob(kfpAntiLambda.GetChi2(), kfpAntiLambda.GetNDF()));
      fHistMassAntiLambda_woMassCut->Fill(kfpAntiLambda.GetMass());
      fHistMassAntiLambda_BeforeSecSel->Fill(kfpAntiLambda.GetMass());

      if (TMath::Abs(kfpAntiLambda.GetRapidity())<0.5) {
        fHistLDeltaLRec_AntiLambda->Fill(nErr_l_AntiLambda);
        if (fIsMC) {
          Int_t lab_AntiLam = MatchToMCAntiLambda(trkN, trkP, mcArray);
          if (lab_AntiLam>-1) fHistLDeltaLRecMC_AntiLambda->Fill(nErr_l_AntiLambda);
          Int_t lab_AntiLamFromXi = MatchToMCAntiLambdaFromXi(trkN, trkP, mcArray);
          if (lab_AntiLamFromXi>-1) fHistLDeltaLRecMC_AntiLambdaFromXi->Fill(nErr_l_AntiLambda);
        }
      }


          // select secondary Anti-Lambda
          /*
          KFParticle kfpAntiLambda_sec = kfpAntiLambda;
          kfpAntiLambda_sec.SetProductionVertex(PV);
          */

//          if ( kfpAntiLambda_sec.GetChi2()/kfpAntiLambda_sec.GetNDF()>=3 ) {

        /*
        fCounterRec_Cuts_AntiLambda->Fill(13);
        f2DCounterRecMC_CutsVsPt_AntiLambda->Fill(kfpAntiLambda.GetPt(), 13);
        if ( fIsMC ) {
          Int_t lab_AntiLam = MatchToMCAntiLambda(trkN, trkP, mcArray);
          if (lab_AntiLam>-1) fCounterRecMC_Cuts_AntiLambda->Fill(13);
        }
        */

        fHistMassAntiLambda->Fill(kfpAntiLambda.GetMass());
        fHistMassLambdaTot->Fill(kfpAntiLambda.GetMass());
        f2DHistMassPtAntiLambda->Fill(kfpAntiLambda.GetPt(), kfpAntiLambda.GetMass());
        f2DHistChi2vsNDF_AntiLambda->Fill(kfpAntiLambda.GetNDF(), kfpAntiLambda.GetChi2());
//        fHistChi2ndfLambdaTot->Fill(kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF());
//        fHistProbLambdaTot->Fill(TMath::Prob(kfpAntiLambda.GetChi2(), kfpAntiLambda.GetNDF()));
        Float_t ptAntiLambda_Rec, err_ptAntiLambda;
        kfpAntiLambda.GetPt(ptAntiLambda_Rec, err_ptAntiLambda);
        fHistMassAntiLambda_V0->Fill(v0->MassAntiLambda());
        fHistMassLambdaTot_V0->Fill(v0->MassAntiLambda());
        Double_t AntiLambda_KF_V0 = kfpAntiLambda.GetMass() - v0->MassAntiLambda();
        fHistMassAntiLambda_KF_V0->Fill(AntiLambda_KF_V0);
        fHistPtAntiLambda->Fill(kfpAntiLambda.GetPt());
        fHistPtLambdaTot->Fill(kfpAntiLambda.GetPt());


        /*
        if (fIsMC && trkP->GetLabel()>0 && trkN->GetLabel()>0) {
          AliAODMCParticle* mcTrkP = static_cast<AliAODMCParticle*>(mcArray->At(trkP->GetLabel()));
          AliAODMCParticle* mcTrkN = static_cast<AliAODMCParticle*>(mcArray->At(trkN->GetLabel()));
          Int_t pdgDgv0[2]={2212, 211};
          Int_t labelLambda = v0->MatchToMC(3122, mcArray, 2, pdgDgv0);                                 

          if ( labelLambda>-1 && mcTrkP->GetPdgCode()==211 && mcTrkN->GetPdgCode()==-2212 ) {
            fCounterRec_Cuts_AntiLambda->Fill(14);
          }
        }
        */

        //======================================== reconstruct Xi+ ============================
        if ( TMath::Abs(massAntiLambda_Rec-massLambda)<=(fAnaCuts->GetProdMassTolLambda()) ) {
          f2DHistArmenterosPodolanski_AntiLam->Fill(alpha_FirstDaugPos, v0->PtArmV0());
          KFParticle kfpAntiLambda_m = kfpAntiLambda;
          kfpAntiLambda_m.SetNonlinearMassConstraint(massLambda);
//          kfpAntiLambda_m.SetMassConstraint(massLambda);
          fHistMassAntiLambda_M->Fill(kfpAntiLambda.GetMass());
          if ( CheckKFParticleCov(kfpAntiLambda_m) && TMath::Abs(kfpAntiLambda.GetE())>TMath::Abs(kfpAntiLambda.GetPz()) ) {

            for (Int_t itrkPion2=0; itrkPion2<flag_trkP; itrkPion2++) { // Loop for bachelor in cascade pion+
              if ( trackP[itrkPion2]->GetID()==trkP->GetID() ) continue;
              fCounterRec_Cuts_XiPlus->Fill(1);
              if ( fIsMC ) {
                Int_t lab_XiPlus = MatchToMCXiPlus(trkN, trkP, trackP[itrkPion2], mcArray);
                if (lab_XiPlus>-1) {
                  fCounterRecMC_Cuts_XiPlus->Fill(1);
                }
              }

              if ( !fAnaCuts->PassedTrackQualityCuts_SecondaryPion(trackP[itrkPion2]) ) continue;

              KFParticle kfpPion2 = CreateKFParticleFromAODtrack(trackP[itrkPion2], 211);

              // select secondary pion tracks
              /*
              KFParticle kfpPion2_sec = kfpPion2;
              kfpPion2_sec.SetProductionVertex(PV);
              if ( kfpPion2_sec.GetChi2()<=18.4 ) continue;
              fCounterRec_Cuts_XiPlus->Fill(2);
              if ( fIsMC ) {
                Int_t lab_XiPlus = MatchToMCXiPlus(trkN, trkP, trackP[itrkPion2], mcArray);
                if (lab_XiPlus>-1) {
                  fCounterRecMC_Cuts_XiPlus->Fill(2);
                }
              }
              */

              KFParticle kfpXiPlus;
//              kfpXiPlus.SetPDG(-3312);
              const KFParticle *vXiDs[2] = {&kfpPion2, &kfpAntiLambda_m};
              kfpXiPlus.Construct(vXiDs, NDaughters);

              // check rapidity of Xi+
              if ( TMath::Abs(kfpXiPlus.GetE())<=TMath::Abs(kfpXiPlus.GetPz()) ) continue;

              if ( fabs(kfpXiPlus.GetRapidity())<0.8 ) {
                fCounterRec_Cuts_XiPlus->Fill(2);
                if ( fIsMC ) {
                  Int_t lab_XiPlus = MatchToMCXiPlus(trkN, trkP, trackP[itrkPion2], mcArray);
                  if (lab_XiPlus>-1) {
                    fCounterRecMC_Cuts_XiPlus->Fill(2);
                    f2DCounterRecMC_CutsVsPt_XiPlus->Fill(2, kfpXiPlus.GetPt());
                  }
                }
              }

              // chi2>0 && NDF>0
              if ( kfpXiPlus.GetNDF()<=0 || kfpXiPlus.GetChi2()<=0 ) continue;

              // Prefilter
              if ( kfpXiPlus.GetChi2()/kfpXiPlus.GetNDF() >= fAnaCuts->GetKFPXi_Chi2geoMax() ) continue;

              // check covariance matrix
              if ( !CheckKFParticleCov(kfpXiPlus) ) continue;

              // cut on mass window
              Float_t massXiPlus_Rec, err_massXiPlus;
              kfpXiPlus.GetMass(massXiPlus_Rec, err_massXiPlus);
              if ( err_massXiPlus<=0 ) continue;

              if ( fabs(kfpXiPlus.GetRapidity())<0.8 ) {
                fCounterRec_Cuts_XiPlus->Fill(3);
                if ( fIsMC ) {
                  Int_t lab_XiPlus = MatchToMCXiPlus(trkN, trkP, trackP[itrkPion2], mcArray);
                  if (lab_XiPlus>-1) {
                    fCounterRecMC_Cuts_XiPlus->Fill(3);
                    f2DCounterRecMC_CutsVsPt_XiPlus->Fill(3, kfpXiPlus.GetPt());
                  }
                }
              }
              fHistChi2ndfXiPlus->Fill(kfpXiPlus.GetChi2()/kfpXiPlus.GetNDF());
              fHistProbXiPlus->Fill(TMath::Prob(kfpXiPlus.GetChi2(), kfpXiPlus.GetNDF()));
              /*
              // 0<chi2/NDF<cut
              if ( kfpXiPlus.GetChi2()/kfpXiPlus.GetNDF()<=0 || (kfpXiPlus.GetChi2()/kfpXiPlus.GetNDF())>=(fAnaCuts->GetKFPXiPlusChi2Cut()) ) continue;
              fCounterRec_Cuts_XiPlus->Fill(5);
              if ( fIsMC ) {
                Int_t lab_XiPlus = MatchToMCXiPlus(trkN, trkP, trackP[itrkPion2], mcArray);
                if (lab_XiPlus>-1) {
                  fCounterRecMC_Cuts_XiPlus->Fill(5);
                  f2DCounterRecMC_CutsVsPt_XiPlus->Fill(5, kfpXiPlus.GetPt());
                }
              }
              fHistProbXiPlus_chi2cut->Fill(TMath::Prob(kfpXiPlus.GetChi2(), kfpXiPlus.GetNDF()));
              */
              // set production vertext to PV
//              KFParticle kfpXiPlus_pv = kfpXiPlus;
//              kfpXiPlus_pv.SetProductionVertex(PV);
              /*
              if ( kfpXiPlus_pv.GetChi2()/kfpXiPlus_pv.GetNDF()>=3 ) continue;
              fCounterRec_Cuts_XiPlus->Fill(6);
              if ( fIsMC ) {
                Int_t lab_XiPlus = MatchToMCXiPlus(trkN, trkP, trackP[itrkPion2], mcArray);
                if (lab_XiPlus>-1) {
                  fCounterRecMC_Cuts_XiPlus->Fill(6);
                  f2DCounterRecMC_CutsVsPt_XiPlus->Fill(6, kfpXiPlus_pv.GetPt());
                }
              }
              */
              // store mass of Xi+
              fHistMassXiPlus->Fill(kfpXiPlus.GetMass());
              fHistMassXiTot->Fill(kfpXiPlus.GetMass());
              f2DHistMassPtXiPlus->Fill(kfpXiPlus.GetPt(), kfpXiPlus.GetMass());

              if ( fabs(kfpXiPlus.GetRapidity())<0.8 ) {
                fCounterRec_Cuts_XiPlus->Fill(4);
                if (TMath::Abs(kfpXiPlus.GetRapidity())<0.5) fCounterRec_Cuts_XiPlus->Fill(5);
                if ( fIsMC ) {
                  Int_t lab_XiPlus = MatchToMCXiPlus(trkN, trkP, trackP[itrkPion2], mcArray);
                  if (lab_XiPlus>-1) {
                    fCounterRecMC_Cuts_XiPlus->Fill(4);
                    f2DCounterRecMC_CutsVsPt_XiPlus->Fill(4, kfpXiPlus.GetPt());
                    if (TMath::Abs(kfpXiPlus.GetRapidity())<0.5) {
                      fCounterRecMC_Cuts_XiPlus->Fill(5);
                      f2DCounterRecMC_CutsVsPt_XiPlus->Fill(5, kfpXiPlus.GetPt());
                    }
                  }
                }
              }

              /*
              // calculate proper decay length of Xi+ and Anti-lambda
              if ( kfpAntiLambda.GetP()>0 || kfpXiPlus_pv.GetP()>0 ) {
                Double_t properDecayL_XiPlus = kfpXiPlus_pv.GetDecayLength() * kfpXiPlus_pv.GetMass() / kfpXiPlus_pv.GetP();
                fHistDecayLXiPlus->Fill(properDecayL_XiPlus);
                f2DHistXiPlusXY_DV->Fill(kfpXiPlus.GetX(), kfpXiPlus.GetY());
                f2DHistXiPlusXY_PV->Fill(kfpXiPlus_pv.GetX(), kfpXiPlus_pv.GetY());
                // calculate proper decay length of Anti-lambda
                kfpAntiLambda.SetProductionVertex(kfpXiPlus);
                Double_t properDecayL_AntiLambda = kfpAntiLambda.GetDecayLength() * kfpAntiLambda.GetMass() / kfpAntiLambda.GetP();
                fHistDecayLAntiLambda->Fill(properDecayL_AntiLambda);
              }
              */
              



              //======================================== reconstruct Anti-XicZero ============================
              if ( TMath::Abs(massXiPlus_Rec-massXi)<=(fAnaCuts->GetProdMassTolXi()) ) {
              // set nonlinear mass constraint
              KFParticle kfpXiPlus_m = kfpXiPlus;
              kfpXiPlus_m.SetNonlinearMassConstraint(massXi);
//              kfpXiPlus_m.SetMassConstraint(massXi);
              fHistMassXiPlus_M->Fill(kfpXiPlus_m.GetMass());
              if ( CheckKFParticleCov(kfpXiPlus_m) && TMath::Abs(kfpXiPlus.GetE())>TMath::Abs(kfpXiPlus.GetPz()) ) {
              for (Int_t itrkBP=0; itrkBP<flag_trkN; itrkBP++) { // Loop for first bachelor pion-
                if ( trackN[itrkBP]->GetID()==trkN->GetID() ) continue;

//            if ( !fAnaCuts->SingleTrkCuts(trackN[itrkBP]) ) continue;
                if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP]) ) continue;

            // DCA of Cascade-bachelor to PV (cm)
//            Double_t d0z0bach[2],covd0z0bach[3];
//            if (!trackP[itrkPion2]->PropagateToDCA(fpVtx,fBzkG,kVeryBig,d0z0bach,covd0z0bach)) continue;
//            if ( d0z0bach[0] <= fAnaCuts->GetProdDcaBachToPrimVertexMin() ) continue;

                KFParticle kfpBP    = CreateKFParticleFromAODtrack(trackN[itrkBP], -211);

                // select primary pion tracks
                /*
                KFParticle kfpBP_prm = kfpBP;
                kfpBP_prm.SetProductionVertex(PV);
                if ( kfpBP_prm.GetChi2()>18.4 ) continue;
                */

//            if (fIsMC) {
//              Int_t labelXiMinus = MatchToMCXiMinus(trkP, trkN, trackN[itrkPion2], mcArray);
//              Int_t labelXiPlus = MatchToMCXiPlus(trkN, trkP, trackP[itrkBP], mcArray);
//              if (labelXiMinus>-1) {
//                fHistMassXiMinus_Match->Fill(kfpXiMinus.GetMass());
//              }
//            }

                // reconstruct Anti-XicZero
                KFParticle kfpAntiXicZero;
                const KFParticle *vXicZeroDs[2] = {&kfpBP, &kfpXiPlus_m};
                kfpAntiXicZero.Construct(vXicZeroDs, NDaughters);
                fHistProbAntiXicZero->Fill(TMath::Prob(kfpAntiXicZero.GetChi2(), kfpAntiXicZero.GetNDF()));
                // check rapidity of Anti-XicZero
                if ( TMath::Abs(kfpAntiXicZero.GetE())<=TMath::Abs(kfpAntiXicZero.GetPz()) ) continue;
                // chi2>0 && NDF>0
                if ( kfpAntiXicZero.GetNDF()<=0 || kfpAntiXicZero.GetChi2()<=0 ) continue;

                // Prefilter
                if ( kfpAntiXicZero.GetChi2()/kfpAntiXicZero.GetNDF() >= fAnaCuts->GetKFPXic0_Chi2geoMax() ) continue;
                if ( kfpAntiXicZero.GetPt() < fAnaCuts->GetPtMinXic0() ) continue;

                // check covariance matrix
                if ( !CheckKFParticleCov(kfpAntiXicZero) ) continue;
                // err_massAntiXicZero > 0
                Float_t massAntiXicZero_Rec, err_massAntiXicZero;
                kfpAntiXicZero.GetMass(massAntiXicZero_Rec, err_massAntiXicZero);
                if ( err_massAntiXicZero<=0 ) continue;

                /*
                // 0<chi2/NDF<cut
                if ( kfpAntiXicZero.GetChi2()/kfpAntiXicZero.GetNDF()<=0 || (kfpAntiXicZero.GetChi2()/kfpAntiXicZero.GetNDF())>=(fAnaCuts->GetKFPAntiXicZeroChi2Cut()) ) continue;
                fHistProbAntiXicZero_chi2cut->Fill(TMath::Prob(kfpAntiXicZero.GetChi2(), kfpAntiXicZero.GetNDF()));
                */
                // set production vertext to PV
//                KFParticle kfpAntiXicZero_pv = kfpAntiXicZero;
//                kfpAntiXicZero_pv.SetProductionVertex(PV);
                /*
                if ( kfpAntiXicZero_pv.GetChi2()/kfpAntiXicZero_pv.GetNDF()>=3 ) continue;
                */
                /*
                // store mass of XicZero
                fHistMassAntiXicZero->Fill(kfpAntiXicZero_pv.GetMass());
                fHistMassXicZeroTot->Fill(kfpAntiXicZero_pv.GetMass());
                f2DHistMassPtAntiXicZero->Fill(kfpAntiXicZero_pv.GetPt(), kfpAntiXicZero_pv.GetMass());
                f2DHistMassPtXicZeroTot->Fill(kfpAntiXicZero_pv.GetPt(), kfpAntiXicZero_pv.GetMass());


                TVector3 momXiPlus(kfpXiPlus.GetPx(), kfpXiPlus.GetPy(), kfpXiPlus.GetPz());
//                TVector3 momBP(kfpBP_prm.GetPx(), kfpBP_prm.GetPy(), kfpBP_prm.GetPz());
                TVector3 momBP(kfpBP.GetPx(), kfpBP.GetPy(), kfpBP.GetPz());
                TVector3 momAntiXicZero(kfpAntiXicZero_pv.GetPx(), kfpAntiXicZero_pv.GetPy(), kfpAntiXicZero_pv.GetPz());
                Float_t QtDiffPionXiPlus = momBP.Perp(momAntiXicZero) - momXiPlus.Perp(momAntiXicZero);
                fHistQtDiffPionXiPlus->Fill(QtDiffPionXiPlus);
                */

//            if (TMath::Abs(QtDiffPionXiPlus)>1.e-5) continue;

                if (fWriteXic0Tree) {
                  Int_t lab_AntiXic0 = -9999.;
                  if (fIsMC) {
                    lab_AntiXic0 = MatchToMCAntiXic0(trkN, trkP, trackP[itrkPion2], trackN[itrkBP], mcArray);
                  }
                  FillTreeRecXic0FromV0(kfpAntiXicZero, trackN[itrkBP], kfpBP, kfpXiPlus, kfpXiPlus_m, trackP[itrkPion2], v0, kfpK0Short, kfpAntiLambda, kfpAntiLambda_m, trkN, trkP, PV, mcArray, lab_AntiXic0);
                }

                kfpAntiXicZero.Clear();
                kfpBP.Clear();
              } // loop for Anti-XicZero
              }
              kfpXiPlus_m.Clear();
              } // Xi+ mass cut
            kfpXiPlus.Clear();
            kfpPion2.Clear();
            } // loop for Xi+
          } // check cov. of Anti-lambda && |E|>|Pz|
          kfpAntiLambda_m.Clear();
        } // Anti-lambda mass cut
//          } // select secondary anti-lambda
    } // separate lambda & anti-lambda and check cov. of Anti-lambda

/*
    // check MC of Lambda
    if (fIsMC) {
      if ( trkP->GetLabel()<0 || trkN->GetLabel()<0 ) continue;
      AliAODMCParticle* mcTrkP = static_cast<AliAODMCParticle*>(mcArray->At(trkP->GetLabel()));  
      AliAODMCParticle* mcTrkN = static_cast<AliAODMCParticle*>(mcArray->At(trkN->GetLabel()));  
      Double_t xTrkP_Rec_MC = xyzP[0] - mcTrkP->Xv();
      Double_t yTrkP_Rec_MC = xyzP[1] - mcTrkP->Yv();
      Double_t zTrkP_Rec_MC = xyzP[2] - mcTrkP->Zv();
      Double_t xTrkP_Rec_MC_Xv = xvyvzvP[0] - mcTrkP->Xv();
      Double_t yTrkP_Rec_MC_Yv = xvyvzvP[1] - mcTrkP->Yv();
      Double_t zTrkP_Rec_MC_Zv = xvyvzvP[2] - mcTrkP->Zv();
//        cout << "x==" << xyzP[0] << endl;
//        cout << "xv==" << xvyvzvP[0] << endl;
//        cout << "xatDCA==" << trk->XAtDCA() << endl;
      fHistXtrkP_Rec_MC->Fill(xTrkP_Rec_MC);
      fHistYtrkP_Rec_MC->Fill(yTrkP_Rec_MC);
      fHistZtrkP_Rec_MC->Fill(zTrkP_Rec_MC);
      fHistXtrkP_Rec_MC_XYZv->Fill(xTrkP_Rec_MC_Xv);
      fHistYtrkP_Rec_MC_XYZv->Fill(yTrkP_Rec_MC_Yv);
      fHistZtrkP_Rec_MC_XYZv->Fill(zTrkP_Rec_MC_Zv);
      Double_t xTrkN_Rec_MC = xyzN[0] - mcTrkN->Xv();
      Double_t yTrkN_Rec_MC = xyzN[1] - mcTrkN->Yv();
      Double_t zTrkN_Rec_MC = xyzN[2] - mcTrkN->Zv();
      Double_t xTrkN_Rec_MC_Xv = xvyvzvN[0] - mcTrkN->Xv();
      Double_t yTrkN_Rec_MC_Yv = xvyvzvN[1] - mcTrkN->Yv();
      Double_t zTrkN_Rec_MC_Zv = xvyvzvN[2] - mcTrkN->Zv();
      fHistXtrkN_Rec_MC->Fill(xTrkN_Rec_MC);
      fHistYtrkN_Rec_MC->Fill(yTrkN_Rec_MC);
      fHistZtrkN_Rec_MC->Fill(zTrkN_Rec_MC);
      fHistXtrkN_Rec_MC_XYZv->Fill(xTrkN_Rec_MC_Xv);
      fHistYtrkN_Rec_MC_XYZv->Fill(yTrkN_Rec_MC_Yv);
      fHistZtrkN_Rec_MC_XYZv->Fill(zTrkN_Rec_MC_Zv);

      Int_t labelLambda = MatchToMCLambda(trkP, trkN, mcArray);
      Int_t labelAntiLambda = MatchToMCAntiLambda(trkN, trkP, mcArray);


      // Lambda
      if (labelLambda>-1) {

        fHistMassLambda_Match->Fill(kfpLambda.GetMass());
        fHistMassLambdaTot_Match->Fill(kfpLambda.GetMass());
        f2DHistChi2vsNDF_Lambda_Match->Fill(kfpLambda.GetNDF(), kfpLambda.GetChi2());
        fHistChi2ndfLambda_Match->Fill(kfpLambda.GetChi2()/kfpLambda.GetNDF());
        fHistProbLambda_Match->Fill(TMath::Prob(kfpLambda.GetChi2(), kfpLambda.GetNDF()));

        AliAODMCParticle* mcLambda = static_cast<AliAODMCParticle*>(mcArray->At(labelLambda));
        AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(mcLambda->GetDaughterFirst()));
        AliAODMCParticle* mcPion = static_cast<AliAODMCParticle*>(mcArray->At(mcLambda->GetDaughterLast()));
//        cout << "PDG Proton ==" << mcProton->GetPdgCode() << endl;
//        cout << "PDG Pion ==" << mcPion->GetPdgCode() << endl;

        //====================== check daughters =============================
        Double_t par_Proton[3], cov_Proton[6]={0.};
        Double_t par_Pion[3], cov_Pion[6]={0.};
        mcProton->XvYvZv(par_Proton);
        mcPion->XvYvZv(par_Pion);
        KFParticle kfpMCproton = CreateKFVertex(par_Proton, cov_Proton);
        KFParticle kfpMCpion = CreateKFVertex(par_Pion, cov_Pion);
        kfpProton.SetProductionVertex(kfpMCproton);
        kfpPionMinus.SetProductionVertex(kfpMCpion);
//        if ( kfpProton.GetNDF()<=0 || kfpProton.GetChi2()<0 ) continue;
//        if ( kfpProton.GetChi2()/kfpProton.GetNDF()<=0 || kfpProton.GetChi2()/kfpProton.GetNDF()>=3 ) continue;

        Double_t xProton_KF_MC = kfpProton.GetX() - mcProton->Xv();
        Double_t yProton_KF_MC = kfpProton.GetY() - mcProton->Yv();
        Double_t zProton_KF_MC = kfpProton.GetZ() - mcProton->Zv();
        Double_t xPion_KF_MC = kfpPionMinus.GetX() - mcPion->Xv();
        Double_t yPion_KF_MC = kfpPionMinus.GetY() - mcPion->Yv();
        Double_t zPion_KF_MC = kfpPionMinus.GetZ() - mcPion->Zv();
        fHistXProton_KF_MC->Fill(xProton_KF_MC);
        fHistYProton_KF_MC->Fill(yProton_KF_MC);
        fHistZProton_KF_MC->Fill(zProton_KF_MC);
        fHistXPion_KF_MC->Fill(xPion_KF_MC);
        fHistYPion_KF_MC->Fill(yPion_KF_MC);
        fHistZPion_KF_MC->Fill(zPion_KF_MC);
        f2DHistXKFMCvsPt_Proton->Fill(kfpProton.GetPt(), xProton_KF_MC);
        f2DHistYKFMCvsPt_Proton->Fill(kfpProton.GetPt(), yProton_KF_MC);
        f2DHistZKFMCvsPt_Proton->Fill(kfpProton.GetPt(), zProton_KF_MC);
        f2DHistXKFMCvsPt_Pion->Fill(kfpPionMinus.GetPt(), xPion_KF_MC);
        f2DHistYKFMCvsPt_Pion->Fill(kfpPionMinus.GetPt(), yPion_KF_MC);
        f2DHistZKFMCvsPt_Pion->Fill(kfpPionMinus.GetPt(), zPion_KF_MC);

        fHistChi2ndfProton->Fill(kfpProton.GetChi2()/kfpProton.GetNDF());
        fHistProbProton->Fill(TMath::Prob(kfpProton.GetChi2(), kfpProton.GetNDF()));
        fHistChi2ndfPion->Fill(kfpPionMinus.GetChi2()/kfpPionMinus.GetNDF());
        fHistProbPion->Fill(TMath::Prob(kfpPionMinus.GetChi2(), kfpPionMinus.GetNDF()));

        if (kfpProton.GetCovariance(0,0)>0) {
          Double_t err_xProton = TMath::Sqrt(kfpProton.GetCovariance(0,0));
          Double_t xProton_PULL = xProton_KF_MC/err_xProton;
          fHistXProton_PULL->Fill(xProton_PULL);
          f2DHistXPULLvsPt_Proton->Fill(kfpProton.GetPt(), xProton_PULL);
        }
        if (kfpProton.GetCovariance(1,1)>0) {
          Double_t err_yProton = TMath::Sqrt(kfpProton.GetCovariance(1,1));
          Double_t yProton_PULL = yProton_KF_MC/err_yProton;
          fHistYProton_PULL->Fill(yProton_PULL);
          f2DHistYPULLvsPt_Proton->Fill(kfpProton.GetPt(), yProton_PULL);
        }
        if (kfpProton.GetCovariance(2,2)>0) {
          Double_t err_zProton = TMath::Sqrt(kfpProton.GetCovariance(2,2));
          Double_t zProton_PULL = zProton_KF_MC/err_zProton;
          fHistZProton_PULL->Fill(zProton_PULL);
          f2DHistZPULLvsPt_Proton->Fill(kfpProton.GetPt(), zProton_PULL);
        }

        if (kfpPionMinus.GetCovariance(0,0)>0) {
          Double_t err_xPion = TMath::Sqrt(kfpPionMinus.GetCovariance(0,0));
          Double_t xPion_PULL = xPion_KF_MC/err_xPion;
          fHistXPion_PULL->Fill(xPion_PULL);
          f2DHistXPULLvsPt_Pion->Fill(kfpPionMinus.GetPt(), xPion_PULL);
        }
        if (kfpPionMinus.GetCovariance(1,1)>0) {
          Double_t err_yPion = TMath::Sqrt(kfpPionMinus.GetCovariance(1,1));
          Double_t yPion_PULL = yPion_KF_MC/err_yPion;
          fHistYPion_PULL->Fill(yPion_PULL);
          f2DHistYPULLvsPt_Pion->Fill(kfpPionMinus.GetPt(), yPion_PULL);
        }
        if (kfpPionMinus.GetCovariance(2,2)>0) {
          Double_t err_zPion = TMath::Sqrt(kfpPionMinus.GetCovariance(2,2));
          Double_t zPion_PULL = zPion_KF_MC/err_zPion;
          fHistZPion_PULL->Fill(zPion_PULL);
          f2DHistZPULLvsPt_Pion->Fill(kfpPionMinus.GetPt(), zPion_PULL);
        }
        //=====================================================================

        Double_t xLambda_KF_MC = kfpLambda.GetX() - mcProton->Xv();
        Double_t yLambda_KF_MC = kfpLambda.GetY() - mcProton->Yv();
        Double_t zLambda_KF_MC = kfpLambda.GetZ() - mcProton->Zv();

        Double_t xLambda_V0_MC = v0->DecayVertexV0X() - mcProton->Xv();

        Double_t ptLambda_KF_MC = ptLambda_Rec - mcLambda->Pt();
        Double_t ptLambda_KF_MC_R = ptLambda_KF_MC/ptLambda_Rec;

        Double_t ptLambda_V0_MC = v0->Pt() - mcLambda->Pt();
        Double_t ptLambda_V0_MC_R = ptLambda_V0_MC/(v0->Pt());

        fHistXLambda_KF_MC->Fill(xLambda_KF_MC);
        fHistYLambda_KF_MC->Fill(yLambda_KF_MC);
        fHistZLambda_KF_MC->Fill(zLambda_KF_MC);
        f2DHistXRecMCvsPt_Lambda->Fill(ptLambda_Rec, xLambda_KF_MC);
        f2DHistPtRecMCvsPt_Lambda->Fill(ptLambda_Rec, ptLambda_KF_MC_R);

        fHistXLambda_V0_MC->Fill(xLambda_V0_MC);
        f2DHistXV0MCvsPt_Lambda->Fill(v0->Pt(), xLambda_V0_MC);
        f2DHistPtV0MCvsPt_Lambda->Fill(v0->Pt(), ptLambda_V0_MC_R);


        if (err_ptLambda>0) {
          Double_t ptLambda_PULL = ptLambda_KF_MC/err_ptLambda;
          f2DHistPtPULLvsPt_Lambda->Fill(ptLambda_Rec, ptLambda_PULL);
        }

        if (err_massLambda>0) {
          Double_t massLambda_KF_MC = massLambda_Rec - massLambda;
          Double_t massLambda_V0_MC = v0->MassLambda() - massLambda;
          fHistMassLambda_KF_MC->Fill(massLambda_KF_MC);
          fHistMassLambda_V0_MC->Fill(massLambda_V0_MC);
          f2DHistMassRecMCvsPt_Lambda->Fill(ptLambda_Rec, massLambda_KF_MC);
          f2DHistMassV0MCvsPt_Lambda->Fill(v0->Pt(), massLambda_V0_MC);
          Double_t massLambda_PULL = massLambda_KF_MC/err_massLambda;
          fHistMassLambda_PULL_KF->Fill(massLambda_PULL);
          f2DHistMassPULLvsPt_Lambda->Fill(ptLambda_Rec, massLambda_PULL);
          Double_t radius_Lambda = TMath::Sqrt(kfpLambda.GetX()*kfpLambda.GetX()+kfpLambda.GetY()*kfpLambda.GetY());
          f2DHistMassPULLvsRadius_Lambda->Fill(radius_Lambda, massLambda_PULL);
        }
        if (kfpLambda.GetCovariance(0,0)>0) {
          Double_t err_xLambda = TMath::Sqrt(kfpLambda.GetCovariance(0,0));
          Double_t xLambda_PULL = xLambda_KF_MC/err_xLambda;
          fHistXLambda_PULL->Fill(xLambda_PULL);
          f2DHistXPULLvsPt_Lambda->Fill(ptLambda_Rec, xLambda_PULL);
        }
        if (kfpLambda.GetCovariance(1,1)>0) {
          Double_t err_yLambda = TMath::Sqrt(kfpLambda.GetCovariance(1,1));
          Double_t yLambda_PULL = yLambda_KF_MC/err_yLambda;
          fHistYLambda_PULL->Fill(yLambda_PULL);
        }
        if (kfpLambda.GetCovariance(2,2)>0) {
          Double_t err_zLambda = TMath::Sqrt(kfpLambda.GetCovariance(2,2));
          Double_t zLambda_PULL = zLambda_KF_MC/err_zLambda;
          fHistZLambda_PULL->Fill(zLambda_PULL);
        }



//        Double_t dz_Lambda[2], cov[3];
//        AliAODVertex *vtx_v0 = reinterpret_cast<AliAODVertex*>(v0);
//        trkP->PropagateToDCA(vtx_v0, fBzkG, kVeryBig, dz_Lambda, cov);
//        Double_t xyz_new[3]={-999,-999,-999};
//        trkP->GetXYZ(xyz_new);
//        Double_t xProton_V0_MC = xyz_new[0] - mcProton->Xv();
//        fHistXLambda_V0_MC->Fill(xProton_V0_MC);
//        f2DHistXV0MCvsPt_Lambda->Fill(trkP->Pt(), xProton_V0_MC);
//        if (cov[0]>0) {
//          Double_t err_xProton_V0 = TMath::Sqrt(cov[0]);
//          Double_t xProton_PULL_V0 = xProton_V0_MC/err_xProton_V0;
//          f2DHistXPULLvsPt_Lambda_V0->Fill(trkP->Pt(), xProton_PULL_V0);
//        }
        




      }

      // Anti-Lambda
      if (labelAntiLambda>-1) {
        fHistMassAntiLambda_Match->Fill(kfpAntiLambda.GetMass());
        fHistMassLambdaTot_Match->Fill(kfpAntiLambda.GetMass());
        fHistChi2ndfAntiLambda_Match->Fill(kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF());
        fHistProbAntiLambda_Match->Fill(TMath::Prob(kfpAntiLambda.GetChi2(), kfpAntiLambda.GetNDF()));

        AliAODMCParticle* mcAntiLambda = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiLambda));
        AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(mcAntiLambda->GetDaughterFirst()));
        AliAODMCParticle* mcAntiPion   = static_cast<AliAODMCParticle*>(mcArray->At(mcAntiLambda->GetDaughterLast()));
//        cout << "PDG antiProton ==" << mcAntiProton->GetPdgCode() << endl;
//        cout << "PDG antiPion ==" << mcAntiPion->GetPdgCode() << endl;

        Double_t xAntiLambda_Rec_MC = kfpAntiLambda.GetX() - mcAntiProton->Xv();
        Double_t yAntiLambda_Rec_MC = kfpAntiLambda.GetY() - mcAntiProton->Yv();
        Double_t zAntiLambda_Rec_MC = kfpAntiLambda.GetZ() - mcAntiProton->Zv();
        Double_t ptAntiLmabda_Rec_MC = ptAntiLambda_Rec - mcAntiLambda->Pt();
        Double_t ptAntiLmabda_Rec_MC_R = ptAntiLmabda_Rec_MC/ptAntiLambda_Rec;

        fHistXAntiLambda_Rec_MC->Fill(xAntiLambda_Rec_MC);
        fHistYAntiLambda_Rec_MC->Fill(yAntiLambda_Rec_MC);
        fHistZAntiLambda_Rec_MC->Fill(zAntiLambda_Rec_MC);
        f2DHistXRecMCvsPt_AntiLambda->Fill(ptAntiLambda_Rec, xAntiLambda_Rec_MC);
        f2DHistPtRecMCvsPt_AntiLambda->Fill(ptAntiLambda_Rec, ptAntiLmabda_Rec_MC_R);

        if (err_ptAntiLambda>0) {
          Double_t ptAntiLambda_PULL = ptAntiLmabda_Rec_MC/err_ptAntiLambda;
          f2DHistPtPULLvsPt_AntiLambda->Fill(ptAntiLambda_Rec, ptAntiLambda_PULL);
        }
        if (err_massAntiLambda>0) {
          Double_t massAntiLambda_KF_MC = massAntiLambda_Rec - massLambda;
          fHistMassAntiLambda_KF_MC->Fill(massAntiLambda_KF_MC);
          f2DHistMassRecMCvsPt_AntiLambda->Fill(ptAntiLambda_Rec, massAntiLambda_KF_MC);
          Double_t massAntiLambda_PULL = massAntiLambda_KF_MC/err_massAntiLambda;
          fHistMassAntiLambda_PULL_KF->Fill(massAntiLambda_PULL);
          f2DHistMassPULLvsPt_AntiLambda->Fill(ptAntiLambda_Rec, massAntiLambda_PULL);
          Double_t radius_AntiLambda = TMath::Sqrt(kfpAntiLambda.GetX()*kfpAntiLambda.GetX()+kfpAntiLambda.GetY()*kfpAntiLambda.GetY());
          f2DHistMassPULLvsRadius_AntiLambda->Fill(radius_AntiLambda, massAntiLambda_PULL);
        }
        if (kfpAntiLambda.GetCovariance(0,0)>0) {
          Double_t err_xAntiLambda = TMath::Sqrt(kfpAntiLambda.GetCovariance(0,0));
          Double_t xAntiLambda_PULL = xAntiLambda_Rec_MC/err_xAntiLambda;
          fHistXAntiLambda_PULL->Fill(xAntiLambda_PULL);
          f2DHistXPULLvsPt_AntiLambda->Fill(ptAntiLambda_Rec, xAntiLambda_PULL);
        }
        if (kfpAntiLambda.GetCovariance(1,1)>0) {
          Double_t err_yAntiLambda = TMath::Sqrt(kfpAntiLambda.GetCovariance(1,1));
          Double_t yAntiLambda_PULL = yAntiLambda_Rec_MC/err_yAntiLambda;
          fHistYAntiLambda_PULL->Fill(yAntiLambda_PULL);
        }
        if (kfpAntiLambda.GetCovariance(2,2)>0) {
          Double_t err_zAntiLambda = TMath::Sqrt(kfpAntiLambda.GetCovariance(2,2));
          Double_t zAntiLambda_PULL = zAntiLambda_Rec_MC/err_zAntiLambda;
          fHistZAntiLambda_PULL->Fill(zAntiLambda_PULL);
        }
      }
    }
*/






/*
    if (kfpLambda.GetChi2()<0) {
      cout << "==========" << endl;
      cout << "Positive track (daughter1)" << endl;
      cout << covP[0] << endl;
      cout << covP[1] << "  " << covP[2] << endl;
      cout << covP[3] << "  " << covP[4] << "  " << covP[5] << endl;
      cout << covP[6] << "  " << covP[7] << "  " << covP[8] << "  " << covP[9] << endl;
      cout << covP[10] << "  " << covP[11] << "  " << covP[12] << "  " << covP[13] << "  " << covP[14] << endl;
      cout << covP[15] << "  " << covP[16] << "  " << covP[17] << "  " << covP[18] << "  " << covP[19] << "  " << covP[20] << endl;
      cout << "Negative track (daughter2)" << endl;
      cout << covN[0] << endl;
      cout << covN[1] << "  " << covN[2] << endl;
      cout << covN[3] << "  " << covN[4] << "  " << covN[5] << endl;
      cout << covN[6] << "  " << covN[7] << "  " << covN[8] << "  " << covN[9] << endl;
      cout << covN[10] << "  " << covN[11] << "  " << covN[12] << "  " << covN[13] << "  " << covN[14] << endl;
      cout << covN[15] << "  " << covN[16] << "  " << covN[17] << "  " << covN[18] << "  " << covN[19] << "  " << covN[20] << endl;
      std::cout << "Lambda at DP" << std::endl << "  " << kfpLambda << std::endl;
    }

    if (kfpLambda.GetNDF()<0) {
      cout << "==========" << endl;
      cout << "Positive track (daughter1)" << endl;
      cout << kfpProton << endl;
      cout << "Negative track (daughter2)" << endl;
      cout << kfpPionMinus << endl;
      std::cout << "Lambda at DP" << std::endl << "  " << kfpLambda << std::endl;
    }
*/





/*
    f2DHistChi2vsNDF_Lambda->Fill(kfpLambda.GetNDF(), kfpLambda.GetChi2());
    f2DHistChi2vsNDF_AntiLambda->Fill(kfpAntiLambda.GetNDF(), kfpAntiLambda.GetChi2());
    fHistChi2ndfLambda->Fill(kfpLambda.GetChi2()/kfpLambda.GetNDF());
    fHistChi2ndfAntiLambda->Fill(kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF());
    fHistChi2ndfLambdaTot->Fill(kfpLambda.GetChi2()/kfpLambda.GetNDF());
    fHistChi2ndfLambdaTot->Fill(kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF());

    fHistProbLambda->Fill(TMath::Prob(kfpLambda.GetChi2(), kfpLambda.GetNDF()));
    fHistProbAntiLambda->Fill(TMath::Prob(kfpAntiLambda.GetChi2(), kfpAntiLambda.GetNDF()));
    fHistProbLambdaTot->Fill(TMath::Prob(kfpLambda.GetChi2(), kfpLambda.GetNDF()));
    fHistProbLambdaTot->Fill(TMath::Prob(kfpAntiLambda.GetChi2(), kfpAntiLambda.GetNDF()));

    Float_t massLambda_Rec, err_massLambda;
    Float_t massAntiLambda_Rec, err_massAntiLambda;
    Float_t ptLambda_Rec, err_ptLambda;
    Float_t ptAntiLambda_Rec, err_ptAntiLambda;

    kfpLambda.GetMass(massLambda_Rec, err_massLambda);
    kfpAntiLambda.GetMass(massAntiLambda_Rec, err_massAntiLambda);
    kfpLambda.GetPt(ptLambda_Rec, err_ptLambda);
    kfpAntiLambda.GetPt(ptAntiLambda_Rec, err_ptAntiLambda);

//    fHistMassLambda->Fill(kfpLambda.GetMass());
    fHistMassAntiLambda->Fill(kfpAntiLambda.GetMass());
    fHistMassLambdaTot->Fill(kfpLambda.GetMass());
    fHistMassLambdaTot->Fill(kfpAntiLambda.GetMass());

    fHistMassLambda_V0->Fill(v0->MassLambda());
    fHistMassAntiLambda_V0->Fill(v0->MassAntiLambda());
    fHistMassLambdaTot_V0->Fill(v0->MassLambda());
    fHistMassLambdaTot_V0->Fill(v0->MassAntiLambda());

    Double_t Lambda_KF_V0 = kfpLambda.GetMass() - v0->MassLambda();
    Double_t AntiLambda_KF_V0 = kfpAntiLambda.GetMass() - v0->MassAntiLambda();

    fHistMassLambda_KF_V0->Fill(Lambda_KF_V0);
    fHistMassAntiLambda_KF_V0->Fill(AntiLambda_KF_V0);

    fHistPtLambda->Fill(kfpLambda.GetPt());
    fHistPtAntiLambda->Fill(kfpAntiLambda.GetPt());
    fHistPtLambdaTot->Fill(kfpLambda.GetPt());
    fHistPtLambdaTot->Fill(kfpAntiLambda.GetPt());
*/



//      std::cout << "Proton at DP" << std::endl << "  " << kfpProton << std::endl;
//      std::cout << "Pion3 at DP" << std::endl << "  " << kfpPion3 << std::endl;
//      std::cout << "Lambda at DP" << std::endl << "  " << kfpLambda << std::endl;


      // set linearised mass constraint
//      kfpLambda.SetMassConstraint(massLambda); 
//      kfpAntiLambda.SetMassConstraint(massLambda); 
      // set the exact nonlinear mass constraint

/*
      if ( kfpLambda.GetChi2()/kfpLambda.GetNDF()<3 && kfpLambda.GetMass()>0 ) {
        kfpLambda.SetNonlinearMassConstraint(massLambda);
        fHistMassLambda_M->Fill(kfpLambda.GetMass());
        fHistMassLambdaTot_M->Fill(kfpLambda.GetMass());
      }
      if ( kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF()<3 && kfpAntiLambda.GetMass()>0 ) {
        kfpAntiLambda.SetNonlinearMassConstraint(massLambda);
        fHistMassAntiLambda_M->Fill(kfpAntiLambda.GetMass());
        fHistMassLambdaTot_M->Fill(kfpAntiLambda.GetMass());
      }
*/

    kfpK0Short.Clear();
    kfpAntiLambda.Clear();
    kfpLambda.Clear();
    kfpPionPlus.Clear();
    kfpPionMinus.Clear();
    kfpAntiProton.Clear();
    kfpProton.Clear();
  }

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::MakeAnaXicZeroFromCasc(AliAODEvent *AODEvent, TClonesArray *mcArray, KFParticle PV)
{
  // Main analysis called from "UserExec"

//  std::cout.setf(std::ios::fixed);
//  std::cout.setf(std::ios::showpoint);
//  std::cout.precision(3);

  // set the magnetic field
  KFParticle::SetField(fBzkG);

  const UInt_t nCasc = AODEvent->GetNumberOfCascades();

  Double_t xyzP[3], xyzN[3];
  Double_t xvyvzvP[3], xvyvzvN[3];
  Double_t covP[21], covN[21], covB[21];
  const Int_t NDaughters = 2;
  const Float_t massLambda = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  const Float_t massXi     = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  const Float_t massK0S    = TDatabasePDG::Instance()->GetParticle(310)->Mass();

  // select good candidates for pion
  const UInt_t nTracks = AODEvent->GetNumberOfTracks();
  AliAODTrack *trackP[nTracks], *trackN[nTracks];
  Int_t flag_trkP = 0, flag_trkN = 0;

  for (UInt_t itrk=0; itrk<nTracks; itrk++) {
    AliAODTrack *trk = static_cast<AliAODTrack*>(AODEvent->GetTrack(itrk));
    Double_t covtest[21];
    if ( !trk || trk->GetID()<0 || !trk->GetCovarianceXYZPxPyPz(covtest) || !CheckTrackCov(trk) ) continue;

    if  (trk->Charge() > 0 ) {
      trackP[flag_trkP] = trk;
      flag_trkP++;
    }
    if  (trk->Charge() < 0 ) {
      trackN[flag_trkN] = trk;
      flag_trkN++;
    }
  }

  for (UInt_t iCasc=0; iCasc<nCasc; iCasc++) {
    AliAODcascade *casc = AODEvent->GetCascade(iCasc);
    // cascade cut
    if ( !fAnaCuts->SingleCascCuts(casc) ) continue;

    AliAODTrack *ptrack = (AliAODTrack*) (casc->GetDaughter(0));
    AliAODTrack *ntrack = (AliAODTrack*) (casc->GetDaughter(1));
    AliAODTrack *btrack = (AliAODTrack*) (casc->GetDecayVertexXi()->GetDaughter(0));

    if ( !ptrack||!ntrack||!btrack ) continue;

    // check charge of the first daughter, if negative, define it as the second one
    if ( ptrack->Charge()<0 ) {
      ptrack = (AliAODTrack*) (casc->GetDaughter(1));
      ntrack = (AliAODTrack*) (casc->GetDaughter(0));
    }

    if ( !ptrack->GetCovarianceXYZPxPyPz(covP) || !ntrack->GetCovarianceXYZPxPyPz(covN) || !btrack->GetCovarianceXYZPxPyPz(covB) ) continue;

    if ( !CheckTrackCov(ptrack) || !CheckTrackCov(ntrack) || !CheckTrackCov(btrack) ) continue;

    KFParticle kfpProton     = CreateKFParticleFromAODtrack(ptrack, 2212);
    KFParticle kfpPionMinus  = CreateKFParticleFromAODtrack(ntrack, -211);
    KFParticle kfpAntiProton = CreateKFParticleFromAODtrack(ntrack, -2212);
    KFParticle kfpPionPlus   = CreateKFParticleFromAODtrack(ptrack, 211);

    KFParticle kfpElePlus    = CreateKFParticleFromAODtrack(ptrack, -11);
    KFParticle kfpEleMinus   = CreateKFParticleFromAODtrack(ntrack, 11);

    // === K0S ===
    KFParticle kfpK0Short;
    const KFParticle *vk0sDaughters[2]  = {&kfpPionPlus, &kfpPionMinus};
    kfpK0Short.Construct(vk0sDaughters, NDaughters);
    // ============
    // === Gamma ===
    KFParticle kfpGamma;
    const KFParticle *vGammaDaughters[2]  = {&kfpElePlus, &kfpEleMinus};
    kfpGamma.Construct(vGammaDaughters, NDaughters);
    // =============

    if ( btrack->Charge()<0 ) {

      const KFParticle *vDaughters[2] = {&kfpProton, &kfpPionMinus};

      KFParticle kfpLambda;
      kfpLambda.Construct(vDaughters, NDaughters);
      Float_t massLambda_Rec, err_massLambda;
      kfpLambda.GetMass(massLambda_Rec, err_massLambda);

      // check rapidity of lambda
      if ( TMath::Abs(kfpLambda.GetE())<=TMath::Abs(kfpLambda.GetPz()) ) continue;

      // chi2>0 && NDF>0 for selecting Lambda
      if ( (kfpLambda.GetNDF()<=0 || kfpLambda.GetChi2()<=0) ) continue;

      // check cov. of Lambda
      if ( !CheckKFParticleCov(kfpLambda) ) continue;

      // err_mass>0 of Lambda
      if ( err_massLambda<=0 ) continue;

      // Chi2geo cut of Lambda
      if ( (kfpLambda.GetChi2()/kfpLambda.GetNDF()) >= fAnaCuts->GetKFPLam_Chi2geoMax() ) continue; 

      //************************** calculate l/Δl for Lambda *************************************
      Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
      Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
      Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
      Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
      Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
      if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
      dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
      Double_t nErr_l_Lambda = l_Lambda/dl_Lambda;
      //***************************************************************************************

      // l/Deltal cut of Lambda
      if ( nErr_l_Lambda <= fAnaCuts->GetKFPLam_lDeltalMin() ) continue;

      // mass window cut of Lambda
      if ( TMath::Abs(massLambda_Rec-massLambda) > (fAnaCuts->GetProdMassTolLambda()) ) continue;

      KFParticle kfpLambda_m = kfpLambda;
      kfpLambda_m.SetNonlinearMassConstraint(massLambda);

      if ( !CheckKFParticleCov(kfpLambda_m) || TMath::Abs(kfpLambda.GetE()) <= TMath::Abs(kfpLambda.GetPz()) ) continue;

      KFParticle kfpPion2 = CreateKFParticleFromAODtrack(btrack, -211);
      KFParticle kfpXiMinus;
      const KFParticle *vXiDs[2] = {&kfpPion2, &kfpLambda_m};
      kfpXiMinus.Construct(vXiDs, NDaughters);

      // check rapidity of Xi-
      if ( TMath::Abs(kfpXiMinus.GetE())<=TMath::Abs(kfpXiMinus.GetPz()) ) continue;

      // err_massXi > 0
      Float_t massXiMinus_Rec, err_massXiMinus;
      kfpXiMinus.GetMass(massXiMinus_Rec, err_massXiMinus);
      if ( err_massXiMinus<=0 ) continue;

      // chi2>0 && NDF>0
      if ( kfpXiMinus.GetNDF()<=0 || kfpXiMinus.GetChi2()<=0 ) continue;

      // Prefilter
      if ( kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF() >= fAnaCuts->GetKFPXi_Chi2geoMax() ) continue;

      // check covariance matrix
      if ( !CheckKFParticleCov(kfpXiMinus) ) continue;

      // mass window cut of Xi-
      if ( TMath::Abs(massXiMinus_Rec-massXi) > (fAnaCuts->GetProdMassTolXi()) ) continue;

      KFParticle kfpXiMinus_m = kfpXiMinus;
      kfpXiMinus_m.SetNonlinearMassConstraint(massXi);

      if ( !CheckKFParticleCov(kfpXiMinus_m) || TMath::Abs(kfpXiMinus.GetE()) <= TMath::Abs(kfpXiMinus.GetPz()) ) continue;

      for (Int_t itrkBP=0; itrkBP<flag_trkP; itrkBP++) { // Loop for first bachelor pion+

        if ( trackP[itrkBP]->GetID() == ptrack->GetID() ) continue;

        if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackP[itrkBP]) ) continue;

        KFParticle kfpBP = CreateKFParticleFromAODtrack(trackP[itrkBP], 211);

        // reconstruct Xic0
        KFParticle kfpXic0;
        const KFParticle *vXicZeroDs[2] = {&kfpBP, &kfpXiMinus_m};
        kfpXic0.Construct(vXicZeroDs, NDaughters);

        // chi2>0 && NDF>0
        if ( kfpXic0.GetNDF()<=0 || kfpXic0.GetChi2()<=0 ) continue;

        // Prefilter
        if ( kfpXic0.GetChi2()/kfpXic0.GetNDF() >= fAnaCuts->GetKFPXic0_Chi2geoMax() ) continue;
        if ( kfpXic0.GetPt() < fAnaCuts->GetPtMinXic0() ) continue;

        // check rapidity of Xic0
        if ( TMath::Abs(kfpXic0.GetE())<=TMath::Abs(kfpXic0.GetPz()) ) continue;

        // check covariance matrix
        if ( !CheckKFParticleCov(kfpXic0) ) continue;

        // err_massXic0 > 0
        Float_t massXic0_Rec, err_massXic0;
        kfpXic0.GetMass(massXic0_Rec, err_massXic0);
        if ( err_massXic0<=0 ) continue;

        if (fWriteXic0Tree) {
           Int_t lab_Xic0 = -9999;
           if (fIsMC) {
             lab_Xic0 = MatchToMCXic0(ptrack, ntrack, btrack, trackP[itrkBP], mcArray);
             if (lab_Xic0>=0) FillTreeRecXic0FromCasc(kfpXic0, trackP[itrkBP], kfpBP, kfpXiMinus, kfpXiMinus_m, btrack, casc, kfpK0Short, kfpGamma, kfpLambda, kfpLambda_m, ptrack, ntrack, PV, mcArray, lab_Xic0);
           }
           if (!fIsMC) {
             FillTreeRecXic0FromCasc(kfpXic0, trackP[itrkBP], kfpBP, kfpXiMinus, kfpXiMinus_m, btrack, casc, kfpK0Short, kfpGamma, kfpLambda, kfpLambda_m, ptrack, ntrack, PV, mcArray, lab_Xic0);
           }
         }
         kfpXic0.Clear();
         kfpBP.Clear();
      } // Loop for first bachelor pion+
      kfpXiMinus_m.Clear();
      kfpXiMinus.Clear();
      kfpPion2.Clear();
      kfpLambda_m.Clear();
      kfpLambda.Clear();
    }

    if ( btrack->Charge()>0 ) {

      const KFParticle *vAntiDaughters[2] = {&kfpPionPlus, &kfpAntiProton};

      KFParticle kfpAntiLambda;
      kfpAntiLambda.Construct(vAntiDaughters, NDaughters);
      Float_t massAntiLambda_Rec, err_massAntiLambda;
      kfpAntiLambda.GetMass(massAntiLambda_Rec, err_massAntiLambda);

      // check rapidity of Anti-Lambda
      if ( TMath::Abs(kfpAntiLambda.GetE())<=TMath::Abs(kfpAntiLambda.GetPz()) ) continue;

      // chi2>0 && NDF>0 for selecting Anti-Lambda
      if ( kfpAntiLambda.GetNDF()<=0 || kfpAntiLambda.GetChi2()<=0 ) continue;

      // check cov. of Anti-Lambda
      if ( !CheckKFParticleCov(kfpAntiLambda) ) continue;

      // err_mass>0 of Anti-Lambda
      if ( err_massAntiLambda<=0 ) continue;

      // Chi2geo cut of Anti-Lambda
      if ( (kfpAntiLambda.GetChi2()/kfpAntiLambda.GetNDF()) >= fAnaCuts->GetKFPLam_Chi2geoMax() ) continue;

      //************************** calculate l/Δl for Anti-Lambda *************************************
      Double_t dx_AntiLambda = PV.GetX()-kfpAntiLambda.GetX();
      Double_t dy_AntiLambda = PV.GetY()-kfpAntiLambda.GetY();
      Double_t dz_AntiLambda = PV.GetZ()-kfpAntiLambda.GetZ();
      Double_t l_AntiLambda = TMath::Sqrt(dx_AntiLambda*dx_AntiLambda + dy_AntiLambda*dy_AntiLambda + dz_AntiLambda*dz_AntiLambda);
      Double_t dl_AntiLambda = (PV.GetCovariance(0)+kfpAntiLambda.GetCovariance(0))*dx_AntiLambda*dx_AntiLambda + (PV.GetCovariance(2)+kfpAntiLambda.GetCovariance(2))*dy_AntiLambda*dy_AntiLambda + (PV.GetCovariance(5)+kfpAntiLambda.GetCovariance(5))*dz_AntiLambda*dz_AntiLambda + 2*( (PV.GetCovariance(1)+kfpAntiLambda.GetCovariance(1))*dx_AntiLambda*dy_AntiLambda + (PV.GetCovariance(3)+kfpAntiLambda.GetCovariance(3))*dx_AntiLambda*dz_AntiLambda + (PV.GetCovariance(4)+kfpAntiLambda.GetCovariance(4))*dy_AntiLambda*dz_AntiLambda );
      if ( fabs(l_AntiLambda)<1.e-8f ) l_AntiLambda = 1.e-8f;
      dl_AntiLambda = dl_AntiLambda<0. ? 1.e8f : sqrt(dl_AntiLambda)/l_AntiLambda;
      Double_t nErr_l_AntiLambda = l_AntiLambda/dl_AntiLambda;
      //***************************************************************************************

      // l/Deltal cut of Anti-Lambda
      if ( nErr_l_AntiLambda <= fAnaCuts->GetKFPLam_lDeltalMin() ) continue;

      // mass window cut of Anti-Lambda
      if ( TMath::Abs(massAntiLambda_Rec-massLambda) > (fAnaCuts->GetProdMassTolLambda()) ) continue;

      KFParticle kfpAntiLambda_m = kfpAntiLambda;
      kfpAntiLambda_m.SetNonlinearMassConstraint(massLambda);

      if ( !CheckKFParticleCov(kfpAntiLambda_m) || TMath::Abs(kfpAntiLambda.GetE()) <= TMath::Abs(kfpAntiLambda.GetPz()) ) continue;

      KFParticle kfpPion2 = CreateKFParticleFromAODtrack(btrack, 211);
      KFParticle kfpXiPlus;
      const KFParticle *vXiDs[2] = {&kfpPion2, &kfpAntiLambda_m};
      kfpXiPlus.Construct(vXiDs, NDaughters);

      // check rapidity of Xi+
      if ( TMath::Abs(kfpXiPlus.GetE())<=TMath::Abs(kfpXiPlus.GetPz()) ) continue;

      // err_massXi > 0
      Float_t massXiPlus_Rec, err_massXiPlus;
      kfpXiPlus.GetMass(massXiPlus_Rec, err_massXiPlus);
      if ( err_massXiPlus<=0 ) continue;

      // chi2>0 && NDF>0
      if ( kfpXiPlus.GetNDF()<=0 || kfpXiPlus.GetChi2()<=0 ) continue;

      // Prefilter
      if ( kfpXiPlus.GetChi2()/kfpXiPlus.GetNDF() >= fAnaCuts->GetKFPXi_Chi2geoMax() ) continue;

      // check covariance matrix
      if ( !CheckKFParticleCov(kfpXiPlus) ) continue;

      // mass window cut of Xi+
      if ( TMath::Abs(massXiPlus_Rec-massXi) > (fAnaCuts->GetProdMassTolXi()) ) continue;

      KFParticle kfpXiPlus_m = kfpXiPlus;
      kfpXiPlus_m.SetNonlinearMassConstraint(massXi);

      if ( !CheckKFParticleCov(kfpXiPlus_m) || TMath::Abs(kfpXiPlus.GetE()) <= TMath::Abs(kfpXiPlus.GetPz()) ) continue;

      for (Int_t itrkBP=0; itrkBP<flag_trkN; itrkBP++) { // Loop for first bachelor pion-

        if ( trackN[itrkBP]->GetID() == ntrack->GetID() ) continue;

        if ( !fAnaCuts->PassedTrackQualityCuts_PrimaryPion(trackN[itrkBP]) ) continue;

        KFParticle kfpBP = CreateKFParticleFromAODtrack(trackN[itrkBP], -211);

        // reconstruct Anti-Xic0
        KFParticle kfpAntiXic0;
        const KFParticle *vXic0Ds[2] = {&kfpBP, &kfpXiPlus_m};
        kfpAntiXic0.Construct(vXic0Ds, NDaughters);

        // chi2>0 && NDF>0
        if ( kfpAntiXic0.GetNDF()<=0 || kfpAntiXic0.GetChi2()<=0 ) continue;

        // Prefilter
        if ( kfpAntiXic0.GetChi2()/kfpAntiXic0.GetNDF() >= fAnaCuts->GetKFPXic0_Chi2geoMax() ) continue;
        if ( kfpAntiXic0.GetPt() < fAnaCuts->GetPtMinXic0() ) continue;

        // check rapidity of Anti-Xic0
        if ( TMath::Abs(kfpAntiXic0.GetE())<=TMath::Abs(kfpAntiXic0.GetPz()) ) continue;

        // check covariance matrix
        if ( !CheckKFParticleCov(kfpAntiXic0) ) continue;

        // err_massAntiXic0 > 0
        Float_t massAntiXic0_Rec, err_massAntiXic0;
        kfpAntiXic0.GetMass(massAntiXic0_Rec, err_massAntiXic0);
        if ( err_massAntiXic0<=0 ) continue;

        if (fWriteXic0Tree) {
          Int_t lab_AntiXic0 = -9999.;
          if (fIsMC) {
            lab_AntiXic0 = MatchToMCAntiXic0(ntrack, ptrack, btrack, trackN[itrkBP], mcArray);
            if (lab_AntiXic0>=0) FillTreeRecXic0FromCasc(kfpAntiXic0, trackN[itrkBP], kfpBP, kfpXiPlus, kfpXiPlus_m, btrack, casc, kfpK0Short, kfpGamma, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, PV, mcArray, lab_AntiXic0);
          }
          if (!fIsMC) {
            FillTreeRecXic0FromCasc(kfpAntiXic0, trackN[itrkBP], kfpBP, kfpXiPlus, kfpXiPlus_m, btrack, casc, kfpK0Short, kfpGamma, kfpAntiLambda, kfpAntiLambda_m, ntrack, ptrack, PV, mcArray, lab_AntiXic0);
          }
        }
        kfpAntiXic0.Clear();
        kfpBP.Clear();
      } // Loop for first bachelor pion-
      kfpXiPlus_m.Clear();
      kfpXiPlus.Clear();
      kfpPion2.Clear();
      kfpAntiLambda_m.Clear();
      kfpAntiLambda.Clear();
    }

    kfpK0Short.Clear();
    kfpPionPlus.Clear();
    kfpAntiProton.Clear();
    kfpPionMinus.Clear();
    kfpProton.Clear();

  }

  return;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicZero2XiPifromKFP::MatchToMCXic0(AliAODTrack *trackProton, AliAODTrack *trackPion3, AliAODTrack *trackPion2, AliAODTrack *trackAntiPion1, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle
  
  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
  Int_t labelPion3  = fabs(trackPion3->GetLabel());
  if (labelPion3<0) return -1;
  AliAODMCParticle* mcPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion3));
  Int_t labelPion2  = fabs(trackPion2->GetLabel());
  if (labelPion2<0) return -1;
  AliAODMCParticle* mcPion2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion2));
  Int_t labelAntiPion1  = fabs(trackAntiPion1->GetLabel());
  if (labelAntiPion1<0) return -1;
  AliAODMCParticle* mcAntiPion1 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion1));

  if ( mcProton->GetPdgCode() != 2212 || mcPion3->GetPdgCode() != -211 || mcPion2->GetPdgCode() != -211 || mcAntiPion1->GetPdgCode() != 211) return -1; // check pdg

  Int_t IndexMother[2];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPion3->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is lambda

  IndexMother[0] = mcMother->GetMother(); // mother of lambda
  IndexMother[1] = mcPion2->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3312 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xi-

  IndexMother[0] = mcMother->GetMother(); // mother of Xi-
  IndexMother[1] = mcAntiPion1->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 4132 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xic0

//  if ( mcMother->IsPrimary() ) return 1;
//  if ( mcMother->IsPhysicalPrimary() ) return 2;
//  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
//  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
//  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

  Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcMother,kTRUE);
  return CheckOrigin;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicZero2XiPifromKFP::MatchToMCAntiXic0(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, AliAODTrack *trackAntiPion2, AliAODTrack *trackPion1, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle
  
  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
  Int_t labelAntiPion3  = fabs(trackAntiPion3->GetLabel());
  if (labelAntiPion3<0) return -1;
  AliAODMCParticle* mcAntiPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion3));
  Int_t labelAntiPion2  = fabs(trackAntiPion2->GetLabel());
  if (labelAntiPion2<0) return -1;
  AliAODMCParticle* mcAntiPion2 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion2));
  Int_t labelPion1  = fabs(trackPion1->GetLabel());
  if (labelPion1<0) return -1;
  AliAODMCParticle* mcPion1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion1));
  if ( mcAntiProton->GetPdgCode() != -2212 || mcAntiPion3->GetPdgCode() != 211 || mcAntiPion2->GetPdgCode() != 211 || mcPion1->GetPdgCode() != -211) return -1; // check pdg

  Int_t IndexMother[2];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcAntiPion3->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is anti-lambda

  IndexMother[0] = mcMother->GetMother(); // mother of lambda
  IndexMother[1] = mcAntiPion2->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3312 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xi+

  IndexMother[0] = mcMother->GetMother(); // mother of Xi+
  IndexMother[1] = mcPion1->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -4132 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is anti-Xic0

//  if ( mcMother->IsPrimary() ) return 1;
//  if ( mcMother->IsPhysicalPrimary() ) return 2;
//  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
//  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
//  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

  Int_t CheckOrigin = AliVertexingHFUtils::CheckOrigin(mcArray,mcMother,kTRUE);
  return CheckOrigin;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::SelectTrack(AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks, Bool_t *seleFlags)
{
  // Select good tracks using fAnaCuts (AliRDHFCuts object)
  if(trkEntries==0) return;

  nSeleTrks=0;                                                                                                 
  for(Int_t i=0; i<trkEntries; i++) {
    seleFlags[i] = kFALSE;
    
    AliVTrack *track;
    track = (AliVTrack*)event->GetTrack(i);
    
//    if(track->GetID()<0) continue;
    Double_t covtest[21];
    if(!track->GetCovarianceXYZPxPyPz(covtest)) continue;
    
//    AliAODTrack *aodt = (AliAODTrack*)track;

/*
    if(!fAnaCuts) continue;
    if(fAnaCuts->SingleTrkCuts(aodt)){
      seleFlags[i]=kTRUE;
      nSeleTrks++;
//      fHistoPiPtRef->Fill(aodt->Pt());
    }
*/
  } // end loop on tracks
}
//_____________________________________________________________________________
KFParticle AliAnalysisTaskSEXicZero2XiPifromKFP::CreateKFTrack(Double_t *param, Double_t *cov, Float_t Chi2perNDF, Int_t charge, Int_t pdg)
{
  KFParticle::SetField(fBzkG);

  // Interface to KFParticle
  KFPTrack kfpTrk;
  // Set the values
  kfpTrk.SetParameters((Float_t) param[0],(Float_t) param[1],(Float_t) param[2],
                       (Float_t) param[3],(Float_t) param[4],(Float_t) param[5]);
  kfpTrk.SetCharge(charge);
  Float_t covF[21];
  for (Int_t i = 0; i<21;i++) { covF[i] = (Float_t) cov[i]; }
  kfpTrk.SetCovarianceMatrix(covF);
  kfpTrk.SetNDF(1);
  kfpTrk.SetChi2(Chi2perNDF);

  // Build KFParticle
  KFParticle kfp(kfpTrk, pdg);
  return kfp;
}
//_____________________________________________________________________________
KFVertex AliAnalysisTaskSEXicZero2XiPifromKFP::CreateKFVertex(Double_t *param, Double_t *cov)
{
  KFParticle::SetField(fBzkG);

  KFPVertex kfpVtx;
  // Set the values
  Float_t paramF[3] = {(Float_t) param[0],(Float_t) param[1],(Float_t) param[2]};
  kfpVtx.SetXYZ(paramF);
  Float_t covF[6] = {(Float_t) cov[0],(Float_t) cov[1],(Float_t) cov[2],
                     (Float_t) cov[3],(Float_t) cov[4],(Float_t) cov[5]};
  kfpVtx.SetCovarianceMatrix(covF);
  KFVertex kfv(kfpVtx);
  return kfv;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEXicZero2XiPifromKFP::CheckVertexCov(AliAODVertex *vtx)
{
  Double_t covMatrix[6];
  vtx->GetCovarianceMatrix(covMatrix);
  Double_t cov[3][3]={0.};
  cov[0][0] = covMatrix[0];
  cov[1][0] = covMatrix[1];
  cov[1][1] = covMatrix[2];
  cov[2][0] = covMatrix[3];
  cov[2][1] = covMatrix[4];
  cov[2][2] = covMatrix[5];
  if ( cov[0][0]<0 || cov[1][1]<0 || cov[2][2]<0 ) return kFALSE;
  for (Int_t i=0; i<3; i++) {
    for (Int_t j=0; j<3; j++) {
      if (i<=j) continue;
      if ( fabs(cov[i][j]) > TMath::Sqrt(cov[i][i]*cov[j][j]) ) return kFALSE;
    }
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEXicZero2XiPifromKFP::CheckTrackCov(AliAODTrack *track)
{
  Double_t covMatrix[21];
  track->GetCovarianceXYZPxPyPz(covMatrix);
  Double_t cov[6][6]={0.};
  cov[0][0] = covMatrix[0];
  cov[1][0] = covMatrix[1];
  cov[1][1] = covMatrix[2];
  cov[2][0] = covMatrix[3];
  cov[2][1] = covMatrix[4];
  cov[2][2] = covMatrix[5];
  cov[3][0] = covMatrix[6];
  cov[3][1] = covMatrix[7];
  cov[3][2] = covMatrix[8];
  cov[3][3] = covMatrix[9];
  cov[4][0] = covMatrix[10];
  cov[4][1] = covMatrix[11];
  cov[4][2] = covMatrix[12];
  cov[4][3] = covMatrix[13];
  cov[4][4] = covMatrix[14];
  cov[5][0] = covMatrix[15];
  cov[5][1] = covMatrix[16];
  cov[5][2] = covMatrix[17];
  cov[5][3] = covMatrix[18];
  cov[5][4] = covMatrix[19];
  cov[5][5] = covMatrix[20];
  if ( cov[0][0]<0 || cov[1][1]<0 || cov[2][2]<0 || cov[3][3]<0 || cov[4][4]<0 || cov[5][5]<0 ) return kFALSE;
  for (Int_t i=0; i<6; i++) {
    for (Int_t j=0; j<6; j++) {
      if (i<=j) continue;
      if ( fabs(cov[i][j]) > TMath::Sqrt(cov[i][i]*cov[j][j]) ) return kFALSE;
    }
  }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEXicZero2XiPifromKFP::CheckKFParticleCov(KFParticle kfp)
{
  if ( kfp.GetCovariance(0,0)<0 || kfp.GetCovariance(1,1)<0 || kfp.GetCovariance(2,2)<0 || kfp.GetCovariance(3,3)<0 || kfp.GetCovariance(4,4)<0 || kfp.GetCovariance(5,5)<0 ) return kFALSE;
  for (Int_t i=0; i<6; i++) {
    for (Int_t j=0; j<6; j++) {
      if (i<=j) continue;
      if ( fabs(kfp.GetCovariance(i,j)) > TMath::Sqrt(kfp.GetCovariance(i,i)*kfp.GetCovariance(j,j)) ) return kFALSE;
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
KFParticle AliAnalysisTaskSEXicZero2XiPifromKFP::CreateKFParticleFromAODtrack(AliAODTrack *track, Int_t pdg)
{
  KFParticle::SetField(fBzkG);

  Double_t trackParam[6];
  Double_t covMatrix[21];

  Bool_t IsDCA = track->GetXYZ(trackParam);
//  if (IsDCA) cout << "track position is at DCA" << endl;
//  if (!IsDCA) cout << "track position is at first point" << endl;
//  track->XvYvZv(trackParam);
  track->GetPxPyPz(&trackParam[3]);
  track->GetCovarianceXYZPxPyPz(covMatrix);
  
  KFPTrack kfpt;
  kfpt.SetParameters(trackParam);
  kfpt.SetCovarianceMatrix(covMatrix);
  kfpt.SetCharge(track->Charge());
  kfpt.SetNDF(1);
  kfpt.SetChi2(track->Chi2perNDF());

  KFParticle kfp(kfpt, pdg);

  return kfp;
}
//_____________________________________________________________________________
KFParticle AliAnalysisTaskSEXicZero2XiPifromKFP::CreateKFMotherParticle(AliAODTrack *track1, AliAODTrack *track2, Int_t pdg1, Int_t pdg2)
{
  KFParticle::SetField(fBzkG);

  Double_t trackParam[6];
  Double_t covMatrix[21];

  track1->GetXYZ(trackParam);
  track1->GetPxPyPz(&trackParam[3]);
  track1->GetCovarianceXYZPxPyPz(covMatrix);

  KFPTrack kfpt1;
  kfpt1.SetParameters(trackParam);
  kfpt1.SetCovarianceMatrix(covMatrix);
  kfpt1.SetCharge(track1->Charge());
  kfpt1.SetNDF(1);
  kfpt1.SetChi2(track1->Chi2perNDF());

  track2->GetXYZ(trackParam);
  track2->GetPxPyPz(&trackParam[3]);
  track2->GetCovarianceXYZPxPyPz(covMatrix);

  KFPTrack kfpt2;
  kfpt2.SetParameters(trackParam);
  kfpt2.SetCovarianceMatrix(covMatrix);
  kfpt2.SetCharge(track2->Charge());
  kfpt2.SetNDF(1);
  kfpt2.SetChi2(track2->Chi2perNDF());

  // now we have all info to create the KFParticle version of the daughters
  KFParticle kfpDaughter1(kfpt1, pdg1);
  KFParticle kfpDaughter2(kfpt2, pdg2);
                       
  KFParticle kfpMother(kfpDaughter1, kfpDaughter2);

  return kfpMother;
}
//_____________________________________________________________________________
KFParticle AliAnalysisTaskSEXicZero2XiPifromKFP::CreateSecKFParticle(KFParticle kfp1, AliAODTrack *track2, Int_t pdg1, Int_t pdg2)
{
  KFParticle::SetField(fBzkG);

  Double_t trackParam[6];
  Double_t covMatrix[21];

  track2->GetXYZ(trackParam);
//  track2->XvYvZv(trackParam);
  track2->GetPxPyPz(&trackParam[3]);
  track2->GetCovarianceXYZPxPyPz(covMatrix);

  KFPTrack kfpt2;
  kfpt2.SetParameters(trackParam);
  kfpt2.SetCovarianceMatrix(covMatrix);
  kfpt2.SetCharge(track2->Charge());
  kfpt2.SetNDF(1);
  kfpt2.SetChi2(track2->Chi2perNDF());

  KFParticle kfpDaughter2(kfpt2, pdg2);

  kfp1.SetPDG(pdg1);
//  KFParticle kfpMother(kfp1, kfpDaughter2);
  KFParticle kfpMother;
//  const KFParticle *vDaughters[2] = {&kfp1, &kfpDaughter2};
  const KFParticle *vDaughters[2] = {&kfpDaughter2, &kfp1}; // the order is important
  kfpMother.Construct(vDaughters, 2);

  return kfpMother;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskSEXicZero2XiPifromKFP::CosPointingAngleKF(KFParticle kfp, KFParticle kfpmother)
{
  Double_t v[3];
  v[0] = kfp.GetX() - kfpmother.GetX();
  v[1] = kfp.GetY() - kfpmother.GetY();
  v[2] = kfp.GetZ() - kfpmother.GetZ();

  Double_t p[3];
  p[0] = kfp.GetPx();
  p[1] = kfp.GetPy();
  p[2] = kfp.GetPz();

  Double_t ptimesv2 = (p[0]*p[0]+p[1]*p[1]+p[2]*p[2])*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  if ( ptimesv2<=0 ) return 0.0;
  else {
    Double_t cos = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2]) / TMath::Sqrt(ptimesv2);
    if(cos >  1.0) cos =  1.0;
    if(cos < -1.0) cos = -1.0;
    return cos;
  }
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskSEXicZero2XiPifromKFP::CosThetaStarKF(Int_t ip, UInt_t pdgvtx, UInt_t pdgprong0, UInt_t pdgprong1, KFParticle kfpvtx, KFParticle kfpprong0, KFParticle kfpprong1)
{
  Double_t massvtx = TDatabasePDG::Instance()->GetParticle(pdgvtx)->Mass();
  Double_t massp[2];

  massp[0] = TDatabasePDG::Instance()->GetParticle(pdgprong0)->Mass();
  massp[1] = TDatabasePDG::Instance()->GetParticle(pdgprong1)->Mass();

  Double_t pStar = TMath::Sqrt((massvtx*massvtx-massp[0]*massp[0]-massp[1]*massp[1])*(massvtx*massvtx-massp[0]*massp[0]-massp[1]*massp[1])-4.*massp[0]*massp[0]*massp[1]*massp[1])/(2.*massvtx);

  Double_t e = kfpvtx.GetE();
  Double_t beta = kfpvtx.GetP()/e;
  Double_t gamma = e/massvtx;

  TVector3 mom;
  TVector3 momTot(kfpvtx.GetPx(), kfpvtx.GetPy(), kfpvtx.GetPz());

  if (ip==0) {
    mom.SetXYZ(kfpprong0.GetPx(), kfpprong0.GetPy(), kfpprong0.GetPz());
  }
  if (ip==1) {
    mom.SetXYZ(kfpprong1.GetPx(), kfpprong1.GetPy(), kfpprong1.GetPz());
  }

  Double_t cts = ( (mom.Dot(momTot)/momTot.Mag()) /gamma-beta*TMath::Sqrt(pStar*pStar+massp[ip]*massp[ip]) ) / pStar;

  return cts;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicZero2XiPifromKFP::MatchToMCXiMinus(AliAODTrack *trackProton, AliAODTrack *trackPion3, AliAODTrack *trackPion2, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
  Int_t labelPion3  = fabs(trackPion3->GetLabel());
  if (labelPion3<0) return -1;
  AliAODMCParticle* mcPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion3));
  Int_t labelPion2  = fabs(trackPion2->GetLabel());
  if (labelPion2<0) return -1;
  AliAODMCParticle* mcPion2 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion2));

  if ( mcProton->GetPdgCode() != 2212 || mcPion3->GetPdgCode() != -211 || mcPion2->GetPdgCode() != -211 ) return -1; // check pdg

  Int_t IndexMother[4];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPion3->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is lambda

  IndexMother[2] = mcMother->GetMother(); // mother of lambda
  IndexMother[3] = mcPion2->GetMother();
  if ( IndexMother[2]<0 || IndexMother[3]<0 ) return -1; // check mother exist
  if ( IndexMother[2] != IndexMother[3] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[2]));
  if ( mcMother->GetPdgCode() != 3312 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xi-

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[2];
  return 1;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicZero2XiPifromKFP::MatchToMCXiPlus(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, AliAODTrack *trackAntiPion2, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
  Int_t labelAntiPion3  = fabs(trackAntiPion3->GetLabel());
  if (labelAntiPion3<0) return -1;
  AliAODMCParticle* mcAntiPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion3));
  Int_t labelAntiPion2  = fabs(trackAntiPion2->GetLabel());
  if (labelAntiPion2<0) return -1;
  AliAODMCParticle* mcAntiPion2 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion2));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcAntiPion3->GetPdgCode() != 211 || mcAntiPion2->GetPdgCode() != 211 ) return -1; // check pdg

  Int_t IndexMother[4];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcAntiPion3->GetMother();
  if ( IndexMother[0]<0 || IndexMother[1]<0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Anti-lambda

  IndexMother[2] = mcMother->GetMother(); // mother of lambda
  IndexMother[3] = mcAntiPion2->GetMother();
  if ( IndexMother[2]<0 || IndexMother[3]<0 ) return -1; // check mother exist
  if ( IndexMother[2] != IndexMother[3] ) return -1; // check the same mother
  mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[2]));
  if ( mcMother->GetPdgCode() != -3312 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Xi+

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[2];
  return 1;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicZero2XiPifromKFP::MatchToMCLambda(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<=0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
  Int_t labelPion3  = fabs(trackPion3->GetLabel());
  if (labelPion3<=0) return -1;
  AliAODMCParticle* mcPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion3));

  if ( mcProton->GetPdgCode() != 2212 || mcPion3->GetPdgCode() != -211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is lambda and only have two daughters

//  AliAODMCParticle* mcMother2 = static_cast<AliAODMCParticle*>(mcArray->At(fabs(mcMother->GetMother())));

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[0];
  return 1;

}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicZero2XiPifromKFP::MatchToMCLambdaFromXi(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelProton = fabs(trackProton->GetLabel());
  if (labelProton<=0) return -1;
  AliAODMCParticle* mcProton = static_cast<AliAODMCParticle*>(mcArray->At(labelProton));
  Int_t labelPion3  = fabs(trackPion3->GetLabel());
  if (labelPion3<=0) return -1;
  AliAODMCParticle* mcPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion3));

  if ( mcProton->GetPdgCode() != 2212 || mcPion3->GetPdgCode() != -211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcProton->GetMother();
  IndexMother[1] = mcPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != 3122 ) return -1; // check mother is lambda
  
  Int_t Index_Mother_Xi = mcMother->GetMother();
  if ( Index_Mother_Xi<=0 ) return -1;
  AliAODMCParticle* mcMother_Xi = static_cast<AliAODMCParticle*>(mcArray->At(Index_Mother_Xi));
  if ( mcMother_Xi->GetPdgCode() != 3312 ) return -1;

//  if ( !mcMother_Xi->IsPhysicalPrimary() ) return -1; // check IsPhysicalPrimary()

  return IndexMother[0];

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicZero2XiPifromKFP::MatchToMCAntiLambda(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<=0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
  Int_t labelAntiPion3  = fabs(trackAntiPion3->GetLabel());
  if (labelAntiPion3<=0) return -1;
  AliAODMCParticle* mcAntiPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion3));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcAntiPion3->GetPdgCode() != 211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcAntiPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 || mcMother->GetNDaughters()!=2 ) return -1; // check mother is Anti-lambda and only have two daughters

//  if ( mcMother->IsPrimary() ) return 1;
  if ( mcMother->IsPhysicalPrimary() ) return 2;
  if ( mcMother->IsSecondaryFromWeakDecay() ) return 3;
  if ( mcMother->IsSecondaryFromMaterial() ) return 4;
  if ( mcMother->IsFromSubsidiaryEvent() ) return 5;

//  return IndexMother[0];
  return 1;

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicZero2XiPifromKFP::MatchToMCAntiLambdaFromXi(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray)
{
  // Check if all of the tracks is matched to a MC signal
  // If no, return -1;
  // If yes, return label (>=0) of the AliAODMCParticle

  Int_t labelAntiProton = fabs(trackAntiProton->GetLabel());
  if (labelAntiProton<=0) return -1;
  AliAODMCParticle* mcAntiProton = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiProton));
  Int_t labelAntiPion3  = fabs(trackAntiPion3->GetLabel());
  if (labelAntiPion3<=0) return -1;
  AliAODMCParticle* mcAntiPion3 = static_cast<AliAODMCParticle*>(mcArray->At(labelAntiPion3));

  if ( mcAntiProton->GetPdgCode() != -2212 || mcAntiPion3->GetPdgCode() != 211 ) return -1; // check PDG

  Int_t IndexMother[2];
  IndexMother[0] = mcAntiProton->GetMother();
  IndexMother[1] = mcAntiPion3->GetMother();
  if ( IndexMother[0]<=0 || IndexMother[1]<=0 ) return -1; // check mother exist
  if ( IndexMother[0] != IndexMother[1] ) return -1; // check the same mother
  AliAODMCParticle* mcMother = static_cast<AliAODMCParticle*>(mcArray->At(IndexMother[0]));
  if ( mcMother->GetPdgCode() != -3122 ) return -1; // check mother is Anti-lambda

  Int_t Index_Mother_Xi = mcMother->GetMother();
  if ( Index_Mother_Xi<=0 ) return -1;
  AliAODMCParticle* mcMother_Xi = static_cast<AliAODMCParticle*>(mcArray->At(Index_Mother_Xi));
  if ( mcMother_Xi->GetPdgCode() != -3312 ) return -1;

//  if ( !mcMother_Xi->IsPhysicalPrimary() ) return -1; // check IsPhysicalPrimary()

  return IndexMother[0];

}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicZero2XiPifromKFP::MatchToMCPion(AliAODTrack *track, TClonesArray *mcArray)
{
  Int_t labelPion = fabs(track->GetLabel());
  if (labelPion<=0) return -1;
  AliAODMCParticle* mcPion = static_cast<AliAODMCParticle*>(mcArray->At(labelPion));
  if ( TMath::Abs(mcPion->GetPdgCode()) != 211 ) return -1;

  return labelPion;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskSEXicZero2XiPifromKFP::MatchToXicZeroMC(TClonesArray *mcArray, Int_t PDGXicZero, const Int_t nDaughters, const Int_t *daughterIndex, const Int_t *daughterPDG)
{
  ///
  /// Check if this candidate is matched to a MC signal XicZero
  /// If yes, return Index (>=0) of the AliAODMCParticle
  /// If no, return -1
  ///

  Int_t IndexMom[10] = {0};
  Int_t IndexMother=-1;
  Bool_t pdgUsed[10] = {0};

  // loop on daughter Index
  for(Int_t i=0; i<nDaughters; i++) {
    IndexMom[i]=-1;
    Int_t Index = daughterIndex[i];
    if(Index<0) {
      printf("daughter with negative index %d\n", Index);
      return -1;
    }
    AliAODMCParticle *part = (AliAODMCParticle*)mcArray->At(Index);
    if(!part) { 
      printf("no MC particle\n");
      return -1;
    }

    // check the PDG of the daughter
    Int_t pdgPart = part->GetPdgCode();
    for(Int_t j=0; j<nDaughters; j++) {
      if(!pdgUsed[j] && pdgPart==daughterPDG[j]) {
        pdgUsed[j]=kTRUE;
        break;
      }
    }

    AliAODMCParticle *mother = part;
    while ( mother->GetMother()>=0 ) {
      IndexMother = mother->GetMother();
      mother = (AliAODMCParticle*)mcArray->At(IndexMother);
      if (!mother) {
        printf("no MC mother particle\n");
        break;
      }
      Int_t pdgMother = mother->GetPdgCode();
      if ( pdgMother==PDGXicZero ) { // check mother is XicZero
        IndexMom[i]=IndexMother;
        break;
      } else if( pdgMother>PDGXicZero || pdgMother<10 ) {
        break;
      }
    }

    if( IndexMom[i]==-1 ) return -1; // mother PDG not ok for this daughter

  } // end loop on daughters

  IndexMother=IndexMom[0];
  for(Int_t i=0; i<nDaughters; i++) {
    // all Index have to be the same and !=-1
    if(IndexMom[i]==-1)          return -1;
    if(IndexMom[i]!=IndexMother) return -1;
    // check that all daughter PDGs are matched
    if(pdgUsed[i]==kFALSE)       return -1;
  }

  return IndexMother;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::DefineEvent()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(3)->GetContainer()->GetName();
  fTree_Event = new TTree(nameoutput, "Event");
  Int_t nVar = 7;
  fVar_Event = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  fVarNames[0]  = "centrality";
  fVarNames[1]  = "z_vtx_reco";
  fVarNames[2]  = "n_vtx_contributors";
  fVarNames[3]  = "n_tracks";
  fVarNames[4]  = "is_ev_rej";
  fVarNames[5]  = "run_number";
  fVarNames[6]  = "ev_id";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_Event->Branch(fVarNames[ivar].Data(), &fVar_Event[ivar], Form("%s/f", fVarNames[ivar].Data()));
  }

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::DefineTreeRecXic0()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  fTree_Xic0 = new TTree(nameoutput, "Xic0 variables tree");
  Int_t nVar = 37;
  fVar_Xic0 = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  fVarNames[0]  = "nSigmaTPC_PiFromXic0"; // TPC nsigma for pion coming from Xic0
  fVarNames[1]  = "nSigmaTOF_PiFromXic0"; // TOF nsigma for pion coming from Xic0
  fVarNames[2]  = "nSigmaTPC_PiFromXi"; // TPC nsigma for pion coming from Xi
  fVarNames[3]  = "nSigmaTOF_PiFromXi"; // TOF nsigma for pion coming from Xi
  fVarNames[4]  = "nSigmaTPC_PiFromLam"; // TPC nsigma for pion coming from Lambda
  fVarNames[5]  = "nSigmaTPC_PrFromLam"; // TPC nsigma for proton coming from Lambda

  fVarNames[6]  = "DCA_LamDau"; // Distance between proton and pion from Lambda decays (calculated from AOD v0)

  fVarNames[7]  = "chi2geo_Lam"; // chi2_geometry of Lambda
  fVarNames[8]  = "ldl_Lam"; // l/dl of Lambda
  fVarNames[9]  = "chi2topo_LamToPV"; // chi2_topo of Lambda to PV
  fVarNames[10] = "chi2geo_Xi"; // chi2_geometry of Xi (with Lambda mass const.)
  fVarNames[11] = "ldl_Xi"; // l/dl of Xi (with Lambda mass const.)
  fVarNames[12] = "chi2topo_XiToPV"; // chi2_topo of Xi to PV

  fVarNames[13] = "DecayLxy_Lam"; // decay length of Lambda in x-y plane
  fVarNames[14] = "ct_Lam"; // life time of Lambda
  fVarNames[15] = "DecayLxy_Xi"; // decay length of Xi in x-y plane
  fVarNames[16] = "ct_Xi"; // life time of Xi

  fVarNames[17] = "PA_Lam"; // pointing angle of Lmabda (pointing back to Xi)
  fVarNames[18] = "PA_LamToPV"; // pointing angle of Lambda (pointing back to PV)
  fVarNames[19] = "PA_Xi"; // pointing angle of Xi (pointing back to PV)
  fVarNames[20] = "mass_Lam"; // mass of Lambda (without mass const.)
  fVarNames[21] = "mass_Xi"; // mass of Xi (without mass const.)
  fVarNames[22] = "pt_PiFromXic0"; // pt of pion coming from Xic0
  fVarNames[23] = "pt_Xic0"; // pt of Xic0
  fVarNames[24] = "rap_Xic0"; // rapidity of Xic0
  fVarNames[25] = "mass_Xic0"; // mass of Xic0
  fVarNames[26] = "CosThetaStar_PiFromXic0"; // CosThetaStar of pion coming from Xic0
  fVarNames[27] = "CosThetaStar_Xi"; // CosThetaStar of Xi coming from Xic0
  fVarNames[28] = "chi2prim_PiFromXic0"; // DCA of pion coming from Xic0 in x-y plane
  fVarNames[29] = "DCAxy_PiFromXic0"; // DCA of pion coming from Xic0 in x-y plane
  fVarNames[30] = "Source_Xic0"; // flag for Xic0 MC truth (“4” prompt, "5" feed-down, “<0” background)
  fVarNames[31] = "mass_K0S"; // mass of Ks0
  fVarNames[32] = "mass_Gamma"; // mass of e+e-
  fVarNames[33] = "chi2topo_LamToXi"; // chi2_topo of Lambda to Xi
  fVarNames[34] = "chi2topo_XiToXic0"; // chi2_topo of Xi to Xic0
  fVarNames[35] = "DecayLxy_Xic0"; // decay length of Xic0 in x-y plane
  fVarNames[36] = "ct_Xic0"; // life time of Xic0

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_Xic0->Branch(fVarNames[ivar].Data(), &fVar_Xic0[ivar], Form("%s/f", fVarNames[ivar].Data()));
  }

  return;
}

/*
//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::DefineVarTreePiPlus()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(6)->GetContainer()->GetName();
  fVarTree_PiPlus = new TTree(nameoutput, "PiPlus variables tree");
  Int_t nVar = 5;
  fVar_PiPlus = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  fVarNames[0]  = "charge_PiPlus";
  fVarNames[1]  = "pt_PiPlus";
  fVarNames[2]  = "pz_PiPlus";
  fVarNames[3]  = "eta_PiPlus";
  fVarNames[4]  = "MatchToMC";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVarTree_PiPlus->Branch(fVarNames[ivar].Data(), &fVar_PiPlus[ivar], Form("%s/f", fVarNames[ivar].Data()));
  }

  return;
}
*/

/*
//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::DefineVarTreePiMinus()
{
  // This is to define tree variables

  const char* nameoutput = GetOutputSlot(9)->GetContainer()->GetName();
  fVarTree_PiMinus = new TTree(nameoutput, "PiMinus variables tree");
  Int_t nVar = 5;
  fVar_PiMinus = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];

  fVarNames[0]  = "charge_PiMinus";
  fVarNames[1]  = "pt_PiMinus";
  fVarNames[2]  = "pz_PiMinus";
  fVarNames[3]  = "eta_PiMinus";
  fVarNames[4]  = "MatchToMC";

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fVarTree_PiMinus->Branch(fVarNames[ivar].Data(), &fVar_PiMinus[ivar], Form("%s/f", fVarNames[ivar].Data()));
  }

  return;
}
*/

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::DefineTreeGenXic0()
{
  const char* nameoutput = GetOutputSlot(5)->GetContainer()->GetName();
  fTree_Xic0MCGen = new TTree(nameoutput,"Xic0 MC variables tree");
  Int_t nVar = 3;
  fVar_Xic0MCGen = new Float_t[nVar];
  TString *fVarNames = new TString[nVar];
  
  fVarNames[0] = "XicY";
  fVarNames[1] = "XicPt";
  fVarNames[2] = "Source_Xic0";
  /*
  fVarNames[ 0]="Centrality";
  fVarNames[ 1]="DecayType";
  fVarNames[ 2]="XicSource";
  fVarNames[ 3]="XicEta";
  fVarNames[ 4]="XicY";
  fVarNames[ 5]="XicPx";
  fVarNames[ 6]="XicPy";
  fVarNames[ 7]="XicPz";
  fVarNames[ 8]="PiPx";
  fVarNames[ 9]="PiPy";
  fVarNames[10]="PiPz";
  fVarNames[11]="CascPx";
  fVarNames[12]="CascPy";
  fVarNames[13]="CascPz";
  fVarNames[14]="XicPdgCode";
  fVarNames[15]="PiPdgCode";
  fVarNames[16]="CascPdgCode";
  fVarNames[17]="RunNumber";
  fVarNames[18]="EvNumber";
  fVarNames[19]="IsPrompt";
  */

  for (Int_t ivar=0; ivar<nVar; ivar++) {
    fTree_Xic0MCGen->Branch(fVarNames[ivar].Data(),&fVar_Xic0MCGen[ivar],Form("%s/f",fVarNames[ivar].Data()));
  }

  return;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskSEXicZero2XiPifromKFP::InvMassV0atPV(AliAODTrack *trk1, AliAODTrack *trk2, Int_t pdg1, Int_t pdg2)
{
  
  Double_t mass1 = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
  Double_t mass2 = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();
  Double_t E1 = TMath::Sqrt(mass1*mass1 + trk1->P()*trk1->P());
  Double_t E2 = TMath::Sqrt(mass2*mass2 + trk2->P()*trk2->P());
  Double_t mass = TMath::Sqrt( (E1+E2)*(E1+E2) - (trk1->Px()+trk2->Px())*(trk1->Px()+trk2->Px()) - (trk1->Py()+trk2->Py())*(trk1->Py()+trk2->Py()) - (trk1->Pz()+trk2->Pz())*(trk1->Pz()+trk2->Pz()) );

  return mass;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::DefineAnaHist()
{
  // Define analysis histograms
  
  Int_t bins_xicmcgen[3]    = {20, 20, 10};
  Double_t xmin_xicmcgen[3] = {0., -1., 0.};
  Double_t xmax_xicmcgen[3] = {20., 1., 100.};
  fHistMCGen_XicZeroTot = new THnSparseF("fHistMCGen_XicZeroTot","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fHistMCGen_XicZero = new THnSparseF("fHistMCGen_XicZero","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fHistMCGen_AntiXicZero = new THnSparseF("fHistMCGen_AntiXicZero","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fHistMCGen_PionTot = new THnSparseF("fHistMCGen_PionTot","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fHistMCGen_PionPlus = new THnSparseF("fHistMCGen_PionPlus","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fHistMCGen_PionMinus = new THnSparseF("fHistMCGen_PionMinus","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fHistMCGen_XiTot = new THnSparseF("fHistMCGen_XiTot","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fHistMCGen_XiMinus = new THnSparseF("fHistMCGen_XiMinus","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fHistMCGen_XiPlus = new THnSparseF("fHistMCGen_XiPlus","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fHistMCGen_Lambda = new THnSparseF("fHistMCGen_Lambda","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);
  fHistMCGen_AntiLambda = new THnSparseF("fHistMCGen_AntiLambda","",3,bins_xicmcgen,xmin_xicmcgen,xmax_xicmcgen);

  Int_t bins_masspt[3] = {300, 20, 10};
  Double_t xmin_masspt[3] = {1., 0.,  0.};
  Double_t xmax_masspt[3] = {4., 20., 100.};
  fHistMCGen_PiXiInvMass = new THnSparseF("fHistMCGen_PiXiInvMass","",3,bins_masspt,xmin_masspt,xmax_masspt);
  fHistMCGen_PiXiMassvsPiPt = new THnSparseF("fHistMCGen_PiXiMassvsPiPt","",3,bins_masspt,xmin_masspt,xmax_masspt);
  fHistMCGen_PiXiMassvsPiPt_PionPlus = new THnSparseF("fHistMCGen_PiXiMassvsPiPt_PionPlus","",3,bins_masspt,xmin_masspt,xmax_masspt);
  fHistMCGen_PiXiMassvsPiPt_PionMinus = new THnSparseF("fHistMCGen_PiXiMassvsPiPt_PionMinus","",3,bins_masspt,xmin_masspt,xmax_masspt);

  fHistMCXicZeroDecayType = new TH1F("fHistMCXicZeroDecayType","",4,0.5,4.5);
  fHistMCXiDecayType = new TH1F("fHistMCXiDecayType","",4,0.5,4.5);

  fHistMCpdg_All = new TH1F("fHistMCpdg_All", "PDG", 20000, -10000, 10000);
  fHistMCpdg_Dau_XicZero = new TH1F("fHistMCpdg_Dau_XicZero", "PDG of daughters of Xic0 & anti-Xic0", 20000, -10000, 10000);
  fHistMCpdg_Dau_XicPM   = new TH1F("fHistMCpdg_Dau_XicPM", "PDG of daughters of Xic+ & Xic-", 20000, -10000, 10000);

  fOutputList->Add(fHistMCGen_XicZeroTot);
  fOutputList->Add(fHistMCGen_XicZero);
  fOutputList->Add(fHistMCGen_AntiXicZero);
  fOutputList->Add(fHistMCGen_PionTot);
  fOutputList->Add(fHistMCGen_PionPlus);
  fOutputList->Add(fHistMCGen_PionMinus);
  fOutputList->Add(fHistMCGen_XiTot);
  fOutputList->Add(fHistMCGen_XiMinus);
  fOutputList->Add(fHistMCGen_XiPlus);
  fOutputList->Add(fHistMCGen_Lambda);
  fOutputList->Add(fHistMCGen_AntiLambda);
  fOutputList->Add(fHistMCGen_PiXiInvMass);
  fOutputList->Add(fHistMCGen_PiXiMassvsPiPt);
  fOutputList->Add(fHistMCGen_PiXiMassvsPiPt_PionPlus);
  fOutputList->Add(fHistMCGen_PiXiMassvsPiPt_PionMinus);

  fOutputList->Add(fHistMCXicZeroDecayType);
  fOutputList->Add(fHistMCXiDecayType);

  fOutputList->Add(fHistMCpdg_All);
  fOutputList->Add(fHistMCpdg_Dau_XicZero);
  fOutputList->Add(fHistMCpdg_Dau_XicPM);

}

/*
//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::FillPiPlusROOTObjects(AliAODTrack *track, TClonesArray *mcArray)
{
  for (Int_t i=0; i<5; i++) {
    fVar_PiPlus[i] = -9999.;
  }

  fVar_PiPlus[0] = track->Charge();
  fVar_PiPlus[1] = track->Pt();
  fVar_PiPlus[2] = track->Pz();
  fVar_PiPlus[3] = track->Eta();

  if (fIsMC) {
    Int_t lab_Pion = MatchToMCPion(track, mcArray);
    if (lab_Pion>-1) fVar_PiPlus[4] = 1;
  }

  fVarTree_PiPlus->Fill();

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::FillPiMinusROOTObjects(AliAODTrack *track, TClonesArray *mcArray)
{
  for (Int_t i=0; i<5; i++) {
    fVar_PiMinus[i] = -9999.;
  }

  fVar_PiMinus[0] = track->Charge();
  fVar_PiMinus[1] = track->Pt();
  fVar_PiMinus[2] = track->Pz();
  fVar_PiMinus[3] = track->Eta();

  if (fIsMC) {
    Int_t lab_Pion = MatchToMCPion(track, mcArray);
    if (lab_Pion>-1) fVar_PiMinus[4] = 1;
  }

  fVarTree_PiMinus->Fill();

  return;

}
*/

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::FillEventROOTObjects()
{

  for (Int_t i=0; i<7; i++) {
    fVar_Event[i] = 0.;
  }

  Double_t pos[3];
  fpVtx->GetXYZ(pos);

  fVar_Event[1] = pos[2];

  fTree_Event->Fill();

  return;

}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::FillTreeRecXic0FromV0(KFParticle kfpXic0, AliAODTrack *trackPiFromXic0, KFParticle kfpBP, KFParticle kfpXiMinus, KFParticle kfpXiMinus_m, AliAODTrack *trackPiFromXi, AliAODv0 *v0, KFParticle kfpK0Short, KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkProton, AliAODTrack *trkPion, KFParticle PV, TClonesArray *mcArray, Int_t lab_Xic0)
{
  
  for (Int_t i=0; i<32; i++) {
    fVar_Xic0[i] = -9999.;
  }

//  Float_t nSigmaTOF_PiFromXic0 = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trackPiFromXic0,AliPID::kPion);
//  Float_t nSigmaTOF_PiFromXi   = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trackPiFromXi,AliPID::kPion);

//  Float_t nSigmaTPC_PiFromXic0 = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trackPiFromXic0,AliPID::kPion);
//  Float_t nSigmaTPC_PiFromXi   = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trackPiFromXi,AliPID::kPion);
//  Float_t nSigmaTPC_PrFromLam  = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkProton,AliPID::kProton);
//  Float_t nSigmaTPC_PiFromLam  = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkPion,AliPID::kPion);


  Float_t nSigmaTOF_PiFromXic0 = fPID->NumberOfSigmasTOF(trackPiFromXic0,AliPID::kPion);
  Float_t nSigmaTOF_PiFromXi   = fPID->NumberOfSigmasTOF(trackPiFromXi,AliPID::kPion);
  Float_t nSigmaTPC_PiFromXic0 = fPID->NumberOfSigmasTPC(trackPiFromXic0,AliPID::kPion);
  Float_t nSigmaTPC_PiFromXi   = fPID->NumberOfSigmasTPC(trackPiFromXi,AliPID::kPion);
  Float_t nSigmaTPC_PrFromLam  = fPID->NumberOfSigmasTPC(trkProton,AliPID::kProton);
  Float_t nSigmaTPC_PiFromLam  = fPID->NumberOfSigmasTPC(trkPion,AliPID::kPion);
  
  if ( fabs(nSigmaTPC_PiFromXic0)>=4. || fabs(nSigmaTPC_PiFromXi)>=4. || fabs(nSigmaTPC_PrFromLam)>=4. || fabs(nSigmaTPC_PiFromLam)>=4. ) return;

//  AliAODTrack *trk0 = (AliAODTrack*) (v0->GetDaughter(0));

//  Double_t alpha_FirstDaugPos = (v0->MomPosAlongV0() - v0->MomNegAlongV0())/(v0->MomPosAlongV0() + v0->MomNegAlongV0());

//  if (trk0->Charge()<0) {
//    alpha_FirstDaugPos = (v0->MomNegAlongV0() - v0->MomPosAlongV0())/(v0->MomPosAlongV0() + v0->MomNegAlongV0());
//  }

  KFParticle kfpLambda_pv = kfpLambda;
  kfpLambda_pv.SetProductionVertex(PV);
//  KFParticle kfpXiMinus_pv = kfpXiMinus;
//  kfpXiMinus_pv.SetProductionVertex(PV);
  KFParticle kfpXic0_pv = kfpXic0;
  kfpXic0_pv.SetProductionVertex(PV);

  KFParticle kfpXiMinus_Xic0 = kfpXiMinus;
  kfpXiMinus_Xic0.SetProductionVertex(kfpXic0_pv);

  KFParticle kfpBP_Xic0 = kfpBP;
  kfpBP_Xic0.SetProductionVertex(kfpXic0);

  //************************** calculate l/Δl for Lambda *************************************
  Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
  Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
  Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
  Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
  Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
  if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
  dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
  if ( dl_Lambda<=0 ) return;
  //***************************************************************************************
  //************************** calculate l/Δl for Xi- *************************************
  Double_t dx_Xi = PV.GetX()-kfpXiMinus.GetX();
  Double_t dy_Xi = PV.GetY()-kfpXiMinus.GetY();
  Double_t dz_Xi = PV.GetZ()-kfpXiMinus.GetZ();
  Double_t l_Xi = TMath::Sqrt(dx_Xi*dx_Xi + dy_Xi*dy_Xi + dz_Xi*dz_Xi);
  Double_t dl_Xi = (PV.GetCovariance(0)+kfpXiMinus.GetCovariance(0))*dx_Xi*dx_Xi + (PV.GetCovariance(2)+kfpXiMinus.GetCovariance(2))*dy_Xi*dy_Xi + (PV.GetCovariance(5)+kfpXiMinus.GetCovariance(5))*dz_Xi*dz_Xi + 2*( (PV.GetCovariance(1)+kfpXiMinus.GetCovariance(1))*dx_Xi*dy_Xi + (PV.GetCovariance(3)+kfpXiMinus.GetCovariance(3))*dx_Xi*dz_Xi + (PV.GetCovariance(4)+kfpXiMinus.GetCovariance(4))*dy_Xi*dz_Xi );
  if ( fabs(l_Xi)<1.e-8f ) l_Xi = 1.e-8f;
  dl_Xi = dl_Xi<0. ? 1.e8f : sqrt(dl_Xi)/l_Xi;
  if ( dl_Xi<=0 ) return;
  //***************************************************************************************

  // calculate CosPointingAngle
  Double_t cosPA_v0      = CosPointingAngleKF(kfpLambda_m, kfpXiMinus);
  KFParticle kfpLambda_Xi = kfpLambda_m;
  kfpLambda_Xi.SetProductionVertex(kfpXiMinus);

  Double_t cosPA_XiToPV  = CosPointingAngleKF(kfpXiMinus_m, PV);
  Double_t cosPA_XiToXic = CosPointingAngleKF(kfpXiMinus_m, kfpXic0);
  KFParticle kfpXiMinus_Xic = kfpXiMinus_m;
  kfpXiMinus_Xic.SetProductionVertex(kfpXic0);

//  if ( kfpLambda_Xi.GetChi2()/kfpLambda_Xi.GetNDF() >= fAnaCuts->GetKFPLam_Chi2topoMax() ) return;
  if ( kfpXiMinus_Xic.GetChi2()/kfpXiMinus_Xic.GetNDF() >= fAnaCuts->GetKFPXi_Chi2topoMax() ) return;
  if ( l_Xi/dl_Xi <= fAnaCuts->GetKFPXi_lDeltalMin() ) return;

  const Float_t PDGmassXic0 = TDatabasePDG::Instance()->GetParticle(4132)->Mass();
  Float_t mass_Xic_PV, err_mass_Xic_PV;
  kfpXic0_pv.GetMass(mass_Xic_PV, err_mass_Xic_PV);
  fVar_Xic0[29] = mass_Xic_PV;

  if ( fabs(mass_Xic_PV-PDGmassXic0) > fAnaCuts->GetProdMassTolXic0() ) return;

  fVar_Xic0[0]  = nSigmaTPC_PiFromXic0;
  fVar_Xic0[1]  = nSigmaTOF_PiFromXic0;
  fVar_Xic0[2]  = nSigmaTPC_PiFromXi;
  fVar_Xic0[3]  = nSigmaTOF_PiFromXi;
  fVar_Xic0[4]  = nSigmaTPC_PiFromLam; 
  fVar_Xic0[5]  = nSigmaTPC_PrFromLam;
  fVar_Xic0[6]  = v0->DcaV0Daughters(); // DCA_LamDau
  fVar_Xic0[7]  = kfpLambda.GetChi2()/kfpLambda.GetNDF(); // chi2geo_Lam
  fVar_Xic0[8]  = kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF();
  fVar_Xic0[9]  = kfpLambda_Xi.GetChi2()/kfpLambda_Xi.GetNDF();
  fVar_Xic0[10] = kfpXiMinus_Xic.GetChi2()/kfpXiMinus_Xic.GetNDF();

  fVar_Xic0[11] = l_Lambda/dl_Lambda;
  fVar_Xic0[12] = l_Xi/dl_Xi;

  Float_t DecayL_Lam, err_DecayL_Lam;
  Float_t DecayL_Xi, err_DecayL_Xi;
  Float_t DecayL_XicZero, err_DecayL_XicZero;
  Float_t DecayLxy_Lam, err_DecayLxy_Lam;
  Float_t DecayLxy_Xi, err_DecayLxy_Xi;
  Float_t DecayLxy_XicZero, err_DecayLxy_XicZero;
  kfpLambda_Xi.GetDecayLength(DecayL_Lam, err_DecayL_Lam);
  kfpXiMinus_Xic.GetDecayLength(DecayL_Xi, err_DecayL_Xi);
//  kfpXic0_pv.GetDecayLength(DecayL_XicZero, err_DecayL_XicZero);
  kfpLambda_Xi.GetDecayLengthXY(DecayLxy_Lam, err_DecayLxy_Lam);
  kfpXiMinus_Xic.GetDecayLengthXY(DecayLxy_Xi, err_DecayLxy_Xi);
//  kfpXic0_pv.GetDecayLengthXY(DecayLxy_XicZero, err_DecayLxy_XicZero);
  fVar_Xic0[13] = DecayL_Lam;
  fVar_Xic0[14] = DecayL_Xi;
//  fVar_Xic0[34] = DecayL_XicZero;
  fVar_Xic0[15] = DecayLxy_Lam;
  fVar_Xic0[16] = DecayLxy_Xi;
//  fVar_Xic0[37] = DecayLxy_XicZero;

  fVar_Xic0[17] = TMath::ACos(cosPA_v0); // PA_Lam
  fVar_Xic0[18] = TMath::ACos(cosPA_XiToPV); // PA_Xi

  fVar_Xic0[19] = kfpLambda.GetPt();
  fVar_Xic0[20] = kfpLambda.GetRapidity();
  Float_t mass_Lam, err_mass_Lam;
  kfpLambda.GetMass(mass_Lam, err_mass_Lam);
  fVar_Xic0[21] = mass_Lam;

  fVar_Xic0[22] = kfpXiMinus.GetPt();
  fVar_Xic0[23] = kfpXiMinus.GetRapidity();
  Float_t mass_Xi, err_mass_Xi;
  kfpXiMinus.GetMass(mass_Xi, err_mass_Xi);
  fVar_Xic0[24] = mass_Xi;

  fVar_Xic0[25] = trackPiFromXic0->Pt();
  fVar_Xic0[26] = trackPiFromXic0->Eta();

  fVar_Xic0[27] = kfpXic0_pv.GetPt();
  if ( TMath::Abs(kfpXic0_pv.GetE())>TMath::Abs(kfpXic0_pv.GetPz()) ) {
    fVar_Xic0[28] = kfpXic0_pv.GetRapidity();
  }

  Float_t massK0S_Rec, err_massK0S;
  kfpK0Short.GetMass(massK0S_Rec, err_massK0S);
  fVar_Xic0[30] = massK0S_Rec;
  if (fIsMC) {
    fVar_Xic0[31] = lab_Xic0;
  }

  fVar_Xic0[32] = CosThetaStarKF(0, 4132, 211, 3312, kfpXic0, kfpBP_Xic0, kfpXiMinus_Xic0);
  fVar_Xic0[33] = CosThetaStarKF(1, 4132, 211, 3312, kfpXic0, kfpBP_Xic0, kfpXiMinus_Xic0);
//  fVar_Xic0[32] = CosThetaStarKF(0, 4132, 211, 3312, kfpXic0, kfpBP, kfpXiMinus_m);
//  fVar_Xic0[33] = CosThetaStarKF(1, 4132, 211, 3312, kfpXic0, kfpBP, kfpXiMinus_m);

//  fVar_Xic0[3]  = TMath::Prob(kfpLambda.GetChi2(), kfpLambda.GetNDF());
//  fVar_Xic0[7]  = kfpLambda.GetPz();
//  fVar_Xic0[8]  = kfpLambda.GetX();
//  fVar_Xic0[9]  = kfpLambda.GetY();
//  fVar_Xic0[10] = kfpLambda.GetZ();
//  fVar_Xic0[8] = kfpLambda_pv.GetChi2()/kfpLambda_pv.GetNDF();
//  fVar_Xic0[9] = TMath::Prob(kfpLambda_pv.GetChi2(), kfpLambda_pv.GetNDF());
//  fVar_Xic0[15] = alpha_FirstDaugPos;
//  fVar_Xic0[16] = v0->PtArmV0();

//  fVar_Xic0[17] = trackPiFromXi->Charge();
//  fVar_Xic0[20] = TMath::Prob(kfpXiMinus.GetChi2(), kfpXiMinus.GetNDF());
//  fVar_Xic0[21] = kfpXiMinus_pv.GetChi2()/kfpXiMinus_pv.GetNDF();
//  fVar_Xic0[22] = TMath::Prob(kfpXiMinus_pv.GetChi2(), kfpXiMinus_pv.GetNDF());
//  Float_t mass_Xi_PV, err_mass_Xi_PV;
//  kfpXiMinus_pv.GetMass(mass_Xi_PV, err_mass_Xi_PV);
//  fVar_Xic0[23] = mass_Xi_PV;
//  if ( TMath::Abs(kfpXiMinus_pv.GetE())>TMath::Abs(kfpXiMinus_pv.GetPz()) ) {
//    fVar_Xic0[24] = kfpXiMinus_pv.GetRapidity();
//  }
//  fVar_Xic0[25] = kfpXiMinus_pv.GetPt();
//  fVar_Xic0[26] = kfpXiMinus_pv.GetPz();
//  fVar_Xic0[27] = kfpXiMinus_pv.GetX();
//  fVar_Xic0[28] = kfpXiMinus_pv.GetY();
//  fVar_Xic0[29] = kfpXiMinus_pv.GetZ();
//  fVar_Xic0[33] = kfpXiMinus.GetPz();
//  fVar_Xic0[34] = kfpXiMinus.GetX();
//  fVar_Xic0[35] = kfpXiMinus.GetY();
//  fVar_Xic0[36] = kfpXiMinus.GetZ();

//  fVar_Xic0[39] = trackPiFromXic0->Charge();
//  fVar_Xic0[41] = trackPiFromXic0->Pz();

//  fVar_Xic0[20] = kfpXic0.GetChi2()/kfpXic0.GetNDF();
//  fVar_Xic0[45] = TMath::Prob(kfpXic0.GetChi2(), kfpXic0.GetNDF());
//  fVar_Xic0[21] = kfpXic0_pv.GetChi2()/kfpXic0_pv.GetNDF();
//  fVar_Xic0[47] = TMath::Prob(kfpXic0_pv.GetChi2(), kfpXic0_pv.GetNDF());
//  fVar_Xic0[51] = kfpXic0_pv.GetPz();
//  fVar_Xic0[52] = kfpXic0_pv.GetX();
//  fVar_Xic0[53] = kfpXic0_pv.GetY();
//  fVar_Xic0[54] = kfpXic0_pv.GetZ();
//  Float_t mass_Xic, err_mass_Xic;
//  kfpXic0.GetMass(mass_Xic, err_mass_Xic);
//  fVar_Xic0[55] = mass_Xic;
//  fVar_Xic0[56] = kfpXic0.GetRapidity();
//  fVar_Xic0[57] = kfpXic0.GetPt();
//  fVar_Xic0[58] = kfpXic0.GetPz();
//  fVar_Xic0[59] = kfpXic0.GetX();
//  fVar_Xic0[60] = kfpXic0.GetY();
//  fVar_Xic0[61] = kfpXic0.GetZ();

//  fVar_Xic0[27] = cosPA_XiToXic;


//  Float_t CT_Lam, err_CT_Lam;
//  Float_t CT_Xi, err_CT_Xi;
//  Float_t CT_XicZero, err_CT_XicZero;
//  kfpLambda_Xi.GetLifeTime(CT_Lam, err_CT_Lam);
//  kfpXiMinus_pv.GetLifeTime(CT_Xi, err_CT_Xi);
//  kfpXic0_pv.GetLifeTime(CT_XicZero, err_CT_XicZero);
//  fVar_Xic0[38] = CT_Lam;
//  fVar_Xic0[39] = CT_Xi;
//  fVar_Xic0[40] = CT_XicZero;

//  if (fIsMC) {
//    fVar_Xic0[78] = lab_Xic0;
//  }


  fTree_Xic0->Fill();

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEXicZero2XiPifromKFP::FillTreeRecXic0FromCasc(KFParticle kfpXic0, AliAODTrack *trackPiFromXic0, KFParticle kfpBP, KFParticle kfpXiMinus, KFParticle kfpXiMinus_m, AliAODTrack *trackPiFromXi, AliAODcascade *casc, KFParticle kfpK0Short, KFParticle kfpGamma, KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkProton, AliAODTrack *trkPion, KFParticle PV, TClonesArray *mcArray, Int_t lab_Xic0)
{

  for (Int_t i=0; i<37; i++) {
    fVar_Xic0[i] = -9999.;
  }

//  Float_t nSigmaTOF_PiFromXic0 = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trackPiFromXic0,AliPID::kPion);
//  Float_t nSigmaTOF_PiFromXi   = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTOF(trackPiFromXi,AliPID::kPion);

//  Float_t nSigmaTPC_PiFromXic0 = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trackPiFromXic0,AliPID::kPion);
//  Float_t nSigmaTPC_PiFromXi   = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trackPiFromXi,AliPID::kPion);
//  Float_t nSigmaTPC_PrFromLam  = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkProton,AliPID::kProton);
//  Float_t nSigmaTPC_PiFromLam  = fAnaCuts->GetPidHF()->GetPidResponse()->NumberOfSigmasTPC(trkPion,AliPID::kPion);

  Float_t nSigmaTOF_PiFromXic0 = fPID->NumberOfSigmasTOF(trackPiFromXic0,AliPID::kPion);
  Float_t nSigmaTOF_PiFromXi   = fPID->NumberOfSigmasTOF(trackPiFromXi,AliPID::kPion);
  Float_t nSigmaTPC_PiFromXic0 = fPID->NumberOfSigmasTPC(trackPiFromXic0,AliPID::kPion);
  Float_t nSigmaTPC_PiFromXi   = fPID->NumberOfSigmasTPC(trackPiFromXi,AliPID::kPion);
  Float_t nSigmaTPC_PrFromLam  = fPID->NumberOfSigmasTPC(trkProton,AliPID::kProton);
  Float_t nSigmaTPC_PiFromLam  = fPID->NumberOfSigmasTPC(trkPion,AliPID::kPion);

  if ( fabs(nSigmaTPC_PiFromXic0)>=4. || fabs(nSigmaTPC_PiFromXi)>=4. || fabs(nSigmaTPC_PrFromLam)>=4. || fabs(nSigmaTPC_PiFromLam)>=4. ) return;

//  AliAODTrack *trk0 = (AliAODTrack*) (v0->GetDaughter(0));

//  Double_t alpha_FirstDaugPos = (v0->MomPosAlongV0() - v0->MomNegAlongV0())/(v0->MomPosAlongV0() + v0->MomNegAlongV0());

//  if (trk0->Charge()<0) {
//    alpha_FirstDaugPos = (v0->MomNegAlongV0() - v0->MomPosAlongV0())/(v0->MomPosAlongV0() + v0->MomNegAlongV0());
//  }

  KFParticle kfpXic0_PV = kfpXic0;
  kfpXic0_PV.SetProductionVertex(PV);

  KFParticle kfpXiMinus_Xic0 = kfpXiMinus_m;
  kfpXiMinus_Xic0.SetProductionVertex(kfpXic0);
  KFParticle kfpXiMinus_PV = kfpXiMinus_m;
  kfpXiMinus_PV.SetProductionVertex(PV);

  KFParticle kfpLambda_Xi = kfpLambda_m;
  kfpLambda_Xi.SetProductionVertex(kfpXiMinus);
  KFParticle kfpLambda_PV = kfpLambda_m;
  kfpLambda_PV.SetProductionVertex(PV);
//  KFParticle kfpXiMinus_pv = kfpXiMinus;
//  kfpXiMinus_pv.SetProductionVertex(PV);

  KFParticle kfpBP_Xic0 = kfpBP;
  kfpBP_Xic0.SetProductionVertex(kfpXic0);

  //************************** calculate l/Δl for Lambda *************************************
  Double_t dx_Lambda = PV.GetX()-kfpLambda.GetX();
  Double_t dy_Lambda = PV.GetY()-kfpLambda.GetY();
  Double_t dz_Lambda = PV.GetZ()-kfpLambda.GetZ();
  Double_t l_Lambda = TMath::Sqrt(dx_Lambda*dx_Lambda + dy_Lambda*dy_Lambda + dz_Lambda*dz_Lambda);
  Double_t dl_Lambda = (PV.GetCovariance(0)+kfpLambda.GetCovariance(0))*dx_Lambda*dx_Lambda + (PV.GetCovariance(2)+kfpLambda.GetCovariance(2))*dy_Lambda*dy_Lambda + (PV.GetCovariance(5)+kfpLambda.GetCovariance(5))*dz_Lambda*dz_Lambda + 2*( (PV.GetCovariance(1)+kfpLambda.GetCovariance(1))*dx_Lambda*dy_Lambda + (PV.GetCovariance(3)+kfpLambda.GetCovariance(3))*dx_Lambda*dz_Lambda + (PV.GetCovariance(4)+kfpLambda.GetCovariance(4))*dy_Lambda*dz_Lambda );
  if ( fabs(l_Lambda)<1.e-8f ) l_Lambda = 1.e-8f;
  dl_Lambda = dl_Lambda<0. ? 1.e8f : sqrt(dl_Lambda)/l_Lambda;
  if ( dl_Lambda<=0 ) return;
  //***************************************************************************************
  //************************** calculate l/Δl for Xi- *************************************
  Double_t dx_Xi = PV.GetX()-kfpXiMinus.GetX();
  Double_t dy_Xi = PV.GetY()-kfpXiMinus.GetY();
  Double_t dz_Xi = PV.GetZ()-kfpXiMinus.GetZ();
  Double_t l_Xi = TMath::Sqrt(dx_Xi*dx_Xi + dy_Xi*dy_Xi + dz_Xi*dz_Xi);
  Double_t dl_Xi = (PV.GetCovariance(0)+kfpXiMinus.GetCovariance(0))*dx_Xi*dx_Xi + (PV.GetCovariance(2)+kfpXiMinus.GetCovariance(2))*dy_Xi*dy_Xi + (PV.GetCovariance(5)+kfpXiMinus.GetCovariance(5))*dz_Xi*dz_Xi + 2*( (PV.GetCovariance(1)+kfpXiMinus.GetCovariance(1))*dx_Xi*dy_Xi + (PV.GetCovariance(3)+kfpXiMinus.GetCovariance(3))*dx_Xi*dz_Xi + (PV.GetCovariance(4)+kfpXiMinus.GetCovariance(4))*dy_Xi*dz_Xi );
  if ( fabs(l_Xi)<1.e-8f ) l_Xi = 1.e-8f;
  dl_Xi = dl_Xi<0. ? 1.e8f : sqrt(dl_Xi)/l_Xi;
  if ( dl_Xi<=0 ) return;
  //***************************************************************************************

  if ( kfpLambda_PV.GetChi2()/kfpLambda_PV.GetNDF() <= fAnaCuts->GetKFPLam_Chi2topoMin() ) return;
  if ( kfpXiMinus_PV.GetChi2()/kfpXiMinus_PV.GetNDF() >= fAnaCuts->GetKFPXi_Chi2topoMax() ) return;
  if ( l_Xi/dl_Xi <= fAnaCuts->GetKFPXi_lDeltalMin() ) return;

  const Float_t PDGmassXic0 = TDatabasePDG::Instance()->GetParticle(4132)->Mass();
  Float_t mass_Xic0_PV, err_mass_Xic0_PV;
  kfpXic0_PV.GetMass(mass_Xic0_PV, err_mass_Xic0_PV);
  fVar_Xic0[25] = mass_Xic0_PV; // mass of Xic0

  if ( fabs(mass_Xic0_PV-PDGmassXic0) > fAnaCuts->GetProdMassTolXic0() ) return;


  fVar_Xic0[0]  = nSigmaTPC_PiFromXic0;
  fVar_Xic0[1]  = nSigmaTOF_PiFromXic0;
  fVar_Xic0[2]  = nSigmaTPC_PiFromXi;
  fVar_Xic0[3]  = nSigmaTOF_PiFromXi;
  fVar_Xic0[4]  = nSigmaTPC_PiFromLam; 
  fVar_Xic0[5]  = nSigmaTPC_PrFromLam;

  fVar_Xic0[6]  = casc->DcaXiDaughters(); // DCA_LamDau

  fVar_Xic0[7]  = kfpLambda.GetChi2()/kfpLambda.GetNDF(); // chi2geo_Lam
  fVar_Xic0[8]  = l_Lambda/dl_Lambda;
  fVar_Xic0[9]  = kfpLambda_PV.GetChi2()/kfpLambda_PV.GetNDF();
  fVar_Xic0[10] = kfpXiMinus.GetChi2()/kfpXiMinus.GetNDF();
  fVar_Xic0[11] = l_Xi/dl_Xi;
  fVar_Xic0[12] = kfpXiMinus_PV.GetChi2()/kfpXiMinus_PV.GetNDF();

  Float_t DecayLxy_Lam, err_DecayLxy_Lam;
  kfpLambda_Xi.GetDecayLengthXY(DecayLxy_Lam, err_DecayLxy_Lam);
  fVar_Xic0[13] = DecayLxy_Lam;
  Float_t ct_Lam=0., err_ct_Lam=0.;
  kfpLambda_Xi.GetLifeTime(ct_Lam, err_ct_Lam);
  fVar_Xic0[14] = ct_Lam;

  Float_t DecayLxy_Xi, err_DecayLxy_Xi;
  kfpXiMinus_PV.GetDecayLengthXY(DecayLxy_Xi, err_DecayLxy_Xi);
  fVar_Xic0[15] = DecayLxy_Xi;
  Float_t ct_Xi=0., err_ct_Xi=0.;
  kfpXiMinus_PV.GetLifeTime(ct_Xi, err_ct_Xi);
  fVar_Xic0[16] = ct_Xi;

  // calculate CosPointingAngle
  Double_t cosPA_v0toXi = CosPointingAngleKF(kfpLambda_m, kfpXiMinus);
  Double_t cosPA_v0toPV = CosPointingAngleKF(kfpLambda_m, PV);
  Double_t cosPA_XiToPV = CosPointingAngleKF(kfpXiMinus_m, PV);
  fVar_Xic0[17] = TMath::ACos(cosPA_v0toXi); // PA_LamToXi
  fVar_Xic0[18] = TMath::ACos(cosPA_v0toPV); // PA_LamToPV
  fVar_Xic0[19] = TMath::ACos(cosPA_XiToPV); // PA_XiToPV

  Float_t mass_Lam, err_mass_Lam;
  kfpLambda.GetMass(mass_Lam, err_mass_Lam);
  fVar_Xic0[20] = mass_Lam;

  Float_t mass_Xi, err_mass_Xi;
  kfpXiMinus.GetMass(mass_Xi, err_mass_Xi);
  fVar_Xic0[21] = mass_Xi;

  fVar_Xic0[22] = trackPiFromXic0->Pt();

  fVar_Xic0[23] = kfpXic0_PV.GetPt();
  if ( TMath::Abs(kfpXic0_PV.GetE())>TMath::Abs(kfpXic0_PV.GetPz()) ) {
    fVar_Xic0[24] = kfpXic0_PV.GetRapidity();
  }

  fVar_Xic0[26] = CosThetaStarKF(0, 4132, 211, 3312, kfpXic0, kfpBP_Xic0, kfpXiMinus_Xic0);
  fVar_Xic0[27] = CosThetaStarKF(1, 4132, 211, 3312, kfpXic0, kfpBP_Xic0, kfpXiMinus_Xic0);

  // --- chi2_prim of Pion to PV ---
  KFParticle kfpBP_PV = kfpBP;
  kfpBP_PV.SetProductionVertex(PV);
  fVar_Xic0[28] = kfpBP_PV.GetChi2()/kfpBP_PV.GetNDF();
  // -------------------------------
  // --- DCA of Pion to PV ---
  fVar_Xic0[29] = kfpBP.GetDistanceFromVertexXY(PV); // DCA of pion in x-y
  // -------------------------

  Float_t massK0S_Rec, err_massK0S;
  kfpK0Short.GetMass(massK0S_Rec, err_massK0S);
  fVar_Xic0[31] = massK0S_Rec;
  Float_t massGamma_Rec, err_massGamma;
  kfpGamma.GetMass(massGamma_Rec, err_massGamma);
  fVar_Xic0[32] = massGamma_Rec;
  fVar_Xic0[33] = kfpLambda_Xi.GetChi2()/kfpLambda_Xi.GetNDF();
  fVar_Xic0[34] = kfpXiMinus_Xic0.GetChi2()/kfpXiMinus_Xic0.GetNDF();

  Float_t DecayLxy_Xic0, err_DecayLxy_Xic0;
  kfpXic0_PV.GetDecayLengthXY(DecayLxy_Xic0, err_DecayLxy_Xic0);
  fVar_Xic0[35] = DecayLxy_Xic0;
  Float_t ct_Xic0=0., err_ct_Xic0=0.;
  kfpXic0_PV.GetLifeTime(ct_Xic0, err_ct_Xic0);
  fVar_Xic0[36] = ct_Xic0;

  if (fIsMC) {
    fVar_Xic0[30] = lab_Xic0;
    if (lab_Xic0>0 && fabs(fVar_Xic0[24])<0.8 && (massK0S_Rec<=0.487614 || massK0S_Rec>=0.507614)) {
      Int_t labelPion1 = fabs(trackPiFromXic0->GetLabel());
      AliAODMCParticle* mcPion1 = static_cast<AliAODMCParticle*>(mcArray->At(labelPion1));
      AliAODMCParticle* mcXic0  = static_cast<AliAODMCParticle*>(mcArray->At(mcPion1->GetMother()));
      f2DHistMCRec_Xic0Pt_weight->Fill(fVar_Xic0[27], fVar_Xic0[25], fWeight->Eval(mcXic0->Pt()));
    }
  }

  // pt(Xic0)>=1
  if (fVar_Xic0[23]>0.9999) fTree_Xic0->Fill();



//  fVar_Xic0[32] = CosThetaStarKF(0, 4132, 211, 3312, kfpXic0, kfpBP, kfpXiMinus_m);
//  fVar_Xic0[33] = CosThetaStarKF(1, 4132, 211, 3312, kfpXic0, kfpBP, kfpXiMinus_m);

//  fVar_Xic0[3]  = TMath::Prob(kfpLambda.GetChi2(), kfpLambda.GetNDF());
//  fVar_Xic0[7]  = kfpLambda.GetPz();
//  fVar_Xic0[8]  = kfpLambda.GetX();
//  fVar_Xic0[9]  = kfpLambda.GetY();
//  fVar_Xic0[10] = kfpLambda.GetZ();
//  fVar_Xic0[8] = kfpLambda_pv.GetChi2()/kfpLambda_pv.GetNDF();
//  fVar_Xic0[9] = TMath::Prob(kfpLambda_pv.GetChi2(), kfpLambda_pv.GetNDF());
//  fVar_Xic0[15] = alpha_FirstDaugPos;
//  fVar_Xic0[16] = v0->PtArmV0();

//  fVar_Xic0[17] = trackPiFromXi->Charge();
//  fVar_Xic0[20] = TMath::Prob(kfpXiMinus.GetChi2(), kfpXiMinus.GetNDF());
//  fVar_Xic0[21] = kfpXiMinus_pv.GetChi2()/kfpXiMinus_pv.GetNDF();
//  fVar_Xic0[22] = TMath::Prob(kfpXiMinus_pv.GetChi2(), kfpXiMinus_pv.GetNDF());
//  Float_t mass_Xi_PV, err_mass_Xi_PV;
//  kfpXiMinus_pv.GetMass(mass_Xi_PV, err_mass_Xi_PV);
//  fVar_Xic0[23] = mass_Xi_PV;
//  if ( TMath::Abs(kfpXiMinus_pv.GetE())>TMath::Abs(kfpXiMinus_pv.GetPz()) ) {
//    fVar_Xic0[24] = kfpXiMinus_pv.GetRapidity();
//  }
//  fVar_Xic0[25] = kfpXiMinus_pv.GetPt();
//  fVar_Xic0[26] = kfpXiMinus_pv.GetPz();
//  fVar_Xic0[27] = kfpXiMinus_pv.GetX();
//  fVar_Xic0[28] = kfpXiMinus_pv.GetY();
//  fVar_Xic0[29] = kfpXiMinus_pv.GetZ();
//  fVar_Xic0[33] = kfpXiMinus.GetPz();
//  fVar_Xic0[34] = kfpXiMinus.GetX();
//  fVar_Xic0[35] = kfpXiMinus.GetY();
//  fVar_Xic0[36] = kfpXiMinus.GetZ();

//  fVar_Xic0[39] = trackPiFromXic0->Charge();
//  fVar_Xic0[41] = trackPiFromXic0->Pz();

//  fVar_Xic0[20] = kfpXic0.GetChi2()/kfpXic0.GetNDF();
//  fVar_Xic0[45] = TMath::Prob(kfpXic0.GetChi2(), kfpXic0.GetNDF());
//  fVar_Xic0[21] = kfpXic0_PV.GetChi2()/kfpXic0_PV.GetNDF();
//  fVar_Xic0[47] = TMath::Prob(kfpXic0_PV.GetChi2(), kfpXic0_PV.GetNDF());
//  fVar_Xic0[51] = kfpXic0_PV.GetPz();
//  fVar_Xic0[52] = kfpXic0_PV.GetX();
//  fVar_Xic0[53] = kfpXic0_PV.GetY();
//  fVar_Xic0[54] = kfpXic0_PV.GetZ();
//  Float_t mass_Xic, err_mass_Xic;
//  kfpXic0.GetMass(mass_Xic, err_mass_Xic);
//  fVar_Xic0[55] = mass_Xic;
//  fVar_Xic0[56] = kfpXic0.GetRapidity();
//  fVar_Xic0[57] = kfpXic0.GetPt();
//  fVar_Xic0[58] = kfpXic0.GetPz();
//  fVar_Xic0[59] = kfpXic0.GetX();
//  fVar_Xic0[60] = kfpXic0.GetY();
//  fVar_Xic0[61] = kfpXic0.GetZ();

//  fVar_Xic0[27] = cosPA_XiToXic;


//  Float_t CT_Lam, err_CT_Lam;
//  Float_t CT_Xi, err_CT_Xi;
//  Float_t CT_XicZero, err_CT_XicZero;
//  kfpLambda_Xi.GetLifeTime(CT_Lam, err_CT_Lam);
//  kfpXiMinus_pv.GetLifeTime(CT_Xi, err_CT_Xi);
//  kfpXic0_PV.GetLifeTime(CT_XicZero, err_CT_XicZero);
//  fVar_Xic0[38] = CT_Lam;
//  fVar_Xic0[39] = CT_Xi;
//  fVar_Xic0[40] = CT_XicZero;

//  if (fIsMC) {
//    fVar_Xic0[78] = lab_Xic0;
//  }

  return;
}
