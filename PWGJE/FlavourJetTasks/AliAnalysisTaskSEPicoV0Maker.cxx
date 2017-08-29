/**************************************************************************
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

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for RecoDecay object (K0 short, Lambda,
// D mesons ...) filtering
//
// Author: X-M. Zhang, xmzhang@lbl.gov
/////////////////////////////////////////////////////////////

#include <TH2D.h>
#include <TMath.h>
#include <THnSparse.h>
#include <TClonesArray.h>
#include <TParticle.h>

#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"

#include "AliHeader.h"
#include "AliInputEventHandler.h"
#include "AliGenDPMjetEventHeader.h"

#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliV0vertexer.h"

#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"

#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"

#include "AliCentrality.h"
#include "AliMultSelection.h"

#include "AliPicoBase.h"
#include "AliPicoV0RD.h"
#include "AliPicoV0MC.h"
#include "AliAnalysisTaskSEPicoV0Maker.h"

ClassImp(AliAnalysisTaskSEPicoV0Maker)
//=============================================================================

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Maker::AliAnalysisTaskSEPicoV0Maker() :
AliAnalysisTaskSE(),
fTriggerMask(0),
fCollisionType(0),
fIsAnaUseMC(kFALSE),
fIsDPMjetMC(kFALSE),
fUseMultOld(kFALSE),
fUseAnaUtils(kFALSE),
fIsSkipFastOnly(kFALSE),
fIsRefitV0sESD(kFALSE),
fRapidityShift(0.),
fMultEstDef(""),
fCutMinMult(0.),
fCutMaxMult(0.),
fCutMinV0Pt(0.),
fCutMaxV0Pt(0.),
fCutMinV0Rap(0.),
fCutMaxV0Rap(0.),
fCutMinDauPt(0.),
fCutMinDauEta(0.),
fCutMaxDauEta(0.),
fCutMaxV0Chi2(0.),
fCutMinV0Radius(0.),
fCutMaxV0Radius(0.),
fCutMaxDausDCA(0.),
fCutMinDauDCAtoPV(0.),
fCutMinDauXrowsTPC(0.),
fCutMinDauXrowsOverFindableClusTPC(0.),
fCutMaxKshortSigmaTPC(0.),
fCutMinKshortCosPA(0.),
fCutMaxKshortCtau(0.),
fCutMaxKshortArmFrac(0.),
fCutMinKshortDeltaM(0.),
fCutMaxLambdaSigmaTPC(0.),
fCutMinLambdaCosPA(0.),
fCutMaxLambdaCtau(0.),
fCutMaxLambdaArmFrac(0.),
fCutMinLambdaDeletaM(0.),
fEventAOD(nullptr),
fEventESD(nullptr),
fRespoPID(nullptr),
fEventAcptMask(0),
fMultEsti(),
fPicoV0sClArr(nullptr),
fOutputListEH(nullptr),
fOutputListMC(nullptr)
{
//
// Default constructor
//

  for (auto &d : fPrimaryVtx) d = -999.;
}

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Maker::AliAnalysisTaskSEPicoV0Maker(const char *name, Bool_t bIsMC) :
AliAnalysisTaskSE(name),
fTriggerMask(AliVEvent::kAny),
fCollisionType(AliPicoBase::kPP),
fIsAnaUseMC(bIsMC),
fIsDPMjetMC(kFALSE),
fUseMultOld(kFALSE),
fUseAnaUtils(kFALSE),
fIsSkipFastOnly(kTRUE),
fIsRefitV0sESD(kFALSE),
fRapidityShift(0.),
fMultEstDef(""),
fCutMinMult(-99999.),
fCutMaxMult(999999.),
fCutMinV0Pt(0.),
fCutMaxV0Pt(100.),
fCutMinV0Rap(-10.),
fCutMaxV0Rap(10.),
fCutMinDauPt(0.),
fCutMinDauEta(-10.),
fCutMaxDauEta(10.),
fCutMaxV0Chi2(33.),
fCutMinV0Radius(0.3),
fCutMaxV0Radius(200.),
fCutMaxDausDCA(1.5),
fCutMinDauDCAtoPV(0.05),
fCutMinDauXrowsTPC(70.),
fCutMinDauXrowsOverFindableClusTPC(0.8),
fCutMaxKshortSigmaTPC(-1.),
fCutMinKshortCosPA(0.95),
fCutMaxKshortCtau(30.),
fCutMaxKshortArmFrac(-1.),
fCutMinKshortDeltaM(0.003),
fCutMaxLambdaSigmaTPC(7.),
fCutMinLambdaCosPA(0.993),
fCutMaxLambdaCtau(40.),
fCutMaxLambdaArmFrac(-1.),
fCutMinLambdaDeletaM(-1.),
fEventAOD(nullptr),
fEventESD(nullptr),
fRespoPID(nullptr),
fEventAcptMask(0),
fMultEsti(),
fPicoV0sClArr(nullptr),
fOutputListEH(nullptr),
fOutputListMC(nullptr)
{
//
// Constructor
//

  for (auto &d : fPrimaryVtx) d = -999.;

  DefineOutput(1, TList::Class());
  if (fIsAnaUseMC) DefineOutput(2, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0Maker::~AliAnalysisTaskSEPicoV0Maker()
{
//
// Default destructor
//

  if (fEventAOD) { delete fEventAOD; fEventAOD = nullptr; }
  if (fEventESD) { delete fEventESD; fEventESD = nullptr; }
  if (fRespoPID) { delete fRespoPID; fRespoPID = nullptr; }

  if (fPicoV0sClArr) { delete fPicoV0sClArr; fPicoV0sClArr = nullptr; }
  if (fOutputListEH) { delete fOutputListEH; fOutputListEH = nullptr; }
  if (fOutputListMC) { delete fOutputListMC; fOutputListMC = nullptr; }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::Init()
{
//
//  AliAnalysisTaskSEPicoV0Maker::Init
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::UserCreateOutputObjects()
{
//
//  AliAnalysisTaskSEPicoV0Maker::UserCreateOutputObjects
//

  InitAnalysis();
//=============================================================================

  if (fPicoV0sClArr) {
    delete fPicoV0sClArr;
    fPicoV0sClArr = nullptr;
  }

  if (fIsAnaUseMC) {
    fPicoV0sClArr = new TClonesArray("AliPicoV0MC");
    fPicoV0sClArr->SetName("PicoV0s");
  } else {
    fPicoV0sClArr = new TClonesArray("AliPicoV0RD");
    fPicoV0sClArr->SetName("PicoV0s");
  }
//=============================================================================

  if (fOutputListEH) {
    delete fOutputListEH;
    fOutputListEH = nullptr;
  }

  fOutputListEH = new TList();
  fOutputListEH->SetOwner();

  CreateHistogramsEH();
  PostData(1, fOutputListEH);
//=============================================================================

  if (fIsAnaUseMC) {
    if (fOutputListMC) {
      delete fOutputListMC;
      fOutputListMC = nullptr;
    }

    fOutputListMC = new TList();
    fOutputListMC->SetOwner();

    CreateHistogramsMC();
    PostData(2, fOutputListMC);
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::UserExec(Option_t */*opt*/)
{
//
//  AliAnalysisTaskSEPicoV0Maker::UserExec
//

  fPicoV0sClArr->Delete();
  if (!(InputEvent()->FindListObject("PicoV0s"))) InputEvent()->AddObject(fPicoV0sClArr);
//=============================================================================

  if (IsEventNotAcpt()) return;
  FillHistogramsEH();
//=============================================================================

  if (IsEventNotINEL()) return;
  if (fIsAnaUseMC) FillHistogramsMC();
//=============================================================================

  if (IsEventNotMBsa()) return;
//=============================================================================

  FillPicoV0s();
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::Terminate(Option_t */*opt*/)
{
//
//  AliAnalysisTaskSEPicoV0Maker::Terminate
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::NotifyRun()
{
//
//  AliAnalysisTaskSEPicoV0Maker::NotifyRun
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::FillPicoV0s()
{
//
//  AliAnalysisTaskSEPicoV0Maker::FillPicoRecoV0s
//

  const auto nV0s(fEventAOD ? fEventAOD->GetNumberOfV0s() :
                              fEventESD->GetNumberOfV0s());

  if (nV0s<=0) return;
//=============================================================================

  auto nAt(fPicoV0sClArr->GetEntriesFast());
  auto hKshortPtInvM(static_cast<TH2D*>(fOutputListEH->FindObject("hKshortPtInvM")));
  auto hLambdaPtInvM(static_cast<TH2D*>(fOutputListEH->FindObject("hLambdaPtInvM")));
  auto hAntiLaPtInvM(static_cast<TH2D*>(fOutputListEH->FindObject("hAntiLaPtInvM")));
//=============================================================================

  for (auto iV0=0; iV0<nV0s; ++iV0) {
    AliPicoV0RD *pV0RD(nullptr);
    AliPicoV0MC *pV0MC(nullptr);

    if (fEventAOD) {
      auto pV0(fEventAOD->GetV0(iV0));
      if (!pV0) continue;

      if (fIsAnaUseMC) {
        pV0MC = SelectV0CandidateMC(pV0);
      } else {
        pV0RD = SelectV0CandidateRD(pV0);
      }
    }

    if (fEventESD) {
      auto pV0(fEventESD->GetV0(iV0));
      if (!pV0) continue;

      if (fIsAnaUseMC) {
        pV0MC = SelectV0CandidateMC(pV0);
      } else {
        pV0RD = SelectV0CandidateRD(pV0);
      }
    }
//=============================================================================

    if (pV0RD) {
      pV0RD->FillKshortPtInvM(hKshortPtInvM);
      pV0RD->FillLambdaPtInvM(hLambdaPtInvM);
      pV0RD->FillAntiLaPtInvM(hAntiLaPtInvM);
      new ((*fPicoV0sClArr)[nAt++]) AliPicoV0RD(*pV0RD);
      delete pV0RD; pV0RD=nullptr;
    }

    if (pV0MC) {
      pV0MC->FillKshortPtInvM(hKshortPtInvM);
      pV0MC->FillLambdaPtInvM(hLambdaPtInvM);
      pV0MC->FillAntiLaPtInvM(hAntiLaPtInvM);
      new ((*fPicoV0sClArr)[nAt++]) AliPicoV0MC(*pV0MC);
      delete pV0MC; pV0MC=nullptr;
    }
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
AliPicoV0RD *AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateRD(AliAODv0 const *pV0)
{
//
//  AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateRD
//

  if  (pV0->GetOnFlyStatus()) return nullptr;
  if ((pV0->Chi2V0())>fCutMaxV0Chi2) return nullptr;
//=============================================================================

  const auto dV0Pt(pV0->Pt()); if ((dV0Pt<fCutMinV0Pt) || (dV0Pt>fCutMaxV0Pt)) return nullptr;
  const auto dKaRap(pV0->RapK0Short()); if ((dKaRap<fCutMinV0Rap) || (dKaRap>fCutMaxV0Rap)) return nullptr;
  const auto dLaRap(pV0->RapLambda());  if ((dLaRap<fCutMinV0Rap) || (dLaRap>fCutMaxV0Rap)) return nullptr;
//=============================================================================

  Double_t dV0Vtx[3]; pV0->GetXYZ(dV0Vtx);
  const auto dV0Radius(TMath::Sqrt(dV0Vtx[0]*dV0Vtx[0] + dV0Vtx[1]*dV0Vtx[1]));
  if ((dV0Radius<fCutMinV0Radius) || (dV0Radius>fCutMaxV0Radius)) return nullptr;

  const auto dDausDCA(pV0->DcaV0Daughters()); if (dDausDCA>fCutMaxDausDCA) return nullptr;
  const auto dPosDCAtoPV(pV0->DcaPosToPrimVertex()); if (dPosDCAtoPV<fCutMinDauDCAtoPV) return nullptr;
  const auto dNegDCAtoPV(pV0->DcaNegToPrimVertex()); if (dNegDCAtoPV<fCutMinDauDCAtoPV) return nullptr;
//=============================================================================

  auto pDauPos(static_cast<AliAODTrack*>(pV0->GetDaughter(0))); if (!pDauPos) return nullptr;
  auto pDauNeg(static_cast<AliAODTrack*>(pV0->GetDaughter(1))); if (!pDauNeg) return nullptr;

  if (!(pDauPos->GetStatus() & AliESDtrack::kTPCrefit)) return nullptr;
  if (!(pDauNeg->GetStatus() & AliESDtrack::kTPCrefit)) return nullptr;

  if ((pDauPos->GetProdVertex()->GetType())==((Char_t)AliAODVertex::kKink)) return nullptr;
  if ((pDauNeg->GetProdVertex()->GetType())==((Char_t)AliAODVertex::kKink)) return nullptr;

  const auto dPosXrowsTPC(pDauPos->GetTPCClusterInfo(2,1));
  const auto dNegXrowsTPC(pDauNeg->GetTPCClusterInfo(2,1));
  const auto dDauXrowsTPC((dPosXrowsTPC<dNegXrowsTPC) ? dPosXrowsTPC : dNegXrowsTPC);
  if (dDauXrowsTPC<fCutMinDauXrowsTPC) return nullptr;

  const auto wPosTPCNClsF(pDauPos->GetTPCNclsF()); if (wPosTPCNClsF<=0) return nullptr;
  const auto wNegTPCNClsF(pDauNeg->GetTPCNclsF()); if (wNegTPCNClsF<=0) return nullptr;
  const auto dPosXrowsOverFindableClusTPC( ((Double_t)dPosXrowsTPC) / ((Double_t)wPosTPCNClsF) );
  const auto dNegXrowsOverFindableClusTPC( ((Double_t)dNegXrowsTPC) / ((Double_t)wNegTPCNClsF) );

  const auto dDauXrowsOverFindableClusTPC((dPosXrowsOverFindableClusTPC<dNegXrowsOverFindableClusTPC) ?
                                           dPosXrowsOverFindableClusTPC :
                                           dNegXrowsOverFindableClusTPC);
  if (dDauXrowsOverFindableClusTPC<fCutMinDauXrowsOverFindableClusTPC) return nullptr;
//=============================================================================

  const auto nPosCharge(pDauPos->Charge());
  const auto nNegCharge(pDauNeg->Charge());
  if ((nPosCharge==0) || (nNegCharge==0) || (nPosCharge==nNegCharge)) return nullptr;

  Double_t dPosPxPyPz[3] = { 0., 0., 0. };
  Double_t dNegPxPyPz[3] = { 0., 0., 0. };
  if ((nPosCharge<0) && (nNegCharge>0)) {
    pDauPos = (AliAODTrack*)pV0->GetDaughter(1);
    pDauNeg = (AliAODTrack*)pV0->GetDaughter(0);

    dPosPxPyPz[0] = pV0->MomNegX(); dPosPxPyPz[1] = pV0->MomNegY(); dPosPxPyPz[2] = pV0->MomNegZ();
    dNegPxPyPz[0] = pV0->MomPosX(); dNegPxPyPz[1] = pV0->MomPosY(); dNegPxPyPz[2] = pV0->MomPosZ();
  } else {
    dPosPxPyPz[0] = pV0->MomPosX(); dPosPxPyPz[1] = pV0->MomPosY(); dPosPxPyPz[2] = pV0->MomPosZ();
    dNegPxPyPz[0] = pV0->MomNegX(); dNegPxPyPz[1] = pV0->MomNegY(); dNegPxPyPz[2] = pV0->MomNegZ();
  }

  const TVector3 v3Pos(dPosPxPyPz);
  const TVector3 v3Neg(dNegPxPyPz);
  if ((v3Pos.Pt()<fCutMinDauPt) || (v3Neg.Pt()<fCutMinDauPt)) return nullptr;
  const auto dPosEta(v3Pos.Eta()); if ((dPosEta<fCutMinDauEta) || (dPosEta>fCutMaxDauEta)) return nullptr;
  const auto dNegEta(v3Neg.Eta()); if ((dNegEta<fCutMinDauEta) || (dNegEta>fCutMaxDauEta)) return nullptr;
//=============================================================================

  auto bIsKshort(kTRUE);
  auto bIsLambda(kTRUE);
  auto bIsAntiLa(kTRUE);
  const auto dPosPionSigmaTPC(fRespoPID->NumberOfSigmasTPC(pDauPos,AliPID::kPion));
  const auto dNegPionSigmaTPC(fRespoPID->NumberOfSigmasTPC(pDauNeg,AliPID::kPion));

  const auto dPosProtonSigmaTPC(fRespoPID->NumberOfSigmasTPC(pDauPos,AliPID::kProton));
  const auto dNegProtonSigmaTPC(fRespoPID->NumberOfSigmasTPC(pDauNeg,AliPID::kProton));

  if (fCutMaxKshortSigmaTPC>0.) {
    bIsKshort = ((TMath::Abs(dPosPionSigmaTPC)<fCutMaxKshortSigmaTPC) &&
                 (TMath::Abs(dNegPionSigmaTPC)<fCutMaxKshortSigmaTPC));
  }

  if (fCutMaxLambdaSigmaTPC>0.) {
    bIsLambda = ((TMath::Abs(dPosProtonSigmaTPC)<fCutMaxLambdaSigmaTPC) &&
                 (TMath::Abs(dNegPionSigmaTPC)  <fCutMaxLambdaSigmaTPC));

    bIsAntiLa = ((TMath::Abs(dPosPionSigmaTPC)  <fCutMaxLambdaSigmaTPC) &&
                 (TMath::Abs(dNegProtonSigmaTPC)<fCutMaxLambdaSigmaTPC));
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  const auto dV0CosPA(pV0->CosPointingAngle(fPrimaryVtx));

  if (bIsKshort) if (dV0CosPA<fCutMinKshortCosPA) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if (dV0CosPA<fCutMinLambdaCosPA) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  auto dV0DistToPV(0.);
  for (auto i=0; i<3; ++i) dV0DistToPV +=  ((dV0Vtx[i]-fPrimaryVtx[i]) * (dV0Vtx[i]-fPrimaryVtx[i]));
  const auto dV0DistToPVoverP(TMath::Sqrt(dV0DistToPV) / (pV0->P()+1e-10));

  if (bIsKshort) if ((dV0DistToPVoverP*AliPicoBase::MassKshort())>fCutMaxKshortCtau) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if ((dV0DistToPVoverP*AliPicoBase::MassLambda())>fCutMaxLambdaCtau) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  const auto dV0ArmFrac(pV0->PtArmV0() / (TMath::Abs(pV0->AlphaV0())+1e-12));
  if (bIsKshort && (fCutMaxKshortArmFrac>0.)) if (dV0ArmFrac>fCutMaxKshortArmFrac) {
    bIsKshort = kFALSE;
  }

  if ((bIsLambda || bIsAntiLa) && fCutMaxLambdaArmFrac>0.) if (dV0ArmFrac>fCutMaxLambdaArmFrac) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, AliPicoBase::MassPion());
  TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, AliPicoBase::MassPion());

  TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, AliPicoBase::MassProton());
  TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, AliPicoBase::MassProton());

  const auto dKshortInvM((vPosPion  +vNegPion).M());
  const auto dLambdaInvM((vPosProton+vNegPion).M());
  const auto dAntiLaInvM((vNegProton+vPosPion).M());
  if (bIsKshort) if ((dKshortInvM<(0.430006 - 0.0110029*dV0Pt)) ||
                     (dKshortInvM>(0.563707 + 0.0114979*dV0Pt))) bIsKshort = kFALSE;

  if (bIsLambda || bIsAntiLa) {
    const auto dLower(1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt));
    const auto dUpper(1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt));
    if (bIsLambda) if ((dLambdaInvM<dLower) || (dLambdaInvM>dUpper)) bIsLambda = kFALSE;
    if (bIsAntiLa) if ((dAntiLaInvM<dLower) || (dAntiLaInvM>dUpper)) bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  if (bIsKshort && (fCutMinKshortDeltaM>0.)) {
    if ((TMath::Abs(dLambdaInvM-AliPicoBase::MassLambda())<fCutMinKshortDeltaM) ||
        (TMath::Abs(dAntiLaInvM-AliPicoBase::MassLambda())<fCutMinKshortDeltaM)) bIsKshort = kFALSE;
  }

  if ((bIsLambda || bIsAntiLa) && (fCutMinLambdaDeletaM>0.)) {
    if ((TMath::Abs(dKshortInvM-AliPicoBase::MassKshort())<fCutMinLambdaDeletaM)) {
      bIsLambda = kFALSE;
      bIsAntiLa = kFALSE;
    }
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  UInt_t wMask(0);
  if (bIsKshort) wMask |= AliPicoBase::kKshort;
  if (bIsLambda) wMask |= AliPicoBase::kLambda;
  if (bIsAntiLa) wMask |= AliPicoBase::kAntiLambda;
//=============================================================================

  auto bPosInJC(kFALSE);
  auto bNegInJC(kFALSE);
  return (new AliPicoV0RD(wMask,
                          dV0Radius,
                          dV0CosPA,
                          dV0DistToPVoverP,
                          dDausDCA,
                          dPosDCAtoPV,
                          dNegDCAtoPV,
                          dDauXrowsTPC,
                          dDauXrowsOverFindableClusTPC,
                          v3Pos.Px(), v3Pos.Py(), v3Pos.Pz(),
                          v3Neg.Px(), v3Neg.Py(), v3Neg.Pz(),
                          bPosInJC, bNegInJC,
                          dPosPionSigmaTPC, dPosProtonSigmaTPC,
                          dNegPionSigmaTPC, dNegProtonSigmaTPC));
}

//_____________________________________________________________________________
AliPicoV0RD *AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateRD(AliESDv0 const *pV0)
{
//
//  AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateRD
//

  if (pV0->GetOnFlyStatus()) return nullptr;
  if (pV0->GetChi2V0()>fCutMaxV0Chi2) return nullptr;
//=============================================================================

  const auto dV0Pt(pV0->Pt()); if ((dV0Pt<fCutMinV0Pt) || (dV0Pt>fCutMaxV0Pt)) return nullptr;
  const auto dKaRap(pV0->RapK0Short()); if ((dKaRap<fCutMinV0Rap) || (dKaRap>fCutMaxV0Rap)) return nullptr;
  const auto dLaRap(pV0->RapLambda());  if ((dLaRap<fCutMinV0Rap) || (dLaRap>fCutMaxV0Rap)) return nullptr;
//=============================================================================

  Double_t dV0Vtx[3];  pV0->GetXYZ(dV0Vtx[0], dV0Vtx[1], dV0Vtx[2]);
  const auto dV0Radius(TMath::Sqrt(dV0Vtx[0]*dV0Vtx[0] + dV0Vtx[1]*dV0Vtx[1]));
  if ((dV0Radius<fCutMinV0Radius) || (dV0Radius>fCutMaxV0Radius)) return nullptr;

  const auto dDausDCA(pV0->GetDcaV0Daughters()); if (dDausDCA>fCutMaxDausDCA) return nullptr;
//=============================================================================

  const auto nPosIndex(pV0->GetPindex()); if (nPosIndex<0) return nullptr;
  const auto nNegIndex(pV0->GetNindex()); if (nNegIndex<0) return nullptr;

  auto pDauPos(fEventESD->GetTrack(nPosIndex)); if (!pDauPos) return nullptr;
  auto pDauNeg(fEventESD->GetTrack(nNegIndex)); if (!pDauNeg) return nullptr;

  const auto dMegField(fEventESD->GetMagneticField());
  const auto dPosDCAtoPV(TMath::Abs(pDauPos->GetD(fPrimaryVtx[0],fPrimaryVtx[1],dMegField)));
  if (dPosDCAtoPV<fCutMinDauDCAtoPV) return nullptr;

  const auto dNegDCAtoPV(TMath::Abs(pDauNeg->GetD(fPrimaryVtx[0],fPrimaryVtx[1],dMegField)));
  if (dNegDCAtoPV<fCutMinDauDCAtoPV) return nullptr;
//=============================================================================

  if (!(pDauPos->GetStatus() & AliESDtrack::kTPCrefit)) return nullptr;
  if (!(pDauNeg->GetStatus() & AliESDtrack::kTPCrefit)) return nullptr;
  if ((pDauPos->GetKinkIndex(0)>0) || (pDauNeg->GetKinkIndex(0)>0)) return nullptr;

  const auto dPosXrowsTPC(pDauPos->GetTPCClusterInfo(2,1));
  const auto dNegXrowsTPC(pDauNeg->GetTPCClusterInfo(2,1));
  const auto dDauXrowsTPC((dPosXrowsTPC<dNegXrowsTPC) ? dPosXrowsTPC : dNegXrowsTPC);
  if (dDauXrowsTPC<fCutMinDauXrowsTPC) return nullptr;

  const auto wPosTPCNClsF(pDauPos->GetTPCNclsF()); if (wPosTPCNClsF<=0) return nullptr;
  const auto wNegTPCNClsF(pDauNeg->GetTPCNclsF()); if (wNegTPCNClsF<=0) return nullptr;
  const auto dPosXrowsOverFindableClusTPC( ((Double_t)dPosXrowsTPC) / ((Double_t)wPosTPCNClsF) );
  const auto dNegXrowsOverFindableClusTPC( ((Double_t)dNegXrowsTPC) / ((Double_t)wNegTPCNClsF) );

  const auto dDauXrowsOverFindableClusTPC((dPosXrowsOverFindableClusTPC<dNegXrowsOverFindableClusTPC) ?
                                           dPosXrowsOverFindableClusTPC :
                                           dNegXrowsOverFindableClusTPC);
  if (dDauXrowsOverFindableClusTPC<fCutMinDauXrowsOverFindableClusTPC) return nullptr;
//=============================================================================

  const auto nPosCharge(pDauPos->Charge());
  const auto nNegCharge(pDauNeg->Charge());
  if ((nPosCharge==0) || (nNegCharge==0) || (nPosCharge==nNegCharge)) return nullptr;

  Double_t dPosPxPyPz[3] = { 0., 0., 0. };
  Double_t dNegPxPyPz[3] = { 0., 0., 0. };
  if ((nPosCharge<0) && (nNegCharge>0)) {
    pDauPos = fEventESD->GetTrack(nNegIndex);
    pDauNeg = fEventESD->GetTrack(nPosIndex);

    pV0->GetNPxPyPz(dPosPxPyPz[0], dPosPxPyPz[1], dPosPxPyPz[2]);
    pV0->GetPPxPyPz(dNegPxPyPz[0], dNegPxPyPz[1], dNegPxPyPz[2]);
  } else {
    pV0->GetPPxPyPz(dPosPxPyPz[0], dPosPxPyPz[1], dPosPxPyPz[2]);
    pV0->GetNPxPyPz(dNegPxPyPz[0], dNegPxPyPz[1], dNegPxPyPz[2]);
  }

  const TVector3 v3Pos(dPosPxPyPz);
  const TVector3 v3Neg(dNegPxPyPz);
  if ((v3Pos.Pt()<fCutMinDauPt) || (v3Neg.Pt()<fCutMinDauPt)) return nullptr;
  const auto dPosEta(v3Pos.Eta()); if ((dPosEta<fCutMinDauEta) || (dPosEta>fCutMaxDauEta)) return nullptr;
  const auto dNegEta(v3Neg.Eta()); if ((dNegEta<fCutMinDauEta) || (dNegEta>fCutMaxDauEta)) return nullptr;
//=============================================================================

  auto bIsKshort(kTRUE);
  auto bIsLambda(kTRUE);
  auto bIsAntiLa(kTRUE);
  const auto dPosPionSigmaTPC(fRespoPID->NumberOfSigmasTPC(pDauPos,AliPID::kPion));
  const auto dNegPionSigmaTPC(fRespoPID->NumberOfSigmasTPC(pDauNeg,AliPID::kPion));

  const auto dPosProtonSigmaTPC(fRespoPID->NumberOfSigmasTPC(pDauPos,AliPID::kProton));
  const auto dNegProtonSigmaTPC(fRespoPID->NumberOfSigmasTPC(pDauNeg,AliPID::kProton));

  if (fCutMaxKshortSigmaTPC>0.) {
    bIsKshort = ((TMath::Abs(dPosPionSigmaTPC)<fCutMaxKshortSigmaTPC) &&
                 (TMath::Abs(dNegPionSigmaTPC)<fCutMaxKshortSigmaTPC));
  }

  if (fCutMaxLambdaSigmaTPC>0.) {
    bIsLambda = ((TMath::Abs(dPosProtonSigmaTPC)<fCutMaxLambdaSigmaTPC) &&
                 (TMath::Abs(dNegPionSigmaTPC)  <fCutMaxLambdaSigmaTPC));

    bIsAntiLa = ((TMath::Abs(dPosPionSigmaTPC)  <fCutMaxLambdaSigmaTPC) &&
                 (TMath::Abs(dNegProtonSigmaTPC)<fCutMaxLambdaSigmaTPC));
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  const auto dV0CosPA(pV0->GetV0CosineOfPointingAngle(fPrimaryVtx[0],fPrimaryVtx[1],fPrimaryVtx[2]));

  if (bIsKshort) if (dV0CosPA<fCutMinKshortCosPA) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if (dV0CosPA<fCutMinLambdaCosPA) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  auto dV0DistToPV(0.);
  for (auto i=0; i<3; ++i) dV0DistToPV +=  ((dV0Vtx[i]-fPrimaryVtx[i]) * (dV0Vtx[i]-fPrimaryVtx[i]));
  const auto dV0DistToPVoverP(TMath::Sqrt(dV0DistToPV) / (pV0->P()+1e-10));

  if (bIsKshort) if ((dV0DistToPVoverP*AliPicoBase::MassKshort())>fCutMaxKshortCtau) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if ((dV0DistToPVoverP*AliPicoBase::MassLambda())>fCutMaxLambdaCtau) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  const auto dV0ArmFrac(pV0->PtArmV0() / (TMath::Abs(pV0->AlphaV0())+1e-12));
  if (bIsKshort && (fCutMaxKshortArmFrac>0.)) if (dV0ArmFrac>fCutMaxKshortArmFrac) {
    bIsKshort = kFALSE;
  }

  if ((bIsLambda && bIsAntiLa) && (fCutMaxLambdaArmFrac>0.)) if (dV0ArmFrac>fCutMaxLambdaArmFrac) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, AliPicoBase::MassPion());
  TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, AliPicoBase::MassPion());

  TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, AliPicoBase::MassProton());
  TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, AliPicoBase::MassProton());

  const auto dKshortInvM((vPosPion  +vNegPion).M());
  const auto dLambdaInvM((vPosProton+vNegPion).M());
  const auto dAntiLaInvM((vNegProton+vPosPion).M());
  if (bIsKshort) if ((dKshortInvM<(0.430006 - 0.0110029*dV0Pt)) ||
                     (dKshortInvM>(0.563707 + 0.0114979*dV0Pt))) bIsKshort = kFALSE;

  if (bIsLambda || bIsAntiLa) {
    const auto dLower(1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt));
    const auto dUpper(1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt));
    if (bIsLambda) if ((dLambdaInvM<dLower) || (dLambdaInvM>dUpper)) bIsLambda = kFALSE;
    if (bIsAntiLa) if ((dAntiLaInvM<dLower) || (dAntiLaInvM>dUpper)) bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  if (bIsKshort && (fCutMinKshortDeltaM>0.)) {
    if ((TMath::Abs(dLambdaInvM-AliPicoBase::MassLambda())<fCutMinKshortDeltaM) ||
        (TMath::Abs(dAntiLaInvM-AliPicoBase::MassLambda())<fCutMinKshortDeltaM)) {
      bIsKshort = kFALSE;
    }
  }

  if ((bIsLambda || bIsAntiLa) && (fCutMinLambdaDeletaM>0.)) {
    if ((TMath::Abs(dKshortInvM-AliPicoBase::MassKshort())<fCutMinLambdaDeletaM)) {
      bIsLambda = kFALSE;
      bIsAntiLa = kFALSE;
    }
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  UInt_t wMask(0);
  if (bIsKshort) wMask |= AliPicoBase::kKshort;
  if (bIsLambda) wMask |= AliPicoBase::kLambda;
  if (bIsAntiLa) wMask |= AliPicoBase::kAntiLambda;
//=============================================================================

  auto bPosInJC(kFALSE);
  auto bNegInJC(kFALSE);
/*if (fJetContisClArr) {
    const auto idPos(pDauPos->GetID());
    const auto idNeg(pDauNeg->GetID());
    for (auto i=0; i<fJetContisClArr->GetEntriesFast(); ++i) {
      auto pTrk(fJetContisClArr->At(i)); if (!pTrk) continue;

      const auto id(pTrkAOD->GetID());
      if (idPos==id) bPosInJC = kTRUE;
      if (idNeg==id) bNegInJC = kTRUE;
      if (bPosInJC && bNegInJC) break;
    }
  }*/
//=============================================================================

  return (new AliPicoV0RD(wMask,
                          dV0Radius,
                          dV0CosPA,
                          dV0DistToPVoverP,
                          dDausDCA,
                          dPosDCAtoPV,
                          dNegDCAtoPV,
                          dDauXrowsTPC,
                          dDauXrowsOverFindableClusTPC,
                          v3Pos.Px(), v3Pos.Py(), v3Pos.Pz(),
                          v3Neg.Px(), v3Neg.Py(), v3Neg.Pz(),
                          bPosInJC, bNegInJC,
                          dPosPionSigmaTPC, dPosProtonSigmaTPC,
                          dNegPionSigmaTPC, dNegProtonSigmaTPC));
}

//_____________________________________________________________________________
AliPicoV0MC *AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateMC(AliAODv0 const *pV0RD)
{
//
//  AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateMC
//

  if (pV0RD->GetOnFlyStatus()) return nullptr;
  if ((pV0RD->Chi2V0())>fCutMaxV0Chi2) return nullptr;

  const auto dV0Pt(pV0RD->Pt()); if ((dV0Pt<fCutMinV0Pt) || (dV0Pt>fCutMaxV0Pt)) return nullptr;
//=============================================================================

  Double_t dV0Vtx[3]; pV0RD->GetXYZ(dV0Vtx);
  const auto dV0Radius(TMath::Sqrt(dV0Vtx[0]*dV0Vtx[0] + dV0Vtx[1]*dV0Vtx[1]));
  if ((dV0Radius<fCutMinV0Radius) || (dV0Radius>fCutMaxV0Radius)) return nullptr;

  const auto dDausDCA(pV0RD->DcaV0Daughters()); if (dDausDCA>fCutMaxDausDCA) return nullptr;
  const auto dPosDCAtoPV(pV0RD->DcaPosToPrimVertex()); if (dPosDCAtoPV<fCutMinDauDCAtoPV) return nullptr;
  const auto dNegDCAtoPV(pV0RD->DcaNegToPrimVertex()); if (dNegDCAtoPV<fCutMinDauDCAtoPV) return nullptr;
//=============================================================================

  auto pDauPosRD(static_cast<AliAODTrack*>(pV0RD->GetDaughter(0))); if (!pDauPosRD) return nullptr;
  auto pDauNegRD(static_cast<AliAODTrack*>(pV0RD->GetDaughter(1))); if (!pDauNegRD) return nullptr;

  if (!(pDauPosRD->GetStatus() & AliESDtrack::kTPCrefit)) return nullptr;
  if (!(pDauNegRD->GetStatus() & AliESDtrack::kTPCrefit)) return nullptr;

  if ((pDauPosRD->GetProdVertex()->GetType())==((Char_t)AliAODVertex::kKink)) return nullptr;
  if ((pDauNegRD->GetProdVertex()->GetType())==((Char_t)AliAODVertex::kKink)) return nullptr;

  const auto dPosXrowsTPC(pDauPosRD->GetTPCClusterInfo(2,1));
  const auto dNegXrowsTPC(pDauNegRD->GetTPCClusterInfo(2,1));
  const auto dDauXrowsTPC(dPosXrowsTPC<dNegXrowsTPC ? dPosXrowsTPC : dNegXrowsTPC);
  if (dDauXrowsTPC<fCutMinDauXrowsTPC) return nullptr;

  const auto wPosTPCNClsF(pDauPosRD->GetTPCNclsF()); if (wPosTPCNClsF<=0) return nullptr;
  const auto wNegTPCNClsF(pDauNegRD->GetTPCNclsF()); if (wNegTPCNClsF<=0) return nullptr;
  const auto dPosXrowsOverFindableClusTPC( ((Double_t)dPosXrowsTPC) / ((Double_t)wPosTPCNClsF) );
  const auto dNegXrowsOverFindableClusTPC( ((Double_t)dNegXrowsTPC) / ((Double_t)wNegTPCNClsF) );

  const auto dDauXrowsOverFindableClusTPC(dPosXrowsOverFindableClusTPC<dNegXrowsOverFindableClusTPC ?
                                          dPosXrowsOverFindableClusTPC :
                                          dNegXrowsOverFindableClusTPC);
  if (dDauXrowsOverFindableClusTPC<fCutMinDauXrowsOverFindableClusTPC) return nullptr;
//=============================================================================

  const auto nPosCharge(pDauPosRD->Charge());
  const auto nNegCharge(pDauNegRD->Charge());
  if ((nPosCharge==0) || (nNegCharge==0) || (nPosCharge==nNegCharge)) return nullptr;

  Double_t dPosPxPyPz[3] = { 0., 0., 0. };
  Double_t dNegPxPyPz[3] = { 0., 0., 0. };
  if ((nPosCharge<0) && (nNegCharge>0)) {
    pDauPosRD = (AliAODTrack*)pV0RD->GetDaughter(1);
    pDauNegRD = (AliAODTrack*)pV0RD->GetDaughter(0);

    dPosPxPyPz[0] = pV0RD->MomNegX(); dPosPxPyPz[1] = pV0RD->MomNegY(); dPosPxPyPz[2] = pV0RD->MomNegZ();
    dNegPxPyPz[0] = pV0RD->MomPosX(); dNegPxPyPz[1] = pV0RD->MomPosY(); dNegPxPyPz[2] = pV0RD->MomPosZ();
  } else {
    dPosPxPyPz[0] = pV0RD->MomPosX(); dPosPxPyPz[1] = pV0RD->MomPosY(); dPosPxPyPz[2] = pV0RD->MomPosZ();
    dNegPxPyPz[0] = pV0RD->MomNegX(); dNegPxPyPz[1] = pV0RD->MomNegY(); dNegPxPyPz[2] = pV0RD->MomNegZ();
  }

  const TVector3 v3Pos(dPosPxPyPz);
  const TVector3 v3Neg(dNegPxPyPz);
  if ((v3Pos.Pt()<fCutMinDauPt) || (v3Neg.Pt()<fCutMinDauPt)) return nullptr;
  const auto dPosEta(v3Pos.Eta()); if ((dPosEta<fCutMinDauEta) || (dPosEta>fCutMaxDauEta)) return nullptr;
  const auto dNegEta(v3Neg.Eta()); if ((dNegEta<fCutMinDauEta) || (dNegEta>fCutMaxDauEta)) return nullptr;
//=============================================================================

  const auto inp(TMath::Abs(pDauPosRD->GetLabel())); if (inp<0) return nullptr;
  auto pDauPosMC(static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(inp))); if (!pDauPosMC) return nullptr;
  const auto imp(pDauPosMC->GetMother()); if (imp<0) return nullptr;

  const auto inn(TMath::Abs(pDauNegRD->GetLabel())); if (inn<0) return nullptr;
  const auto pDauNegMC(static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(inn))); if (!pDauNegMC) return nullptr;
  const auto imn(pDauNegMC->GetMother()); if (imn<0) return nullptr;

  if (imp != imn) return nullptr;
  const auto pV0MC(static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(imp))); if (!pV0MC) return nullptr;
  if (((pV0MC->Y())<fCutMinV0Rap) || ((pV0MC->Y())>fCutMaxV0Rap)) return nullptr;

  const auto idvMC(pV0MC->GetPdgCode());
  const auto idp(pDauPosMC->GetPdgCode());
  const auto idn(pDauNegMC->GetPdgCode());
  auto bIsKshort((idp==211)  && (idn==-211)  && (idvMC== 310));
  auto bIsLambda((idp==2212) && (idn==-211)  && (idvMC== 3122));
  auto bIsAntiLa((idp==211)  && (idn==-2212) && (idvMC==-3122));
  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  UInt_t wsvMC(0);
  if (pV0MC->IsPrimary())                wsvMC |= AliPicoBase::kPrimary;
  if (pV0MC->IsPhysicalPrimary())        wsvMC |= AliPicoBase::kPhysicalPrimary;
  if (pV0MC->IsSecondaryFromWeakDecay()) wsvMC |= AliPicoBase::kSecondaryFromWeakDecay;
  if (pV0MC->IsSecondaryFromMaterial())  wsvMC |= AliPicoBase::kSecondaryFromMaterial;

  auto idmMC(0);
  UInt_t wsmMC(0);
  auto dMotherPt(0.);
  auto dMotherEta(0.);
  auto dMotherRap(0.);
  if (bIsLambda || bIsAntiLa) {
    const auto imv(pV0MC->GetMother());

    if (imv>=0) {
      const auto pMother(static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(imv)));

      if (pMother) {
        idmMC = pMother->GetPdgCode();
        if ((bIsLambda && ((idmMC== 3312) || (idmMC== 3322))) ||
            (bIsAntiLa && ((idmMC==-3312) || (idmMC==-3322)))) {
          dMotherPt  = pMother->Pt();
          dMotherEta = pMother->Eta();
          dMotherRap = pMother->Y();

          if (pMother->IsPrimary())                wsmMC |= AliPicoBase::kPrimary;
          if (pMother->IsPhysicalPrimary())        wsmMC |= AliPicoBase::kPhysicalPrimary;
          if (pMother->IsSecondaryFromWeakDecay()) wsmMC |= AliPicoBase::kSecondaryFromWeakDecay;
          if (pMother->IsSecondaryFromMaterial())  wsmMC |= AliPicoBase::kSecondaryFromMaterial;
        }
      }
    }
  }
//=============================================================================

  const auto dV0CosPA(pV0RD->CosPointingAngle(fPrimaryVtx));

  if (bIsKshort) if (dV0CosPA<fCutMinKshortCosPA) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if (dV0CosPA<fCutMinLambdaCosPA) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  auto dV0DistToPV(0.);
  for (auto i=0; i<3; ++i) dV0DistToPV +=  ((dV0Vtx[i]-fPrimaryVtx[i]) * (dV0Vtx[i]-fPrimaryVtx[i]));
  const auto dV0DistToPVoverP(TMath::Sqrt(dV0DistToPV) / (pV0RD->P()+1e-10));

  if (bIsKshort) if ((dV0DistToPVoverP*AliPicoBase::MassKshort())>fCutMaxKshortCtau) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if ((dV0DistToPVoverP*AliPicoBase::MassLambda())>fCutMaxLambdaCtau) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  const auto dV0ArmFrac(pV0RD->PtArmV0() / (TMath::Abs(pV0RD->AlphaV0())+1e-12));

  if (bIsKshort && (fCutMaxKshortArmFrac>0.)) if (dV0ArmFrac>fCutMaxKshortArmFrac) {
    bIsKshort = kFALSE;
  }

  if ((bIsLambda || bIsAntiLa) && fCutMaxLambdaArmFrac>0.) if (dV0ArmFrac>fCutMaxLambdaArmFrac) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  UInt_t wMask(0);
  if (bIsKshort) {
    TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, AliPicoBase::MassPion());
    TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, AliPicoBase::MassPion());

    const auto dKshortInvM((vPosPion+vNegPion).M());
    if ((dKshortInvM<(0.430006-0.0110029*dV0Pt)) ||
        (dKshortInvM>(0.563707+0.0114979*dV0Pt))) return nullptr;

    if (fCutMinKshortDeltaM>0.) {
      TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, AliPicoBase::MassProton());
      TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, AliPicoBase::MassProton());
      if ((TMath::Abs((vPosProton+vNegPion).M()-AliPicoBase::MassLambda())<fCutMinKshortDeltaM) ||
          (TMath::Abs((vNegProton+vPosPion).M()-AliPicoBase::MassLambda())<fCutMinKshortDeltaM)) return nullptr;
    }

    wMask = AliPicoBase::kKshort;
  }

  if (bIsLambda) {
    TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, AliPicoBase::MassProton());
    TLorentzVector vNegPion;   vNegPion.SetVectM(v3Neg, AliPicoBase::MassPion());

    const auto dLambdaInvM((vPosProton+vNegPion).M());
    if ((dLambdaInvM<(1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt))) ||
        (dLambdaInvM>(1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt)))) return nullptr;

    if (fCutMinLambdaDeletaM>0.) {
      TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, AliPicoBase::MassPion());
      if ((TMath::Abs((vPosPion+vNegPion).M()-AliPicoBase::MassKshort())<fCutMinLambdaDeletaM)) return nullptr;
    }

    wMask = AliPicoBase::kLambda;
  }

  if (bIsAntiLa) {
    TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, AliPicoBase::MassProton());
    TLorentzVector vPosPion;   vPosPion.SetVectM(v3Pos, AliPicoBase::MassPion());

    const auto dAntiLaInvM((vNegProton+vPosPion).M());
    if ((dAntiLaInvM<(1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt))) ||
        (dAntiLaInvM>(1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt)))) return nullptr;

    if (fCutMinLambdaDeletaM>0.) {
      TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, AliPicoBase::MassPion());
      if ((TMath::Abs((vPosPion+vNegPion).M()-AliPicoBase::MassKshort())<fCutMinLambdaDeletaM)) return nullptr;
    }

    wMask = AliPicoBase::kAntiLambda;
  }
//=============================================================================

  auto bPosInJC(kFALSE);
  auto bNegInJC(kFALSE);
  return (new AliPicoV0MC(wMask,
                          dV0Radius,
                          dV0CosPA,
                          dV0DistToPVoverP,
                          dDausDCA,
                          dPosDCAtoPV,
                          dNegDCAtoPV,
                          dDauXrowsTPC,
                          dDauXrowsOverFindableClusTPC,
                          v3Pos.Px(), v3Pos.Py(), v3Pos.Pz(),
                          v3Neg.Px(), v3Neg.Py(), v3Neg.Pz(),
                          bPosInJC, bNegInJC,
                          idvMC, wsvMC, pV0MC->Px(), pV0MC->Py(), pV0MC->Pz(), pV0MC->E(),
                          idmMC, wsmMC, dMotherPt, dMotherEta, dMotherRap));
}

//_____________________________________________________________________________
AliPicoV0MC *AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateMC(AliESDv0 const *pV0RD)
{
//
//  AliAnalysisTaskSEPicoV0Maker::SelectV0CandidateMC
//

  const auto pStack(MCEvent()->Stack()); if (!pStack) return nullptr;
//=============================================================================

  if (pV0RD->GetOnFlyStatus()) return nullptr;
  if (pV0RD->GetChi2V0()>fCutMaxV0Chi2) return nullptr;

  const auto dV0Pt(pV0RD->Pt()); if ((dV0Pt<fCutMinV0Pt) || (dV0Pt>fCutMaxV0Pt)) return nullptr;
//=============================================================================

  Double_t dV0Vtx[3];  pV0RD->GetXYZ(dV0Vtx[0], dV0Vtx[1], dV0Vtx[2]);
  const auto dV0Radius(TMath::Sqrt(dV0Vtx[0]*dV0Vtx[0] + dV0Vtx[1]*dV0Vtx[1]));
  if ((dV0Radius<fCutMinV0Radius) || (dV0Radius>fCutMaxV0Radius)) return nullptr;

  const auto dDausDCA(pV0RD->GetDcaV0Daughters()); if (dDausDCA>fCutMaxDausDCA) return nullptr;
//=============================================================================

  const auto nPosIndex(TMath::Abs(pV0RD->GetPindex())); if (nPosIndex<0) return nullptr;
  const auto nNegIndex(TMath::Abs(pV0RD->GetNindex())); if (nNegIndex<0) return nullptr;

  auto pDauPosRD(fEventESD->GetTrack(nPosIndex)); if (!pDauPosRD) return nullptr;
  auto pDauNegRD(fEventESD->GetTrack(nNegIndex)); if (!pDauNegRD) return nullptr;

  const auto dMegField(fEventESD->GetMagneticField());
  const auto dPosDCAtoPV(TMath::Abs(pDauPosRD->GetD(fPrimaryVtx[0],fPrimaryVtx[1],dMegField)));
  if (dPosDCAtoPV<fCutMinDauDCAtoPV) return nullptr;

  const auto dNegDCAtoPV(TMath::Abs(pDauNegRD->GetD(fPrimaryVtx[0],fPrimaryVtx[1],dMegField)));
  if (dNegDCAtoPV<fCutMinDauDCAtoPV) return nullptr;
//=============================================================================

  if (!(pDauPosRD->GetStatus() & AliESDtrack::kTPCrefit)) return nullptr;
  if (!(pDauNegRD->GetStatus() & AliESDtrack::kTPCrefit)) return nullptr;
  if ((pDauPosRD->GetKinkIndex(0)>0) || (pDauNegRD->GetKinkIndex(0)>0)) return nullptr;

  const auto dPosXrowsTPC(pDauPosRD->GetTPCClusterInfo(2,1));
  const auto dNegXrowsTPC(pDauNegRD->GetTPCClusterInfo(2,1));
  const auto dDauXrowsTPC((dPosXrowsTPC<dNegXrowsTPC) ? dPosXrowsTPC : dNegXrowsTPC);
  if (dDauXrowsTPC<fCutMinDauXrowsTPC) return nullptr;

  const auto wPosTPCNClsF(pDauPosRD->GetTPCNclsF()); if (wPosTPCNClsF<=0) return nullptr;
  const auto wNegTPCNClsF(pDauNegRD->GetTPCNclsF()); if (wNegTPCNClsF<=0) return nullptr;
  const auto dPosXrowsOverFindableClusTPC( ((Double_t)dPosXrowsTPC) / ((Double_t)wPosTPCNClsF) );
  const auto dNegXrowsOverFindableClusTPC( ((Double_t)dNegXrowsTPC) / ((Double_t)wNegTPCNClsF) );

  const auto dDauXrowsOverFindableClusTPC((dPosXrowsOverFindableClusTPC<dNegXrowsOverFindableClusTPC) ?
                                           dPosXrowsOverFindableClusTPC :
                                           dNegXrowsOverFindableClusTPC);
  if (dDauXrowsOverFindableClusTPC<fCutMinDauXrowsOverFindableClusTPC) return nullptr;
//=============================================================================

  const auto nPosCharge(pDauPosRD->Charge());
  const auto nNegCharge(pDauNegRD->Charge());
  if ((nPosCharge==0) || (nNegCharge==0) || (nPosCharge==nNegCharge)) return nullptr;

  Double_t dPosPxPyPz[3] = { 0., 0., 0. };
  Double_t dNegPxPyPz[3] = { 0., 0., 0. };
  if ((nPosCharge<0) && (nNegCharge>0)) {
    pDauPosRD = fEventESD->GetTrack(nNegIndex);
    pDauNegRD = fEventESD->GetTrack(nPosIndex);

    pV0RD->GetNPxPyPz(dPosPxPyPz[0], dPosPxPyPz[1], dPosPxPyPz[2]);
    pV0RD->GetPPxPyPz(dNegPxPyPz[0], dNegPxPyPz[1], dNegPxPyPz[2]);
  } else {
    pV0RD->GetPPxPyPz(dPosPxPyPz[0], dPosPxPyPz[1], dPosPxPyPz[2]);
    pV0RD->GetNPxPyPz(dNegPxPyPz[0], dNegPxPyPz[1], dNegPxPyPz[2]);
  }

  const TVector3 v3Pos(dPosPxPyPz);
  const TVector3 v3Neg(dNegPxPyPz);
  if ((v3Pos.Pt()<fCutMinDauPt) || (v3Neg.Pt()<fCutMinDauPt)) return nullptr;
  const auto dPosEta(v3Pos.Eta()); if ((dPosEta<fCutMinDauEta) || (dPosEta>fCutMaxDauEta)) return nullptr;
  const auto dNegEta(v3Neg.Eta()); if ((dNegEta<fCutMinDauEta) || (dNegEta>fCutMaxDauEta)) return nullptr;
//=============================================================================

  AliMCParticle *pDauTmpMC(nullptr);
  const auto inp(TMath::Abs(pDauPosRD->GetLabel())); if (inp<0) return nullptr;
  pDauTmpMC = static_cast<AliMCParticle*>(MCEvent()->GetTrack(inp)); if (!pDauTmpMC) return nullptr;
  const auto pDauPosMC(pDauTmpMC->Particle()); if (!pDauPosMC) return nullptr;
  const auto imp(pDauPosMC->GetFirstMother()); if (imp<0) return nullptr;

  const auto inn(TMath::Abs(pDauNegRD->GetLabel())); if (inn<0) return nullptr;
  pDauTmpMC = static_cast<AliMCParticle*>(MCEvent()->GetTrack(inn)); if (!pDauTmpMC) return nullptr;
  const auto pDauNegMC(pDauTmpMC->Particle()); if (!pDauNegMC) return nullptr;
  const auto imn(pDauNegMC->GetFirstMother()); if (imn<0) return nullptr;

  if (imp != imn) return nullptr;
  pDauTmpMC = static_cast<AliMCParticle*>(MCEvent()->GetTrack(imp)); if (!pDauTmpMC) return nullptr;
  const auto pV0MC(pDauTmpMC->Particle()); if (!pV0MC) return nullptr;
  if (((pV0MC->Y())<fCutMinV0Rap) || ((pV0MC->Y())>fCutMaxV0Rap)) return nullptr;

  const auto idvMC(pV0MC->GetPdgCode());
  const auto idp(pDauPosMC->GetPdgCode());
  const auto idn(pDauNegMC->GetPdgCode());
  auto bIsKshort((idp==211)  && (idn==-211)  && (idvMC== 310));
  auto bIsLambda((idp==2212) && (idn==-211)  && (idvMC== 3122));
  auto bIsAntiLa((idp==211)  && (idn==-2212) && (idvMC==-3122));
  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  UInt_t wsvMC(0);
  if (imp<pStack->GetNprimary())             wsvMC |= AliPicoBase::kPrimary;
  if (pStack->IsPhysicalPrimary(imp))        wsvMC |= AliPicoBase::kPhysicalPrimary;
  if (pStack->IsSecondaryFromWeakDecay(imp)) wsvMC |= AliPicoBase::kSecondaryFromWeakDecay;
  if (pStack->IsSecondaryFromMaterial(imp))  wsvMC |= AliPicoBase::kSecondaryFromMaterial;

  auto   idmMC(0);
  UInt_t wsmMC(0);
  auto dMotherPt(0.);
  auto dMotherEta(0.);
  auto dMotherRap(0.);
  if (bIsLambda || bIsAntiLa) {
    const auto imv(pV0MC->GetFirstMother());

    if (imv>=0) {
      pDauTmpMC = static_cast<AliMCParticle*>(MCEvent()->GetTrack(imv));

      if (pDauTmpMC) {
        auto pMother(pDauTmpMC->Particle());

        if (pMother) {
          idmMC = pMother->GetPdgCode();
          if ((bIsLambda && ((idmMC== 3312) || (idmMC== 3322))) ||
              (bIsAntiLa && ((idmMC==-3312) || (idmMC==-3322)))) {
            dMotherPt  = pMother->Pt();
            dMotherEta = pMother->Eta();
            dMotherRap = pMother->Y();

            if (imp<pStack->GetNprimary())             wsmMC |= AliPicoBase::kPrimary;
            if (pStack->IsPhysicalPrimary(imv))        wsmMC |= AliPicoBase::kPhysicalPrimary;
            if (pStack->IsSecondaryFromWeakDecay(imv)) wsmMC |= AliPicoBase::kSecondaryFromWeakDecay;
            if (pStack->IsSecondaryFromMaterial(imv))  wsmMC |= AliPicoBase::kSecondaryFromMaterial;
          }
        }
      }
    }
  }
//=============================================================================

  const auto dV0CosPA(pV0RD->GetV0CosineOfPointingAngle(fPrimaryVtx[0],fPrimaryVtx[1],fPrimaryVtx[2]));

  if (bIsKshort) if (dV0CosPA<fCutMinKshortCosPA) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if (dV0CosPA<fCutMinLambdaCosPA) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  auto dV0DistToPV(0.);
  for (auto i=0; i<3; ++i) dV0DistToPV +=  ((dV0Vtx[i]-fPrimaryVtx[i]) * (dV0Vtx[i]-fPrimaryVtx[i]));
  const auto dV0DistToPVoverP(TMath::Sqrt(dV0DistToPV) / (pV0RD->P()+1e-10));

  if (bIsKshort) if ((dV0DistToPVoverP*AliPicoBase::MassKshort())>fCutMaxKshortCtau) {
    bIsKshort = kFALSE;
  }

  if (bIsLambda || bIsAntiLa) if ((dV0DistToPVoverP*AliPicoBase::MassLambda())>fCutMaxLambdaCtau) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  const auto dV0ArmFrac(pV0RD->PtArmV0() / (TMath::Abs(pV0RD->AlphaV0())+1e-12));

  if (bIsKshort && (fCutMaxKshortArmFrac>0.)) if (dV0ArmFrac>fCutMaxKshortArmFrac) {
    bIsKshort = kFALSE;
  }

  if ((bIsLambda && bIsAntiLa) && (fCutMaxLambdaArmFrac>0.)) if (dV0ArmFrac>fCutMaxLambdaArmFrac) {
    bIsLambda = kFALSE;
    bIsAntiLa = kFALSE;
  }

  if (!(bIsKshort || bIsLambda || bIsAntiLa)) return nullptr;
//=============================================================================

  UInt_t wMask(0);
  if (bIsKshort) {
    TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, AliPicoBase::MassPion());
    TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, AliPicoBase::MassPion());

    const auto dKshortInvM((vPosPion+vNegPion).M());
    if ((dKshortInvM<(0.430006 - 0.0110029*dV0Pt)) ||
        (dKshortInvM>(0.563707 + 0.0114979*dV0Pt))) return nullptr;

    if (fCutMinKshortDeltaM>0.) {
      TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, AliPicoBase::MassProton());
      TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, AliPicoBase::MassProton());
      if ((TMath::Abs((vPosProton+vNegPion).M()-AliPicoBase::MassLambda())<fCutMinKshortDeltaM) ||
          (TMath::Abs((vNegProton+vPosPion).M()-AliPicoBase::MassLambda())<fCutMinKshortDeltaM)) return nullptr;
    }

    wMask = AliPicoBase::kKshort;
  }

  if (bIsLambda) {
    TLorentzVector vPosProton; vPosProton.SetVectM(v3Pos, AliPicoBase::MassProton());
    TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, AliPicoBase::MassPion());

    const auto dLambdaInvM((vPosProton+vNegPion).M());
    if ((dLambdaInvM<(1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt))) ||
        (dLambdaInvM>(1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt)))) return nullptr;

    if (fCutMinLambdaDeletaM>0.) {
      TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, AliPicoBase::MassPion());
      if ((TMath::Abs((vPosPion+vNegPion).M()-AliPicoBase::MassKshort())<fCutMinLambdaDeletaM)) return nullptr;
    }

    wMask = AliPicoBase::kLambda;
  }

  if (bIsAntiLa) {
    TLorentzVector vNegProton; vNegProton.SetVectM(v3Neg, AliPicoBase::MassProton());
    TLorentzVector vPosPion; vPosPion.SetVectM(v3Pos, AliPicoBase::MassPion());

    const auto dAntiLaInvM((vNegProton+vPosPion).M());
    if ((dAntiLaInvM<(1.09501 - 0.00523272*dV0Pt - 0.075269*TMath::Exp(-3.46339*dV0Pt))) ||
        (dAntiLaInvM>(1.13688 + 0.00527838*dV0Pt + 0.084222*TMath::Exp(-3.80595*dV0Pt)))) return nullptr;

    if (fCutMinLambdaDeletaM>0.) {
      TLorentzVector vNegPion; vNegPion.SetVectM(v3Neg, AliPicoBase::MassPion());
      if ((TMath::Abs((vPosPion+vNegPion).M()-AliPicoBase::MassKshort())<fCutMinLambdaDeletaM)) return nullptr;
    }

    wMask = AliPicoBase::kAntiLambda;
  }
//=============================================================================

  auto bPosInJC(kFALSE);
  auto bNegInJC(kFALSE);
  return (new AliPicoV0MC(wMask,
                          dV0Radius,
                          dV0CosPA,
                          dV0DistToPVoverP,
                          dDausDCA,
                          dPosDCAtoPV,
                          dNegDCAtoPV,
                          dDauXrowsTPC,
                          dDauXrowsOverFindableClusTPC,
                          v3Pos.Px(), v3Pos.Py(), v3Pos.Pz(),
                          v3Neg.Px(), v3Neg.Py(), v3Neg.Pz(),
                          bPosInJC, bNegInJC,
                          idvMC, wsvMC, pV0MC->Px(), pV0MC->Py(), pV0MC->Pz(), pV0MC->Energy(),
                          idmMC, wsmMC, dMotherPt, dMotherEta, dMotherRap));
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEPicoV0Maker::IsEventNotAcpt()
{
//
//  AliAnalysisTaskSEPicoV0Maker::IsEventNotAcpt
//

  fEventAcptMask = 0;
  if (!InputEvent())  return kTRUE;
  if (!fInputHandler) return kTRUE;

  fEventAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fEventESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if ((!fEventAOD) && (!fEventESD)) return kTRUE;
//=============================================================================

  fRespoPID = fInputHandler->GetPIDResponse();
  if (!fRespoPID) return kTRUE;
//=============================================================================

  if (fIsAnaUseMC) {
    if (MCEvent()) {
      if (MCEvent()->GetNumberOfTracks()<=0) return kTRUE;
    } else {
      return kTRUE;
    }

    const auto pHeader(MCEvent()->Header()); if (!pHeader) return kTRUE;

    if (fIsDPMjetMC) {
      const auto pDPMjetH(dynamic_cast<AliGenDPMjetEventHeader*>(pHeader->GenEventHeader()));

      if (pDPMjetH) {
        auto nd0(0), nd1(0), nd2(0);
        pDPMjetH->GetNDiffractive(nd1, nd2, nd0);
        if ((nd1+nd2) != (pDPMjetH->ProjectileParticipants() + pDPMjetH->TargetParticipants())) return kTRUE;
      }
    }
  }

  fEventAcptMask |= AliPicoBase::kEventCheck;
//=============================================================================

  if ((fMultEsti.GetEntries()>0) || (!fMultEstDef.IsNull())) {
    Float_t dMult(-999.);

    if (fUseMultOld) {
      auto pCentSel(InputEvent()->GetCentrality());
      if (!pCentSel) { fEventAcptMask=0; return kTRUE; }
      if (pCentSel->GetQuality()!=0) return kFALSE;

      if (fMultEsti.GetEntries()>0) {
        TObjString *ps(nullptr);
        const auto next(fMultEsti.MakeIterator());
        while ((ps = static_cast<TObjString*>((*next)()))) {
          const auto s(ps->String());
          const auto p(static_cast<TParameter<Float_t>*>((fMultEsti(s.Data()))));
          if (p) p->SetVal(pCentSel->GetCentralityPercentile(s.Data()));
        }
      }

      if (!fMultEstDef.IsNull()) dMult = pCentSel->GetCentralityPercentile(fMultEstDef.Data());
    } else {
      auto pMultSel(static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection")));
      if (!pMultSel) { fEventAcptMask=0; return kTRUE; }

      if (fMultEsti.GetEntries()>0) {
        TObjString *ps(nullptr);
        const auto next(fMultEsti.MakeIterator());
        while (((ps = static_cast<TObjString*>((*next)())))) {
          const auto s(ps->String());
          const auto p(static_cast<TParameter<Float_t>*>(fMultEsti(s.Data())));
          if (p) p->SetVal(pMultSel->GetMultiplicityPercentile(s.Data()));
        }
      }

      if (!fMultEstDef.IsNull()) dMult = pMultSel->GetMultiplicityPercentile(fMultEstDef.Data());
    }

    if (!fMultEstDef.IsNull()) if ((dMult<fCutMinMult) || (dMult>=fCutMaxMult)) return kFALSE;
  }

  fEventAcptMask |= AliPicoBase::kEventMult;
//=============================================================================

  const UInt_t wMask(fInputHandler->IsEventSelected());
  if ((wMask & fTriggerMask) != fTriggerMask) return kFALSE;
  if (fIsSkipFastOnly) if ((wMask & AliVEvent::kFastOnly) == AliVEvent::kFastOnly) return kFALSE;

  fEventAcptMask |= AliPicoBase::kEventTrigger;
//=============================================================================

  const auto pVertex(InputEvent()->GetPrimaryVertex());
  if (!pVertex) return kFALSE;
  pVertex->GetXYZ(fPrimaryVtx);

  if (fUseAnaUtils) {
    auto pUtils(new AliAnalysisUtils());
    if (!pUtils->IsVertexSelected2013pA(InputEvent()))  return kFALSE;
    if (pUtils->IsSPDClusterVsTrackletBG(InputEvent())) return kFALSE;
    if (pUtils->IsPileUpEvent(InputEvent())) return kFALSE;

    if ((fCollisionType==(AliPicoBase::kPA)) ||
        (fCollisionType==(AliPicoBase::kAP))) {
      if (pUtils->IsFirstEventInChunk(InputEvent())) return kFALSE;
    }
  }

/*if ((fCollisionType==(AliPicoBase::kPA)) || (fCollisionType==(AliPicoBase::kAP))) {
    if (fEventAOD) {
      const AliAODVertex *pVtxSPD = fEventAOD->GetPrimaryVertexSPD();
      const AliAODVertex *pVtxTrk = fEventAOD->GetPrimaryVertex();
      if ((!pVtxSPD) && (!pVtxTrk)) return (fEventAcptMask==0);
    }

    if (fEventESD) {
      Bool_t fHasVertex = kFALSE;
      const AliESDVertex *pVtxESD = fEventESD->GetPrimaryVertexTracks();
      if (pVtxESD->GetNContributors()<1) {
        pVtxESD = fEventESD->GetPrimaryVertexSPD();
        if (pVtxESD->GetNContributors()<1) fHasVertex = kFALSE;
        else fHasVertex = kTRUE;

        Double_t cov[6] = { 0., 0., 0., 0., 0., 0. };
        pVtxESD->GetCovarianceMatrix(cov);
        Double_t zRes = TMath::Sqrt(cov[5]);
        if (pVtxESD->IsFromVertexerZ() && (zRes>0.25)) fHasVertex = kFALSE;
      } else fHasVertex = kTRUE;

      if (!fHasVertex) return (fEventAcptMask==0);
    }
  } else {
    if (fEventAOD) {
      const auto pVtxSPD(fEventAOD->GetPrimaryVertexSPD()); if (!pVtxSPD) return (fEventAcptMask==0);
      const auto pVtxTrk(fEventAOD->GetPrimaryVertex());    if (!pVtxTrk) return (fEventAcptMask==0);
    }

    if (fEventESD) {
      const auto pVtxPri(fEventESD->GetPrimaryVertex());       if (!pVtxPri) return (fEventAcptMask==0);
      const auto pVtxSPD(fEventESD->GetPrimaryVertexSPD());    if (!pVtxSPD) return (fEventAcptMask==0);
      const auto pVtxTrk(fEventESD->GetPrimaryVertexTracks()); if (!pVtxTrk) return (fEventAcptMask==0);
      if ((!(pVtxPri->GetStatus())) &&
          (!(pVtxSPD->GetStatus())) &&
          (!(pVtxTrk->GetStatus()))) return (fEventAcptMask==0);
    }
  }*/

  fEventAcptMask |= AliPicoBase::kEventVertex;
//=============================================================================

  if (fIsRefitV0sESD && fEventESD) {
    Double_t dCuts[7] = { fCutMaxV0Chi2,
                          fCutMinDauDCAtoPV,
                          fCutMinDauDCAtoPV,
                          fCutMaxDausDCA,
                          fCutMinKshortCosPA,
                          fCutMinV0Radius,
                          fCutMaxV0Radius };

    fEventESD->ResetV0s();
    AliV0vertexer aV0vtxer;
    aV0vtxer.SetDefaultCuts(dCuts);
    aV0vtxer.Tracks2V0vertices(fEventESD);
  }
//=============================================================================

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEPicoV0Maker::IsEventNotINEL()
{
//
//  AliAnalysisTaskSEPicoV0Maker::IsEventNotINEL
//

  if ((fEventAcptMask & AliPicoBase::kEventCheck) != AliPicoBase::kEventCheck) return kTRUE;
  if ((fEventAcptMask & AliPicoBase::kEventMult)  != AliPicoBase::kEventMult)  return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEPicoV0Maker::IsEventNotMBsa()
{
//
//  AliAnalysisTaskSEPicoV0Maker::IsEventNotMBsa
//

  if (IsEventNotINEL()) return kTRUE;
  if ((fEventAcptMask & AliPicoBase::kEventTrigger) != AliPicoBase::kEventTrigger) return kTRUE;
  if ((fEventAcptMask & AliPicoBase::kEventVertex)  != AliPicoBase::kEventVertex)  return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::FillHistogramsEH()
{
//
//  AliAnalysisTaskSEPicoV0Maker::FillHistogramsEH
//

  const auto h(static_cast<TH1D*>(fOutputListEH->FindObject("hEventCheck")));

  if (fEventAcptMask==0) h->Fill(0.);
  if ((fEventAcptMask & AliPicoBase::kEventCheck)   == AliPicoBase::kEventCheck)   h->Fill(1.);
  if ((fEventAcptMask & AliPicoBase::kEventMult)    == AliPicoBase::kEventMult)    h->Fill(2.);
  if ((fEventAcptMask & AliPicoBase::kEventTrigger) == AliPicoBase::kEventTrigger) h->Fill(3.);
  if ((fEventAcptMask & AliPicoBase::kEventVertex)  == AliPicoBase::kEventVertex)  h->Fill(4.);

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::FillHistogramsMC()
{
//
//  AliAnalysisTaskSEPicoV0Maker::FillHistogramsMC
//

  auto nPrimary(0);
  AliStack *pStack(nullptr);

  if (fEventESD) {
    pStack   = MCEvent()->Stack(); if (!pStack) return;
    nPrimary = pStack->GetNprimary();
  }
//=============================================================================

  const auto n(6+fMultEsti.GetEntries());
  const auto dEvType(IsEventNotMBsa() ? 1.5 : 0.5);
  auto hsV0(static_cast<THnSparseD*>(fOutputListMC->FindObject("hsV0")));
  auto hsXi(static_cast<THnSparseD*>(fOutputListMC->FindObject("hsXi")));

  if (!(hsV0 && hsXi)) {
    AliFatal("Cannot find hsV0 or hsXi; should not happen");
    return;
  }
//=============================================================================

  for (auto i=0; i<MCEvent()->GetNumberOfTracks(); ++i) {
    TParticle *pESD(nullptr);
    AliAODMCParticle *pAOD(nullptr);

    if (fEventAOD) {
      pAOD = static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(i));
      if (!pAOD) continue;
    }

    if (fEventESD) {
      const auto pMC(static_cast<AliMCParticle*>(MCEvent()->GetTrack(i))); if (!pMC) continue;
      pESD = pMC->Particle(); if (!pESD) continue;
    }
//=============================================================================

    const auto bPri(pAOD ? pAOD->IsPrimary() : (i<nPrimary));
    const auto bPhy(pAOD ? pAOD->IsPhysicalPrimary() : pStack->IsPhysicalPrimary(i));
    if (!(bPri || bPhy)) continue;
//=============================================================================

    const auto id(pAOD ? pAOD->GetPdgCode() : pESD->GetPdgCode());

    const auto bXi(bPri && ((id==3312) || (id==-3312)));
    const auto bV0(bPhy && ((id==3122) || (id==-3122) || (id==310)));
    if (!(bXi || bV0)) continue;
//=============================================================================

    const auto dEta(pAOD ? pAOD->Eta() : pESD->Eta());
    if ((dEta<-5.) || (dEta>=5.)) continue;

    const auto dRapLab(pAOD ? pAOD->Y() : pESD->Y());
    if ((dRapLab<-5.) || (dRapLab>=5.)) continue;

    const auto dRapCMS(dRapLab + fRapidityShift);
    if ((dRapCMS<-5.) || (dRapCMS>=5.)) continue;
//=============================================================================

    Double_t dVar[n];
    dVar[0] = (pAOD ? pAOD->Pt() : pESD->Pt());

    dVar[1] = dEta;
    dVar[2] = dRapLab;
    dVar[3] = dRapCMS;
    dVar[4] = dEvType;

    auto l(6);
    TObjString *ps(nullptr);
    const auto next(fMultEsti.MakeIterator());
    while ((ps = static_cast<TObjString*>((*next)()))) {
      const auto s(ps->String());
      const auto p(static_cast<TParameter<Float_t>*>(fMultEsti(s.Data())));
      if (p) dVar[l++] = p->GetVal();
    }

    if (bV0) {
      if (id== 310 ) dVar[5] = 0.;
      if (id== 3122) dVar[5] = 1.;
      if (id==-3122) dVar[5] = 2.;
      hsV0->Fill(dVar);
    }

    if (bXi) {
      if (id== 3312) dVar[5] = 0.;
      if (id==-3312) dVar[5] = 1.;
      hsXi->Fill(dVar);
    }
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::CreateHistogramsEH()
{
//
//  AliAnalysisTaskSEPicoV0Maker::CreateHistogramsEH
//

  const auto b(TH1::AddDirectoryStatus());
  TH1::AddDirectory(kFALSE);
//=============================================================================

  fOutputListEH->Add(new TH1D("hEventCheck", "", 5, -0.5, 4.5));
  fOutputListEH->Add(new TH2D("hKshortPtInvM", "", 1000, 0., 100., 300, AliPicoBase::MassKshort()-0.15,
                                                                        AliPicoBase::MassKshort()+0.15));

  fOutputListEH->Add(new TH2D("hLambdaPtInvM", "", 1000, 0., 100., 200, AliPicoBase::MassLambda()-0.10,
                                                                        AliPicoBase::MassLambda()+0.10));

  fOutputListEH->Add(new TH2D("hAntiLaPtInvM", "", 1000, 0., 100., 200, AliPicoBase::MassLambda()-0.10,
                                                                        AliPicoBase::MassLambda()+0.10));

  TH1 *h(nullptr);
  TListIter next(fOutputListEH);
  while ((h = static_cast<TH1*>(next()))) h->Sumw2();

  TH1::AddDirectory(b);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::CreateHistogramsMC()
{
//
//  AliAnalysisTaskSEPicoV0Maker::CreateHistogramsMC
//

  const auto b(TH1::AddDirectoryStatus());
  TH1::AddDirectory(kFALSE);
//=============================================================================

  const auto nr(5); // 0: Pt
                    // 1: eta
                    // 2: rap in Lab
                    // 3: rap in CMS
                    // 4: Event type
                    //    ==0.5, INEL
                    //    ==1.5, MB
                    // 5: particle type
                    //   V0
                    //     ==0, Kshort
                    //     ==1, Lambda
                    //     ==2, AntiLa
                    //   Xi
                    //     ==0, XiNeg
                    //     ==1, XiPos
  const Int_t    nBin[nr] = { 1000, 100, 100, 100, 2  };
  const Double_t dMin[nr] = {   0., -5., -5., -5., 0. };
  const Double_t dMax[nr] = { 100.,  5.,  5.,  5., 2. };
  const auto ns(1 + nr + fMultEsti.GetEntries());
//=============================================================================

  Int_t nV0Bin[ns], nXiBin[ns];
  Double_t dV0Min[ns], dV0Max[ns];
  Double_t dXiMin[ns], dXiMax[ns];

  for (auto i=0; i<ns; ++i) {
    if (i<nr) {
      nV0Bin[i] = nXiBin[i] = nBin[i];
      dV0Min[i] = dXiMin[i] = dMin[i];
      dV0Max[i] = dXiMax[i] = dMax[i];
    }

    if (i==nr) {
      nV0Bin[i] = 3; dV0Min[i] = -0.5; dV0Max[i] = 2.5;
      nXiBin[i] = 2; dXiMin[i] = -0.5; dXiMax[i] = 1.5;
    }

    if (i>nr) {
      nV0Bin[i] = nXiBin[i] =  110;
      dV0Min[i] = dXiMin[i] =  -5.;
      dV0Max[i] = dXiMax[i] = 105.;
    }
  }
//=============================================================================

  const TString sa[nr+1] { "pT", "eta", "y_lab", "y_cms", "evt_t", "par_t" };
  fOutputListMC->Add(new THnSparseD("hsV0", "", ns, nV0Bin, dV0Min, dV0Max));
  fOutputListMC->Add(new THnSparseD("hsXi", "", ns, nXiBin, dXiMin, dXiMax));

  TObjString *ps(nullptr);
  THnSparseD *hs(nullptr);
  TListIter next(fOutputListMC);
  const auto pn(fMultEsti.MakeIterator());
  while ((hs = static_cast<THnSparseD*>(next()))) {
    for (auto i=0; i<nr; ++i) hs->GetAxis(i)->SetName(sa[i].Data());

    auto l(nr);
    while ((ps = static_cast<TObjString*>((*pn)()))) {
      const auto s(ps->String());
      const auto p(hs->GetAxis(l++));
      if (p) p->SetName(s.Data());
    }
  }
//=============================================================================

  TH1::AddDirectory(b);

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::InitAnalysis()
{
//
//  AliAnalysisTaskSEPicoV0Maker::InitAnalysis
//

  if (fCollisionType==(AliPicoBase::kPP)) InitParamsPP();
  if (fCollisionType==(AliPicoBase::kPA)) InitParamsPA();
  if (fCollisionType==(AliPicoBase::kAP)) InitParamsAP();
  if (fCollisionType==(AliPicoBase::kAA)) InitParamsAA();

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::InitParamsPP()
{
//
//  AliAnalysisTaskSEPicoV0Maker::InitParametersPP
//

  fRapidityShift = 0.;
//=============================================================================

  fCutMaxV0Chi2   = 33.;
  fCutMinV0Radius = 0.3;   // default: >0.5; uncertainty: >0.3, 0.4, 0.6, 0.7;
  fCutMaxV0Radius = 200.;  // default: not applied in pp
//=============================================================================

  fCutMaxDausDCA                     = 1.5;  // default: <1;    uncertainty: <0.5, 0.75, 1.25, 1.5
  fCutMinDauDCAtoPV                  = 0.05; // default: >0.06; uncertainty: >0.05, 0.055, 0.07, 0.08
  fCutMinDauXrowsTPC                 = 70.;  // default: >70;   uncertainty: >75, 80
  fCutMinDauXrowsOverFindableClusTPC = 0.8;  // default: >0.8;  uncertainty: >0.95
//=============================================================================

  fCutMaxKshortSigmaTPC = -1.;   // default: <5;     uncertainty: w/o cut
  fCutMinKshortCosPA    = 0.95;  // default: >0.97;  uncertainty: >0.95, 0.96, 0.98, 0.99
  fCutMaxKshortCtau     = 30.;   // default: <20;    uncertainty: <12, 30
  fCutMaxKshortArmFrac  = -1.;   // default: not applied in pp
  fCutMinKshortDeltaM   = 0.003; // default: >0.005; uncertainty: >0.003, 0.006
//=============================================================================

  fCutMaxLambdaSigmaTPC  = 7.;    // default: <5;     uncertainty: 4, 6, 7
  fCutMinLambdaCosPA     = 0.993; // default: >0.995; uncertainty: >0.993, 0.994, 0.996, 0.997
  fCutMaxLambdaCtau      = 40.;   // default: <30;    uncertainty: <20, 40
  fCutMaxLambdaArmFrac   = -1.;   // default: not applied in pp
  fCutMinLambdaDeletaM   = -1.;   // default: >0.01;  uncertainty: w/o rejection

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::InitParamsPA()
{
//
//  AliAnalysisTaskSEPicoV0Maker::InitParametersPA
//

  InitParamsPP();

  fRapidityShift = 0.465;

  fCutMaxLambdaSigmaTPC = 6;  // default: <5; uncertaity: <4, 6

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::InitParamsAP()
{
//
//  AliAnalysisTaskSEPicoV0Maker::InitParametersAP
//

  InitParamsPP();

  fRapidityShift = -0.465;

  fCutMaxLambdaSigmaTPC = 6;  // default: <5; uncertaity: <4, 6

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0Maker::InitParamsAA()
{
//
//  AliAnalysisTaskSEPicoV0Maker::InitParametersAA
//

  InitParamsPP();

  fCutMinV0Radius = 0.9;  // default: 5; uncertainty: varying around
  fCutMaxV0Radius = 100.;

//fCutMinDauPt      = 0.16;  // default: 0.16
  fCutMinDauDCAtoPV = 0.08;  // default: 0.1; uncertainty: 0.08, 0.12

  fCutMaxKshortSigmaTPC = 6;      // default: <5;     uncertainty: <4, 6;
  fCutMinKshortCosPA    = 0.997;  // default: >0.998; uncertainty: 0.997, 0.999
  fCutMaxKshortArmFrac  = 0.2;    // default: <0.2

  fCutMaxLambdaSigmaTPC = 6;      // default: <5;     uncertaity: <4, 6
  fCutMinLambdaCosPA    = 0.997;  // default: >0.998; uncertainty: 0.997, 0.999
  fCutMinLambdaDeletaM  = 0.008;  // default: >0.01;  uncertainty: 0.008, 0.012

  return;
}
