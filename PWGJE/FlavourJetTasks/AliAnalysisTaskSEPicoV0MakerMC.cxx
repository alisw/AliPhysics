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
#include <TString.h>
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
#include "AliPicoV0MC.h"
#include "AliAnalysisTaskSEPicoV0MakerMC.h"

ClassImp(AliAnalysisTaskSEPicoV0MakerMC)
//=============================================================================

AliAnalysisTaskSEPicoV0MakerMC::AliAnalysisTaskSEPicoV0MakerMC() :
AliAnalysisTaskSE(),
fTriggerMask(0),
fCollisionType(0),
fUseAnaUtils(kFALSE),
fIsDPMjetMC(kFALSE),
fMultEst(""),
fMultMin(0.),
fMultMax(0.),
fMultOld(kFALSE),
fIsSkipFastOnly(kFALSE),
fIsRefitV0sESD(kFALSE),
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
fEventAOD(0),
fEventESD(0),
fRespoPID(0),
fEventAcptMask(0),
fPicoV0sClArr(0),
fListUserOutputs(0)
{
//
// Default constructor
//

  for (auto &d : fPrimaryVtx) d = -999.;
}

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0MakerMC::AliAnalysisTaskSEPicoV0MakerMC(const char *name) :
AliAnalysisTaskSE(name),
fTriggerMask(AliVEvent::kAny),
fCollisionType(AliPicoBase::kPP),
fUseAnaUtils(kFALSE),
fIsDPMjetMC(kFALSE),
fMultEst(""),
fMultMin(-99999.),
fMultMax(999999.),
fMultOld(kFALSE),
fIsSkipFastOnly(kFALSE),
fIsRefitV0sESD(kFALSE),
fCutMinV0Pt(0.),
fCutMaxV0Pt(100.),
fCutMinV0Rap(-10.),
fCutMaxV0Rap(10.),
fCutMinDauPt(0.),
fCutMinDauEta(-10.),
fCutMaxDauEta(10.),
fCutMaxV0Chi2(33.),
fCutMinV0Radius(0.5),
fCutMaxV0Radius(200.),
fCutMaxDausDCA(1.),
fCutMinDauDCAtoPV(0.06),
fCutMinDauXrowsTPC(70.),
fCutMinDauXrowsOverFindableClusTPC(0.8),
fCutMaxKshortSigmaTPC(5.),
fCutMinKshortCosPA(0.97),
fCutMaxKshortCtau(20.),
fCutMaxKshortArmFrac(-1.),
fCutMinKshortDeltaM(0.005),
fCutMaxLambdaSigmaTPC(5.),
fCutMinLambdaCosPA(0.995),
fCutMaxLambdaCtau(30.),
fCutMaxLambdaArmFrac(-1.),
fCutMinLambdaDeletaM(0.01),
fEventAOD(0),
fEventESD(0),
fRespoPID(0),
fEventAcptMask(0),
fPicoV0sClArr(0),
fListUserOutputs(0)
{
//
// Constructor
//

  for (auto &d : fPrimaryVtx) d = -999.;
  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskSEPicoV0MakerMC::~AliAnalysisTaskSEPicoV0MakerMC()
{
//
// Default destructor
//

  if (fEventAOD) { delete fEventAOD; fEventAOD = 0; }
  if (fEventESD) { delete fEventESD; fEventESD = 0; }
  if (fRespoPID) { delete fRespoPID; fRespoPID = 0; }

  if (fPicoV0sClArr)    { delete fPicoV0sClArr;    fPicoV0sClArr    = 0; }
  if (fListUserOutputs) { delete fListUserOutputs; fListUserOutputs = 0; }
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0MakerMC::Init()
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::Init
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0MakerMC::UserCreateOutputObjects()
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::UserCreateOutputObjects
//

  InitAnalysis();
//=============================================================================

  if (fPicoV0sClArr) {
    delete fPicoV0sClArr;
    fPicoV0sClArr = nullptr;
  }

  fPicoV0sClArr = new TClonesArray("AliPicoV0MC");
  fPicoV0sClArr->SetName("PicoV0s");
//=============================================================================

  if (fListUserOutputs) {
    delete fListUserOutputs;
    fListUserOutputs = nullptr;
  }

  fListUserOutputs = new TList();
  fListUserOutputs->SetOwner();

  CreateHistograms();
  PostData(1, fListUserOutputs);
//=============================================================================

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0MakerMC::UserExec(Option_t */*opt*/)
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::UserExec
//

  fPicoV0sClArr->Delete();
  if (!(InputEvent()->FindListObject("PicoV0s"))) InputEvent()->AddObject(fPicoV0sClArr);
//=============================================================================

  if (IsEventNotAcpt()) return;
//=============================================================================

  FillHistograms(); if (IsEventNotMBsa()) return;
//=============================================================================

  FillPicoV0s();
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0MakerMC::Terminate(Option_t */*opt*/)
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::Terminate
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0MakerMC::NotifyRun()
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::NotifyRun
//

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0MakerMC::FillPicoV0s()
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::FillPicoRecoV0s
//

  const auto nV0s(fEventAOD ? fEventAOD->GetNumberOfV0s() :
                              fEventESD->GetNumberOfV0s());

  if (nV0s<=0) return;
//=============================================================================

  auto nAt(fPicoV0sClArr->GetEntriesFast());
  auto hKshortPtInvM(static_cast<TH2D*>(fListUserOutputs->FindObject("hKshortPtInvM")));
  auto hLambdaPtInvM(static_cast<TH2D*>(fListUserOutputs->FindObject("hLambdaPtInvM")));
  auto hAntiLaPtInvM(static_cast<TH2D*>(fListUserOutputs->FindObject("hAntiLaPtInvM")));
//=============================================================================

  for (auto iV0=0; iV0<nV0s; iV0++) {
    AliPicoV0MC *pV0MC(nullptr);

    if (fEventAOD) {
      auto pV0(fEventAOD->GetV0(iV0));
      if (!pV0) continue;

      pV0MC = SelectV0Candidate(pV0);
    }

    if (fEventESD) {
      auto pV0(fEventESD->GetV0(iV0));
      if (!pV0) continue;

      pV0MC = SelectV0Candidate(pV0);
    }

    if (pV0MC) {
      pV0MC->FillKshortPtInvM(hKshortPtInvM);
      pV0MC->FillLambdaPtInvM(hLambdaPtInvM);
      pV0MC->FillAntiLaPtInvM(hAntiLaPtInvM);
      new ((*fPicoV0sClArr)[nAt++]) AliPicoV0MC(*pV0MC);
      delete pV0MC; pV0MC = nullptr;
    }
  }
//=============================================================================

  return;
}

//_____________________________________________________________________________
AliPicoV0MC *AliAnalysisTaskSEPicoV0MakerMC::SelectV0Candidate(AliAODv0 const *pV0RD)
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::SelectV0Candidate
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

  const auto dPosXrowsTPC = pDauPosRD->GetTPCClusterInfo(2,1);
  const auto dNegXrowsTPC = pDauNegRD->GetTPCClusterInfo(2,1);
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
  Double_t dPosEta = v3Pos.Eta(); if ((dPosEta<fCutMinDauEta) || (dPosEta>fCutMaxDauEta)) return nullptr;
  Double_t dNegEta = v3Neg.Eta(); if ((dNegEta<fCutMinDauEta) || (dNegEta>fCutMaxDauEta)) return nullptr;
//=============================================================================

  const auto inp(TMath::Abs(pDauPosRD->GetLabel())); if (inp<0) return nullptr;
  auto pDauPosMC(static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(inp))); if (!pDauPosMC) return nullptr;
  const auto imp(pDauPosMC->GetMother()); if (imp<0) return nullptr;

  const auto inn(TMath::Abs(pDauNegRD->GetLabel())); if (inn<0) return nullptr;
  auto pDauNegMC(static_cast<AliAODMCParticle*>(MCEvent()->GetTrack(inn))); if (!pDauNegMC) return nullptr;
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
  UInt_t wsmMC = 0;
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

    const auto dKshortInvM((vPosPion + vNegPion).M());
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
AliPicoV0MC *AliAnalysisTaskSEPicoV0MakerMC::SelectV0Candidate(AliESDv0 const *pV0RD)
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::SelectV0Candidate
//

  auto pStack(MCEvent()->Stack()); if (!pStack) return nullptr;
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

  auto  idmMC(0);
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

  Double_t dV0DistToPV(0.);
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
Bool_t AliAnalysisTaskSEPicoV0MakerMC::IsEventNotAcpt()
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::IsEventNotAcpt
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

  fEventAcptMask |= AliPicoBase::kEventCheck;
//=============================================================================

  if (!fMultEst.IsNull()) {
    Float_t dMult(-999.);

    if (fMultOld) {
      auto pCentSel(InputEvent()->GetCentrality());
      if (!pCentSel) { fEventAcptMask=0; return kTRUE; }
      if (pCentSel->GetQuality()!=0) return kFALSE;
      dMult = pCentSel->GetCentralityPercentile(fMultEst.Data());
    } else {
      auto pMultSel(static_cast<AliMultSelection*>(InputEvent()->FindListObject("MultSelection")));
      if (!pMultSel) { fEventAcptMask=0; return kTRUE; }
      dMult = pMultSel->GetMultiplicityPercentile(fMultEst.Data());
    }

    if ((dMult<fMultMin) || (dMult>=fMultMax)) return kFALSE;
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
//  if (pUtils->IsSPDClusterVsTrackletBG(InputEvent())) return kFALSE;
    if (pUtils->IsPileUpEvent(InputEvent())) return kFALSE;

    if ((fCollisionType==(AliPicoBase::kPA)) ||
        (fCollisionType==(AliPicoBase::kAP))) {
      if (pUtils->IsFirstEventInChunk(InputEvent())) return kFALSE;
    }
  }

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
Bool_t AliAnalysisTaskSEPicoV0MakerMC::IsEventNotINEL()
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::IsEventNotINEL
//

  if ((fEventAcptMask & AliPicoBase::kEventCheck) != AliPicoBase::kEventCheck) return kTRUE;
  if ((fEventAcptMask & AliPicoBase::kEventMult)  != AliPicoBase::kEventMult)  return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEPicoV0MakerMC::IsEventNotMBsa()
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::IsEventNotMBsa
//

  if (IsEventNotINEL()) return kTRUE;

  if ((fEventAcptMask & AliPicoBase::kEventTrigger) != AliPicoBase::kEventTrigger) return kTRUE;
  if ((fEventAcptMask & AliPicoBase::kEventVertex)  != AliPicoBase::kEventVertex)  return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0MakerMC::FillHistograms()
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::FillHistograms
//

  const auto h(static_cast<TH1D*>(fListUserOutputs->FindObject("hEventCheck")));

  if (fEventAcptMask==0) h->Fill(0.);
  if ((fEventAcptMask & AliPicoBase::kEventCheck)   == AliPicoBase::kEventCheck)   h->Fill(1.);
  if ((fEventAcptMask & AliPicoBase::kEventMult)    == AliPicoBase::kEventMult)    h->Fill(2.);
  if ((fEventAcptMask & AliPicoBase::kEventTrigger) == AliPicoBase::kEventTrigger) h->Fill(3.);
  if ((fEventAcptMask & AliPicoBase::kEventVertex)  == AliPicoBase::kEventVertex)  h->Fill(4.);

  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0MakerMC::CreateHistograms()
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::CreateHistograms
//

  const auto b(TH1::AddDirectoryStatus());
  TH1::AddDirectory(kFALSE);
//=============================================================================

  fListUserOutputs->Add(new TH1D("hEventCheck", "", 5, -0.5, 4.5));
  fListUserOutputs->Add(new TH2D("hKshortPtInvM", "", 1000, 0., 100., 300, AliPicoBase::MassKshort()-0.15,
                                                                           AliPicoBase::MassKshort()+0.15));

  fListUserOutputs->Add(new TH2D("hLambdaPtInvM", "", 1000, 0., 100., 200, AliPicoBase::MassLambda()-0.10,
                                                                           AliPicoBase::MassLambda()+0.10));

  fListUserOutputs->Add(new TH2D("hAntiLaPtInvM", "", 1000, 0., 100., 200, AliPicoBase::MassLambda()-0.10,
                                                                           AliPicoBase::MassLambda()+0.10));

  TH1 *h(nullptr);
  TListIter nl(fListUserOutputs);
  while ((h = static_cast<TH1*>(nl()))) h->Sumw2();

  TH1::AddDirectory(b);
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0MakerMC::InitAnalysis()
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::InitAnalysis
//
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPicoV0MakerMC::SetV0Cuts(Double_t d[14])
{
//
//  AliAnalysisTaskSEPicoV0MakerMC::SetCuts
//

  fCutMinV0Radius = d[0];
  fCutMaxV0Radius = d[1];

  fCutMaxDausDCA                     = d[2];
  fCutMinDauDCAtoPV                  = d[3];
  fCutMinDauXrowsTPC                 = d[4];
  fCutMinDauXrowsOverFindableClusTPC = d[5];

  fCutMaxKshortSigmaTPC = d[6];
  fCutMinKshortCosPA    = d[7];
  fCutMaxKshortCtau     = d[8];
  fCutMinKshortDeltaM   = d[9];

  fCutMaxLambdaSigmaTPC = d[10];
  fCutMinLambdaCosPA    = d[11];
  fCutMaxLambdaCtau     = d[12];
  fCutMinLambdaDeletaM  = d[13];

  return;
}
