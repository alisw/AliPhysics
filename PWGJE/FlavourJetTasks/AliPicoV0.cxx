#include <TH2D.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "AliPicoV0.h"

ClassImp(AliPicoV0)

//_____________________________________________________________________________
AliPicoV0::AliPicoV0() :
TObject(),
fMask(0),
fV0Radius(0.),
fV0CosPA(0.),
fV0DistToPVoverP(0.),
fDausDCA(0.),
fPosDCAtoPV(0.),
fNegDCAtoPV(0.),
fDauXrowsTPC(0.),
fDauXrowsOverFindableClusTPC(0.),
fP3Pos(),
fP3Neg(),
fIsPosInJC(kFALSE),
fIsNegInJC(kFALSE)
{
//
//  AliPicoV0::AliPicoV0
//
}

//_____________________________________________________________________________
AliPicoV0::AliPicoV0(UInt_t   wMask,
                             Double_t dV0Radius,
                             Double_t dV0CosPA,
                             Double_t dV0DistToPVoverP,
                             Double_t dDausDCA,
                             Double_t dPosDCAtoPV,
                             Double_t dNegDCAtoPV,
                             Float_t  dDauXrowsTPC,
                             Double_t dDauXrowsOverFindableClusTPC,
                             Double_t dPosPx, Double_t dPosPy, Double_t dPosPz,
                             Double_t dNegPx, Double_t dNegPy, Double_t dNegPz,
                             Bool_t bPosInJC,
                             Bool_t bNegInJC) :
TObject(),
fMask(wMask),
fV0Radius(dV0Radius),
fV0CosPA(dV0CosPA),
fV0DistToPVoverP(dV0DistToPVoverP),
fDausDCA(dDausDCA),
fPosDCAtoPV(dPosDCAtoPV),
fNegDCAtoPV(dNegDCAtoPV),
fDauXrowsTPC(dDauXrowsTPC),
fDauXrowsOverFindableClusTPC(dDauXrowsOverFindableClusTPC),
fP3Pos(dPosPx, dPosPy, dPosPz),
fP3Neg(dNegPx, dNegPy, dNegPz),
fIsPosInJC(bPosInJC),
fIsNegInJC(bNegInJC)
{
//
//  AliPicoV0::AliPicoV0
//
}

//_____________________________________________________________________________
AliPicoV0::AliPicoV0(const AliPicoV0 &src) :
TObject(src),
fMask(src.fMask),
fV0Radius(src.fV0Radius),
fV0CosPA(src.fV0CosPA),
fV0DistToPVoverP(src.fV0DistToPVoverP),
fDausDCA(src.fDausDCA),
fPosDCAtoPV(src.fPosDCAtoPV),
fNegDCAtoPV(src.fNegDCAtoPV),
fDauXrowsTPC(src.fDauXrowsTPC),
fDauXrowsOverFindableClusTPC(src.fDauXrowsOverFindableClusTPC),
fP3Pos(src.fP3Pos),
fP3Neg(src.fP3Neg),
fIsPosInJC(src.fIsPosInJC),
fIsNegInJC(src.fIsNegInJC)
{
//
//  AliPicoV0::AliPicoV0
//
}

//_____________________________________________________________________________
AliPicoV0& AliPicoV0::operator=(const AliPicoV0 &src)
{
//
//  AliPicoV0::operator=
//

  if (&src==this) return *this;

  TObject::operator=(src);

  fMask                        = src.fMask;
  fV0Radius                    = src.fV0Radius;
  fV0CosPA                     = src.fV0CosPA;
  fV0DistToPVoverP             = src.fV0DistToPVoverP;
  fDausDCA                     = src.fDausDCA;
  fPosDCAtoPV                  = src.fPosDCAtoPV;
  fNegDCAtoPV                  = src.fNegDCAtoPV;
  fDauXrowsTPC                 = src.fDauXrowsTPC;
  fDauXrowsOverFindableClusTPC = src.fDauXrowsOverFindableClusTPC;
  fP3Pos                       = src.fP3Pos;
  fP3Neg                       = src.fP3Neg;
  fIsPosInJC                   = src.fIsPosInJC;
  fIsNegInJC                   = src.fIsNegInJC;

  return *this;
}

//_____________________________________________________________________________
AliPicoV0::~AliPicoV0()
{
//
//  AliPicoV0::~AliPicoV0
//
}

//_____________________________________________________________________________
Bool_t AliPicoV0::IsKa(const Double_t dCutMinV0Radius,
                       const Double_t dCutMinV0CosPA,
                       const Double_t dCutMaxV0Ctau,
                       const Double_t dCutMaxDausDCA,
                       const Double_t dCutMinPosDCAtoPV,
                       const Double_t dCutMinNegDCAtoPV,
                       const Float_t  dCutMinDauXrowsTPC,
                       const Double_t dCutMinDauXrowsOverFindableClusTPC,
                       const Double_t dCutMinDauDeltaM) const
{
//
//  AliPicoV0::IsKa
//

  if (!IsCandidateSelected(dCutMinV0Radius,
                           dCutMinV0CosPA,
                           dCutMaxDausDCA,
                           dCutMinPosDCAtoPV,
                           dCutMinNegDCAtoPV,
                           dCutMinDauXrowsTPC,
                           dCutMinDauXrowsOverFindableClusTPC)) return kFALSE;

  if (!IsKaSelected(dCutMaxV0Ctau, dCutMinDauDeltaM)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0::IsLa(const Double_t dCutMinV0Radius,
                       const Double_t dCutMinV0CosPA,
                       const Double_t dCutMaxV0Ctau,
                       const Double_t dCutMaxDausDCA,
                       const Double_t dCutMinPosDCAtoPV,
                       const Double_t dCutMinNegDCAtoPV,
                       const Float_t  dCutMinDauXrowsTPC,
                       const Double_t dCutMinDauXrowsOverFindableClusTPC,
                       const Double_t dCutMinDauDeltaM) const
{
//
//  AliPicoV0::IsLa
//

  if (!IsCandidateSelected(dCutMinV0Radius,
                           dCutMinV0CosPA,
                           dCutMaxDausDCA,
                           dCutMinPosDCAtoPV,
                           dCutMinNegDCAtoPV,
                           dCutMinDauXrowsTPC,
                           dCutMinDauXrowsOverFindableClusTPC)) return kFALSE;

  if (!IsLaSelected(dCutMaxV0Ctau, dCutMinDauDeltaM)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0::IsKaSelected(const Double_t dCutMaxV0Ctau, const Double_t dCutMinDauDeltaM) const
{
//
//  AliPicoV0::IsKaSelected
//

  if ((fV0DistToPVoverP*AliPicoBase::MassKshort())>dCutMaxV0Ctau) return kFALSE;

  if (dCutMinDauDeltaM>0.) {
    Double_t dMassLambda = KineLambda().M();
    Double_t dMassAntiLa = KineAntiLa().M();
    if ((TMath::Abs(dMassLambda-AliPicoBase::MassLambda())<dCutMinDauDeltaM) ||
        (TMath::Abs(dMassAntiLa-AliPicoBase::MassLambda())<dCutMinDauDeltaM)) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0::IsLaSelected(const Double_t dCutMaxV0Ctau, const Double_t dCutMinDauDeltaM) const
{
//
//  AliPicoV0::IsLaSelected
//

  if ((fV0DistToPVoverP*AliPicoBase::MassLambda())>dCutMaxV0Ctau) return kFALSE;

  if (dCutMinDauDeltaM>0.) {
    Double_t dMassKshort = KineKshort().M();
    if (TMath::Abs(dMassKshort-AliPicoBase::MassKshort())<dCutMinDauDeltaM) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0::IsCandidateSelected(const Double_t dCutMinV0Radius,
                                      const Double_t dCutMinV0CosPA,
                                      const Double_t dCutMaxDausDCA,
                                      const Double_t dCutMinPosDCAtoPV,
                                      const Double_t dCutMinNegDCAtoPV,
                                      const Float_t  dCutMinDauXrowsTPC,
                                      const Double_t dCutMinDauXrowsOverFindableClusTPC) const
{
//
//  AliPicoV0::IsCandidateSelected
//

  if (fV0Radius<dCutMinV0Radius) return kFALSE;
  if (fV0CosPA <dCutMinV0CosPA)  return kFALSE;

  if (fDausDCA>dCutMaxDausDCA) return kFALSE;
  if (fPosDCAtoPV<dCutMinPosDCAtoPV) return kFALSE;
  if (fNegDCAtoPV<dCutMinNegDCAtoPV) return kFALSE;

  if (fDauXrowsTPC<dCutMinDauXrowsTPC) return kFALSE;
  if (fDauXrowsOverFindableClusTPC<dCutMinDauXrowsOverFindableClusTPC) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliPicoV0::FillKshortPtInvM(TH2D* const h, const Double_t *dCuts) const
{
//
//  AliPicoV0::FillKshortPtInvM
//

  if (!h) return;
  if (!IsKshort(dCuts)) return;
//=============================================================================

  TLorentzVector v = KineKshort();
  Double_t dPt = v.Pt(); if (dPt<h->GetXaxis()->GetBinLowEdge(1)) return;

  h->Fill(dPt, v.M());
  return;
}

//_____________________________________________________________________________
void AliPicoV0::FillLambdaPtInvM(TH2D* const h, const Double_t *dCuts) const
{
//
//  AliPicoV0::FillLambdaPtInvM
//

  if (!h) return;
  if (!IsLambda(dCuts)) return;
//=============================================================================

  TLorentzVector v = KineLambda();
  Double_t dPt = v.Pt(); if (dPt<h->GetXaxis()->GetBinLowEdge(1)) return;

  h->Fill(dPt, v.M());
  return;
}

//_____________________________________________________________________________
void AliPicoV0::FillAntiLaPtInvM(TH2D* const h, const Double_t *dCuts) const
{
//
//  AliPicoV0::FillAntiLaPtInvM
//

  if (!h) return;
  if (!IsAntiLa(dCuts)) return;
//=============================================================================

  TLorentzVector v = KineAntiLa();
  Double_t dPt = v.Pt(); if (dPt<h->GetXaxis()->GetBinLowEdge(1)) return;

  h->Fill(dPt, v.M());
  return;
}

//_____________________________________________________________________________
Bool_t AliPicoV0::IsKaInRapAcc(Double_t dMin, Double_t dMax)
{
//
//  AliPicoV0::IsKaInRapAcc
//

  Double_t dRap = RapidityKa(); return ((dRap>=dMin) && (dRap<dMax));
}

//_____________________________________________________________________________
Bool_t AliPicoV0::IsLaInRapAcc(Double_t dMin, Double_t dMax)
{
//
//  AliPicoV0::IsLaInRapAcc
//

  Double_t dRap = RapidityLa(); return ((dRap>=dMin) && (dRap<dMax));
}

//_____________________________________________________________________________
Bool_t AliPicoV0::IsV0InEtaAcc(Double_t dMin, Double_t dMax)
{
//
//  AliPicoV0::IsV0InEtaAcc
//

  Double_t dEta = KineRD().Eta();
  if ((dEta<dMin) || (dEta>=dMax)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0::IsDausInEtaAcc(Double_t dMin, Double_t dMax)
{
//
//  AliPicoV0::IsDausInEtaAcc
//

  Double_t dPosEta = fP3Pos.Eta(); if ((dPosEta<dMin) || (dPosEta>=dMax)) return kFALSE;
  Double_t dNegEta = fP3Neg.Eta(); if ((dNegEta<dMin) || (dNegEta>=dMax)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
const Double_t AliPicoV0::RapidityKa() const
{
//
//  AliPicoV0::RapidityKa
//

  TLorentzVector v; v.SetVectM(KineRD(), AliPicoBase::MassKshort());
  return v.Rapidity();
}

//_____________________________________________________________________________
const Double_t AliPicoV0::RapidityLa() const
{
//
//  AliPicoV0::RapidityLa
//

  TLorentzVector v; v.SetVectM(KineRD(), AliPicoBase::MassLambda());
  return v.Rapidity();
}

//_____________________________________________________________________________
const TLorentzVector AliPicoV0::KineKshort() const
{
//
//  AliPicoV0::KineKshort
//

  TLorentzVector vPos; vPos.SetVectM(fP3Pos, AliPicoBase::MassPion());
  TLorentzVector vNeg; vNeg.SetVectM(fP3Neg, AliPicoBase::MassPion());

  return (vPos + vNeg);
}

//_____________________________________________________________________________
const TLorentzVector AliPicoV0::KineLambda() const
{
//
//  AliPicoV0::KineLambda
//

  TLorentzVector vPos; vPos.SetVectM(fP3Pos, AliPicoBase::MassProton());
  TLorentzVector vNeg; vNeg.SetVectM(fP3Neg, AliPicoBase::MassPion());

  return (vPos + vNeg);
}

//_____________________________________________________________________________
const TLorentzVector AliPicoV0::KineAntiLa() const
{
//
//  AliPicoV0::KineAntiLa
//

  TLorentzVector vPos; vPos.SetVectM(fP3Pos, AliPicoBase::MassPion());
  TLorentzVector vNeg; vNeg.SetVectM(fP3Neg, AliPicoBase::MassProton());

  return (vPos + vNeg);
}
