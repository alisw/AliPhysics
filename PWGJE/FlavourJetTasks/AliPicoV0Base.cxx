#include <TH2D.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "AliPicoV0Base.h"

ClassImp(AliPicoV0Base)

//_____________________________________________________________________________
const Double_t AliPicoV0Base::fgkMassPion   = 0.13957;
const Double_t AliPicoV0Base::fgkMassKshort = 0.497614;
const Double_t AliPicoV0Base::fgkMassProton = 0.938272;
const Double_t AliPicoV0Base::fgkMassLambda = 1.11568;

//_____________________________________________________________________________
AliPicoV0Base::AliPicoV0Base() :
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
//  AliPicoV0Base::AliPicoV0Base
//
}

//_____________________________________________________________________________
AliPicoV0Base::AliPicoV0Base(UInt_t   wMask,
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
//  AliPicoV0Base::AliPicoV0Base
//
}

//_____________________________________________________________________________
AliPicoV0Base::AliPicoV0Base(const AliPicoV0Base &src) :
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
//  AliPicoV0Base::AliPicoV0Base
//
}

//_____________________________________________________________________________
AliPicoV0Base& AliPicoV0Base::operator=(const AliPicoV0Base &src)
{
//
//  AliPicoV0Base::operator=
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
AliPicoV0Base::~AliPicoV0Base()
{
//
//  AliPicoV0Base::~AliPicoV0Base
//
}

//_____________________________________________________________________________
Bool_t AliPicoV0Base::IsKa(Double_t dCutMinV0Radius,
                           Double_t dCutMinV0CosPA,
                           Double_t dCutMaxV0Ctau,
                           Double_t dCutMaxDausDCA,
                           Double_t dCutMinPosDCAtoPV,
                           Double_t dCutMinNegDCAtoPV,
                           Float_t  dCutMinDauXrowsTPC,
                           Double_t dCutMinDauXrowsOverFindableClusTPC,
                           Double_t dCutMinDauDeltaM)
{
//
//  AliPicoV0Base::IsKa
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
Bool_t AliPicoV0Base::IsLa(Double_t dCutMinV0Radius,
                           Double_t dCutMinV0CosPA,
                           Double_t dCutMaxV0Ctau,
                           Double_t dCutMaxDausDCA,
                           Double_t dCutMinPosDCAtoPV,
                           Double_t dCutMinNegDCAtoPV,
                           Float_t  dCutMinDauXrowsTPC,
                           Double_t dCutMinDauXrowsOverFindableClusTPC,
                           Double_t dCutMinDauDeltaM)
{
//
//  AliPicoV0Base::IsLa
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
Bool_t AliPicoV0Base::IsKaSelected(Double_t dCutMaxV0Ctau, Double_t dCutMinDauDeltaM)
{
//
//  AliPicoV0Base::IsKaSelected
//

  if ((fV0DistToPVoverP*fgkMassKshort)>dCutMaxV0Ctau) return kFALSE;

  if (dCutMinDauDeltaM>0.) {
    Double_t dMassLambda = KineLambda().M();
    Double_t dMassAntiLa = KineAntiLa().M();
    if ((TMath::Abs(dMassLambda-fgkMassLambda)<dCutMinDauDeltaM) ||
        (TMath::Abs(dMassAntiLa-fgkMassLambda)<dCutMinDauDeltaM)) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0Base::IsLaSelected(Double_t dCutMaxV0Ctau, Double_t dCutMinDauDeltaM)
{
//
//  AliPicoV0Base::IsLaSelected
//

  if ((fV0DistToPVoverP*fgkMassLambda)>dCutMaxV0Ctau) return kFALSE;

  if (dCutMinDauDeltaM>0.) {
    Double_t dMassKshort = KineKshort().M();
    if (TMath::Abs(dMassKshort-fgkMassKshort)<dCutMinDauDeltaM) return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0Base::IsCandidateSelected(Double_t dCutMinV0Radius,
                                          Double_t dCutMinV0CosPA,
                                          Double_t dCutMaxDausDCA,
                                          Double_t dCutMinPosDCAtoPV,
                                          Double_t dCutMinNegDCAtoPV,
                                          Float_t  dCutMinDauXrowsTPC,
                                          Double_t dCutMinDauXrowsOverFindableClusTPC)
{
//
//  AliPicoV0Base::IsCandidateSelected
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
void AliPicoV0Base::FillKshortPtInvM(TH2D *h)
{
//
//  AliPicoV0Base::FillKshortPtInvM
//

  if (!h) return;
  if (!IsKshort()) return;
//=============================================================================

  TLorentzVector v = KineKshort();
  Double_t dPt = v.Pt(); if (dPt<h->GetXaxis()->GetBinLowEdge(1)) return;

  h->Fill(dPt, v.M());
  return;
}

//_____________________________________________________________________________
void AliPicoV0Base::FillLambdaPtInvM(TH2D *h)
{
//
//  AliPicoV0Base::FillLambdaPtInvM
//

  if (!h) return;
  if (!IsLambda()) return;
//=============================================================================

  TLorentzVector v = KineLambda();
  Double_t dPt = v.Pt(); if (dPt<h->GetXaxis()->GetBinLowEdge(1)) return;

  h->Fill(dPt, v.M());
  return;
}

//_____________________________________________________________________________
void AliPicoV0Base::FillAntiLaPtInvM(TH2D *h)
{
//
//  AliPicoV0Base::FillAntiLaPtInvM
//

  if (!h) return;
  if (!IsAntiLa()) return;
//=============================================================================

  TLorentzVector v = KineAntiLa();
  Double_t dPt = v.Pt(); if (dPt<h->GetXaxis()->GetBinLowEdge(1)) return;

  h->Fill(dPt, v.M());
  return;
}

//_____________________________________________________________________________
Bool_t AliPicoV0Base::IsKaInRapAcc(Double_t dMin, Double_t dMax)
{
//
//  AliPicoV0Base::IsKaInRapAcc
//

  if (!IsKshort()) return kFALSE;
  Double_t dRap = RapidityKa(); return ((dRap>=dMin) && (dRap<dMax));
}

//_____________________________________________________________________________
Bool_t AliPicoV0Base::IsLaInRapAcc(Double_t dMin, Double_t dMax)
{
//
//  AliPicoV0Base::IsLaInRapAcc
//

  if (!(IsLambda() || IsAntiLa())) return kFALSE;
  Double_t dRap = RapidityLa(); return ((dRap>=dMin) && (dRap<dMax));
}

//_____________________________________________________________________________
Bool_t AliPicoV0Base::IsV0InEtaAcc(Double_t dMin, Double_t dMax)
{
//
//  AliPicoV0Base::IsV0InEtaAcc
//

  Double_t dEta = KineRD().Eta();
  if ((dEta<dMin) || (dEta>=dMax)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0Base::IsDausInEtaAcc(Double_t dMin, Double_t dMax)
{
//
//  AliPicoV0Base::IsDausInEtaAcc
//

  Double_t dPosEta = fP3Pos.Eta(); if ((dPosEta<dMin) || (dPosEta>=dMax)) return kFALSE;
  Double_t dNegEta = fP3Neg.Eta(); if ((dNegEta<dMin) || (dNegEta>=dMax)) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Double_t AliPicoV0Base::RapidityKa()
{
//
//  AliPicoV0Base::RapidityKa
//

  TLorentzVector v; v.SetVectM(KineRD(), fgkMassKshort);
  return v.Rapidity();
}

//_____________________________________________________________________________
Double_t AliPicoV0Base::RapidityLa()
{
//
//  AliPicoV0Base::RapidityLa
//

  TLorentzVector v; v.SetVectM(KineRD(), fgkMassLambda);
  return v.Rapidity();
}

//_____________________________________________________________________________
TLorentzVector AliPicoV0Base::KineKshort()
{
//
//  AliPicoV0Base::KineKshort
//

  TLorentzVector vPos; vPos.SetVectM(fP3Pos, fgkMassPion);
  TLorentzVector vNeg; vNeg.SetVectM(fP3Neg, fgkMassPion);

  return (vPos + vNeg);
}

//_____________________________________________________________________________
TLorentzVector AliPicoV0Base::KineLambda()
{
//
//  AliPicoV0Base::KineLambda
//

  TLorentzVector vPos; vPos.SetVectM(fP3Pos, fgkMassProton);
  TLorentzVector vNeg; vNeg.SetVectM(fP3Neg, fgkMassPion);

  return (vPos + vNeg);
}

//_____________________________________________________________________________
TLorentzVector AliPicoV0Base::KineAntiLa()
{
//
//  AliPicoV0Base::KineAntiLa
//

  TLorentzVector vPos; vPos.SetVectM(fP3Pos, fgkMassPion);
  TLorentzVector vNeg; vNeg.SetVectM(fP3Neg, fgkMassProton);

  return (vPos + vNeg);
}
