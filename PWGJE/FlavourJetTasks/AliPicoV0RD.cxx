#include <TMath.h>
#include <TLorentzVector.h>

#include "AliPicoV0RD.h"

ClassImp(AliPicoV0RD)

//_____________________________________________________________________________
AliPicoV0RD::AliPicoV0RD() :
AliPicoV0(),
fPosPionSigmaTPC(0.),
fNegPionSigmaTPC(0.),
fPosProtonSigmaTPC(0.),
fNegProtonSigmaTPC(0.)
{
//
//  AliPicoV0RD::AliPicoV0RD
//
}

//_____________________________________________________________________________
AliPicoV0RD::AliPicoV0RD(UInt_t   wMask,
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
                         Bool_t   bPosInJC, Bool_t bNegInJC,
                         Float_t  dPosPionSigmaTPC, Float_t dPosProtonSigmaTPC,
                         Float_t  dNegPionSigmaTPC, Float_t dNegProtonSigmaTPC) :
AliPicoV0(wMask,
          dV0Radius,
          dV0CosPA,
          dV0DistToPVoverP,
          dDausDCA,
          dPosDCAtoPV,
          dNegDCAtoPV,
          dDauXrowsTPC,
          dDauXrowsOverFindableClusTPC,
          dPosPx, dPosPy, dPosPz,
          dNegPx, dNegPy, dNegPz,
          bPosInJC, bNegInJC),
fPosPionSigmaTPC(dPosPionSigmaTPC),
fNegPionSigmaTPC(dNegPionSigmaTPC),
fPosProtonSigmaTPC(dPosProtonSigmaTPC),
fNegProtonSigmaTPC(dNegProtonSigmaTPC)
{
//
//  AliPicoV0RD::AliPicoV0RD
//
}

//_____________________________________________________________________________
AliPicoV0RD::AliPicoV0RD(const AliPicoV0RD &src) :
AliPicoV0(src),
fPosPionSigmaTPC(src.fPosPionSigmaTPC),
fNegPionSigmaTPC(src.fNegPionSigmaTPC),
fPosProtonSigmaTPC(src.fPosProtonSigmaTPC),
fNegProtonSigmaTPC(src.fNegProtonSigmaTPC)
{
//
//  AliPicoV0RD::AliPicoV0RD
//
}

//_____________________________________________________________________________
AliPicoV0RD& AliPicoV0RD::operator=(const AliPicoV0RD &src)
{
//
//  AliPicoV0RD::operator=
//

  if (&src==this) return *this;

  AliPicoV0::operator=(src);

  fPosPionSigmaTPC   = src.fPosPionSigmaTPC;
  fNegPionSigmaTPC   = src.fNegPionSigmaTPC;
  fPosProtonSigmaTPC = src.fPosProtonSigmaTPC;
  fNegProtonSigmaTPC = src.fNegProtonSigmaTPC;

  return *this;
}

//_____________________________________________________________________________
AliPicoV0RD::~AliPicoV0RD()
{
//
//  AliPicoV0RD::~AliPicoV0RD
//
}

//_____________________________________________________________________________
Bool_t AliPicoV0RD::IsKshort(Double_t const dCuts[10]) const
{
//
//  Bool_t AliPicoV0RD::IsKshort(Double_t dCuts[10]) const
//

  if (!AliPicoV0::IsKshort()) return kFALSE;
  if (!dCuts) return kTRUE;
//=============================================================================

  if (dCuts[9]>0.) {
    if (!((TMath::Abs(fPosPionSigmaTPC)<dCuts[9]) &&
          (TMath::Abs(fNegPionSigmaTPC)<dCuts[9]))) return kFALSE;
  }

  if (!IsKa(dCuts[0],dCuts[1],dCuts[2],dCuts[3],dCuts[4],dCuts[5],dCuts[6],dCuts[7],dCuts[8])) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0RD::IsLambda(Double_t const dCuts[10]) const
{
//
//  Bool_t AliPicoV0RD::IsLambda(Double_t dCuts[10]) const
//

  if (!AliPicoV0::IsLambda()) return kFALSE;
  if (!dCuts) return kTRUE;
//=============================================================================

  if (dCuts[9]>0.) {
    if (!((TMath::Abs(fPosProtonSigmaTPC)<dCuts[9]) &&
          (TMath::Abs(fNegPionSigmaTPC)  <dCuts[9]))) return kFALSE;
  }

  if (!IsLa(dCuts[0],dCuts[1],dCuts[2],dCuts[3],dCuts[4],dCuts[5],dCuts[6],dCuts[7],dCuts[8])) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0RD::IsAntiLa(Double_t const dCuts[10]) const
{
//
//  Bool_t AliPicoV0RD::IsAntiLa(Double_t dCuts[10]) const
//

  if (!AliPicoV0::IsAntiLa()) return kFALSE;
  if (!dCuts) return kTRUE;
//=============================================================================

  if (dCuts[9]>0.) {
    if (!((TMath::Abs(fPosPionSigmaTPC)  <dCuts[9]) &&
          (TMath::Abs(fNegProtonSigmaTPC)<dCuts[9]))) return kFALSE;
  }

  if (!IsLa(dCuts[0],dCuts[1],dCuts[2],dCuts[3],dCuts[4],dCuts[5],dCuts[6],dCuts[7],dCuts[8])) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
void AliPicoV0RD::GetControlVariables(Float_t d[20]) const
{
//
//  void AliPicoV0RD::GetControlVariables(Float_t d[20]) const
//

  d[ 0] = (Float_t)fV0Radius;
  d[ 1] = (Float_t)fV0CosPA;
  d[ 2] = (Float_t)fV0DistToPVoverP;
  d[ 3] = (Float_t)fDausDCA;
  d[ 4] = (Float_t)fPosDCAtoPV;
  d[ 5] = (Float_t)fNegDCAtoPV;
  d[ 6] = (Float_t)fDauXrowsTPC;
  d[ 7] = (Float_t)fDauXrowsOverFindableClusTPC;
  d[ 8] = (Float_t)fPosPionSigmaTPC;
  d[ 9] = (Float_t)fNegPionSigmaTPC;
  d[10] = (Float_t)fPosProtonSigmaTPC;
  d[11] = (Float_t)fNegProtonSigmaTPC;
  d[12] = (Float_t)KineRD().Pt();
  d[13] = (Float_t)RapidityKa();
  d[14] = (Float_t)RapidityLa();
  d[15] = (Float_t)KineKshort().M();
  d[16] = (Float_t)KineLambda().M();
  d[17] = (Float_t)KineAntiLa().M();
  d[18] = (Float_t)fP3Pos.Eta();
  d[19] = (Float_t)fP3Neg.Eta();

  return;
}
