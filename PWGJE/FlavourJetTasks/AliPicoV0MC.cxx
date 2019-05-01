#include "AliPicoV0MC.h"

ClassImp(AliPicoV0MC)

//_____________________________________________________________________________
AliPicoV0MC::AliPicoV0MC() :
AliPicoV0(),
fV0PDG(0),
fV0Status(0),
fV0Kine(),
fMotherPDG(0),
fMotherStatus(0),
fMotherPt(0.),
fMotherEta(0.),
fMotherRap(0.)
{
//
//  AliPicoV0MC::AliPicoV0MC
//
}

//_____________________________________________________________________________
AliPicoV0MC::AliPicoV0MC(UInt_t   wMask,
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
                         Int_t idV, UInt_t wsV, Double_t dV0Px, Double_t dV0Py, Double_t dV0Pz, Double_t dV0E,
                         Int_t idM, UInt_t wsM, Double_t dPtM,  Double_t dEtaM, Double_t dRapM) :
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
fV0PDG(idV),
fV0Status(wsV),
fV0Kine(dV0Px, dV0Py, dV0Pz, dV0E),
fMotherPDG(idM),
fMotherStatus(wsM),
fMotherPt(dPtM),
fMotherEta(dEtaM),
fMotherRap(dRapM)
{
//
//  AliPicoV0MC::AliPicoV0MC
//
}

//_____________________________________________________________________________
AliPicoV0MC::AliPicoV0MC(const AliPicoV0MC &src) :
AliPicoV0(src),
fV0PDG(src.fV0PDG),
fV0Status(src.fV0Status),
fV0Kine(src.fV0Kine),
fMotherPDG(src.fMotherPDG),
fMotherStatus(src.fMotherStatus),
fMotherPt(src.fMotherPt),
fMotherEta(src.fMotherEta),
fMotherRap(src.fMotherRap)
{
//
//  AliPicoV0MC::AliPicoV0MC
//
}

//_____________________________________________________________________________
AliPicoV0MC& AliPicoV0MC::operator=(const AliPicoV0MC &src)
{
//
//  AliPicoV0RD::operator=
//

  if (&src==this) return *this;

  AliPicoV0::operator=(src);

  fV0PDG        = src.fV0PDG;
  fV0Status     = src.fV0Status;
  fV0Kine       = src.fV0Kine;

  fMotherPDG    = src.fMotherPDG;
  fMotherStatus = src.fMotherStatus;
  fMotherPt     = src.fMotherPt;
  fMotherEta    = src.fMotherEta;
  fMotherRap    = src.fMotherRap;

  return *this;
}

//_____________________________________________________________________________
AliPicoV0MC::~AliPicoV0MC()
{
//
//  AliPicoV0MC::~AliPicoV0MC
//
}

//_____________________________________________________________________________
Bool_t AliPicoV0MC::IsKshort(Double_t const dCuts[9]) const
{
//
//  Bool_t AliPicoV0MC::IsKshort(Double_t dCuts[9]) const
//

  if (!IsKshortMC()) return kFALSE;
  if (!AliPicoV0::IsKshort()) return kFALSE;
  if (!dCuts) return kTRUE;
//=============================================================================

  if (!IsKa(dCuts[0],dCuts[1],dCuts[2],dCuts[3],dCuts[4],dCuts[5],dCuts[6],dCuts[7],dCuts[8])) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0MC::IsLambda(Double_t const dCuts[9]) const
{
//
//  Bool_t AliPicoV0MC::IsLambda(Double_t dCuts[9]) const
//

  if (!IsLambdaMC()) return kFALSE;
  if (!AliPicoV0::IsLambda()) return kFALSE;
  if (!dCuts) return kTRUE;
//=============================================================================

  if (!IsLa(dCuts[0],dCuts[1],dCuts[2],dCuts[3],dCuts[4],dCuts[5],dCuts[6],dCuts[7],dCuts[8])) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0MC::IsAntiLa(Double_t const dCuts[9]) const
{
//
//  Bool_t AliPicoV0MC::IsAntiLa(Double_t dCuts[9]) const
//

  if (!IsAntiLaMC()) return kFALSE;
  if (!AliPicoV0::IsAntiLa()) return kFALSE;
  if (!dCuts) return kTRUE;
//=============================================================================

  if (!IsLa(dCuts[0],dCuts[1],dCuts[2],dCuts[3],dCuts[4],dCuts[5],dCuts[6],dCuts[7],dCuts[8])) return kFALSE;

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliPicoV0MC::IsV0InRapAcc(Double_t dMin, Double_t dMax)
{
//
//  AliPicoV0MC::IsV0InRapAcc
//

  Double_t dRap = fV0Kine.Rapidity();
  return ((dRap>=dMin) && (dRap<dMax));
}

//_____________________________________________________________________________
void AliPicoV0MC::GetControlVariables(Float_t d[18]) const
{
//
//  AliPicoV0MC::GetControlVariables
//

  d[ 0] = (Float_t)fV0Radius;
  d[ 1] = (Float_t)fV0CosPA;
  d[ 2] = (Float_t)fV0DistToPVoverP;
  d[ 3] = (Float_t)fDausDCA;
  d[ 4] = (Float_t)fPosDCAtoPV;
  d[ 5] = (Float_t)fNegDCAtoPV;
  d[ 6] = (Float_t)fDauXrowsTPC;
  d[ 7] = (Float_t)fDauXrowsOverFindableClusTPC;
  d[ 8] = (Float_t)KineRD().Pt();
  d[ 9] = (Float_t)RapidityKa();
  d[10] = (Float_t)RapidityLa();
  d[11] = (Float_t)KineMC().Pt();
  d[12] = (Float_t)KineMC().Rapidity();
  d[13] = (Float_t)KineKshort().M();
  d[14] = (Float_t)KineLambda().M();
  d[15] = (Float_t)KineAntiLa().M();
  d[16] = (Float_t)fP3Pos.Eta();
  d[17] = (Float_t)fP3Neg.Eta();

  return;
}
