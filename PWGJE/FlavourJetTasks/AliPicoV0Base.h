#ifndef ALIPICOV0BASE_H
#define ALIPICOV0BASE_H

#include <TObject.h>
#include <TVector3.h>

#include "AliPicoHeaderCJ.h"

class TH2D;
class TLorentzVector;

class AliPicoV0Base : public TObject {

 public :

  AliPicoV0Base();
  AliPicoV0Base(UInt_t   wMask,
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
                Bool_t   bPosInJC, Bool_t   bNegInJC);
  AliPicoV0Base(const AliPicoV0Base &src);
  AliPicoV0Base& operator=(const AliPicoV0Base &src);
  virtual ~AliPicoV0Base();
//=============================================================================

  TVector3 KinePos() const { return  fP3Pos; }
  TVector3 KineNeg() const { return  fP3Neg; }
  TVector3 KineRD()  const { return (fP3Pos + fP3Neg); }

  TLorentzVector KineKshort();
  TLorentzVector KineLambda();
  TLorentzVector KineAntiLa();

  Double_t RapidityKa();
  Double_t RapidityLa();
//=============================================================================

  Bool_t IsKshort() const { return ((fMask & AliPicoHeaderCJ::kKshort)     == AliPicoHeaderCJ::kKshort);     }
  Bool_t IsLambda() const { return ((fMask & AliPicoHeaderCJ::kLambda)     == AliPicoHeaderCJ::kLambda);     }
  Bool_t IsAntiLa() const { return ((fMask & AliPicoHeaderCJ::kAntiLambda) == AliPicoHeaderCJ::kAntiLambda); }

  Bool_t IsKaInRapAcc(Double_t dMin, Double_t dMax);
  Bool_t IsLaInRapAcc(Double_t dMin, Double_t dMax);
  Bool_t IsV0InEtaAcc(Double_t dMin, Double_t dMax);
  Bool_t IsDausInEtaAcc(Double_t dMin, Double_t dMax);

  Bool_t IsPosInJC() const { return  fIsPosInJC; }
  Bool_t IsNegInJC() const { return  fIsNegInJC; }
  Bool_t IsTwoInJC() const { return (fIsPosInJC && fIsNegInJC); }
  Bool_t IsOneInJC() const { return (fIsPosInJC || fIsNegInJC); }

  void FillKshortPtInvM(TH2D *h);
  void FillLambdaPtInvM(TH2D *h);
  void FillAntiLaPtInvM(TH2D *h);
//=============================================================================

 protected :

  Bool_t IsKa(Double_t dCutMinV0Radius                    = 0.5,
              Double_t dCutMinV0CosPA                     = 0.97,
              Double_t dCutMaxV0Ctau                      = 20.,
              Double_t dCutMaxDausDCA                     = 1.,
              Double_t dCutMinPosDCAtoPV                  = 0.06,
              Double_t dCutMinNegDCAtoPV                  = 0.06,
              Float_t  dCutMinDauXrowsTPC                 = 70.,
              Double_t dCutMinDauXrowsOverFindableClusTPC = 0.8,
              Double_t dCutMinDauDeltaM                   = 0.005);

  Bool_t IsLa(Double_t dCutMinV0Radius                    = 0.5,
              Double_t dCutMinV0CosPA                     = 0.995,
              Double_t dCutMaxV0Ctau                      = 30.,
              Double_t dCutMaxDausDCA                     = 1.,
              Double_t dCutMinPosDCAtoPV                  = 0.06,
              Double_t dCutMinNegDCAtoPV                  = 0.06,
              Float_t  dCutMinDauXrowsTPC                 = 70.,
              Double_t dCutMinDauXrowsOverFindableClusTPC = 0.8,
              Double_t dCutMinDauDeltaM                   = 0.01);

  Bool_t IsCandidateSelected(Double_t dCutMinV0Radius,
                             Double_t dCutMinV0CosPA,
                             Double_t dCutMaxDausDCA,
                             Double_t dCutMinPosDCAtoPV,
                             Double_t dCutMinNegDCAtoPV,
                             Float_t  dCutMinDauXrowsTPC,
                             Double_t dCutMinDauXrowsOverFindableClusTPC);

  Bool_t IsKaSelected(Double_t dCutMaxV0Ctau, Double_t dCutMinDauDeltaM);
  Bool_t IsLaSelected(Double_t dCutMaxV0Ctau, Double_t dCutMinDauDeltaM);
//=============================================================================

  UInt_t fMask;  //

  Double_t fV0Radius;  //
  Double_t fV0CosPA;  //
  Double_t fV0DistToPVoverP;  //

  Double_t fDausDCA;  //

  Double_t fPosDCAtoPV;  //
  Double_t fNegDCAtoPV;  //

  Float_t  fDauXrowsTPC;  //
  Double_t fDauXrowsOverFindableClusTPC;  //

  TVector3 fP3Pos;  //
  TVector3 fP3Neg;  //

  Bool_t fIsPosInJC;  // match w/ jet consti
  Bool_t fIsNegInJC;  // match w/ jet consti

  static const Double_t fgkMassPion;    //
  static const Double_t fgkMassKshort;  //
  static const Double_t fgkMassLambda;  //
  static const Double_t fgkMassProton;  //

  ClassDef(AliPicoV0Base, 3)
};

#endif
