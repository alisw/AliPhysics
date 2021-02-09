#ifndef ALIPICOV0BASE_H
#define ALIPICOV0BASE_H

#include <TObject.h>
#include <TVector3.h>

#include "AliPicoBase.h"

class TH2D;
class TLorentzVector;

class AliPicoV0 : public TObject {

 public :

  AliPicoV0();
  AliPicoV0(UInt_t   wMask,
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
  AliPicoV0(const AliPicoV0 &src);
  AliPicoV0& operator=(const AliPicoV0 &src);
  virtual ~AliPicoV0();
//=============================================================================

  const TVector3 &KinePos() const { return  fP3Pos; }
  const TVector3 &KineNeg() const { return  fP3Neg; }
  const TVector3  KineRD()  const { return (fP3Pos + fP3Neg); }

  const TLorentzVector KineKshort() const;
  const TLorentzVector KineLambda() const;
  const TLorentzVector KineAntiLa() const;

  const Double_t RapidityKa() const;
  const Double_t RapidityLa() const;
//=============================================================================

  virtual Bool_t IsKshort(Double_t const */*dCuts*/=nullptr) const {
    return ((fMask & AliPicoBase::kKshort) == AliPicoBase::kKshort);
  }

  virtual Bool_t IsLambda(Double_t const */*dCuts*/=nullptr) const {
    return ((fMask & AliPicoBase::kLambda) == AliPicoBase::kLambda);
  }

  virtual Bool_t IsAntiLa(Double_t const */*dCuts*/=nullptr) const {
    return ((fMask & AliPicoBase::kAntiLambda) == AliPicoBase::kAntiLambda);
  }

  Bool_t IsKaInRapAcc(Double_t dMin, Double_t dMax);
  Bool_t IsLaInRapAcc(Double_t dMin, Double_t dMax);
  Bool_t IsV0InEtaAcc(Double_t dMin, Double_t dMax);
  Bool_t IsDausInEtaAcc(Double_t dMin, Double_t dMax);

  Bool_t IsPosInJC() const { return  fIsPosInJC; }
  Bool_t IsNegInJC() const { return  fIsNegInJC; }
  Bool_t IsTwoInJC() const { return (fIsPosInJC && fIsNegInJC); }
  Bool_t IsOneInJC() const { return (fIsPosInJC || fIsNegInJC); }

  virtual void GetControlVariables(Float_t */*d*/=nullptr) const = 0;
//=============================================================================

  void FillKshortPtInvM(TH2D* const h, Double_t const *dCuts=nullptr) const;
  void FillLambdaPtInvM(TH2D* const h, Double_t const *dCuts=nullptr) const;
  void FillAntiLaPtInvM(TH2D* const h, Double_t const *dCuts=nullptr) const;
//=============================================================================

 protected :

  Bool_t IsKa(const Double_t dCutMinV0Radius                    = 0.5,
              const Double_t dCutMinV0CosPA                     = 0.97,
              const Double_t dCutMaxV0Ctau                      = 20.,
              const Double_t dCutMaxDausDCA                     = 1.,
              const Double_t dCutMinPosDCAtoPV                  = 0.06,
              const Double_t dCutMinNegDCAtoPV                  = 0.06,
              const Float_t  dCutMinDauXrowsTPC                 = 70.,
              const Double_t dCutMinDauXrowsOverFindableClusTPC = 0.8,
              const Double_t dCutMinDauDeltaM                   = 0.005) const;

  Bool_t IsLa(const Double_t dCutMinV0Radius                    = 0.5,
              const Double_t dCutMinV0CosPA                     = 0.995,
              const Double_t dCutMaxV0Ctau                      = 30.,
              const Double_t dCutMaxDausDCA                     = 1.,
              const Double_t dCutMinPosDCAtoPV                  = 0.06,
              const Double_t dCutMinNegDCAtoPV                  = 0.06,
              const Float_t  dCutMinDauXrowsTPC                 = 70.,
              const Double_t dCutMinDauXrowsOverFindableClusTPC = 0.8,
              const Double_t dCutMinDauDeltaM                   = 0.01) const;

  Bool_t IsCandidateSelected(const Double_t dCutMinV0Radius,
                             const Double_t dCutMinV0CosPA,
                             const Double_t dCutMaxDausDCA,
                             const Double_t dCutMinPosDCAtoPV,
                             const Double_t dCutMinNegDCAtoPV,
                             const Float_t  dCutMinDauXrowsTPC,
                             const Double_t dCutMinDauXrowsOverFindableClusTPC) const;

  Bool_t IsKaSelected(const Double_t dCutMaxV0Ctau, const Double_t dCutMinDauDeltaM) const;
  Bool_t IsLaSelected(const Double_t dCutMaxV0Ctau, const Double_t dCutMinDauDeltaM) const;
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

  ClassDef(AliPicoV0, 5)
};

#endif
