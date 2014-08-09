#ifndef ALIPICOV0MC_H
#define ALIPICOV0MC_H

#include <TLorentzVector.h>

#include "AliPicoV0Base.h"


class AliPicoV0MC : public AliPicoV0Base {

 public :

  AliPicoV0MC();
  AliPicoV0MC(UInt_t   wMask,
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
              Int_t idV, UInt_t wsV, Double_t dV0Px,  Double_t dV0Py, Double_t dV0Pz, Double_t dV0E,
              Int_t idM=0, UInt_t wsM=0, Double_t dPtM=0., Double_t dEtaM=0., Double_t dRapM=0.);
  AliPicoV0MC(const AliPicoV0MC &src);
  AliPicoV0MC& operator=(const AliPicoV0MC &src);
  virtual ~AliPicoV0MC();
//=============================================================================

  TLorentzVector KineMC() const { return fV0Kine; }

  Double_t MotherPt()  const { return fMotherPt;  }
  Double_t MotherEta() const { return fMotherEta; }
  Double_t MotherRap() const { return fMotherRap; }
//=============================================================================

  Bool_t IsKshort(Double_t dCuts[9]);
  Bool_t IsLambda(Double_t dCuts[9]);
  Bool_t IsAntiLa(Double_t dCuts[9]);
  Bool_t IsV0InRapAcc(Double_t dMin, Double_t dMax);

  Bool_t IsKshort() const { return (fV0PDG== 310);  }
  Bool_t IsLambda() const { return (fV0PDG== 3122); }
  Bool_t IsAntiLa() const { return (fV0PDG==-3122); }

  Bool_t IsMotherXiNeg() const { return (fMotherPDG== 3312); }
  Bool_t IsMotherXiPos() const { return (fMotherPDG==-3312); }

  Bool_t IsLambdaFd() const { return (AliPicoV0Base::IsLambda() && IsMotherXiNeg()); }
  Bool_t IsAntiLaFd() const { return (AliPicoV0Base::IsAntiLa() && IsMotherXiPos()); }

  void GetControlVariables(Float_t d[18]);
//=============================================================================

  Bool_t IsV0Primary() const {
    return ((fV0Status & AliPicoHeaderCJ::kPrimary) == AliPicoHeaderCJ::kPrimary);
  }

  Bool_t IsV0PhysicalPrimary() const {
    return ((fV0Status & AliPicoHeaderCJ::kPhysicalPrimary) == AliPicoHeaderCJ::kPhysicalPrimary);
  }

  Bool_t IsV0SecondaryFromWeakDecay() const {
    return ((fV0Status & AliPicoHeaderCJ::kSecondaryFromWeakDecay) == AliPicoHeaderCJ::kSecondaryFromWeakDecay);
  }

  Bool_t IsV0SecondaryFromMaterial() const {
    return ((fV0Status & AliPicoHeaderCJ::kSecondaryFromMaterial) == AliPicoHeaderCJ::kSecondaryFromMaterial);
  }
//=============================================================================

  Bool_t IsMotherPrimary() const {
    return ((fMotherStatus & AliPicoHeaderCJ::kPrimary) == AliPicoHeaderCJ::kPrimary);
  }

  Bool_t IsMotherPhysicalPrimary() const {
    return ((fMotherStatus & AliPicoHeaderCJ::kPhysicalPrimary) == AliPicoHeaderCJ::kPhysicalPrimary);
  }

  Bool_t IsMotherSecondaryFromWeakDecay() const {
    return ((fMotherStatus & AliPicoHeaderCJ::kSecondaryFromWeakDecay) == AliPicoHeaderCJ::kSecondaryFromWeakDecay);
  }

  Bool_t IsMotherSecondaryFromMaterial() const {
    return ((fMotherStatus & AliPicoHeaderCJ::kSecondaryFromMaterial) == AliPicoHeaderCJ::kSecondaryFromMaterial);
  }
//=============================================================================

 private :

  Int_t  fV0PDG;     //
  UInt_t fV0Status;  //
  TLorentzVector fV0Kine;  //

  Int_t  fMotherPDG;     //
  UInt_t fMotherStatus;  //

  Double_t fMotherPt;   //
  Double_t fMotherEta;  //
  Double_t fMotherRap;  //

  ClassDef(AliPicoV0MC, 4);
};

#endif
