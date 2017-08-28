#ifndef ALIPICOV0MC_H
#define ALIPICOV0MC_H

#include <TLorentzVector.h>

#include "AliPicoV0.h"


class AliPicoV0MC : public AliPicoV0 {

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

  virtual Bool_t IsKshort(Double_t const dCuts[9]) const;
  virtual Bool_t IsLambda(Double_t const dCuts[9]) const;
  virtual Bool_t IsAntiLa(Double_t const dCuts[9]) const;

  virtual void GetControlVariables(Float_t d[18]) const;
//=============================================================================

  TLorentzVector KineMC() const { return fV0Kine; }

  Double_t MotherPt()  const { return fMotherPt;  }
  Double_t MotherEta() const { return fMotherEta; }
  Double_t MotherRap() const { return fMotherRap; }
//=============================================================================

  Bool_t IsV0InRapAcc(Double_t dMin, Double_t dMax);

  Bool_t IsKshortMC() const { return (fV0PDG== 310);  }
  Bool_t IsLambdaMC() const { return (fV0PDG== 3122); }
  Bool_t IsAntiLaMC() const { return (fV0PDG==-3122); }

  Bool_t IsMotherXiNeg() const { return (fMotherPDG== 3312); }
  Bool_t IsMotherXiPos() const { return (fMotherPDG==-3312); }

  Bool_t IsLambdaFd() const { return (AliPicoV0::IsLambda() && IsMotherXiNeg()); }
  Bool_t IsAntiLaFd() const { return (AliPicoV0::IsAntiLa() && IsMotherXiPos()); }
//=============================================================================

  Bool_t IsV0Primary() const {
    return ((fV0Status & AliPicoBase::kPrimary) == AliPicoBase::kPrimary);
  }

  Bool_t IsV0PhysicalPrimary() const {
    return ((fV0Status & AliPicoBase::kPhysicalPrimary) == AliPicoBase::kPhysicalPrimary);
  }

  Bool_t IsV0SecondaryFromWeakDecay() const {
    return ((fV0Status & AliPicoBase::kSecondaryFromWeakDecay) == AliPicoBase::kSecondaryFromWeakDecay);
  }

  Bool_t IsV0SecondaryFromMaterial() const {
    return ((fV0Status & AliPicoBase::kSecondaryFromMaterial) == AliPicoBase::kSecondaryFromMaterial);
  }
//=============================================================================

  Bool_t IsMotherPrimary() const {
    return ((fMotherStatus & AliPicoBase::kPrimary) == AliPicoBase::kPrimary);
  }

  Bool_t IsMotherPhysicalPrimary() const {
    return ((fMotherStatus & AliPicoBase::kPhysicalPrimary) == AliPicoBase::kPhysicalPrimary);
  }

  Bool_t IsMotherSecondaryFromWeakDecay() const {
    return ((fMotherStatus & AliPicoBase::kSecondaryFromWeakDecay) == AliPicoBase::kSecondaryFromWeakDecay);
  }

  Bool_t IsMotherSecondaryFromMaterial() const {
    return ((fMotherStatus & AliPicoBase::kSecondaryFromMaterial) == AliPicoBase::kSecondaryFromMaterial);
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

  ClassDef(AliPicoV0MC, 5);
};

#endif
