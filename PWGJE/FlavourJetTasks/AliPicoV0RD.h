#ifndef ALIPICOV0RD_H
#define ALIPICOV0RD_H

#include "AliPicoV0.h"

class AliPicoV0RD : public AliPicoV0 {

 public :

  AliPicoV0RD();
  AliPicoV0RD(UInt_t   wMask,
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
              Float_t  dNegPionSigmaTPC, Float_t dNegProtonSigmaTPC);
  AliPicoV0RD(const AliPicoV0RD &src);
  AliPicoV0RD& operator=(const AliPicoV0RD &src);
  virtual ~AliPicoV0RD();

  virtual Bool_t IsKshort(Double_t const dCuts[10]) const;
  virtual Bool_t IsLambda(Double_t const dCuts[10]) const;
  virtual Bool_t IsAntiLa(Double_t const dCuts[10]) const;

  virtual void GetControlVariables(Float_t d[20]) const;
//=============================================================================

 private :

  Float_t fPosPionSigmaTPC;  //
  Float_t fNegPionSigmaTPC;  //

  Float_t fPosProtonSigmaTPC;  //
  Float_t fNegProtonSigmaTPC;  //

  ClassDef(AliPicoV0RD, 5);
};

#endif
