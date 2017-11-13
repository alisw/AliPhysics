// -*- C++ -*-
#ifndef ALI_AD_CHARGE_EQUALIZATION_H
#define ALI_AD_CHARGE_EQUALIZATION_H

#include <TNamed.h>

class AliVAD;

class AliADChargeEqualization : public TNamed {
public:
  enum Side {
    kCSide     = 0,
    kASide     = 1,
    kBothSides = 2
  };
  AliADChargeEqualization(Float_t *quantiles=NULL);
  virtual ~AliADChargeEqualization() {}

  virtual void Print(Option_t *option="") const;

  Float_t GetMult(const AliVAD*, Side) const;
protected:
private:
  Float_t fNormAll[16];       //
  Float_t fNormPerSide[2][8]; //

  ClassDef(AliADChargeEqualization, 1);
} ;

#endif // ALI_AD_CHARGE_EQUALIZATION_H
