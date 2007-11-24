#ifndef ALI_TPC_TRANSFORM_H
#define ALI_TPC_TRANSFORM_H

#include "TGeoMatrix.h"
#include "AliTransform.h"

class AliTPCTransform:public AliTransform {
public:
  AliTPCTransform();
  virtual ~AliTPCTransform();
  virtual void Transform(Double_t *x,Int_t *i,UInt_t time,
			 Int_t coordinateType);
protected:
  void Local2RotatedGlobal(Int_t sec,  Double_t *x) const;
  void RotatedGlobal2Global(Int_t sector,Double_t *x) const;
  void Global2RotatedGlobal(Int_t sector,Double_t *x) const;
  void GetCosAndSin(Int_t sector,Double_t &cos,Double_t &sin) const;
private:
  Double_t fCoss[18];
  Double_t fSins[18];

  ClassDef(AliTPCTransform,1)
};

#endif
