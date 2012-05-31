#ifndef ALI_TRANSFORM_H
#define ALI_TRANSFORM_H

#include "TObject.h"

class AliTransform:public TObject {
public:
  virtual ~AliTransform() {};
  virtual void Transform(Double_t *x,Int_t *i,UInt_t time,
			 Int_t coordinateType)=0;

  ClassDef(AliTransform,0)
};

#endif
