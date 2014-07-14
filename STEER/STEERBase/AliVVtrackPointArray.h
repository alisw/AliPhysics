#ifndef ALIVVTRACKPOINTARRAY_H
#define ALIVVTRACKPOINTARRAY_H

#include "Rtypes.h"

#include "AliVVtrackPoint.h"

class AliVVtrackPointArray {
  public:
  AliVVtrackPointArray() {}
  virtual ~AliVVtrackPointArray() {}
  virtual Int_t GetNPoints() const {return 0;}
  virtual Bool_t GetPoint(AliVVtrackPoint& /*p*/, Int_t /*i*/) const =0 ;

  ClassDef(AliVVtrackPointArray, 1);
};

#endif
