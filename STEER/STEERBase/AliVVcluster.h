#ifndef ALIVVCLUSTER_H
#define ALIVVCLUSTER_H

#include "Rtypes.h"

class AliVVcluster
{
  public:
  AliVVcluster() {}
  AliVVcluster(Bool_t) {}
  virtual ~AliVVcluster() {}
  virtual void SetX(Float_t /*x*/)             {}
  virtual void SetY(Float_t /*y*/)             {}
  virtual void SetZ(Float_t /*z*/)             {}
  virtual void SetPadRow(Short_t /*padrow*/)   {}
  virtual void SetSigmaY2(Float_t /*sigmaY2*/) {}
  virtual void SetSigmaZ2(Float_t /*sigmaZ2*/) {}
  virtual void SetCharge(UShort_t /*charge*/)  {}
  virtual void SetQMax(UShort_t /*qmax*/)      {}

  virtual Float_t  GetX()       const      {return 0.;}
  virtual Float_t  GetY()       const      {return 0.;}
  virtual Float_t  GetZ()       const      {return 0.;}
  virtual UShort_t GetPadRow()  const      {return 0;}
  virtual Float_t  GetSigmaY2() const      {return 0.;}
  virtual Float_t  GetSigmaZ2() const      {return 0.;}
  virtual UShort_t GetCharge()  const      {return 0;}
  virtual UShort_t GetQMax()    const      {return 0;}
};

#endif
