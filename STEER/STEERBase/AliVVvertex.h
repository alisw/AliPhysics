#ifndef ALIVVVERTEX_H
#define ALIVVVERTEX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Author: Mikolaj Krzewicki

 * >> interface for vertices <<
 */

#include "Rtypes.h"

class AliVVvertex
{
  public:
  AliVVvertex() {} 
  virtual ~AliVVvertex() {}
  virtual Double_t GetX() const { return 0; }
  virtual Double_t GetY() const { return 0; }
  virtual Double_t GetZ() const { return 0; }
  virtual void GetXYZ(Double_t pos[3]) const {if (pos[0]<0) return;}
  virtual void     GetCovarianceMatrix(Double_t covmatrix[6]) const = 0;
  virtual Double_t GetChi2perNDF() const = 0;
  virtual Double_t GetChi2() const = 0;
  virtual Int_t    GetNDF() const = 0;
  virtual void     PrintIndices() const = 0;
  virtual void     Print(Option_t* option = "") const = 0;
  virtual Int_t    GetBC() const = 0;
  virtual void Clear(Option_t* /*option*/) = 0;
  virtual Int_t    GetNContributors() const { return 0; }

  //  ClassDef(AliVVvertex, 1)   

};

#endif
