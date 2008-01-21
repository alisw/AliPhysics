#ifndef ALITPCTRANSFORM_H
#define ALITPCTRANSFORM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
//    Class for tranformation of the coordinate frame
//    Transformation  
//      local coordinate frame (sector, padrow, pad, timebine) ==>
//      rotated global (tracking) cooridnate frame (sector, lx,ly,lz)
//

#include "AliTransform.h"

class AliTPCTransform:public AliTransform {
public:
  AliTPCTransform();
  virtual ~AliTPCTransform();
  virtual void Transform(Double_t *x,Int_t *i,UInt_t time,
			 Int_t coordinateType);
  void SetPrimVertex(Double_t *vtx);
protected:
  void Local2RotatedGlobal(Int_t sec,  Double_t *x) const;
  void RotatedGlobal2Global(Int_t sector,Double_t *x) const;
  void Global2RotatedGlobal(Int_t sector,Double_t *x) const;
  void GetCosAndSin(Int_t sector,Double_t &cos,Double_t &sin) const;
private:
  Double_t fCoss[18];  // cache the transformation
  Double_t fSins[18];  // cache the transformation
  Double_t fPrimVtx[3];// position of the primary vertex - needed for TOF correction
  ClassDef(AliTPCTransform,1)
};

#endif
