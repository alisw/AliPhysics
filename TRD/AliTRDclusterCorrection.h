#ifndef ALITRDCLUSTERCORRECTON_H
#define ALITRDCLUSTERCORRECTON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "TObject.h"

class AliTRDclusterCorrection : public TObject {

 public:  
  AliTRDclusterCorrection();
  Float_t GetCorrection(Int_t plane, Int_t timebin, Float_t angle);
  Float_t GetSigma(Int_t plane, Int_t timebin, Float_t angle);
  Float_t GetOffsetAngle(){return fOffsetAngle;}
  void SetOffsetAngle(Float_t angle){fOffsetAngle=angle;}
  void SetCorrection(Int_t plane,Int_t timebin, Float_t angle, Float_t value,Float_t sigma);
  Float_t GetAngle(Int_t i){return (i-10.)/10.+fOffsetAngle;}
  static  AliTRDclusterCorrection * GetCorrection();
 protected:
  Float_t  fCorrections[6][30][20][2];
  Float_t fOffsetAngle;
				 
  ClassDef(AliTRDclusterCorrection,1) // ClusterCorrection for the TRD
 
};

#endif
