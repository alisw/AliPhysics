#ifndef ALIESDRESOLPARAMS_H
#define ALIESDRESOLPARAMS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//   ESD tracks and V0 resolution  parameterization
//
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"
#include <TVectorD.h>
class TTree;

class AliESDresolParams : public TObject{
 public:
  AliESDresolParams();
  void SetInstance(AliESDresolParams *param){fgInstance = param;}
  //
  Double_t GetResolPrimFast(Int_t param, Float_t onept, Float_t tanth) const;
  Double_t GetResolRFast(Int_t param, Float_t onept, Float_t radius) const;
  //
  static Double_t  SGetResolPrimFast(Int_t sdim, Float_t onept, Float_t tanth){ return fgInstance->GetResolPrimFast(sdim,onept,tanth);}
  static Double_t  SGetResolRFast(Int_t param, Float_t onept, Float_t radius){ return fgInstance->GetResolRFast(param,onept,radius);}
  void SetResolPrimFast(TObjArray* array);
  void SetResolRFast(TObjArray* array);
  // protected:
 public:
  //Resolution at prim vertex
  TVectorD *fResolDCAyy;            // resolution Y parameterization
  TVectorD *fResolDCAzz;            // resolution Z parameterization 
  TVectorD *fResolDCAphi;           // resolution phi parameterization - pt-theta
  TVectorD *fResolDCAth;            // resolution theta parameterization -pt-theta
  TVectorD *fResolDCA1pt;           // resolution 1/pt parameterization - pt-theta
  //
  // Resolution at V0 - radial dependent
  //
  TVectorD *fResolCyy;              // resolution Y parameterization - r-pt
  TVectorD *fResolCzz;              // resolution Z parameterization - r-pt
  TVectorD *fResolCphi;             // resolution phi parameterization - r-pt
  TVectorD *fResolCth;              // resolution theta parameterization - r-pt
  TVectorD *fResolC1pt;             // resolution 1/pt parameterization - r-pt  
  // 
  static AliESDresolParams*   fgInstance; //! Instance of this class (singleton implementation)
  ClassDef(AliESDresolParams,1)      // ESD resolution parametereization
};



#endif
