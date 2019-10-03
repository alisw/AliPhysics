#ifndef ALIESDRESOLMAKERFAST_H
#define ALIESDRESOLMAKERFAST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//   ESD tracks and V0 resolution  parameterization maker
//
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"
#include <TVectorD.h>
class TTree;
class TObjArray; 
class TCut;
//
class AliESDresolMakerFast : public TObject{
 public:
  AliESDresolMakerFast();
  //
  static TObjArray * MakeParamPrimFast(TTree * tree, TCut &cutDCA, Float_t fraction=-1, Int_t entries=100000);
  static TObjArray * MakeParamRFast(TTree * tree, TCut &cutV0, Float_t fraction=-1, Int_t entries=100000);
  // protected:
 public:
  // 
  ClassDef(AliESDresolMakerFast,1)      // ESD resolution parametereization
};



#endif
