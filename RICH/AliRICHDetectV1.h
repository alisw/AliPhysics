#ifndef ALIRICHDETECTV1_H
#define ALIRICHDETECTV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////////////
//  Manager and hits classes for set:RICH default version //
////////////////////////////////////////////////////////////

#include "AliRICHDetect.h"

class AliRICHDetectV1 : public AliRICHDetect {
    
 public:

  AliRICHDetectV1();
  AliRICHDetectV1(const char *name, const char *title);
  virtual       ~AliRICHDetectV1();
  void   Detect(Int_t nev, Int_t type);
  Int_t  ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
  void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh);
    
  ClassDef(AliRICHDetectV1,1)  //Reconstruction module for :RICH version 1
	};


	
	
#endif
	






