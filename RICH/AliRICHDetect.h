#ifndef AliRICHDetect_H
#define AliRICHDetect_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//   Reconstruction classes for set:RICH version 0     //
/////////////////////////////////////////////////////////

#include "AliRICH.h"

class AliRICHDetect;

class AliRICHDetect : public TObject {
    
 public:
  AliRICHDetect();
  AliRICHDetect(const char *name, const char *title);
  virtual       ~AliRICHDetect() {}
  void   Detect();
  float Area(float theta,float OMEGA);
  //virtual void AddRecHit(const AliRICHRecHit);
  //virtual int  ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
  ClassDef(AliRICHDetect,1)  //Reconstruction module for :RICH version 0
      
      
	
	};


	
	
#endif
