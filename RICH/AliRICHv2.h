#ifndef ALIRICHV2_H
#define ALIRICHV2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


///////////////////////////////////////////////////////////
//  Manager and hits classes for set: RICH full version  //
///////////////////////////////////////////////////////////

#include "AliRICHv1.h"

class AliRICHv2 : public AliRICHv1 {
    
 public:
    
  //Int_t fCkov_number;
  //Int_t fFreon_prod;

    AliRICHv2();
    AliRICHv2(const char *name, const char *title);
    virtual void   Init();
    virtual       ~AliRICHv2() {}
   
 private:
    ClassDef(AliRICHv2,1)  //Hits manager for set: RICH full version, CONFIGURABLE 
	
	};
	
	
#endif
	






