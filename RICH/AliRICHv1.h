#ifndef ALIRICHV1_H
#define ALIRICHV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


///////////////////////////////////////////////////////////
//  Manager and hits classes for set: RICH full version  //
///////////////////////////////////////////////////////////

#include "AliRICHv0.h"

class AliRICHv1 : public AliRICHv0 {
    
 public:
    
  //Int_t fCkov_number;
  //Int_t fFreon_prod;

    AliRICHv1();
    AliRICHv1(const char *name, const char *title);
    virtual void   Init();
    virtual       ~AliRICHv1() {}
    virtual void   StepManager();

 private:
    ClassDef(AliRICHv1,1)  //Hits manager for set: RICH full version 
	
	};
	
	
#endif
	






