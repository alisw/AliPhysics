#ifndef ALIRICHV1_H
#define ALIRICHV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


///////////////////////////////////////////////////////////
//  Manager and hits classes for set: RICH full version  //
///////////////////////////////////////////////////////////

#include "AliRICHv0.h"

class AliRICHv1 : public AliRICH {
    
 public:
    
  //Int_t fCkov_number;
  //Int_t fFreon_prod;

    AliRICHv1();
    AliRICHv1(const char *name, const char *title);
    virtual void   Init();
    virtual Int_t  IsVersion() const {return 1;}
    virtual       ~AliRICHv1() {}

 private:
    ClassDef(AliRICHv1,1)  //Hits manager for set: RICH full version 
	
	};
	
	
#endif
	






