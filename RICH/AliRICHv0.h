#ifndef RICHv0_H
#define RICHv0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:RICH version 0    //
/////////////////////////////////////////////////////////

#include "AliRICH.h"

class AliRICHv0 : public AliRICH {
    
 public:
    AliRICHv0();
    AliRICHv0(const char *name, const char *title);
    virtual       ~AliRICHv0() {}
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init();
    virtual Int_t  IsVersion() const {return 0;}
    virtual void   StepManager();
//   virtual void   Trigger(Float_t (*)[4], Float_t (*)[4], Int_t& iflag);
 private:
    ClassDef(AliRICHv0,1)  //Hits manager for set:RICH version 0
	
	};
	
	
#endif
	






