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
    Int_t fCkov_number;
    Int_t fCkov_quarz;
    Int_t fCkov_gap;
    Int_t fCkov_csi;
    Int_t fLost_rfreo;
    Int_t fLost_rquar;
    Int_t fLost_afreo;
    Int_t fLost_aquarz;
    Int_t fLost_ameta;
    Int_t fLost_csi;
    Int_t fLost_wires;
    Int_t fFreon_prod;
    Float_t fMipx;
    Float_t fMipy;
    Int_t fFeedbacks;
    Int_t fLost_fresnel;
    
    AliRICHv0();
    AliRICHv0(const char *name, const char *title);
    virtual       ~AliRICHv0() {}
    virtual void   CreateGeometry();
    virtual void   CreateMaterials();
    virtual void   Init();
    virtual Int_t  IsVersion() const {return 0;}
    virtual void   StepManager();
    Float_t        Fresnel(Float_t ene,Float_t pdoti, Bool_t pola);
    Float_t        AbsoCH4(Float_t x);
//   virtual void   Trigger(Float_t (*)[4], Float_t (*)[4], Int_t& iflag);
 private:
    ClassDef(AliRICHv0,1)  //Hits manager for set:RICH version 0
	
	};
	
	
#endif
	






