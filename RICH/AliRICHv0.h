#ifndef ALIRICHV0_H
#define ALIRICHV0_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


////////////////////////////////////////////////////////////
//  Manager and hits classes for set:RICH default version //
////////////////////////////////////////////////////////////

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
    Float_t        Fresnel(Float_t ene,Float_t pdoti, Bool_t pola);
    Float_t        AbsoCH4(Float_t x);
//   virtual void   Trigger(Float_t (*)[4], Float_t (*)[4], Int_t& iflag);
 protected:
    Int_t fCkovNumber;                   // Number of Cerenkov photons
    Int_t fCkovQuarz;                    // Cerenkovs crossing quartz
    Int_t fCkovGap;                      // Cerenkovs crossing gap
    Int_t fCkovCsi;                      // Cerenkovs crossing csi
    Int_t fLostRfreo;                    // Cerenkovs reflected in freon
    Int_t fLostRquar;                    // Cerenkovs reflected in quartz
    Int_t fLostAfreo;                    // Cerenkovs absorbed in freon 
    Int_t fLostAquarz;                   // Cerenkovs absorbed in quartz
    Int_t fLostAmeta;                    // Cerenkovs absorbed in methane
    Int_t fLostCsi;                      // Cerenkovs below csi quantum efficiency 
    Int_t fLostWires;                    // Cerenkovs lost in wires
    Int_t fFreonProd;                    // Cerenkovs produced in freon
    Float_t fMipx;                       // x coord. of MIP
    Float_t fMipy;                       // y coord. of MIP
    Int_t fFeedbacks;                    // Number of feedback photons
    Int_t fLostFresnel;                  // Cerenkovs lost by Fresnel reflection
    ClassDef(AliRICHv0,1)  //Hits manager for set: RICH default version
	
	};
	
	
#endif
	






