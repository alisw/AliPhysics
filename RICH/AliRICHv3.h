#ifndef AliRICHv3_h
#define AliRICHv3_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRICH.h"

class AliRICHv3 : public AliRICH 
{
    
public:
    
   AliRICHv3():AliRICH()                                {} // Default ctor
   AliRICHv3(const char *pcName, const char *pcTitle);     // Named ctor 
   virtual       ~AliRICHv3();                             // Dtor
   
   virtual Int_t  IsVersion()                     const {return 3;}
      
   virtual void   CreateMaterials(); // Provides material definition for simulation (currently GEANT)   
   virtual void   CreateGeometry();  // Provides geometry structure for simulation (currently GEANT volumes tree)
   virtual void   BuildGeometry();   // Provides geometry structure for event display (ROOT TNode tree)
   virtual void   Init();            // Makes nothing for a while
   virtual void   StepManager();     // Processes the particle transport 
       
private:
   Double_t* RotateXY(const Double_t* r, Double_t a);
   
    ClassDef(AliRICHv3,1)  //RICH full version, configurable with azimuthal rotation	
};// class AliRICHv3

#endif // AliRICHv3_h
