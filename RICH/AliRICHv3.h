#ifndef AliRICHv3_h
#define AliRICHv3_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRICH.h"

class AliRICHv3 : public AliRICH 
{    
public:
    
 inline                AliRICHv3():AliRICH()                                {} 
                       AliRICHv3(const char *pcName, const char *pcTitle);    
        virtual       ~AliRICHv3();                             
 inline virtual Int_t  IsVersion()                                     const{return 3;}
        virtual void   StepManager();     
        virtual void   CreateGeometry();  
        virtual void   BuildGeometry();   
        virtual void   Init();            // Makes nothing for a while          
private:
            Double_t* RotateXY(const Double_t* r, Double_t a);   //Rotation in the X-Y plane in G3 notation
  ClassDef(AliRICHv3,1)  //RICH full version, configurable with azimuthal rotation	
};// class AliRICHv3

#endif // AliRICHv3_h
