#ifndef AliRICHSegmentationV1_h
#define AliRICHSegmentationV1_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliRICHSegmentationV0.h"
#include <Riostream.h>

class AliRICHSegmentationV1 : public AliRICHSegmentationV0 
{    
public:
// ctor $ dtor:      
            AliRICHSegmentationV1();      // default ctor
   virtual ~AliRICHSegmentationV1() {}    // dtor
// The following staff is defined in AliRICHSegmentation.cxx:
   virtual void   Init(Int_t id); // Recalculates all the values after some of them have been changed

   virtual Int_t  Sector(Float_t x, Float_t y);    // calculate sector from x-y coordinates

   virtual void    GetPadI(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);         // Transform from pad to real coordinates
   virtual void    GetPadI(Float_t x, Float_t y , Float_t /*z*/, Int_t &ix, Int_t &iy)  {GetPadI(x, y, ix, iy);}
    
   virtual void    GetPadC(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );    // Transform from real to pad coordinates
   virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z) {z=0; GetPadC(ix, iy, x , y);}
    
   virtual void IntegrationLimits (Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);    // Current integration limits
// inline methods:
   virtual Int_t  ISector() const{return fSector;}   // Get current sector
   
   inline virtual void Print(Option_t *option)const; // Prints debug information
    
private:
    ClassDef(AliRICHSegmentationV1,1)
};

inline void AliRICHSegmentationV1::Print(Option_t *option)const
{
   TObject::Print(option);
   cout<<"Pad width in cm:         "<<fDpx                 <<endl;
   cout<<"Pad heights in cm:       "<<fDpy                 <<endl;
   cout<<"Pad number along x:      "<<fNpx                 <<endl;
   cout<<"Pad number along y:      "<<fNpy                 <<endl;
   cout<<"Sector:                  "<<fSector              <<endl;
   cout<<"Wire pitch:              "<<fWireD               <<endl; 
   cout<<"Dead zone in cm:         "<<fDeadZone            <<endl;
   cout<<"Pad plane width in cm:   "<<fPadPlane_Width      <<endl;
   cout<<"Pad plane heights in cm: "<<fPadPlane_Length     <<endl;
}//void AliRICHSegmentationV1::Print(Option_t *option)const
	
#endif//AliRICHSegmentationV1_h
