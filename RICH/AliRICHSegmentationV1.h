#ifndef AliRICHSegmentationV1_h
#define AliRICHSegmentationV1_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRICHSegmentationV0.h"

class AliRICHSegmentationV1 : public AliRICHSegmentationV0 
{    
public:
            AliRICHSegmentationV1();     
   virtual ~AliRICHSegmentationV1() {}   
   virtual void   Init(Int_t id); 
   virtual Int_t  Sector(Float_t x, Float_t y);    
   virtual void   GetPadI(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);    
   virtual void   GetPadI(Float_t x, Float_t y , Float_t /*z*/, Int_t &ix, Int_t &iy)  {GetPadI(x, y, ix, iy);}
   virtual void   GetPadC(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );    
   virtual void   GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z) {z=0; GetPadC(ix, iy, x , y);}    
   virtual void IntegrationLimits (Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2); 
   virtual Int_t  ISector() const{return fSector;}   
private:
    ClassDef(AliRICHSegmentationV1,1)
};
#endif//AliRICHSegmentationV1_h
