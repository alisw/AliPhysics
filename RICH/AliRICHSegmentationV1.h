#ifndef ALIRICHSEGMENTATIONV1_H
#define ALIRICHSEGMENTATIONV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliRICHSegmentationV0.h"

class AliRICHSegmentationV1 : public AliRICHSegmentationV0 {
    
 public:
    AliRICHSegmentationV1();
    virtual ~AliRICHSegmentationV1();
    // current sector
    virtual Int_t  ISector() {return fSector;}
    // calculate sector from x-y coordinates
    virtual Int_t  Sector(Float_t x, Float_t y);

    // Transform from pad to real coordinates
    virtual void    GetPadI(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);
    virtual void    GetPadI(Float_t x, Float_t y , Float_t z, Int_t &ix, Int_t &iy)  
	{GetPadI(x, y, ix, iy);}
    // Transform from real to pad coordinates
    virtual void    GetPadC(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );
    virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z) 
	{z=0; GetPadC(ix, iy, x , y);}
    // Current integration limits
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
 private:
    Int_t fSector;             //Pad plane sector
    ClassDef(AliRICHSegmentationV1,1)
};
	
#endif








