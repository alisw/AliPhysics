#ifndef RICHSegResV1_H
#define RICHSegResV1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliRICHSegResV0.h"

class AliRICHSegmentationV1 : public AliRICHSegmentationV0 {
    
 public:
    AliRICHSegmentationV1();
    virtual ~AliRICHSegmentationV1();
    
    Int_t fSector;

    // current sector
    virtual Int_t  ISector(){return fSector;}
    // calculate sector from x-y coordinates
    virtual Int_t  Sector(Float_t x, Float_t y);

    // Transform from pad to real coordinates
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);
    // Transform from real to pad coordinates
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );
    // Current integration limits
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    
    ClassDef(AliRICHSegmentationV1,1)
	
	
	};
	
#endif
