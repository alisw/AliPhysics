#ifndef ALIMUONSEGMENTATIONV02_H
#define ALIMUONSEGMENTATIONV02_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation and Response classes version 02   //
/////////////////////////////////////////////////////
 

#include "AliMUONSegmentationV01.h"

class AliMUONSegmentationV02 :
public AliMUONSegmentationV01 {
 public:
    AliMUONSegmentationV02(){};
    virtual ~AliMUONSegmentationV02(){}
    //
    // Pad size Dx*Dy 
    virtual void SetPadSize(Float_t p1, Float_t p2);
    //
    // Get member data
    // Pad size in x
    virtual Float_t Dpx() const {return fDpy;}
    // Pad size in y
    virtual Float_t Dpy() const {return fDpx;}
    // Pad size in x by Sector
    virtual Float_t Dpx(Int_t isec) const;
    // Pad size in y by Sector
    virtual Float_t Dpy(Int_t isec) const;
    // Max number of Pads in x
    virtual Int_t   Npx()  const;
     // max number of Pads in y
    virtual Int_t   Npy()  const;
    // calculate sector from pad coordinates
    virtual Int_t   Sector(Int_t ix, Int_t iy);
    virtual void Draw(const char *opt="") const {}
    //
    // Transform from pad (wire) to real coordinates and vice versa
    // Transform from pad to real coordinates
    virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x ,Float_t &y );
    virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z) 
	{z=fZ; GetPadC(ix, iy, x , y);}
    // Transform from pad to real coordinates
    virtual void    GetPadI(Float_t x ,Float_t y , Int_t &ix, Int_t &iy);
    virtual void    GetPad(Float_t x, Float_t y , Float_t z, Int_t &ix, Int_t &iy) 
	{GetPadI(x, y, ix, iy);}
    // Set pad position
    virtual void    SetPad(Int_t ix,Int_t iy);
    // Stepper
    virtual void    NextPad();
    // Condition
    virtual Int_t   MorePads();
    // Get next neighbours 
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
    // Get next neighbours 
    ClassDef(AliMUONSegmentationV02,1)  // Segmentation approximating circular zones with different pad size
	};
#endif









