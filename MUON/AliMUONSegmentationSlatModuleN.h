#ifndef ALIMUONSEGMENTATIONSLATMODULEN_H
#define ALIMUONSEGMENTATIONSLATMODULEN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation classes for slat modules          //
//  to be used with AluMUONSegmentationSlat        //
/////////////////////////////////////////////////////

#include "AliMUONSegmentationSlatModule.h"

class  AliMUONSegmentationSlatModuleN :
public AliMUONSegmentationSlatModule {
 public:
    AliMUONSegmentationSlatModuleN();
    AliMUONSegmentationSlatModuleN(Int_t nsec);
    virtual ~AliMUONSegmentationSlatModuleN(){}
    // Transform from pad to real coordinates
    virtual void    GetPadI(Float_t x ,Float_t y ,Int_t   &ix,Int_t &iy);
    // Transform from real to pad coordinates
    virtual void    GetPadC(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );
    // Initialisation
    virtual void Init(Int_t chamber);
    //
    // Get member data
    //
    // Pad size in x by Sector
    virtual Float_t Dpx(Int_t isec) const;
    // Pad size in y by Sector
    virtual Float_t Dpy(Int_t isec) const;
    // Iterate over pads
    // Stepper
    virtual void  NextPad();
    virtual Int_t MorePads();
    
    // Get next neighbours
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
 protected:
    Int_t   fNpxPCB;        // Number of strips per PCB board 
    
    ClassDef(AliMUONSegmentationSlatModuleN,1) // Segmentation class for non bending plane slat
};

#endif






