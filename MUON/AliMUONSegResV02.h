#ifndef MUONv02_H
#define MUONv02_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation and Response classes version 01   //
/////////////////////////////////////////////////////
 
#include "AliMUON.h"
#include "TArrayF.h"
#include "TArrayI.h"
#include "AliMUONSegResV01.h"
class AliMUONsegmentationV02 :
public AliMUONsegmentationV01 {
 public:
    AliMUONsegmentationV02(){};
    virtual ~AliMUONsegmentationV02(){}
    //
    virtual void SetPADSIZ(Float_t p1, Float_t p2);
    //
    // Get member data
    // Pad size in x
    virtual Float_t Dpx() {return fDpy;}
    // Pad size in y
    virtual Float_t Dpy() {return fDpx;}
    // Pad size in x by Sector
    virtual Float_t Dpx(Int_t isec);
    // Pad size in y by Sector
    virtual Float_t Dpy(Int_t isec);
    // Max number of Pads in x
    virtual Int_t   Npx();
     // max number of Pads in y
    virtual Int_t   Npy();
    // calculate sector from pad coordinates
    virtual Int_t   Sector(Int_t ix, Int_t iy);
    //
    // Transform from pad (wire) to real coordinates and vice versa
    // Transform from pad to real coordinates
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );
    // Transform from pad to real coordinates
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);
    virtual void    SetPad(Int_t ix,Int_t iy);
    // Stepper
    virtual void    NextPad();
    // Condition
    virtual Int_t   MorePads();
    
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
    // Get next neighbours 
    ClassDef(AliMUONsegmentationV02,1) //Muon chamber segmentation version 02
	};
#endif









