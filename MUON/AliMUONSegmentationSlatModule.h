#ifndef ALIMUONSEGMENTATIONSLATMODULE_H
#define ALIMUONSEGMENTATIONSLATMODULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation classes for slat modules          //
//  to be used with AluMUONSegmentationSlat        //
/////////////////////////////////////////////////////

class TArrayF;
class TArrayI;
class TObjArray;

#include  "AliMUONSegmentationV0.h"

class AliMUONSegmentationSlatModule :
public AliMUONSegmentationV0 {
 public:
    AliMUONSegmentationSlatModule();
    virtual ~AliMUONSegmentationSlatModule(){}
    //    
    // Set Chamber Segmentation Parameters
    // 
    virtual  void    SetPadDivision(Int_t ndiv[4]);
    // Transform from pad to real coordinates
    virtual void    GetPadI(Float_t x ,Float_t y ,Int_t   &ix,Int_t &iy);
    virtual void    GetPadI(Float_t x, Float_t y , Float_t z, Int_t &ix, Int_t &iy)
	{GetPadI(x, y, ix, iy);}
    // Transform from real to pad coordinates
    virtual void    GetPadC(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );
    virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z)
	{z=0; GetPadC(ix, iy, x , y);}
    // Initialisation
    virtual void Init(Int_t chamber);
    // Set Segmentation Zones (PCB Boards)
    virtual void SetPcbBoards(Int_t n[4]);
    //
    // Get member data
    //
    // Pad size in x by Sector
    virtual Float_t Dpx(Int_t isec) const;
    // Pad size in y by Sector
    virtual Float_t Dpy(Int_t isec) const;
    //
    virtual void    SetPad(Int_t ix,Int_t iy);
    virtual void    SetHit(Float_t xhit, Float_t yhit);
    virtual void    SetHit(Float_t xhit, Float_t yhit, Float_t zhit)
	{SetHit(xhit, yhit);}
    //
    // Iterate over pads
    // Initialiser
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t zhit, Float_t dx, Float_t dy)
	{FirstPad(xhit, yhit, dx, dy);}
    // Stepper
    virtual void  NextPad();
    // Condition
    virtual Int_t MorePads();
    // Get next neighbours 
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]); 
    //
    // Current Pad during Integration
    // current sector
    virtual Int_t ISector()  {return fSector;}
    // calculate sector from pad coordinates
    virtual Int_t Sector(Int_t ix, Int_t iy);
    //
    // Integration
    // Current integration limits
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    //
    virtual void SetId(Int_t id) {fId=id;}
    
 protected:
    //
    //  Geometry
    //
    Int_t       fId;             // Id of this module
    Int_t       fNsec;           // Number of sectors
    TArrayI*    fNDiv;           // Pad size division
    TArrayF*    fDpxD;           // y pad width per sector
    
    // Segmentation map
    Int_t      fNpxS[10];       // Number of pads per sector in x
    Int_t      fNpyS[10];       // Number of pads per sector in y    
    Float_t    fCx[10];         // pad-sector contour x vs y  
    // Chamber region consideres during disintegration
    // (lower left and upper right corner)
    //
    Float_t fXmin;           // lower left  x
    Float_t fXmax;           // lower left  y
    Float_t fYmin;           // upper right x
    Float_t fYmax;           // upper right y 

    //
    // Current pad during integration (cursor for disintegration)
    Int_t   fSector;         // Current sector

    Float_t fDxPCB;          // x-size of PCB board
    Float_t fDyPCB;          // y-size of PCB board
    Int_t   fPcbBoards[4];   // number of PCB boards per segmentation region
    
    ClassDef(AliMUONSegmentationSlatModule,1) // Segmenation for bending plate slats
	
};
#endif






