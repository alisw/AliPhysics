#ifndef ALIRICHSEGMENTATIONV0_H
#define ALIRICHSEGMENTATIONV0_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliRICHSegmentation.h"

class AliRICHSegmentationV0 :
public AliRICHSegmentation {
 public:
    AliRICHSegmentationV0(){}
    virtual ~AliRICHSegmentationV0(){}
    //    
    // Set Chamber Segmentation Parameters
    //
    // Pad size Dx*Dy 
    virtual  void    SetPadSize(Float_t p1, Float_t p2);
    // Anod Pitch
    virtual  void    SetDAnod(Float_t D) {fWireD = D;};
        
    //
    // Transform from pad (wire) to real coordinates and vice versa
    //
    // Anod wire coordinate closest to xhit
    virtual Float_t GetAnod(Float_t xhit);
    // Transform from pad to real coordinates
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);
    // Transform from real to pad coordinates
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );
    //
    // Initialisation
    virtual void Init(AliRICHChamber* chamber);
    //
    // Get member data
    //
    // Pad size in x
    virtual Float_t Dpx(){return fDpx;}
    //
    // Pad size in y
    virtual Float_t Dpy(){return fDpy;}
    // Pad size in x by Sector
    virtual Float_t Dpx(Int_t) {return fDpx;}
    // Pad size in y by Sector
    virtual Float_t Dpy(Int_t) {return fDpy;}
    // Max number of Pads in x
    virtual Int_t   Npx(){return fNpx;}
    // Max number of Pads in y
    virtual Int_t   Npy(){return fNpy;}
    

    // set pad position
    virtual void     SetPad(Int_t ix, Int_t iy);
    // set hit position
    virtual void     SetHit(Float_t xhit , Float_t yhit);
    //
    // Iterate over pads
    // Initialiser
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);
    // Stepper
    virtual void  NextPad();
    // Condition
    virtual Int_t MorePads();
    //
    // Distance between 1 pad and a position
    virtual Float_t Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t *
				       dummy);
    // Number of pads read in parallel and offset to add to x 
    // (specific to LYON, but mandatory for display)
    virtual void GetNParallelAndOffset(Int_t iX, Int_t iY,
				       Int_t *Nparallel, Int_t *Offset) {*Nparallel=1;*Offset=0;}
    // Get next neighbours 
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
    //
    // Current Pad during Integration
    // x-coordinate
    virtual Int_t  Ix(){return fix;}
    // y-coordinate
    virtual Int_t  Iy(){return fiy;}
    // current sector
    virtual Int_t  ISector(){return 1;}
    // calculate sector from x-y coordinates
    virtual Int_t  Sector(Float_t x, Float_t y){return 1;}
    //
    // Signal Generation Condition during Stepping
    virtual Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
    // Initialise signal gneration at coord (x,y,z)
    virtual void  SigGenInit(Float_t x, Float_t y, Float_t z);
    // Current integration limits
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    // Test points for auto calibration
    virtual void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y);
    // Debugging utilities
    virtual void Draw();
    // Function for systematic corrections
    virtual void SetCorrFunc(Int_t dum, TF1* func) {fCorr=func;}
    
    virtual TF1* CorrFunc(Int_t) {return fCorr;} 
    ClassDef(AliRICHSegmentationV0,1)
	protected:
    //
    // Implementation of the segmentation data
    // Version 0 models rectangular pads with the same dimensions all
    // over the cathode plane
    //
    //  geometry
    //
    Float_t    fDpx;             // x pad width per sector  
    Float_t    fDpy;             // y pad base width
    Int_t      fNpx;             // Number of pads in x
    Int_t      fNpy;             // Number of pads in y
    Int_t      fSector;          // Current padplane
    Float_t    fWireD;           // wire pitch
    

        
    // Chamber region consideres during disintegration (lower left and upper right corner)
    //
    Int_t fixmin; // lower left  x
    Int_t fixmax; // lower left  y
    Int_t fiymin; // upper right x
    Int_t fiymax; // upper right y 
    //
    // Current pad during integration (cursor for disintegration)
    Int_t fix;  // pad coord. x 
    Int_t fiy;  // pad coord. y 
    Float_t fx; // x
    Float_t fy; // y
    //
    // Current pad and wire during tracking (cursor at hit centre)
    Float_t fxhit;   //x position
    Float_t fyhit;   //y position
    // Reference point to define signal generation condition
    Int_t fixt;     // pad coord. x
    Int_t fiyt;     // pad coord. y
    Int_t fiwt;     // wire number
    Float_t fxt;    // x
    Float_t fyt;    // y
    TF1*    fCorr;  // correction function
};
#endif
