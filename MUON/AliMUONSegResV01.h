#ifndef MUONSegResV01_H
#define MUONSegResV01_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////
//  Segmentation and Response classes version 01   //
/////////////////////////////////////////////////////
 
#include "AliMUON.h"
#include "AliMUONSegResV0.h"
#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"

class AliMUONSegmentationV01 :
public AliMUONSegmentationV0 {
 public:
    AliMUONSegmentationV01();
    virtual ~AliMUONSegmentationV01(){}
    //    
    // Set Chamber Segmentation Parameters
    // 
    virtual  void    SetPadDivision(Int_t ndiv[4]);
    // Radii
    virtual  void    SetSegRadii(Float_t  r[4]);
    //
    // Transform from pad (wire) to real coordinates and vice versa
    //
    // Transform from pad to real coordinates
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);
    // Transform from real to pad coordinates
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );    
    //
    // Initialisation
    virtual void Init(AliMUONChamber*);
    //
    // Get member data
    //
    // Pad size in x by Sector
    virtual Float_t Dpx(Int_t isec);
    // Pad size in y by Sector
    virtual Float_t Dpy(Int_t isec);
    // Max number of Pads in x
    virtual Int_t   Npx(){return fNpxS[fNsec-1][1]+1;}
    //
    virtual void     SetPad(Int_t,Int_t);
    //
    // Iterate over pads
    // Initialiser
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);
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
    // Test points for auto calibration
    void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y);
    //
    // Draw segmentation zones
    virtual void Draw();
    // Function for systematic corrections
    // Set the correction function
    virtual void SetCorrFunc(Int_t dum, TF1* func);
    // Get the correction function
    virtual TF1* CorrFunc(Int_t);
    ClassDef(AliMUONSegmentationV01,1) // Segmentation approximating circular zones with different pad size
 protected:
    //  Geometry
    //
    Int_t      fNsec;           // Number of sectors
    TArrayF    fRSec;           // Sector outer radia
    TArrayI    fNDiv;           // Pad size division
    TArrayF    fDpxD;           // y pad width per sector
    // Segmentation map
    Int_t      fNpxS[10][1000]; // Number of pads per sector in x
    Float_t    fCx[10][1000];   // pad-sector contour x vs y  
    // Chamber region consideres during disintegration
    // (lower left and upper right corner)
    //
    Float_t fxmin; // lower left  x
    Float_t fxmax; // lower left  y
    Float_t fymin; // upper right x
    Float_t fymax; // upper right y 

    //
    // Current pad during integration (cursor for disintegration)
    Int_t   fSector; // Current sector
    //
    TObjArray *fCorr; // Correction functions
};
#endif






