#ifndef MUONSegResV01_H
#define MUONSegResV01_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation and Response classes version 01   //
/////////////////////////////////////////////////////
 
#include "AliMUON.h"
#include "AliMUONSegResV0.h"
#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"

class AliMUONsegmentationV01 :
public AliMUONsegmentationV0 {
 public:
    AliMUONsegmentationV01();
    virtual ~AliMUONsegmentationV01(){}
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
    virtual void Init(AliMUONchamber*);
    //
    // Get member data
    //
    // Pad size in x
    virtual Float_t Dpx(){return AliMUONsegmentationV0::Dpx();}
    // Pad size in y
    virtual Float_t Dpy(){return AliMUONsegmentationV0::Dpy();}
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
    // Debugging utilities
    virtual void Draw(Option_t *);
    // Function for systematic corrections
    virtual void SetCorrFunc(Int_t dum, TF1* func);
    virtual TF1* CorrFunc(Int_t);
    

    ClassDef(AliMUONsegmentationV01,1)
 protected:
    //
    // Implementation of the segmentation data
    // Version 0 models rectangular pads with the same dimensions all
    // over the cathode plane
    //
    //  geometry
    //
    Int_t      fNsec;           // Number of sectors
    TArrayF    fRSec;           // sector outer radia
    TArrayI    fNDiv;           // pad size division
    TArrayF    fDpxD;           // y pad width per sector
    Int_t      fNpxS[10][1000]; // Number of pads per sector in x
    Float_t    fCx[10][1000];   // pad-sector contour x vs y  
    // Chamber region consideres during disintegration
    // (lower left and upper right corner)
    //
    Float_t fxmin;
    Float_t fxmax;
    Float_t fymin;
    Float_t fymax;

    //
    // Current pad during integration (cursor for disintegration)
    Int_t   fSector;
    TObjArray *fCorr;
};
#endif






