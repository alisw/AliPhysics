#ifndef RICHSegRes_H
#define RICHSegRes_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"
#include "TClonesArray.h"
#include "TF1.h"
class AliRICHChamber;

class AliRICHSegmentation :
public TObject {
    
 public:
    
    // Set Chamber Segmentation Parameters
    //
    // Pad size Dx*Dy 
    virtual void    SetPadSize(Float_t p1, Float_t p2)  =0;
    // Anod Pitch
    virtual void    SetDAnod(Float_t D)                 =0;
    
    //
    // Anod wire coordinate closest to xhit
    virtual Float_t GetAnod(Float_t xhit)               =0;
    // Transform from pad (wire) to real coordinates
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy)=0;
    // Transform from real to pad coordinates
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y )=0;
    //
    // Initialisation
    virtual void Init(AliRICHChamber*)                  =0;
    //
    // Get member data
    //
    // Pad size in x
    virtual Float_t Dpx()                               =0;
    // Pad size in y
    virtual Float_t Dpy()                               =0;
    // Pad size in x by Sector 
    virtual Float_t Dpx(Int_t)                          =0;
    // Pad size in y by Sector 
    virtual Float_t Dpy(Int_t)                          =0;
    // Max number of Pads in x
    virtual Int_t Npx()                                 =0;
    // Max number of Pads in y
    virtual Int_t Npy()                                 =0;
    

    // set pad position
    virtual void     SetPad(Int_t, Int_t)               =0;
    // set hit position
    virtual void     SetHit(Float_t, Float_t)           =0;
    //
    // Iterate over pads
    // Initialiser
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy) =0;
    // Stepper
    virtual void  NextPad()=0;
    // Condition
    virtual Int_t MorePads()                           =0;
    //
    // Distance between 1 pad and a position
    virtual Float_t Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t *dummy) =0;
    // Number of pads read in parallel and offset to add to x 
    // (specific to LYON, but mandatory for display)
    virtual void GetNParallelAndOffset(Int_t iX, Int_t iY,
				       Int_t *Nparallel, Int_t *Offset) =0;
    // Get next neighbours 
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10])     =0;
    //
    // Current pad cursor during disintegration
    // x-coordinate
    virtual Int_t  Ix()                                =0;
    // y-coordinate
    virtual Int_t  Iy()                                =0;
    // current Sector
    virtual Int_t  ISector()                           =0;
    // calculate sector from pad coordinates
    virtual Int_t  Sector(Float_t ix, Float_t iy)          =0;
    //
    // Signal Generation Condition during Stepping
    virtual Int_t SigGenCond(Float_t x, Float_t y, Float_t z) = 0;
    // Initialise signal gneration at coord (x,y,z)
    virtual void  SigGenInit(Float_t x, Float_t y, Float_t z) = 0;
    // Current integration limits 
    virtual void  IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2)  = 0;
    // Test points for auto calibration
    virtual void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y) = 0;
    // Debug utilities
    virtual void Draw()                                           = 0;
    // Function for systematic corrections
    virtual void SetCorrFunc(Int_t, TF1*)                         = 0;
    virtual TF1* CorrFunc(Int_t)                                  = 0;
    ClassDef(AliRICHSegmentation,1)
	};
//----------------------------------------------
//
// Chamber response virtual base class
//
class AliRICHResponse :
public TObject {
 public:
    //
    // Configuration methods
    //
    // Number of sigmas over which cluster didintegration is performed
    virtual void    SetSigmaIntegration(Float_t p1)           =0;
    virtual Float_t SigmaIntegration()                        =0;
    // charge slope in ADC/e
    virtual void    SetChargeSlope(Float_t p1)                =0;
    virtual Float_t ChargeSlope()                             =0;
    // sigma of the charge spread function
    virtual void    SetChargeSpread(Float_t p1, Float_t p2)   =0;
    virtual Float_t ChargeSpreadX()                           =0;
    virtual Float_t ChargeSpreadY()                           =0;
    // Adc-count saturation value
    virtual void    SetMaxAdc(Float_t p1)                     =0;
    virtual Float_t MaxAdc()                                  =0;
    // anode cathode Pitch
    virtual void    SetPitch(Float_t)                         =0;
    virtual Float_t Pitch()                                   =0;
    // alpha feedback
    virtual void    SetAlphaFeedback(Float_t)                 =0;
    virtual Float_t AlphaFeedback()                           =0;
    // ionisation enrgy
    virtual void    SetEIonisation(Float_t)                   =0;
    virtual Float_t EIonisation()                             =0;
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t eloss)                       =0;
    virtual Float_t IntPH()                                    =0;
    // Charge disintegration
    virtual Float_t IntXY(AliRICHSegmentation *)                 =0;
    virtual Int_t   FeedBackPhotons(Float_t *source, Float_t qtot) =0;
    //
    // Mathieson parameters
    virtual void   SetSqrtKx3(Float_t p1)                        =0;
    virtual void   SetKx2(Float_t p1)                            =0;
    virtual void   SetKx4(Float_t p1)                            =0;
    virtual void   SetSqrtKy3(Float_t p1)                        =0;
    virtual void   SetKy2(Float_t p1)                            =0;
    virtual void   SetKy4(Float_t p1)                            =0;
    ClassDef(AliRICHResponse,1)
};

//----------------------------------------------
//
// Chamber geometry virtual base class
//
class AliRICHGeometry :
public TObject {
 public:
    //
    // Configuration methods
    //
   // Radiator Thickness 
    virtual void   SetGapThickness(Float_t t)      =0;
    // Proximity Gap Thickness
    virtual void   SetProximityGapThickness(Float_t t)           =0;
    // Quartz Length
    virtual void   SetQuartzLength(Float_t t)           =0;
    // Quartz Width
    virtual void   SetQuartzWidth(Float_t t)            =0;
    // Quartz Thickness
    virtual void   SetQuartzThickness(Float_t t)        =0;
    // Freon Length
    virtual void   SetOuterFreonLength(Float_t t)           =0;
    // Freon Width
    virtual void   SetOuterFreonWidth(Float_t t)            =0;
    // Freon Length
    virtual void   SetInnerFreonLength(Float_t t)           =0;
    // Freon Width
    virtual void   SetInnerFreonWidth(Float_t t)            =0;
    // Quartz Thickness
    virtual void   SetFreonThickness(Float_t t)        =0;
    // Distance between radiator and pads
    virtual void   SetRadiatorToPads(Float_t)                =0;

    // Radiator thickness
    virtual Float_t  GetGapThickness()             =0;
    // Proximity Gap thickness
    virtual Float_t  GetProximityGapThickness()                  =0;
    // Quartz Length
    virtual Float_t  GetQuartzLength()                  =0;
    // Quartz Width
    virtual Float_t  GetQuartzWidth()                   =0;
    // Quartz Thickness
    virtual Float_t  GetQuartzThickness()               =0;
    // Freon Length
    virtual Float_t  GetOuterFreonLength()              =0;
    // Freon Width
    virtual Float_t  GetOuterFreonWidth()               =0;
    // Freon Length
    virtual Float_t  GetInnerFreonLength()              =0;
    // Freon Width
    virtual Float_t  GetInnerFreonWidth()               =0;
    // Freon Thickness
    virtual Float_t  GetFreonThickness()                =0;
    // Get distance between radiator and pads
    virtual Float_t  GetRadiatorToPads()                =0;

    ClassDef(AliRICHGeometry,1)
};

#endif
