#ifndef RICHSegResV0_H
#define RICHSegResV0_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include "AliRICH.h"
#include "AliRICHv0.h"
#include "AliRICHSegRes.h"

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
    virtual void Init(AliRICHChamber*);
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
    virtual void     SetPad(Int_t, Int_t);
    // set hit position
    virtual void     SetHit(Float_t, Float_t);
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
    Float_t fxhit;
    Float_t fyhit;
    // Reference point to define signal generation condition
    Int_t fixt;     // pad coord. x
    Int_t fiyt;     // pad coord. y
    Int_t fiwt;     // wire number
    Float_t fxt;    // x
    Float_t fyt;    // y
    TF1*    fCorr;  // correction function
};

class AliRICHResponseV0 : //Mathieson response
public AliRICHResponse {
 public:
    AliRICHResponseV0(){}
    virtual ~AliRICHResponseV0(){}
    //
    // Configuration methods
    // 
    // Number of sigmas over which cluster didintegration is performed
    virtual void    SetSigmaIntegration(Float_t p1) {fSigmaIntegration=p1;}
    virtual Float_t SigmaIntegration() {return fSigmaIntegration;}    
    // charge slope in ADC/e
    virtual void    SetChargeSlope(Float_t p1) {fChargeSlope=p1;}
    virtual Float_t ChargeSlope()      {return fChargeSlope;}
    // sigma of the charge spread function
    virtual void    SetChargeSpread(Float_t p1, Float_t p2)
	{fChargeSpreadX=p1; fChargeSpreadY=p2;}
    virtual Float_t ChargeSpreadX()    {return fChargeSpreadX;}    
    virtual Float_t ChargeSpreadY()    {return fChargeSpreadY;}        
    // Adc-count saturation value
    virtual void    SetMaxAdc(Float_t p1) {fMaxAdc=p1;}
    virtual Float_t MaxAdc()           {return fMaxAdc;}
    // anode cathode Pitch
    virtual Float_t Pitch()            {return fPitch;}
    virtual void    SetPitch(Float_t p1) {fPitch=p1;};
    // alpha feedback
    virtual void    SetAlphaFeedback(Float_t alpha) {fAlphaFeedback=alpha;}
    virtual Float_t AlphaFeedback()  {return fAlphaFeedback;}
    // ionisation enrgy
    virtual void    SetEIonisation(Float_t e) {fEIonisation=e;}
    virtual Float_t EIonisation() {return fEIonisation;}                            
    // Mathieson parameters
    virtual void   SetSqrtKx3(Float_t p1) {fSqrtKx3=p1;};
    virtual void   SetKx2(Float_t p1) {fKx2=p1;};
    virtual void   SetKx4(Float_t p1) {fKx4=p1;};
    virtual void   SetSqrtKy3(Float_t p1) {fSqrtKy3=p1;};
    virtual void   SetKy2(Float_t p1) {fKy2=p1;};
    virtual void   SetKy4(Float_t p1) {fKy4=p1;};
    //  
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t eloss);
    virtual Float_t IntPH();
    // Charge disintegration
    virtual Float_t IntXY(AliRICHSegmentation * segmentation);
    virtual Int_t   FeedBackPhotons(Float_t *source, Float_t qtot);
	protected:
    Float_t fChargeSlope;              // Slope of the charge distribution
    Float_t fChargeSpreadX;            // Width of the charge distribution in x
    Float_t fChargeSpreadY;            // Width of the charge distribution in y
    Float_t fSigmaIntegration;         // Number of sigma's used for charge distribution
    Float_t fAlphaFeedback;            // Feedback photons coefficient
    Float_t fEIonisation;              // Mean ionisation energy
    Float_t fMaxAdc;                   // Maximum ADC channel
    Float_t fSqrtKx3;         // Mathieson parameters for x
    Float_t fKx2;
    Float_t fKx4;
    Float_t fSqrtKy3;         // Mathieson parameters for y
    Float_t fKy2;
    Float_t fKy4;
    Float_t fPitch;           //anode-cathode pitch
    ClassDef(AliRICHResponseV0,1)
};

class AliRICHGeometryV0 : //Chamber geometry
public AliRICHGeometry {
 public:
    // Radiator Thickness
    virtual void   SetGapThickness(Float_t thickness)  {fGapThickness=thickness;} 
    // Proximity Gap Thickness
    virtual void   SetProximityGapThickness(Float_t thickness)       {fProximityGapThickness=thickness;}
    // Quartz Length
    virtual void   SetQuartzLength(Float_t length)          {fQuartzLength=length;}
    // Quartz Width
    virtual void   SetQuartzWidth(Float_t width)            {fQuartzWidth=width;}
    // Quartz Thickness
    virtual void   SetQuartzThickness(Float_t thickness)    {fQuartzThickness=thickness;}
    // Freon Length
    virtual void   SetOuterFreonLength(Float_t length)          {fOuterFreonLength=length;}
    // Freon Width
    virtual void   SetOuterFreonWidth(Float_t width)            {fOuterFreonWidth=width;}
    // Freon Length
    virtual void   SetInnerFreonLength(Float_t length)          {fInnerFreonLength=length;}
    // Freon Width
    virtual void   SetInnerFreonWidth(Float_t width)            {fInnerFreonWidth=width;}
    // Freon Thickness
    virtual void   SetFreonThickness(Float_t thickness)    {fFreonThickness=thickness;}
    // Freon Thickness
    virtual void   SetRadiatorToPads(Float_t distance)   {fRadiatorToPads=distance;}

    // Radiator thickness
    virtual Float_t  GetGapThickness(){return fGapThickness;}
    // Proximity Gap thickness
    virtual Float_t  GetProximityGapThickness(){return fProximityGapThickness;}
    // Quartz Length
    virtual Float_t  GetQuartzLength(){return fQuartzLength;}
    // Quartz Width
    virtual Float_t  GetQuartzWidth(){return fQuartzWidth;}
    // Quartz Thickness
    virtual Float_t  GetQuartzThickness(){return fQuartzThickness;}
    // Freon Length
    virtual Float_t  GetOuterFreonLength(){return fOuterFreonLength;}
    // Freon Width
    virtual Float_t  GetOuterFreonWidth(){return fOuterFreonWidth;}
    // Freon Length
    virtual Float_t  GetInnerFreonLength(){return fInnerFreonLength;}
    // Freon Width
    virtual Float_t  GetInnerFreonWidth(){return fInnerFreonWidth;}
    // Freon Thickness
    virtual Float_t  GetFreonThickness(){return fFreonThickness;}
    // Get distance between radiator and pads
    virtual Float_t  GetRadiatorToPads(){return fRadiatorToPads;}   
    

    //Self-explainable:

    Float_t fGapThickness;
    Float_t fProximityGapThickness;
    Float_t fQuartzLength;
    Float_t fQuartzWidth;
    Float_t fQuartzThickness;
    Float_t fOuterFreonLength;
    Float_t fOuterFreonWidth;
    Float_t fInnerFreonLength;
    Float_t fInnerFreonWidth;
    Float_t fFreonThickness;
    
    Float_t fRadiatorToPads;      // Distance from radiator to pads

    ClassDef(AliRICHGeometryV0,1)
};

#endif
