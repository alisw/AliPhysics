#ifndef RICHSegResV0_H
#define RICHSegResV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliRICH.h"
#include "AliRICHv0.h"

class AliRICHsegmentationV0 :
public AliRICHsegmentation {
 public:
    AliRICHsegmentationV0(){}
    virtual ~AliRICHsegmentationV0(){}
    //    
    // Set Chamber Segmentation Parameters
    virtual  void    SetPADSIZ(Float_t p1, Float_t p2);
    virtual  void    SetDAnod(Float_t D) {fWireD = D;};
    //
    // Transform from pad (wire) to real coordinates and vice versa  
    virtual Float_t GetAnod(Float_t xhit);
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy);
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y );
    //
    // Initialisation
    virtual void Init(AliRICHchamber*);
    //
    // Get member data
    virtual Float_t Dpx(){return fDpx;}
    virtual Float_t Dpy(){return fDpy;}
    virtual Int_t   Npx(){return fNpx;}
    virtual Int_t   Npy(){return fNpy;}
    //
    // Iterate over pads
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);
    virtual void  NextPad();
    virtual Int_t MorePads();
    // Get next neighbours 
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
    // Provisory RecCluster coordinates reconstructor
    virtual void FitXY(AliRICHRecCluster* Cluster,TClonesArray* RICHdigits);
    //
    // Current Pad during Integration
    virtual Int_t  Ix(){return fix;}
    virtual Int_t  Iy(){return fiy;}    
    virtual Int_t  ISector(){return 1;}
    //
    // Signal Generation Condition during Stepping
    virtual Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
    virtual void  SigGenInit(Float_t x, Float_t y, Float_t z);
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    //
    // Identification
    virtual char* YourName(){return fName;}
    ClassDef(AliRICHsegmentationV0,1)
	protected:
    //
    // Implementation of the segmentation data
    // Version 0 models rectangular pads with the same dimensions all
    // over the cathode plane
    //
    //  geometry
    //
    Float_t    fDpx;           // x pad width per sector  
    Float_t    fDpy;           // y pad base width
    Int_t      fNpx;
    Int_t      fNpy;           // Number of pads in y
    Float_t    fWireD;         // wire pitch
    
    // Chamber region consideres during disintegration (lower left and upper right corner)
    //
    Int_t fixmin;
    Int_t fixmax;
    Int_t fiymin;
    Int_t fiymax;
    //
    // Current pad during integration (cursor for disintegration)
    Int_t fix;
    Int_t fiy;
    Float_t fx;
    Float_t fy;
    //
    // Current pad and wire during tracking (cursor at hit centre)
    Int_t fixt;
    Int_t fiyt;
    Int_t fiwt;
    Float_t fxt;
    Float_t fyt;
    //
    char    *fName;       //! Version Identifier
};

class AliRICHresponseV0 : //Mathieson response
public AliRICHresponse {
 public:
    AliRICHresponseV0(){}
    virtual ~AliRICHresponseV0(){}
    
    
    
    //
    // Configuration methods
    // 
    virtual void   SetRSIGM(Float_t p1) {fNsigma=p1;} 
    virtual void   SetMUCHSP(Float_t p1) {fChslope=p1;}
    virtual void   SetMUSIGM(Float_t p1, Float_t p2) {fChwX=p1; fChwY=p2;}
    virtual void   SetMAXADC(Float_t p1) {fadc_satm=p1;}
    // Mathieson parameters
    virtual void   SetSqrtKx3(Float_t p1) {fSqrtKx3=p1;};
    virtual void   SetKx2(Float_t p1) {fKx2=p1;};
    virtual void   SetKx4(Float_t p1) {fKx4=p1;};
    virtual void   SetSqrtKy3(Float_t p1) {fSqrtKy3=p1;};
    virtual void   SetKy2(Float_t p1) {fKy2=p1;};
    virtual void   SetKy4(Float_t p1) {fKy4=p1;};
    virtual void   SetPitch(Float_t p1) {fPitch=p1;};
    
    //
    // Get member data
    virtual Float_t Chslope() {return fChslope;}
    virtual Float_t ChwX() {return fChwX;}    
    virtual Float_t ChwY() {return fChwY;}        
    virtual Float_t Nsigma() {return fNsigma;}    
    virtual Float_t adc_satm() {return fadc_satm;}
    //  
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t eloss=0);
    virtual Int_t FeedBackPhotons(Float_t *source, Float_t qtot);
    
    // Charge disintegration
    virtual Float_t IntXY(AliRICHsegmentation * segmentation);
    // Identification
    //
    virtual char* YourName() {return fName;}
    
    ClassDef(AliRICHresponseV0,1)
	protected:
    Float_t fChslope;         // Slope of the charge distribution
    Float_t fChwX;            // Width of the charge distribution in x
    Float_t fChwY;            // Width of the charge distribution in y
    Float_t fNsigma;          // Number of sigma's used for charge distribution
    Float_t fadc_satm;        // Maximum ADC channel
    Float_t fSqrtKx3;         // Mathieson parameters for x
    Float_t fKx2;
    Float_t fKx4;
    Float_t fSqrtKy3;         // Mathieson parameters for y
    Float_t fKy2;
    Float_t fKy4;
    Float_t fPitch;           //anode-cathode pitch
    char    *fName;           //! Version Identifier
};

#endif
