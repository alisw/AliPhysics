#ifndef MUONSegRes_H
#define MUONSegRes_H
#include "TObject.h"
#include "TClonesArray.h"
#include "TF1.h"
class AliMUONchamber;

//----------------------------------------------
//
// Chamber segmentation virtual base class
//
class AliMUONsegmentation :
public TObject {
 public:
    // Set Chamber Segmentation Parameters
    //
    // Pad size Dx*Dy 
    virtual void    SetPADSIZ(Float_t p1, Float_t p2)  =0;
    // Anod Pitch
    virtual void    SetDAnod(Float_t D)                =0;
    // Transform from pad (wire) to real coordinates and vice versa
    //
    // Anod wire coordinate closest to xhit
    virtual Float_t GetAnod(Float_t xhit)              =0;
    // Transform from pad to real coordinates
    virtual void    GetPadIxy(Float_t x ,Float_t y ,Int_t   &ix,Int_t   &iy)=0;
    // Transform from real to pad coordinates
    virtual void    GetPadCxy(Int_t   ix,Int_t   iy,Float_t &x ,Float_t &y )=0;
    //
    // Initialisation
    virtual void Init(AliMUONchamber*)                 =0;
    //
    // Get member data
    //
    // Pad size in x
    virtual Float_t Dpx()                              =0;
    // Pad size in y 
    virtual Float_t Dpy()                              =0;
    // Pad size in x by Sector 
    virtual Float_t Dpx(Int_t)                         =0;
    // Pad size in y by Sector 
    virtual Float_t Dpy(Int_t)                         =0;
    // Max number of Pads in x
    virtual Int_t    Npx()                             =0;
    // max number of Pads in y
    virtual Int_t    Npy()                             =0;
    // set pad position
    virtual void     SetPad(Int_t, Int_t)              =0;
    // set hit position
    virtual void     SetHit(Float_t, Float_t)          =0;
    
    //
    // Iterate over pads
    // Initialiser
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy) =0;
    // Stepper
    virtual void  NextPad()                            =0;
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
    // Current pad cursor during disintegration
    // x-coordinate
    virtual Int_t  Ix()                                =0;
    // y-coordinate
    virtual Int_t  Iy()                                =0;
    // current sector
    virtual Int_t  ISector()                           =0;
    // calculate sector from pad coordinates
    virtual Int_t  Sector(Int_t ix, Int_t iy)          =0;
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
    virtual void Draw(Option_t *)                             = 0;
    // Function for systematic corrections
    virtual void SetCorrFunc(Int_t, TF1*)                         = 0;
    virtual TF1* CorrFunc(Int_t)                                  = 0;
	    
    ClassDef(AliMUONsegmentation,1) //Segmentation class for homogeneous segmentation
};
//----------------------------------------------
//
// Chamber response virtual base class
//
class AliMUONresponse :
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
    //  
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t eloss)                      =0;
    // Charge disintegration 
    virtual Float_t IntXY(AliMUONsegmentation *)              =0;

    ClassDef(AliMUONresponse,1) // Implementation of Mathieson CPC response 
};
#endif



