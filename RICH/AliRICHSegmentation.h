#ifndef AliRICHSegmentation_h
#define ALIRICHSegmentation_h


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "TObject.h"
class TF1;

class AliRICHSegmentation : public TObject 
{
public:
   AliRICHSegmentation();
   virtual ~AliRICHSegmentation(){}
// pad size getter&setter   
   void    SetPadSize(Float_t p1, Float_t p2) {fDpx=p1;fDpy=p2;Init();}
   Float_t Dpx()                         const{return fDpx;}        // Get pad width in cm
   Float_t Dpy()                         const{return fDpy;}        // Get pad heights in cm
   
   Int_t   Npx()                         const{return fNpx;}        // Get number of Pads in x direction
   Int_t   Npy()                         const{return fNpy;}        // Get number of Pads in y direction
   Float_t DeadZone()                    const{return fDeadZone;}   // Dead zone width in cm

   void    SetDAnod(Float_t D) {fWireD = D; Init();}
   Float_t GetAnod(Float_t xhit) const;
        
   Int_t  Ix()                          const{return fIx;}     // Current pad during integration: x-coordinate
   Int_t  Iy()                          const{return fIy;}     // Current pad during integration: y-coordinate
   Int_t  ISector()                     const{return fSector;} // Get current sector
   
   void  GetPadI(Float_t x, Float_t y , Int_t &ix, Int_t &iy);
//    void  GetPadI(Float_t x, Float_t y , Float_t z, Int_t &ix, Int_t &iy)  {GetPadI(x, y, ix, iy);}
    // Transform from real to pad coordinates
    void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y);
//    void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z) {z=0; GetPadC(ix, iy, x , y);}
    //
   void Init();    // Recalculation after changing the paremeters
    Float_t Dpx(Int_t) const {return fDpx;}
    Float_t Dpy(Int_t) const {return fDpy;}

    void SetCorrFunc(Int_t dum, TF1* func) {fCorr=func;}// Set function for systematic corrections
    TF1*    CorrFunc(Int_t)          const {return fCorr;}// Set function for systematic corrections
    
    Float_t GetPadPlaneWidth()  const {return fPadPlane_Width;}
    Float_t GetPadPlaneLength() const {return fPadPlane_Length;}
    
    void     SetPad(Int_t ix, Int_t iy);     // set pad position    
    void     SetHit(Float_t xhit , Float_t yhit);// set hit position
    
    
//    void     SetHit(Float_t xhit, Float_t yhit, Float_t zhit){SetHit(xhit, yhit);}
    //
    // Iterate over pads
    // Initialiser
    void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);
    void  FirstPad(Float_t xhit, Float_t yhit, Float_t zhit, Float_t dx, Float_t dy)
	{FirstPad(xhit, yhit, dx, dy);}
    // Stepper
    void  NextPad();
    // Condition
    Int_t MorePads();
    //
    // Distance between 1 pad and a position
    Float_t Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t *dummy);
    // Number of pads read in parallel and offset to add to x 
    // (specific to LYON, but mandatory for display)
    virtual void GetNParallelAndOffset(Int_t iX, Int_t iY,Int_t *Nparallel, Int_t *Offset) {*Nparallel=1;*Offset=0;}
     
    void Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
    
    
    Int_t  Sector(Float_t x,Float_t y);     // calculate sector from x-y coordinates

    // Signal Generation Condition during Stepping
    Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
    
    void  SigGenInit(Float_t x, Float_t y, Float_t z);// Initialise signal gneration at coord (x,y,z)
    void IntegrationLimits(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);             // Current integration limits
    void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y) const{n=1;x[0]=0.;y[0]=x[0];}     // Test points for auto calibration

    
     
protected:
   Float_t    fDpx;             // x pad width per sector  
   Float_t    fDpy;             // y pad base width
   Int_t      fNpx;             // Number of pads in x
   Int_t      fNpy;             // Number of pads in y
   Int_t      fSector;          // Current padplane
   Float_t    fWireD;           // wire pitch
    
   Float_t fDeadZone;               //width of deadzones beteween CsI padplanes
   Float_t fPadPlane_Width;         //width of CsI padplanes
   Float_t fPadPlane_Length;        //length of CsI padplanes

        
// Chamber region consideres during disintegration (lower left and upper right corner)
   Int_t fIxmin; // lower left  x
   Int_t fIxmax; // lower left  y
   Int_t fIymin; // upper right x
   Int_t fIymax; // upper right y 
    
   Int_t fIx;  // Current x pad during integration (cursor for disintegration)
   Int_t fIy;  // Current y pad during integration (cursor for disintegration)
   Float_t fX; // x
   Float_t fY; // y
    //
    // Current pad and wire during tracking (cursor at hit centre)
   Float_t fXhit;   //x position
   Float_t fYhit;   //y position
    // Reference point to define signal generation condition
   Int_t fIxt;     // pad coord. x
   Int_t fIyt;     // pad coord. y
   Int_t fIwt;     // wire number
   Float_t fXt;    // x
   Float_t fYt;    // y
   TF1*    fCorr;  // correction function
    
   ClassDef(AliRICHSegmentation,1)
};
#endif //AliRICHSegmentation_h
