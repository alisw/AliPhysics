#ifndef ALIRICHSEGMENTATIONV0_H
#define ALIRICHSEGMENTATIONV0_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliSegmentation.h"

class AliRICHSegmentationV0 :public AliSegmentation 
{
 public:
  AliRICHSegmentationV0();
    virtual ~AliRICHSegmentationV0(){}
    virtual  void    SetPadSize(Float_t p1, Float_t p2);
    // Anod Pitch
    virtual  void    SetDAnod(Float_t D) {fWireD = D;};
        
    virtual Float_t GetAnod(Float_t xhit) const;
    virtual void    GetPadI(Float_t x, Float_t y , Int_t &ix, Int_t &iy);
    virtual void    GetPadI(Float_t x, Float_t y , Float_t /*z*/, Int_t &ix, Int_t &iy) {GetPadI(x, y, ix, iy);}
    virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y);
    virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z) 
	{z=0; GetPadC(ix, iy, x , y);}
    virtual void Init(Int_t i);
    virtual Float_t Dpx() const {return fDpx;}
    virtual Float_t Dpy() const {return fDpy;}
    virtual Float_t Dpx(Int_t) const {return fDpx;}
    virtual Float_t Dpy(Int_t) const {return fDpy;}
    virtual Int_t   Npx() const {return fNpx;}
    virtual Int_t   Npy() const {return fNpy;}
    virtual Float_t   DeadZone() const {return fDeadZone;}
    virtual Float_t GetPadPlaneWidth()  const {return fPadPlane_Width;}
    virtual Float_t GetPadPlaneLength() const {return fPadPlane_Length;}
    virtual void     SetPad(Int_t ix, Int_t iy);
    virtual void     SetHit(Float_t xhit , Float_t yhit);
    virtual void     SetHit(Float_t xhit, Float_t yhit, Float_t /*zhit*/){SetHit(xhit, yhit);}
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t /*zhit*/, Float_t dx, Float_t dy){FirstPad(xhit, yhit, dx, dy);}
    virtual void  NextPad();
    virtual Int_t MorePads();
    virtual Float_t Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t *dummy);
    virtual void GetNParallelAndOffset(Int_t /*iX*/, Int_t /*iY*/,Int_t *Nparallel, Int_t *Offset) {*Nparallel=1;*Offset=0;}
    virtual void Neighbours	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
    virtual Int_t  Ix() {return fIx;}
    virtual Int_t  Iy() {return fIy;}
    virtual Int_t  ISector() {return 1;}
    virtual Int_t  Sector(Int_t,Int_t) {return 1;}
    virtual Int_t  Sector(Float_t,Float_t) {return 1;}
    virtual Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
    virtual void  SigGenInit(Float_t x, Float_t y, Float_t z);
    virtual void IntegrationLimits(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    virtual void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y) const;
    virtual void Draw(const char* = "") const; 
    virtual void SetCorrFunc(Int_t /*dum*/, TF1* func) {fCorr=func;}
    virtual TF1* CorrFunc(Int_t) const {return fCorr;} 
    ClassDef(AliRICHSegmentationV0,1)
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
    // Current pad during integration (cursor for disintegration)
    Int_t fIx;  // pad coord. x 
    Int_t fIy;  // pad coord. y 
    Float_t fX; // x
    Float_t fY; // y
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
};
#endif


