#ifndef ALIMUONSEGMENTATIONSLAT_H
#define ALIMUONSEGMENTATIONSLAT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation classes for slat modules          //
//  to be used with AluMUONSegmentationSlat        //
/////////////////////////////////////////////////////

#include "AliSegmentation.h"

class TArrayI;
class TObjArray;
class AliMUONSegmentationSlatModule;
class AliMUONChamber;


class AliMUONSegmentationSlat :
public AliSegmentation {
 public:
    AliMUONSegmentationSlat();
    virtual ~AliMUONSegmentationSlat(){}
    //    
    // Set Chamber Segmentation Parameters
    //
    // Pad size Dx*Dy 
    virtual  void    SetPadSize(Float_t p1, Float_t p2);
    // Anod Pitch
    virtual  void    SetDAnod(Float_t D) {fWireD = D;};
    
    // Anod wire coordinate closest to xhit
    virtual Float_t GetAnod(Float_t xhit) const;
    
    void SetPadDivision(Int_t ndiv[4]);
    
    // Transform from pad to real coordinates
    virtual void    GetPadI(Float_t x, Float_t y , Float_t z, Int_t &ix, Int_t &iy);
    // Transform from real to pad coordinates
    virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z);
     // Initialisation
    virtual void Init(Int_t chamber);
    //
    // Get member data
    //
    //
    // Pad size in x
    virtual Float_t Dpx() const {return fDpx;}
    
    // Pad size in y 
    virtual Float_t Dpy() const {return fDpy;}
    // Pad size in x by Sector
    virtual Float_t Dpx(Int_t isec) const;
    // Pad size in y by Sector
    virtual Float_t Dpy(Int_t isec) const;
    // Maximum number of Pads in x
    virtual Int_t    Npx() const {return fNpx;}
    // Maximum number of Pads in y
    virtual Int_t    Npy() const {return fNpy;}
    //
    virtual void     SetPad(Int_t ix,Int_t iy);
    
    virtual void     SetHit(Float_t xhit, Float_t yhit, Float_t zhit);
    //
    // Iterate over pads
    // Initialiser
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t zhit, Float_t dx, Float_t dy);
    // Stepper
    virtual void  NextPad();
    // Condition
    virtual Int_t MorePads();
    // Get next neighbours 
    virtual void Neighbours
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
    virtual Float_t Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t *dummy) {return 0.;}
    virtual void GetNParallelAndOffset(Int_t iX, Int_t iY,
				       Int_t *Nparallel, Int_t *Offset) {*Nparallel=1;*Offset=0;}
    //
    // Current Pad during Integration
    // x-coordinate
    virtual Int_t  Ix();
    // y-coordinate
    virtual Int_t  Iy();
    // current sector
    virtual Int_t ISector();
    // calculate sector from pad coordinates
    virtual Int_t Sector(Int_t ix, Int_t iy);
    //
    // Signal Generation Condition during Stepping
    virtual Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
    // Initialise signal generation at coord (x,y,z)
    virtual void  SigGenInit(Float_t x, Float_t y, Float_t z);

    //
    // Integration
    // Current integration limits
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    //
    // Class specific methods
    virtual void SetNSlats(Int_t nslats) {fNSlats = nslats;}
    virtual void SetShift(Float_t shift) {fShift = shift;}
    virtual void SetNPCBperSector(Int_t *npcb);
    virtual void SetSlatXPositions(Float_t *xpos);
    virtual AliMUONSegmentationSlatModule* Slat(Int_t index) const;
    
// Not used
    // Test points for auto calibration
    virtual void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y)  const {;}
    // Draw the segmentation zones
    virtual void Draw(const char *opt = "") const {;}
    
    // Function for systematic corrections
    // Set the correction function
    virtual void SetCorrFunc(Int_t, TF1*) {;}
    
    // Get the correction Function
    virtual TF1* CorrFunc(Int_t) const {return NULL;}
 protected:
    
    virtual void GlobalToLocal(
	Float_t x, Float_t y, Float_t z,  Int_t &islat, Float_t &xlocal, Float_t &ylocal);
    virtual void GlobalToLocal(
	Int_t ix, Int_t iy, Int_t &islat, Int_t &ixlocal, Int_t &iylocal);

    virtual void LocalToGlobal(
	 Int_t islat, Float_t xlocal, Float_t ylocal, Float_t &x, Float_t  &y, Float_t &z);
    virtual void LocalToGlobal(
	 Int_t islat, Int_t ixlocal, Int_t iylocal, Int_t &ix, Int_t &iy);
    virtual void SetSymmetry(Int_t   ix,   Int_t iy);
    virtual void SetSymmetry(Float_t  x, Float_t  y);
   // Factory method for associated slat module class  
    virtual AliMUONSegmentationSlatModule* CreateSlatModule();
    
 protected:

    AliMUONChamber*      fChamber;               // Parent Chamber
    
    //
    //  Geometry
    //
    Float_t    fWireD;                            // Wire Pitch
    Int_t      fNSlats;                           // Number of slats
    Int_t      fPcb[10][4];                       // PcbSegmentation
    Float_t    fXPosition[10];                    // x-position of slats
    Float_t    fYPosition[10];                    // y-position of slats
    Float_t    fSlatX[10];                        // Slat x-dimension
    Float_t    fSlatY;                            // Slat y-dimension
    Float_t    fDpx;                              // Pad size x
    Float_t    fDpy;                              // Pad size y
    Int_t      fNpx;                              // maximum number of pads in x
    Int_t      fNpy;                              // maximum number of pads in y
    Int_t      fSym[2];                           // signs for symmetry trafo
    Float_t    fShift;                            // Half overlap of pad planes
    Float_t    fDz;                               // Half distance between slat planes
    
    TArrayI*    fNDiv;                             // Pad size division
    // Slats
    TObjArray*  fSlats;                           // Array of Slats
    // Proxy data
    AliMUONSegmentationSlatModule* fCurrentSlat;  // Pointer to current slat
    Int_t       fSlatIndex;                       // Current slat index
    ClassDef(AliMUONSegmentationSlat,1)           // Segmentation for Muon Chamber built from Slat Modules
};
	

#endif









