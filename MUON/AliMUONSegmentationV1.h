#ifndef ALIMUONSEGMENTATIONV1_H
#define ALIMUONSEGMENTATIONV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliSegmentation.h"

const Int_t kNzone = 3;                // Specific for chamber with equal pads
const Int_t kNzonem1 = 2;              // kNzone - 1
const Int_t kNzoneCUT = 30;            

class AliMUONSegmentationV1 :
public AliSegmentation {
 public:
    AliMUONSegmentationV1();
    AliMUONSegmentationV1(const AliMUONSegmentationV1 & segmentation);
    virtual ~AliMUONSegmentationV1(){}
    //    
    // Set Chamber Segmentation Parameters
    // Set Number of zones
    void SetNzone(Int_t N) {fNzone = N;};
    // Pad size Dx*Dy 
    virtual  void    SetPadSize(Float_t p1, Float_t p2);
    // Set Sense wire offset
    void SetSensOffset(Float_t Offset) {fSensOffset = Offset;};
    // Anod Pitch
    void SetDAnod(Float_t D) {fDAnod = D;};
    // max x and y for the zone in number of pads units 
    //(WARNING : first pad is labelled 0 !!) 
    virtual void AddCut(Int_t Zone, Int_t nX, Int_t nY);
    // Apply default cut
    virtual void DefaultCut(void);
    //
    // Initialisation
    virtual void Init(Int_t id);
    //
    // Get member data
    //
    // Pad size in x
    virtual Float_t Dpx(){return fDpx;}
     // Pad size in y 
    virtual Float_t Dpy(){return fDpy;}
    // Pad size in x by Sector
    virtual Float_t Dpx(Int_t i){return fDpx;}
    // Pad size in y by Sector 
    virtual Float_t Dpy(Int_t i){return fDpy;}
    // Maximum number of Pads in x
    virtual Int_t   Npx(){return fNpx;}
    // Maximum number of Pads in y
    virtual Int_t   Npy(){return fNpy;}
    //
    // Get  the zone of segmentation
    virtual Int_t GetZone(Float_t X, Float_t Y);
    virtual Int_t GetZone(Int_t X, Int_t Y);
    //
    // Transform from pad (wire) to real coordinates and vice versa  
    virtual Int_t GetiAnod(Float_t xhit);
    // Anod wire coordinate closest to xhit
    virtual Float_t GetAnod(Float_t xhit);
    // Transform from pad to real coordinates
    virtual void    GetPadI(Float_t x ,Float_t y , Int_t &ix, Int_t &iy);
    virtual void    GetPadI(Float_t x, Float_t y , Float_t z, Int_t &ix, Int_t &iy)
	{GetPadI(x, y, ix, iy);}
    // Transform from real to pad coordinates
    virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y);
    virtual void    GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z)
	{z=0; GetPadC(ix, iy, x , y);}
    // Set pad position
    virtual void     SetPad(Int_t ix, Int_t iy);
    // Set hit position
    virtual void     SetHit(Float_t xhit, Float_t yhit);
    virtual void     SetHit(Float_t xhit, Float_t yhit, Float_t zhit)
	 {SetHit(xhit, yhit);}
    //
    // Iterate over pads
    // Set Pad coordinates
    virtual void SetPadCoord(Int_t iX, Int_t iY);
    // Initialiser
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t zhit, Float_t dx, Float_t dy)
	{FirstPad(xhit, yhit, dx, dy);}
    // Stepper
    virtual void  NextPad();
    // Condition
    virtual Int_t MorePads();
    // Get next neighbours 
    virtual void Neighbours // implementation Neighbours function
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t *Xlist, Int_t *Ylist);
    virtual void NeighboursDiag // with diagonal elements
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t *Xlist, Int_t *Ylist);
    virtual void NeighboursNonDiag // without diagonal elements
	(Int_t iX, Int_t iY, Int_t* Nlist, Int_t *Xlist, Int_t *Ylist);
    void CleanNeighbours(Int_t* Nlist, Int_t *Xlist, Int_t *Ylist);
    //
    // Current pad cursor during disintegration
    // x-coordinate
    virtual Int_t Ix(Int_t trueX, Int_t trueY);
    virtual Int_t Ix();
    // y-coordinate
    virtual Int_t Iy(){return fiy;}
    // current sector
    virtual Int_t ISector();
    // calculate sector from pad coordinates
    virtual Int_t Sector(Int_t ix, Int_t iy) {return 1;}
    // Position of pad in perellel read-out
    virtual Int_t IsParallel2(Int_t iX, Int_t iY);
    virtual Int_t IsParallel3(Int_t iX, Int_t iY);
    // Number of pads read in parallel
    virtual Int_t NParallel2(Int_t iX, Int_t iY);
    virtual Int_t NParallel3(Int_t iX, Int_t iY);
    //
    // Number of pads read in parallel and offset to add to x
    virtual void GetNParallelAndOffset(Int_t iX, Int_t iY,
    	Int_t *Nparallel, Int_t *Offset);
    // Minimum distance between 1 pad and a position
    virtual Float_t Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t *Offset);
    //
    // Signal Generation Condition during Stepping
    Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
    // Initialise signal generation at coord (x,y,z)
    void  SigGenInit(Float_t x, Float_t y, Float_t z);
    // Test points for auto calibration
    void  GiveTestPoints(Int_t &n, Float_t *x, Float_t *y);
    // Current integration limits 
    virtual void IntegrationLimits
	(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2);
    // Draw the segmentation zones
    virtual void Draw(){;}
    // Function for systematic corrections
    // Set the correction function
    virtual void SetCorrFunc(Int_t dum, TF1* func) {fCorr=func;}
    // Get the correction function
    virtual TF1* CorrFunc(Int_t) {return fCorr;}
    //
    AliMUONSegmentationV1& operator=(const AliMUONSegmentationV1& rhs);
    ClassDef(AliMUONSegmentationV1,1) // Implementation of the Lyon type chamber segmentation with parallel read-out
 protected:
    //
    // Implementation of the segmentation data
    // Version This models rectangular pads with the same dimensions all
    // over the cathode plane but let the possibilit for different design.
    //
    //  geometry
    Int_t fNzone; // Number of differents sensitive zones
    Float_t fDpx;         // X pad width
    Float_t fDpy;         // Y pad width
    Int_t   fNZoneCut[kNzonem1];    // Number of cuts for given zone 
    Int_t fZoneX[kNzonem1][kNzoneCUT]; // X descriptor of zone segmentations
    Int_t fZoneY[kNzonem1][kNzoneCUT]; // Y descriptor of zone segmentations
    Float_t frSensMax2; // square of maximum sensitive radius
    Float_t frSensMin2; // square of minimum sensitive radius
    Int_t   fNpx;         // Maximum number of pads along x
    Int_t   fNpy;         // Maximum number of pads along y
    Float_t fDAnod;       // Anod gap
    Float_t fSensOffset;  // Offset of sensitive zone with respect to quadrant (positive)
    
    // Chamber region consideres during disintegration (lower left and upper right corner)
    //
    Int_t fixmin; // lower left  x
    Int_t fixmax; // lower left  y
    Int_t fiymin; // upper right x
    Int_t fiymax; // upper right y 
    //
    // Current pad during integration (cursor for disintegration)
    Int_t fix;   // pad coord.  x
    Int_t fiy;   // pad coord.  y
    Float_t fx;  // real coord. x
    Float_t fy;  // real ccord. y
    //
    // Current pad and wire during tracking (cursor at hit centre)
    Int_t fixt;  // x-position of hit
    Int_t fiyt;  // y-position of hit
    // Reference point to define signal generation condition
    Int_t fiwt;     // wire number
    Float_t fxt;    // x
    Float_t fyt;    // y
    Float_t fxhit;  // x-position of hit
    Float_t fyhit;  // y-position of hit
    
    TF1* fCorr;     // correction function
};

#endif


