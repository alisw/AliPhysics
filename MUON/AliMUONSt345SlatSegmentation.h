#ifndef ALIMUONST345SLATSEGMENTATION_H
#define ALIMUONST345SLATSEGMENTATION_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: */

//*********************************************************
//  Segmentation classes for slat modules          
//  This class works with local coordinates
//  of the slats via the class AliMUONGeometrySegmentation
//  This class contains the size of the slats and the
//  and the differents PCB densities. 
//*********************************************************

#include  "AliMUONVGeometryDESegmentation.h"

class TArrayF;
class TArrayI;

class AliMUONSt345SlatSegmentation : public AliMUONVGeometryDESegmentation 
{
 public:
    AliMUONSt345SlatSegmentation();
    AliMUONSt345SlatSegmentation(Bool_t bending);
    virtual ~AliMUONSt345SlatSegmentation();
      
    virtual Float_t  Distance2AndOffset(Int_t iX, Int_t iY, Float_t X, Float_t Y, Int_t * dummy);  // Distance between 1 pad and a position
    virtual Float_t  Dpx() const {return fDpx;}  // Pad size in x   
    virtual Float_t  Dpy() const {return fDpy;}  // Pad size in y   
    virtual Float_t  Dpx(Int_t isec) const;       // Pad size in x by Sector
    virtual Float_t  Dpy(Int_t isec) const;       // Pad size in y by Sector
    virtual void     Draw(const char */*opt*/ = "") {}  // Not implemented
    virtual void     FirstPad(Float_t xhit, Float_t yhit, Float_t dx, Float_t dy);  // Initialisation for pad iteration
    virtual void     FirstPad(Float_t xhit, Float_t yhit, Float_t zhit, Float_t dx, Float_t dy);

    virtual Bool_t   HasPad(Float_t /*x*/, Float_t /*y*/, Float_t /*z*/) { return true; }
    virtual Bool_t   HasPad(Int_t /*ix*/, Int_t /*iy*/) { return true; }
    virtual AliMUONGeometryDirection  GetDirection() { return kDirUndefined; } 

    virtual Float_t  GetAnod(Float_t xhit) const;  // Anod wire coordinate closest to xhit
    virtual void     GetPadI(Float_t x ,Float_t y ,Int_t   &ix,Int_t &iy);  // Transform from pad to real coordinates
    virtual void     GetPadI(Float_t x, Float_t y , Float_t z, Int_t &ix, Int_t &iy);
    virtual void     GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y);
    virtual void     GetPadC(Int_t ix, Int_t iy, Float_t &x, Float_t &y, Float_t &z) {z=0; GetPadC(ix, iy, x , y);}
 
 
    virtual void     IntegrationLimits(Float_t& x1, Float_t& x2, Float_t& y1, Float_t& y2); //Current integration limits
    virtual Int_t    ISector()  {return fSector;} // Current Pad during Integration (current sector)
    virtual Int_t    Ix() {return fIx;} // x-coordinate
    virtual Int_t    Iy() {return fIy;} // y-coordinate
  
    virtual Int_t    MorePads();  // Condition
 
    virtual void     Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);  // Get next neighbours
    virtual void     NextPad(); // Stepper
    virtual Int_t    Npx() const {return fNpx;} // Maximum number of Pads in x
    virtual Int_t    Npy() const {return fNpy;} // Maximum number of Pads in y

    virtual void     SetDAnod(Float_t D) {fWireD = D;};  // Anod pitch
    virtual Int_t    Sector(Int_t ix, Int_t iy);         // Calculate sector from pad coordinates
    virtual void     SetHit(Float_t xhit, Float_t yhit); // Set hit position
    virtual void     SetHit(Float_t xhit, Float_t yhit, Float_t zhit);
    virtual void     SetId(Int_t id) {fId=id;}  // Setting detection element
    virtual void     SetPad(Int_t ix, Int_t iy);         // Set pad position
    virtual void     SetPadDivision(Int_t ndiv[4]);      // Set Slat Segmentation Parameters
    virtual void     SetPadSize(Float_t p1, Float_t p2); // Pad size Dx*Dy 
    virtual void     SetPcbBoards(Int_t n[4]);           // Set Segmentation Zones (PCB Boards)
 
    // The following function could be obsolet for this class, but they are pure virtual in AliSegmentation
    virtual void     GetNParallelAndOffset(Int_t /*iX*/, Int_t /*iY*/, Int_t */*Nparallel*/, Int_t */*Offset*/) {};
    virtual Int_t    SigGenCond(Float_t /*x*/, Float_t /*y*/, Float_t /*z*/){return 0;} ;  // Signal Generation Condition during Stepping
    virtual void     SigGenInit(Float_t /*x*/, Float_t /*y*/, Float_t /*z*/){};  // Initialise signal gneration at coord (x,y,z)
    virtual void     GiveTestPoints(Int_t &/*n*/, Float_t * /*x*/, Float_t */*y*/) const{};   // Test points for auto calibration
    virtual void     SetCorrFunc(Int_t /*dum*/, TF1* /*func*/){}; // Function for systematic corrections, Set the correction function
    virtual TF1*     CorrFunc(Int_t) const {return 0x0;} // Get the correction Function
    virtual Int_t    Sector(Float_t /*x*/, Float_t /*y*/) {return 1;}

    virtual void     Init(Int_t detectionElementId); // Initialisation
    // Current integration limits

 protected:

    AliMUONSt345SlatSegmentation(const AliMUONSt345SlatSegmentation& rhs);
    AliMUONSt345SlatSegmentation& operator=(const AliMUONSt345SlatSegmentation& rhs);
    
    //  Internal geometry of the slat 
    Bool_t      fBending;        // 0: Bending or 1:Non Bending segmentation
    Int_t       fId;             // Identifier of detection element
    Int_t       fNsec;           // Number of density sectors (should be 4, if not not warranty about the output
    TArrayI*    fNDiv;           // Densities (d1, d2, d3, d4). It should be (4, 4, 2, 1) which goes from beam to out-beam
    TArrayF*    fDpxD;           // x pad width per density sector
    TArrayF*    fDpyD;           // x pad width per density sector
    Float_t     fDpx;            // x pad base width  
    Float_t     fDpy;            // y pad base width
    Int_t       fNpx;            // Number of pads in x
    Int_t       fNpy;            // Number of pads in y
    Float_t     fWireD;          // wire pitch
    // 
    Int_t       fSector;         // Current density sector
    Float_t     fDxPCB;          // x-size of PCB board
    Float_t     fDyPCB;          // y-size of PCB board
    Int_t       fPcbBoards[4];   // number of PCB boards per density sector n1,n2,n3,n4 
    // n1 PcbBoard with density d1, n2 PcbBoards with density d2, etc ....
   
    // Segmentation map
    Int_t       fNpxS[10];       // Number of pads per sector in x
    Int_t       fNpyS[10];       // Number of pads per sector in y    
    Float_t     fCx[10];         // pad-sector contour x vs y      
    Float_t     fCy;             // y offset      

    // Current pad and wire during tracking (cursor at hit centre)
    Float_t     fXhit;  // ! x-position of hit
    Float_t     fYhit;  // ! y-position of hit

    // Current pad and wire during tracking (cursor at hit centre)
    Int_t       fIx;   // ! pad coord.  x 
    Int_t       fIy;   // ! pad coord.  y 
    Float_t     fX;    // ! real coord. x
    Float_t     fY;    // ! real ccord. y
    
    // Chamber region consideres during disintegration   
    Int_t       fIxmin; // ! lower left  x
    Int_t       fIxmax; // ! lower left  y
    Int_t       fIymin; // ! upper right x
    Int_t       fIymax; // ! upper right y 

    // Chamber region consideres during disintegration  (lower left and upper right corner)
    Float_t     fXmin;           // lower left  x
    Float_t     fXmax;           // lower left  y
    Float_t     fYmin;           // upper right x
    Float_t     fYmax;           // upper right y 

    ClassDef(AliMUONSt345SlatSegmentation,1) 
};
#endif






