#ifndef ALI_MUON_ST1_SEGMENTATION_H
#define ALI_MUON_ST1_SEGMENTATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONSt1Segmentation
// -----------------------------
// Segmentation for MUON station 1 using the external 
// mapping package

#include "AliMpPlaneType.h"

#include "AliSegmentation.h"

class TObjArray;

class AliMpPlane;
class AliMpPlaneSegmentation;
class AliMpVPadIterator;
class AliMpPad;

class AliMUONChamber;

class AliMUONSt1Segmentation : public AliSegmentation 
{
  public:
    AliMUONSt1Segmentation(const AliMpPlaneType planeType);
    AliMUONSt1Segmentation();
    
    virtual ~AliMUONSt1Segmentation();
    
    //    
    // methods derived from base class    
    // 

    // Set Chamber Segmentation Parameters
    //
    virtual void SetPadSize(Float_t p1, Float_t p2);
                       // Pad size Dx*Dy 
    virtual void SetDAnod(Float_t D);
                       // Anode Pitch

    // Transform from pad (wire) to real coordinates and vice versa
    //
    virtual Float_t GetAnod(Float_t xhit) const;
                       // Anode wire coordinate closest to xhit
    virtual void  GetPadI(Float_t x, Float_t y, Float_t  z, Int_t& ix, Int_t& iy);
    virtual void  GetPadI(Float_t x, Float_t y , Int_t &ix, Int_t &iy) ;
                       // Transform from pad to real coordinates
    virtual void  GetPadC(Int_t ix, Int_t iy, Float_t& x, Float_t& y, Float_t& z);
    virtual void  GetPadC(Int_t ix, Int_t iy, Float_t& x, Float_t& y);
                       // Transform from real to pad coordinates
                      

    // Initialisation
    //
    virtual void Init(Int_t chamber);
 
    // Get member data
    //
    virtual Float_t Dpx() const;
    virtual Float_t Dpy() const;
                      // Pad size in x, y 
    virtual Float_t Dpx(Int_t isector) const;
    virtual Float_t Dpy(Int_t isector) const;
                      // Pad size in x, y by Sector 
    virtual Int_t   Npx() const;
    virtual Int_t   Npy() const;
                      // Maximum number of Pads in y

    virtual void  SetPad(Int_t ix, Int_t iy);
                      // Set pad position
    virtual void  SetHit(Float_t xhit, Float_t yhit, Float_t zhit);
                      // Set hit position
    
    // Iterate over pads
    //
    virtual void  FirstPad(Float_t xhit, Float_t yhit, Float_t zhit, 
                           Float_t dx, Float_t dy);
    virtual void  NextPad();
    virtual Int_t MorePads();

    virtual Float_t Distance2AndOffset(Int_t iX, Int_t iY, 
                                       Float_t X, Float_t Y, Int_t* dummy) ;
                      // Distance between 1 pad and a position
    virtual void GetNParallelAndOffset(Int_t iX, Int_t iY,
				       Int_t* Nparallel, Int_t* Offset);
                      // Number of pads read in parallel and offset to add to x 
                      // (specific to LYON, but mandatory for display)
    virtual void Neighbours(Int_t iX, Int_t iY, 
                            Int_t* Nlist, Int_t Xlist[10], Int_t Ylist[10]);
                      // Get next neighbours 

    // Current values
    //
    virtual Int_t  Ix();
    virtual Int_t  Iy();
                     // Current pad cursor during disintegration
                     // x, y-coordinate
    virtual Int_t  ISector();
                    // current sector

    virtual Int_t  Sector(Int_t ix, Int_t iy);
    virtual Int_t  Sector(Float_t x, Float_t y);
                    // calculate sector from pad coordinates

    virtual void  IntegrationLimits(Float_t& x1, Float_t& x2,
                                    Float_t& y1, Float_t& y2);
                   // Current integration limits 

    // Signal Generation
    //
    virtual Int_t SigGenCond(Float_t x, Float_t y, Float_t z);
                    // Signal Generation Condition during Stepping
    virtual void  SigGenInit(Float_t x, Float_t y, Float_t z);
                    // Initialise signal generation at coord (x,y,z)
		    
    
    virtual void GiveTestPoints(Int_t& n, Float_t* x, Float_t* y) const;
                   // Test points for auto calibration
    virtual void Draw(const char *opt = "") const;
                   // Draw the segmentation zones

    // Function for systematic corrections
    //
    virtual void SetCorrFunc(Int_t isec,  TF1* func);
                   // Set the correction function
    virtual TF1* CorrFunc(Int_t isec)  const;
                   // Get the correction Function
 
  protected:
    AliMUONSt1Segmentation(const AliMUONSt1Segmentation& rhs);
  
    // operators
    AliMUONSt1Segmentation& operator=(const AliMUONSt1Segmentation & rhs);

  private:
    // methods
    void UpdateCurrentPadValues(const AliMpPad& pad);
  
    // constants
    static const Float_t  fgkWireD;     // default wire pitch
    static const Float_t  fgkLengthUnit;// conversion between length units
                                        // from mapping (mm) to AliRoot (cm)
  
    // data members

    // From mapping
    //
    AliMpPlane*             fPlane;            //  plane (from mapping)
    AliMpPlaneSegmentation* fPlaneSegmentation;//  plane segmantation (from mapping)
    AliMpVPadIterator*      fPlaneIterator;    // ! iterator over pads

    // Wire pitch
    //
    Float_t         fWireD;  // wire pitch
                             // (smaller distance between anode wires)
    
    // Reference to mother chamber
    //
    AliMUONChamber* fChamber; // ! Reference to mother chamber
    Int_t           fId;      // Identifier
    Float_t         fRmin;    // inner radius
    Float_t         fRmax;    // outer radius
    Float_t         fZ;       // z-position of chamber

    // Current pad during integration (cursor for disintegration)
    //
    Int_t   fIx;     // ! pad coord.  x 
    Int_t   fIy;     // ! pad coord.  y 
    Float_t fX;      // ! real coord. x
    Float_t fY;      // ! real ccord. y
    Int_t   fSector;    // ! Current sector
    
    // Current pad and wire during tracking (cursor at hit centre)
    //
    Float_t fXhit;  // ! x-position of hit
    Float_t fYhit;  // ! y-position of hit

    // Reference point to define signal generation condition
    //
    Int_t   fIxt;   // ! pad coord. x
    Int_t   fIyt;   // ! pad coord. y
    Int_t   fIwt;   // ! wire number
    Float_t fXt;    // ! x
    Float_t fYt;    // ! y

    TObjArray* fCorrA; // ! Array of correction functions

  ClassDef(AliMUONSt1Segmentation,1) // Station1 segmentation
};

#endif //ALI_MUON_ST1_SEGMENTATION_H








