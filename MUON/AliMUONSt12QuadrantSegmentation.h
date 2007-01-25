#ifndef ALI_MUON_ST12_QUADRANT_SEGMENTATION_H
#define ALI_MUON_ST12_QUADRANT_SEGMENTATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONSt12QuadrantSegmentation
/// \brief Segmentation for MUON quadrants of stations 1 and 2 using 
/// the mapping package
///
/// \author Ivana Hrivnacova, IPN Orsay

#include "AliMpStationType.h"
#include "AliMpPlaneType.h"

#include "AliMUONVGeometryDESegmentation.h"

class TObjArray;

class AliMpSector;
class AliMpSectorSegmentation;
class AliMpVPadIterator;
class AliMpPad;

class AliMUONChamber;

class AliMUONSt12QuadrantSegmentation : public AliMUONVGeometryDESegmentation 
{
  public:
    AliMUONSt12QuadrantSegmentation(AliMpVSegmentation* segmentation,
                                    AliMp::StationType stationType, 
                                    AliMp::PlaneType planeType);
    AliMUONSt12QuadrantSegmentation();
    
    virtual ~AliMUONSt12QuadrantSegmentation();
    
    //    
    // methods derived from base class    
    // 

    // Set Chamber Segmentation Parameters
    //
    virtual void SetPadSize(Float_t p1, Float_t p2);
                       // Pad size Dx*Dy 
    virtual void SetDAnod(Float_t D);
                       // Anode Pitch

    // Check if pad exists
    //
    virtual Bool_t  HasPad(Float_t x, Float_t y, Float_t z); 
                       // Returns true if a pad exists in the given position
    virtual Bool_t  HasPad(Int_t ix, Int_t iy);
                       // Returns true if a pad with given indices exists

    // Quadrant type
    //
    virtual AliMUONGeometryDirection  GetDirection();
                       // Returns the direction with a constant pad size
    // Access to mapping
    virtual const AliMpVSegmentation* GetMpSegmentation() const; 		       

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
    virtual void Draw(const char *opt = "");
                   // Draw the segmentation zones

    // Function for systematic corrections
    //
    virtual void SetCorrFunc(Int_t isec,  TF1* func);
                   // Set the correction function
    virtual TF1* CorrFunc(Int_t isec)  const;
                   // Get the correction Function
 
  protected:
    AliMUONSt12QuadrantSegmentation(const AliMUONSt12QuadrantSegmentation& rhs);
  
    // operators
    AliMUONSt12QuadrantSegmentation& operator=(const AliMUONSt12QuadrantSegmentation & rhs);

  private:
    // methods
    void UpdateCurrentPadValues(const AliMpPad& pad);
  
    // data members

    // From mapping
    //
    AliMp::StationType       fStationType;       ///< Station type
    AliMp::PlaneType         fPlaneType;         ///< Plane type
    const AliMpSector*       fSector;            ///< Sector (from mapping)
    AliMpSectorSegmentation* fSectorSegmentation;///< Sector segmentation (from mapping)
    AliMpVPadIterator*       fSectorIterator;    //!< Iterator over pads

    // Wire pitch
    //
    Float_t         fWireD;  ///< \ brief Wire pitch
                             ///< (smaller distance between anode wires)
    
    // Reference to mother chamber
    //
    AliMUONChamber* fChamber; //!< Reference to mother chamber
    Int_t           fId;      ///< Identifier
    Float_t         fRmin;    ///< Inner radius
    Float_t         fRmax;    ///< Outer radius
    Float_t         fZ;       ///< Z-position of chamber

    // Current pad during integration (cursor for disintegration)
    //
    Int_t   fIx;     //!< Pad coord.  x 
    Int_t   fIy;     //!< Pad coord.  y 
    Float_t fX;      //!< Real coord. x
    Float_t fY;      //!< Real ccord. y
    Int_t   fZone;   //!< Current zone (sector in AliSegmentation naming)
    
    // Current pad and wire during tracking (cursor at hit centre)
    //
    Float_t fXhit;  //!< X-position of hit
    Float_t fYhit;  //!< Y-position of hit

    // Reference point to define signal generation condition
    //
    Int_t   fIxt;   //!< Pad coord. x
    Int_t   fIyt;   //!< Pad coord. y
    Int_t   fIwt;   //!< Wire number
    Float_t fXt;    //!< X
    Float_t fYt;    //!< Y

    TObjArray* fCorrA; //!< Array of correction functions

  ClassDef(AliMUONSt12QuadrantSegmentation,2) // Station1 segmentation
};

#endif //ALI_MUON_ST12_QUADRANT_SEGMENTATION_H








