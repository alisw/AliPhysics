/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Class AliMUONGeometrySegmentation
// ----------------------------------
// New class for the geometry segmentation 
// composed of the segmentations of detection elements.
// Applies transformations defined in geometry.
//
// Author:Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_GEOMETRY_SEGMENTATION_H
#define ALI_MUON_GEOMETRY_SEGMENTATION_H

#include <TObject.h>

#include "AliMUONGeometryDirection.h"

class TObjArray;
class TF1;

class AliMUONGeometryModule;
class AliMUONGeometryStore;
class AliMUONGeometryDetElement;
class AliMUONVGeometryDESegmentation;

class AliMUONGeometrySegmentation : public TObject
{
  public:
    AliMUONGeometrySegmentation(AliMUONGeometryModule* geometry);
    AliMUONGeometrySegmentation();
    virtual ~AliMUONGeometrySegmentation();

    // methods
    void Add(Int_t detElemId, 
             AliMUONVGeometryDESegmentation* segmentation); 
 
    // get methods
    AliMUONGeometryModule* GetGeometry() const;	      
    virtual AliMUONGeometryDirection GetDirection(Int_t detElemId) const;
                       // Direction with a constant pad size  (Direction or coordinate where the resolution is the best)
    
    //    
    // redefined methods from AliSegmentation interface
    // 

    // Set Chamber Segmentation Parameters
    //
    virtual void SetPadSize(Float_t p1, Float_t p2);
                       // Pad size Dx*Dy 
    virtual void SetDAnod(Float_t D);
                       // Anode Pitch

    // Transform from pad (wire) to real coordinates and vice versa
    //
    virtual Float_t GetAnod(Int_t detElemId, Float_t xlhit) const;
                       // Anode wire coordinate closest to xhit
    virtual Bool_t  GetPadI(Int_t detElemId,
                          Float_t xg, Float_t yg, Float_t  zg, 
                          Int_t& ix, Int_t& iy);
                       // Transform from pad to real coordinates
    virtual Bool_t  GetPadC(Int_t detElemId,
                          Int_t ix, Int_t iy,
                          Float_t& x, Float_t& y, Float_t& z);
                       // Transform from real to pad coordinates

    // Initialisation
    //
    virtual void Init(Int_t chamber);
 
    // Get member data
    //
    virtual Float_t Dpx(Int_t detElemId) const;
    virtual Float_t Dpy(Int_t detElemId) const ;
                      // Pad size in x, y 
    virtual Float_t Dpx(Int_t detElemId, Int_t isector) const;
    virtual Float_t Dpy(Int_t detElemId, Int_t isector) const;
                      // Pad size in x, y by Sector 
    virtual Int_t   Npx(Int_t detElemId) const;
    virtual Int_t   Npy(Int_t detElemId) const;
                      // Maximum number of Pads in y

    virtual void  SetPad(Int_t detElemId, Int_t ix, Int_t iy);
                      // Set pad position
    virtual void  SetHit(Int_t detElemId, Float_t xghit, Float_t yghit, Float_t zghit);
                      // Set hit position
    
    // Iterate over pads
    //    
    virtual void  FirstPad(Int_t detElemId, 
                           Float_t xghit, Float_t yghit, Float_t zghit, 
                           Float_t dx, Float_t dy);
    virtual void  NextPad(Int_t detElemId);
    virtual Int_t MorePads(Int_t detElemId);

    virtual Float_t Distance2AndOffset(Int_t detElemId,
                           Int_t ix, Int_t iy, 
                           Float_t xg, Float_t yg, Float_t zg, 
			   Int_t* dummy);
                      // Distance between 1 pad and a position
    virtual void GetNParallelAndOffset(Int_t detElemId,
                           Int_t ix, Int_t iy,
		           Int_t* nparallel, Int_t* offset);
                      // Number of pads read in parallel and offset to add to x 
                      // (specific to LYON, but mandatory for display)
    virtual void Neighbours(Int_t detElemId,
                            Int_t ix, Int_t iy,
                            Int_t* nlist, Int_t xlist[10], Int_t ylist[10]);
                      // Get next neighbours 

    // Current values
    //
    virtual Int_t  Ix();
    virtual Int_t  Iy();
    virtual Int_t  DetElemId();
                     // Current pad cursor during disintegration
                     // x, y-coordinate
    virtual Int_t  ISector();
                    // current sector

    virtual Int_t  Sector(Int_t detElemId,
                          Int_t ix, Int_t iy);
    virtual Int_t  Sector(Int_t detElemId,
                          Float_t xg, Float_t yg, Float_t zg);
                    // calculate sector from pad coordinates

    virtual void  IntegrationLimits(Int_t detElemId,
                          Float_t& x1, Float_t& x2,
                          Float_t& y1, Float_t& y2);
                   // Current integration limits 

    // Signal Generation
    //
    virtual Int_t SigGenCond(Int_t detElemId,
                          Float_t xg, Float_t yg, Float_t zg);
                    // Signal Generation Condition during Stepping
    virtual void  SigGenInit(Int_t detElemId,
                          Float_t xg, Float_t yg, Float_t zg);
                    // Initialise signal generation at coord (x,y,z)
		    
    
    virtual void GiveTestPoints(Int_t detElemId,
                          Int_t& n, Float_t* xg, Float_t* yg) const;
                   // Test points for auto calibration
    virtual void Draw(Int_t detElemId,
                          const char *opt = "") const;
                   // Draw the segmentation zones

    // Function for systematic corrections
    //
    virtual void SetCorrFunc(Int_t detElemId,
                          Int_t isec,  TF1* func);
                   // Set the correction function
    virtual TF1* CorrFunc(Int_t detElemId, Int_t isec) const;
                   // Get the correction Function
 
  protected:
    AliMUONGeometrySegmentation(const AliMUONGeometrySegmentation& rhs);
  
    // operators
    AliMUONGeometrySegmentation& operator=(const AliMUONGeometrySegmentation & rhs);

  private:
    // methods
    Bool_t Notify(Int_t detElemId) const;
  
    // data members
    mutable  Int_t                           fCurrentDetElemId;  
    mutable  AliMUONGeometryDetElement*      fCurrentDetElement;  
    mutable  AliMUONVGeometryDESegmentation* fCurrentSegmentation;
    AliMUONGeometryModule*    fGeometryModule;
    AliMUONGeometryStore*     fDESegmentations;
 
   ClassDef(AliMUONGeometrySegmentation,2) // Geometry segmentation
};

// inline functions

inline AliMUONGeometryModule* AliMUONGeometrySegmentation::GetGeometry() const
{ return fGeometryModule; }	      

#endif //ALI_MUON_GEOMETRY_SEGMENTATION_H








