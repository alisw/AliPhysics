/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 
/* $Id$ */

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONSt1Segmentation
// ------------------------------
// Segmentation for MUON station 1 using the external 
// mapping package
// Included in AliRoot 2003/01/28

#include <TError.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TVector2.h>

#include "AliMpPad.h"
#include "AliMpPlane.h"
#include "AliMpPlaneAreaPadIterator.h"
#include "AliMpPlaneSegmentation.h"

#include "AliMUONSt1Segmentation.h"
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"

ClassImp(AliMUONSt1Segmentation)

const Float_t  AliMUONSt1Segmentation::fgkWireD = 0.20; 
const Float_t  AliMUONSt1Segmentation::fgkLengthUnit = 0.1; 

//______________________________________________________________________________
AliMUONSt1Segmentation::AliMUONSt1Segmentation(const AliMpPlaneType planeType) 
: AliSegmentation(),
  fPlane(0),
  fPlaneSegmentation(0),
  fPlaneIterator(0),
  fWireD(fgkWireD),
  fChamber(0),
  fId(0),
  fRmin(0.),
  fRmax(0.),
  fZ(0.),
  fIx(0),
  fIy(0),
  fX(0.),
  fY(0.),
  fSector(0),
  fXhit(0.),
  fYhit(0.),
  fIxt(0),
  fIyt(0),
  fIwt(0),
  fXt(0.),
  fYt(0.),
  fCorrA(0)
{
// Normal constructor
  fPlane = AliMpPlane::Create(kStation1, planeType);
  fPlaneSegmentation = new AliMpPlaneSegmentation(fPlane);

  fCorrA = new TObjArray(3);
  fCorrA->AddAt(0,0);
  fCorrA->AddAt(0,1);
  fCorrA->AddAt(0,2);
}

//______________________________________________________________________________
AliMUONSt1Segmentation::AliMUONSt1Segmentation() 
: AliSegmentation(),
  fPlane(0),
  fPlaneSegmentation(0),
  fPlaneIterator(0),
  fWireD(fgkWireD),
  fChamber(0),
  fId(0),
  fRmin(0.),
  fRmax(0.),
  fZ(0.),
  fIx(0),
  fIy(0),
  fX(0.),
  fY(0.),
  fSector(0),
  fXhit(0.),
  fYhit(0.),
  fIxt(0),
  fIyt(0),
  fIwt(0),
  fXt(0.),
  fYt(0.),
  fCorrA(0) {
// Default Constructor
}

//______________________________________________________________________________
AliMUONSt1Segmentation::AliMUONSt1Segmentation(const AliMUONSt1Segmentation& rhs) :AliSegmentation(rhs)
{
// Copy constructor
  Fatal("Copy constructor", 
        "Copy constructor is not implemented.");
}

//______________________________________________________________________________
AliMUONSt1Segmentation::~AliMUONSt1Segmentation() {
// Destructor

  delete fPlane;
  delete fPlaneSegmentation;  
  delete fPlaneIterator;  
} 

//
// operators
//

//______________________________________________________________________________
AliMUONSt1Segmentation& 
AliMUONSt1Segmentation::operator=(const AliMUONSt1Segmentation& rhs)
{
// Copy operator 

  // check assignement to self
  if (this == &rhs) return *this;

  Fatal("operator=", 
        "Assignment operator is not implemented.");
    
  return *this;  
}

//
// private methods
//

//______________________________________________________________________________
void AliMUONSt1Segmentation::UpdateCurrentPadValues(const AliMpPad& pad)
{
// Updates current pad values.
// ---

  fIx = pad.GetIndices().GetFirst();
  fIy = pad.GetIndices().GetSecond();
  fX = pad.Position().X() * fgkLengthUnit;
  fY = pad.Position().Y() * fgkLengthUnit;
  fSector = fPlaneSegmentation->Zone(pad);
}  

//
// public methods
//

//______________________________________________________________________________
void AliMUONSt1Segmentation::SetPadSize(Float_t /*p1*/, Float_t /*p2*/)
{
// Set pad size Dx*Dy 
// ---

  Fatal("SetPadSize", "Not uniform pad size.");
}

//______________________________________________________________________________
void AliMUONSt1Segmentation::SetDAnod(Float_t d)
{
// Set anod pitch
// ---

  fWireD = d;
}

//______________________________________________________________________________
Float_t AliMUONSt1Segmentation::GetAnod(Float_t xhit) const
{
// Anod wire coordinate closest to xhit
// Returns for a hit position xhit the position of the nearest anode wire    
// From AliMUONSegmentationV01.
// ---

  Float_t wire= (xhit>0) ? Int_t(xhit/fWireD) + 0.5
                         : Int_t(xhit/fWireD) - 0.5;
  return fWireD*wire;

}

//______________________________________________________________________________
void  AliMUONSt1Segmentation::GetPadI(Float_t x, Float_t y, Float_t /*z*/, 
                                      Int_t& ix, Int_t& iy)
{					
//  Returns pad coordinates (ix,iy) for given real coordinates (x,y)
// ---

  GetPadI(x, y, ix, iy);
}
		       
//______________________________________________________________________________
void  AliMUONSt1Segmentation::GetPadI(Float_t x, Float_t y,
                                      Int_t& ix, Int_t& iy)
{					
// Returns pad coordinates (ix,iy) for given real coordinates (x,y)
// If there is no pad, ix = 0, iy = 0 are returned.
// ---

  AliMpPad pad = fPlaneSegmentation
               ->PadByPosition(TVector2(x/fgkLengthUnit, y/fgkLengthUnit), false);

  ix = pad.GetIndices().GetFirst();
  iy = pad.GetIndices().GetSecond();
}
		       
//______________________________________________________________________________
void  AliMUONSt1Segmentation::GetPadC(Int_t ix, Int_t iy, 
                                      Float_t& x, Float_t& y, Float_t& z)
{					
// Transform from pad to real coordinates
// ---

  z = fZ;
  GetPadC(ix, iy, x , y);
}

//______________________________________________________________________________
void  AliMUONSt1Segmentation::GetPadC(Int_t ix, Int_t iy, 
                                      Float_t& x, Float_t& y)
{					
// Transform from pad to real coordinates
// If there is no pad, x = 0., y = 0. are returned.
// ---

  AliMpPad pad = fPlaneSegmentation->PadByIndices(AliMpIntPair(ix,iy));

  x = pad.Position().X() * fgkLengthUnit;
  y = pad.Position().Y() * fgkLengthUnit;
}


//______________________________________________________________________________
void AliMUONSt1Segmentation::Init(Int_t chamber)
{
// Initialize segmentation
// ---

  // find Npx, Npy and save this info
  
  // reference to chamber
 AliMUON *pMUON  = (AliMUON *) gAlice->GetModule("MUON");
 fChamber = &(pMUON->Chamber(chamber));
 fRmin=fChamber->RInner();
 fRmax=fChamber->ROuter();  
 fZ = fChamber->Z();
 fId=chamber;
}
 
//______________________________________________________________________________
Float_t AliMUONSt1Segmentation::Dpx() const
{
// Get pad size in x
// ---

  Fatal("Dpx", "Not uniform pad size.");
  return 0.;
}

//______________________________________________________________________________
Float_t AliMUONSt1Segmentation::Dpy() const
{
// Get pad size in y
// ---

  Fatal("Dpy", "Not uniform pad size.");
  return 0.;
}
 
//______________________________________________________________________________
Float_t AliMUONSt1Segmentation::Dpx(Int_t isector) const
{
// Pad size in x by sector
// ---

  return fPlaneSegmentation->PadDimensions(isector).X()*2.*fgkLengthUnit;
} 

//______________________________________________________________________________
Float_t AliMUONSt1Segmentation::Dpy(Int_t isector) const
{
// Pad size in x, y by Sector 
// ---

  return fPlaneSegmentation->PadDimensions(isector).Y()*2.*fgkLengthUnit;
}

//______________________________________________________________________________
Int_t AliMUONSt1Segmentation::Npx() const
{
// Maximum number of Pads in x
// ---

  //Fatal("Npx", "Not yet implemented.");
  return 142; //hard coded for the time being
}

//______________________________________________________________________________
Int_t AliMUONSt1Segmentation::Npy() const
{
// Maximum number of Pads in y
// ---

  //Fatal("Npy", "Not yet implemented.");
  return 213; //hard coded for the time being
}

//______________________________________________________________________________
void  AliMUONSt1Segmentation::SetPad(Int_t ix, Int_t iy)
{
// Set pad position.
// Sets virtual pad coordinates, needed for evaluating pad response 
// outside the tracking program.
// From AliMUONSegmentationV01.
// ---

  GetPadC(ix, iy, fX, fY);
  fSector = Sector(ix, iy);
}

//______________________________________________________________________________
void  AliMUONSt1Segmentation::SetHit(Float_t xhit, Float_t yhit, Float_t /*zhit*/)
{
// Set hit position
// Sets virtual hit position, needed for evaluating pad response 
// outside the tracking program 
// From AliMUONSegmentationV01.

  fXhit = xhit;
  fYhit = yhit;
}
    
//______________________________________________________________________________
void  AliMUONSt1Segmentation::FirstPad(Float_t xhit, Float_t yhit, Float_t /*zhit*/, 
                                       Float_t dx, Float_t dy) 
{					 
// Iterate over pads - initialiser
// ---

  // Sets the current pad to that located in the hit position
 
  SetHit(GetAnod(xhit), yhit, 0.); 
  
  // Enable iterating in one dimension
  if (dx == 0.)  dx = 0.01;
  if (dy == 0.)  dy = 0.01;
  
  fPlaneIterator 
    = fPlaneSegmentation
        ->CreateIterator(AliMpArea(TVector2(fXhit/fgkLengthUnit, fYhit/fgkLengthUnit), 
	                           TVector2(dx/fgkLengthUnit, dy/fgkLengthUnit)));

  fPlaneIterator->First();		

  if (! fPlaneIterator->IsDone())
    UpdateCurrentPadValues(fPlaneIterator->CurrentItem());
}
 
//______________________________________________________________________________
void  AliMUONSt1Segmentation::NextPad()
{
// Iterate over pads - stepper
// ---

  fPlaneIterator->Next();				

  if (! fPlaneIterator->IsDone())
    UpdateCurrentPadValues(fPlaneIterator->CurrentItem());
}

//______________________________________________________________________________
Int_t AliMUONSt1Segmentation::MorePads()
{
// Iterate over pads - condition
// ---

  if (fPlaneIterator->IsDone())
    return 0;
  else
    return 1;  
}

//______________________________________________________________________________
Float_t AliMUONSt1Segmentation::Distance2AndOffset(Int_t iX, Int_t iY, 
						   Float_t x, Float_t y, Int_t* /*dummy*/)
{					   
// Returns the square of the distance between 1 pad
// labelled by its channel numbers and a coordinate
// ---

  AliMpPad pad = fPlaneSegmentation->PadByIndices(AliMpIntPair(iX, iY));
  
  if (!pad.IsValid())
    Fatal("Distance2AndOffset", "Cannot locate pad.");

  return (pad.Position()*fgkLengthUnit - TVector2(x, y)).Mod2();
}

//______________________________________________________________________________
void AliMUONSt1Segmentation::GetNParallelAndOffset(Int_t /*iX*/, Int_t /*iY*/,
						   Int_t* /*Nparallel*/, Int_t* /*Offset*/)
{					   
// Number of pads read in parallel and offset to add to x 
// (specific to LYON, but mandatory for display)
// ---

  Fatal("GetNParallelAndOffset", "Not yet implemented.");
}


//______________________________________________________________________________
void AliMUONSt1Segmentation::Neighbours(Int_t iX, Int_t iY, 
                                        Int_t* Nlist, 
					Int_t Xlist[10], Int_t Ylist[10])
{					  
// Get next neighbours 
// ---

  AliMpPad pad = fPlaneSegmentation->PadByIndices(AliMpIntPair(iX,iY));
  Int_t &i = *Nlist;
  i=0;
  AliMpVPadIterator* iter
    = fPlaneSegmentation
      ->CreateIterator(AliMpArea(pad.Position(),2.*pad.Dimensions()*1.1));

  for( iter->First(); !iter->IsDone() && i<10; iter->Next()) {
    Xlist[i] = iter->CurrentItem().GetIndices().GetFirst();
    Ylist[i] = iter->CurrentItem().GetIndices().GetSecond();
    i++;
  }
  
  delete iter;
}

//______________________________________________________________________________
Int_t  AliMUONSt1Segmentation::Ix()
{
// Current pad cursor during disintegration
// x, y-coordinate
// ---

  return fPlaneIterator->CurrentItem().GetIndices().GetFirst();
}

//______________________________________________________________________________
Int_t  AliMUONSt1Segmentation::Iy()
{
// Current pad cursor during disintegration
// x, y-coordinate
// ---

  return fPlaneIterator->CurrentItem().GetIndices().GetSecond();
}

//______________________________________________________________________________
Int_t  AliMUONSt1Segmentation::ISector()
{
// Current sector
// ---

  return fSector;
}

//______________________________________________________________________________
Int_t AliMUONSt1Segmentation::Sector(Int_t ix, Int_t iy)
{
// Calculate sector from pad coordinates
// ---

  return fPlaneSegmentation
           ->Zone(fPlaneSegmentation->PadByIndices(AliMpIntPair(ix, iy)));
}

//______________________________________________________________________________
Int_t AliMUONSt1Segmentation::Sector(Float_t x, Float_t y)
{
// Calculate sector from pad coordinates
// ---

  return fPlaneSegmentation
           ->Zone(fPlaneSegmentation
	            ->PadByPosition(TVector2(x/fgkLengthUnit, y/fgkLengthUnit)));
}

//______________________________________________________________________________
void  AliMUONSt1Segmentation::IntegrationLimits(Float_t& x1, Float_t& x2,
                                                Float_t& y1, Float_t& y2)
{						  
// Current integration limits 
// ---
 
  x1 = fXhit - fX - Dpx(fSector)/2.;
  x2 = x1 + Dpx(fSector);

  y1 = fYhit - fY - Dpy(fSector)/2.;
  y2 = y1 + Dpy(fSector);    
}

//______________________________________________________________________________
Int_t AliMUONSt1Segmentation::SigGenCond(Float_t x, Float_t y, Float_t /*z*/)
{
// Signal Generation Condition during Stepping
//  0: don't generate signal
//  1: generate signal 
//  Comments: 
//
//  Crossing of pad boundary and mid plane between neighbouring wires is checked.
//  To correctly simulate the dependence of the spatial resolution on the angle 
//  of incidence signal must be generated for constant steps on 
//  the projection of the trajectory along the anode wire.
//
//  Signal will be generated if particle crosses pad boundary or
//  boundary between two wires. 
//
// From AliMUONSegmentationV01
// ---

  Int_t ixt, iyt;
  GetPadI(x, y, ixt, iyt);
  Int_t iwt=(x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1;
  
  if ((ixt != fIxt) || (iyt !=fIyt) || (iwt != fIwt)) {
    return 1;
  } 
  else {
    return 0;
  }
}


//______________________________________________________________________________
void  AliMUONSt1Segmentation::SigGenInit(Float_t x, Float_t y, Float_t /*z*/)
{
// Initialise signal generation at coord (x,y,z)
// Initialises pad and wire position during stepping.
// From AliMUONSegmentationV01
// ---

  fXt = x;
  fYt = y;
  GetPadI(x, y, fIxt, fIyt);
  fIwt = (x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1 ;
}		    
    
//______________________________________________________________________________
void AliMUONSt1Segmentation::GiveTestPoints(Int_t& n, Float_t* x, Float_t* y) const
{					      
// Test points for auto calibration
// Returns test point on the pad plane.
// Used during determination of the segmoid correction of the COG-method
// From AliMUONSegmentationV01
// ---

  n=1;
  x[0] = (fRmax+fRmin)/2/TMath::Sqrt(2.);
  y[0] = x[0];
}

//______________________________________________________________________________
void AliMUONSt1Segmentation::Draw(const char * /*opt*/) const
{
// Draw the segmentation zones.
// (Called from AliMUON::BuildGeometry)
// ---

  Warning("Draw", "Not yet implemented.");
}

//______________________________________________________________________________
void AliMUONSt1Segmentation::SetCorrFunc(Int_t isec, TF1* func)
{
// Set the correction function.
// From AliMUONSegmentationV01
// ---

  fCorrA->AddAt(func, isec);
}

//______________________________________________________________________________
TF1* AliMUONSt1Segmentation::CorrFunc(Int_t isec) const
{
// Get the correction Function.
// From AliMUONSegmentationV01
// ---

  return (TF1*) fCorrA->At(isec);
} 

