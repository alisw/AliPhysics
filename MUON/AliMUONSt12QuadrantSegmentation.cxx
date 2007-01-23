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

// -------------------------------------
// Class AliMUONSt12QuadrantSegmentation
// -------------------------------------
// Segmentation for MUON quadrants of stations 1 and 2 using 
// the mapping package
// Author: Ivana Hrivnacova, IPN Orsay
 
#include <TError.h>
#include <TF1.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TVector2.h>
#include <TSystem.h>

#include "AliLog.h"

#include "AliMpPad.h"
#include "AliMpArea.h"
#include "AliMpSectorReader.h"
#include "AliMpSector.h"
#include "AliMpVPadIterator.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpFiles.h"
#include "AliMpNeighboursPadIterator.h"

#include "AliMUONSt12QuadrantSegmentation.h"
#include "AliMUONConstants.h"

/// \cond CLASSIMP
ClassImp(AliMUONSt12QuadrantSegmentation)
/// \endcond

//______________________________________________________________________________
AliMUONSt12QuadrantSegmentation::AliMUONSt12QuadrantSegmentation(
                                       AliMpVSegmentation* segmentation,
                                       AliMpStationType stationType,
				       AliMpPlaneType planeType) 
: AliMUONVGeometryDESegmentation(),
  fStationType(stationType),
  fPlaneType(planeType),
  fSector(0),
  fSectorSegmentation(0),
  fSectorIterator(0),
  fWireD(0),
  fChamber(0),
  fId(0),
  fRmin(0.),
  fRmax(0.),
  fZ(0.),
  fIx(0),
  fIy(0),
  fX(0.),
  fY(0.),
  fZone(0),
  fXhit(0.),
  fYhit(0.),
  fIxt(0),
  fIyt(0),
  fIwt(0),
  fXt(0.),
  fYt(0.),
  fCorrA(0)
{
/// Standard constructor

  fSectorSegmentation = dynamic_cast<AliMpSectorSegmentation*>(segmentation);
  if (fSectorSegmentation)
    fSector = fSectorSegmentation->GetSector();
  else 
    AliFatal("Wrong mapping segmentation type");
    
  // Anod pitch
  if ( stationType == kStation1 )  
    fWireD = AliMUONConstants::PitchSt1();
  else if ( stationType == kStation2 )   
    fWireD = AliMUONConstants::Pitch();
  else  
    AliFatal("Wrong station type");

  fCorrA = new TObjArray(3);
  fCorrA->AddAt(0,0);
  fCorrA->AddAt(0,1);
  fCorrA->AddAt(0,2);

  AliDebug(1, Form("ctor this = %p", this) ); 
}

//______________________________________________________________________________
AliMUONSt12QuadrantSegmentation::AliMUONSt12QuadrantSegmentation() 
: AliMUONVGeometryDESegmentation(),
  fStationType(kStation1),
  fPlaneType(kBendingPlane),
  fSector(0),
  fSectorSegmentation(0),
  fSectorIterator(0),
  fWireD(0),
  fChamber(0),
  fId(0),
  fRmin(0.),
  fRmax(0.),
  fZ(0.),
  fIx(0),
  fIy(0),
  fX(0.),
  fY(0.),
  fZone(0),
  fXhit(0.),
  fYhit(0.),
  fIxt(0),
  fIyt(0),
  fIwt(0),
  fXt(0.),
  fYt(0.),
  fCorrA(0) 
{
/// Default Constructor

  AliDebug(1, Form("default (empty) ctor this = %p", this));
}

//______________________________________________________________________________
AliMUONSt12QuadrantSegmentation::~AliMUONSt12QuadrantSegmentation() 
{
/// Destructor

  AliDebug(1, Form("dtor this = %p", this));

  delete fSectorIterator;  
} 

//
// private methods
//

//______________________________________________________________________________
void AliMUONSt12QuadrantSegmentation::UpdateCurrentPadValues(const AliMpPad& pad)
{
/// Updates current pad values.

  fIx = pad.GetIndices().GetFirst();
  fIy = pad.GetIndices().GetSecond();
  fX = pad.Position().X();
  fY = pad.Position().Y();
  fZone = fSectorSegmentation->Zone(pad);
}  


//
// public methods
//

//______________________________________________________________________________
void AliMUONSt12QuadrantSegmentation::SetPadSize(Float_t /*p1*/, Float_t /*p2*/)
{
/// Set pad size Dx*Dy 

  AliFatal("Not uniform pad size.");
}

//______________________________________________________________________________
void AliMUONSt12QuadrantSegmentation::SetDAnod(Float_t d)
{
/// Set anod pitch

  fWireD = d;
}

#include "AliMpMotifMap.h"
//______________________________________________________________________________
Bool_t  AliMUONSt12QuadrantSegmentation::HasPad(Float_t x, Float_t y, Float_t /*z*/)
{ 
/// Returns true if a pad exists in the given position

  // fSector->GetMotifMap()->Print();

  AliMpPad pad = fSectorSegmentation
               ->PadByPosition(TVector2(x,y), false);

  return pad.IsValid();
}  


//______________________________________________________________________________
Bool_t  AliMUONSt12QuadrantSegmentation::HasPad(Int_t ix, Int_t iy)
{
/// Returns true if a pad with given indices exists

  AliMpPad pad = fSectorSegmentation->PadByIndices(AliMpIntPair(ix,iy), false);

  return pad.IsValid();
}  


//______________________________________________________________________________
AliMUONGeometryDirection  AliMUONSt12QuadrantSegmentation::GetDirection()
{
/// Returns the direction with a constant pad size  
/// (Direction or coordinate where the resolution is the best)

  switch ( fSector->GetDirection() ) {
    case kX: return kDirX;
    case kY: return kDirY;
    default: return kDirUndefined;
  }  
}  

//______________________________________________________________________________
const AliMpVSegmentation*  
AliMUONSt12QuadrantSegmentation::GetMpSegmentation() const
{
/// Returns the mapping segmentation
/// (provides access to electronics info)

  return fSectorSegmentation;
}  

//______________________________________________________________________________
Float_t AliMUONSt12QuadrantSegmentation::GetAnod(Float_t xhit) const
{
/// Anod wire coordinate closest to xhit
/// Returns for a hit position xhit the position of the nearest anode wire    
/// From AliMUONSegmentationV01.

  Float_t wire= (xhit>0) ? Int_t(xhit/fWireD) + 0.5
                         : Int_t(xhit/fWireD) - 0.5;
  return fWireD*wire;

}

//______________________________________________________________________________
void  AliMUONSt12QuadrantSegmentation::GetPadI(Float_t x, Float_t y, Float_t /*z*/, 
                                      Int_t& ix, Int_t& iy)
{					
///  Returns pad coordinates (ix,iy) for given real coordinates (x,y)

  GetPadI(x, y, ix, iy);
}
		       
//______________________________________________________________________________
void  AliMUONSt12QuadrantSegmentation::GetPadI(Float_t x, Float_t y,
                                      Int_t& ix, Int_t& iy)
{					
/// Returns pad coordinates (ix,iy) for given real coordinates (x,y)
/// If there is no pad, ix = 0, iy = 0 are returned.

  AliMpPad pad = fSectorSegmentation->PadByPosition(TVector2(x,y), true);

  ix = pad.GetIndices().GetFirst();
  iy = pad.GetIndices().GetSecond();
}
		       
//______________________________________________________________________________
void  AliMUONSt12QuadrantSegmentation::GetPadC(Int_t ix, Int_t iy, 
                                      Float_t& x, Float_t& y, Float_t& z)
{					
/// Transform from pad to real coordinates

  z = fZ;
  GetPadC(ix, iy, x , y);
}

//______________________________________________________________________________
void  AliMUONSt12QuadrantSegmentation::GetPadC(Int_t ix, Int_t iy, 
                                      Float_t& x, Float_t& y)
{					
/// Transform from pad to real coordinates
/// If there is no pad, x = 0., y = 0. are returned.

  AliMpPad pad = fSectorSegmentation->PadByIndices(AliMpIntPair(ix,iy), true);

  x = pad.Position().X();
  y = pad.Position().Y();
}


//______________________________________________________________________________
void AliMUONSt12QuadrantSegmentation::Init(Int_t chamber)
{
/// Initialize segmentation

 // find Npx, Npy and save this info
  
  // reference to chamber
 fRmin=AliMUONConstants::Rmin(0);
 fRmax=AliMUONConstants::Rmax(0);  
 fZ = 0;
 fId=chamber;
}
 
//______________________________________________________________________________
Float_t AliMUONSt12QuadrantSegmentation::Dpx() const
{
/// Get pad size in x

  AliFatal( "Not uniform pad size.");
  return 0.;
}

//______________________________________________________________________________
Float_t AliMUONSt12QuadrantSegmentation::Dpy() const
{
/// Get pad size in y

  AliFatal("Not uniform pad size.");
  return 0.;
}
 
//______________________________________________________________________________
Float_t AliMUONSt12QuadrantSegmentation::Dpx(Int_t isector) const
{
/// Pad size in x by sector

  return fSectorSegmentation->PadDimensions(isector).X()*2.0;
} 

//______________________________________________________________________________
Float_t AliMUONSt12QuadrantSegmentation::Dpy(Int_t isector) const
{
/// Pad size in x, y by Sector 

  return fSectorSegmentation->PadDimensions(isector).Y()*2.0;
}

//______________________________________________________________________________
Int_t AliMUONSt12QuadrantSegmentation::Npx() const
{
/// Maximum number of Pads in x
/// hard coded for the time being

  return fSectorSegmentation->MaxPadIndexX();
}

//______________________________________________________________________________
Int_t AliMUONSt12QuadrantSegmentation::Npy() const
{
/// Maximum number of Pads in y
/// hard coded for the time being

  return fSectorSegmentation->MaxPadIndexY();
}

//______________________________________________________________________________
void  AliMUONSt12QuadrantSegmentation::SetPad(Int_t ix, Int_t iy)
{
/// Set pad position.
/// Sets virtual pad coordinates, needed for evaluating pad response 
/// outside the tracking program.
/// From AliMUONSegmentationV01.

  GetPadC(ix, iy, fX, fY);
  fZone = Sector(ix, iy);
}

//______________________________________________________________________________
void  AliMUONSt12QuadrantSegmentation::SetHit(Float_t xhit, Float_t yhit, Float_t /*zhit*/)
{
/// Set hit position
/// Sets virtual hit position, needed for evaluating pad response 
/// outside the tracking program 
/// From AliMUONSegmentationV01.

  fXhit = xhit;
  fYhit = yhit;
}
    
//______________________________________________________________________________
void  AliMUONSt12QuadrantSegmentation::FirstPad(Float_t xhit, Float_t yhit, Float_t /*zhit*/, 
                                       Float_t dx, Float_t dy) 
{					 
/// Iterate over pads - initialiser

  // Sets the current pad to that located in the hit position
 
  SetHit(GetAnod(xhit), yhit, 0.); 
  
  // Enable iterating in one dimension
  if (dx == 0.)  dx = 0.01;
  if (dy == 0.)  dy = 0.01;
  
  // Delete previous iterator
  delete  fSectorIterator;
  
  fSectorIterator 
    = fSectorSegmentation
        ->CreateIterator(AliMpArea(TVector2(fXhit,fYhit),TVector2(dx,dy)));

  AliDebug(1,Form("CreateIterator area=%e,%e +- %e,%e %s",
                  fXhit,fYhit,dx,dy,PlaneTypeName(fPlaneType).Data()));
  
  fSectorIterator->First();		

  if (! fSectorIterator->IsDone())
    UpdateCurrentPadValues(fSectorIterator->CurrentItem());
}
 
//______________________________________________________________________________
void  AliMUONSt12QuadrantSegmentation::NextPad()
{
/// Iterate over pads - stepper

  fSectorIterator->Next();				

  if (! fSectorIterator->IsDone())
    UpdateCurrentPadValues(fSectorIterator->CurrentItem());
}

//______________________________________________________________________________
Int_t AliMUONSt12QuadrantSegmentation::MorePads()
{
/// Iterate over pads - condition

  if (fSectorIterator->IsDone())
    return 0;
  else
    return 1;  
}

//______________________________________________________________________________
Float_t AliMUONSt12QuadrantSegmentation::Distance2AndOffset(Int_t iX, Int_t iY, 
						   Float_t x, Float_t y, Int_t* /*dummy*/)
{					   
/// Returns the square of the distance between 1 pad
/// labelled by its channel numbers and a coordinate

  AliMpPad pad = fSectorSegmentation->PadByIndices(AliMpIntPair(iX, iY));
  
  if (!pad.IsValid())
    AliFatal("Cannot locate pad.");

  return (pad.Position() - TVector2(x, y)).Mod2();
}

//______________________________________________________________________________
void AliMUONSt12QuadrantSegmentation::GetNParallelAndOffset(Int_t /*iX*/, Int_t /*iY*/,
						   Int_t* /*Nparallel*/, Int_t* /*Offset*/)
{					   
/// Number of pads read in parallel and offset to add to x 
/// (specific to LYON, but mandatory for display)

  AliFatal( "Not yet implemented.");
}


//______________________________________________________________________________
void AliMUONSt12QuadrantSegmentation::Neighbours(Int_t iX, Int_t iY, 
                                        Int_t* Nlist, 
					Int_t Xlist[10], Int_t Ylist[10])
{					  
/// Get next neighbours 

  AliMpPad pad = fSectorSegmentation->PadByIndices(AliMpIntPair(iX,iY));
  Int_t &i = *Nlist;
  i=0;
  AliMpNeighboursPadIterator iter(fSectorSegmentation, pad, kFALSE);

  for( iter.First(); !iter.IsDone() && i<10; iter.Next()) {
    Xlist[i] = iter.CurrentItem().GetIndices().GetFirst();
    Ylist[i] = iter.CurrentItem().GetIndices().GetSecond();
    i++;
  }
}
 
//______________________________________________________________________________
Int_t  AliMUONSt12QuadrantSegmentation::Ix()
{
/// Current pad cursor during disintegration
/// x, y-coordinate

  return fSectorIterator->CurrentItem().GetIndices().GetFirst();
}

//______________________________________________________________________________
Int_t  AliMUONSt12QuadrantSegmentation::Iy()
{
/// Current pad cursor during disintegration
/// x, y-coordinate

  return fSectorIterator->CurrentItem().GetIndices().GetSecond();
}

//______________________________________________________________________________
Int_t  AliMUONSt12QuadrantSegmentation::ISector()
{
/// Current sector

  return fZone;
}

//______________________________________________________________________________
Int_t AliMUONSt12QuadrantSegmentation::Sector(Int_t ix, Int_t iy)
{
/// Calculate sector from pad coordinates

  return fSectorSegmentation
           ->Zone(fSectorSegmentation->PadByIndices(AliMpIntPair(ix, iy)));
}

//______________________________________________________________________________
Int_t AliMUONSt12QuadrantSegmentation::Sector(Float_t x, Float_t y)
{
/// Calculate sector from pad coordinates

  return fSectorSegmentation
           ->Zone(fSectorSegmentation
	            ->PadByPosition(TVector2(x,y)));
}

//______________________________________________________________________________
void  AliMUONSt12QuadrantSegmentation::IntegrationLimits(Float_t& x1, Float_t& x2,
                                                Float_t& y1, Float_t& y2)
{						  
/// Current integration limits 
 
  x1 = fXhit - fX - Dpx(fZone)/2.;
  x2 = x1 + Dpx(fZone);

  y1 = fYhit - fY - Dpy(fZone)/2.;
  y2 = y1 + Dpy(fZone);    
}

//______________________________________________________________________________
Int_t AliMUONSt12QuadrantSegmentation::SigGenCond(Float_t x, Float_t y, Float_t /*z*/)
{
/// Signal Generation Condition during Stepping
/// -  0: don't generate signal
/// -  1: generate signal 
///
///  Comments:                                                                  \n 
///
///  Crossing of pad boundary and mid plane between neighbouring wires is checked.
///  To correctly simulate the dependence of the spatial resolution on the angle 
///  of incidence signal must be generated for constant steps on 
///  the projection of the trajectory along the anode wire.
///
///  Signal will be generated if particle crosses pad boundary or
///  boundary between two wires. 
///
/// From AliMUONSegmentationV01

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
void  AliMUONSt12QuadrantSegmentation::SigGenInit(Float_t x, Float_t y, Float_t /*z*/)
{
/// Initialise signal generation at coord (x,y,z)
/// Initialises pad and wire position during stepping.
/// From AliMUONSegmentationV01

  fXt = x;
  fYt = y;
  GetPadI(x, y, fIxt, fIyt);
  fIwt = (x>0) ? Int_t(x/fWireD)+1 : Int_t(x/fWireD)-1 ;
}		    
    
//______________________________________________________________________________
void AliMUONSt12QuadrantSegmentation::GiveTestPoints(Int_t& n, Float_t* x, Float_t* y) const
{					      
/// Test points for auto calibration
/// Returns test point on the pad plane.
/// Used during determination of the segmoid correction of the COG-method
/// From AliMUONSegmentationV01

  n=1;
  x[0] = (fRmax+fRmin)/2/TMath::Sqrt(2.);
  y[0] = x[0];
}

//______________________________________________________________________________
void AliMUONSt12QuadrantSegmentation::Draw(const char * /*opt*/)
{
/// Draw the segmentation zones.
/// (Called from AliMUON::BuildGeometry)

  AliWarning("Not yet implemented.");
}

//______________________________________________________________________________
void AliMUONSt12QuadrantSegmentation::SetCorrFunc(Int_t isec, TF1* func)
{
/// Set the correction function.
/// From AliMUONSegmentationV01

  fCorrA->AddAt(func, isec);
}

//______________________________________________________________________________
TF1* AliMUONSt12QuadrantSegmentation::CorrFunc(Int_t isec) const
{
/// Get the correction Function.
/// From AliMUONSegmentationV01

  return (TF1*) fCorrA->At(isec);
} 
