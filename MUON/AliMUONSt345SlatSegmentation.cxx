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

#include "AliMUONSt345SlatSegmentation.h"
#include "AliMUONConstants.h"

#include "AliMpArea.h"
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpPCB.h" 
// FIXME: we should not need access to PCB at this level
// if the Dpx, Dpy, Sector interface would not be used.
// Investigate instead to direct use of AliMpVSegmentation and AliMpPad
// from the clusterFinder ? // :aphecetc:20050722 
// or, or ... have the AliMpSlat handles the notion of zone, finally
// (where "zone" in mapping is equivalent to "sector" here), and thus
// let AliMpSlat offer the interface to get information similar to what
// Dpx, Dpy, Sector offer now to the clustering. That indeed should be
// handled directly at the level of AliMpVSegmentation...
#include "AliMpVPadIterator.h"

#include "AliLog.h"

#include "Riostream.h"

/// \cond CLASSIMP
ClassImp(AliMUONSt345SlatSegmentation)
/// \endcond

namespace
{
  Float_t FMAX(1E9);
}

//_____________________________________________________________________________
AliMUONSt345SlatSegmentation::AliMUONSt345SlatSegmentation()
: AliMUONVGeometryDESegmentation(),
fDetElemId(-1),
fPlaneType(kBendingPlane),
fSlat(0),
fSlatSegmentation(0),
fPadIterator(0),
fCurrentPad(),
fXhit(FMAX),
fYhit(FMAX)
{
/// Default ctor

	AliDebug(1,Form("this=%p default (empty) ctor",this));
}

//_____________________________________________________________________________
AliMUONSt345SlatSegmentation::AliMUONSt345SlatSegmentation(
                                   AliMpVSegmentation* segmentation,
                                   Int_t detElemId, AliMpPlaneType bendingOrNonBending)
: AliMUONVGeometryDESegmentation(),
fDetElemId(detElemId),
fPlaneType(bendingOrNonBending),
fSlat(0),
fSlatSegmentation(0),
fPadIterator(0),
fCurrentPad(),
fXhit(FMAX),
fYhit(FMAX)
{ 
/// Normal ctor.

  fSlatSegmentation = dynamic_cast<AliMpSlatSegmentation*>(segmentation);
  if (fSlatSegmentation)
    fSlat = fSlatSegmentation->Slat();
  else 
    AliFatal("Wrong mapping segmentation type");
		
  AliDebug(1,Form("this=%p detElemId=%3d %s fSlatSegmentation=%p",this,detElemId,
									( (bendingOrNonBending==kBendingPlane)?"Bending":"NonBending" ),
									fSlatSegmentation));
}

//_____________________________________________________________________________
AliMUONSt345SlatSegmentation::~AliMUONSt345SlatSegmentation()
{
/// Destructor

	AliDebug(1,Form("dtor this=%p",this));
  delete fPadIterator;
}

//_____________________________________________________________________________
TF1*
AliMUONSt345SlatSegmentation::CorrFunc(Int_t /*isec*/) const
{
/// Not implemented

  AliFatal("Not Implemented");
  return 0x0;
}

//_____________________________________________________________________________
Float_t 
AliMUONSt345SlatSegmentation::Distance2AndOffset(Int_t /*iX*/, Int_t /*iY*/, 
			          Float_t /*x*/, Float_t /*y*/, Int_t* /*dummy*/){
/// Not implemented

  AliFatal("Not Implemented");
  return 0.0;
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::Draw(Option_t* /*opt*/)
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Float_t
AliMUONSt345SlatSegmentation::Dpx() const
{
/// Not implemented

  AliFatal("Not Implemented");
  return 0.0;
}

//_____________________________________________________________________________
Float_t
AliMUONSt345SlatSegmentation::Dpy() const
{
/// Not implemented

  AliFatal("Not Implemented");
  return 0.0;
}

//_____________________________________________________________________________
Float_t
AliMUONSt345SlatSegmentation::Dpx(int ipcb) const
{
/// Get pad size in x

	AliMpPCB* pcb = fSlat->GetPCB(ipcb);
	if (!pcb) 
  {
    AliFatal("pcb is null!");
  }
	return pcb->PadSizeX();
}

//_____________________________________________________________________________
Float_t
AliMUONSt345SlatSegmentation::Dpy(int ipcb) const
{
/// Get pad size in y

	AliMpPCB* pcb = fSlat->GetPCB(ipcb);
	if (!pcb) 
  {
    AliFatal("pcb is null!");
  }
	return pcb->PadSizeY();
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::FirstPad(Float_t xhit, Float_t yhit,Float_t /*zhit*/,
                                         Float_t dx, Float_t dy)
{
/// OK. We will do it in 2 steps. First get the area over which to
/// iterate, based on hit coordinates and (dx,dy). This first step
/// has nothing to do with segmentation in the first place, but with
/// how we simulate the fact that at some point the charge is shared
/// amongst several pads.
/// The second step is the actual pad iteration and is handled by 
/// a specific class (which has to do only with iteration...)
///
/// \todo FIXME: this method should not be here in the first place, IMHO.
	
  // Find the wire position (center of charge distribution)
  Float_t xwire = GetAnod(xhit);
  fXhit = xwire;
  fYhit = yhit;
	
  Double_t x01 = xwire - dx;
  Double_t x02 = xwire + dx;
  Double_t y01 = yhit - dy;
  Double_t y02 = yhit + dy;
	
  Double_t xext = x02 - x01;
  Double_t yext = y02 - y01;
	
  // we do not check area here, as we assume the iterator
  // will do it, and will possibly truncate it if needed.
  // Note that we convert the area position to a reference frame
  // located in the lower-left corner of the slat, instead of its
  // center.
//  AliMpArea area(TVector2((x01+x02)/2.0+fSlat->DX(),
//													(y01+y02)/2.0+fSlat->DY()),
//								 TVector2(xext/2.0,yext/2.0));
  AliMpArea area(TVector2((x01+x02)/2.0,(y01+y02)/2.0),
								 TVector2(xext/2.0,yext/2.0));
	
  delete fPadIterator;
	
  fPadIterator = fSlatSegmentation->CreateIterator(area);
	
  fPadIterator->First();
	
  fCurrentPad = fPadIterator->CurrentItem();
	
  if ( !fCurrentPad.IsValid() ) 
	{
		AliError(Form("Cannot get a valid pad for (xhit,yhit,dx,dy)=(%e,%e,%e,%e)",xhit,yhit,dx,dy));
	}
	
  AliDebug(4,Form("xhit,yhit,dx,dy=%e,%e,%e,%e ix,iy=%3d,%3d slat=%s",
									xhit,yhit,dx,dy,
									fCurrentPad.GetIndices().GetFirst(),
									fCurrentPad.GetIndices().GetSecond(),fSlat->GetID()));		  
}

//_____________________________________________________________________________
Float_t
AliMUONSt345SlatSegmentation::GetAnod(Float_t xhit) const
{
/// Gets the x-coordinate of the wire which is the closest to xhit.
	
  Int_t n = Int_t(xhit/AliMUONConstants::Pitch());
  Float_t wire = (xhit>0) ? n+0.5 : n-0.5;
  return AliMUONConstants::Pitch()*wire;
}

//_____________________________________________________________________________
AliMUONGeometryDirection
AliMUONSt345SlatSegmentation::GetDirection()
{
/// Not implemented

  //AliWarning("Not Implemented");
  return kDirUndefined;
}

//______________________________________________________________________________
const AliMpVSegmentation*  
AliMUONSt345SlatSegmentation::GetMpSegmentation() const
{
/// Returns the mapping segmentation
/// (provides access to electronics info)

  return fSlatSegmentation;
}  


//_____________________________________________________________________________
void 
AliMUONSt345SlatSegmentation::GetNParallelAndOffset(Int_t /*iX*/, Int_t /*iY*/,
				         Int_t* /*Nparallel*/, Int_t* /*Offset*/)
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::GetPadC(Int_t ix, Int_t iy, 
                                        Float_t& x, Float_t& y, Float_t& z)
{					 
/// Transform from pad to real coordinates

  z = 0;
  GetPadC(ix,iy,x,y);
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::GetPadC(Int_t ix, Int_t iy, 
                                        Float_t& x, Float_t& y)
{
/// Transform from pad to real coordinates

  AliMpPad pad = 
  fSlatSegmentation->PadByIndices(AliMpIntPair(ix,iy),kTRUE);
  x = pad.Position().X();
  y = pad.Position().Y();
}


//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::GetPadI(Float_t x, Float_t y, Float_t /*z*/,
Int_t& ix, Int_t& iy)
{
///  Returns pad coordinates (ix,iy) for given real coordinates (x,y)

  GetPadI(x,y,ix,iy);
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::GetPadI(Float_t x, Float_t y,
                                        Int_t& ix, Int_t& iy)
{
///  Returns pad coordinates (ix,iy) for given real coordinates (x,y)

  AliMpPad pad = fSlatSegmentation->PadByPosition(TVector2(x,y), kTRUE);
	
  if ( pad != AliMpPad::Invalid() )
	{
		ix = pad.GetIndices().GetFirst();
		iy = pad.GetIndices().GetSecond();
	}
  else
	{
		ix = iy = -1;
	}
}

//_____________________________________________________________________________
void 
AliMUONSt345SlatSegmentation::GiveTestPoints(Int_t& /*n*/, 
                                    Float_t* /*x*/, Float_t* /*y*/) const
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Bool_t
AliMUONSt345SlatSegmentation::HasPad(Float_t x, Float_t y, Float_t z)
{
/// Returns true if a pad exists in the given position

  Int_t ix, iy;
  GetPadI(x,y,z,ix,iy);
  return HasPad(ix,iy);
}

//_____________________________________________________________________________
Bool_t
AliMUONSt345SlatSegmentation::HasPad(Int_t ix, Int_t iy)
{
/// Returns true if a pad with given indices exists

  return fSlatSegmentation->HasPad(AliMpIntPair(ix,iy));
}

//_____________________________________________________________________________
void  
AliMUONSt345SlatSegmentation::IntegrationLimits(Float_t& x1, Float_t& x2,
                                                  Float_t& y1, Float_t& y2)
{
///  Returns integration limits for current pad
	
	//   x1 = fXhit - fX - Dpx(fSector)/2.;
	//   x2 = x1 + Dpx(fSector);
	//   y1 = fYhit - fY - Dpy(fSector)/2.;
	//   y2 = y1 + Dpy(fSector);    
	
  Float_t x = fCurrentPad.Position().X();
  Float_t y = fCurrentPad.Position().Y();
	
  Float_t padsizex = fCurrentPad.Dimensions().X() * 2.0;
  Float_t padsizey = fCurrentPad.Dimensions().Y() * 2.0;
	
  x1 = fXhit - x - padsizex/2.0;
  x2 = x1 + padsizex;
  y1 = fYhit - y - padsizey/2.0;
  y2 = y1 + padsizey;
	
  AliDebug(4,Form("xhit,yhit=%e,%e x,y=%e,%e, x1,x2,y1,y2=%e,%e,%e,%e",fXhit,fYhit,x,y,x1,x2,y1,y2));
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentation::ISector()
{
/// \todo FIXME: remove the usage of ISector from all the code.

  return -10;
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentation::Ix()
{
/// Current pad cursor during disintegration
/// x, y-coordinate

  if ( fPadIterator )
	{
		return fPadIterator->CurrentItem().GetIndices().GetFirst();
	}
  else
	{
		return -1;
	}
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentation::Iy()
{
/// Current pad cursor during disintegration
/// x, y-coordinate

  if ( fPadIterator ) 
	{
		return fPadIterator->CurrentItem().GetIndices().GetSecond();
	}
  else
	{
		return -1;
	}
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentation::MorePads()
{
/// Iterate over pads - condition

  return (fPadIterator && !fPadIterator->IsDone());
}

//_____________________________________________________________________________
void 
AliMUONSt345SlatSegmentation::Neighbours(Int_t iX, Int_t iY, Int_t* Nlist,
                                           Int_t Xlist[10], Int_t Ylist[10])
{
/// Find pad at (ix,iy) for which we'll search neighbours.

  AliMpPad pad = 
	fSlatSegmentation->PadByIndices(AliMpIntPair(iX,iY),kTRUE);
	
  // Define the region to look into : a region slightly bigger
  // than the pad itself (5% bigger), in order to catch first neighbours.

  AliMpArea area(pad.Position(),pad.Dimensions()*1.05); 
		
  AliMpVPadIterator* it = fSlatSegmentation->CreateIterator(area);
  it->First();
  Int_t n = 0;
  while ( !it->IsDone() && n < 10 )
	{
		AliMpPad p = it->CurrentItem();
		if ( p != pad ) // skip self
		{
			Xlist[n] = p.GetIndices().GetFirst();
			Ylist[n] = p.GetIndices().GetSecond();
			++n;
		}
		it->Next();
	}
	delete it;
  *Nlist = n;
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::NextPad()
{
/// Iterate over pads - stepper

  if ( fPadIterator )
	{
		fPadIterator->Next();
		fCurrentPad = fPadIterator->CurrentItem();
	}
  else
	{
		AliError("PadIterator not initialized. Please use First() first ;-)");
	}
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentation::Npx() const
{
/// Maximum number of Pads in x

  return fSlatSegmentation->MaxPadIndexX()+1;
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentation::Npy() const
{
/// Maximum number of Pads in y

  return fSlatSegmentation->MaxPadIndexY()+1;
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::Print(Option_t* /*opt*/) const
{
/// Printing

  cout << "DetElemId=" << fDetElemId << " PlaneType=" 
  << fPlaneType << " Npx,Npy=" << Npx() << "," << Npy() << " fSlat=" << fSlat 
  << " fSlatSegmentation=" << fSlatSegmentation
  << " fSlatSegmentation->Slat()=" << fSlatSegmentation->Slat() << endl;
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentation::Sector(Int_t ix, Int_t /*iy*/)
{
/// Calculate sector from pad coordinates

  return fSlat->FindPCBIndex(ix);
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentation::Sector(Float_t x, Float_t y)
{
/// Calculate sector from pad coordinates

  return fSlat->FindPCBIndex(x,y);
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::SetCorrFunc(Int_t /*isec*/,TF1* /*func*/)
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::SetDAnod(float /*d*/)
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::SetHit(Float_t x, Float_t y, Float_t)
{
/// Set hit position
/// Sets virtual hit position, needed for evaluating pad response 
/// outside the tracking program 

  fXhit = x;
  fYhit = y;
	
  //
  // insure we're within the slat limits. If not, issue a simple
  // warning, as this might be correct (due to clustering/fitting algorithm
  // that is allowed to go a little bit outside limits).
  // That should only be considered an error in case the point is way
  // outside (but by how much is the question you'll have to determine
  // by yourself...)
  //
  if ( fXhit < -fSlat->DX() || fXhit > fSlat->DX() ||
       fYhit < -fSlat->DY() || fYhit > fSlat->DY() )
	{
    Double_t dx = - fSlat->DX() + TMath::Abs(fXhit);
    Double_t dy = - fSlat->DY() + TMath::Abs(fYhit);
    dx = ( dx > 0 ? dx : 0);
    dy = ( dy > 0 ? dy : 0);
		AliWarning(Form("Hit outside slat %s limits (x,y)hit = (%e,%e). "
                    "Outside by (%e,%e) cm. Might be ok though.",
                    fSlat->GetID(),fXhit,fYhit,dx,dy));
	}  
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::SetPad(Int_t ix, Int_t iy)
{
/// Set pad position.
/// Sets virtual pad coordinates, needed for evaluating pad response 
/// outside the tracking program.

  fCurrentPad = 
	fSlatSegmentation->PadByIndices(AliMpIntPair(ix,iy),kTRUE);
  if ( !fCurrentPad.IsValid() )
	{
		AliError(Form("Setting current pad to invalid ! (ix,iy)=(%4d,%4d)",ix,iy));
	}
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentation::SetPadSize(float /*p1*/,float /*p2*/)
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Int_t 
AliMUONSt345SlatSegmentation::SigGenCond(Float_t /*x*/, Float_t /*y*/, Float_t /*z*/)
{
/// Not implemented

  AliFatal("Not Implemented");
  return 0;
}

//_____________________________________________________________________________
void 
AliMUONSt345SlatSegmentation::SigGenInit(Float_t,Float_t,Float_t)
{
/// Not implemented

  AliFatal("Not Implemented");
}
