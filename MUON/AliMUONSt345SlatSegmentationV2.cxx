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

#include "AliMUONSt345SlatSegmentationV2.h"

#include "AliLog.h"
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
#include "AliMUONSegmentationManager.h"
#include "AliMUONConstants.h"

#include "Riostream.h"

#include "TClass.h"

ClassImp(AliMUONSt345SlatSegmentationV2)

namespace
{
  Int_t fgIntOffset = 1;
  Float_t FMAX(1E9);
}

//_____________________________________________________________________________
AliMUONSt345SlatSegmentationV2::AliMUONSt345SlatSegmentationV2()
: AliMUONVGeometryDESegmentation(),
fDetElemId(-1),
fPlaneType(kBendingPlane),
fSlat(0),
fSlatSegmentation(0),
fPadIterator(0),
fXhit(FMAX),
fYhit(FMAX)
{
	AliDebug(1,Form("this=%p default (empty) ctor",this));
}

//_____________________________________________________________________________
AliMUONSt345SlatSegmentationV2::AliMUONSt345SlatSegmentationV2(Int_t detElemId, AliMpPlaneType bendingOrNonBending)
: AliMUONVGeometryDESegmentation(),
fDetElemId(detElemId),
fPlaneType(bendingOrNonBending),
fSlat(0),
fSlatSegmentation(0),
fPadIterator(0),
fXhit(FMAX),
fYhit(FMAX)
{ 
  //
  // Normal ctor.
  //
	
	ReadMappingData();
		
  AliDebug(1,Form("this=%p detElemId=%3d %s fSlatSegmentation=%p",this,detElemId,
									( (bendingOrNonBending==kBendingPlane)?"Bending":"NonBending" ),
									fSlatSegmentation));
}

//_____________________________________________________________________________
AliMUONSt345SlatSegmentationV2::~AliMUONSt345SlatSegmentationV2()
{
	AliDebug(1,Form("dtor this=%p",this));
  delete fPadIterator;
}

//_____________________________________________________________________________
TF1*
AliMUONSt345SlatSegmentationV2::CorrFunc(Int_t) const
{
  AliFatal("Not Implemented");
  return 0x0;
}

//_____________________________________________________________________________
Float_t 
AliMUONSt345SlatSegmentationV2::Distance2AndOffset(Int_t,Int_t,
																									 Float_t,Float_t,Int_t*)
{
  AliFatal("Not Implemented");
  return 0.0;
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::Draw(Option_t*)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Float_t
AliMUONSt345SlatSegmentationV2::Dpx() const
{
  AliFatal("Not Implemented");
  return 0.0;
}

//_____________________________________________________________________________
Float_t
AliMUONSt345SlatSegmentationV2::Dpy() const
{
  AliFatal("Not Implemented");
  return 0.0;
}

//_____________________________________________________________________________
Float_t
AliMUONSt345SlatSegmentationV2::Dpx(int ipcb) const
{
	AliMpPCB* pcb = fSlat->GetPCB(ipcb);
	if (!pcb) 
  {
    AliFatal("pcb is null!");
  }
	return pcb->PadSizeX();
}

//_____________________________________________________________________________
Float_t
AliMUONSt345SlatSegmentationV2::Dpy(int ipcb) const
{
	AliMpPCB* pcb = fSlat->GetPCB(ipcb);
	if (!pcb) 
  {
    AliFatal("pcb is null!");
  }
	return pcb->PadSizeY();
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::FirstPad(Float_t xhit, Float_t yhit, 
																				 Float_t /*zhit*/, 
																				 Float_t dx, Float_t dy)
{
  // OK. We will do it in 2 steps. First get the area over which to
  // iterate, based on hit coordinates and (dx,dy). This first step
  // has nothing to do with segmentation in the first place, but with
  // how we simulate the fact that at some point the charge is shared
  // amongst several pads.
  // The second step is the actual pad iteration and is handled by 
  // a specific class (which has to do only with iteration...)
  //
  // FIXME: this method should not be here in the first place, IMHO.
  //
	
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
									fCurrentPad.GetIndices().GetFirst()+1,
									fCurrentPad.GetIndices().GetSecond()+1,fSlat->GetID()));		  
}

//_____________________________________________________________________________
Float_t
AliMUONSt345SlatSegmentationV2::GetAnod(Float_t xhit) const
{
  // Gets the x-coordinate of the wire which is the closest to xhit.
	
  Int_t n = Int_t(xhit/AliMUONConstants::Pitch());
  Float_t wire = (xhit>0) ? n+0.5 : n-0.5;
  return AliMUONConstants::Pitch()*wire;
}

//_____________________________________________________________________________
AliMUONGeometryDirection
AliMUONSt345SlatSegmentationV2::GetDirection()
{
  //AliWarning("Not Implemented");
  return kDirUndefined;
}

//______________________________________________________________________________
const AliMpVSegmentation*  
AliMUONSt345SlatSegmentationV2::GetMpSegmentation() const
{
// Returns the mapping segmentation
// (provides access to electronics info)

  return fSlatSegmentation;
}  


//_____________________________________________________________________________
void 
AliMUONSt345SlatSegmentationV2::GetNParallelAndOffset(Int_t,Int_t,Int_t*,Int_t*)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::GetPadC(Int_t ix, Int_t iy, 
																				Float_t& x, Float_t& y, Float_t& z)
{
  z = 0;
  GetPadC(ix,iy,x,y);
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::GetPadC(Int_t ix, Int_t iy, 
																				Float_t& x, Float_t& y)
{
  AliMpPad pad = 
  fSlatSegmentation->PadByIndices(AliMpIntPair(ix-fgIntOffset,iy-fgIntOffset),
                                  kTRUE);
  x = pad.Position().X();
  y = pad.Position().Y();
}


//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::GetPadI(Float_t x, Float_t y, Float_t /*z*/,
Int_t& ix, Int_t& iy)
{
  GetPadI(x,y,ix,iy);
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::GetPadI(Float_t x, Float_t y,
																				Int_t& ix, Int_t& iy)
{
//  Double_t slatx = fSlat->DX();
//  Double_t slaty = fSlat->DY();
  AliMpPad pad = 
    fSlatSegmentation->PadByPosition(TVector2(x,y), kTRUE);
//  fSlatSegmentation->PadByPosition(TVector2(x+slatx,y+slaty), kTRUE);
	
  if ( pad != AliMpPad::Invalid() )
	{
		ix = pad.GetIndices().GetFirst() + fgIntOffset;
		iy = pad.GetIndices().GetSecond() + fgIntOffset;
	}
  else
	{
		ix = iy = -1;
	}
}

//_____________________________________________________________________________
void 
AliMUONSt345SlatSegmentationV2::GiveTestPoints(Int_t&,Float_t*,Float_t*) const
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Bool_t
AliMUONSt345SlatSegmentationV2::HasPad(Float_t x, Float_t y, Float_t z)
{
  Int_t ix, iy;
  GetPadI(x,y,z,ix,iy);
  return HasPad(ix,iy);
}

//_____________________________________________________________________________
Bool_t
AliMUONSt345SlatSegmentationV2::HasPad(Int_t ix, Int_t iy)
{
  return fSlatSegmentation->HasPad(AliMpIntPair(ix-fgIntOffset,iy-fgIntOffset));
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::Init(int)
{
  AliWarning("Not Implemented because not needed ;-)");
}

//_____________________________________________________________________________
void  
AliMUONSt345SlatSegmentationV2::IntegrationLimits(Float_t& x1, Float_t& x2, 
																									Float_t& y1, Float_t& y2)
{
  //
  //  Returns integration limits for current pad
  //
	
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
AliMUONSt345SlatSegmentationV2::ISector()
{
  // FIXME: remove the usage of ISector from all the code.
  return -10;
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentationV2::Ix()
{
  if ( fPadIterator )
	{
		return fPadIterator->CurrentItem().GetIndices().GetFirst() + fgIntOffset;
	}
  else
	{
		return -1;
	}
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentationV2::Iy()
{
  if ( fPadIterator ) 
	{
		return fPadIterator->CurrentItem().GetIndices().GetSecond() + fgIntOffset;
	}
  else
	{
		return -1;
	}
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentationV2::MorePads()
{
  return (fPadIterator && !fPadIterator->IsDone());
}

//_____________________________________________________________________________
void 
AliMUONSt345SlatSegmentationV2::Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, 
																					 Int_t Xlist[10], Int_t Ylist[10])
{
  // Find pad at (ix,iy) for which we'll search neighbours.
  AliMpPad pad = 
	fSlatSegmentation->PadByIndices(AliMpIntPair(iX-fgIntOffset,iY-fgIntOffset),kTRUE);
	
  // Define the region to look into : a region slightly bigger
  // than the pad itself (10% bigger), in order to catch first neighbours.
//  AliMpArea area(pad.Position()+fSlat->Dimensions(),pad.Dimensions()*2.1); 
  AliMpArea area(pad.Position(),pad.Dimensions()*2.1); 
		
  AliMpVPadIterator* it = fSlatSegmentation->CreateIterator(area);
  it->First();
  Int_t n = 0;
  while ( !it->IsDone() && n < 10 )
	{
		AliMpPad p = it->CurrentItem();
		if ( p != pad ) // skip self
		{
			Xlist[n] = p.GetIndices().GetFirst() + fgIntOffset;
			Ylist[n] = p.GetIndices().GetSecond() + fgIntOffset;
			++n;
		}
		it->Next();
	}
	delete it;
  *Nlist = n;
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::NextPad()
{
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
AliMUONSt345SlatSegmentationV2::Npx() const
{
  return fSlatSegmentation->MaxPadIndexX()+1;
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentationV2::Npy() const
{
  return fSlatSegmentation->MaxPadIndexY()+1;
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::Print(Option_t*) const
{
  cout << "DetElemId=" << fDetElemId << " PlaneType=" 
  << fPlaneType << " Npx,Npy=" << Npx() << "," << Npy() << " fSlat=" << fSlat 
  << " fSlatSegmentation=" << fSlatSegmentation
  << " fSlatSegmentation->Slat()=" << fSlatSegmentation->Slat() << endl;
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::ReadMappingData()
{
	fSlatSegmentation = dynamic_cast<AliMpSlatSegmentation*>
  (AliMUONSegmentationManager::Segmentation(fDetElemId,fPlaneType));
	
  if (!fSlatSegmentation)
	{
		AliFatal("Wrong segmentation type encountered");
	}
  fSlat = fSlatSegmentation->Slat();
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentationV2::Sector(Int_t ix, Int_t)
{
  return fSlat->FindPCBIndex(ix - fgIntOffset);
}

//_____________________________________________________________________________
Int_t
AliMUONSt345SlatSegmentationV2::Sector(Float_t x, Float_t y)
{
  return fSlat->FindPCBIndex(x,y);
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::SetCorrFunc(Int_t,TF1*)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::SetDAnod(float)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::SetHit(Float_t x, Float_t y, Float_t)
{
  fXhit = x;
  fYhit = y;
	
  //
  // insure we're within the slat limits. If not, issue an error and sets
  // the current hit to slat center.
  // FIXME: this should probably a) not happen at all b) be a fatal error
  //
  if ( fXhit < -fSlat->DX() || fXhit > fSlat->DX() ||
       fYhit < -fSlat->DY() || fYhit > fSlat->DY() )
	{
		AliError(Form("Hit outside slat %s limits (x,y)hit = (%e,%e)."
                  " Forcing to (0,0)",fSlat->GetID(),fXhit,fYhit));
		fXhit = 0.0;
		fYhit = 0.0;
	}
  
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::SetPad(Int_t ix, Int_t iy)
{
  fCurrentPad = 
	fSlatSegmentation->PadByIndices(AliMpIntPair(ix-fgIntOffset,iy-fgIntOffset),kTRUE);
  if ( !fCurrentPad.IsValid() )
	{
		AliError(Form("Setting current pad to invalid ! (ix,iy)=(%4d,%4d)",ix,iy));
	}
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::SetPadSize(float,float)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Int_t 
AliMUONSt345SlatSegmentationV2::SigGenCond(Float_t,Float_t,Float_t)
{
  AliFatal("Not Implemented");
  return 0;
}

//_____________________________________________________________________________
void 
AliMUONSt345SlatSegmentationV2::SigGenInit(Float_t,Float_t,Float_t)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONSt345SlatSegmentationV2::Streamer(TBuffer &R__b)
{
  if (R__b.IsReading()) 
	{
    AliMUONSt345SlatSegmentationV2::Class()->ReadBuffer(R__b, this);
    ReadMappingData();
  } 
  else 
	{
    AliMUONSt345SlatSegmentationV2::Class()->WriteBuffer(R__b, this);
  }
}

