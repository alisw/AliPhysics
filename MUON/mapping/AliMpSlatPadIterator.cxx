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

// $Id$
// $MpId: AliMpSlatPadIterator.cxx,v 1.6 2006/05/24 13:58:50 ivana Exp $

#include "AliMpSlatPadIterator.h"

#include "AliLog.h"
#include "AliMpArea.h"
#include "AliMpPCB.h"
#include "AliMpSlat.h"
#include "AliMpSlatZonePadIterator.h"

#include <algorithm>
#include <limits>
#include <cassert>

///
/// \class AliMpSlatPadIterator
///
/// Implementation of AliMpVPadIterator for slats.
///
/// This class first split the input area (upon which to iterate)
/// into a set of areas of constant pad size.
/// Then each of those areas is iterated over, using
/// AliMpSlatZonePadIterator objects.
///
/// \author L. Aphecetche
///

/// \cond CLASSIMP
ClassImp(AliMpSlatPadIterator)
/// \endcond

//const Double_t
//AliMpSlatPadIterator::fgkDmax = std::numeric_limits<Double_t>::max();

//_____________________________________________________________________________
AliMpSlatPadIterator::AliMpSlatPadIterator()
: AliMpVPadIterator(),
fkSlat(0),
fDelegates(),
fCurrentDelegate(0),
fCurrentDelegateIndex(0)
{
  //
  // Empty (default) ctor.
  //
}

//_____________________________________________________________________________
AliMpSlatPadIterator::AliMpSlatPadIterator(const AliMpSlat* slat,
																					 const AliMpArea& area)
: AliMpVPadIterator(),
fkSlat(slat),
fDelegates(),
fCurrentDelegate(0),
fCurrentDelegateIndex(0)
{
  //
  // Normal ctor.
  // The iteration will occur on the given slat over the specified area.
  //
  AliDebug(1,Form("this=%p ctor",this));
  if (!Prepare(area)) 
	{
		AliError("Iterator invalidated by improper initialization (e.g. incorrect area given ?)");
	}
}

//_____________________________________________________________________________
AliMpSlatPadIterator::~AliMpSlatPadIterator()
{ 
  //
  // Dtor.
  //
  AliDebug(1,Form("this=%p dtor",this));
  Invalidate();
}

//_____________________________________________________________________________
AliMpArea
AliMpSlatPadIterator::Intersect(const AliMpArea& a, const AliMpArea& b) const
{ 
  //
  // Returns the common part of a and b.
  //
  AliDebug(4,Form("a=(%7.2f,%7.2f;%7.2f,%7.2f) b=(%7.2f,%7.2f;%7.2f,%7.2f)",
									a.LeftBorder(),a.DownBorder(),a.RightBorder(),a.UpBorder(),
									b.LeftBorder(),b.DownBorder(),b.RightBorder(),b.UpBorder()));
	
	Double_t xmin = std::max(a.LeftBorder(),b.LeftBorder());
  Double_t xmax = std::min(a.RightBorder(),b.RightBorder());
  Double_t ymin = std::max(a.DownBorder(),b.DownBorder());
  Double_t ymax = std::min(a.UpBorder(),b.UpBorder());
  AliMpArea c( TVector2( (xmin+xmax)/2.0, (ymin+ymax)/2.0 ),
							 TVector2( (xmax-xmin)/2.0, (ymax-ymin)/2.0 ) );
	
  AliDebug(4,Form("a intersect b = (%7.2f,%7.2f;%7.2f,%7.2f)",
									c.LeftBorder(),c.DownBorder(),c.RightBorder(),c.UpBorder()));
  return c;
}

//_____________________________________________________________________________
Bool_t
AliMpSlatPadIterator::Prepare(const AliMpArea& area)
{
  //
  // Split area into smaller area intersecting pcbs,
  // and allocate the corresponding delegate iterators.
	
  for ( AliMpSlat::Size_t i = 0; i < fkSlat->GetSize(); ++i )
	{
		const AliMpPCB* pcb = fkSlat->GetPCB(i);
		AliMpArea pcbArea( TVector2( (pcb->Xmin()+pcb->Xmax())/2.0,fkSlat->DY()),
											 TVector2( pcb->DX(), pcb->DY() ) );
		AliMpArea zone = Intersect(pcbArea,area);
		AliDebug(3,Form("i=%2d zone is %7.2f,%7.2f->%7.2f,%7.2f %d",i,
										zone.LeftBorder(),zone.DownBorder(),
										zone.RightBorder(),zone.UpBorder(),
										zone.IsValid()));
		if ( zone.IsValid() )
		{
			fDelegates.push_back(new AliMpSlatZonePadIterator(fkSlat,zone));
		}
	}
  return !fDelegates.empty();
}

//_____________________________________________________________________________
AliMpPad
AliMpSlatPadIterator::CurrentItem() const
{
  //
  // Returns the current pad of the iteration.
  //
  if ( fCurrentDelegate )
	{
		return fCurrentDelegate->CurrentItem();
	}
  else
	{
		return AliMpPad::Invalid();
	}
}

//_____________________________________________________________________________
void
AliMpSlatPadIterator::First()
{
  //
  // (Re)starts the iteration.
  //
  if ( fDelegates.empty() )
	{
		AliError("Iterator is not valid, as it gets no delegates at all !");
	}
  else
	{
		fCurrentDelegateIndex = 0;
		fCurrentDelegate = fDelegates[0];
		fCurrentDelegate->First();
	}
}

//_____________________________________________________________________________
void
AliMpSlatPadIterator::Invalidate()
{
  //
  // Make the iterator invalid.
  //
  for ( size_t i = 0; i < fDelegates.size(); ++i )
	{
		delete fDelegates[i];
		fDelegates[i] = 0;
	}
  fDelegates.clear();
  fCurrentDelegate = 0;
  fCurrentDelegateIndex = 0;
}

//_____________________________________________________________________________
Bool_t
AliMpSlatPadIterator::IsDone() const
{
  //
  // Returns whether the iteration is ended or not.
  //
  return ( !fCurrentDelegate ||
					 ( fCurrentDelegateIndex >= fDelegates.size() && 
						 fCurrentDelegate->IsDone() ) );
}

//_____________________________________________________________________________
void
AliMpSlatPadIterator::Next()
{
  //
  // Next step of the iteration.
  //
  if (IsDone()) return;
	
  fCurrentDelegate->Next();
	
  if ( fCurrentDelegate->IsDone() )
	{
		AliDebug(3,"Moving to next delegate");
		++fCurrentDelegateIndex;
		if ( fCurrentDelegateIndex < fDelegates.size() )
		{
			fCurrentDelegate = fDelegates[fCurrentDelegateIndex];
			fCurrentDelegate->First();
		}
	}
}
