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
// $MpId: AliMpTriggerPadIterator.cxx,v 1.6 2006/05/24 13:58:50 ivana Exp $

#include "AliMpTriggerPadIterator.h"

#include "TArrayI.h"
#include "AliLog.h"
#include "AliMpTrigger.h"
#include "AliMpSlat.h"
#include "AliMpVSegmentation.h"


//-----------------------------------------------------------------------------
/// \class AliMpTriggerPadIterator
///
/// Implementation of AliMpVPadIterator for trigger slats.
///
/// The class iterate over the pads in a trigger slats
///
/// \author D. Stocco
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMpTriggerPadIterator)
/// \endcond

//_____________________________________________________________________________
AliMpTriggerPadIterator::AliMpTriggerPadIterator()
: AliMpVPadIterator(),
fkTriggerSlat(0),
fPadList(),
fCurrentPadIndex(0)
{
  ///
  /// Empty (default) ctor.
  ///
}

//_____________________________________________________________________________
AliMpTriggerPadIterator::AliMpTriggerPadIterator(const AliMpTrigger* slat )
: AliMpVPadIterator(),
fkTriggerSlat(slat),
fPadList(),
fCurrentPadIndex(0)
{
  ///
  /// Normal ctor.
  /// The iteration will occur on the given slat
  ///
  Prepare();
  fPadList.SetOwner();
}

//_____________________________________________________________________________
AliMpTriggerPadIterator::~AliMpTriggerPadIterator()
{ 
  ///
  /// Dtor.
  ///
  AliDebug(1,Form("this=%p dtor",this));
  Invalidate();
}

//_____________________________________________________________________________
Bool_t
AliMpTriggerPadIterator::Prepare()
{
  ///
  /// Allocate the corresponding delegate iterators.

  TArrayI boards;
  fkTriggerSlat->GetLayer(0)->GetAllMotifPositionsIDs(boards);
  AliMpVSegmentation* seg = fkTriggerSlat->GetLayerSegmentation(0);
  for ( Int_t iboard=0; iboard<boards.GetSize(); iboard++ ) {
    for (Int_t ibitxy=0; ibitxy<16; ++ibitxy) {
      AliMpPad pad = seg->PadByLocation(boards[iboard],ibitxy,kFALSE);
      if ( ! pad.IsValid() ) continue;
      fPadList.Add(pad.Clone());
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
AliMpPad
AliMpTriggerPadIterator::CurrentItem() const
{
  ///
  /// Returns the current pad of the iteration.
  ///
  if ( fCurrentPadIndex == fPadList.GetEntriesFast() ) return AliMpPad::Invalid();
  return *(static_cast<AliMpPad*>(fPadList.At(fCurrentPadIndex)));
}

//_____________________________________________________________________________
void
AliMpTriggerPadIterator::First()
{
  ///
  /// (Re)starts the iteration.
  ///
  if ( fPadList.GetEntriesFast() == 0 )
	{
		AliError("Iterator is not valid, as it gets no delegates at all !");
	}
  else
	{
		fCurrentPadIndex = 0;
	}
}

//_____________________________________________________________________________
void
AliMpTriggerPadIterator::Invalidate()
{
  ///
  /// Make the iterator invalid.
  ///
  fPadList.Delete();
  fCurrentPadIndex = 0;
}

//_____________________________________________________________________________
Bool_t
AliMpTriggerPadIterator::IsDone() const
{
  ///
  /// Returns whether the iteration is ended or not.
  ///
  return ( fCurrentPadIndex >= fPadList.GetEntriesFast() );
}

//_____________________________________________________________________________
void
AliMpTriggerPadIterator::Next()
{
  ///
  /// Next step of the iteration.
  ///
  if (IsDone()) return;

  fCurrentPadIndex++;
}
