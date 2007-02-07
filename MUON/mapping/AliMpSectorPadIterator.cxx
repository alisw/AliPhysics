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
// $MpId: AliMpSectorPadIterator.cxx,v 1.6 2006/05/24 13:58:46 ivana Exp $
// Category: sector
//
// Class AliMpSectorPadIterator
// ----------------------------
// Class, which defines an iterator over the pads of a sector
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpSectorPadIterator.h"
#include "AliMpIntPair.h"
#include "AliMpSector.h"
#include "AliMpMotifType.h"

#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"

/// \cond CLASSIMP
ClassImp(AliMpSectorPadIterator)
/// \endcond

//______________________________________________________________________________
AliMpSectorPadIterator::AliMpSectorPadIterator()
  : AliMpVPadIterator(),
    fkSector(0),
    fCurrentIndex(0),
    fMotifPos(0),
    fIterator()
{
/// Default constructor, set the current position to "invalid"
}

//______________________________________________________________________________
AliMpSectorPadIterator::AliMpSectorPadIterator(const AliMpSector* const sector)
  : AliMpVPadIterator(),
    fkSector(sector),
    fCurrentIndex(0),
    fMotifPos(0),
    fIterator()
{
/// Standard constructor, set *this to invalid position  
}

//______________________________________________________________________________
AliMpSectorPadIterator::AliMpSectorPadIterator(const AliMpSectorPadIterator& right)
  : AliMpVPadIterator(right),
    fkSector(0),
    fCurrentIndex(0),
    fMotifPos(0),
    fIterator()
{
/// Copy constructor
 
  *this = right;
}

//______________________________________________________________________________
AliMpSectorPadIterator::~AliMpSectorPadIterator()
{
/// Destructor
}

//
// operators
//

//______________________________________________________________________________
AliMpSectorPadIterator& 
AliMpSectorPadIterator::operator = (const AliMpSectorPadIterator& right)
{
/// Assignment operator

  // check assignment to self
  if (this == &right) return *this;

  // base class assignment
  AliMpVPadIterator::operator=(right);

  fkSector      = right.fkSector;
  fCurrentIndex = right.fCurrentIndex,
  fMotifPos     = right.fMotifPos;
  fIterator     = right.fIterator;

  return *this;
} 

//private methods

//______________________________________________________________________________
AliMpMotifPosition* AliMpSectorPadIterator::ResetToCurrentMotifPosition()
{
/// Find the AliMpMotifType object associated with the triplet
/// (fCurrentRow, fCurrentSeg, fCurrentMotif),
/// place it in the private fMotifType member and return it.

  if ( fCurrentIndex == fkSector->GetMotifMap()->GetNofMotifPositions() ) {
    Invalidate();
    return 0;
  }  
    
  fMotifPos = fkSector->GetMotifMap()->GetMotifPosition(fCurrentIndex);
  fIterator = AliMpMotifPositionPadIterator(fMotifPos);
  fIterator.First();

  return fMotifPos;
}

//______________________________________________________________________________
Bool_t AliMpSectorPadIterator::IsValid() const
{
/// Is the iterator in a valid position?

    return (fkSector!=0) && (fMotifPos!=0);
} 

//
//public methods
//

//______________________________________________________________________________
void AliMpSectorPadIterator::First()
{
/// Reset the iterator, so that it points to the first available
/// pad in the sector

    if (!fkSector) {
        Invalidate();
        return;
    }
    fCurrentIndex =0;
    ResetToCurrentMotifPosition();

    return;
}

//______________________________________________________________________________
void AliMpSectorPadIterator::Next()
{
/// Move the iterator to the next valid pad.

  if (!IsValid()) return;

  fIterator.Next();
  
  if (!fIterator.IsDone()) return;
  

  // Go to the next motif, in the current segment
  ++fCurrentIndex;
  if (ResetToCurrentMotifPosition()) return;

  Invalidate();
  return;

}

//______________________________________________________________________________
Bool_t AliMpSectorPadIterator::IsDone() const
{
/// Is the iterator in the end? 

  return ! IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpSectorPadIterator::CurrentItem () const 
{
/// Return current pad.

  if (!IsValid())
    return AliMpPad::Invalid();
      

  // no more verification, since IsValid() is TRUE here.
  return fIterator.CurrentItem();
}

//______________________________________________________________________________
void AliMpSectorPadIterator::Invalidate()
{
/// Let the iterator point to the invalid position
    fMotifPos = 0;
    fIterator.Invalidate();
} 

