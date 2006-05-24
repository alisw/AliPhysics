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
// $MpId: AliMpMotifPositionPadIterator.cxx,v 1.6 2006/05/24 13:58:41 ivana Exp $
// Category: motif
//
// Class AliMpMotifPositionPadIterator
// -----------------------------------
// Class, which defines an iterator over the pads of a given motif type
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpMotifPositionPadIterator.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotifType.h"
#include "AliMpConnection.h"

/// \cond CLASSIMP
ClassImp(AliMpMotifPositionPadIterator)
/// \endcond

//______________________________________________________________________________
AliMpMotifPositionPadIterator::AliMpMotifPositionPadIterator():
    AliMpVPadIterator(),
    fMotifPos(0),
    fIterator()
{
/// Default constructor, set the current position to "invalid"
}

//______________________________________________________________________________

AliMpMotifPositionPadIterator::AliMpMotifPositionPadIterator(
                                    const AliMpMotifPosition* motifPos)
  : AliMpVPadIterator(),
    fMotifPos(motifPos),
    fIterator(motifPos->GetMotif()->GetMotifType())
{
/// Standard constructor, let *this to invalid position
}

//______________________________________________________________________________
AliMpMotifPositionPadIterator::AliMpMotifPositionPadIterator(
                                    const AliMpMotifPositionPadIterator& right)
  : AliMpVPadIterator(right),
    fMotifPos(right.fMotifPos),
    fIterator(right.fIterator)
    
{
/// Copy constructor
}

//______________________________________________________________________________
AliMpMotifPositionPadIterator::~AliMpMotifPositionPadIterator()
{
/// Destructor
}

// operators

//______________________________________________________________________________
AliMpMotifPositionPadIterator& 
AliMpMotifPositionPadIterator::operator = (const AliMpMotifPositionPadIterator& right)
{
/// Assignment operator

// if the right hand iterator isn't of good type
// the current operator is invalidated

  // check assignment to self
  if (this == &right) return *this;

  // base class assignment
  AliMpVPadIterator::operator=(right);

  fMotifPos = right.fMotifPos;
  fIterator = right.fIterator;

  return *this;
}  

//private methods


//______________________________________________________________________________
Bool_t AliMpMotifPositionPadIterator::IsValid() const
{
/// Is the iterator in a valid position?

    return (fMotifPos!=0) && (!fIterator.IsDone());
} 

//
// public methods
//

//______________________________________________________________________________
void AliMpMotifPositionPadIterator::First()
{
/// Reset the iterator, so that it points to the first available
/// pad in the motif type

    if (!fMotifPos) {
        Invalidate();
        return ;
    }

    fIterator.First();
    return;
}

//______________________________________________________________________________
void AliMpMotifPositionPadIterator::Next()
{
/// Move the iterator to the next valid pad.

  fIterator.Next();
}

//______________________________________________________________________________
Bool_t AliMpMotifPositionPadIterator::IsDone() const
{
/// Is the iterator in the end? 

  return !IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpMotifPositionPadIterator::CurrentItem() const 
{
/// Return current pad.

    if (!fMotifPos)
        return AliMpPad::Invalid();
    else {
      AliMpIntPair ind = fIterator.CurrentItem().GetIndices();
      AliMpMotifType* mt = fMotifPos->GetMotif()->GetMotifType();
      AliMpConnection* connect = 
        mt->FindConnectionByLocalIndices(ind);
      return AliMpPad(AliMpIntPair(fMotifPos->GetID(),connect->GetGassiNum()),
                  fMotifPos->GlobalIndices(ind),
                  fMotifPos->Position()+fMotifPos->GetMotif()->PadPositionLocal(ind),
                  fMotifPos->GetMotif()->GetPadDimensions(ind));
    }
}

//______________________________________________________________________________
void AliMpMotifPositionPadIterator::Invalidate()
{
/// Let the iterator point to the invalid position

  fIterator.Invalidate();
} 

