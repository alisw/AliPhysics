// $Id$
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

ClassImp(AliMpMotifPositionPadIterator)

//______________________________________________________________________________
AliMpMotifPositionPadIterator::AliMpMotifPositionPadIterator():
    AliMpVPadIterator(),
    fMotifPos(0),
    fIterator()
{
// default constructor, set the current position to "invalid"
}

//______________________________________________________________________________

AliMpMotifPositionPadIterator::AliMpMotifPositionPadIterator(
                                    const AliMpMotifPosition* motifPos)
  : AliMpVPadIterator(),
    fMotifPos(motifPos),
    fIterator(motifPos->GetMotif()->GetMotifType())
{
// normal constructor, let *this to invalid position
}

//______________________________________________________________________________
AliMpMotifPositionPadIterator::AliMpMotifPositionPadIterator(
                                    const AliMpMotifPositionPadIterator& right)
  : AliMpVPadIterator(right),
    fMotifPos(right.fMotifPos),
    fIterator(right.fIterator)
    
{
// copy constructor
}

//______________________________________________________________________________
AliMpMotifPositionPadIterator::~AliMpMotifPositionPadIterator()
{
// destructor
}

// operators

//______________________________________________________________________________
AliMpMotifPositionPadIterator& 
AliMpMotifPositionPadIterator::operator = (const AliMpMotifPositionPadIterator& right)
{
// assignement operator
// if the right hand iterator isn't of good type
// the current operator is invalidated

  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  AliMpVPadIterator::operator=(right);

  fMotifPos = right.fMotifPos;
  fIterator = right.fIterator;

  return *this;
}  

//private methods


//______________________________________________________________________________
Bool_t AliMpMotifPositionPadIterator::IsValid() const
{
// Is the iterator in a valid position?

    return (fMotifPos!=0) && (!fIterator.IsDone());
} 

//public methods

//______________________________________________________________________________
void AliMpMotifPositionPadIterator::First()
{
// Reset the iterator, so that it points to the first available
// pad in the motif type

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
// Move the iterator to the next valid pad.
  fIterator.Next();
}

//______________________________________________________________________________
Bool_t AliMpMotifPositionPadIterator::IsDone() const
{
// 
  return !IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpMotifPositionPadIterator::CurrentItem() const 
{
// Returns current pad.

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
// Let the iterator points to the invalid position
  fIterator.Invalidate();
} 

