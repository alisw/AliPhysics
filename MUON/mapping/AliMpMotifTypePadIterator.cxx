// $Id$
// Category: motif
//
// Class AliMpMotifTypePadIterator
// -------------------------------
// Class, which defines an iterator over the pads of a given motif type
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpMotifTypePadIterator.h"
#include "AliMpMotifType.h"

ClassImp(AliMpMotifTypePadIterator)

//______________________________________________________________________________
AliMpMotifTypePadIterator::AliMpMotifTypePadIterator():
    AliMpVPadIterator(),
    fMotifType(0),
    fCurrentPosition(AliMpIntPair::Invalid())
{
// default constructor, set the current position to "invalid"
}

//______________________________________________________________________________
AliMpMotifTypePadIterator::AliMpMotifTypePadIterator( 
                                const AliMpMotifType* motifType)
  : AliMpVPadIterator(),
    fMotifType(motifType),
    fCurrentPosition(AliMpIntPair::Invalid())
{
// normal constructor, let *this to invalid position
}

//______________________________________________________________________________
AliMpMotifTypePadIterator::AliMpMotifTypePadIterator(
                                const AliMpMotifTypePadIterator& right)
  : AliMpVPadIterator(right),
    fMotifType(right.fMotifType),
    fCurrentPosition(right.fCurrentPosition)
    
{
// copy constructor
}

//______________________________________________________________________________
AliMpMotifTypePadIterator::~AliMpMotifTypePadIterator()
{
// destructor
}

// operators

//______________________________________________________________________________
AliMpMotifTypePadIterator& 
AliMpMotifTypePadIterator::operator = (const AliMpMotifTypePadIterator& right)
{
// assignement operator
// if the right hand iterator isn't of good type
// the current operator is invalidated

  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  AliMpVPadIterator::operator=(right);

  fMotifType = right.fMotifType;
  fCurrentPosition = right.fCurrentPosition;

  return *this;
}  

//
//private methods
//

//______________________________________________________________________________
AliMpIntPair 
AliMpMotifTypePadIterator::FindFirstPadInLine(AliMpIntPair indices) const
{
// Find the indices of the first pad in the same line
// as the <indices>, and in column, at least equal, to the
// one of <indices>

    if (!fMotifType) return AliMpIntPair::Invalid();

    while (indices.GetFirst() < fMotifType->GetNofPadsX()) {
        if (fMotifType->HasPad(indices)) return indices;
        indices += AliMpIntPair(1,0);
    }
    return AliMpIntPair::Invalid();
} 

//______________________________________________________________________________
Bool_t AliMpMotifTypePadIterator::IsValid() const
{
// Is the iterator in a valid position?

    return fMotifType!=0 && fCurrentPosition.IsValid();
} 

//
//public methods
//

//______________________________________________________________________________
void AliMpMotifTypePadIterator::First()
{
// Reset the iterator, so that it points to the first available
// pad in the motif type

    if (!fMotifType) {
        Invalidate();
        return ;
    }
    fCurrentPosition = AliMpIntPair(0,0);
    if (fMotifType->HasPad(fCurrentPosition)) return;
    
    
    // if (0,0) is not available
    // place itself to the first avalable motif after (0,0) (if exists)
    // ++(*this);
    Next();
    
    return;
}

//______________________________________________________________________________
void AliMpMotifTypePadIterator::Next()
{
// Move the iterator to the next valid pad.

    //if (!IsValid()) return *this;
    if (!IsValid()) return;

    while (fCurrentPosition.GetSecond() < fMotifType->GetNofPadsY()){
        AliMpIntPair nextTry 
	  = FindFirstPadInLine(fCurrentPosition + AliMpIntPair(1,0));
        if (nextTry.IsValid()){
            fCurrentPosition = nextTry;
            return;
        }
        fCurrentPosition.SetFirst(-1);
        fCurrentPosition.SetSecond(fCurrentPosition.GetSecond()+1);
    }
    
    // if the loop is finished, there's not available pads at all...
    Invalidate();
    return;
}

//______________________________________________________________________________
Bool_t AliMpMotifTypePadIterator::IsDone() const
{
// 
  return !IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpMotifTypePadIterator::CurrentItem() const 
{
// Returns current pad.

    if (!fMotifType)
        return AliMpPad::Invalid();
    else
        return AliMpPad(AliMpIntPair::Invalid(),
	                fCurrentPosition,TVector2(),TVector2());
}

//______________________________________________________________________________
void AliMpMotifTypePadIterator::Invalidate()
{
// Let the iterator points to the invalid position
    fCurrentPosition = AliMpIntPair::Invalid();

} 

