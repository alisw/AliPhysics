// $Id$
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

ClassImp(AliMpSectorPadIterator)

//______________________________________________________________________________
AliMpSectorPadIterator::AliMpSectorPadIterator()
  : AliMpVPadIterator(),
    fkSector(0),
    fCurrentRow(0),
    fCurrentSeg(0),
    fCurrentMotif(0),
    fMotifPos(0),
    fIterator()
{
// default constructor, set the current position to "invalid"
}

//______________________________________________________________________________
AliMpSectorPadIterator::AliMpSectorPadIterator(const AliMpSector* const sector)
  : AliMpVPadIterator(),
    fkSector(sector),
    fCurrentRow(0),
    fCurrentSeg(0),
    fCurrentMotif(0),
    fMotifPos(0),
    fIterator()
{
// normal constructor, set *this to invalid position  
}

//______________________________________________________________________________
AliMpSectorPadIterator::AliMpSectorPadIterator(const AliMpSectorPadIterator& right)
  : AliMpVPadIterator(right)
{
// copy constructor
 
  *this = right;
}

//______________________________________________________________________________
AliMpSectorPadIterator::~AliMpSectorPadIterator()
{
// destructor
}

// operators

//______________________________________________________________________________
AliMpSectorPadIterator& 
AliMpSectorPadIterator::operator = (const AliMpSectorPadIterator& right)
{
// assignement operator

  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  AliMpVPadIterator::operator=(right);

  fkSector      = right.fkSector;
  fCurrentRow   = right.fCurrentRow;
  fCurrentSeg   = right.fCurrentSeg;
  fCurrentMotif = right.fCurrentMotif;
  fIterator     = right.fIterator;

  return *this;
} 

//private methods

//______________________________________________________________________________
AliMpMotifPosition* AliMpSectorPadIterator::ResetToCurrentMotifPosition()
{
  // Find the AliMpMotifType object associated with the triplet
  // (fCurrentRow, fCurrentSeg, fCurrentMotif)
  // place it in the private fMotifType member and returns it.
  
  fMotifPos =0;
  
  if (fkSector){
    AliMpRow* row;
    if ((fCurrentRow >= 0) && (fCurrentRow < fkSector->GetNofRows())){
      row= fkSector->GetRow(fCurrentRow);

      AliMpVRowSegment* seg;
      if (fCurrentSeg<row->GetNofRowSegments()){
        seg = row->GetRowSegment(fCurrentSeg);

        if (fCurrentMotif<seg->GetNofMotifs()){
          fMotifPos = 
           fkSector->GetMotifMap()->FindMotifPosition(
                seg->GetMotifPositionId(fCurrentMotif));
        }
      }
    }
  }
  
  if (fMotifPos) {
    fIterator = AliMpMotifPositionPadIterator(fMotifPos);
    fIterator.First();
  }
  else
    Invalidate();

  return fMotifPos;
}

//______________________________________________________________________________
Bool_t AliMpSectorPadIterator::IsValid() const
{
// Is the iterator in a valid position?
    return (fkSector!=0) && (fMotifPos!=0);
} 

//public methods

//______________________________________________________________________________
void AliMpSectorPadIterator::First()
{
// Reset the iterator, so that it points to the first available
// pad in the sector

    if (!fkSector) {
        Invalidate();
        return;
    }
    fCurrentRow =0;
    fCurrentSeg=0;
    fCurrentMotif=0;
    
    ResetToCurrentMotifPosition();

    return;
}

//______________________________________________________________________________
void AliMpSectorPadIterator::Next()
{
// Move the iterator to the next valid pad.


  //if (!IsValid()) return *this;
  if (!IsValid()) return;

  fIterator.Next();
  
  if (!fIterator.IsDone()) return;
  

  // Go to ne next motif, in the current segment
  ++fCurrentMotif;
  if (ResetToCurrentMotifPosition()) return;


  // if motif number is too big, set it to 0 and pass to the next row segment
  fCurrentMotif=0;
  ++fCurrentSeg;
  if (ResetToCurrentMotifPosition()) return;


  // if row segment number is too big, pass to the next row
  fCurrentSeg=0;
  ++fCurrentRow;
  if (ResetToCurrentMotifPosition()) return;
  
  // if row number is too big, the invalidate the iterator (==End())
  Invalidate();
  return;

}

//______________________________________________________________________________
Bool_t AliMpSectorPadIterator::IsDone() const
{
// 
  return !IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpSectorPadIterator::CurrentItem () const 
{
// Returns current pad.

  if (!IsValid())
    return AliMpPad::Invalid();
      

  // no more verification, since IsValid() is TRUE here.

  return fIterator.CurrentItem();
}

//______________________________________________________________________________
void AliMpSectorPadIterator::Invalidate()
{
// Let the iterator points to the invalid position
    fMotifPos = 0;
    fIterator.Invalidate();
} 

