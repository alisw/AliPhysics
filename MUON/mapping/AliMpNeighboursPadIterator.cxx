// $Id$
// Category: sector
//
// Class AliMpNeighboursPadIterator
// --------------------------------
// Class, which defines an iterator over the pads surrounding a given pad
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <set>

#include <TVector2.h>

#include "AliMpNeighboursPadIterator.h"
#include "AliMpIntPair.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpRow.h"
#include "AliMpConstants.h"

ClassImp(AliMpNeighboursPadIterator)

const UInt_t AliMpNeighboursPadIterator::fgkInvalidIndex = 9999; 
                                                   //never so much neighbours...

//______________________________________________________________________________
AliMpNeighboursPadIterator::AliMpNeighboursPadIterator()
  : AliMpVPadIterator(),
    fkSegmentation(0),
    fCenterPad(AliMpPad::Invalid()),
    fPads(),
    fIndex(fgkInvalidIndex)
{
// default constructor, set the current position to "invalid"
}

//______________________________________________________________________________
AliMpNeighboursPadIterator::AliMpNeighboursPadIterator(
                                 const AliMpSectorSegmentation* segmentation,
                                 const AliMpPad& centerPad,
                                 Bool_t includeCenter)
  : AliMpVPadIterator(),
    fkSegmentation(segmentation),
    fCenterPad(centerPad),
    fIndex(fgkInvalidIndex)
{
// normal constructor, set *this to invalid position

    FillPadsVector(includeCenter);
}

//______________________________________________________________________________
AliMpNeighboursPadIterator::AliMpNeighboursPadIterator(
                                 const AliMpNeighboursPadIterator& right)
  : AliMpVPadIterator(right)
{
// copy constructor

  *this = right;
}

//______________________________________________________________________________
AliMpNeighboursPadIterator::~AliMpNeighboursPadIterator()
{
// destructor
}

// operators

//______________________________________________________________________________
AliMpNeighboursPadIterator& 
AliMpNeighboursPadIterator::operator = (const AliMpNeighboursPadIterator& right)
{
// assignement operator
// if the right hand iterator isn't of good type
// the current operator is invalidated

  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  AliMpVPadIterator::operator=(right);

  fkSegmentation = right.fkSegmentation;
  fCenterPad     = right.fCenterPad;
  fPads          = right.fPads;
  fIndex         = right.fIndex;
    

  return *this;
} 

//private methods

//______________________________________________________________________________
Bool_t AliMpNeighboursPadIterator::IsNeighbours(const AliMpPad& pad) const
{
// true if the pad located by <padIndice> is a neighbours of those
// located at <fCenterPad>

    
    TVector2 relPos  = pad.Position()   - fCenterPad.Position();
    TVector2 bounds  = pad.Dimensions() + fCenterPad.Dimensions();
    return (TMath::Abs(relPos.X())- bounds.X()<AliMpConstants::LengthTolerance()) && 
           (TMath::Abs(relPos.Y())- bounds.Y()<AliMpConstants::LengthTolerance());

}

//______________________________________________________________________________
PadVector AliMpNeighboursPadIterator::PadVectorLine(const AliMpPad& from,
                                           const AliMpIntPair& direction) const
{
// Fill  a new vector with all pads which have common
// parts with the pad located at <fCenterPad>, in a given line
// starting from <from> and moving by <direction>

    AliMpPad current = from;
    PadVector ans;
    Bool_t cont=kTRUE;
    do {
        if (IsNeighbours(current))
            ans.push_back(current);
        else
            cont=kFALSE;
        TVector2 nextPos = current.Position() + TVector2(
          current.Dimensions().X()*(AliMpConstants::LengthStep()+1.)*direction.GetFirst(),
          current.Dimensions().Y()*(AliMpConstants::LengthStep()+1.)*direction.GetSecond());
        current = fkSegmentation->PadByPosition(nextPos);
    } while (cont);
    return ans;
}

//______________________________________________________________________________
void AliMpNeighboursPadIterator::FillPadsVector(Bool_t includeCenter)
{
// Fill the indices vector with all indices of pads which have common
// parts with the pad located at <fCenterPad>

    if (!fkSegmentation || !fCenterPad.IsValid()) return;
    
    
    AliMpPad from;
    AliMpIntPair direction;
    PadVector found;
    
    // repare a unique simple associative container
    // --> no doublons, rapid insersion
    PadSet setTotal;

  /////////////  Left side
  
  ////////////////// up direction
    
    from = fkSegmentation->PadsLeft(fCenterPad).GetFirst();
    direction = AliMpIntPair(0,1);
    found = PadVectorLine(from,direction);
    setTotal.insert(found.begin(),found.end());


  ////////////////// down direction

    from = fkSegmentation->PadsDown(from).GetFirst(); // the Pad down is already added
    direction = AliMpIntPair(0,-1);
    found = PadVectorLine(from,direction);
    setTotal.insert(found.begin(),found.end());
    
  /////////////  Up side
  
  ////////////////// right direction

    from = fkSegmentation->PadsUp(fCenterPad).GetFirst();
    direction = AliMpIntPair(1,0);
    found = PadVectorLine(from,direction);
    setTotal.insert(found.begin(),found.end());
    
  ////////////////// left direction

    from = fkSegmentation->PadsLeft(from).GetFirst(); // the pad up is already added
    direction = AliMpIntPair(-1,0);
    found = PadVectorLine(from,direction);
    setTotal.insert(found.begin(),found.end());
    
  /////////////  Right side
  
  ////////////////// Up direction
    
    from = fkSegmentation->PadsRight(fCenterPad).GetFirst();
    direction = AliMpIntPair(0,1);
    found = PadVectorLine(from,direction);
    setTotal.insert(found.begin(),found.end());
    
  ////////////////// down direction

    from = fkSegmentation->PadsDown(from).GetFirst(); // the pad right is already added
    direction = AliMpIntPair(0,-1);
    found = PadVectorLine(from,direction);
    setTotal.insert(found.begin(),found.end());
    
  /////////////  Down side
  
  ////////////////// Right direction

    from = fkSegmentation->PadsDown(fCenterPad).GetFirst();
    direction = AliMpIntPair(1,0);
    found = PadVectorLine(from,direction);
    setTotal.insert(found.begin(),found.end());
    
  ////////////////// left direction
    
    from = fkSegmentation->PadsLeft(from).GetFirst(); // the pad down is already added
    direction = AliMpIntPair(-1,0);
    found = PadVectorLine(from,direction);
    setTotal.insert(found.begin(),found.end());
    

    // fill the fIndices vector with the set (-->pass from a rapid insertion,
    // to rapid and indexed access, for the rest of the job)

    fPads.clear();
    // include the center pad if requiered
    if (includeCenter) fPads.push_back(fCenterPad);
    //fPads.insert(fPads.end(),setTotal.begin(),setTotal.end());
    
    PadSetIterator it;
    for (it = setTotal.begin(); it != setTotal.end(); it++)
      fPads.push_back((*it));
}

//______________________________________________________________________________
Bool_t AliMpNeighboursPadIterator::IsValid() const
{
// Is the iterator in a valid position?
    return (fkSegmentation!=0 && fIndex!=fgkInvalidIndex);
} 

//public methods

//______________________________________________________________________________
void AliMpNeighboursPadIterator::First()
{
// Reset the iterator, so that it points to the first available
// pad in the sector

    if ((fkSegmentation != 0) && (fPads.size() != 0)) 
      fIndex=0; 
    else 
      fIndex=fgkInvalidIndex;

}

//______________________________________________________________________________
void AliMpNeighboursPadIterator::Next()
{
// pre-increment operator. Should be used by default for iterating over
// pads


  if (!IsValid()) return;
  
  if (fIndex < fPads.size()-1) 
    fIndex++; 
  else 
    Invalidate();
}

//______________________________________________________________________________
Bool_t AliMpNeighboursPadIterator::IsDone() const
{
// 
  return !IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpNeighboursPadIterator::CurrentItem() const 
{
// dereferencement operator
  if (!IsValid())
    return AliMpPad::Invalid();
  else
    return fPads[fIndex];
}

//______________________________________________________________________________
void AliMpNeighboursPadIterator::Invalidate()
{
// Let the iterator points to the invalid position
    fIndex=fgkInvalidIndex;
}

