// $Id$
// Category: sector
//
// Class AliMpSectorAreaVPadIterator
// ---------------------------------
// Class, which defines an iterator over the pads 
// inside a given area in a sector in vertical direction.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>

#include "AliMpSectorAreaVPadIterator.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpConstants.h"

ClassImp(AliMpSectorAreaVPadIterator)

//______________________________________________________________________________
AliMpSectorAreaVPadIterator::AliMpSectorAreaVPadIterator(
                                const AliMpSectorSegmentation* segmentation,
                                const AliMpArea& area) 
 : AliMpVPadIterator(),
   fkSegmentation(segmentation),
   fkArea(area),
   fCurrentPad(AliMpPad::Invalid()),
   fCurrentColumnPosition(0.)
{
// normal constructor, start in invalid position
}

//______________________________________________________________________________
AliMpSectorAreaVPadIterator::AliMpSectorAreaVPadIterator(
                                const AliMpSectorAreaVPadIterator& right)
  : AliMpVPadIterator(right)
{
// copy constructor
 
  Fatal("Copy constructor", "Not implemented");
}

//______________________________________________________________________________
AliMpSectorAreaVPadIterator::AliMpSectorAreaVPadIterator()
 : AliMpVPadIterator(),
   fkSegmentation(0),
   fkArea(AliMpArea()),
   fCurrentPad(AliMpPad::Invalid()),
   fCurrentColumnPosition(0.)
{
// Dummy default constructor.
}

//______________________________________________________________________________
AliMpSectorAreaVPadIterator::~AliMpSectorAreaVPadIterator()
{
// destructor
}

//
// operators
//

//______________________________________________________________________________
AliMpSectorAreaVPadIterator& 
AliMpSectorAreaVPadIterator::operator = (const AliMpSectorAreaVPadIterator& right)
{
// Assignement operator

  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  AliMpVPadIterator::operator=(right);

  fkSegmentation = right.fkSegmentation;
  fkArea         = right.fkArea;
  fCurrentPad    = right.fCurrentPad;
  fCurrentColumnPosition = right.fCurrentColumnPosition;

  return *this;
} 

// 
// private methods
//

//______________________________________________________________________________
Bool_t AliMpSectorAreaVPadIterator::IsValid() const
{
// Is the iterator in a valid position?
// ---

  return fCurrentPad.IsValid() ;
}

//______________________________________________________________________________
void AliMpSectorAreaVPadIterator::MoveRight()
{
// Increase the current row position and searches the first valid pad.
// ---

  Double_t step = 2.* fkSegmentation->GetMinPadDimensions().X();

  while ( !fCurrentPad.IsValid() && 
          fCurrentColumnPosition + step < fkArea.RightBorder())
  {
    //cout << "#########  Move right ##########" << endl;
  
    fCurrentColumnPosition += step;
    TVector2 position = TVector2(fCurrentColumnPosition, fkArea.DownBorder());
    
    fCurrentPad = fkSegmentation->PadByDirection(position, fkArea.UpBorder());
  } 
}

//
// public methods
//

//______________________________________________________________________________
void AliMpSectorAreaVPadIterator::First()
{
// Reset the iterator, so that it points to the first available
// pad in the area
// ---

  if (!fkSegmentation) {
    Fatal("First", "Segmentation is not defined");
    return;
  }  

  // Start position = left down corner of the area
  //
  fCurrentColumnPosition = fkArea.DownBorder();
  TVector2 position(fkArea.LeftDownCorner()); 
  
  fCurrentPad = fkSegmentation->PadByDirection(position, fkArea.UpBorder());
  
  MoveRight();

  // Set the column position to the center of pad
  //
  if (fCurrentPad.IsValid()) fCurrentColumnPosition = fCurrentPad.Position().X();
}

//______________________________________________________________________________
void AliMpSectorAreaVPadIterator::Next()
{
// Move the iterator to the next valid pad.
// ---

  if (!IsValid()) return;
  
  // Start position = up board of current pad + little step
  //
  TVector2 position 
    = fCurrentPad.Position() 
    + TVector2(0., fCurrentPad.Dimensions().Y() + AliMpConstants::LengthStep());

  fCurrentPad = fkSegmentation->PadByDirection(position, fkArea.UpBorder());  
  
  if (fCurrentPad.IsValid()) return;

  MoveRight();
}

//______________________________________________________________________________
Bool_t AliMpSectorAreaVPadIterator::IsDone() const
{
// 
  return !IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpSectorAreaVPadIterator::CurrentItem () const 
{
// Returns current pad.
// ---

  return fCurrentPad;
}
//______________________________________________________________________________
void AliMpSectorAreaVPadIterator::Invalidate()
{
// 
  fCurrentPad = AliMpPad::Invalid();
  fCurrentColumnPosition = 0;
}



