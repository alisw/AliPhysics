// $Id$
// Category: sector
//
// Class AliMpSectorAreaHPadIterator
// ---------------------------------
// Class, which defines an iterator over the pads 
// inside a given area in a sector in horizontal direction.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <Riostream.h>

#include "AliMpSectorAreaHPadIterator.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpConstants.h"

ClassImp(AliMpSectorAreaHPadIterator)

//______________________________________________________________________________
AliMpSectorAreaHPadIterator::AliMpSectorAreaHPadIterator(
                                const AliMpSectorSegmentation* segmentation,
                                const AliMpArea& area) 
 : AliMpVPadIterator(),
   fkSegmentation(segmentation),
   fkArea(area),
   fCurrentPad(AliMpPad::Invalid()),
   fCurrentRowPosition(0.)
{
// normal constructor, start in invalid position
}

//______________________________________________________________________________
AliMpSectorAreaHPadIterator::AliMpSectorAreaHPadIterator(
                                const AliMpSectorAreaHPadIterator& right)
  : AliMpVPadIterator(right)
{
// copy constructor
 
  Fatal("Copy constructor", "Not implemented");
}

//______________________________________________________________________________
AliMpSectorAreaHPadIterator::AliMpSectorAreaHPadIterator()
 : AliMpVPadIterator(),
   fkSegmentation(0),
   fkArea(AliMpArea()),
   fCurrentPad(AliMpPad::Invalid()),
   fCurrentRowPosition(0.)
{
// Dummy default constructor.
}

//______________________________________________________________________________
AliMpSectorAreaHPadIterator::~AliMpSectorAreaHPadIterator()
{
// destructor
}

//
// operators
//

//______________________________________________________________________________
AliMpSectorAreaHPadIterator& 
AliMpSectorAreaHPadIterator::operator = (const AliMpSectorAreaHPadIterator& right)
{
// Assignement operator

  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  AliMpVPadIterator::operator=(right);

  fkSegmentation = right.fkSegmentation;
  fkArea         = right.fkArea;
  fCurrentPad    = right.fCurrentPad;
  fCurrentRowPosition = right.fCurrentRowPosition;

  return *this;
} 

// 
// private methods
//

//______________________________________________________________________________
Bool_t AliMpSectorAreaHPadIterator::IsValid() const
{
// Is the iterator in a valid position?
// ---

  return fCurrentPad.IsValid() ;
}

//______________________________________________________________________________
void AliMpSectorAreaHPadIterator::MoveUp()
{
// Increase the current row position and searches the first valid pad.
// ---

  Double_t step = 2.* fkSegmentation->GetMinPadDimensions().Y();

  while ( !fCurrentPad.IsValid() && 
          fCurrentRowPosition + step < fkArea.UpBorder())
  {
    //cout << "#########  Move up ##########" << endl;
  
    fCurrentRowPosition += step;
    TVector2 position = TVector2(fkArea.LeftBorder(), fCurrentRowPosition);
    
    fCurrentPad = fkSegmentation->PadByDirection(position, fkArea.RightBorder());
  } 
}

//
// public methods
//

//______________________________________________________________________________
void AliMpSectorAreaHPadIterator::First()
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
  fCurrentRowPosition = fkArea.DownBorder();
  TVector2 position(fkArea.LeftDownCorner()); 
  
  fCurrentPad = fkSegmentation->PadByDirection(position, fkArea.RightBorder());

  MoveUp();
  
  // Set the row position to the center of pad
  //
  if (fCurrentPad.IsValid()) fCurrentRowPosition = fCurrentPad.Position().Y();
}

//______________________________________________________________________________
void AliMpSectorAreaHPadIterator::Next()
{
// Move the iterator to the next valid pad.
// ---

  if (!IsValid()) return;
  
  // Start position = right board of current pad + little step
  //
  TVector2 position 
    = fCurrentPad.Position() 
    + TVector2(fCurrentPad.Dimensions().X() + AliMpConstants::LengthStep(), 0.);

  fCurrentPad = fkSegmentation->PadByDirection(position, fkArea.RightBorder());  

  if (fCurrentPad.IsValid()) return;

  MoveUp();
}

//______________________________________________________________________________
Bool_t AliMpSectorAreaHPadIterator::IsDone() const
{
// 
  return !IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpSectorAreaHPadIterator::CurrentItem () const 
{
// Returns current pad.
// ---

  return fCurrentPad;
}
//______________________________________________________________________________
void AliMpSectorAreaHPadIterator::Invalidate()
{
// 
  fCurrentPad = AliMpPad::Invalid();
  fCurrentRowPosition = 0;
}



