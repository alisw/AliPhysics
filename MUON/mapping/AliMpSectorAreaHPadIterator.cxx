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
// $MpId: AliMpSectorAreaHPadIterator.cxx,v 1.7 2006/05/24 13:58:46 ivana Exp $
// Category: sector
//
// Class AliMpSectorAreaHPadIterator
// ---------------------------------
// Class, which defines an iterator over the pads 
// inside a given area in a sector in horizontal direction.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpSectorAreaHPadIterator.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpConstants.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMpSectorAreaHPadIterator)
/// \endcond

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
/// Standard constructor, start in invalid position
}

//______________________________________________________________________________
AliMpSectorAreaHPadIterator::AliMpSectorAreaHPadIterator(
                                const AliMpSectorAreaHPadIterator& right)
  : AliMpVPadIterator(right),
    fkSegmentation(0),
    fkArea(AliMpArea()),
    fCurrentPad(AliMpPad::Invalid()),
    fCurrentRowPosition(0.)
{
/// Copy constructor
 
  *this = right;
}

//______________________________________________________________________________
AliMpSectorAreaHPadIterator::AliMpSectorAreaHPadIterator()
 : AliMpVPadIterator(),
   fkSegmentation(0),
   fkArea(AliMpArea()),
   fCurrentPad(AliMpPad::Invalid()),
   fCurrentRowPosition(0.)
{
/// Default constructor.
}

//______________________________________________________________________________
AliMpSectorAreaHPadIterator::~AliMpSectorAreaHPadIterator()
{
/// Destructor
}

//
// operators
//

//______________________________________________________________________________
AliMpSectorAreaHPadIterator& 
AliMpSectorAreaHPadIterator::operator = (const AliMpSectorAreaHPadIterator& right)
{
/// Assignment operator

  // check assignment to self
  if (this == &right) return *this;

  // base class assignment
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
/// Is the iterator in a valid position?

  return fCurrentPad.IsValid() ;
}

//______________________________________________________________________________
void AliMpSectorAreaHPadIterator::MoveUp()
{
/// Increase the current row position and searches the first valid pad.

  Double_t dy = fkSegmentation->GetMinPadDimensions().Y();

  while ( !fCurrentPad.IsValid() && 
          fCurrentRowPosition + dy < fkArea.UpBorder())
  {
    fCurrentRowPosition += 2.*dy;
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
/// Reset the iterator, so that it points to the first available
/// pad in the area

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
/// Move the iterator to the next valid pad.

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
/// Is the iterator in the end ?
 
  return !IsValid();
}

//______________________________________________________________________________
AliMpPad AliMpSectorAreaHPadIterator::CurrentItem () const 
{
/// Return current pad.

  return fCurrentPad;
}
//______________________________________________________________________________
void AliMpSectorAreaHPadIterator::Invalidate()
{
/// Let the iterator point to the invalid position

  fCurrentPad = AliMpPad::Invalid();
  fCurrentRowPosition = 0;
}



