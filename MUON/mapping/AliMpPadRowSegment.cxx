// $Id$
// Category: sector
//
// Class AliMpPadRowSegment
// --------------------
// Class describing a pad row segment composed of the 
// the identic pads.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpPadRowSegment.h"
#include "AliMpPadRow.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"

ClassImp(AliMpPadRowSegment)

//_____________________________________________________________________________
AliMpPadRowSegment::AliMpPadRowSegment(AliMpPadRow* padRow, AliMpMotif* motif, 
                                       Int_t motifPositionId, Int_t nofPads)
  : TObject(),
    fNofPads(nofPads),
    fOffsetX(0.),
    fPadRow(padRow),
    fMotif(motif),
    fMotifPositionId(motifPositionId)
{
// 
}

//_____________________________________________________________________________
AliMpPadRowSegment::AliMpPadRowSegment() 
  : TObject(),
    fNofPads(0),
    fOffsetX(0.),
    fPadRow(0),
    fMotif(0),
    fMotifPositionId(0)
{
//
}

//_____________________________________________________________________________
AliMpPadRowSegment::~AliMpPadRowSegment() {
//  
}

//
// private methods  
//

//_____________________________________________________________________________
Double_t AliMpPadRowSegment::FirstPadCenterX() const
{
// Returns the x coordinate of the first (the most right) pad center
// in global coordinate system.
// ---

  return fOffsetX - fMotif->GetPadDimensions().X();
}  

//_____________________________________________________________________________
Double_t AliMpPadRowSegment::LastPadCenterX() const
{
// Returns the x coordinate of the last (the most left) pad center
// in global coordinate system.
// !! numbering of pads is in (-x) direction
// ---

  return fOffsetX - (2.*fNofPads - 1)*fMotif->GetPadDimensions().X();
}

//_____________________________________________________________________________
Double_t AliMpPadRowSegment::FirstPadBorderX() const
{
// Returns the x coordinate of the right border of the first (the most right) 
// pad in global coordinate system.
// ---

  return fOffsetX;
         // Also could be
         // return FirstPadCenterX() + fMotif->GetPadDimensions().X();
}  

//_____________________________________________________________________________
Double_t AliMpPadRowSegment::LastPadBorderX() const
{
// Returns the x coordinate of the left border of the last (the most left)
// pad in global coordinate system.
// ---

  return LastPadCenterX() - fMotif->GetPadDimensions().X();
}  

//
// public methods  
//

//_____________________________________________________________________________
Double_t  AliMpPadRowSegment::LeftBorderX() const
{
// Returns the x coordinate of the left row segment border
// in global coordinate system.
// ---

  return LastPadBorderX();
}

//_____________________________________________________________________________
Double_t  AliMpPadRowSegment::RightBorderX() const
{
// Returns the x coordinate of the right row segment border
// in global coordinate system.
// ---

  return FirstPadBorderX();
}

//_____________________________________________________________________________
Double_t  AliMpPadRowSegment::HalfSizeY() const
{
// Returns the size in y of this row segment.
// ---

  return fMotif->GetPadDimensions().Y();
}

//_____________________________________________________________________________
AliMpPadRow*  AliMpPadRowSegment::GetPadRow() const
{
// Returns the pad row.which this pad row segment belongs to.
// ---

  return fPadRow;
}  

//_____________________________________________________________________________
AliMpMotif*  AliMpPadRowSegment::GetMotif() const
{
// Returns the motif of this pad row segment. 
// ---

  return fMotif;
}  

//_____________________________________________________________________________
Int_t  AliMpPadRowSegment::GetMotifPositionId() const
{
// Returns the motif of this pad row segment. 
// ---

  return fMotifPositionId;
}  

//_____________________________________________________________________________
void  AliMpPadRowSegment::SetOffsetX(Double_t offsetX)
{
// Sets the x offset.
// ---

  fOffsetX = offsetX;
}    

