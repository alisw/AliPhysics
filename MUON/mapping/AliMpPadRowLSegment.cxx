// $Id$
// Category: sector
//
// Class AliMpPadRowLSegment
// -------------------------
// Class describing a pad row segment composed of the 
// the identic pads;
// the pads are placed from the offset (defined in the base class)
// to the left.
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpPadRowLSegment.h"
#include "AliMpPadRow.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"

ClassImp(AliMpPadRowLSegment)

//_____________________________________________________________________________
AliMpPadRowLSegment::AliMpPadRowLSegment(
                          AliMpPadRow* padRow, AliMpMotif* motif, 
                          Int_t motifPositionId, Int_t nofPads)
  : AliMpVPadRowSegment(padRow, motif, motifPositionId, nofPads)
{
// 
}

//_____________________________________________________________________________
AliMpPadRowLSegment::AliMpPadRowLSegment() 
  : AliMpVPadRowSegment()
{
//
}

//_____________________________________________________________________________
AliMpPadRowLSegment::~AliMpPadRowLSegment() {
//  
}

//
// private methods  
//

//_____________________________________________________________________________
Double_t AliMpPadRowLSegment::FirstPadCenterX() const
{
// Returns the x coordinate of the first (the most right) pad center
// in global coordinate system.
// ---

  return GetOffsetX() - GetMotif()->GetPadDimensions().X();
}  

//_____________________________________________________________________________
Double_t AliMpPadRowLSegment::LastPadCenterX() const
{
// Returns the x coordinate of the last (the most left) pad center
// in global coordinate system.
// !! numbering of pads is in (-x) direction
// ---

  return GetOffsetX() - (2.*GetNofPads() - 1)*GetMotif()->GetPadDimensions().X();
}

//_____________________________________________________________________________
Double_t AliMpPadRowLSegment::FirstPadBorderX() const
{
// Returns the x coordinate of the right border of the first (the most right) 
// pad in global coordinate system.
// ---

  return GetOffsetX();
         // Also could be
         // return FirstPadCenterX() + GetMotif()->GetPadDimensions().X();
}  

//_____________________________________________________________________________
Double_t AliMpPadRowLSegment::LastPadBorderX() const
{
// Returns the x coordinate of the left border of the last (the most left)
// pad in global coordinate system.
// ---

  return LastPadCenterX() - GetMotif()->GetPadDimensions().X();
}  

//
// public methods  
//

//_____________________________________________________________________________
Double_t  AliMpPadRowLSegment::LeftBorderX() const
{
// Returns the x coordinate of the left row segment border
// in global coordinate system.
// ---

  return LastPadBorderX();
}

//_____________________________________________________________________________
Double_t  AliMpPadRowLSegment::RightBorderX() const
{
// Returns the x coordinate of the right row segment border
// in global coordinate system.
// ---

  return FirstPadBorderX();
}
