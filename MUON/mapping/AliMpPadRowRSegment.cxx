// $Id$
// Category: sector
//
// Class AliMpPadRowRSegment
// -------------------------
// Class describing a pad row segment composed of the 
// the identic pads;
// the pads are placed from the offset (defined in the base class)
// to the right.
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpPadRowRSegment.h"
#include "AliMpPadRow.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"

ClassImp(AliMpPadRowRSegment)

//______________________________________________________________________________
AliMpPadRowRSegment::AliMpPadRowRSegment(AliMpPadRow* padRow, AliMpMotif* motif, 
                                         Int_t motifPositionId, Int_t nofPads)
  : AliMpVPadRowSegment(padRow, motif, motifPositionId, nofPads)
{
// 
}

//______________________________________________________________________________
AliMpPadRowRSegment::AliMpPadRowRSegment() 
  : AliMpVPadRowSegment()
{
//
}

//______________________________________________________________________________
AliMpPadRowRSegment::~AliMpPadRowRSegment() {
//  
}

//
// private methods  
//

//______________________________________________________________________________
Double_t AliMpPadRowRSegment::FirstPadCenterX() const
{
// Returns the x coordinate of the first (the most left) pad center
// in global coordinate system.
// ---

  return GetOffsetX() + GetMotif()->GetPadDimensions().X();
}  

//______________________________________________________________________________
Double_t AliMpPadRowRSegment::LastPadCenterX() const
{
// Returns the x coordinate of the last (the most right) pad center
// in global coordinate system.
// !! numbering of pads is in (-x) direction
// ---

  return GetOffsetX() + (2.*GetNofPads() - 1)*GetMotif()->GetPadDimensions().X();
}

//______________________________________________________________________________
Double_t AliMpPadRowRSegment::FirstPadBorderX() const
{
// Returns the x coordinate of the left border of the first (the most left) 
// pad in global coordinate system.
// ---

  return GetOffsetX();
         // Also could be
         // return FirstPadCenterX() + GetMotif()->GetPadDimensions().X();
}  

//______________________________________________________________________________
Double_t AliMpPadRowRSegment::LastPadBorderX() const
{
// Returns the x coordinate of the right border of the last (the most right)
// pad in global coordinate system.
// ---

  return LastPadCenterX() + GetMotif()->GetPadDimensions().X();
}  

//
// public methods  
//

//______________________________________________________________________________
Double_t  AliMpPadRowRSegment::LeftBorderX() const
{
// Returns the x coordinate of the left row segment border
// in global coordinate system.
// ---

  return FirstPadBorderX();
}

//______________________________________________________________________________
Double_t  AliMpPadRowRSegment::RightBorderX() const
{
// Returns the x coordinate of the right row segment border
// in global coordinate system.
// ---

  return LastPadBorderX();
}

