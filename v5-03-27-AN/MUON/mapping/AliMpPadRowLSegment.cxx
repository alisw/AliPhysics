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
// $MpId: AliMpPadRowLSegment.cxx,v 1.6 2006/05/24 13:58:46 ivana Exp $
// Category: sector

//-----------------------------------------------------------------------------
// Class AliMpPadRowLSegment
// -------------------------
// Class describing a pad row segment composed of the 
// the identic pads;
// the pads are placed from the offset (defined in the base class)
// to the left.
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpPadRowLSegment.h"
#include "AliMpPadRow.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"

/// \cond CLASSIMP
ClassImp(AliMpPadRowLSegment)
/// \endcond

//_____________________________________________________________________________
AliMpPadRowLSegment::AliMpPadRowLSegment(
                          AliMpPadRow* padRow, AliMpMotif* motif, 
                          Int_t motifPositionId, Int_t nofPads)
  : AliMpVPadRowSegment(padRow, motif, motifPositionId, nofPads)
{
/// Standard constructor 
}

//_____________________________________________________________________________
AliMpPadRowLSegment::AliMpPadRowLSegment() 
  : AliMpVPadRowSegment()
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpPadRowLSegment::~AliMpPadRowLSegment() 
{
/// Destructor  
}

//
// private methods  
//

//_____________________________________________________________________________
Double_t AliMpPadRowLSegment::FirstPadCenterX() const
{
/// Return the x coordinate of the first (the most right) pad center
/// in the global coordinate system.

  return GetOffsetX() - GetMotif()->GetPadDimensionX();
}  

//_____________________________________________________________________________
Double_t AliMpPadRowLSegment::LastPadCenterX() const
{
/// Return the x coordinate of the last (the most left) pad center
/// in the global coordinate system.                                         \n       
/// !! numbering of pads is in (-x) direction

  return GetOffsetX() - (2.*GetNofPads() - 1)*GetMotif()->GetPadDimensionX();
}

//_____________________________________________________________________________
Double_t AliMpPadRowLSegment::FirstPadBorderX() const
{
/// Return the x coordinate of the right border of the first (the most right) 
/// pad in the global coordinate system.

  return GetOffsetX();
         // Also could be
         // return FirstPadCenterX() + GetMotif()->GetPadDimensionX();
}  

//_____________________________________________________________________________
Double_t AliMpPadRowLSegment::LastPadBorderX() const
{
/// Return the x coordinate of the left border of the last (the most left)
/// pad in the global coordinate system.

  return LastPadCenterX() - GetMotif()->GetPadDimensionX();
}  

//
// public methods  
//

//_____________________________________________________________________________
Double_t  AliMpPadRowLSegment::LeftBorderX() const
{
/// Return the x coordinate of the left row segment border
/// in the global coordinate system.

  return LastPadBorderX();
}

//_____________________________________________________________________________
Double_t  AliMpPadRowLSegment::RightBorderX() const
{
/// Return the x coordinate of the right row segment border
/// in the global coordinate system.

  return FirstPadBorderX();
}
