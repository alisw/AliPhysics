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
// $MpId: AliMpPadRowRSegment.cxx,v 1.6 2006/05/24 13:58:46 ivana Exp $
// Category: sector

//-----------------------------------------------------------------------------
// Class AliMpPadRowRSegment
// -------------------------
// Class describing a pad row segment composed of the 
// the identic pads;
// the pads are placed from the offset (defined in the base class)
// to the right.
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpPadRowRSegment.h"
#include "AliMpPadRow.h"
#include "AliMpMotif.h"
#include "AliMpMotifType.h"

/// \cond CLASSIMP
ClassImp(AliMpPadRowRSegment)
/// \endcond

//______________________________________________________________________________
AliMpPadRowRSegment::AliMpPadRowRSegment(AliMpPadRow* padRow, AliMpMotif* motif, 
                                         Int_t motifPositionId, Int_t nofPads)
  : AliMpVPadRowSegment(padRow, motif, motifPositionId, nofPads)
{
/// Standard constructor 
}

//______________________________________________________________________________
AliMpPadRowRSegment::AliMpPadRowRSegment() 
  : AliMpVPadRowSegment()
{
/// Default constructor
}

//______________________________________________________________________________
AliMpPadRowRSegment::~AliMpPadRowRSegment() 
{
/// Destructor  
}

//
// private methods  
//

//______________________________________________________________________________
Double_t AliMpPadRowRSegment::FirstPadCenterX() const
{
/// Return the x coordinate of the first (the most left) pad center
/// in the global coordinate system.

  return GetOffsetX() + GetMotif()->GetPadDimensionX();
}  

//______________________________________________________________________________
Double_t AliMpPadRowRSegment::LastPadCenterX() const
{
/// Return the x coordinate of the last (the most right) pad center
/// in the global coordinate system.                                             \n
/// !! numbering of pads is in (-x) direction

  return GetOffsetX() + (2.*GetNofPads() - 1)*GetMotif()->GetPadDimensionX();
}

//______________________________________________________________________________
Double_t AliMpPadRowRSegment::FirstPadBorderX() const
{
/// Return the x coordinate of the left border of the first (the most left) 
/// pad in the global coordinate system.

  return GetOffsetX();
         // Also could be
         // return FirstPadCenterX() + GetMotif()->GetPadDimensionX();
}  

//______________________________________________________________________________
Double_t AliMpPadRowRSegment::LastPadBorderX() const
{
/// Return the x coordinate of the right border of the last (the most right)
/// pad in the global coordinate system.

  return LastPadCenterX() + GetMotif()->GetPadDimensionX();
}  

//
// public methods  
//

//______________________________________________________________________________
Double_t  AliMpPadRowRSegment::LeftBorderX() const
{
/// Return the x coordinate of the left row segment border
/// in the global coordinate system.

  return FirstPadBorderX();
}

//______________________________________________________________________________
Double_t  AliMpPadRowRSegment::RightBorderX() const
{
/// Return the x coordinate of the right row segment border
/// in the global coordinate system.

  return LastPadBorderX();
}

