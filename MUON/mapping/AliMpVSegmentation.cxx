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
// $MpId: AliMpVSegmentation.cxx,v 1.5 2006/05/24 13:58:29 ivana Exp $
// Category: basic
//
// Class AliMpVSegmentation
// ------------------------
// The abstract base class for the segmentation.
// Provides methods related to pads:
// conversion between pad indices, pad location, pad position;
// finding pad neighbour.
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpVSegmentation.h"
#include "AliMpConstants.h"

/// \cond CLASSIMP
ClassImp(AliMpVSegmentation)
/// \endcond

//_____________________________________________________________________________
AliMpVSegmentation::AliMpVSegmentation() 
  : TObject()
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpVSegmentation::~AliMpVSegmentation() 
{
/// Destructor 
}

//
// private methods
//

//_____________________________________________________________________________
AliMpPadPair AliMpVSegmentation::FindPads(const TVector2& position1, 
                                          const TVector2& position2) const
{
/// Return a pair of pads with specified position.
/// If both pads are identical, the second pad in pair is set to invalid.

  AliMpPad pad1 = PadByPosition(position1, false);
  AliMpPad pad2 = PadByPosition(position2, false);
				     
  if (pad1 == pad2) pad2 = AliMpPad::Invalid();				     

  return AliMpPadPair(pad1, pad2);
}  

//
// public methods
//

//_____________________________________________________________________________
AliMpPadPair AliMpVSegmentation::PadsUp(const AliMpPad& pad) const
{
/// Return a pair of pads neighbouring up to the specified pad.
/// If there is only one neighbouring pad,
/// the second pad in pair is invalid.

  TVector2 position1 
    = pad.Position()+ TVector2((-1.)*AliMpConstants::LengthStep(), 
	                       pad.Dimensions().Y()+ AliMpConstants::LengthStep());
  TVector2 position2 
    = pad.Position()+ TVector2(AliMpConstants::LengthStep(), 
	                       pad.Dimensions().Y()+ AliMpConstants::LengthStep());
			       
  return FindPads(position1, position2);
}

//_____________________________________________________________________________
AliMpPadPair AliMpVSegmentation::PadsDown(const AliMpPad& pad) const
{
/// Return a pair of pads neighbouring down to the specified pad.
/// If there is only one neighbouring pad,
/// the second pad in pair is invalid.

  TVector2 position1 
    = pad.Position()- TVector2(AliMpConstants::LengthStep(), 
	                       pad.Dimensions().Y()+ AliMpConstants::LengthStep());

  TVector2 position2
    = pad.Position()- TVector2((-1.)*AliMpConstants::LengthStep(), 
	                       pad.Dimensions().Y()+ AliMpConstants::LengthStep());
				     
  return FindPads(position1, position2);
}

//_____________________________________________________________________________
AliMpPadPair AliMpVSegmentation::PadsLeft(const AliMpPad& pad) const
{
/// Return a pair of pads neighbouring left to the specified pad.
/// If there is only one neighbouring pad,
/// the second in pair is invalid.

  TVector2 position1 
    = pad.Position() - TVector2(pad.Dimensions().X() + AliMpConstants::LengthStep(),
			        AliMpConstants::LengthStep()); 
  TVector2 position2
    = pad.Position() - TVector2(pad.Dimensions().X() + AliMpConstants::LengthStep(),
			        (-1.)*AliMpConstants::LengthStep()); 

  return FindPads(position1, position2);
}

//_____________________________________________________________________________
AliMpPadPair AliMpVSegmentation::PadsRight(const AliMpPad& pad) const
{
/// Return a pair of pads neighbouring right to the specified pad.
/// If there is only one neighbouring pad,
/// the second in pair is invalid.

  TVector2 position1 
    = pad.Position() + TVector2(pad.Dimensions().X() + AliMpConstants::LengthStep(),
			        (-1.)*AliMpConstants::LengthStep()); 
  TVector2 position2
    = pad.Position() + TVector2(pad.Dimensions().X() + AliMpConstants::LengthStep(),
		                AliMpConstants::LengthStep()); 
				     
  return FindPads(position1, position2);
}

