// $Id$
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

ClassImp(AliMpVSegmentation)

//_____________________________________________________________________________
AliMpVSegmentation::AliMpVSegmentation() 
  : TObject()
{
//
}

//_____________________________________________________________________________
AliMpVSegmentation::~AliMpVSegmentation() {
// 
}

//
// private methods
//

//_____________________________________________________________________________
AliMpPadPair AliMpVSegmentation::FindPads(const TVector2& position1, 
                                          const TVector2& position2) const
{
// Returns a pair of pads with specified position.
// If both pads are identical, the second pad in pair is set to invalid.
// ---

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
// Returns a pair of pads neighbouring up to the specified pad.
// If there is only one neighbouring pad,
// the second pad in pair is invalid.
// ---

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
// Returns a pair of pads neighbouring down to the specified pad.
// If there is only one neighbouring pad,
// the second pad in pair is invalid.
// ---

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
// Returns a pair of pads neighbouring left to the specified pad.
// If there is only one neighbouring pad,
// the second in pair is invalid.
// ---

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
// Returns a pair of pads neighbouring right to the specified pad.
// If there is only one neighbouring pad,
// the second in pair is invalid.
// ---

  TVector2 position1 
    = pad.Position() + TVector2(pad.Dimensions().X() + AliMpConstants::LengthStep(),
			        (-1.)*AliMpConstants::LengthStep()); 
  TVector2 position2
    = pad.Position() + TVector2(pad.Dimensions().X() + AliMpConstants::LengthStep(),
		                AliMpConstants::LengthStep()); 
				     
  return FindPads(position1, position2);
}

