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

//-----------------------------------------------------------------------------
// Class AliMpVSegmentation
// ------------------------
// The abstract base class for the segmentation.
// Provides methods related to pads:
// conversion between pad indices, pad location, pad position;
// finding pad neighbour.
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//         Laurent Aphecetche, SUBATECH
//-----------------------------------------------------------------------------


#include "AliMpVSegmentation.h"
#include "AliMpArea.h"
#include "AliMpConstants.h"

#include "AliLog.h"

#include "TObjArray.h"

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

//_____________________________________________________________________________
Int_t 
AliMpVSegmentation::GetNeighbours(const AliMpPad& pad, 
                                  TObjArray& neighbours,
                                  Bool_t includeSelf,
                                  Bool_t includeVoid) const
{
  /// Returns the list of neighbours of pad
  /// testPositions are the positions (L,T,R,B) relative to pad's center (O)
  /// were we'll try to get a neighbouring pad, by getting a little
  /// bit outside the pad itself.
  /// The pad density can only decrease when going from left to right except
  /// for round slates where it is the opposite.
  /// The pad density can only decrease when going from bottom to top but
  /// to be symmetric we also consider the opposite.
  /// The order in which we actually test the positions has some importance,
  /// i.e. when using this information to compute status map later on. Here's
  /// the sequence :
  /// <pre>
  /// 4- 5- 6-7
  /// |       |
  /// 3       8
  /// |   0   |
  /// 2       9
  /// |       |
  /// 1-12-11-10
  /// </pre>

  static const Int_t kNofTestPositions(12);
  static Double_t shiftx[12] = {-1., -1., -1., -1., -1./3., 1./3., 1., 1., 1., 1., 1./3., -1./3.};
  static Double_t shifty[12] = {-1., -1./3., 1./3., 1., 1., 1., 1., 1./3., -1./3., -1., -1., -1.};

  static const Double_t kEpsilon(AliMpConstants::LengthTolerance()*2.0);
  
  neighbours.Delete();
  neighbours.SetOwner(kTRUE);
  
  if ( ! pad.IsValid() ) return 0;
  
  AliMpPad invalid(AliMpPad::Invalid());
  AliMpPad *previous = &invalid;
  Int_t n(0);
  
  // consider adding the pad itself
  if ( includeSelf )
  {
    neighbours.Add(new AliMpPad(pad));
    ++n;
  }
  else if ( includeVoid )
  {
    neighbours.Add(new AliMpPad(invalid));
    ++n;
  }
  
  // add the neighbours (only once)
  for ( Int_t i = 0; i < kNofTestPositions; ++i )
  {
    
    AliMpPad p
      = PadByPosition(pad.GetPositionX() + ( pad.GetDimensionX() + kEpsilon )*shiftx[i], 
                      pad.GetPositionY() + ( pad.GetDimensionY() + kEpsilon )*shifty[i],
                      kFALSE);
    
    if ( p.IsValid() && p != *previous )
    {
      previous = new AliMpPad(p);
      neighbours.Add(previous);
      ++n;
    }
    else if ( includeVoid )
    {
      neighbours.Add(new AliMpPad(invalid));
      ++n;
    }
    
  }
  
  return n;
  
}

//
// public methods
//

//_____________________________________________________________________________
Bool_t 
AliMpVSegmentation::HasPadByIndices(Int_t ix, Int_t iy) const
{
  /// Default implementation. Must be overwritten if can be made more
  /// efficient in the child class
  
  return ( PadByIndices(ix, iy, kFALSE) != AliMpPad::Invalid() );
}

//_____________________________________________________________________________
Bool_t 
AliMpVSegmentation::HasPadByLocation(Int_t manuId, Int_t manuChannel) const
{
  /// Default implementation. Must be overwritten if can be made more
  /// efficient in the child class
  
  return (PadByLocation(manuId, manuChannel, kFALSE) != AliMpPad::Invalid());
}

//_____________________________________________________________________________
Bool_t 
AliMpVSegmentation::HasMotifPosition(Int_t manuId) const
{
  /// Default implementation to know if we hold a given manu
  return ( MotifPosition(manuId) != 0x0 );
}


