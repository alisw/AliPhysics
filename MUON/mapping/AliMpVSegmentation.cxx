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

#include "TVector2.h"
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
  static TVector2* testPositions(0x0);
  static const Int_t kNofTestPositions(11);
  static const Double_t kEpsilon(AliMpConstants::LengthTolerance()*2.0);
  static Int_t centerIndex(-1);
  
  // testPositions are the positions (L,T,R,B) relative to pad's center (O)
  // were we'll try to get a neighbouring pad, by getting a little
  // bit outside the pad itself.
  // Note that it's not symmetric as we assume that pad density
  // can always decrease when going from left to right (or from bottom to top)
  //
  // L--T--R
  // |     |
  // L     |
  // |  O  R
  // L     |
  // |     |
  // L-B-B-R
  //
  // The order in which we actually test the positions has some importance,
  // i.e. when using this information to compute status map later on. Here's
  // the sequence :
  //
  // 4-- 5-- 6
  // |       |
  // 3       |
  // |   0   7
  // 2       |
  // |       |
  // 1-10- 9-8
  
  neighbours.Delete();
  neighbours.SetOwner(kTRUE);
  
  if (!pad.IsValid()) return 0;
  
  if (!testPositions)
  {
    testPositions = new TVector2[kNofTestPositions];
    Int_t n(0); 
    // starting center
    centerIndex = 0;
    testPositions[n++] = TVector2(0,0); // O (pad center)
    // then left column (L), starting from bottom
    testPositions[n++] = TVector2(-1,-1); // 1
    testPositions[n++] = TVector2(-1,-1/3.0); // 2 
    testPositions[n++] = TVector2(-1,1/3.0); // 3
    testPositions[n++] = TVector2(-1,1); // 4
    // top (T)
    testPositions[n++] = TVector2(0,1); // 5
    // right column (R), starting from top
    testPositions[n++] = TVector2(1,1); // 6
    testPositions[n++] = TVector2(1,0); // 7
    testPositions[n++] = TVector2(1,-1); // 8
    // bottom (B), starting from right
    testPositions[n++] = TVector2(1/3.0,-1); // 9 
    testPositions[n++] = TVector2(-1/3.0,-1); // 10
    // pad center
    if ( n != kNofTestPositions ) {
      AliError("Test on number of test positions failed.");
    }  
  }
  
  Int_t n(0);
  
  AliMpPad previous(AliMpPad::Invalid());
  
  for ( Int_t i = 0; i < kNofTestPositions; ++i ) 
  {
    if ( i == centerIndex && !includeSelf )
    {
      if ( includeVoid ) 
      {
        previous = AliMpPad::Invalid();
        neighbours.Add(new AliMpPad(previous));
        ++n;
      }
      continue;
    }
    
    TVector2 shift = testPositions[i];
    TVector2 pos = pad.Position();
    pos += TVector2((pad.Dimensions().X()+kEpsilon)*shift.X(),
                    (pad.Dimensions().Y()+kEpsilon)*shift.Y());

    
    AliMpPad p = PadByPosition(pos,kFALSE);
    
    if ( !p.IsValid() && !includeVoid ) continue;
    
    if ( p != previous || !previous.IsValid() ) 
    {
      previous = p;
      neighbours.Add(new AliMpPad(p));
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


