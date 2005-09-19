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
// $MpId: AliMpVIndexed.cxx,v 1.6 2005/08/26 15:43:36 ivana Exp $
// Category: basic
//
// Class AliMpVIndexed
// -------------------
// Class that defines the limits of global pad indices.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpVIndexed.h"

ClassImp(AliMpVIndexed)

//_____________________________________________________________________________
AliMpVIndexed::AliMpVIndexed(const AliMpIntPair& lowLimit, 
                             const AliMpIntPair& highLimit)
  : TObject(),
    fLowIndicesLimit(lowLimit),
    fHighIndicesLimit(highLimit) 
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMpVIndexed::AliMpVIndexed()
  : TObject(),
    fLowIndicesLimit(AliMpIntPair::Invalid()),
    fHighIndicesLimit(AliMpIntPair::Invalid()) 
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpVIndexed::~AliMpVIndexed()
{
/// Destructor 
}


//_____________________________________________________________________________
AliMpIntPair AliMpVIndexed::GlobalIndices(const AliMpIntPair& localIndices) const
{
/// Return the global indices corresponding to the given local indices.

  return GetLowIndicesLimit()+localIndices;

}

//_____________________________________________________________________________
Bool_t AliMpVIndexed::HasIndices(const AliMpIntPair& indices) const
{
/// Return true in the specified indices are within the limits.
  
  return (indices.GetFirst()  >= fLowIndicesLimit.GetFirst() && 
          indices.GetSecond() >= fLowIndicesLimit.GetSecond() && 
          indices.GetFirst()  <= fHighIndicesLimit.GetFirst() && 
          indices.GetSecond() <= fHighIndicesLimit.GetSecond() );
}  

//_____________________________________________________________________________
Bool_t AliMpVIndexed::HasValidIndices() const
{
/// Returns true if both indices limits have valid values.
  
  return (fLowIndicesLimit.IsValid() && fHighIndicesLimit.IsValid() );
}	   


  
