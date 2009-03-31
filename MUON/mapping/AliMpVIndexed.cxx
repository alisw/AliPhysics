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
// $MpId: AliMpVIndexed.cxx,v 1.7 2006/05/24 13:58:29 ivana Exp $
// Category: basic

//-----------------------------------------------------------------------------
// Class AliMpVIndexed
// -------------------
// Class that defines the limits of global pad indices.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpVIndexed.h"

/// \cond CLASSIMP
ClassImp(AliMpVIndexed)
/// \endcond

//_____________________________________________________________________________
AliMpVIndexed::AliMpVIndexed()
  : TObject(),
    fLowLimit(0),
    fHighLimit(0),
    fLowValid(false),
    fHighValid(false)
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpVIndexed::~AliMpVIndexed()
{
/// Destructor 
}

//_____________________________________________________________________________
MpPair_t AliMpVIndexed::GlobalIndices(MpPair_t localIndices) const
{
/// Return the global indices corresponding to the given local indices.

  return fLowLimit + localIndices;
}

//_____________________________________________________________________________
Int_t  AliMpVIndexed::GlobalIx(Int_t localIx) const
{
/// Return the global indices ix corresponding to the given local indices

  return GetLowLimitIx() + localIx;
}  


//_____________________________________________________________________________
Int_t  AliMpVIndexed::GlobalIy(Int_t localIy) const
{
/// Return the global indices iy corresponding to the given local indices

  return GetLowLimitIy() + localIy;
}  

//_____________________________________________________________________________
void AliMpVIndexed::SetLowIndicesLimit(MpPair_t limit, Bool_t valid)
{ 
/// Set low indices limit

  fLowLimit = limit; 
  fLowValid = valid ; 
}
  
//_____________________________________________________________________________
void AliMpVIndexed::SetLowIndicesLimit(Int_t ix, Int_t iy, Bool_t valid)
{ 
/// Set low indices limit

  fLowLimit = AliMp::Pair(ix, iy); 
  fLowValid = valid; 
}
  
//_____________________________________________________________________________
void AliMpVIndexed::SetHighIndicesLimit(MpPair_t limit, Bool_t valid)
{ 
/// Set high indices limit

  fHighLimit = limit; 
  fHighValid = valid ; 
}
  
//_____________________________________________________________________________
void AliMpVIndexed::SetHighIndicesLimit(Int_t ix, Int_t iy, Bool_t valid)
{ 
/// Set high indices limit

  fHighLimit = AliMp::Pair(ix, iy); 
  fHighValid = valid; 
}

//_____________________________________________________________________________
Bool_t AliMpVIndexed::HasIndices(MpPair_t indices) const
{
/// Return true in the specified indices are within the limits.
  
  return ( AliMp::PairFirst(indices)  >= GetLowLimitIx() && 
           AliMp::PairSecond(indices) >= GetLowLimitIy() && 
           AliMp::PairFirst(indices)  <= GetHighLimitIx() && 
           AliMp::PairSecond(indices) <= GetHighLimitIy() );
}  

//_____________________________________________________________________________
Bool_t AliMpVIndexed::HasIndices(Int_t ix, Int_t iy) const
{
/// Return true in the specified indices are within the limits.
  
  return (ix  >= GetLowLimitIx() && 
          iy  >= GetLowLimitIy() && 
          ix  <= GetHighLimitIx() && 
          iy  <= GetHighLimitIy() );
}  

//_____________________________________________________________________________
Bool_t AliMpVIndexed::HasValidIndices() const
{
/// Returns true if both indices limits have valid values.
  
  return ( fLowValid && fHighValid );
}	   

//_____________________________________________________________________________
MpPair_t AliMpVIndexed::GetLowIndicesLimit() const
{ 
/// Return low indices limit

  // if ( ! fLowValid )  return 0;

  return fLowLimit; 
}

//_____________________________________________________________________________
Int_t  AliMpVIndexed::GetLowLimitIx() const
{ 
/// Return low indices ix limit

  // if ( ! fLowValid )  return 0;

  return AliMp::PairFirst(fLowLimit); 
}

//_____________________________________________________________________________
Int_t  AliMpVIndexed::GetLowLimitIy() const
{ 
/// Return low indices iy limit

  // if ( ! fLowValid )  return 0;

  return AliMp::PairSecond(fLowLimit); 
}

//_____________________________________________________________________________
Bool_t AliMpVIndexed::IsLowLimitValid() const  
{
/// Return true, if low indices limit is set 

  return fLowValid; 
}

//_____________________________________________________________________________
MpPair_t AliMpVIndexed::GetHighIndicesLimit() const
{ 
/// Return high indices limit

  // if ( ! fHighValid )  return 0;

  return fHighLimit; 
}

//_____________________________________________________________________________
Int_t  AliMpVIndexed::GetHighLimitIx() const
{ 
/// Return high indices ix limit

  // if ( ! fHighValid )  return 0;

  return AliMp::PairFirst(fHighLimit); 
}

//_____________________________________________________________________________
Int_t  AliMpVIndexed::GetHighLimitIy() const
{ 
/// Return high indices iy limit

  // if ( ! fHighValid )  return 0;

  return AliMp::PairSecond(fHighLimit); 
}

//_____________________________________________________________________________
Bool_t AliMpVIndexed::IsHighLimitValid() const  
{
/// Return true, if high indices limit is set 

  return fHighValid; 
}



  
