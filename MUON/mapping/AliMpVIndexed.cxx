// $Id$
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
    fHighIndicesLimit(highLimit) {
//
}

//_____________________________________________________________________________
AliMpVIndexed::AliMpVIndexed()
  : TObject(),
    fLowIndicesLimit(AliMpIntPair::Invalid()),
    fHighIndicesLimit(AliMpIntPair::Invalid()) {
//
}

//_____________________________________________________________________________
AliMpVIndexed::~AliMpVIndexed(){
// 
}


//_____________________________________________________________________________
AliMpIntPair AliMpVIndexed::GlobalIndices(const AliMpIntPair& localIndices) const
{
// Returns the global indices corresponding to the given local indices.
// ---

  return GetLowIndicesLimit()+localIndices;

}

//_____________________________________________________________________________
Bool_t AliMpVIndexed::HasIndices(const AliMpIntPair& indices) const
{
// Returns true in the specified indices are within the limits.
// ---
  
  return (indices.GetFirst()  >= fLowIndicesLimit.GetFirst() && 
          indices.GetSecond() >= fLowIndicesLimit.GetSecond() && 
          indices.GetFirst()  <= fHighIndicesLimit.GetFirst() && 
          indices.GetSecond() <= fHighIndicesLimit.GetSecond() );
}  

//_____________________________________________________________________________
Bool_t AliMpVIndexed::HasValidIndices() const
{
// Returns true if both indices limits have valid values.
// ---
  
  return (fLowIndicesLimit.IsValid() && fHighIndicesLimit.IsValid() );
}	   


  
