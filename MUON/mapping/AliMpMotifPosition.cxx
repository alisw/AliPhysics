// $Id$
//
// Class AliMpMotifPosition
// ------------------------
// Class that represents a placed motif.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TError.h>

#include "AliMpMotifPosition.h"
#include "AliMpMotifPositionPadIterator.h"
#include "AliMpMotifType.h"

ClassImp(AliMpMotifPosition)

//______________________________________________________________________________
AliMpMotifPosition::AliMpMotifPosition(Int_t id, AliMpVMotif* motif, 
                                       TVector2 position)
  : AliMpVIndexed(),
    fID(id),
    fMotif(motif),
    fPosition(position) {
//
}

//______________________________________________________________________________
AliMpMotifPosition::AliMpMotifPosition()
  : AliMpVIndexed(), 
    fID(0),
    fMotif(0),
    fPosition(TVector2(0.,0.)) {
//
}

//_____________________________________________________________________________
AliMpMotifPosition::AliMpMotifPosition(const AliMpMotifPosition& right) 
  : AliMpVIndexed(right) {
// 
  Fatal("AliMpMotifPosition", "Copy constructor not provided.");
}

//______________________________________________________________________________
AliMpMotifPosition::~AliMpMotifPosition(){
// 
}

// operators

//_____________________________________________________________________________
AliMpMotifPosition& 
AliMpMotifPosition::operator=(const AliMpMotifPosition& right)
{
  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//______________________________________________________________________________
AliMpVPadIterator* AliMpMotifPosition::CreateIterator() const
{
// Iterator is not yet implemented.
//

  return new AliMpMotifPositionPadIterator(this);
}  

//______________________________________________________________________________
Bool_t AliMpMotifPosition::HasPad(const AliMpIntPair& indices) const
{
// Returns true if pad with the specified indices exists in 
// this motif position.
// ---

  if (!HasIndices(indices)) return kFALSE;
  
  if (fMotif->GetMotifType()->IsFull()) return kTRUE;
  
  return fMotif->GetMotifType()->HasPad(indices-GetLowIndicesLimit());
}

