// $Id$
//
// Author: Ivana Hrivnacova, IPN Orsay

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

//______________________________________________________________________________
AliMpMotifPosition::~AliMpMotifPosition(){
// 
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

