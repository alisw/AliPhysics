// $Id$
// Category: sector
//
// Class AliMpZone
// ---------------
// Class describing a zone composed of the zone segments.
// The zone contains pads of the same dimensions.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpZone.h"
#include "AliMpSubZone.h"

ClassImp(AliMpZone)

//_____________________________________________________________________________
AliMpZone::AliMpZone(Int_t id) 
  : TObject(),
    fID(id),
    fPadDimensions(TVector2())
{
//
}

//_____________________________________________________________________________
AliMpZone::AliMpZone() 
  : TObject(),
    fID(0),
    fPadDimensions(TVector2())
{
//
}

//_____________________________________________________________________________
AliMpZone::~AliMpZone() {
//

  for (Int_t i=0; i<GetNofSubZones(); i++)
    delete fSubZones[i];  
}

//
// public methods
//

//_____________________________________________________________________________
void AliMpZone::AddSubZone(AliMpSubZone* subZone)
{
// Adds row segment.
// ---

  fSubZones.push_back(subZone);
}  
  
//_____________________________________________________________________________
AliMpSubZone* AliMpZone::FindSubZone(AliMpVMotif* motif) const
{
// Finds a subzone with a specified motif;
// returns 0 if not found.
// ---

  for (Int_t i=0; i<GetNofSubZones(); i++) {
    AliMpSubZone* subZone = GetSubZone(i);
    if (subZone->GetMotif() == motif) return subZone;
  }
  
  return 0;  
}

//_____________________________________________________________________________
Int_t AliMpZone::GetNofSubZones() const 
{
// Returns number of row segments.

  return fSubZones.size();
}  

//_____________________________________________________________________________
AliMpSubZone* AliMpZone::GetSubZone(Int_t i) const 
{
  if (i<0 || i>=GetNofSubZones()) {
    Warning("GetSubZone", "Index outside range");
    return 0;
  }
  
  return fSubZones[i];  
}
