// $Id$
// Category: sector
//
// Class AliMpZone
// ---------------
// Class describing a zone composed of the zone segments.
// The zone contains pads of the same dimensions.
// Included in AliRoot: 2003/05/02
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

#ifdef WITH_STL
  fSubZones.push_back(subZone);
#endif

#ifdef WITH_ROOT
  fSubZones.Add(subZone);
#endif
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
// ---

#ifdef WITH_STL
  return fSubZones.size();
#endif

#ifdef WITH_ROOT
  return fSubZones.GetEntriesFast();
#endif
}  

//_____________________________________________________________________________
AliMpSubZone* AliMpZone::GetSubZone(Int_t i) const 
{
// Returns i-th sub zone.
// ---

  if (i<0 || i>=GetNofSubZones()) {
    Warning("GetSubZone", "Index outside range");
    return 0;
  }
  
#ifdef WITH_STL
  return fSubZones[i];  
#endif

#ifdef WITH_ROOT
  return (AliMpSubZone*)fSubZones[i];  
#endif
}
