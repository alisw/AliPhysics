// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4SpecialCutsForChargedHadron
// ------------------------------------
// See the class description in the header file.

#include "TG4SpecialCutsForChargedHadron.h"
#include "TG4Limits.h"


//_____________________________________________________________________________
TG4SpecialCutsForChargedHadron::TG4SpecialCutsForChargedHadron(
					  const G4String& processName)
  : TG4VSpecialCuts(processName) {
//
}

//_____________________________________________________________________________
TG4SpecialCutsForChargedHadron::TG4SpecialCutsForChargedHadron() {
//
}

//_____________________________________________________________________________
TG4SpecialCutsForChargedHadron::~TG4SpecialCutsForChargedHadron() {
//
}

// public methods

//_____________________________________________________________________________
G4double TG4SpecialCutsForChargedHadron::GetMinEkine(const TG4Limits& limits,
                                                     const G4Track& track) const
{					     
// Returns the min kinetic energy cut from limits.
// --- 

  return limits.GetMinEkineForChargedHadron(track);
}  

