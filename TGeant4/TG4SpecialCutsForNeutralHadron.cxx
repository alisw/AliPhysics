// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4SpecialCutsForNeutralHadron
// ------------------------------------
// See the class description in the header file.

#include "TG4SpecialCutsForNeutralHadron.h"
#include "TG4Limits.h"


//_____________________________________________________________________________
TG4SpecialCutsForNeutralHadron::TG4SpecialCutsForNeutralHadron(
                                              const G4String& processName)
  : TG4VSpecialCuts(processName) {
//
}

//_____________________________________________________________________________
TG4SpecialCutsForNeutralHadron::TG4SpecialCutsForNeutralHadron() {
//
}

//_____________________________________________________________________________
TG4SpecialCutsForNeutralHadron::~TG4SpecialCutsForNeutralHadron() {
//
}

// public methods

//_____________________________________________________________________________
G4double TG4SpecialCutsForNeutralHadron::GetMinEkine(const TG4Limits& limits,
                                                     const G4Track& track) const
{					     
// Returns the min kinetic energy cut from limits.
// --- 

  return limits.GetMinEkineForNeutralHadron(track);
}  

