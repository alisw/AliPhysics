// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4SpecialCutsForOther
// ----------------------------
// See the class description in the header file.

#include "TG4SpecialCutsForOther.h"
#include "TG4Limits.h"


//_____________________________________________________________________________
TG4SpecialCutsForOther::TG4SpecialCutsForOther(const G4String& processName)
  : TG4VSpecialCuts(processName) {
//
}

//_____________________________________________________________________________
TG4SpecialCutsForOther::TG4SpecialCutsForOther() {
//
}

//_____________________________________________________________________________
TG4SpecialCutsForOther::~TG4SpecialCutsForOther() {
//
}

// public methods

//_____________________________________________________________________________
G4double TG4SpecialCutsForOther::GetMinEkine(const TG4Limits& limits,
                                            const G4Track& track) const
{					     
// Returns the min kinetic energy cut from limits.
// --- 

  return limits.GetMinEkineForOther(track);
}  

