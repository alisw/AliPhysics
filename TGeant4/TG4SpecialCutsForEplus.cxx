// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4SpecialCutsForEplus
// ----------------------------
// See the class description in the header file.

#include "TG4SpecialCutsForEplus.h"
#include "TG4Limits.h"


//_____________________________________________________________________________
TG4SpecialCutsForEplus::TG4SpecialCutsForEplus(const G4String& processName)
  : TG4VSpecialCuts(processName) {
//
}

//_____________________________________________________________________________
TG4SpecialCutsForEplus::TG4SpecialCutsForEplus() {
//
}

//_____________________________________________________________________________
TG4SpecialCutsForEplus::~TG4SpecialCutsForEplus() {
//
}

// public methods

//_____________________________________________________________________________
G4double TG4SpecialCutsForEplus::GetMinEkine(const TG4Limits& limits,
                                             const G4Track& track) const
{					     
// Returns the min kinetic energy cut from limits.
// --- 

  return limits.GetMinEkineForEplus(track);
}  

