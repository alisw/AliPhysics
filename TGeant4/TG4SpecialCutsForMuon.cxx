// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4SpecialCutsForMuon
// ---------------------------
// See the class description in the header file.

#include "TG4SpecialCutsForMuon.h"
#include "TG4Limits.h"


//_____________________________________________________________________________
TG4SpecialCutsForMuon::TG4SpecialCutsForMuon(const G4String& processName)
  : TG4VSpecialCuts(processName) {
//
}

//_____________________________________________________________________________
TG4SpecialCutsForMuon::TG4SpecialCutsForMuon() {
//
}

//_____________________________________________________________________________
TG4SpecialCutsForMuon::~TG4SpecialCutsForMuon() {
//
}

// public methods

//_____________________________________________________________________________
G4double TG4SpecialCutsForMuon::GetMinEkine(const TG4Limits& limits,
                                            const G4Track& track) const
{					     
// Returns the min kinetic energy cut from limits.
// --- 

  return limits.GetMinEkineForMuon(track);
}  

