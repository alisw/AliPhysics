// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliTrackInformation
// -------------------------
// See the class description in the header file.

#include "AliTrackInformation.h"

G4Allocator<AliTrackInformation> gAliTrackInfoAllocator;

//_____________________________________________________________________________
AliTrackInformation::AliTrackInformation()
  : G4VUserTrackInformation(),
    fTrackParticleID(-1),
    fParentParticleID(-1) 
{
//
}

//_____________________________________________________________________________
AliTrackInformation::AliTrackInformation(G4int trackParticleID) 
  : G4VUserTrackInformation(),
    fTrackParticleID(trackParticleID),
    fParentParticleID(-1)
{
//
}    

//_____________________________________________________________________________
AliTrackInformation::AliTrackInformation(G4int trackParticleID, 
                                         G4int parentParticleID)
  : G4VUserTrackInformation(),
    fTrackParticleID(trackParticleID),
    fParentParticleID(parentParticleID)
{
//
}    

//_____________________________________________________________________________
AliTrackInformation::~AliTrackInformation(){
//
}    

// public methods

//_____________________________________________________________________________
void AliTrackInformation::Print() const
{
// Prints track information.
// ---

  G4cout << "TrackParticleID: " << fTrackParticleID << "   "
         << "ParentParticleID: " << fParentParticleID << G4endl;
}
