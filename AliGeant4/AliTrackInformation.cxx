// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliTrackInformation.h"

G4Allocator<AliTrackInformation> trackInfoAllocator;

AliTrackInformation::AliTrackInformation()
  : G4VUserTrackInformation(),
    fTrackParticleID(-1),
    fParentParticleID(-1) 
{
//
}

AliTrackInformation::AliTrackInformation(G4int trackParticleID) 
  : G4VUserTrackInformation(),
    fTrackParticleID(trackParticleID),
    fParentParticleID(-1)
{
//
}    

AliTrackInformation::AliTrackInformation(G4int trackParticleID, 
                                         G4int parentParticleID)
  : G4VUserTrackInformation(),
    fTrackParticleID(trackParticleID),
    fParentParticleID(parentParticleID)
{
//
}    

AliTrackInformation::~AliTrackInformation(){
//
}    

// public methods

void AliTrackInformation::Print() const
{
// Prints track information.
// ---

  G4cout << "TrackParticleID: " << fTrackParticleID << "   "
         << "ParentParticleID: " << fParentParticleID << G4endl;
}
