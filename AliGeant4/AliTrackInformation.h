// $Id$
// Category: event
//
// Class with additional track information.

#ifndef ALI_TRACK_INFORMATION_H
#define ALI_TRACK_INFORMATION_H

#include <G4VUserTrackInformation.hh>
#include <G4Allocator.hh>
#include <globals.hh>

class AliTrackInformation : public G4VUserTrackInformation
{
  public:
    AliTrackInformation();
    AliTrackInformation(G4int trackParticleID);
    AliTrackInformation(G4int trackParticleID, G4int parentParticleID);
    virtual ~AliTrackInformation();
   
    // operators
    // required by G4
    inline void *operator new(size_t);
      // Override "new" for "G4Allocator".
    inline void operator delete(void *trackInformation);
      // Override "delete" for "G4Allocator".
      
    // methods
    virtual void Print() const;  

    // set methods
    void SetTrackParticleID(G4int trackParticleID);
    void SetParentParticleID(G4int parentParticleID);

    // get methods
    G4int GetTrackParticleID() const;
    G4int GetParentParticleID() const;

  private:
    // data members
    G4int  fTrackParticleID; //the index of track particle in AliRoot stack
    G4int  fParentParticleID;//the index of parent track in AliRoot stack
};

// inline methods
#include "AliTrackInformation.icc"

#endif //ALI_TRACK_INFORMATION_H
