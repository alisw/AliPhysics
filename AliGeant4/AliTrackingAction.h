// $Id$
// Category: event
//
// Class that takes care of storing kinematics.

#ifndef ALI_TRACKING_ACTION_H
#define ALI_TRACKING_ACTION_H

#include "TG4TrackingAction.h"

#include <globals.hh>

class AliTrackingActionMessenger;

class G4Track;

class TClonesArray;

class AliTrackingAction : public TG4TrackingAction 
{
  public:
    AliTrackingAction();
    // --> protected
    // AliTrackingAction(const AliTrackingAction& right);
    virtual ~AliTrackingAction();
   
    // static get method
    static AliTrackingAction* Instance();

    // methods
    void PrepareNewEvent();
    virtual void PreTrackingAction(const G4Track* aTrack);
    virtual void PostTrackingAction(const G4Track* aTrack);
    void SaveParticle(const G4Track* track, G4String processName);
    void SaveAndDestroyTrack();

    // set methods
    void SetVerboseLevel(G4int level);
    void SetSavePrimaries(G4bool savePrimaries);

    // get methods
    G4int GetVerboseLevel() const;
    G4bool GetSavePrimaries() const;
    G4int GetNofTracks() const;
    G4int GetNofPrimaryTracks() const;
    G4int GetNofSavedTracks() const;

  protected:
    AliTrackingAction(const AliTrackingAction& right);

    // operators
    AliTrackingAction& operator=(const AliTrackingAction& right);

  private:
    // methods
    G4int GetParticleIndex(G4int trackID);
  
    // static data members
    static AliTrackingAction*    fgInstance; //this instance

    // data members
    TClonesArray*  fParticles;         //AliRun::fParticles
    G4int          fPrimaryTrackID;    //primary track ID 
    G4bool         fSavePrimaries;     //control of saving primaries
    G4int          fVerboseLevel;      //verbose level
    G4int          fPrimariesCounter;  //primary particles counter
    G4int          fParticlesCounter;  //particles counter
    G4int          fTrackCounter;      //tracks counter
    G4int          fLastParticleIndex; //index of the last particle in fParticles
    AliTrackingActionMessenger*  fMessenger; //messenger
};

// inline methods

inline void AliTrackingAction::SetVerboseLevel(G4int level)
{ fVerboseLevel = level; }

inline void AliTrackingAction::SetSavePrimaries(G4bool savePrimaries)
{ fSavePrimaries = savePrimaries; }

inline G4int AliTrackingAction::GetVerboseLevel() const
{ return fVerboseLevel; }

inline G4bool AliTrackingAction::GetSavePrimaries() const
{ return fSavePrimaries; }

inline G4int AliTrackingAction::GetNofTracks() const
{ return fTrackCounter; }

inline G4int AliTrackingAction::GetNofPrimaryTracks() const
{ return fPrimariesCounter; }

inline G4int AliTrackingAction::GetNofSavedTracks() const
{ return fParticlesCounter; }

#endif //ALI_TRACKING_ACTION_H
