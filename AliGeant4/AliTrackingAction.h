// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliTrackingAction
// -----------------------
// Class that takes care of storing kinematics.

#ifndef ALI_TRACKING_ACTION_H
#define ALI_TRACKING_ACTION_H

#include "AliVerbose.h"
#include "AliTrackingActionMessenger.h"

#include "TG4TrackingAction.h"

#include <globals.hh>

class AliTrackInformation;

class G4Track;

class TClonesArray;

class AliTrackingAction : public TG4TrackingAction,
                          public AliVerbose
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
    void FinishPrimaryTrack();
    void SaveTrack(const G4Track* track);

    // set methods
    void SetNewVerboseLevel(G4int level);
    void SetNewVerboseTrackID(G4int trackID);
    void SetSavePrimaries(G4bool savePrimaries);

    // get methods
    G4bool GetSavePrimaries() const;
    G4int GetNofTracks() const;

  protected:
    AliTrackingAction(const AliTrackingAction& right);

    // operators
    AliTrackingAction& operator=(const AliTrackingAction& right);

  private:
    // methods
    G4int GetParticleIndex(G4int trackID);
    AliTrackInformation* GetTrackInformation(const G4Track* track,
                                             const G4String& method) const;
  
    // static data members
    static AliTrackingAction*   fgInstance; //this instance

    // data members
    AliTrackingActionMessenger  fMessenger; //messenger
    G4int   fPrimaryTrackID;    //current primary track ID 
    G4bool  fSavePrimaries;     //control of saving primaries
    G4int   fNewVerboseLevel;   //new /tracking/verbose level
    G4int   fNewVerboseTrackID; //track ID for which new /tracking/verbose level
                                // is applied
    G4int   fTrackCounter;      //tracks counter
};

// inline methods

inline AliTrackingAction* AliTrackingAction::Instance() 
{ return fgInstance; }

inline void AliTrackingAction::SetSavePrimaries(G4bool savePrimaries)
{ fSavePrimaries = savePrimaries; }

inline G4bool AliTrackingAction::GetSavePrimaries() const
{ return fSavePrimaries; }

inline G4int AliTrackingAction::GetNofTracks() const
{ return fTrackCounter; }

#endif //ALI_TRACKING_ACTION_H
