// $Id$
// Category: event
//
// Class that defines Alice stacking mechanism.

#ifndef ALI_STACKING_ACTION_H
#define ALI_STACKING_ACTION_H

#include <G4UserStackingAction.hh>
#include <globals.hh>

class AliStackingActionMessenger;
class AliTrackingAction;
class G4Track;

class AliStackingAction : public G4UserStackingAction
{
  public:
    AliStackingAction();
    // --> protected
    // AliStackingAction(const AliStackingAction& right);
    virtual ~AliStackingAction();

    // methods
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track);
    void NewStage();
    void ClearPrimaryStack();
    void PrepareNewEvent();

    // set methods
    void SetVerboseLevel(G4int level);

    // get methods
    G4int GetVerboseLevel() const;

  protected:
    AliStackingAction(const AliStackingAction& right);

    // operators
    AliStackingAction& operator=(const AliStackingAction& right);

  private:
    // data members
    G4int                        fStage;          //stage number
    G4int                        fVerboseLevel;   //verbose level
    G4bool                       fSavePrimaries;  //control of saving primaries
    G4TrackStack*                fPrimaryStack;   //stack of primary tracks
    AliTrackingAction*           fTrackingAction; //AliTrackingAction
    AliStackingActionMessenger*  fMessenger;      //messenger
};


// inline methods

inline void AliStackingAction::SetVerboseLevel(G4int level)
{ fVerboseLevel = level; }

inline G4int AliStackingAction::GetVerboseLevel() const
{ return fVerboseLevel; }

#endif //ALI_STACKING_ACTION_H

