// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliStackingAction
// -----------------------
// Class that defines AliRoot specific stacking mechanism.

#ifndef ALI_STACKING_ACTION_H
#define ALI_STACKING_ACTION_H

#include "AliVerbose.h"

#include <G4UserStackingAction.hh>
#include <globals.hh>

class AliTrackingAction;

class G4Track;
class G4TrackStack;

class AliStackingAction : public G4UserStackingAction,
                          public AliVerbose
{
  public:
    AliStackingAction();
    // --> protected
    // AliStackingAction(const AliStackingAction& right);
    virtual ~AliStackingAction();

    // methods
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track);
    void NewStage();
    void PrepareNewEvent();

  protected:
    AliStackingAction(const AliStackingAction& right);

    // operators
    AliStackingAction& operator=(const AliStackingAction& right);

  private:
    // data members
    G4int                        fStage;          //stage number
    G4bool                       fSavePrimaries;  //control of saving primaries
    G4TrackStack*                fPrimaryStack;   //stack of primary tracks
    AliTrackingAction*           fTrackingAction; //AliTrackingAction
};

#endif //ALI_STACKING_ACTION_H

