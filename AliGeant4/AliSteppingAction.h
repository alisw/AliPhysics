// $Id$
// Category: event
//
// Author: I. Hrivnacova
//
// Class AliSteppingAction
// -----------------------
// Class takes care of stopping particles
// if they get outside of user defined tracking region
// (in AliRun).

#ifndef ALI_STEPPING_ACTION_H
#define ALI_STEPPING_ACTION_H

#include "AliSteppingActionMessenger.h"

#include "TG4SteppingAction.h"


class AliSteppingAction : public TG4SteppingAction
{
  public:
    AliSteppingAction();
    virtual ~AliSteppingAction();

    // methods
    virtual void SteppingAction(const G4Step* step);
    
  private:
    AliSteppingActionMessenger  fMessenger; //messenger
};

#endif //ALI_STEPPING_ACTION_H

