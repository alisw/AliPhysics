// $Id$
// Category: event
//
// Class for detecting and stopping looping particles
// or particles that reached maximal number of steps.

#ifndef ALI_STEPPING_ACTION_H
#define ALI_STEPPING_ACTION_H

#include "TG4SteppingAction.h"

#include <G4ThreeVector.hh>
#include <globals.hh>

class AliSteppingActionMessenger;

class AliSteppingAction : public TG4SteppingAction
{
  enum { 
    kCheckNofSteps = 100,
  };

  public:
    AliSteppingAction();
    // protected
    // AliSteppingAction(const AliSteppingAction& right);
    virtual ~AliSteppingAction();

    // methods
    virtual void SteppingAction(const G4Step* step);
    
  protected:
    AliSteppingAction(const AliSteppingAction& right);

    // operators
    AliSteppingAction& operator=(const AliSteppingAction& right);

  private:
    // static data members
    static const G4double fgkTolerance; //tolerance used in detecting 
                                        //of looping particles

    // data members
    G4ThreeVector  fKeptStepPoint;           //kept step point
    AliSteppingActionMessenger*  fMessenger; //messenger
};

#endif //ALI_STEPPING_ACTION_H

