// $Id$
// Category: event
//
// Class that ensures additional call to sensitive detector
// when track crosses geometrical boundary.

#ifndef TG4_STEPPING_ACTION_H
#define TG4_STEPPING_ACTION_H

#include <G4UserSteppingAction.hh>
#include <globals.hh>

class G4Step;

class TG4SteppingAction : public G4UserSteppingAction 
{
  public:
    TG4SteppingAction();
    // --> protected
    // TG4SteppingAction(const TG4SteppingAction& right);
    virtual ~TG4SteppingAction();
   
    // methods
    virtual void SteppingAction(const G4Step* step) {;}
                  // the following method should not
		  // be overwritten in a derived class
    virtual void UserSteppingAction(const G4Step* step);


  protected:
    TG4SteppingAction(const TG4SteppingAction& right);

    // operators
    TG4SteppingAction& operator=(const TG4SteppingAction& right);
};

#endif //TG4_STEPPING_ACTION_H
