// $Id$
// Category: event
//
// Class that ensures additional call to sensitive detector
// when track crosses geometrical boundary.

#ifndef TG4_STEPPING_ACTION_H
#define TG4_STEPPING_ACTION_H

#include <G4UserSteppingAction.hh>

#include <globals.hh>

class G4Track;
class G4Step;

class TG4SteppingAction : public G4UserSteppingAction 
{
  enum { 
    kMaxNofSteps = 10000,
    kMaxNofLoopSteps = 5
  };

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

    // set methods
    void SetLoopVerboseLevel(G4int level);
    void SetMaxNofSteps(G4int number);

    // get methods
    G4int GetLoopVerboseLevel() const;
    G4int GetMaxNofSteps() const;

  protected:
    TG4SteppingAction(const TG4SteppingAction& right);

    // operators
    TG4SteppingAction& operator=(const TG4SteppingAction& right);
    
    // methods
    void PrintTrackInfo(const G4Track* track) const;

    // data members
    G4int  fMaxNofSteps;          //max number of steps allowed
    G4int  fStandardVerboseLevel; //standard tracking verbose level
    G4int  fLoopVerboseLevel;     //tracking verbose level                                           //for looping particles
                                  //for looping particles
    G4int  fLoopStepCounter;      //loop steps counter    
};

// inline methods

inline void TG4SteppingAction::SetLoopVerboseLevel(G4int level)
{ fLoopVerboseLevel = level; }

inline void TG4SteppingAction::SetMaxNofSteps(G4int number)
{ fMaxNofSteps = number; }

inline G4int TG4SteppingAction::GetMaxNofSteps() const
{ return fMaxNofSteps; }

inline G4int TG4SteppingAction::GetLoopVerboseLevel() const
{ return fLoopVerboseLevel; }

#endif //TG4_STEPPING_ACTION_H
