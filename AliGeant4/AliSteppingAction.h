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

class G4Track;

class AliSteppingAction : public TG4SteppingAction
{
  enum { 
    kCheckNofSteps = 100,
    kMaxNofLoopSteps = 5,
    kMaxNofSteps = 5000
  };

  public:
    AliSteppingAction();
    // protected
    // AliSteppingAction(const AliSteppingAction& right);
    virtual ~AliSteppingAction();

    // methods
    virtual void SteppingAction(const G4Step* step);
    
    // set methods
    void SetLoopVerboseLevel(G4int level);
    
    // get methods
    G4int GetLoopVerboseLevel() const;

  protected:
    AliSteppingAction(const AliSteppingAction& right);

    // operators
    AliSteppingAction& operator=(const AliSteppingAction& right);

  private:
    // methods
    void PrintTrackInfo(const G4Track* track) const;
  
    // static data members
    static const G4double fgkTolerance;    //tolerance used in detecting 
                                           //of looping particles

    // data members
    G4ThreeVector   fKeptStepPoint;        //kept step point
    G4int           fLoopVerboseLevel;     //tracking verbose level
                                           //for looping particles
    G4int           fStandardVerboseLevel; //standard tracking verbose level
    G4int           fLoopStepCounter;      //loop steps counter
    AliSteppingActionMessenger*  fMessenger;  //messenger
};

// inline methods

inline void AliSteppingAction::SetLoopVerboseLevel(G4int level)
{ fLoopVerboseLevel = level; }

inline G4int AliSteppingAction::GetLoopVerboseLevel() const
{ return fLoopVerboseLevel; }

#endif //ALI_STEPPING_ACTION_H

