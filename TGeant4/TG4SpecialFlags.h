// $Id$
// Category: physics
//
// Special process that applies process control flags 

#ifndef TG4_SPECIAL_FLAGS_H
#define TG4_SPECIAL_FLAGS_H

#include "TG4Globals.h"

#include <G4VProcess.hh>
#include <G4ProcessVector.hh>
#include <globals.hh>

class TG4SpecialFlags : public G4VProcess 
{
  enum Switch { kSwitch, kReswitch, kUnswitch };

  public:     
    TG4SpecialFlags(const G4String& processName ="specialFlag" );
    // --> protected
    // TG4SpecialFlags(const TG4SpecialFlags& right);
    virtual ~TG4SpecialFlags();

    // methods

    virtual G4double PostStepGetPhysicalInteractionLength(
                         const G4Track& track, G4double previousStepSize, 
		         G4ForceCondition* condition);

    virtual G4VParticleChange* PostStepDoIt(const G4Track& track, 
                                   const G4Step& step);

    virtual G4double AlongStepGetPhysicalInteractionLength(
                         const G4Track& track, G4double previousStepSize ,
			 G4double currentMinimumStep, G4double& proposedSafety,
                         G4GPILSelection* selection)
			 { return -1.0; }

    virtual G4VParticleChange* AlongStepDoIt(const G4Track& ,
			           const G4Step& step)
	                 { return 0; }

    virtual G4double AtRestGetPhysicalInteractionLength(
                         const G4Track& track, G4ForceCondition* condition)
                         { return -1.0; }

    virtual G4VParticleChange* AtRestDoIt(const G4Track& track,
			           const G4Step& step)
                         { return 0; }

  protected:
    TG4SpecialFlags(const TG4SpecialFlags& right);
    
    //operators
    TG4SpecialFlags& operator = (const TG4SpecialFlags& right);
    
  private:  
    // data members
    Switch           fSwitchFlags;       //directive passed from PostStepGetPIL
                                         //to PostStepDoIt
    G4ProcessVector  fSwitchedProcesses; //vector of the processes activation of
                                         //which is changed by this process
    TG4boolVector    fSwitchedFlags;     //vector for storing the current values of 
                                         //the processes activation   
};

#endif //TG4_SPECIAL_FLAGS_H

