// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4SpecialControls
// ----------------------------
// Special process that applies the G3 process controls.

#ifndef TG4_SPECIAL_CONTROLS_H
#define TG4_SPECIAL_CONTROLS_H

#include "TG4Globals.h"
#include "TG4G3Control.h"

#include <G4VProcess.hh>
#include <G4ProcessVector.hh>
#include <globals.hh>
#include <g4std/vector>

class TG4SpecialControls : public G4VProcess 
{
  enum Switch { kSwitch, kReswitch, kUnswitch };

  public:     
    TG4SpecialControls(const G4String& processName ="specialControl" );
    // --> protected
    // TG4SpecialControls(const TG4SpecialControls& right);
    virtual ~TG4SpecialControls();

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
    TG4SpecialControls(const TG4SpecialControls& right);
    
    //operators
    TG4SpecialControls& operator = (const TG4SpecialControls& right);
    
  private:  
    // methods
    virtual void Reset();    			 
    
    // data members
    Switch           fSwitchControls;    //directive passed from PostStepGetPIL
                                         //to PostStepDoIt
    G4ProcessVector  fSwitchedProcesses; //vector of the processes the activation of
                                         //which is changed by this process
    TG4boolVector    fSwitchedControls;  //vector for storing the current values of 
                                         //the processes activation   
    G4int            fLastTrackID;       //the last processed track ID
};

#endif //TG4_SPECIAL_CONTROLS_H

