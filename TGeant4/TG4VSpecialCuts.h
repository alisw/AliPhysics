// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4VSpecialCuts
// ---------------------
// Abstract base class for a special process that activates 
// kinetic energy cuts.
// The pure virtual functions GetMinEkine have to be implemented
// by derived classes specific for each particle type
// (see TG4G3ParticleWSP.h).

#ifndef TG4_V_SPECIAL_CUTS_H
#define TG4_V_SPECIAL_CUTS_H

#include <G4VProcess.hh>

class TG4G3CutVector;
class TG4Limits;

class G4Track;

class TG4VSpecialCuts: public G4VProcess
{
  public:
    TG4VSpecialCuts(const G4String& processName);
    // --> protected
    // TG4VSpecialCuts();		   
    virtual ~TG4VSpecialCuts();

    // methods
    virtual G4double GetMinEkine(const TG4Limits& limits,
                                 const G4Track& track) const = 0;
    
    virtual G4double PostStepGetPhysicalInteractionLength(
                         const G4Track& track, G4double previousStepSize,
                         G4ForceCondition* condition);

    virtual G4VParticleChange* PostStepDoIt(const G4Track& track, 
                         const G4Step& step);
			    
    virtual G4double AlongStepGetPhysicalInteractionLength(
                         const G4Track&, G4double, G4double, G4double&,
                         G4GPILSelection*)
			 { return -1.0; }

    virtual G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&)
			 { return 0; }

    virtual G4double AtRestGetPhysicalInteractionLength(const G4Track&,
		         G4ForceCondition* )
			 { return -1.0; }
			    
    virtual G4VParticleChange* AtRestDoIt(
			 const G4Track&, const G4Step&)
			 { return 0; }

  protected:
    TG4VSpecialCuts();		   
    TG4VSpecialCuts(const TG4VSpecialCuts& right);
    
    // operators
    TG4VSpecialCuts& operator = (const TG4VSpecialCuts& right);
};

#endif //TG4_SPECIAL_CUTS_H



