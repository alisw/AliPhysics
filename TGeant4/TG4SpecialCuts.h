// $Id$
// Category: physics
//
// Special process that activates kinetic energy cuts

#ifndef TG4_SPECIAL_CUTS_H
#define TG4_SPECIAL_CUTS_H

#include <G4UserSpecialCuts.hh>
#include "TG3ParticleWSP.h"

class TG4CutVector;
class TG4Limits;

typedef G4double(TG4CutVector::*PtrMinEkineInCutVector)(const G4Track&) const;
typedef G4double(TG4Limits::*PtrMinEkineInLimits)(const G4Track&) const;

class TG4SpecialCuts: public G4UserSpecialCuts
// to do: change to inheritance from G4VProcess
{
  public:
    TG4SpecialCuts(TG3ParticleWSP particle, TG4CutVector* cutVector, 
                   const G4String& processName ="specialCut");
    // --> protected
    // TG4SpecialCuts();		   
    // TG4SpecialCuts(const TG4SpecialCuts& right);
    virtual ~TG4SpecialCuts();

    // methods
    virtual G4double PostStepGetPhysicalInteractionLength(
                       const G4Track& track, G4double previousStepSize,
                       G4ForceCondition* condition);
    virtual G4VParticleChange* PostStepDoIt(const G4Track& , const G4Step&);
			    
  protected:
    TG4SpecialCuts();		   
    TG4SpecialCuts(const TG4SpecialCuts& right);
    
    // operators
    TG4SpecialCuts& operator = (const TG4SpecialCuts& right);
    
  private:
    // data members
    TG4CutVector*           fCutVector;              //TG4CutVector 
    PtrMinEkineInCutVector  fPtrMinEkineInCutVector; //pointer to 
                                //TG4CutVector::GetMinEKineForXX() method for 
				//the particle XX that this process is applied to 
    PtrMinEkineInLimits     fPtrMinEkineInLimits;    //pointer to 
                                //TG4Limits::GetMinEKineForXX() method for 
				//the particle XX that this process is applied to 
};

#endif //TG4_SPECIAL_CUTS_H



