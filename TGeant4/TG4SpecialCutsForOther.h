// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4SpecialCutsForOther
// ----------------------------
// Special process that activates the kinetic energy cuts
// for other (than specified in TParticleWSP) particles.
 
#ifndef TG4_SPECIAL_CUTS_FOR_OTHER_H
#define TG4_SPECIAL_CUTS_FOR_OTHER_H

#include "TG4VSpecialCuts.h"

class TG4Limits;

class G4Track;

class TG4SpecialCutsForOther: public TG4VSpecialCuts
{
  public:
    TG4SpecialCutsForOther(const G4String& processName = "specialCutForOther");
    // --> protected
    // TG4SpecialCutsForOther();		   
    virtual ~TG4SpecialCutsForOther();

    // methods
    virtual G4double GetMinEkine(const TG4Limits& limits,
                                 const G4Track& track) const;
			    
  protected:
    TG4SpecialCutsForOther();		   
};

#endif //TG4_SPECIAL_CUTS_FOR_OTHER_H



