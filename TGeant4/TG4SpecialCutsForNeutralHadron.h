// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4SpecialCutsForNeutralHadron
// ------------------------------------
// Special process that activates kinetic energy cuts
// for neutral hadrons.

#ifndef TG4_SPECIAL_CUTS_FOR_NEUTRAL_HADRON_H
#define TG4_SPECIAL_CUTS_FOR_NEUTRAL_HADRON_H

#include "TG4VSpecialCuts.h"

class TG4Limits;

class G4Track;

class TG4SpecialCutsForNeutralHadron: public TG4VSpecialCuts
{
  public:
    TG4SpecialCutsForNeutralHadron(const G4String& processName 
                                                   = "specialCutForNeuHadron");
    // --> protected
    // TG4SpecialCutsForNeutralHadron();		   
    virtual ~TG4SpecialCutsForNeutralHadron();

    // methods
    virtual G4double GetMinEkine(const TG4Limits& limits,
                                 const G4Track& track) const;
			    
  protected:
    TG4SpecialCutsForNeutralHadron();		   
};

#endif //TG4_SPECIAL_CUTS_FOR_NEUTRAL_HADRON_H



