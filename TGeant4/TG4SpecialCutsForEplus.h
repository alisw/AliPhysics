// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4SpecialCutsForEplus
// ----------------------------
// Special process that activates kinetic energy cuts
// for e+.

#ifndef TG4_SPECIAL_CUTS_FOR_EPLUS_H
#define TG4_SPECIAL_CUTS_FOR_EPLUS_H

#include "TG4VSpecialCuts.h"

class TG4Limits;

class G4Track;

class TG4SpecialCutsForEplus: public TG4VSpecialCuts
{
  public:
    TG4SpecialCutsForEplus(const G4String& processName = "specialCutForEplus");
    // --> protected
    // TG4SpecialCutsForEplus();		   
    virtual ~TG4SpecialCutsForEplus();

    // methods
    virtual G4double GetMinEkine(const TG4Limits& limits,
                                 const G4Track& track) const;
			    
  protected:
    TG4SpecialCutsForEplus();		   
};

#endif //TG4_SPECIAL_CUTS_FOR_EPLUS_H



