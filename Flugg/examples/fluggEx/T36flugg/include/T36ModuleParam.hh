// $Id$
// Flugg tag $Name$

//
//  Parameterisation for t36 module
//

#ifndef T36ModuleParam_H
#define T36ModuleParam_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;
class G4Box;

class T36ModuleParam : public G4VPVParameterisation
{ 
  public:
    T36ModuleParam(G4double ModuleSizeZ, G4double NbOfModules); 
    ~T36ModuleParam();
    void ComputeTransformation
    (const G4int copyNo,G4VPhysicalVolume *physVol) const;

  private:

    G4double fModuleSizeZ;  
    G4double fNbOfModules;      
};

#endif


