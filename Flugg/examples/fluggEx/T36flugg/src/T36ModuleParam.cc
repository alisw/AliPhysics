// $Id$
// Flugg tag $Name$

//
// T36ModuleParam.cc, 3/III/99, Sara Vanini
// parametrisation for t36 modules
//

#include "T36ModuleParam.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4RotationMatrix.hh"


T36ModuleParam::T36ModuleParam(G4double ModuleSizeZ, G4double NbOfModules)
{
  fModuleSizeZ = ModuleSizeZ; 
  fNbOfModules = NbOfModules;
}

T36ModuleParam::~T36ModuleParam()
{}

void T36ModuleParam::ComputeTransformation
(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  G4double      Zposition= - fModuleSizeZ/2*(fNbOfModules-1) + copyNo*fModuleSizeZ;
  G4ThreeVector origin(0,0,Zposition);
  physVol->SetTranslation(origin);
}
