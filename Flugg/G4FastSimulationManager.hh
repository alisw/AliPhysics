//---------------------------------------------------------------
//
//  dummy G4FastSimulationManager.hh
//
//---------------------------------------------------------------


#ifndef G4FastSimulationManager_h
#define G4FastSimulationManager_h 1

#include "G4LogicalVolume.hh"

class G4FastSimulationManager
{
public:  
  G4FastSimulationManager(G4LogicalVolume* logVol);
  ~G4FastSimulationManager();
};

#endif
