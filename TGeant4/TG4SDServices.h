// $Id$
// Category: digits+hits
//
// Sensitive detectors services
// The class provides service methods for accessing to Geant4 geometry,
// namely using AliMC volumes identifiers
// (implemented via TG4VSensitiveDetector instances).

#ifndef TG4_SD_SERVICES_H
#define TG4_SD_SERVICES_H

#include <globals.hh>

#include <Rtypes.h>

class G4LogicalVolume;

class TG4SDServices
{
  public:
    TG4SDServices();
    virtual ~TG4SDServices();

    // static methods
    static TG4SDServices* Instance();

    // get methods
    G4int GetVolumeID(const G4String& volumeName) const;
    G4int GetVolumeID(G4LogicalVolume* volume) const;
    G4String GetVolumeName(G4int volumeId) const;
    G4LogicalVolume* GetLogicalVolume(G4int volumeId) const;   
    Int_t NofSensitiveDetectors() const; 
    G4int GetMediumId(G4int volumeId) const;    

  protected:
    // static data members
    static TG4SDServices* fgInstance;   //this instance
};

// inline methods
inline TG4SDServices* TG4SDServices::Instance() 
{ return fgInstance; }

#endif //TG4_SD_SERVICES_H

