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

class TG4VSensitiveDetector;

class G4LogicalVolume;
class G4VSensitiveDetector;

class TG4SDServices
{
  public:
    TG4SDServices();
    // --> protected
    // TG4SDServices(const TG4SDServices& right);
    virtual ~TG4SDServices();

    // static methods
    static TG4SDServices* Instance();

    // methods
    void PrintStatistics(G4bool open, G4bool close) const;

    // get methods
          // volume IDs conversions
    G4int GetVolumeID(const G4String& volumeName) const;
    G4int GetVolumeID(G4LogicalVolume* volume) const;
    G4String         GetVolumeName(G4int volumeId) const;
    G4LogicalVolume* GetLogicalVolume(G4int volumeId) const;   
    G4int            GetMediumId(G4int volumeId) const; 
          // SDs
    Int_t NofSensitiveDetectors() const; 
    TG4VSensitiveDetector* GetSensitiveDetector(G4VSensitiveDetector* sd) const;  

  protected:
    TG4SDServices(const TG4SDServices& right);

    // operators
    TG4SDServices& operator=(const TG4SDServices& right);
  
    // static data members
    static TG4SDServices* fgInstance;   //this instance
};

// inline methods
inline TG4SDServices* TG4SDServices::Instance() 
{ return fgInstance; }

#endif //TG4_SD_SERVICES_H

