// $Id$
// Category: geometry
//
// Geometry services
// The class provides service methods for accessing to Geant4 geometry,
// namely using AliMC volumes and materials identifiers. 

#ifndef TG4_GEOMETRY_SERVICES_H
#define TG4_GEOMETRY_SERVICES_H

#include "TG4NameMap.h"
#include "TG4Globals.h"

#include <globals.hh>

#include <Rtypes.h>

class TG4CutVector;
class TG4FlagVector;
class TG4GeometryOutputManager;

class G4Material;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4UserLimits;

class TG4GeometryServices
{
  public:
    TG4GeometryServices(TG4intVector* mediumIdVector, TG4NameMap* nameMap);
    // --> protected
    // TG4GeometryServices();
    // TG4GeometryServices(const TG4GeometryServices& right);
    virtual ~TG4GeometryServices();

    // static access method
    static TG4GeometryServices* Instance();

    // methods
    G4double* CreateG4doubleArray(Float_t* array, G4int size) const;
    G4String  CutName(const char* name) const;
    void G4ToG3VolumeName(G4String& name) const;
    G4int SetUserLimits(G4UserLimits* userLimits, G4LogicalVolume* lv);		     
    G4Material* MixMaterials(G4String name, G4double density,
                     TG4StringVector* matNames, TG4doubleVector* matWeights);

    // get methods
    Int_t NofG3Volumes() const; 
    Int_t NofG4LogicalVolumes() const; 
    Int_t NofG4PhysicalVolumes() const; 
    Int_t NofSensitiveDetectors() const; 

    G4int GetVolumeID(const G4String& volumeName) const;
    G4int GetVolumeID(G4LogicalVolume* volume) const;
    G4String GetVolumeName(G4int volumeId) const;
    G4LogicalVolume* GetLogicalVolume(G4int volumeId) const;
    G4bool IsG3Volume(const G4String& lvName) const;
    const G4String& GetMapSecond(const G4String& name);

          // materials
    G4int GetMediumId(G4Material* material) const;    
    G4int GetMediumId(G4int volumeId) const;    
    G4double GetEffA(G4Material* material) const;
    G4double GetEffZ(G4Material* material) const;

  protected:
    TG4GeometryServices();
    TG4GeometryServices(const TG4GeometryServices& right);

    // operators
    TG4GeometryServices& operator=(const TG4GeometryServices& right);

  private:        
    // static data members
    static TG4GeometryServices*  fgInstance;    //this instance

    // data members
    TG4intVector*  fMediumIdVector;  //vector of second indexes for materials
    TG4NameMap*    fNameMap;         //map of volumes names to modules names
};

// inline methods
inline TG4GeometryServices* TG4GeometryServices::Instance()
{ return fgInstance; }

#endif //TG4_GEOMETRY_SERVICES_H

