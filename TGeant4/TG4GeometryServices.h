// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class TG4GeometryServices
// -------------------------
// The class provides service methods for accessing to Geant4 geometry,
// namely using AliMC volumes and materials identifiers. 

#ifndef TG4_GEOMETRY_SERVICES_H
#define TG4_GEOMETRY_SERVICES_H

#include "TG4Verbose.h"
#include "TG4Globals.h"

#include <globals.hh>

#include <Rtypes.h>

class TG4IntMap;
class TG4NameMap;
class TG4Limits;
class TG4G3ControlVector;

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UserLimits;

class TG4GeometryServices : public TG4Verbose
{
  public:
    TG4GeometryServices(TG4IntMap* mediumMap, TG4NameMap* nameMap);
    // --> protected
    // TG4GeometryServices();
    // TG4GeometryServices(const TG4GeometryServices& right);
    virtual ~TG4GeometryServices();

    // static access method
    static TG4GeometryServices* Instance();

    // methods
           // utilities  
    G4double* CreateG4doubleArray(Float_t* array, G4int size) const;
    G4String  CutName(const char* name) const;
    G4String  CutMaterialName(const char* name) const;
    G4String  G4ToG3VolumeName(const G4String& name) const;
    G4String  GenerateLimitsName(G4int id, const G4String& medName,
                                           const G4String& matName) const;

    G4Material* MixMaterials(G4String name, G4double density,
                             const TG4StringVector& matNames, 
			     const TG4doubleVector& matWeights);
           // printing 
    void PrintNameMap() const;
    void PrintLimits(const G4String& name) const;
    void PrintVolumeLimits(const G4String& volumeName) const;
    void PrintStatistics(G4bool open, G4bool close) const;
    void PrintLogicalVolumeStore() const;

    // set methods
    void SetWorld(G4VPhysicalVolume* world);

    // get methods
           // volumes
    Int_t NofG3Volumes() const; 
    Int_t NofG4LogicalVolumes() const; 
    Int_t NofG4PhysicalVolumes() const; 
    G4bool IsSpecialControls() const;
    G4VPhysicalVolume* GetWorld() const;

    TG4Limits* GetLimits(G4UserLimits* limits) const;
    const G4String& GetMapSecond(const G4String& name);

    G4LogicalVolume* FindLogicalVolume(const G4String& name, 
                                       G4bool silent = false) const;
    TG4Limits*       FindLimits(const G4String& name, 
                                       G4bool silent = false) const;

          // materials
    G4int    GetMediumId(G4LogicalVolume* lv) const;    
    G4double GetEffA(G4Material* material) const;
    G4double GetEffZ(G4Material* material) const;
    G4Material* FindMaterial(G4double a, G4double z, G4double density) const;
    G4Material* FindMaterial(G4double* a, G4double* z, G4double density, 
                             G4int nmat, G4double* wmat) const;

  protected:
    TG4GeometryServices();
    TG4GeometryServices(const TG4GeometryServices& right);

    // operators
    TG4GeometryServices& operator=(const TG4GeometryServices& right);

  private:
    // methods        
    G4bool IsG3Volume(const G4String& lvName) const;
    G4bool CompareElement(G4double a, G4double z, const G4Element* elem) const;
    G4bool CompareMaterial(G4int nofElements, G4double density, 
                           const G4Material* material) const;
    G4double* ConvertAtomWeight(G4int nmat, G4double* a, G4double* wmat) const;

    // static data members
    static TG4GeometryServices*  fgInstance;   //this instance
    static const G4double  fgkAZTolerance;     //A,Z tolerance
    static const G4double  fgkDensityTolerance;//density tolerance (percentual)
 
    // data members
    TG4IntMap*         fMediumMap; //map of volumes names to medias IDs
    TG4NameMap*        fNameMap;   //map of volumes names to modules names
    G4VPhysicalVolume* fWorld;     //top pgysical volume (world)
};

// inline methods
inline TG4GeometryServices* TG4GeometryServices::Instance()
{ return fgInstance; }

inline void TG4GeometryServices::SetWorld(G4VPhysicalVolume* world)
{ fWorld = world; }

inline G4VPhysicalVolume* TG4GeometryServices::GetWorld() const
{ return fWorld; }

#endif //TG4_GEOMETRY_SERVICES_H

