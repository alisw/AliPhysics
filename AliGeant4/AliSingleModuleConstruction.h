// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliSingleModuleConstruction
// ---------------------------------
// Class for geometry construction of a standalone
// AliRoot module. 

#ifndef ALI_SINGLE_MODULE_CONSTRUCTION_H
#define ALI_SINGLE_MODULE_CONSTRUCTION_H

#include "AliModuleConstruction.h"
#include "AliModuleType.h"

#include <globals.hh>

class AliFiles;

class G4VPhysicalVolume;
class G4LogicalVolume;

class AliSingleModuleConstruction : public AliModuleConstruction
{  
  public:
    AliSingleModuleConstruction(G4String moduleName, G4int version,
                                AliModuleType moduleType = kDetector);
    AliSingleModuleConstruction(const AliSingleModuleConstruction& right);
    // --> protected
    // AliSingleModuleConstruction();
    virtual ~AliSingleModuleConstruction();
    
    // operators
    AliSingleModuleConstruction& operator=(
                                const AliSingleModuleConstruction &right);

    // static set/get methods 
    static void SetWorld(G4VPhysicalVolume* world);
    static G4VPhysicalVolume* GetWorld();

    // methods
    void Configure(const AliFiles& files);    
    virtual void Construct();

    // set methods
    void SetProcessConfig(G4bool processConfig);
    void SetModuleFrameName(G4String moduleFrameName);
    void SetModuleType(AliModuleType type);

    // get methods
    G4int  GetVersion() const;
    AliModuleType GetType() const;
    G4bool GetProcessConfig() const;
    
  protected:
    AliSingleModuleConstruction();

  private:    
    // static data members
    static G4VPhysicalVolume* fgWorld;       //top (world) physical volume

    // data members
    G4int           fVersion;        //module version
    AliModuleType   fType;           //module type (detector/structure)
    G4bool          fProcessConfig;  //control for processing Config.C
};             

// inline methods

inline G4VPhysicalVolume* AliSingleModuleConstruction::GetWorld()
{ return fgWorld; }

inline void AliSingleModuleConstruction::SetWorld(G4VPhysicalVolume* world)
{ fgWorld = world; }

inline void AliSingleModuleConstruction::SetProcessConfig(G4bool processConfig)
{ fProcessConfig = processConfig; }

inline void AliSingleModuleConstruction::SetModuleFrameName(G4String name)
{ fModuleFrameName = name; }

inline void AliSingleModuleConstruction::SetModuleType(AliModuleType type)
{ fType = type; }

inline G4int AliSingleModuleConstruction::GetVersion() const
{ return fVersion; }

inline AliModuleType AliSingleModuleConstruction::GetType() const
{ return fType; }

inline G4bool AliSingleModuleConstruction::GetProcessConfig() const
{ return fProcessConfig; }

#endif //ALI_SINGLE_MODULE_CONSTRUCTION_H
