// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliModulesComposition
// ---------------------------
// Detector construction base class for building geometry
// composed from modules.

#ifndef ALI_MODULES_COMPOSITION_H
#define ALI_MODULES_COMPOSITION_H

#include "AliModulesCompositionMessenger.h"
#include "AliModuleType.h"
#include "AliVerbose.h"

#include <TG4MagneticFieldType.h>

#include <G4VUserDetectorConstruction.hh>
#include <globals.hh>
#include <g4std/vector>

class AliModuleConstruction;
class G4MagneticField;

class G4VPhysicalVolume;

class AliModulesComposition : public G4VUserDetectorConstruction,
                              public AliVerbose
{
  typedef G4std::vector<AliModuleConstruction*> AliModuleConstructionVector;
  
  public:
    AliModulesComposition();
    // --> protected
    // AliModulesComposition(const AliModulesComposition& right);
    virtual ~AliModulesComposition();

    // methods
    virtual G4VPhysicalVolume* Construct() = 0;
    virtual void GenerateXMLGeometry() const = 0;
    virtual void PrintMaterials() const;

    // set methods
    void SetFieldType(TG4MagneticFieldType fieldType);
    void SetUniformFieldValue(G4double fieldValue);
    void SetReadGeometry(G4bool readGeometry);
    void SetWriteGeometry(G4bool writeGeometry);
    
  protected:
    AliModulesComposition(const AliModulesComposition& right);

    // operators
    AliModulesComposition& operator=(const AliModulesComposition& right);

    // methods  
    void AddModule(const G4String& name, 
                   G4int version,
                   AliModuleType moduleType = kDetector);
    void ConstructModules();
    void SetProcessConfigToModules(G4bool processConfig);
    
  private:    
    // methods
    void CreateMagneticField();
    void Configure();
    void CreateG4Geometry();
    void SetReadGeometryToModules(G4bool readGeometry);
    void SetWriteGeometryToModules(G4bool writeGeometry);

    // data members
    AliModulesCompositionMessenger  fMessenger;      //messenger
    AliModuleConstructionVector     fModuleConstructionVector; //..
                                        //vector of AliModuleConstruction
    TG4MagneticFieldType            fMagneticFieldType;//magnetic field type
    G4MagneticField*   fMagneticField;  //magnetic field
    G4bool             fReadGeometry;   //option applied to all modules
    G4bool             fWriteGeometry;  //option applied to all modules     
};

// inline methods

inline void AliModulesComposition::SetReadGeometry(G4bool readGeometry)
{ fReadGeometry = readGeometry; }

inline void AliModulesComposition::SetWriteGeometry(G4bool writeGeometry)
{ fWriteGeometry = writeGeometry; }

#endif //ALI_MODULES_COMPOSITION_H

