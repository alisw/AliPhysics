// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliModulesComposition
// ---------------------------
// Detector construction base class for building geometry
// composed from independent modules with availability of interactive modules
// setup.

#ifndef ALI_MODULES_COMPOSITION_H
#define ALI_MODULES_COMPOSITION_H

#include "AliModulesCompositionMessenger.h"
#include "AliModuleType.h"

#include <G4VUserDetectorConstruction.hh>
#include <globals.hh>
#include <g4std/vector>

class AliSingleModuleConstruction;
class AliDetSwitch;
class AliMoreModulesConstruction;
class AliMagneticField;

class G4VPhysicalVolume;

class AliModulesComposition : public G4VUserDetectorConstruction
{
  typedef G4std::vector<AliDetSwitch*>    DetSwitchVector;
  typedef DetSwitchVector::iterator       DetSwitchIterator;
  typedef DetSwitchVector::const_iterator DetSwitchConstIterator;

  typedef G4std::vector<AliSingleModuleConstruction*> SingleModuleVector;
  typedef SingleModuleVector::iterator                SingleModuleIterator;    

  public:
    AliModulesComposition();
    // --> protected
    // AliModulesComposition(const AliModulesComposition& right);
    virtual ~AliModulesComposition();

    // methods
    virtual G4VPhysicalVolume* Construct() = 0;
    void SwitchDetOn(const G4String& moduleNameVer);
    void SwitchDetOn(const G4String& moduleName, G4int version);
    void SwitchDetOnDefault(const G4String& moduleName);
    void SwitchDetOff(const G4String& moduleName);
    void PrintSwitchedDets() const;
    void PrintAvailableDets() const;
    void PrintMaterials() const;
    void GenerateXMLGeometry() const;

    // set methods
    void SetMagField(G4double fieldValue);
    void SetReadGeometry(G4bool readGeometry);
    void SetWriteGeometry(G4bool writeGeometry);
    void SetProcessConfigToModules(G4bool processConfig);
    
    // get methods
    G4String GetSwitchedDetsList() const;
    G4String GetAvailableDetsList() const;
    G4String GetAvailableDetsListWithCommas() const;
    G4String GetDetNamesList() const;
    G4String GetDetNamesListWithCommas() const;
    
  protected:
    AliModulesComposition(const AliModulesComposition& right);

    // operators
    AliModulesComposition& operator=(const AliModulesComposition& right);

    // methods  
    void AddDetSwitch(AliDetSwitch* detSwitch);
    void AddSingleModuleConstruction(const G4String& name, G4int version,
                                     AliModuleType moduleType = kDetector);
    void AddMoreModuleConstruction(const G4String& name, G4int version,
                                     AliModuleType moduleType = kDetector);
    void ConstructModules();

    // get methods
    AliDetSwitch* GetDetSwitch(const G4String& moduleName) const;

    // data members
    DetSwitchVector  fDetSwitchVector; //vector of AliDetSwitch
    
  private:    
    // methods
    void SetReadGeometryToModules(G4bool readGeometry);
    void SetWriteGeometryToModules(G4bool writeGeometry);

    // data members
    SingleModuleVector           fModuleConstructionVector; //vector of 
					  //single module constructions 
    AliMoreModulesConstruction*  fMoreModulesConstruction;  //..
                                          //dependent modules construction

    AliMagneticField*               fMagneticField;  //magnetic field
    AliModulesCompositionMessenger  fMessenger;      //messenger
    G4bool  fReadGeometry;        //option applied to all modules
    G4bool  fWriteGeometry;       //option applied to all modules     
};

// inline methods

inline void AliModulesComposition::SetReadGeometry(G4bool readGeometry)
{ fReadGeometry = readGeometry; }

inline void AliModulesComposition::SetWriteGeometry(G4bool writeGeometry)
{ fWriteGeometry = writeGeometry; }

#endif //ALI_MODULES_COMPOSITION_H

