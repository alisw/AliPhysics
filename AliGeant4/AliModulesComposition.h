// $Id$
// Category: geometry
//
// Detector construction base class for building geometry
// composed from independent modules with availability of interactive modules
// setup.

#ifndef ALI_MODULES_COMPOSITION_H
#define ALI_MODULES_COMPOSITION_H

#include "AliSingleModuleConstruction.h"
#include "AliDetSwitch.h"
#include "AliModuleType.h"

#include <G4VUserDetectorConstruction.hh>
#include <globals.hh>

#include <g4rw/tpordvec.h>

class AliModulesCompositionMessenger;
class AliMoreModulesConstruction;
class AliMagneticField;

class G4VPhysicalVolume;

class AliModulesComposition : public G4VUserDetectorConstruction
{
  typedef G4RWTPtrOrderedVector<AliDetSwitch>  AliDetSwitchRWVector;
  typedef G4RWTPtrOrderedVector<AliSingleModuleConstruction>
                                AliSingleModuleConstructionRWVector; 

  public:
    AliModulesComposition();
    // --> protected
    // AliModulesComposition(const AliModulesComposition& right);
    virtual ~AliModulesComposition();

    // methods
    virtual G4VPhysicalVolume* Construct() = 0;
    void SwitchDetOn(G4String moduleNameVer);
    void SwitchDetOn(G4String moduleName, G4int version);
    void SwitchDetOnDefault(G4String moduleName);
    void SwitchDetOff(G4String moduleName);
    void PrintSwitchedDets() const;
    void PrintAvailableDets() const;

    // set methods
    void SetMagField(G4double fieldValue);
    void SetAllLVSensitive(G4bool allLVSensitive);
    void SetReadGeometry(G4bool readGeometry);
    void SetWriteGeometry(G4bool writeGeometry);
    void SetProcessConfigToModules(G4bool processConfig);
    
    // get methods
    const G4RWTPtrOrderedVector<AliDetSwitch>& GetDetSwitchVector() const;
    G4String GetSwitchedDetsList() const;
    G4String GetAvailableDetsList() const;
    G4String GetAvailableDetsListWithCommas() const;
    G4String GetDetNamesList() const;
    G4String GetDetNamesListWithCommas() const;
    //G4ThreeVector GetMagField() const;
    
  protected:
    AliModulesComposition(const AliModulesComposition& right);

    // operators
    AliModulesComposition& operator=(const AliModulesComposition& right);

    // methods  
    void AddDetSwitch(AliDetSwitch* detSwitch);
    void AddSingleModuleConstruction(G4String moduleName, G4int version,
                                     AliModuleType moduleType = kDetector);
    void AddMoreModuleConstruction(G4String moduleName, G4int version,
                                     AliModuleType moduleType = kDetector);
    void ConstructModules();

  private:    
    // methods
    void SetReadGeometryToModules(G4bool readGeometry);
    void SetWriteGeometryToModules(G4bool writeGeometry);
    void SetAllLVSensitiveToModules(G4bool allSensitive);

    // data members
    AliDetSwitchRWVector                fDetSwitchVector;          //..         
                                          //vector of AliDetSwitch
    AliSingleModuleConstructionRWVector fModuleConstructionVector; //..
				          //vector of 
					  //AliSingleModuleConstruction 
    AliMoreModulesConstruction*         fMoreModulesConstruction;  //..
                                          //AliMoreModulesConstruction

    AliMagneticField*                fMagneticField;  //magnetic field
    AliModulesCompositionMessenger*  fMessenger;      //messenger
    G4bool  fAllLVSensitive; //option applied to all modules
    G4bool  fReadGeometry;   //option applied to all modules
    G4bool  fWriteGeometry;  //option applied to all modules  
};

// inline methods

inline void AliModulesComposition::SetAllLVSensitive(G4bool allLVSensitive)
{ fAllLVSensitive = allLVSensitive; }

inline void AliModulesComposition::SetReadGeometry(G4bool readGeometry)
{ fReadGeometry = readGeometry; }

inline void AliModulesComposition::SetWriteGeometry(G4bool writeGeometry)
{ fWriteGeometry = writeGeometry; }

#endif //ALI_MODULES_COMPOSITION_H

