// $Id$
// Category: geometry
//
// Author: I. Hrivnacova
//
// Class AliModulesCompositionMessenger
// ------------------------------------
// Messenger class that defines commands for AliModulesComposition.

#ifndef ALI_MODULES_COMPOSITION_MESSENGER_H
#define ALI_MODULES_COMPOSITION_MESSENGER_H

#include <G4UImessenger.hh>
#include <globals.hh>

class AliModulesComposition;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class AliModulesCompositionMessenger: public G4UImessenger
{
  public:
    AliModulesCompositionMessenger(AliModulesComposition* modulesComposition);
    // --> protected
    // AliModulesCompositionMessenger();
    // AliModulesCompositionMessenger(
    //       const AliModulesCompositionMessenger& right);
    virtual ~AliModulesCompositionMessenger();

    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String newValues);
    
  protected:
    AliModulesCompositionMessenger();
    AliModulesCompositionMessenger(
             const AliModulesCompositionMessenger& right);

    // operators
    AliModulesCompositionMessenger& operator=(
             const AliModulesCompositionMessenger &right);
             
  private:
    AliModulesComposition*  fModulesComposition; //associated class
    G4UIdirectory*          fDirectory;          //command directory
    
    // commands data members
    G4UIcmdWithAString*         fFieldTypeCmd;        //command: fieldType
    G4UIcmdWithADoubleAndUnit*  fUniformFieldValueCmd;//command: uniformFieldValue
    G4UIcmdWithABool*           fSetReadGeometryCmd;  //command: readGeometry   
    G4UIcmdWithABool*           fSetWriteGeometryCmd; //command: writeGeometry    
    G4UIcmdWithoutParameter*    fPrintMaterialsCmd;   //command: printMatrials     
    G4UIcmdWithoutParameter*    fGenerateXMLCmd;      //command: generateXML
};

#endif //ALI_MODULES_COMPOSITION_MESSENGER_H

