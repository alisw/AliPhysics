// $Id$
// Category: geometry
//
// Messenger class that defines commands for AliModulesComposition.

#ifndef ALI_MODULES_COMPOSITION_MESSENGER_H
#define ALI_MODULES_COMPOSITION_MESSENGER_H

#include <G4UImessenger.hh>
#include <globals.hh>

class AliModulesComposition;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
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
    void SetCandidates();
    
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
    G4UIcmdWithAString*         fSwitchOnCmd;         //command: switchOn 
    G4UIcmdWithAString*         fSwitchOffCmd;        //command: switchOn
    G4UIcmdWithoutParameter*    fListCmd;             //command: list
    G4UIcmdWithoutParameter*    fListAvailableCmd;    //command: listAvailable
    G4UIcmdWithADoubleAndUnit*  fFieldValueCmd;       //command: fieldValue
    G4UIcmdWithABool*           fSetAllSensitiveCmd;  //command: setAllSensitive   
    G4UIcmdWithABool*           fSetReadGeometryCmd;  //command: readGeometry   
    G4UIcmdWithABool*           fSetWriteGeometryCmd; //command: writeGeometry    
};

#endif //ALI_MODULES_COMPOSITION_MESSENGER_H

