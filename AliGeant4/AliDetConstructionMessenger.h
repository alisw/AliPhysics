// $Id$
// Category: geometry
//
// Messenger class that defines commands for AliDetConstruction.

#ifndef ALI_DET_CONSTRUCTION_MESSENGER_H
#define ALI_DET_CONSTRUCTION_MESSENGER_H

#include <G4UImessenger.hh>
#include <globals.hh>

class AliDetConstruction;

class G4UIcommand;
class G4UIcmdWithABool;

class AliDetConstructionMessenger: public G4UImessenger
{
  public:
    AliDetConstructionMessenger(AliDetConstruction* detConstruction);
    // --> protected
    // AliDetConstructionMessenger();  
    // AliDetConstructionMessenger(const AliDetConstructionMessenger& right);
    virtual ~AliDetConstructionMessenger();

    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String newValues);

  protected:
    AliDetConstructionMessenger();
    AliDetConstructionMessenger(const AliDetConstructionMessenger& right);

    // operators
    AliDetConstructionMessenger& operator=(
                                const AliDetConstructionMessenger& right);
  
  private:
    // data members
    AliDetConstruction*  fDetConstruction;     //associated class
    G4UIcmdWithABool*    fSetAllSensitiveCmd;  //command: setAllSensitive   
    G4UIcmdWithABool*    fSetReadGeometryCmd;  //command: setReadGeometry   
    G4UIcmdWithABool*    fSetWriteGeometryCmd; //command: setWriteGeometry   
};

#endif //ALI_DET_CONSTRUCTION_MESSENGER_H

