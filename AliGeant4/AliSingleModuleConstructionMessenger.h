// $Id$
// Category: geometry
//
// Messenger class that defines command directory for the
// AliSingleModuleConstruction instances.

#ifndef ALI_SINGLE_MODULE_CONSTRUCTION_MESSENGER_H
#define ALI_SINGLE_MODULE_CONSTRUCTION_MESSENGER_H

#include <G4UImessenger.hh>
#include <globals.hh>

class AliSingleModuleConstruction;

class G4UIcmdWithABool;

class AliSingleModuleConstructionMessenger: public G4UImessenger
{
  public:
    AliSingleModuleConstructionMessenger(
       AliSingleModuleConstruction* moduleConstruction, G4String moduleName);
    // --> protected
    // AliSingleModuleConstructionMessenger();
    // AliSingleModuleConstructionMessenger(
    //                 const AliSingleModuleConstructionMessenger& right); 
    virtual ~AliSingleModuleConstructionMessenger();

    // methods
    virtual void SetNewValue(G4UIcommand* command, G4String newValues);

  protected:
    AliSingleModuleConstructionMessenger();
    AliSingleModuleConstructionMessenger(
                       const AliSingleModuleConstructionMessenger& right);

    // operators
    AliSingleModuleConstructionMessenger& operator=(
                       const AliSingleModuleConstructionMessenger& right);
  
  private:
    // data members
    AliSingleModuleConstruction*  fModuleConstruction; //associated class

    G4UIcmdWithABool*  fSetAllSensitiveCmd; //command: setAllSensitive   
};

#endif //ALI_MODULE_CONSTRUCTION_MESSENGER_H

