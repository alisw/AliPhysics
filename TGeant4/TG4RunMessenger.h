// $Id$
// Category: run
//
// Messenger class that defines commands for TG4RunManager

#ifndef TG4_RUN_MESSENGER_H
#define TG4_RUN_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>

class TG4RunManager;
class TG4UICmdWithAComplexString;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;

class TG4RunMessenger: public G4UImessenger
{
  public:
    TG4RunMessenger(TG4RunManager* runManager);
    // protected
    // TG4RunMessenger();
    // TG4RunMessenger(const TG4RunMessenger& right);
    virtual ~TG4RunMessenger();
   
    // methods 
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    TG4RunMessenger();  
    TG4RunMessenger(const TG4RunMessenger& right);

    // operators
    TG4RunMessenger& operator=(const TG4RunMessenger& right);

  private:
    // data members
    TG4RunManager*  fRunManager; //associated class   
    G4UIdirectory*  fDirectory;  //command directory

    G4UIcmdWithoutParameter*    fRootCmd;        //command: root
    G4UIcmdWithAString*         fRootMacroCmd;   //command: rootMacro 
    TG4UICmdWithAComplexString* fRootCommandCmd; //command: rootCmd 
    G4UIcmdWithoutParameter*    fG3DefaultsCmd;  //command: g3Defaults   
};

#endif //TG4_RUN_MESSENGER_H
