// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4VerboseMessenger
// ------------------
// Messenger class that defines commands for 
// the verbose classes.

#ifndef TG4_VERBOSE_MESSENGER_H
#define TG4_VERBOSE_MESSENGER_H 

#include <G4UImessenger.hh>
#include <globals.hh>
#include <g4std/vector>

class TG4VVerbose;

class G4UIdirectory;
class G4UIcmdWithAnInteger;

class TG4VerboseMessenger: public G4UImessenger
{
  typedef G4std::vector<TG4VVerbose*>           VerboseVector;
  typedef G4std::vector<G4UIcmdWithAnInteger*>  CommandVector;

  public:
    TG4VerboseMessenger(const G4String& directoryName);
    // --> protected   
    // TG4VerboseMessenger();
    // TG4VerboseMessenger(const TG4VerboseMessenger& right);
    virtual ~TG4VerboseMessenger();
   
    // methods 
    virtual void AddCommand(TG4VVerbose* verbose, const G4String& cmdName);
    virtual void SetNewValue(G4UIcommand* command, G4String string);
    
  protected:
    TG4VerboseMessenger();  
    TG4VerboseMessenger(const TG4VerboseMessenger& right);

    // operators
    TG4VerboseMessenger& operator=(const TG4VerboseMessenger& right);

  private:
    // methods
    void SetNewValueToAll(const G4String value) const;
  
    // data members
    const G4String        fkDirectoryName;  //command directory name
    G4UIdirectory*        fDirectory;       //command directory
    G4UIcmdWithAnInteger* fGlobalVerboseCmd;//global verbose command 
    VerboseVector         fVerboseVector;   //associated verbose instances
    CommandVector         fCommandVector;   //verbose commands
};

#endif //TG4_VERBOSE_MESSENGER_H
