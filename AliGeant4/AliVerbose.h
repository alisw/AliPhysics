// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class AliVerbose
// ----------------
// Class defines the verbose level
// and the static messenger (common for all instances).
// Used as a base class for all AliGeant4 verbose classes;
// enables to handle the standard output in a common way.

#ifndef ALI_VERBOSE_H
#define ALI_VERBOSE_H

#include "TG4VVerbose.h"
#include "TG4VerboseMessenger.h"

#include <globals.hh>

class AliVerbose : public TG4VVerbose
{
  public:
    // AliVerbose(); --> private      
    AliVerbose(const G4String& cmdName);
    AliVerbose(const G4String& cmdName, G4int verboseLevel);      
    virtual ~AliVerbose();

  private:
    AliVerbose();

    // methods
    virtual TG4VerboseMessenger* CreateMessenger();    

    // static data members
    static const G4String        fgkDirectoryName;// directory name
    static TG4VerboseMessenger*  fgMessenger;     // messenger
};     

#endif //TG4_VERBOSE_H
