// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4Verbose
// ----------------
// Class defines the verbose level
// and the static messenger (common for all instances).
// Used as a base class for all TGeant4 verbose classes;
// enables to handle the standard output in a common way.

#ifndef TG4_VERBOSE_H
#define TG4_VERBOSE_H

#include "TG4VVerbose.h"
#include "TG4VerboseMessenger.h"

#include <globals.hh>

class TG4Verbose : public TG4VVerbose
{
  public:
    // TG4Verbose(); --> private      
    TG4Verbose(const G4String& cmdName);
    TG4Verbose(const G4String& cmdName, G4int verboseLevel);      
    virtual ~TG4Verbose();

  private:
    TG4Verbose();

    // methods
    virtual TG4VerboseMessenger* CreateMessenger();    

    // static data members
    static const G4String        fgkDirectoryName;// directory name
    static TG4VerboseMessenger*  fgMessenger;     // messenger
};     

#endif //TG4_VERBOSE_H
