// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class AliVerbose
// -------------------
// See the class description in the header file.

#include "AliVerbose.h"
#include "TG4Globals.h"

#include <math.h>

// static data members
const G4String       AliVerbose::fgkDirectoryName = "/aliVerbose/";
TG4VerboseMessenger* AliVerbose::fgMessenger = 0;

//_____________________________________________________________________________
AliVerbose::AliVerbose(const G4String& cmdName)
  : TG4VVerbose() {
//
  CreateMessenger();
  
  fgMessenger->AddCommand(this, cmdName);  
}
  
//_____________________________________________________________________________
AliVerbose::AliVerbose(const G4String& cmdName, G4int verboseLevel) 
  : TG4VVerbose(verboseLevel) {
// 
  CreateMessenger();

  fgMessenger->AddCommand(this, cmdName);  
}

//_____________________________________________________________________________
AliVerbose::AliVerbose()
  : TG4VVerbose() {
//
}
  
//_____________________________________________________________________________
AliVerbose::~AliVerbose() {
//
  if (fgMessenger) {
    delete fgMessenger;
    fgMessenger = 0;
  }  
}

//
// private methods
//

TG4VerboseMessenger* AliVerbose::CreateMessenger() 
{
// Creates static messenger if it does not yet exists.
// ---

  if (!fgMessenger)
    fgMessenger = new TG4VerboseMessenger(fgkDirectoryName);
    
  return fgMessenger;  
}
    
