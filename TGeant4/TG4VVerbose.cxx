// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Class TG4Verbose
// -------------------
// See the class description in the header file.

#include "TG4VVerbose.h"
#include "TG4Globals.h"

#include <math.h>

// static data members
const G4int TG4VVerbose::fgkDefaultVerboseLevel = 1;

//_____________________________________________________________________________
TG4VVerbose::TG4VVerbose()
  : fVerboseLevel(fgkDefaultVerboseLevel) {
//
}
  
//_____________________________________________________________________________
TG4VVerbose::TG4VVerbose(G4int verboseLevel) 
  : fVerboseLevel(verboseLevel) {
// 
}

//_____________________________________________________________________________
TG4VVerbose::~TG4VVerbose() {
//
}

