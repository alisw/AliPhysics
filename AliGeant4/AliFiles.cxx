// $Id$
// Category: global
//
// See the class description in the header file.

#include "AliFiles.h"

#include <stdlib.h>

// static data members

const G4String AliFiles::fgkTop =    getenv("AG4_INSTALL");    
const G4String AliFiles::fgkConfig = fgkTop + "/../macros/Config.C";    
const G4String AliFiles::fgkDetConfig1 = fgkTop + "/macro/";    
const G4String AliFiles::fgkDetConfig2 = "/Config";    
const G4String AliFiles::fgkDetConfig3 = ".C";    
const G4String AliFiles::fgkDetConfig4 = ".in";    
const G4String AliFiles::fgkDetConfigName1 = "Config(";    
const G4String AliFiles::fgkDetConfigName2 = ")";    
const G4String AliFiles::fgkDetData1 = fgkTop + "/macro/";
const G4String AliFiles::fgkDetData2 = "/g3calls_";    
const G4String AliFiles::fgkDetData3 = ".dat";    
const G4String AliFiles::fgkSTRUCT = "STRUCT/";

AliFiles::AliFiles() {
//
}
  
AliFiles::~AliFiles() {
//
}
  
