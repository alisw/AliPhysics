// $Id$
// Category: global
//
// See the class description in the header file.

#include "AliFiles.h"

#include <stdlib.h>

// static data members

const G4String AliFiles::fgTop =    getenv("AG4_INSTALL");    
const G4String AliFiles::fgConfig = fgTop + "/../macros/Config.C";    
const G4String AliFiles::fgDetConfig1 = fgTop + "/macro/";    
const G4String AliFiles::fgDetConfig2 = "/Config";    
const G4String AliFiles::fgDetConfig3 = ".C";    
const G4String AliFiles::fgDetConfig4 = ".in";    
const G4String AliFiles::fgDetConfigName1 = "Config(";    
const G4String AliFiles::fgDetConfigName2 = ")";    
const G4String AliFiles::fgDetData1 = fgTop + "/macro/";
const G4String AliFiles::fgDetData2 = "/g3calls_";    
const G4String AliFiles::fgDetData3 = ".dat";    
const G4String AliFiles::fgSTRUCT = "STRUCT/";

AliFiles::AliFiles() {
//
}
  
AliFiles::~AliFiles() {
//
}
  
