// $Id$
// Flugg tag $Name$

#include "AliDetConstructionHelper.h"
#include "AliDetConstruction.h"
#include "AliFiles.h"
#include "AliBODY.h"
#include "AliHALL.h"

/*
#include <TROOT.h>
#include <TCint.h>
*/

AliDetConstructionHelper::AliDetConstructionHelper()
{
/*
  // initialize the Alice setup
  gROOT->LoadMacro("/home/ivana/AliRoot/macros/Config.C");
  gInterpreter->ProcessLine("Config()");
*/
  G4cout << "Creating world" << G4endl;
  AliBODY *BODY = new AliBODY("BODY","Alice envelop");
  G4cout << "Creating hall" << G4endl;
  AliHALL *HALL = new AliHALL("HALL","Alice Hall");

  fFiles = new AliFiles();
  fDetConstruction = new AliDetConstruction();
}

AliDetConstructionHelper::~AliDetConstructionHelper()
{
  delete fDetConstruction;
  delete fFiles;
}

