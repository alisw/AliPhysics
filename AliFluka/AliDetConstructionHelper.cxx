/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
//
// Author: I. Hrivnacova
//
// Class AliDetConstructionHelper
// ------------------------------
// See the class description in the header file.

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

