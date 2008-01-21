// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <TInterpreter.h>
#include <TRint.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TError.h>

#include <TEveUtil.h>
#include <TEveManager.h>

#include <Getline.h>

int main(int argc, char **argv)
{
  static const TEveException eH("alieve::main");

  if (gSystem->Getenv("REVESYS") == 0) {
    if (gSystem->Getenv("ALICE_ROOT") != 0) {
      Info(eH.Data(), "setting REVESYS from ALICE_ROOT.");
      gSystem->Setenv("REVESYS", Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));
    } else {
      Error(eH.Data(), "REVESYS not defined, neither is ALICE_ROOT.");
      gSystem->Exit(1);
    }
  }
  if (gSystem->AccessPathName(gSystem->Getenv("REVESYS")) == kTRUE) {
    Error(eH.Data(), "REVESYS '%s' does not exist.", gSystem->Getenv("REVESYS"));
    gSystem->Exit(1);
  }

  TString macPath(gROOT->GetMacroPath());
  macPath += Form(":%s/macros", gSystem->Getenv("REVESYS"));
  gInterpreter->AddIncludePath(gSystem->Getenv("REVESYS"));
  if (gSystem->Getenv("ALICE_ROOT") != 0) {
    macPath += Form(":%s/alice-macros", gSystem->Getenv("REVESYS"));
    gInterpreter->AddIncludePath(Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/include", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(gSystem->Getenv("ALICE_ROOT"));
  }
  gROOT->SetMacroPath(macPath);

  TRint app("App", &argc, argv);

  TEveManager::Create();

  app.Run(); // Never returns.

  return 0;
}
