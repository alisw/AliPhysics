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
#include <TPRegexp.h>
#include <TSystem.h>
#include <TError.h>

#include <AliLog.h>

#include <TEveUtil.h>
#include <TEveManager.h>
#include <TEveSelection.h>

#include <Getline.h>

int main(int argc, char **argv)
{
  static const TEveException kEH("alieve::main");

  if (gSystem->Getenv("ALICE_ROOT") == 0)
  {
    Error(kEH.Data(), "ALICE_ROOT is not defined, aborting.");
    gSystem->Exit(1);
  }

  TString evedir(Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));

  if (gSystem->AccessPathName(evedir) == kTRUE)
  {
    Error(kEH.Data(), "Directory $ALICE_ROOT/EVE does not exist.");
    gSystem->Exit(1);
  }

  TString macPath(gROOT->GetMacroPath());
  macPath += Form(":%s/macros", evedir.Data());
  gInterpreter->AddIncludePath(evedir);
  if (gSystem->Getenv("ALICE_ROOT") != 0)
  {
    macPath += Form(":%s/alice-macros", evedir.Data());
    gInterpreter->AddIncludePath(Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/include", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(gSystem->Getenv("ALICE_ROOT"));
  }
  {
    // TabCom fails on double-colon in macro-path.
    TPMERegexp doubleColon(":{2,}", "og");
    doubleColon.Substitute(macPath, ":");
  }
  gROOT->SetMacroPath(macPath);

  // How to hadle AliLog properly?
  AliLog* log = new AliLog;

  TRint app("App", &argc, argv);

  TEveManager::Create();
  gEve->GetSelection()->SetPickToSelect(TEveSelection::kPS_Projectable);
  gEve->GetHighlight()->SetPickToSelect(TEveSelection::kPS_Projectable);

  gEve->RegisterGeometryAlias("Default", Form("%s/alice-data/default_geo.root", evedir.Data()));

  app.Run(kTRUE);

  TEveManager::Terminate();

  delete log;

  return 0;
}
