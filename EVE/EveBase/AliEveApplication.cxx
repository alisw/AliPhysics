// Author: Mihai Niculescu 2012

/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, all rights reserved. 					 *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          										 *
 * full copyright notice.                                                																						 *
 **************************************************************************/

#include <TInterpreter.h>
#include <TSystem.h>
#include <TString.h>
#include <TROOT.h>

#include <AliLog.h>

#include <AliEveApplication.h>
#include <AliEveManager.h>


ClassImp(AliEveApplication)

AliEveApplication::AliEveApplication(const char* appClassName, int* argc, char** argv, void* options, int numOptions, Bool_t noLogo)
	: TRint(appClassName, argc, argv, options, numOptions, noLogo)
{
  Init();
}

AliEveApplication::~AliEveApplication()
{}

void  AliEveApplication::Init()
{

	static const TEveException kEH("alieve::main");

  TString evedir(Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));

  TString macPath(gROOT->GetMacroPath());
  macPath += Form(":%s/macros", evedir.Data());
  gInterpreter->AddIncludePath(evedir);

  macPath += Form(":%s/alice-macros", evedir.Data());
  gInterpreter->AddIncludePath(Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));
  gInterpreter->AddIncludePath(Form("%s/PWG0", gSystem->Getenv("ALICE_ROOT")));
  gInterpreter->AddIncludePath(Form("%s/include", gSystem->Getenv("ALICE_ROOT")));
  gInterpreter->AddIncludePath(gSystem->Getenv("ALICE_ROOT"));
 
  gROOT->SetMacroPath(macPath);

  // make sure logger is instantiated
  AliLog::GetRootLogger();


}
