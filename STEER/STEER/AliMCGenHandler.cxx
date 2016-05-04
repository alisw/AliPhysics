/************************************************************************** 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved  *
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

/* $Id: AliMCGenHandler.cxx 64675 2013-10-23 12:21:37Z hristov $ */

//-------------------------------------------------------------------------
//                          Class AliMCGenHandler
// This class can be used with the analysis framework to generate event on
// the fly and analyse them.
//      
// Origin: Andrei Gheata, Jan Fiete Grosse-Oetringhaus
//-------------------------------------------------------------------------

#include "AliMCGenHandler.h"
#include "AliMCEvent.h"
#include "AliPDG.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliGenerator.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliAnalysisManager.h"

#include "TMacro.h"
#include "TInterpreter.h"

ClassImp(AliMCGenHandler)

AliMCGenHandler::AliMCGenHandler() :
    AliInputEventHandler(),
    fMCEvent(0),
    fEventNumber(0),
    fStack(0),
    fHeader(0),
    fGenerator(0),
    fSeedMode(0),
    fSeed(0),
    fGeneratorMacroPath(),
    fGeneratorMacroParameters(),
    fGeneratorCustomization(0)
{
  //
  // Default constructor
  //
  // Be sure to add all particles to the PDG database
  AliPDG::AddParticlesToPdgDataBase();
}

AliMCGenHandler::AliMCGenHandler(const char* name, const char* title) :
    AliInputEventHandler(name, title),
    fMCEvent(0),
    fEventNumber(0),
    fStack(0),
    fHeader(0),
    fGenerator(0),
    fSeedMode(0),
    fSeed(0),
    fGeneratorMacroPath(),
    fGeneratorMacroParameters(),
    fGeneratorCustomization(0)
{
  //
  // Constructor
  //
  // Be sure to add all particles to the PDG database
  AliPDG::AddParticlesToPdgDataBase();
}

AliMCGenHandler::~AliMCGenHandler()
{ 
    // Destructor
    delete fMCEvent;
    delete fGenerator;
}

Bool_t AliMCGenHandler::Init(Option_t* /*opt*/)
{ 
  // Initialize input
    //
    
    if (!fGenerator)
      CreateGenerator();

    
    if (fSeedMode == 0)
      Printf("AliMCGenHandler::Init: Not setting any seed. Seed needs to be set externally!");
    else
    {
      if (fSeedMode == 1)
      {
	Printf("AliMCGenHandler::Init: Using manually set seed");
      }
      else if (fSeedMode == 2)
      {
	Printf("AliMCGenHandler::Init: Taking seed from current time");
	fSeed = time(0);
      }
      else if (fSeedMode == 3)
      {
	Printf("AliMCGenHandler::Init: Taking seed from AliEn job id");
	TString tmp(gSystem->Getenv("ALIEN_PROC_ID"));
	fSeed = tmp.Atoi();
	if (tmp.Length() == 0 || fSeed == 0)
	  AliFatal(Form("Could not retrieve AliEn job id for seed. The variable ALIEN_PROC_ID contains %s", tmp.Data()));
      }
      else
	AliFatal(Form("Seed mode %d unknown", fSeedMode));

      Printf("AliMCGenHandler::Init: Using seed: %d", fSeed);
      gRandom->SetSeed(fSeed);
      fGenerator->SetSeed(fSeed);
    }

    AliRunLoader* rl = AliRunLoader::Open("galice.root","FASTRUN","recreate");
    rl->MakeTree("E");
    gAlice->SetRunLoader(rl);
    rl->MakeStack();
    fStack = rl->Stack();
    fHeader = rl->GetHeader();
    
    fGenerator->SetStack(fStack);
    fGenerator->Init();
    
    fMCEvent = new AliMCEvent;

    return kTRUE;
}

void AliMCGenHandler::CreateGenerator()
{
      if (fGeneratorMacroPath.Length() == 0)
	AliFatal("fGeneratorMacroPath empty!");
      
      TString macroPath(fGeneratorMacroPath);
      
      if (!fGeneratorMacroPath.BeginsWith("$"))
        macroPath.Form("$ALICE_ROOT/%s", fGeneratorMacroPath.Data());
      
      macroPath = gSystem->ExpandPathName(macroPath.Data());

      if (gSystem->AccessPathName(macroPath))
	  AliFatal(Form("Cannot find macro %s", macroPath.Data()));

      TMacro m(macroPath);
      Int_t error = 0;
      Long64_t retval = m.Exec(fGeneratorMacroParameters, &error);
      if (error != TInterpreter::kNoError)
	  AliFatal(Form("Macro interpretation %s failed", macroPath.Data()));

      if (retval<0)
	AliFatal(Form("The macro %s did not return a valid generator (1)", macroPath.Data()));

      fGenerator = reinterpret_cast<AliGenerator*>(retval);
      if (!fGenerator)
	AliFatal(Form("The macro %s did not return a valid generator (2)", macroPath.Data()));

      // customization from LEGO train
      if (fGeneratorCustomization) {
	fGeneratorCustomization->Exec(Form("(AliGenerator*) %p", fGenerator), &error);
	if (error != TInterpreter::kNoError)
	  AliFatal("Execution of generator customization failed");
      }
}

Bool_t AliMCGenHandler::BeginEvent(Long64_t /*entry*/)
{ 
    // Begin event

    fHeader->Reset(0, fEventNumber++);
    
    if (AliAnalysisManager::GetAnalysisManager()->GetDebugLevel() > 1)
      Printf("AliMCGenHandler::BeginEvent: Generating Event number %lld", fEventNumber);
      
    fGenerator->Generate();

    fHeader->SetStack(fStack);
    fHeader->SetNprimary(fStack->GetNprimary());
    fHeader->SetNtrack(fStack->GetNtrack());  

    fMCEvent->ConnectHeaderAndStack(fHeader);

    return kTRUE;
}

Bool_t AliMCGenHandler::FinishEvent()
{
    // Clean-up after each event
   
    fMCEvent->Stack()->Reset();

    return kTRUE;
}
