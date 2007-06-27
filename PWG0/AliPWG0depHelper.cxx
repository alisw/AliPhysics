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

/* $Id$ */

#include <TList.h>
#include <TParticle.h>

#include <AliPWG0depHelper.h>

#include <AliHeader.h>
#include <AliStack.h>
#include <AliLog.h>

#include <AliGenEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliGenCocktailEventHeader.h>

//____________________________________________________________________
ClassImp(AliPWG0depHelper)

//____________________________________________________________________
Int_t AliPWG0depHelper::GetPythiaEventProcessType(AliHeader* aHeader, Bool_t adebug) {
  //
  // get the process type of the event.
  //

  // can only read pythia headers, either directly or from cocktalil header
  AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(aHeader->GenEventHeader());

  if (!pythiaGenHeader) {

    AliGenCocktailEventHeader* genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader*>(aHeader->GenEventHeader());
    if (!genCocktailHeader) {
      printf("AliPWG0depHelper::GetProcessType : Unknown header type (not Pythia or Cocktail). \n");
      return -1;
    }

    TList* headerList = genCocktailHeader->GetHeaders();
    if (!headerList) {
      return -1;
    }

    for (Int_t i=0; i<headerList->GetEntries(); i++) {
      pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(headerList->At(i));
      if (pythiaGenHeader)
        break;
    }

    if (!pythiaGenHeader) {
      printf("AliPWG0depHelper::GetProcessType : Could not find Pythia header. \n");
      return -1;
    }
  }

  if (adebug) {
    printf("AliPWG0depHelper::GetProcessType : Pythia process type found: %d \n",pythiaGenHeader->ProcessType());
  }

  return pythiaGenHeader->ProcessType();
}

//____________________________________________________________________
TParticle* AliPWG0depHelper::FindPrimaryMother(AliStack* stack, Int_t label)
{
  //
  // Finds the first mother among the primary particles of the particle identified by <label>,
  // i.e. the primary that "caused" this particle
  //

  Int_t motherLabel = FindPrimaryMotherLabel(stack, label);
  if (motherLabel < 0)
    return 0;

  return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliPWG0depHelper::FindPrimaryMotherLabel(AliStack* stack, Int_t label)
{
  //
  // Finds the first mother among the primary particles of the particle identified by <label>,
  // i.e. the primary that "caused" this particle
  //
  // returns its label
  //

  Int_t nPrim  = stack->GetNprimary();

  while (label >= nPrim)
  {
    //printf("Particle %d (pdg %d) is not a primary. Let's check its mother %d\n", label, mother->GetPdgCode(), mother->GetMother(0));

    TParticle* particle = stack->Particle(label);
    if (!particle)
    {
      AliDebugGeneral("FindPrimaryMother", AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack.", label));
      return -1;
    }
 
    // find mother
    if (particle->GetMother(0) < 0)
    {
      AliDebugGeneral("FindPrimaryMother", AliLog::kError, Form("UNEXPECTED: Could not find mother of secondary particle %d.", label));
      return -1;
    }

    label = particle->GetMother(0);
  }

  return label;
}


