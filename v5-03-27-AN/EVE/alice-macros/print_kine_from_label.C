// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TParticle.h>

#include <AliRunLoader.h>
#include <AliStack.h>
#include <AliEveEventManager.h>
#endif

void print_kine_from_label(Int_t label)
{
  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();

  printf("Number primaries %d, all particles %d, label %d\n",
	 stack->GetNprimary(), stack->GetNtrack(), label);
  if (label < 0 || label >= stack->GetNtrack()) {
    printf("  Label exceeds available range.\n");
    return;
  }

  TParticle* part = stack->Particle(label);
  if(part != 0) {
    part->Print();
    while(part->GetMother(0) >= 0) {
      part = stack->Particle(part->GetMother(0));
      part->Print();
    }
  }
}
