// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2008

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Import tracks from kinematics-tree / particle-stack.
// Preliminary/minimal solution.
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TParticle.h>
#include <TParticlePDG.h>

#include <AliStack.h>
#include <AliRunLoader.h>
#include <AliEveEventManager.h>
#endif
void
kine_print(Double_t min_pt = 0, Double_t min_p = 0)
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();
  if (!stack) {
    Error("kine_tracks.C", "can not get kinematics.");
    return;
  }

  printf("\n");
  printf("%4s | %-11s | %3s | %4s | %9s | %s %s\n",
         "id", "name", "sts",
         "mth", "dghtrs", "p", "P");
  printf("------------------------------------------------------------\n");
  Int_t N = stack->GetNtrack();
  for (Int_t i=0; i<N; ++i)
  {
    TParticle* p = stack->Particle(i);
    printf("%4d | %-11s | %3d | %4d | %4d %4d | %d %d\n",
           i, p->GetName(), p->GetStatusCode(),
           p->GetMother(0), p->GetDaughter(0), p->GetDaughter(1),
           p->IsPrimary(), stack->IsPhysicalPrimary(i));
  }
}
