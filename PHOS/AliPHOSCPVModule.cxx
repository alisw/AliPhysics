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

/*
  $Log$
*/

////////////////////////////////////////////////
//  Manager class for one CPV module          //
//                                            //
//  Author: Yuri Kharlov, IHEP, Protvino      //
//  e-mail: Yuri.Kharlov@cern.ch              //
//  Last modified: 2 November 2000            //
////////////////////////////////////////////////
 
// --- ROOT system ---
#include <TTree.h>

// --- Standard library ---
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// --- galice header files ---
#include "AliRun.h"
#include "AliPHOSCPVModule.h"
#include "AliPHOSCPVHit.h"

//==============================================================================
//                              AliPHOSCPVModule
//==============================================================================

ClassImp(AliPHOSCPVModule)

//______________________________________________________________________________

AliPHOSCPVModule::AliPHOSCPVModule(void) {
  //
  // Allocate an array of hits
  //

  if ( NULL==(fHits=new TClonesArray("AliPHOSCPVHit",100)) ) {
    Error("CPV","Can not create array of hits");
    exit(1);
  }
}

//______________________________________________________________________________

AliPHOSCPVModule::~AliPHOSCPVModule(void)
{
  Clear();
}

//______________________________________________________________________________

void AliPHOSCPVModule::Clear(Option_t *opt="")
{
// Clear hit information

  fHits  -> Clear(opt);
}

//______________________________________________________________________________

void AliPHOSCPVModule::AddHit(TLorentzVector p, Float_t *xy, Int_t ipart)
{
  //
  // Add this hit to the hit list in CPV detector.
  //

  TClonesArray &lhits = *(TClonesArray *)fHits;
  new(lhits[fHits->GetEntriesFast()]) AliPHOSCPVHit(p,xy,ipart);
}

//______________________________________________________________________________

void AliPHOSCPVModule::Print(Option_t *opt)
{
  //
  // Print AliPHOSCPVModule information.
  // options:  'p' - print hits in the module
  //

  Int_t nhits,hit;
  if (strcmp(opt,"p")==0) {
    printf ("CPV module has %d hits\n",nhits=fHits->GetEntriesFast());
    for (hit=0;hit<nhits;hit++) {
      AliPHOSCPVHit *cpvHit = (AliPHOSCPVHit*)fHits->UncheckedAt(hit);
      cpvHit->Print();
    }
  }
}

//______________________________________________________________________________

void AliPHOSCPVModule::MakeBranch(Int_t i)
{
  //
  // Create a new branch for a CPV module #i in the current Root Tree
  //

  char branchname[10];
  sprintf(branchname,"CPV%d",i);
  gAlice->TreeH()->Branch(branchname,&fHits, 1000);
}

//_____________________________________________________________________________
void AliPHOSCPVModule::SetTreeAddress(Int_t i)
{
  //
  // Set branch address for the Hits Tree for a CPV module #i
  //

  TBranch *branch;
  char branchname[20];
  TTree *treeH = gAlice->TreeH();
  if (treeH){
    sprintf(branchname,"CPV%d",i);
    branch = treeH->GetBranch(branchname);
    if (branch) branch->SetAddress(&fHits);
  }
}
