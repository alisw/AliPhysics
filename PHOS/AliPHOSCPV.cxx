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
//  Manager and hits classes for set          //
//  Charged Particle Veto (CPV)               //
//                                            //
//  Author: Yuri Kharlov, IHEP, Protvino      //
//  e-mail: Yuri.Kharlov@cern.ch              //
//  Last modified: 28 September 2000          //
////////////////////////////////////////////////
 
// --- ROOT system ---
#include <TTree.h>

// --- Standard library ---
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream.h>

// --- galice header files ---
#include "AliPHOSCPV.h"
#include "AliRun.h"

//==============================================================================
//                              CPVModule
//==============================================================================

ClassImp(CPVModule)

//______________________________________________________________________________

CPVModule::CPVModule(void) {
  //
  // Allocate an array of hits
  //

  if ( NULL==(fHits=new TClonesArray("CPVHit",100)) ) {
    Error("CPV","Can not create array of hits");
    exit(1);
  }
}

//______________________________________________________________________________

CPVModule::~CPVModule(void)
{
  Clear();
}

//______________________________________________________________________________

void CPVModule::Clear(Option_t *opt="")
{
// Clear hit information

  fHits  -> Clear(opt);
}

//______________________________________________________________________________

void CPVModule::AddHit(TLorentzVector p, Float_t *xy, Int_t ipart)
{
  //
  // Add this hit to the hit list in CPV detector.
  //

  TClonesArray &lhits = *(TClonesArray *)fHits;
  new(lhits[fHits->GetEntriesFast()]) CPVHit(p,xy,ipart);
}

//______________________________________________________________________________

void CPVModule::Print(Option_t *opt)
{
  //
  // Print CPVModule information.
  // options:  'p' - print hits in the module
  //

  Int_t nhits,hit;
  if (strcmp(opt,"p")==0) {
    printf ("CPV module has %d hits\n",nhits=fHits->GetEntriesFast());
    for (hit=0;hit<nhits;hit++) {
      CPVHit *cpvHit = (CPVHit*)fHits->UncheckedAt(hit);
      cpvHit->Print();
    }
  }
}

//______________________________________________________________________________

void CPVModule::MakeBranch(Int_t i)
{
  //
  // Create a new branch for a CPV module #i in the current Root Tree
  //

  char branchname[10];
  sprintf(branchname,"CPV%d",i);
  gAlice->TreeH()->Branch(branchname,&fHits, 1000);
}

//_____________________________________________________________________________
void CPVModule::SetTreeAddress(Int_t i)
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

//==============================================================================
//                              CPVHit
//==============================================================================

ClassImp(CPVHit)

//______________________________________________________________________________

CPVHit::CPVHit(TLorentzVector p, Float_t *xy, Int_t ipart)
{
  //
  // Create a CPV hit object
  //

  fMomentum  = p;
  fXhit      = xy[0];
  fYhit      = xy[1];
  fIpart     = ipart;
}

//______________________________________________________________________________
void CPVHit::Print()
{
  //
  // Print CPV hit
  //

  printf("CPV hit: p  = (% .4f, % .4f, % .4f, % .4f) GeV,\n",
	GetMomentum().Px(),GetMomentum().Py(),GetMomentum().Pz(),GetMomentum().E());
  printf("         xy = (%8.4f, %8.4f) cm, ipart = %d\n",
	 fXhit,fYhit,fIpart);
}

//==============================================================================
//                              CPVDigit
//==============================================================================

ClassImp(CPVDigit)

//______________________________________________________________________________

CPVDigit::CPVDigit(Int_t x, Int_t y, Float_t q)
{
  //
  // Create a CPV digit object
  //

  fXpad = x;
  fYpad = y;
  fQpad = q;
}

//______________________________________________________________________________
