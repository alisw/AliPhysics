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
Revision 1.4  1999/11/12 15:04:00  fca
Modifications from A.Maevskaya

Revision 1.3  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  START (T-Zero) Detector                                            //
//  This class contains the base procedures for the START     //
//  detector                                                                 //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliSTARTClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:Alla.Maevskaia@cern.ch">Alla Maevskaia</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>

#include <TTUBE.h>
#include <TNode.h>
#include <TGeometry.h>
#include "AliRun.h"
#include "AliSTART.h"
#include <iostream.h>
#include <fstream.h>
#include "AliMC.h"
 
ClassImp(AliSTART)
 
//_____________________________________________________________________________
AliSTART::AliSTART()
{
  //
  // Default constructor for class AliSTART
  //
  fIshunt   = 0;
}
 
//_____________________________________________________________________________
AliSTART::AliSTART(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // Standard constructor for START Detector
  //
 
  AliModule *fmd = gAlice->GetModule("FMD");
  if(fmd) {
    Int_t fmdversion = fmd->IsVersion();
    if(fmdversion==0 || fmdversion==1) {
      Error("ctor","Versions 0 and 1 of FMD incompatible with START\n");
      exit(1);
    }
  }
  //
  // Initialise Hit array
  fHits   = new TClonesArray("AliSTARThit",  405);
  gAlice->AddHitList(fHits);
  
  fIshunt     =  0;
  fIdSens1    =  0;

  SetMarkerColor(kRed);
}
 
//_____________________________________________________________________________
void AliSTART::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a START hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliSTARThit(fIshunt,track,vol,hits);
}
 
//_____________________________________________________________________________
void AliSTART::BuildGeometry()
{
  //
  // Build simple ROOT TNode geometry for event display
  //
  TNode *Node, *Top;
  const int kColorSTART  = 19;
  //
  Top=gAlice->GetGeometry()->GetNode("alice");

  // START define the different volumes
  new TRotMatrix("rot999","rot999",  90,0,90,90,180,0);

  new TTUBE("S_STR1","START  volume 1","void",5.,10.7,5.3);
  Top->cd();
  Node = new TNode("STR1","STR1","S_STR1",0,0,75.,"");
  Node->SetLineColor(kColorSTART);
  fNodes->Add(Node);

  new TTUBE("S_STR2","START volume 2","void",5.,10.7,5.3);
  Top->cd();
  Node = new TNode("STR2","STR2","S_STR2",0,0,-75,"rot999");
  Node->SetLineColor(kColorSTART);
  fNodes->Add(Node);
}
 
//_____________________________________________________________________________
Int_t AliSTART::DistanceToPrimitive(Int_t px, Int_t py)
{
  //
  // Calculate the distance from the mouse to the START on the screen
  // Dummy routine
  //
  return 9999;
}
 
//_____________________________________________________________________________

//-------------------------------------------------------------------------
void AliSTART::Init()
{
  //
  // Initialis the START after it has been built
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" START_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  // Here the START initialisation code (if any!)
  for(i=0;i<80;i++) printf("*");
  printf("\n");
 //
 //
 fIdSens1=gMC->VolId("PTOP");

}

//---------------------------------------------------------------------------
void AliSTART::MakeBranch(Option_t* option)
{
  
  // Create Tree branches for the START.
  Int_t buffersize = 4000;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option);

  char *D = strstr(option,"D");

  if (fDigits   && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fDigits, buffersize);
    printf("Making Branch %s for digits\n",branchname);
  }
  
}    
 
ClassImp(AliSTARThit)
 
//_____________________________________________________________________________
AliSTARThit::AliSTARThit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Add a START hit
  //
  
//  Int_t i;
  fVolume = vol[0];
  fPmt=vol[1];
//printf("fvolume %d\n",fVolume);
//printf("fpmt %d\n",fPmt);

  fX=hits[0];
  fY=hits[1];
  fZ=hits[2];
  fEdep=hits[3];
  fEtot=hits[4];
  fParticle=Int_t (hits[5]);
  fTime=hits[6];

//for (i=0; i<=6; i++) {printf("Hits up %f\n",hits[i]);} 
}

// ClassImp(AliSTARTdigit)
