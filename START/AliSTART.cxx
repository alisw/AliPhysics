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

#include <Riostream.h>

#include <TFile.h>
#include <TGeometry.h>
#include <TMath.h>
#include <TNode.h>
#include <TParticle.h>
#include <TRandom.h>
#include <TTUBE.h>
#include <TVirtualMC.h>
#include <AliESD.h>

#include "AliLoader.h"
#include "AliRun.h"
#include "AliSTART.h"
#include "AliSTARTLoader.h"
#include "AliSTARTdigit.h"
#include "AliSTARThit.h"
#include "AliSTARThitPhoton.h"
#include "AliSTARTvertex.h"
#include "AliMC.h"
#include "AliSTARTDigitizer.h"

ClassImp(AliSTART)

static  AliSTARTdigit *digits; 

//_____________________________________________________________________________
AliSTART::AliSTART()
{
  //
  // Default constructor for class AliSTART
  //
  fIshunt   = 1;
  fHits     = 0;
  fDigits   = 0;
  fPhotons  = 0;
}
 
//_____________________________________________________________________________
AliSTART::AliSTART(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // Standard constructor for START Detector
  //

  
  //
  // Initialise Hit array
  fHits       = new TClonesArray("AliSTARThit",  405);
  gAlice->GetMCApp()->AddHitList(fHits);

  fPhotons  = new TClonesArray("AliSTARThitPhoton", 10000);
  gAlice->GetMCApp()->AddHitList (fPhotons);
  if (GetDebug()>2) cout<<" Debug "<<endl;
  fIshunt     =  1;
  fIdSens   =  0;
  fNPhotons =  0;
  SetMarkerColor(kRed);
}

//_____________________________________________________________________________
AliSTART::~AliSTART() {
  if (fHits) {
    fHits->Delete();
    delete fHits;
  }
  if (fPhotons) {
    fPhotons->Delete();
    delete fPhotons;
  }
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
void AliSTART::AddHitPhoton(Int_t track, Int_t *vol, Float_t *hits)
{
  //  Add a START hit of photons
  
  TClonesArray &lhits = *fPhotons;
  new(lhits[fNPhotons++]) AliSTARThitPhoton(fIshunt,track,vol,hits);
}

//_____________________________________________________________________________

void AliSTART::AddDigit(Int_t * /*tracks*/, Int_t * /*digits*/)
{
  
  //  Add a START digit to the list. Dummy function.
  
}

//_____________________________________________________________________________
void AliSTART::BuildGeometry()
{
  //
  // Build simple ROOT TNode geometry for event display
  //
  TNode *node, *top;
  const int kColorSTART  = 19;

  top=gAlice->GetGeometry()->GetNode("alice");

  // START define the different volumes
  new TRotMatrix("rotx999","rot999",  90,0,90,90,180,0);

  new TTUBE("S_0ST1","START  volume 1","void",5.,10.7,5.3);
  top->cd();
  node = new TNode("0ST1","0ST01","S_0ST1",0,0,-69.7,"");
  node->SetLineColor(kColorSTART);
  fNodes->Add(node);

  new TTUBE("S_0ST2","START volume 2","void",5.,10.7,5.3);
  top->cd();
  node = new TNode("0ST2","0ST2","S_0ST2",0,0,350,"rotx999");
  node->SetLineColor(kColorSTART);
  fNodes->Add(node);
}
 
//_____________________________________________________________________________
Int_t AliSTART::DistanceToPrimitive(Int_t /*px*/, Int_t /*py*/)
{
  //
  // Calculate the distance from the mouse to the START on the screen
  // Dummy routine
  //
  return 9999;
}
 
//-------------------------------------------------------------------------
void AliSTART::Init()
{
  //
  // Initialis the START after it has been built
  Int_t i;
  //
  if(fDebug) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" START_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the START initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}

//---------------------------------------------------------------------------
void AliSTART::MakeBranch(Option_t* option)
{
  //
  // Specific START branches
  //
  // Create Tree branches for the START.
  Int_t buffersize = 4000;
  char branchname[20];
  sprintf(branchname,"%s",GetName());


  const char *cD = strstr(option,"D");
  const char *cH = strstr(option,"H");
  
  if (cH && fLoader->TreeH())
  {
     if (fPhotons == 0x0) fPhotons  = new TClonesArray("AliSTARThitPhoton", 10000);
     sprintf (branchname, "%shitPhoton", GetName());
     MakeBranchInTree (fLoader->TreeH(), branchname, &fPhotons, 50000, 0);
     if (fHits == 0x0) fHits  = new TClonesArray("AliSTARThit",  405);
  } 
  
  AliDetector::MakeBranch(option);

  if (cD) {
    digits = new AliSTARTdigit();
    MakeBranchInTree(fLoader->TreeD(), branchname, "AliSTARTdigit", digits, buffersize, 1, 0);
  } 
}    

//_____________________________________________________________________________
void AliSTART::ResetHits()
{
  AliDetector::ResetHits();
  
  fNPhotons = 0;
  if (fPhotons)  fPhotons->Clear();
}

//_____________________________________________________________________________
void AliSTART::SetTreeAddress()
{
  TBranch  *branch;
  TTree    *treeH;
 
  
  treeH = TreeH();
  
  if (treeH)
    {
      if (fPhotons == 0x0) fPhotons  = new TClonesArray("AliSTARThitPhoton", 10000);
      branch = treeH->GetBranch("STARThitPhoton");
      if (branch)  branch->SetAddress(&fPhotons);
      if (fHits == 0x0) fHits  = new TClonesArray("AliSTARThit",  405);
    }
    
  AliDetector::SetTreeAddress();
  
}

//______________________________________________________________________
AliLoader* AliSTART::MakeLoader(const char* topfoldername)
{ 
  Info("MakeLoader", "Creating AliSTARTLoader. Top folder is %s.", topfoldername);
  fLoader = new AliSTARTLoader(GetName(), topfoldername);
  return fLoader;
}

//_____________________________________________________________________________
AliDigitizer* AliSTART::CreateDigitizer(AliRunDigitizer* manager) const
{
  return new AliSTARTDigitizer(manager);
}

//_____________________________________________________________________________
void AliSTART::FillESD(AliESD* pESD)  const
{
  AliSTARTvertex reco;
  reco.Reconstruct(fLoader->GetRunLoader(), pESD);
}

