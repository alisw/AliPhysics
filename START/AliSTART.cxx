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

#include "AliLog.h"
#include "AliMC.h"
#include "AliLoader.h"
#include "AliRun.h"

#include "AliSTART.h"
#include "AliSTARTLoader.h"
#include "AliSTARTdigit.h"
#include "AliSTARThit.h"
#include "AliSTARTDigitizer.h"
#include "AliSTARTRawData.h"

ClassImp(AliSTART)

  //static  AliSTARTdigit *digits; 

//_____________________________________________________________________________
AliSTART::AliSTART()
{
  //
  // Default constructor for class AliSTART
  //
  fIshunt   = 1;
  fHits     = 0;
  fDigits   = 0;
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
  fIshunt     =  1;
  fIdSens   =  0;
  SetMarkerColor(kRed);
}

//_____________________________________________________________________________
AliSTART::~AliSTART() {
  if (fHits) {
    fHits->Delete();
    delete fHits;
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
  char branchname[20];
  sprintf(branchname,"%s",GetName());

  const char *cH = strstr(option,"H");
  
  if (cH && fLoader->TreeH())
  {
     if (fHits == 0x0) fHits  = new TClonesArray("AliSTARThit",  405);
  } 
  
  AliDetector::MakeBranch(option);
}    

//_____________________________________________________________________________
void AliSTART::ResetHits()
{
  AliDetector::ResetHits();
  
}

//_____________________________________________________________________________
void AliSTART::SetTreeAddress()
{

  TTree    *treeH;
  treeH = TreeH();
  
  if (treeH)
    {
      if (fHits == 0x0) fHits  = new TClonesArray("AliSTARThit",  405);
    }
    
  AliDetector::SetTreeAddress();
  
}

//______________________________________________________________________
AliLoader* AliSTART::MakeLoader(const char* topfoldername)
{ 

  AliDebug(2,Form(" Creating AliSTARTLoader "));
  fLoader = new AliSTARTLoader(GetName(), topfoldername);
  return fLoader;
}

//_____________________________________________________________________________
AliDigitizer* AliSTART::CreateDigitizer(AliRunDigitizer* manager) const
{
  return new AliSTARTDigitizer(manager);
}
//____________________________________________________________________________
void AliSTART::Digits2Raw()
{
//
// Starting from the START digits, writes the Raw Data objects
//
  AliSTARTLoader* pStartLoader = (AliSTARTLoader*)fLoader;
  pStartLoader ->LoadDigits();
  AliSTARTdigit* fDigits=pStartLoader->Digits();
  AliSTARTRawData rawWriter;
  rawWriter.SetVerbose(0);

  AliDebug(2,Form(" Formatting raw data for START "));
  
  rawWriter.RawDataSTART (fDigits);

   pStartLoader->UnloadDigits();

}
