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

 //////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Forward Multiplicity Detector based on Silicon plates                    //
//  This class contains the base procedures for the Forward Multiplicity     //
//  detector                                                                 //
//  Detector consists of 5 Si volumes covered pseudorapidity interval         //
//  from 1.7 to 5.1.                                                         //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliFMDClass.gif">
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

#define DEBUG

#include <Riostream.h>
#include <stdlib.h>

#include <TClonesArray.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TTree.h>
#include <TVirtualMC.h>

#include "AliDetector.h"
#include "AliFMDdigit.h"
#include "AliFMDhit.h"
#include "AliFMDv1.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliFMDDigitizer.h"

ClassImp (AliFMD)
  //_____________________________________________________________________________
AliFMD::AliFMD ():AliDetector ()
{
  //
  // Default constructor for class AliFMD
  //
  fIshunt = 0;
  fHits     = 0;
  fDigits   = 0;
}

//_____________________________________________________________________________
AliFMD::AliFMD (const char *name, const char *title):
AliDetector (name, title)
{
  //
  // Standard constructor for Forward Multiplicity Detector
  //

  //
  // Initialise Hit array
  fHits = new TClonesArray ("AliFMDhit", 1000);
  // Digits for each Si disk
  fDigits = new TClonesArray ("AliFMDdigit", 1000);
  gAlice->GetMCApp()->AddHitList (fHits);

  fIshunt = 0;
  //  fMerger = 0;
  SetMarkerColor (kRed);
}

//-----------------------------------------------------------------------------
AliFMD::~AliFMD ()
{
  //destructor for base class AliFMD
  if (fHits)
    {
      fHits->Delete ();
      delete fHits;
      fHits = 0;
    }
  if (fDigits)
    {
      fDigits->Delete ();
      delete fDigits;
      fDigits = 0;
    }

}

//_____________________________________________________________________________
void AliFMD::AddHit (Int_t track, Int_t * vol, Float_t * hits)
{
  //
  // Add a hit to the list
  //
  TClonesArray & lhits = *fHits;
  new (lhits[fNhits++]) AliFMDhit (fIshunt, track, vol, hits);
}

//_____________________________________________________________________________
void AliFMD::AddDigit (Int_t * digits)
{
  // add a real digit - as coming from data

  if (fDigits == 0x0) fDigits = new TClonesArray ("AliFMDdigit", 1000);  
  TClonesArray & ldigits = *fDigits;
  new (ldigits[fNdigits++]) AliFMDdigit (digits);
}

//_____________________________________________________________________________
void AliFMD::BuildGeometry ()
{
  //
  // Build simple ROOT TNode geometry for event display
  //
  TNode *node, *top;
  const int kColorFMD = 5;
  //
  top = gAlice->GetGeometry ()->GetNode ("alice");

  // FMD define the different volumes
  new TRotMatrix ("rot901", "rot901", 90, 0, 90, 90, 180, 0);

  new TTUBE ("S_FMD0", "FMD  volume 0", "void", 4.2, 17.2, 1.5);
  top->cd ();
  node = new TNode ("FMD0", "FMD0", "S_FMD0", 0, 0, -62.8, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD1", "FMD  volume 1", "void", 15.4, 28.4, 1.5);
  top->cd ();
  node = new TNode ("FMD1", "FMD1", "S_FMD1", 0, 0, -75.2, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD2", "FMD  volume 2", "void", 4.2, 17.2, 1.5);
  top->cd ();
  node = new TNode ("FMD2", "FMD2", "S_FMD2", 0, 0, 83.2, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD3", "FMD  volume 3", "void", 15.4, 28.4, 1.5);
  top->cd ();
  node = new TNode ("FMD3", "FMD3", "S_FMD3", 0, 0, 75.2, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD4", "FMD  volume 4", "void", 4.2, 17.2, 1.5);
  top->cd ();
  node = new TNode ("FMD4", "FMD4", "S_FMD4", 0, 0, 340, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);
}

//_____________________________________________________________________________
const Int_t AliFMD::DistanceToPrimitive (Int_t /*px*/, Int_t /*py*/)
{
  //
  // Calculate the distance from the mouse to the FMD on the screen
  // Dummy routine
  //
  return 9999;
}

//___________________________________________
void AliFMD::ResetHits ()
{
  // Reset number of clusters and the cluster array for this detector
  AliDetector::ResetHits ();
}

//____________________________________________
void AliFMD::ResetDigits ()
{
  //
  // Reset number of digits and the digits array for this detector
  AliDetector::ResetDigits ();
  //
}

//-------------------------------------------------------------------------
void  AliFMD::Init ()
{
  //
  // Initialis the FMD after it has been built
  Int_t i;
  //
  if (fDebug)
    {
      printf ("\n%s: ", ClassName ());
      for (i = 0; i < 35; i++)
	printf ("*");
      printf (" FMD_INIT ");
      for (i = 0; i < 35; i++)
	printf ("*");
      printf ("\n%s: ", ClassName ());
      //
      // Here the FMD initialisation code (if any!)
      for (i = 0; i < 80; i++)
	printf ("*");
      printf ("\n");
    }
  //
  //
 
}
//---------------------------------------------------------------------
void AliFMD::MakeBranch (Option_t * option)
{
  // Create Tree branches for the FMD.
  char branchname[10];
  const Int_t kBufferSize = 16000;
  sprintf (branchname, "%s", GetName ());
  
  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  
  if (cH && (fHits == 0x0)) fHits = new TClonesArray ("AliFMDhit", 1000);

  AliDetector::MakeBranch (option);
  
  if (cD){
    if (fDigits == 0x0) fDigits = new TClonesArray ("AliFMDdigit", 1000);  
    MakeBranchInTree(fLoader->TreeD(), branchname,&fDigits, kBufferSize, 0);
  }

}

//_____________________________________________________________________________
void AliFMD::SetTreeAddress ()
{
  // Set branch address for the Hits and Digits Tree.

  if (fLoader->TreeH() && (fHits == 0x0)) 
    fHits = new TClonesArray ("AliFMDhit", 1000);  

  AliDetector::SetTreeAddress ();

  TBranch *branch;
  TTree *treeD = fLoader->TreeD();

  if (treeD)
    {
      if (fDigits == 0x0) fDigits = new TClonesArray ("AliFMDdigit", 1000);
      branch = treeD->GetBranch ("FMD");
      if (branch)
       branch->SetAddress (&fDigits);
    }
}



//-----------------------------------------------------------------------

void AliFMD::MakeBranchInTreeD(TTree *treeD, const char *file)
{
    //
    // Create TreeD branches for the FMD
    //
    const Int_t kBufferSize = 4000;
    char branchname[20];
    sprintf(branchname,"%s",GetName());	
    if(treeD)
     {
       MakeBranchInTree(treeD,  branchname,&fDigits, kBufferSize, file);
     }
}

//____________________________________________________________________________
AliDigitizer* AliFMD::CreateDigitizer(AliRunDigitizer* manager) const
{
  return new AliFMDDigitizer(manager);
}
