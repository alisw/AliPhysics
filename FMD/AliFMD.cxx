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
 //////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Forward Multiplicity Detector based on Silicon plates                    //
//  This class contains the base procedures for the Forward Multiplicity     //
//  detector                                                                 //
//  Detector consists of 6 Si volumes covered pseudorapidity interval         //
//  from 1.6 to 6.0.                                                         //
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
#include <TMath.h>
#include <TGeometry.h>
#include <TTUBE.h>
#include <TTree.h>
#include <TNode.h>
#include <TFile.h>

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "AliFMDv1.h"
#include "AliRun.h"
#include "AliDetector.h"
#include <Riostream.h>
#include "AliMagF.h"
#include "AliFMDhit.h"
#include "AliFMDdigit.h"
#include "AliFMDReconstruction.h"
#include "AliFMDReconstParticles.h"
#include <stdlib.h>

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
  fSDigits  = 0;
  fReconParticles=0; 
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
  fSDigits = new TClonesArray ("AliFMDdigit", 1000);
  fReconParticles=new TClonesArray("AliFMDReconstParticles",1000); 
  gAlice->AddHitList (fHits);

  fIshunt = 0;
  fIdSens1 = 0;
  fIdSens2 = 0;
  fIdSens3 = 0;
  fIdSens4 = 0;
  fIdSens5 = 0;
  //  fMerger = 0;
  SetMarkerColor (kRed);
}

//-----------------------------------------------------------------------------
AliFMD::~AliFMD ()
{
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
  if (fSDigits)
    {
      fSDigits->Delete ();
      delete fSDigits;
      fSDigits = 0;
    }
  if (fReconParticles)
    {
      fReconParticles->Delete ();
      delete fReconParticles;
      fReconParticles = 0;
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


  TClonesArray & ldigits = *fDigits;
  new (ldigits[fNdigits++]) AliFMDdigit (digits);

}
//_____________________________________________________________________________
void AliFMD::AddSDigit (Int_t * digits)
{
  // add a real digit - as coming from data

  TClonesArray & ldigits = *fSDigits;
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
  node = new TNode ("FMD0", "FMD0", "S_FMD0", 0, 0, 62.8, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD1", "FMD  volume 1", "void", 15.4, 28.4, 1.5);
  top->cd ();
  node = new TNode ("FMD1", "FMD1", "S_FMD1", 0, 0, 75.2, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD2", "FMD  volume 2", "void", 4.2, 17.2, 1.5);
  top->cd ();
  node = new TNode ("FMD2", "FMD2", "S_FMD2", 0, 0, -83.2, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD3", "FMD  volume 3", "void", 15.4, 28.4, 1.5);
  top->cd ();
  node = new TNode ("FMD3", "FMD3", "S_FMD3", 0, 0, -75.2, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD4", "FMD  volume 4", "void", 4.2, 17.2, 1.5);
  top->cd ();
  //  node = new TNode("FMD4","FMD4","S_FMD4",0,0,-270,"");
  node = new TNode ("FMD4", "FMD4", "S_FMD4", 0, 0, -340, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);
}

//_____________________________________________________________________________
Int_t AliFMD::DistanceToPrimitive (Int_t px, Int_t py)
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
    fIdSens1 = gMC->VolId ("GRN1");	//Si sensetive volume
    fIdSens2 = gMC->VolId ("GRN2");	//Si sensetive volume
    fIdSens3 = gMC->VolId ("GRN3");	//Si sensetive volume
    fIdSens4 = gMC->VolId ("GRN4");	//Si sensetive volume
    fIdSens5 = gMC->VolId ("GRN5");	//Si sensetive volume

}

//---------------------------------------------------------------------
void AliFMD::MakeBranch (Option_t * option, const char *file)
{
  // Create Tree branches for the FMD.
  char branchname[10];
  const Int_t kBufferSize = 16000;
  sprintf (branchname, "%s", GetName ());
  
  AliDetector::MakeBranch (option, file);
  const char *cD = strstr(option,"D");
  const char *cR = strstr(option,"R");
  const char *cS = strstr(option,"S");
  
  if (cS){

    MakeBranchInTree(gAlice->TreeS(), 
		     branchname,&fSDigits, 
    		     kBufferSize, file);
  }
  if (cD){

    MakeBranchInTree(gAlice->TreeD(), 
		     branchname,&fDigits, 
    		     kBufferSize, file);
    cout<<" tree "<<gAlice->TreeD()<<" "<<branchname<<" "<<&fDigits<<endl;
  }
  if (cR){
    MakeBranchInTree(gAlice->TreeR(), 
		     branchname,&fReconParticles,
    		     kBufferSize, file);
  }
  
}

//_____________________________________________________________________________
void AliFMD::SetTreeAddress ()
{
  // Set branch address for the Hits and Digits Tree.
  char branchname[30];
  AliDetector::SetTreeAddress ();

  TBranch *branch;
  TTree *treeD = gAlice->TreeD ();


  if (treeD)
    {
      if (fDigits)
	{
	  branch = treeD->GetBranch (branchname);
	  if (branch)
	    branch->SetAddress (&fDigits);
	}

    }
  if (fSDigits)
    //  fSDigits->Clear ();

  if (gAlice->TreeS () && fSDigits)
    {
      branch = gAlice->TreeS ()->GetBranch ("FMD");
      if (branch)
	branch->SetAddress (&fSDigits);
    }

  if (gAlice->TreeR() && fReconParticles) 
    {
      branch = gAlice->TreeR()->GetBranch("FMD"); 
      if (branch) branch->SetAddress(&fReconParticles) ;
    }   
}

//---------------------------------------------------------------------

void AliFMD::SetRingsSi1(Int_t ringsSi1)
{
  //  fRingsSi1=ringsSi1;
  fRingsSi1=768;
}
void AliFMD::SetSectorsSi1(Int_t sectorsSi1)
{
  fSectorsSi1=20;
}
void AliFMD::SetRingsSi2(Int_t ringsSi2)
{
  fRingsSi2=384;
}
void AliFMD::SetSectorsSi2(Int_t sectorsSi2)
{
  fSectorsSi2=40;
}

//---------------------------------------------------------------------
/*
void AliFMD::SDigits2Digits() 
{
  cout<<"AliFMD::SDigits2Digits"<<endl; 
    if (!fMerger) {
      fMerger = new AliFMDMerger();
    }
    
    fMerger ->SetRingsSi1(fRingsSi1);
    fMerger->SetRingsSi2(fRingsSi2);
    fMerger ->SetSectorsSi1(fSectorsSi1);
    fMerger ->SetSectorsSi2(fSectorsSi2);
     
    fMerger->Init();
    cout<<"AliFMD::SDigits2Digits Init"<<endl; 
    fMerger->Digitise();
    cout<<"AliFMD::SDigits2Digits Digitise() "<<endl; 
 

    }

    //---------------------------------------------------------------------
void   AliFMD::SetMerger(AliFMDMerger* merger)
{
// Set pointer to merger
    fMerger = merger;
}

AliFMDMerger*  AliFMD::Merger()
{
// Return pointer to merger
    return fMerger;
}
*/
//---------------------------------------------------------------------



void
AliFMD::Eta2Radius (Float_t eta, Float_t zDisk, Float_t * radius)
{
  Float_t expEta = TMath::Exp (-eta);
  Float_t theta = TMath::ATan (expEta);
  theta = 2. * theta;
  Float_t rad = zDisk * (TMath::Tan (theta));
  *radius = rad;

  if (fDebug)
    printf ("%s: eta %f radius %f\n", ClassName (), eta, rad);
}

//---------------------------------------------------------------------

void AliFMD::Hits2SDigits ()
{

  //#ifdef DEBUG
  cout<<"ALiFMD::Hits2SDigits> start...\n";
  //#endif
  
  char * fileSDigits = "FMD.SDigits.root";
  char * fileHeader = 0;
  AliFMDSDigitizer * sd = new AliFMDSDigitizer(fileHeader,fileSDigits) ;
  sd->SetRingsSi1(fRingsSi1);
  sd->SetRingsSi2(fRingsSi2);
  sd->SetSectorsSi1(fSectorsSi1);
  sd->SetSectorsSi2(fSectorsSi2);
  //  sd->SetEventNumber(fEvNrSig);
  sd->Exec("") ;
  
  delete sd ;
  
}
//-----------------------------------------------------------------------

void AliFMD::Digits2Reco()
{
  char * fileReconParticles=0;
  char * fileHeader=0;
  AliFMDReconstruction * reconstruction =
    new AliFMDReconstruction(fileHeader,fileReconParticles) ;
  reconstruction->Exec("");
  delete  reconstruction;
}
//-----------------------------------------------------------------------

void AliFMD::MakeBranchInTreeD(TTree *treeD, const char *file)
{
    //
    // Create TreeD branches for the MUON.
    //

    const Int_t kBufferSize = 4000;
    char branchname[20];
    

    sprintf(branchname,"%s",GetName());	
    if(treeD){
    MakeBranchInTree(treeD, 
		     branchname,&fDigits, 
    		     kBufferSize, file);
    }
}

