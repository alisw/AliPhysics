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
//                                                                           //
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
#include "AliMC.h"
#include "AliDetector.h"
#include <iostream.h>
#include <fstream.h>
#include "AliMagF.h"
#include "AliFMDhit.h"
#include "AliFMDdigit.h"
#include "AliFMDReconstruction.h"
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
  gAlice->AddHitList (fHits);

  fIshunt = 0;
  fIdSens1 = 0;

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

  // printf("AddDigit\n");

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
  const int kColorFMD = 7;
  //
  top = gAlice->GetGeometry ()->GetNode ("alice");

  // FMD define the different volumes
  new TRotMatrix ("rot901", "rot901", 90, 0, 90, 90, 180, 0);

  new TTUBE ("S_FMD0", "FMD  volume 0", "void", 4.73, 17.7, 1.5);
  top->cd ();
  node = new TNode ("FMD0", "FMD0", "S_FMD0", 0, 0, 64, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD1", "FMD  volume 1", "void", 23.4, 36., 1.5);
  top->cd ();
  node = new TNode ("FMD1", "FMD1", "S_FMD1", 0, 0, 85, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD2", "FMD  volume 2", "void", 4.73, 17.7, 1.5);
  top->cd ();
  node = new TNode ("FMD2", "FMD2", "S_FMD2", 0, 0, -64, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD3", "FMD  volume 3", "void", 23.4, 36., 1.5);
  top->cd ();
  node = new TNode ("FMD3", "FMD3", "S_FMD3", 0, 0, -85, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  new TTUBE ("S_FMD4", "FMD  volume 4", "void", 5, 15, 0.015);
  top->cd ();
  //  node = new TNode("FMD4","FMD4","S_FMD4",0,0,-270,"");
  node = new TNode ("FMD4", "FMD4", "S_FMD4", 0, 0, -270, "");
  node->SetLineColor (kColorFMD);
  fNodes->Add (node);

  /* 
     new TTUBE("S_FMD5","FMD  volume 5","void",5,14,0.015);
     top->cd();
     node = new TNode("FMD5","FMD5","S_FMD5",0,0,-630,"");
     node->SetLineColor(kColorFMD);
     fNodes->Add(node);
   */
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
  AliDetector::ResetHits ();
  //
}

//-------------------------------------------------------------------------
void  AliFMD::Init ()
{
  //
  // Initialis the FMD after it has been built
  Int_t i;
  AliMC *pMC = AliMC::GetMC ();
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
  if (IsVersion () != 0)
    fIdSens1 = pMC->VolId ("GRIN");	//Si sensetive volume
  else
    fIdSens1 = pMC->VolId ("GFSI");	//Si sensetive volume

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
    fSDigits->Clear ();

  if (gAlice->TreeS () && fSDigits)
    {
      branch = gAlice->TreeS ()->GetBranch ("FMD");
      if (branch)
	branch->SetAddress (&fSDigits);
    }

  if(fReconParticles)
    fReconParticles->Clear();
  if (gAlice->TreeR()) 
    {
      branch = gAlice->TreeR()->GetBranch("FMD"); 
      if (branch) branch->SetAddress(&fReconParticles) ;
    } 

}

//---------------------------------------------------------------------

void AliFMD::SDigits2Digits() 
{
  cout<<"AliFMD::SDigits2Digits"<<endl; 
    if (fMerger) {
      cout<<"AliFMD::SDigits2Digits fMerger"<<fMerger<<endl;
      fMerger->Init();
      cout<<"AliFMD::SDigits2Digits Init"<<endl; 
  
      fMerger->Digitise();
    }

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

  AliFMD *FMD = (AliFMD *) gAlice->GetDetector ("FMD");

  if (fNevents == 0)
    fNevents = (Int_t) gAlice->TreeE ()->GetEntries ();

  for (Int_t ievent = 0; ievent < fNevents; ievent++)
    {
      gAlice->GetEvent (ievent);
      if (gAlice->TreeH () == 0)
	return;
      if (gAlice->TreeS () == 0)
	gAlice->MakeTree ("S");



      Int_t nSdigits = 0;
      
            //Make branches
      char branchname[20];
       sprintf (branchname, "%s", FMD->GetName ());
      //Make branch for digits
      FMD->MakeBranch ("S");
          
    
       //Now made SDigits from hits, for PHOS it is the same
      Int_t volume, sector, ring, charge;
      Float_t e;
      Float_t de[10][20][150];
      Int_t ivol, isec, iring;
      Int_t hit, nbytes;
      TParticle *particle;
      AliFMDhit *fmdHit;
      TClonesArray *FMDhits = FMD->Hits ();

      // Event ------------------------- LOOP  

      for (ivol = 1; ivol <= 5; ivol++)
	for (isec = 1; isec <= 16; isec++)
	  for (iring = 1; iring <= 128; iring++)
	    de[ivol][isec][iring] = 0;

      if (FMD)
	{
	  FMDhits = FMD->Hits ();
	  TTree *TH = gAlice->TreeH ();
	  Stat_t ntracks = TH->GetEntries ();
	  for (Int_t track = 0; track < ntracks; track++)
	    {
	      gAlice->ResetHits ();
	      nbytes += TH->GetEvent (track);
	      particle = gAlice->Particle (track);
	      Int_t nhits = FMDhits->GetEntriesFast ();

	      for (hit = 0; hit < nhits; hit++)
		{
		  fmdHit = (AliFMDhit *) FMDhits->UncheckedAt (hit);

		  volume = fmdHit->Volume ();
		  sector = fmdHit->NumberOfSector ();
		  ring = fmdHit->NumberOfRing ();
		  e = fmdHit->Edep ();
		  de[volume][sector][ring] = de[volume][sector][ring] + e;
		}		//hit loop
	    }			//track loop
	}			//if FMD


      Int_t digit[5];
      Float_t I = 1.664 * 0.04 * 2.33 / 22400;	// = 0.69e-6;
      for (ivol = 1; ivol <= 5; ivol++)
	{
	  for (isec = 1; isec <= 16; isec++)
	    {
	      for (iring = 1; iring <= 128; iring++)
		{
		      digit[0] = ivol;
		      digit[1] = isec;
		      digit[2] = iring;
		      charge = Int_t (de[ivol][isec][iring] / I);
		      digit[3] = charge;
		      //		      if (charge!=0) cout<<" charge "<<charge<<endl;
		      //dinamic diapason from MIP(0.155MeV) to 30MIP(4.65MeV)
		      //1024 ADC channels 
		      Float_t channelWidth = (22400 * 30) / 1024;
		      digit[4] = Int_t (digit[3] / channelWidth);

		      new ((*fSDigits)[nSdigits++]) AliFMDdigit (digit);

		}		// iring loop
	    }			//sector loop
	}			// volume loop
      
      gAlice->TreeS ()->Fill ();
      gAlice->TreeS ()->Print ();

    }				//event loop

}
//-----------------------------------------------------------------------

void AliFMD::Digits2Reco()
{
#ifdef DEBUG
  cout<<"ALiFMD::Digits2Reco> start...";
#endif
  char * fileReconParticles=0;
  char * fileHeader=0;
  AliFMDReconstruction * reconstruction =
    new AliFMDReconstruction(fileHeader,fileReconParticles) ;
  fReconParticles=new TClonesArray("AliFMDReconstParticles",1000);
  reconstruction->Exec(fReconParticles,"");
  delete  reconstruction;
}

