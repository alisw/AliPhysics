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

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "AliFMDv1.h"
#include "AliRun.h"
#include "AliMC.h"
#include <iostream.h>
#include <fstream.h>
#include "AliMagF.h"
#include "AliFMDhit.h"
#include "AliFMDdigit.h"
#include <stdlib.h>


ClassImp(AliFMD)
 
  //_____________________________________________________________________________
  AliFMD::AliFMD(): AliDetector()
{
  //
  // Default constructor for class AliFMD
  //
  fIshunt   = 0;
}
 
//_____________________________________________________________________________
AliFMD::AliFMD(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Standard constructor for Forward Multiplicity Detector
  //
 
  //
  // Initialise Hit array
  fHits     = new TClonesArray("AliFMDhit", 1000);
  // Digits for each Si disk
  fDigits   = new TClonesArray("AliFMDdigit", 1000);
  fSDigits  = new TClonesArray("AliFMDdigit", 1000) ; 
  gAlice->AddHitList(fHits);

  fIshunt     =  0;
  fIdSens1    =  0;

  SetMarkerColor(kRed);
}
 

AliFMD::~AliFMD()
{
  if (fHits) {
      fHits->Delete();
      delete fHits;
  }
  delete fDigits ;
  delete fSDigits ;
   
}
//_____________________________________________________________________________
void AliFMD::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a hit to the list
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliFMDhit(fIshunt,track,vol,hits);
}
//_____________________________________________________________________________
void AliFMD::AddDigit( Int_t *digits) 
{
  // add a real digit - as coming from data

  // printf("AddDigit\n");

   TClonesArray &ldigits = *fDigits;
   // new((*fDigits)[fNdigits++]) AliFMDdigit(digits);
  new(ldigits[fNdigits++]) AliFMDdigit(digits);

}
//_____________________________________________________________________________
void AliFMD::BuildGeometry()
{
  //
  // Build simple ROOT TNode geometry for event display
  //
  TNode *node, *top;
  const int kColorFMD  = 7;
  //
  top=gAlice->GetGeometry()->GetNode("alice");

  // FMD define the different volumes
  new TRotMatrix("rot901","rot901",  90, 0, 90, 90, 180, 0);

  new TTUBE("S_FMD0","FMD  volume 0","void",4.73,17.7,1.5);
  top->cd();
  node = new TNode("FMD0","FMD0","S_FMD0",0,0,64,"");
  node->SetLineColor(kColorFMD);
  fNodes->Add(node);

  new TTUBE("S_FMD1","FMD  volume 1","void",23.4,36.,1.5);
  top->cd();
  node = new TNode("FMD1","FMD1","S_FMD1",0,0,85,"");
  node->SetLineColor(kColorFMD);
  fNodes->Add(node);
  
  new TTUBE("S_FMD2","FMD  volume 2","void",4.73,17.7,1.5);
  top->cd();
  node = new TNode("FMD2","FMD2","S_FMD2",0,0,-64,"");
  node->SetLineColor(kColorFMD);
  fNodes->Add(node);

  new TTUBE("S_FMD3","FMD  volume 3","void",23.4,36.,1.5);
  top->cd();
  node = new TNode("FMD3","FMD3","S_FMD3",0,0,-85,"");
  node->SetLineColor(kColorFMD);
  fNodes->Add(node);

  new TTUBE("S_FMD4","FMD  volume 4","void",5,15,0.015);
  top->cd();
  node = new TNode("FMD4","FMD4","S_FMD4",0,0,-270,"");
  node->SetLineColor(kColorFMD);
  fNodes->Add(node);
  
  
  new TTUBE("S_FMD5","FMD  volume 5","void",5,14,0.015);
  top->cd();
  node = new TNode("FMD5","FMD5","S_FMD5",0,0,-630,"");
  node->SetLineColor(kColorFMD);
  fNodes->Add(node);

}
 
//_____________________________________________________________________________
Int_t AliFMD::DistanceToPrimitive(Int_t px, Int_t py)
{
  //
  // Calculate the distance from the mouse to the FMD on the screen
  // Dummy routine
  //
  return 9999;
}
//___________________________________________
void AliFMD::ResetHits()
{
  // Reset number of clusters and the cluster array for this detector
  AliDetector::ResetHits();
}

//____________________________________________
void AliFMD::ResetDigits()
{
    //
    // Reset number of digits and the digits array for this detector
   AliDetector::ResetHits();
   //
}

//-------------------------------------------------------------------------
void AliFMD::Init()
{
  //
  // Initialis the FMD after it has been built
  Int_t i;
  AliMC* pMC = AliMC::GetMC();
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" FMD_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  // Here the FMD initialisation code (if any!)
  for(i=0;i<80;i++) printf("*");
  printf("\n");
  //
  //
  fIdSens1=pMC->VolId("GRIN"); //Si sensetive volume

}
//---------------------------------------------------------------------
void AliFMD::MakeBranch(Option_t* option, char *file)
{
  // Create Tree branches for the FMD.
  const Int_t kBufferSize = 4000;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option,file);

  const char *cD = strstr(option,"D");
 
  if (cD) {
    
     gAlice->MakeBranchInTree(gAlice->TreeD(), 
                                   branchname, &fDigits, kBufferSize, file) ; 	  

    printf("Making Branch %s for digits\n",branchname);
    gAlice->TreeD()->Print();
  }
}
 
//_____________________________________________________________________________
void AliFMD::SetTreeAddress()
{
  // Set branch address for the Hits and Digits Tree.
  char branchname[30];
  AliDetector::SetTreeAddress();

  TBranch *branch;
  TTree *treeD = gAlice->TreeD();


  if (treeD) {
	if (fDigits) {
	      branch = treeD->GetBranch(branchname);
	      if (branch) branch->SetAddress(&fDigits);
        }
      
  }
  if(fSDigits)
    fSDigits->Clear();

  if (gAlice->TreeS()  && fSDigits ) {
    branch = gAlice->TreeS()->GetBranch("FMD");
    if (branch) branch->SetAddress(&fSDigits) ;
  } 


}

//---------------------------------------------------------------------
/*
void AliFMD::SDigits2Digits() 
{
     Int_t  ibg=0, bgr=0;
    Int_t nbgr_ev= 0;
       
    if (ibg) {
         printf("nbgr_ev %d\n",nbgr_ev);
 //        Hit2Digits(nbgr_ev,"Add"," ","galice_bgr.root");
    } 
       else {
   //      Hit2Digits(nbgr_ev,"rien","",""); 
    }

}
*/
//---------------------------------------------------------------------



void AliFMD::Eta2Radius(Float_t eta, Float_t zDisk, Float_t *radius)
{
  Float_t expEta=TMath::Exp(-eta);
  Float_t theta=TMath::ATan(expEta);
  theta=2.*theta;
  Float_t rad=zDisk*(TMath::Tan(theta));
  *radius=rad;
   
  printf(" eta %f radius %f\n", eta, rad);
}

//---------------------------------------------------------------------

void AliFMD::Hits2SDigits(){
#ifdef DEBUG
  cout<<"ALiFMD::Hits2SDigits> start...\n";
#endif
  char * fileSDigits = 0 ;
  char * fileHeader=0 ;
  AliFMDSDigitizer * sd = new AliFMDSDigitizer("mgalice.root","FMD.SDigit.root") ;
  sd->Exec("") ;
  sd->Print("");
  
  delete sd ;
}



