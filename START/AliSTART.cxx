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
Revision 1.5  2000/01/19 17:17:15  fca
Introducing a list of lists of hits -- more hits allowed for detector now

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

#include <iostream.h>
#include <fstream.h>

#include "TMath.h"
#include "TTUBE.h"
#include "TNode.h"
#include "TGeometry.h"
#include "AliRun.h"
#include "AliSTART.h"
#include "AliSTARTdigit.h"
#include "AliMC.h"
#include "AliSTARThit.h"

ClassImp(AliSTART)
 
//_____________________________________________________________________________
AliSTART::AliSTART()
{
  //
  // Default constructor for class AliSTART
  //
  fIshunt   = 0;
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
  fHits       = new TClonesArray("AliSTARThit",  405);
  gAlice->AddHitList(fHits);
  fDigits     = new TClonesArray("AliSTARTdigit",500);
  
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
void AliSTART::AddDigit(Int_t *tracks,Int_t *digits)
{
  //
  // Add a START digit to the list
  //
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliSTARTdigit(tracks,digits);
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

//_____________________________________________________________________________
void AliSTART::Hit2digit(Int_t evnum) 
{
  
  Float_t x,y,e;
  Int_t nbytes = 0;
  Int_t hit,i;
  Int_t nhits;
  Int_t volume,pmt;
  char nameTH[8];
  Float_t timediff,timeright,timeleft,t1,t2,timeav;
  Float_t besttimeright,besttimeleft;
  Float_t pp_bunch=25;
  Int_t channel_width=10;
  Int_t digits[3];
  Int_t tracks[2];

  TParticle *particle;

  AliSTARThit  *startHit;


  // Event ------------------------- LOOP  
 //   for (evnum=0; evnum<=9; evnum++){

    besttimeright=9999.;
    besttimeleft=9999.;

    Int_t nparticles = gAlice->GetEvent(evnum);
    if (nparticles <= 0) return;
    printf("\n nparticles %d\n",nparticles);
    
    TClonesArray *Particles = gAlice->Particles();
    
    sprintf(nameTH,"TreeH%d",evnum);
    printf("%s\n",nameTH);
    TTree *TH = gAlice->TreeH();
    Int_t ntracks    = (Int_t) TH->GetEntries();
    if (ntracks<=0) return;
    // Start loop on tracks in the hits containers
    for (Int_t track=0; track<ntracks;track++) {
      gAlice->ResetHits();
      nbytes += TH->GetEvent(track);
      particle=(TParticle*)Particles->UncheckedAt(track);
      nhits = fHits->GetEntriesFast();
           
      for (hit=0;hit<nhits;hit++) {
	startHit   = (AliSTARThit*)fHits->UncheckedAt(hit);
	pmt=startHit->fPmt;
	e=startHit->fEtot;
	x=startHit->fX;
	y=startHit->fY;
	volume = startHit->fVolume;
	if(volume==1){
	  timeright = startHit->fTime;
	  if(timeright<besttimeright) {
	    besttimeright=timeright;
	    tracks[0]=track;
	  } //timeright
	}//time for right shoulder
	if(volume==2){            
	  timeleft = startHit->fTime;
	  //                printf("timeleft %f\n",timeleft);
	  if(timeleft<besttimeleft) {
	    besttimeleft=timeleft;
	    tracks[1]=track;
	  } //timeleftbest
	}//time for left shoulder
      } //hit loop
    } //track loop
    printf("\n----time1stright %f \n",besttimeright);     
    printf("----time1stleft %f \n",besttimeleft);     
    timediff=besttimeright-besttimeleft;
    if (timediff!=0 && TMath::Abs(timediff)<100) {
      //we assume centre of bunch is 5ns after TTS signal
      //TOF values are relative of the end of bunch
      pp_bunch=pp_bunch-10/2;
      t1=besttimeleft+pp_bunch;
      t2=besttimeright+pp_bunch;
      t1=1000*t1/channel_width; //time in ps to channel_width
      t2=1000*t2/channel_width; //time in ps to channel_width
      printf(" t1= %f t2= %f\n",t1,t2);

      timeav=(t1+t2)/2.;
      printf("timediff= %f timeav= %f\n",timediff,timeav);

      // Time to TDC signal
      // 1024 channels for timediff, range 1ns
      
     timediff=512+1000*timediff/channel_width; // time in ps
     printf("timediff= %f timeav= %f\n",timediff,timeav);


     digits[0]=evnum;
     digits[1]=(Int_t)(timeav);   // time in ps
     digits[2]=(Int_t)(timediff); // time in ps
     //  new(ldigits[fNdigits++]) AliSTARTdigit(track,digits);
 
    
     for (i=0; i<3; i++){
       printf(" DIGITS on START  %d\n",digits[i]); } 
     for (i=0; i<=1; i++) { printf("START track %d\n",tracks[i]);}
     AddDigit(tracks,digits);
     //     sprintf(nameTD,"TreeD%d",evnum);
     //    gAlice->TreeD()->Fill();
     //gAlice->TreeD()->Write();
     //printf("%s\n",nameTD);
     MakeTree(evnum);
     if (fTreeD!=0) fTreeD->Fill();    
     if (fTreeD!=0) fTreeD->Write();    
    } // if timediff !=0
    
    //   } // event loop
    
} // end of mcro
 
//_____________________________________________________________________________
Bool_t  AliSTART::SetTree(Int_t nevent, TDirectory *dir )
{
  char treeName[100];
  // Get Hits Tree header from file
  sprintf(treeName,"TreeD%d",nevent);
  fTreeD = (TTree*)dir->Get(treeName);
  if (fTreeD == 0) return kFALSE;
  //set Digit branch 
  TBranch *b = fTreeD->GetBranch("Digits");
  if (b==0) return kFALSE;
  b->SetAddress(&fDigits);
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t  AliSTART::MakeTree(Int_t nevent)
{
  char treeName[100];
  // Get Hits Tree header from file
  sprintf(treeName,"TreeD%d",nevent);
  fTreeD =  new TTree(treeName,treeName);
  if (fTreeD == 0) return kFALSE;
  //set Digit branch 
  TBranch *b = fTreeD->Branch("Digits",&fDigits,40000);
  if (b==0) return kFALSE;
  b->SetAddress(&fDigits);
 
  return kTRUE;
}
