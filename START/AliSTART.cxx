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
Revision 1.6  2000/01/21 15:45:23  fca
New Version from Alla

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
#include "TRandom.h"
#include "TGeometry.h"
#include "AliRun.h"
#include "AliSTART.h"
#include "AliSTARTdigit.h"
#include "AliMC.h"
#include "AliSTARThit.h"
#include "AliSTARTvertex.h"

ClassImp(AliSTART)
  AliSTARTdigit *digits; 
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

  
  //
  // Initialise Hit array
  fHits       = new TClonesArray("AliSTARThit",  405);
  //  gAlice->AddHitList(fHits);
  //  fDigits     = new TClonesArray("AliSTARTdigit",500);
  
  fIshunt     =  0;
  fIdSens   =  0;
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
  
  //  Add a START digit to the list
  
//  printf (" AddDigit*******");
    // TClonesArray &ldigits = *fDigits;
    // new(ldigits[fNdigits++]) AliSTARTdigit(tracks,digits);
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
  //  fIdSensRad=gMC->VolId("PTOP");
  //  fIdSensPC =gMC->VolId("T0PC");

}

//---------------------------------------------------------------------------
void AliSTART::MakeBranch(Option_t* option)
{
  
  AliSTARTdigit *digits; 
  // Create Tree branches for the START.
  Int_t buffersize = 400;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option);

  TTree *TD = gAlice->TreeD();
  digits = new AliSTARTdigit();
  TD->Branch(branchname,"AliSTARTdigit",&digits, buffersize);
  printf("Making Branch %s for digits\n",branchname);
    
/*
  gAlice->TreeR()->Branch(branchname,"Int_t",&fZposit, buffersize);
  printf("Making Branch %s for vertex position %d\n",branchname);
  */
}    

//_____________________________________________________________________________

void AliSTART::Hit2digit(Int_t evnum) 
{
  
  Float_t x,y,z,e;
  Int_t nbytes = 0;
  Int_t j,hit;
  Int_t nhits;
  Int_t volume,pmt;
  char nameTH[8],nameTD[8];
  Float_t timediff,timeright,timeleft,timeav;
  Float_t besttimeright,besttimeleft,meanTime;
  Int_t channel_width=10;

  TParticle *particle;
  AliSTARThit  *startHit;

  Int_t buffersize=256;
  Int_t split=1;

  digits= new AliSTARTdigit();
  TBranch *bDig=0;

  /*    
  // Create histograms
  
   TH1F *hTimediff = new TH1F("hTimediff","Time different",100,-2,2);
   TH1F *hMeanTime = new TH1F("hMeanTime","Mean Time",100,2.2,2.8);
  
   TH1F *hTime1stright = new TH1F("hTime1stright","Time flight of 1st  particle right", 100,1.5,3.2);
   TH1F *hTime1stleft = new  TH1F("hTime1sleft","Time flight of 1st particle left",100,1.5,3.2);
  
  */ 
   //   AliSTART *START  = (AliSTART*) gAlice->GetDetector("START");
  
 // Event ------------------------- LOOP  
 
    sprintf(nameTD,"TreeD%d",evnum);
    TTree *TD = new TTree(nameTD,"START");
    bDig = TD->Branch("START","AliSTARTdigit",&digits,buffersize,split);

    besttimeright=9999.;
    besttimeleft=9999.;
    Int_t Timediff=0;
    Int_t Timeav=0;

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
      nhits =fHits->GetEntriesFast();
      
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
	  } //timeright
	}//time for right shoulder
	if(volume==2){            
	  timeleft = startHit->fTime;
	  //                printf("timeleft %f\n",timeleft);
	  if(timeleft<besttimeleft) {
	    besttimeleft=timeleft;
	  } //timeleftbest
	}//time for left shoulder
      } //hit loop
    } //track loop

    //folding with experimental time distribution
   Float_t besttimerightGaus=gRandom->Gaus(besttimeright,0.05);
   Float_t besttimeleftGaus=gRandom->Gaus(besttimeleft,0.05);
   timediff=besttimerightGaus-besttimeleftGaus;
   meanTime=(besttimerightGaus+besttimeleftGaus)/2.;
  if ( TMath::Abs(timediff)<2. && meanTime<3.) 
     {
     //we assume centre of bunch is 5ns after TTS signal
     //TOF values are relative of the end of bunch
       //      hTimediff->Fill(timediff);
       //hMeanTime->Fill(meanTime);
       Float_t pp_bunch=25;
    
       pp_bunch=pp_bunch-10/2;
       Float_t t1=1000.*besttimeleftGaus;
       Float_t t2=1000.*besttimerightGaus;
       t1=t1/channel_width+pp_bunch; //time in ps to channel_width
       t2=t2/channel_width+pp_bunch; //time in ps to channel_width
     
       timeav=(t1+t2)/2.;
     
       // Time to TDC signal
       // 256 channels for timediff, range 1ns
       
       timediff=128+1000*timediff/channel_width; // time in ps
 

       Timeav = (Int_t)(timeav);   // time in ps
       Timediff = (Int_t)(timediff); // time in ps
       digits->Set(Timeav,Timediff);
       TD->Fill();
       digits->MyDump();
       TD->Write();
     } //timediff
   

} // end macro











 

