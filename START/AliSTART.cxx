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
Revision 1.15  2001/01/01 13:10:42  hristov
Local definition of digits removed

Revision 1.14  2000/12/22 16:17:15  hristov
Updated  START code from Alla

Revision 1.13  2000/12/18 11:39:41  alibrary
Quick fix to avoid crash in display waiting for new version

Revision 1.12  2000/12/04 08:48:19  alibrary
Fixing problems in the HEAD

Revision 1.11  2000/10/13 13:14:08  hristov
Bug fixes and code cleaning

Revision 1.10  2000/10/02 21:28:13  fca
Removal of useless dependecies via forward declarations

Revision 1.9  2000/07/13 16:41:29  fca
New START corrected for coding conventions

Revision 1.8  2000/03/27 17:24:25  alla
Modifing geometry

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
#include "TParticle.h"
#include "AliRun.h"
#include "AliSTART.h"
#include "AliSTARTdigit.h"
#include "AliMC.h"
#include "AliSTARThit.h"
#include "AliSTARTvertex.h"

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
  gAlice->AddHitList(fHits);
  
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

void AliSTART::AddDigit(Int_t *tracks,Int_t *digits)
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
  new TRotMatrix("rot999","rot999",  90,0,90,90,180,0);

  new TTUBE("S_STR1","START  volume 1","void",5.,10.7,5.3);
  top->cd();
  node = new TNode("STR1","STR1","S_STR1",0,0,75.,"");
  node->SetLineColor(kColorSTART);
  fNodes->Add(node);

  new TTUBE("S_STR2","START volume 2","void",5.,10.7,5.3);
  top->cd();
  node = new TNode("STR2","STR2","S_STR2",0,0,-75,"rot999");
  node->SetLineColor(kColorSTART);
  fNodes->Add(node);
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
  //
  // Specific START branches
  //
  // Create Tree branches for the START.
  Int_t buffersize = 400;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option);

  TTree *td = gAlice->TreeD();
  digits = new AliSTARTdigit();
  td->Branch(branchname,"AliSTARTdigit",&digits, buffersize);
  printf("Making Branch %s for digits\n",branchname);
    
/*
  gAlice->TreeR()->Branch(branchname,"Int_t",&fZposit, buffersize);
  printf("Making Branch %s for vertex position %d\n",branchname);
  */
}    

//_____________________________________________________________________________

void AliSTART::Hit2digit(Int_t evnum) 
{
  //
  // From hits to digits
  //
  Float_t x,y,e;
  Int_t nbytes = 0;
  Int_t hit;
  Int_t nhits;
  Int_t volume,pmt;
  char nameTH[8],nameTD[8];
  Float_t timediff,timeright,timeleft,timeav;
  Float_t besttimeright,besttimeleft,meanTime;
  Int_t channelWidth=10;

  TParticle *particle;
  AliSTARThit  *startHit;

  Int_t buffersize=256;
  Int_t split=1;

  digits= new AliSTARTdigit();
  TBranch *bDig=0;

  
 // Event ------------------------- LOOP  
 
    sprintf(nameTD,"TreeD%d",evnum);
    TTree *td = new TTree(nameTD,"START");
    bDig = td->Branch("START","AliSTARTdigit",&digits,buffersize,split);

    besttimeright=9999.;
    besttimeleft=9999.;
    Int_t timeDiff=0;
    Int_t timeAv=0;

    Int_t nparticles = gAlice->GetEvent(evnum);
    if (nparticles <= 0) return;
    printf("\n nparticles %d\n",nparticles);
    
    TClonesArray *particles = gAlice->Particles();
   
    sprintf(nameTH,"TreeH%d",evnum);
    printf("%s\n",nameTH);
    TTree *th = gAlice->TreeH();
    Int_t ntracks    = (Int_t) th->GetEntries();
    if (ntracks<=0) return;
    // Start loop on tracks in the hits containers
    for (Int_t track=0; track<ntracks;track++) {
      gAlice->ResetHits();
      nbytes += th->GetEvent(track);
      particle=(TParticle*)particles->UncheckedAt(track);
      nhits =fHits->GetEntriesFast();
      
      for (hit=0;hit<nhits;hit++) {
	startHit   = (AliSTARThit*)fHits->UncheckedAt(hit);
	pmt=startHit->fPmt;
	e=startHit->fEtot;
	x=startHit->X();
	y=startHit->Y();
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
       Float_t ppBunch=25;
    
       ppBunch=ppBunch-10/2;
       Float_t t1=1000.*besttimeleftGaus;
       Float_t t2=1000.*besttimerightGaus;
       t1=t1/channelWidth+ppBunch; //time in ps to channelWidth
       t2=t2/channelWidth+ppBunch; //time in ps to channelWidth
     
       timeav=(t1+t2)/2.;
     
       // Time to TDC signal
       // 256 channels for timediff, range 1ns
       
       timediff=128+1000*timediff/channelWidth; // time in ps 

       timeAv = (Int_t)(timeav);   // time (ps) channel numbres
       timeDiff = (Int_t)(timediff); // time ( ps) channel numbres
       digits->Set(timeAv,timeDiff);
     }
    else
      {timeAv=999999; timeDiff=99999;}
      
    td->Fill();
    printf("digits-> %d \n",digits->GetTime());
    td->Write();

}












 

