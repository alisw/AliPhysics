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
Revision 1.14.4.17  2000/06/14 14:36:46  morsch
- add TriggerCircuit (PC)
- add GlobalTrigger and LocalTrigger and specific methods (PC)

Revision 1.14.4.16  2000/06/09 21:20:28  morsch
Most coding rule violations corrected

Revision 1.14.4.15  2000/05/02 09:54:32  morsch
RULE RN17 violations corrected

Revision 1.14.4.12  2000/04/26 12:25:02  morsch
Code revised by P. Crochet:
- Z position of TriggerChamber changed according to A.Tournaire Priv.Comm.
- ToF included in the method MakePadHits
- inner radius of flange between beam shielding and trigger corrected
- Trigger global volume updated (according to the new geometry)

Revision 1.14.4.11  2000/04/19 19:42:08  morsch
Some changes of variable names curing viols and methods concerning
correlated clusters removed.

Revision 1.14.4.10  2000/03/22 16:44:07  gosset
Memory leak suppressed in function Digitise:
p_adr->Delete() instead of Clear (I.Chevrot and A.Baldisseri)

Revision 1.14.4.9  2000/03/20 18:15:25  morsch
Positions of trigger chambers corrected (P.C.)

Revision 1.14.4.8  2000/02/21 15:38:01  morsch
Call to AddHitList introduced to make this version compatible with head.

Revision 1.14.4.7  2000/02/20 07:45:53  morsch
Bugs in Trigger part of BuildGeomemetry corrected (P.C)

Revision 1.14.4.6  2000/02/17 14:28:54  morsch
Trigger included into initialization and digitization

Revision 1.14.4.5  2000/02/15 10:02:58  morsch
Log messages of previous revisions added

Revision 1.14.4.2  2000/02/04 10:57:34  gosset
Z position of the chambers:
it was the Z position of the stations;
it is now really the Z position of the chambers.
   !!!! WARNING: THE CALLS TO "AliMUONChamber::SetZPOS"
   !!!!                   AND "AliMUONChamber::ZPosition"
   !!!! HAVE TO BE CHANGED TO "AliMUONChamber::"SetZ"
   !!!!                   AND "AliMUONChamber::Z"

Revision 1.14.4.3  2000/02/04 16:19:04  gosset
Correction for mis-spelling of NCH 

Revision 1.14.4.4  2000/02/15 09:43:38  morsch
Log message added

*/


///////////////////////////////////////////////
//  Manager and hits classes for set:MUON     //
////////////////////////////////////////////////

#include <TTUBE.h>
#include <TBRIK.h>
#include <TRotMatrix.h>
#include <TNode.h> 
#include <TTree.h> 
#include <TRandom.h> 
#include <TObject.h>
#include <TVector.h>
#include <TObjArray.h>
#include <TMinuit.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TDirectory.h>
#include <TObjectTable.h>
#include <AliPDG.h>
#include <TTUBE.h>

#include "AliMUON.h"
#include "AliMUONHit.h"
#include "AliMUONPadHit.h"
#include "AliMUONDigit.h"
#include "AliMUONTransientDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONHitMap.h"
#include "AliMUONHitMapA1.h"
#include "AliMUONChamberTrigger.h"

#include "AliMUONClusterFinder.h"
#include "AliMUONTriggerDecision.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"
#include "AliCallf77.h" 
#include "AliConst.h" 

// Defaults parameters for Z positions of chambers
// taken from values for "stations" in AliMUON::AliMUON
//     const Float_t zch[7]={528, 690., 975., 1249., 1449., 1610, 1710.};
// and from array "dstation" in AliMUONv1::CreateGeometry
//          Float_t dstation[5]={20., 20., 20, 20., 20.};
//     for tracking chambers,
//          according to (Z1 = zch - dstation) and  (Z2 = zch + dstation)
//          for the first and second chambers in the station, respectively,
// and from "DTPLANES" in AliMUONv1::CreateGeometry
//           const Float_t DTPLANES = 15.;
//     for trigger chambers,
//          according to (Z1 = zch) and  (Z2 = zch + DTPLANES)
//          for the first and second chambers in the station, respectively
static const Float_t kDefaultChambersZ[kNCH] =
{518., 538., 680., 700., 965., 985., 1239., 1259., 1439., 1459.,
 1603.5, 1618.5, 1703.5, 1718.5}; 

ClassImp(AliMUON)
//___________________________________________
AliMUON::AliMUON()
{
   fIshunt       = 0;
   fHits         = 0;
   fPadHits      = 0;
   fNPadHits     = 0;
   fDchambers    = 0;
   fTriggerCircuits = 0;     // cp new design of AliMUONTriggerDecision
   fNdch         = 0;
   fRawClusters  = 0;
   fNrawch       = 0;
   fGlobalTrigger   = 0;
   fNLocalTrigger   = 0;
   fLocalTrigger    = 0;
   fNLocalTrigger   = 0;
}
 
//___________________________________________
AliMUON::AliMUON(const char *name, const char *title)
       : AliDetector(name,title)
{
//Begin_Html
/*
<img src="gif/alimuon.gif">
*/
//End_Html
 
   fHits     = new TClonesArray("AliMUONHit",1000);
   gAlice->AddHitList(fHits);
   fPadHits = new TClonesArray("AliMUONPadHit",10000);
   fNPadHits  =  0;
   fIshunt     =  0;

   fNdch      = new Int_t[kNCH];

   fDchambers = new TObjArray(kNCH);

   Int_t i;
   
   for (i=0; i<kNCH ;i++) {
       (*fDchambers)[i] = new TClonesArray("AliMUONDigit",10000); 
       fNdch[i]=0;
   }

   fNrawch      = new Int_t[kNTrackingCh];

   fRawClusters = new TObjArray(kNTrackingCh);

   for (i=0; i<kNTrackingCh;i++) {
       (*fRawClusters)[i] = new TClonesArray("AliMUONRawCluster",10000); 
       fNrawch[i]=0;
   }
   cout << " here " << "\n";

   fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1);    
   fNGlobalTrigger = 0;
   fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);    
   fNLocalTrigger = 0;

//   
// Transport angular cut
   fAccCut=0;
   fAccMin=2;
   fAccMax=9;

   SetMarkerColor(kRed);
//
//
//
//
//  inner diameter
   const Float_t kDmin[7]={ 35.,  47.,  66.,   80.,  80., 100., 100.};
//
//  outer diameter
   const Float_t kDmax[7]={183., 245., 316.6,  520.,  520., 830., 880.};
//
    Int_t ch;

    fChambers = new TObjArray(kNCH);

    // Loop over stations
    for (Int_t st = 0; st < kNCH / 2; st++) {
      // Loop over 2 chambers in the station
	for (Int_t stCH = 0; stCH < 2; stCH++) {
//
//    
//    Default Parameters for Muon Tracking Stations


	    ch = 2 * st + stCH;
//
	    if (ch < kNTrackingCh) {
		(*fChambers)[ch] = new AliMUONChamber();
	    } else {
		(*fChambers)[ch] = new AliMUONChamberTrigger();
	    }
	    
	    AliMUONChamber* chamber = (AliMUONChamber*) (*fChambers)[ch];
	    
	    chamber->SetGid(0);
	    // Default values for Z of chambers
	    chamber->SetZ(kDefaultChambersZ[ch]);
//
	    chamber->InitGeo(kDefaultChambersZ[ch]);
	    chamber->SetRInner(kDmin[st]/2);
	    chamber->SetROuter(kDmax[st]/2);
//
	} // Chamber stCH (0, 1) in 
    }     // Station st (0...)
    fMaxStepGas=0.01; 
    fMaxStepAlu=0.1; 
    fMaxDestepGas=-1;
    fMaxDestepAlu=-1;
//
   fMaxIterPad   = 0;
   fCurIterPad   = 0;

   // cp new design of AliMUONTriggerDecision
   fTriggerCircuits = new TObjArray(kNTriggerCircuit);
   for (Int_t circ=0; circ<kNTriggerCircuit; circ++) {
     (*fTriggerCircuits)[circ] = new AliMUONTriggerCircuit();     
   }
   // cp new design of AliMUONTriggerDecision

}
 
//___________________________________________
AliMUON::AliMUON(const AliMUON& rMUON)
{
// Dummy copy constructor
    ;
    
}

AliMUON::~AliMUON()
{

    printf("Calling AliMUON destructor !!!\n");
    
  Int_t i;
  fIshunt  = 0;
  delete fHits;
  delete fPadHits;

  delete fGlobalTrigger;
  fNGlobalTrigger = 0;

  delete fLocalTrigger;
  fNLocalTrigger = 0;

  for (i=0;i<kNCH;i++) {
      delete (*fDchambers)[i];
      fNdch[i]=0;
  }
  delete fDchambers;

  for (i=0;i<kNTrackingCh;i++) {
      delete (*fRawClusters)[i];
      fNrawch[i]=0;
  }
  delete fRawClusters;

  for (Int_t circ=0; circ<kNTriggerCircuit; circ++) {
    delete (*fTriggerCircuits)[circ];
  }
  delete fTriggerCircuits;
}
 
//___________________________________________
void AliMUON::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliMUONHit(fIshunt,track,vol,hits);
}
//___________________________________________
void AliMUON::AddPadHit(Int_t *clhits)
{
   TClonesArray &lclusters = *fPadHits;
   new(lclusters[fNPadHits++]) AliMUONPadHit(clhits);
}
//_____________________________________________________________________________
void AliMUON::AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits)
{
    //
    // Add a MUON digit to the list
    //

    TClonesArray &ldigits = *((TClonesArray*)(*fDchambers)[id]);
    new(ldigits[fNdch[id]++]) AliMUONDigit(tracks,charges,digits);
}

//_____________________________________________________________________________
void AliMUON::AddRawCluster(Int_t id, const AliMUONRawCluster& c)
{
    //
    // Add a MUON digit to the list
    //

    TClonesArray &lrawcl = *((TClonesArray*)(*fRawClusters)[id]);
    new(lrawcl[fNrawch[id]++]) AliMUONRawCluster(c);
}

//___________________________________________
void AliMUON::AddGlobalTrigger(Int_t *singlePlus, Int_t *singleMinus,
			       Int_t *singleUndef,
			       Int_t *pairUnlike, Int_t *pairLike)
{
// add a MUON Global Trigger to the list (only one GlobalTrigger per event !)
  TClonesArray &globalTrigger = *fGlobalTrigger;
  new(globalTrigger[fNGlobalTrigger++]) 
    AliMUONGlobalTrigger(singlePlus, singleMinus,  singleUndef, pairUnlike, 
			 pairLike);
}
//___________________________________________
void AliMUON::AddLocalTrigger(Int_t *localtr)
{
// add a MUON Local Trigger to the list
  TClonesArray &localTrigger = *fLocalTrigger;
  new(localTrigger[fNLocalTrigger++]) AliMUONLocalTrigger(localtr);
}

//___________________________________________
void AliMUON::BuildGeometry()
{
    TNode *node, *nodeF, *top, *nodeS;
    const int kColorMUON  = kBlue;
    const int kColorMUON2 = kGreen; 
    //
    top=gAlice->GetGeometry()->GetNode("alice");
// MUON
//
//  z-Positions of Chambers
    const Float_t kCz[7]={511., 686., 971., 1245., 1445., 1600, 1700.};
//  inner diameter (Xlenght for trigger chamber -> active area)
    const Float_t kDmin[7]={ 35.,  47.,  67.,   86.,  100., 544., 544.};
//  outer diameter (Ylenght for trigger chamber -> active area)
    const Float_t kDmax[7]={183., 245., 346.,  520.,  520., 612., 612.};

    TRotMatrix* rot000 = new TRotMatrix("Rot000"," ", 90,  0, 90, 90, 0, 0);
    TRotMatrix* rot090 = new TRotMatrix("Rot090"," ", 90, 90, 90,180, 0, 0);
    TRotMatrix* rot180 = new TRotMatrix("Rot180"," ", 90,180, 90,270, 0, 0);
    TRotMatrix* rot270 = new TRotMatrix("Rot270"," ", 90,270, 90,  0, 0, 0);
    
    Float_t rmin, rmax, dx, dy, dz, dr, xpos, ypos, zpos;
    Float_t dzc1=4.;           // tracking chambers
    Float_t dzc2=15.;          // trigger chambers
    Float_t hole=102.;          // x-y hole around beam pipe for trig. chambers
    Float_t zscale;            // scaling parameter trigger chambers
    Float_t halfx, halfy;   
    char nameChamber[9], nameSense[9], nameFrame[9], nameNode[8];
    char nameSense1[9], nameSense2[9];    
    for (Int_t i=0; i<7; i++) {
	for (Int_t j=0; j<2; j++) {
	    Int_t id=2*i+j+1;
	    if (i<5) {               // tracking chambers
		if (j==0) {
		    zpos=kCz[i]-dzc1;
		} else {
		    zpos=kCz[i]+dzc1;
		}
	    } else {
		if (j==0) {
		    zpos=kCz[i];
		} else {            
		    zpos=kCz[i]+dzc2;
		}
	    }
	    sprintf(nameChamber,"C_MUON%d",id);
	    sprintf(nameSense,"S_MUON%d",id);
	    sprintf(nameSense1,"S1_MUON%d",id);
	    sprintf(nameSense2,"S2_MUON%d",id);
	    sprintf(nameFrame,"F_MUON%d",id);	
	    if (i<5) {	                      // tracking chambers
		rmin = kDmin[i]/2.-3;
		rmax = kDmax[i]/2.+3;
		new TTUBE(nameChamber,"Mother","void",rmin,rmax,0.25,1.);
		rmin = kDmin[i]/2.;
		rmax = kDmax[i]/2.;
		new TTUBE(nameSense,"Sens. region","void",rmin,rmax,0.25, 1.);
		dx=(rmax-rmin)/2;
		dy=3.;
		dz=0.25;
		TBRIK* frMUON = new TBRIK(nameFrame,"Frame","void",dx,dy,dz);
		top->cd();
		sprintf(nameNode,"MUON%d",100+id);
		node = new TNode(nameNode,"ChamberNode",nameChamber,0,0,zpos,"");
		node->SetLineColor(kColorMUON);
		fNodes->Add(node);
		node->cd();
		sprintf(nameNode,"MUON%d",200+id);
		node = new TNode(nameNode,"Sens. Region Node",nameSense,0,0,0,"");
		node->SetLineColor(kColorMUON);
		node->cd();
		dr=dx+rmin;
		sprintf(nameNode,"MUON%d",300+id);
		nodeF = new TNode(nameNode,"Frame0",frMUON,dr, 0, 0,rot000,"");
		nodeF->SetLineColor(kColorMUON);
		node->cd();
		sprintf(nameNode,"MUON%d",400+id);
		nodeF = new TNode(nameNode,"Frame1",frMUON,0 ,dr,0,rot090,"");
		nodeF->SetLineColor(kColorMUON);
		node->cd();
		sprintf(nameNode,"MUON%d",500+id);
		nodeF = new TNode(nameNode,"Frame2",frMUON,-dr,0,0,rot180,"");
		nodeF->SetLineColor(kColorMUON);
		node  ->cd();
		sprintf(nameNode,"MUON%d",600+id);   
		nodeF = new TNode(nameNode,"Frame3",frMUON,0,-dr,0,rot270,"");
		nodeF->SetLineColor(kColorMUON);
	    } else { 
		zscale=zpos/kCz[5];
		Float_t xsize=kDmin[i]*zscale;
		Float_t ysize=kDmax[i]*zscale;
		Float_t holeScaled=hole*zscale;
		
		halfx=xsize/2.+3.;
		halfy=ysize/2.+3.;	    
		new TBRIK(nameChamber,"Mother","void",halfx,halfy,0.25);
		top->cd();
		sprintf(nameNode,"MUON%d",100+id);
		node = new TNode(nameNode,"Chambernode",nameChamber,0,0,zpos,"");
		node->SetLineColor(kColorMUON2);
		fNodes->Add(node);
		
// up/down of beam pipe
		halfx=xsize/2.;
		halfy=(ysize/2.-holeScaled/2.)/2.;	    
		new TBRIK(nameSense,"Sens. region","void",halfx,halfy,0.25);
		
		node->cd();
		ypos=holeScaled/2.+((ysize/2.-holeScaled/2.)/2.);
		sprintf(nameNode,"MUON%d",200+id);
		nodeS = new TNode(nameNode,"Sens. Region Node",nameSense,0,ypos,0,"");
		nodeS->SetLineColor(kColorMUON2);
		
		node->cd();
		ypos=-1.*ypos;
		sprintf(nameNode,"MUON%d",300+id);
		nodeS = new TNode(nameNode,"Sens. Region Node",nameSense,0,ypos,0,"");
		nodeS->SetLineColor(kColorMUON2);
		
// left/right of beam pipe
		halfx=(xsize/2.-holeScaled/2.)/2.;
		halfy=holeScaled/2.;	
		new TBRIK(nameSense1,"Sens. region","void",halfx,halfy,0.25);
		
		node->cd();
		xpos=holeScaled/2.+((xsize/2.-holeScaled/2.)/2.);	    
		sprintf(nameNode,"MUON%d",400+id);
		nodeS = new TNode(nameNode,"Sens. Region Node",nameSense1,xpos,0,0,"");
		nodeS->SetLineColor(kColorMUON2);
		
		node->cd();
		xpos=-1.*xpos;
		sprintf(nameNode,"MUON%d",500+id);
		nodeS = new TNode(nameNode,"Sens. Region Node",nameSense1,xpos,0,0,"");
		nodeS->SetLineColor(kColorMUON2);
		
// missing corners
		halfx=17.*zscale/2.;
		halfy=halfx;
		new TBRIK(nameSense2,"Sens. region","void",halfx,halfy,0.25);
		
		node->cd();
		xpos=holeScaled/2.-halfx;
		ypos=xpos;
		sprintf(nameNode,"MUON%d",600+id);
		nodeS = new TNode(nameNode,"Sens. Region Node",nameSense2,xpos,ypos,0,"");
		nodeS->SetLineColor(kColorMUON2);
		
		node->cd();
		ypos=-1.*xpos;
		sprintf(nameNode,"MUON%d",700+id);
		nodeS = new TNode(nameNode,"Sens. Region Node",nameSense2,xpos,ypos,0,"");
		nodeS->SetLineColor(kColorMUON2);
		
		node->cd();
		xpos=-1.*xpos;
		sprintf(nameNode,"MUON%d",800+id);
		nodeS = new TNode(nameNode,"Sens. Region Node",nameSense2,xpos,ypos,0,"");
		nodeS->SetLineColor(kColorMUON2);
		
		node->cd();
		ypos=-1.*xpos;
		sprintf(nameNode,"MUON%d",900+id);
		nodeS = new TNode(nameNode,"Sens. Region Node",nameSense2,xpos,ypos,0,"");
		nodeS->SetLineColor(kColorMUON2);
	    } 
	}
    }
}


//___________________________________________
Int_t AliMUON::DistancetoPrimitive(Int_t , Int_t )
{
   return 9999;
}

//___________________________________________
void AliMUON::MakeBranch(Option_t* option)
{
  // Create Tree branches for the MUON.
  
  const Int_t kBufferSize = 4000;
  char branchname[30];
  sprintf(branchname,"%sCluster",GetName());

  AliDetector::MakeBranch(option);

  if (fPadHits   && gAlice->TreeH()) {
    gAlice->TreeH()->Branch(branchname,&fPadHits, kBufferSize);
    printf("Making Branch %s for clusters\n",branchname);
  }

// one branch for digits per chamber
  Int_t i;
  
  for (i=0; i<kNCH ;i++) {
      sprintf(branchname,"%sDigits%d",GetName(),i+1);
      
      if (fDchambers   && gAlice->TreeD()) {
	  gAlice->TreeD()->Branch(branchname,&((*fDchambers)[i]), kBufferSize);
	  printf("Making Branch %s for digits in chamber %d\n",branchname,i+1);
      }	
  }

  printf("Make Branch - TreeR address %p\n",gAlice->TreeR());

// one branch for raw clusters per chamber
  for (i=0; i<kNTrackingCh ;i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
      
      if (fRawClusters   && gAlice->TreeR()) {
	 gAlice->TreeR()->Branch(branchname,&((*fRawClusters)[i]), kBufferSize);
	 printf("Making Branch %s for raw clusters in chamber %d\n",branchname,i+1);
      }	
  }

// one branch for global trigger
  sprintf(branchname,"%sGlobalTrigger",GetName());
  if (fGlobalTrigger && gAlice->TreeR()) {  
    gAlice->TreeR()->Branch(branchname,&fGlobalTrigger,kBufferSize);
    printf("Making Branch %s for Global Trigger\n",branchname);
  }
// one branch for local trigger
  sprintf(branchname,"%sLocalTrigger",GetName());
  if (fLocalTrigger && gAlice->TreeR()) {  
    gAlice->TreeR()->Branch(branchname,&fLocalTrigger,kBufferSize);
    printf("Making Branch %s for Local Trigger\n",branchname);
  }

}

//___________________________________________
void AliMUON::SetTreeAddress()
{
  // Set branch address for the Hits and Digits Tree.
  char branchname[30];
  AliDetector::SetTreeAddress();

  TBranch *branch;
  TTree *treeH = gAlice->TreeH();
  TTree *treeD = gAlice->TreeD();
  TTree *treeR = gAlice->TreeR();

  if (treeH) {
    if (fPadHits) {
      branch = treeH->GetBranch("MUONCluster");
      if (branch) branch->SetAddress(&fPadHits);
    }
  }

  if (treeD) {
      for (int i=0; i<kNCH; i++) {
	  sprintf(branchname,"%sDigits%d",GetName(),i+1);
	  if (fDchambers) {
	      branch = treeD->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fDchambers)[i]));
	  }
      }
  }

  // printf("SetTreeAddress --- treeR address  %p \n",treeR);

  if (treeR) {
      for (int i=0; i<kNTrackingCh; i++) {
	  sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
	  if (fRawClusters) {
	      branch = treeR->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fRawClusters)[i]));
	  }
      }

      if (fLocalTrigger) {
	branch = treeR->GetBranch("MUONLocalTrigger");
	if (branch) branch->SetAddress(&fLocalTrigger);
      }
      if (fGlobalTrigger) {
	branch = treeR->GetBranch("MUONGlobalTrigger");
	if (branch) branch->SetAddress(&fGlobalTrigger);
      }
  }
}
//___________________________________________
void AliMUON::ResetHits()
{
  // Reset number of clusters and the cluster array for this detector
  AliDetector::ResetHits();
  fNPadHits = 0;
  if (fPadHits) fPadHits->Clear();
}

//____________________________________________
void AliMUON::ResetDigits()
{
    //
    // Reset number of digits and the digits array for this detector
    //
    for ( int i=0;i<kNCH;i++ ) {
	if ((*fDchambers)[i])    ((TClonesArray*)(*fDchambers)[i])->Clear();
	if (fNdch)  fNdch[i]=0;
    }
}
//____________________________________________
void AliMUON::ResetRawClusters()
{
    //
    // Reset number of raw clusters and the raw clust array for this detector
    //
    for ( int i=0;i<kNTrackingCh;i++ ) {
	if ((*fRawClusters)[i])    ((TClonesArray*)(*fRawClusters)[i])->Clear();
	if (fNrawch)  fNrawch[i]=0;
    }
}

//____________________________________________
void AliMUON::ResetTrigger()
{
  //  Reset Local and Global Trigger 
  fNGlobalTrigger = 0;
  if (fGlobalTrigger) fGlobalTrigger->Clear();
  fNLocalTrigger = 0;
  if (fLocalTrigger) fLocalTrigger->Clear();
}

//____________________________________________
void AliMUON::SetPadSize(Int_t id, Int_t isec, Float_t p1, Float_t p2)
{
    Int_t i=2*(id-1);
    ((AliMUONChamber*) (*fChambers)[i])  ->SetPadSize(isec,p1,p2);
    ((AliMUONChamber*) (*fChambers)[i+1])->SetPadSize(isec,p1,p2);
}

//___________________________________________
void AliMUON::SetChambersZ(const Float_t *Z)
{
  // Set Z values for all chambers (tracking and trigger)
  // from the array pointed to by "Z"
  for (Int_t ch = 0; ch < kNCH; ch++)
    ((AliMUONChamber*) ((*fChambers)[ch]))->SetZ(Z[ch]);
  return;
}

//___________________________________________
void AliMUON::SetChambersZToDefault()
{
  // Set Z values for all chambers (tracking and trigger)
  // to default values
  SetChambersZ(kDefaultChambersZ);
  return;
}

//___________________________________________
void AliMUON::SetChargeSlope(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONChamber*) (*fChambers)[i])->SetChargeSlope(p1);
    ((AliMUONChamber*) (*fChambers)[i+1])->SetChargeSlope(p1);
}

//___________________________________________
void AliMUON::SetChargeSpread(Int_t id, Float_t p1, Float_t p2)
{
    Int_t i=2*(id-1);
    ((AliMUONChamber*) (*fChambers)[i])->SetChargeSpread(p1,p2);
    ((AliMUONChamber*) (*fChambers)[i+1])->SetChargeSpread(p1,p2);
}

//___________________________________________
void AliMUON::SetSigmaIntegration(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONChamber*) (*fChambers)[i])->SetSigmaIntegration(p1);
    ((AliMUONChamber*) (*fChambers)[i+1])->SetSigmaIntegration(p1);
}

//___________________________________________
void AliMUON::SetMaxAdc(Int_t id, Float_t p1)
{
    Int_t i=2*(id-1);
    ((AliMUONChamber*) (*fChambers)[i])->SetMaxAdc(p1);
    ((AliMUONChamber*) (*fChambers)[i+1])->SetMaxAdc(p1);
}

//___________________________________________
void AliMUON::SetMaxStepGas(Float_t p1)
{
     fMaxStepGas=p1;
}

//___________________________________________
void AliMUON::SetMaxStepAlu(Float_t p1)
{
    fMaxStepAlu=p1;
}

//___________________________________________
void AliMUON::SetMaxDestepGas(Float_t p1)
{
    fMaxDestepGas=p1;
}

//___________________________________________
void AliMUON::SetMaxDestepAlu(Float_t p1)
{
    fMaxDestepAlu=p1;
}
//___________________________________________
void AliMUON::SetMuonAcc(Bool_t acc, Float_t angmin, Float_t angmax)
{
   fAccCut=acc;
   fAccMin=angmin;
   fAccMax=angmax;
}
//___________________________________________
void   AliMUON::SetSegmentationModel(Int_t id, Int_t isec, AliMUONSegmentation *segmentation)
{
    ((AliMUONChamber*) (*fChambers)[id])->SetSegmentationModel(isec, segmentation);

}
//___________________________________________
void   AliMUON::SetResponseModel(Int_t id, AliMUONResponse *response)
{
    ((AliMUONChamber*) (*fChambers)[id])->SetResponseModel(response);
}

void   AliMUON::SetReconstructionModel(Int_t id, AliMUONClusterFinder *reconst)
{
    ((AliMUONChamber*) (*fChambers)[id])->SetReconstructionModel(reconst);
}

void   AliMUON::SetNsec(Int_t id, Int_t nsec)
{
    ((AliMUONChamber*) (*fChambers)[id])->SetNsec(nsec);
}


//___________________________________________



void AliMUON::MakePadHits(Float_t xhit,Float_t yhit,
			  Float_t eloss, Float_t tof,  Int_t idvol)
{
//
//  Calls the charge disintegration method of the current chamber and adds
//  the simulated cluster to the root treee 
//
    Int_t clhits[7];
    Float_t newclust[6][500];
    Int_t nnew;
    
    
//
//  Integrated pulse height on chamber

    
    clhits[0]=fNhits+1;
//
//
    ((AliMUONChamber*) (*fChambers)[idvol])
	->DisIntegration(eloss, tof, xhit, yhit, nnew, newclust);
    Int_t ic=0;
    
//
//  Add new clusters
    for (Int_t i=0; i<nnew; i++) {
	if (Int_t(newclust[3][i]) > 0) {
	    ic++;
// Cathode plane
	    clhits[1] = Int_t(newclust[5][i]);
//  Cluster Charge
	    clhits[2] = Int_t(newclust[0][i]);
//  Pad: ix
	    clhits[3] = Int_t(newclust[1][i]);
//  Pad: iy 
	    clhits[4] = Int_t(newclust[2][i]);
//  Pad: charge
	    clhits[5] = Int_t(newclust[3][i]);
//  Pad: chamber sector
	    clhits[6] = Int_t(newclust[4][i]);
	    
	    AddPadHit(clhits);
	}
    }
}

//----------------------------------------------------------------------

void AliMUON::Digitise(Int_t nev,Int_t bgrEvent,Option_t *option,Option_t *opt,Text_t *filename)
{
    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
  
    static Bool_t first=kTRUE;
    static TFile *file;
    char *addBackground = strstr(option,"Add");

    AliMUONChamber*  iChamber;
    AliMUONSegmentation*  segmentation;

    
    Int_t trk[50];
    Int_t chtrk[50];  
    TObjArray *list=new TObjArray;
    static TClonesArray *pAddress=0;
    if(!pAddress) pAddress=new TClonesArray("TVector",1000);
    Int_t digits[5]; 

    AliMUON *pMUON  = (AliMUON *) gAlice->GetModule("MUON");
    AliMUONHitMap * hitMap[kNCH];
    for (Int_t i=0; i<kNCH; i++) {hitMap[i]=0;}
    if (addBackground ) {
	if(first) {
	    fFileName=filename;
	    cout<<"filename"<<fFileName<<endl;
	    file=new TFile(fFileName);
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
	    fHits2     = new TClonesArray("AliMUONHit",1000  );
	    fPadHits2 = new TClonesArray("AliMUONPadHit",10000);
	}	    
	first=kFALSE;
	file->cd();
	//file->ls();
	// Get Hits Tree header from file
	if(fHits2) fHits2->Clear();
	if(fPadHits2) fPadHits2->Clear();
	if(fTrH1) delete fTrH1;
	fTrH1=0;
	
	char treeName[20];
	sprintf(treeName,"TreeH%d",bgrEvent);
	fTrH1 = (TTree*)gDirectory->Get(treeName);
        //printf("fTrH1 %p of treename %s for event %d \n",fTrH1,treeName,bgrEvent);
	
	if (!fTrH1) {
	    printf("ERROR: cannot find Hits Tree for event:%d\n",bgrEvent);
	}
	// Set branch addresses
	TBranch *branch;
	char branchname[20];
	sprintf(branchname,"%s",GetName());
	if (fTrH1 && fHits2) {
	    branch = fTrH1->GetBranch(branchname);
	    if (branch) branch->SetAddress(&fHits2);
	}
	if (fTrH1 && fPadHits2) {
	    branch = fTrH1->GetBranch("MUONCluster");
	    if (branch) branch->SetAddress(&fPadHits2);
	}
// test
	//Int_t ntracks1 =(Int_t)fTrH1->GetEntries();
	//printf("background - ntracks1 - %d\n",ntracks1);
    }
    //
    // loop over cathodes
    //
    AliMUONHitMap* hm;
    Int_t countadr=0;
    for (int icat=0; icat<2; icat++) { 
	Int_t counter=0;
	for (Int_t i =0; i<kNCH; i++) {
	    iChamber=(AliMUONChamber*) (*fChambers)[i];
	    if (iChamber->Nsec()==1 && icat==1) {
		continue;
	    } else {
		segmentation=iChamber->SegmentationModel(icat+1);
	    }
	    hitMap[i] = new AliMUONHitMapA1(segmentation, list);
	}
	//printf("Start loop over tracks \n");     
//
//   Loop over tracks
//

	TTree *treeH = gAlice->TreeH();
	Int_t ntracks =(Int_t) treeH->GetEntries();
	Int_t nmuon[kNCH]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Float_t xhit[kNCH][2];
	Float_t yhit[kNCH][2];
	
	for (Int_t track=0; track<ntracks; track++) {
	    gAlice->ResetHits();
	    treeH->GetEvent(track);
	    
//
//   Loop over hits
	    for(AliMUONHit* mHit=(AliMUONHit*)pMUON->FirstHit(-1); 
		mHit;
		mHit=(AliMUONHit*)pMUON->NextHit()) 
	    {
		Int_t   nch   = mHit->fChamber-1;  // chamber number
		if (nch > kNCH-1) continue;
		iChamber = &(pMUON->Chamber(nch));
                // new 17.07.99
		if (addBackground) {

		    if (mHit->fParticle == kMuonPlus 
			|| mHit->fParticle == kMuonMinus) {
			xhit[nch][nmuon[nch]]=mHit->fX;
			yhit[nch][nmuon[nch]]=mHit->fY;
			nmuon[nch]++;
			if (nmuon[nch] >2) printf("nmuon %d\n",nmuon[nch]);
		    }
		}
		


		
//
// Loop over pad hits
		for (AliMUONPadHit* mPad=
			 (AliMUONPadHit*)pMUON->FirstPad(mHit,fPadHits);
		     mPad;
		     mPad=(AliMUONPadHit*)pMUON->NextPad(fPadHits))
		{
		    Int_t cathode  = mPad->fCathode;    // cathode number
		    Int_t ipx      = mPad->fPadX;       // pad number on X
		    Int_t ipy      = mPad->fPadY;       // pad number on Y
		    Int_t iqpad    = Int_t(mPad->fQpad);// charge per pad
//
//
		    
		    if (cathode != (icat+1)) continue;
		    // fill the info array
		    Float_t thex, they;
		    segmentation=iChamber->SegmentationModel(cathode);
		    segmentation->GetPadCxy(ipx,ipy,thex,they);
//		    Float_t rpad=TMath::Sqrt(thex*thex+they*they);
//		    if (rpad < rmin || iqpad ==0 || rpad > rmax) continue;

		    new((*pAddress)[countadr++]) TVector(2);
		    TVector &trinfo=*((TVector*) (*pAddress)[countadr-1]);
		    trinfo(0)=(Float_t)track;
		    trinfo(1)=(Float_t)iqpad;

		    digits[0]=ipx;
		    digits[1]=ipy;
		    digits[2]=iqpad;
		    digits[3]=iqpad;
		    if (mHit->fParticle == kMuonPlus ||
			mHit->fParticle == kMuonMinus) {
			digits[4]=mPad->fHitNumber;
		    } else digits[4]=-1;

		    AliMUONTransientDigit* pdigit;
		    // build the list of fired pads and update the info
		    if (!hitMap[nch]->TestHit(ipx, ipy)) {

			list->AddAtAndExpand(
			    new AliMUONTransientDigit(nch,digits),counter);
			
			hitMap[nch]->SetHit(ipx, ipy, counter);
			counter++;
			pdigit=(AliMUONTransientDigit*)list->At(list->GetLast());
			// list of tracks
			TObjArray *trlist=(TObjArray*)pdigit->TrackList();
			trlist->Add(&trinfo);
		    } else {
			pdigit=(AliMUONTransientDigit*) hitMap[nch]->GetHit(ipx, ipy);
			// update charge
			(*pdigit).fSignal+=iqpad;
			(*pdigit).fPhysics+=iqpad;			
			// update list of tracks
			TObjArray* trlist=(TObjArray*)pdigit->TrackList();
			Int_t lastEntry=trlist->GetLast();
			TVector *pTrack=(TVector*)trlist->At(lastEntry);
			TVector &ptrk=*pTrack;
			Int_t lastTrack=Int_t(ptrk(0));
			Int_t lastCharge=Int_t(ptrk(1));
			if (lastTrack==track) {
			    lastCharge+=iqpad;
			    trlist->RemoveAt(lastEntry);
			    trinfo(0)=lastTrack;
			    trinfo(1)=lastCharge;
			    trlist->AddAt(&trinfo,lastEntry);
			} else {
			    trlist->Add(&trinfo);
			}
			// check the track list
			Int_t nptracks=trlist->GetEntriesFast();
			if (nptracks > 2) {
			    for (Int_t tr=0;tr<nptracks;tr++) {
				TVector *ppTrack=(TVector*)trlist->At(tr);
				TVector &pptrk=*ppTrack;
				trk[tr]=Int_t(pptrk(0));
				chtrk[tr]=Int_t(pptrk(1));
			    }
			} // end if nptracks
		    } //  end if pdigit
		} //end loop over clusters
	    } // hit loop
	} // track loop

	// open the file with background
       
	if (addBackground) {
	    ntracks =(Int_t)fTrH1->GetEntries();
//
//   Loop over tracks
//
	    for (Int_t track=0; track<ntracks; track++) {

		if (fHits2)       fHits2->Clear();
		if (fPadHits2)   fPadHits2->Clear();

		fTrH1->GetEvent(track);
//
//   Loop over hits
		AliMUONHit* mHit;
		for(int i=0;i<fHits2->GetEntriesFast();++i) 
		{	
		    mHit=(AliMUONHit*) (*fHits2)[i];
		    Int_t   nch   = mHit->fChamber-1;  // chamber number
		    if (nch >9) continue;
		    iChamber = &(pMUON->Chamber(nch));
		    Int_t rmin = (Int_t)iChamber->RInner();
		    Int_t rmax = (Int_t)iChamber->ROuter();
                    Float_t xbgr=mHit->fX;
		    Float_t ybgr=mHit->fY;
		    Bool_t cond=kFALSE;
		    
		    for (Int_t imuon =0; imuon < nmuon[nch]; imuon++) {
			Float_t dist= (xbgr-xhit[nch][imuon])*(xbgr-xhit[nch][imuon])
			    +(ybgr-yhit[nch][imuon])*(ybgr-yhit[nch][imuon]);
			if (dist<100) cond=kTRUE;
		    }
		    if (!cond) continue;
		    
//
// Loop over pad hits
		    for (AliMUONPadHit* mPad=
			     (AliMUONPadHit*)pMUON->FirstPad(mHit,fPadHits2);
			 mPad;
			 mPad=(AliMUONPadHit*)pMUON->NextPad(fPadHits2))
		    {
			//		    mPad = (AliMUONPadHit*) (*fPadHits2)[j];
			Int_t cathode  = mPad->fCathode;    // cathode number
			Int_t ipx      = mPad->fPadX;       // pad number on X
			Int_t ipy      = mPad->fPadY;       // pad number on Y
			Int_t iqpad    = Int_t(mPad->fQpad);// charge per pad

			if (cathode != (icat+1)) continue;
			Float_t thex, they;
			segmentation=iChamber->SegmentationModel(cathode);
			segmentation->GetPadCxy(ipx,ipy,thex,they);
			Float_t rpad=TMath::Sqrt(thex*thex+they*they);
			if (rpad < rmin || iqpad ==0 || rpad > rmax) continue;
			new((*pAddress)[countadr++]) TVector(2);
			TVector &trinfo=*((TVector*) (*pAddress)[countadr-1]);
			trinfo(0)=-1;  // tag background
			trinfo(1)=-1;
			
			digits[0]=ipx;
			digits[1]=ipy;
			digits[2]=iqpad;
			digits[3]=0;
			digits[4]=-1;
			
			AliMUONTransientDigit* pdigit;
			// build the list of fired pads and update the info
			if (!hitMap[nch]->TestHit(ipx, ipy)) {
			    list->AddAtAndExpand(new AliMUONTransientDigit(nch,digits),counter);
			    
			    hitMap[nch]->SetHit(ipx, ipy, counter);
			    counter++;
			    
			    pdigit=(AliMUONTransientDigit*)list->At(list->GetLast());
			    // list of tracks
			    TObjArray *trlist=(TObjArray*)pdigit->
				TrackList();
			    trlist->Add(&trinfo);
			} else {
			    pdigit=(AliMUONTransientDigit*) hitMap[nch]->GetHit(ipx, ipy);
			    // update charge
			    (*pdigit).fSignal+=iqpad;
			    
			    // update list of tracks
			    TObjArray* trlist=(TObjArray*)pdigit->
				TrackList();
			    Int_t lastEntry=trlist->GetLast();
			    TVector *pTrack=(TVector*)trlist->
				At(lastEntry);
			    TVector &ptrk=*pTrack;
			    Int_t lastTrack=Int_t(ptrk(0));
			    if (lastTrack==-1) {
				continue;
			    } else {
				trlist->Add(&trinfo);
			    }
				// check the track list
			    Int_t nptracks=trlist->GetEntriesFast();
			    if (nptracks > 0) {
				for (Int_t tr=0;tr<nptracks;tr++) {
				    TVector *ppTrack=(TVector*)trlist->At(tr);
				    TVector &pptrk=*ppTrack;
				    trk[tr]=Int_t(pptrk(0));
				    chtrk[tr]=Int_t(pptrk(1));
				}
			    } // end if nptracks
			} //  end if pdigit
		    } //end loop over clusters
		} // hit loop
	    } // track loop
	    //Int_t nentr2=list->GetEntriesFast();
	    //printf(" \n counter2, nentr2 %d %d \n",counter,nentr2);
	    TTree *fAli=gAlice->TreeK();
            TFile *file=NULL;
	    
	    if (fAli) file =fAli->GetCurrentFile();
	    file->cd();
	} // if addBackground	
	
	Int_t tracks[10];
	Int_t charges[10];
	Int_t nentries=list->GetEntriesFast();
	
	for (Int_t nent=0;nent<nentries;nent++) {
	    AliMUONTransientDigit *address=(AliMUONTransientDigit*)list->At(nent);
	    if (address==0) continue; 
	    Int_t ich=address->fChamber;
	    Int_t q=address->fSignal; 
	    iChamber=(AliMUONChamber*) (*fChambers)[ich];
//
//  Digit Response (noise, threshold, saturation, ...)
//		if (address->fPhysics !=0 ) address->fPhysics+=(Int_t)Noise; 
	    AliMUONResponse * response=iChamber->ResponseModel();
	    q=response->DigitResponse(q);
	    
	    if (!q) continue;
	    
	    digits[0]=address->fPadX;
	    digits[1]=address->fPadY;
	    digits[2]=q;
	    digits[3]=address->fPhysics;
	    digits[4]=address->fHit;
	    
	    TObjArray* trlist=(TObjArray*)address->TrackList();
	    Int_t nptracks=trlist->GetEntriesFast();
	    //printf("nptracks, trlist   %d  %p\n",nptracks,trlist);

	    // this was changed to accomodate the real number of tracks
	    if (nptracks > 10) {
		cout<<"Attention - nptracks > 10 "<<nptracks<<endl;
		nptracks=10;
	    }
	    if (nptracks > 2) {
		printf("Attention - nptracks > 2  %d \n",nptracks);
		printf("cat,ich,ix,iy,q %d %d %d %d %d \n",icat,ich,digits[0],digits[1],q);
	    }
	    for (Int_t tr=0;tr<nptracks;tr++) {
		TVector *ppP=(TVector*)trlist->At(tr);
		if(!ppP ) printf("ppP - %p\n",ppP);
		TVector &pp  =*ppP;
		tracks[tr]=Int_t(pp(0));
		charges[tr]=Int_t(pp(1));
                //printf("tracks, charges - %d %d\n",tracks[tr],charges[tr]);
	    }      //end loop over list of tracks for one pad
            // Sort list of tracks according to charge
	    if (nptracks > 1) {
		SortTracks(tracks,charges,nptracks);
	    }
	    if (nptracks < 10 ) {
		for (Int_t i=nptracks; i<10; i++) {
		    tracks[i]=0;
		    charges[i]=0;
		}
	    }
	    
	    // fill digits
	    pMUON->AddDigits(ich,tracks,charges,digits);
	    // delete trlist;
	}
	//cout<<"I'm out of the loops for digitisation"<<endl;
	//	gAlice->GetEvent(nev);
	gAlice->TreeD()->Fill();
	pMUON->ResetDigits();
	list->Delete();
	
	for(Int_t ii=0;ii<kNCH;++ii) {
	    if (hitMap[ii]) {
		hm=hitMap[ii];
		delete hm;
		hitMap[ii]=0;
	    }
	}
    } //end loop over cathodes
    
    char hname[30];
    sprintf(hname,"TreeD%d",nev);
    gAlice->TreeD()->Write(hname);
    // reset tree
    gAlice->TreeD()->Reset();
    delete list;
    
    pAddress->Delete();
    // gObjectTable->Print();
}

void AliMUON::SortTracks(Int_t *tracks,Int_t *charges,Int_t ntr)
{
  //
  // Sort the list of tracks contributing to a given digit
  // Only the 3 most significant tracks are acctually sorted
  //
  
  //
  //  Loop over signals, only 3 times
  //
  
  Int_t qmax;
  Int_t jmax;
  Int_t idx[3] = {-2,-2,-2};
  Int_t jch[3] = {-2,-2,-2};
  Int_t jtr[3] = {-2,-2,-2};
  Int_t i,j,imax;
  
  if (ntr<3) imax=ntr;
  else imax=3;
  for(i=0;i<imax;i++){
    qmax=0;
    jmax=0;
    
    for(j=0;j<ntr;j++){
      
      if((i == 1 && j == idx[i-1]) 
	 ||(i == 2 && (j == idx[i-1] || j == idx[i-2]))) continue;
      
      if(charges[j] > qmax) {
	qmax = charges[j];
	jmax=j;
      }       
    } 
    
    if(qmax > 0) {
      idx[i]=jmax;
      jch[i]=charges[jmax]; 
      jtr[i]=tracks[jmax]; 
    }
    
  } 
  
  for(i=0;i<3;i++){
    if (jtr[i] == -2) {
         charges[i]=0;
         tracks[i]=0;
    } else {
         charges[i]=jch[i];
         tracks[i]=jtr[i];
    }
  }

}

//___________________________________________
void AliMUON::Trigger(Int_t nev){
// call the Trigger Algorithm and fill TreeR

  Int_t singlePlus[3]  = {0,0,0}; 
  Int_t singleMinus[3] = {0,0,0}; 
  Int_t singleUndef[3] = {0,0,0};
  Int_t pairUnlike[3]  = {0,0,0}; 
  Int_t pairLike[3]    = {0,0,0};

  ResetTrigger();

  AliMUONTriggerDecision* decision= new AliMUONTriggerDecision(1);
  decision->Trigger();   
  decision->GetGlobalTrigger(singlePlus, singleMinus, singleUndef,
			     pairUnlike, pairLike);
// add a local trigger in the list 
  AddGlobalTrigger(singlePlus, singleMinus, singleUndef, pairUnlike, pairLike);
  
  for (Int_t icirc=0; icirc<kNTriggerCircuit; icirc++) { 
    if(decision->GetITrigger(icirc)==1) {
      Int_t localtr[7]={0,0,0,0,0,0,0};      
      Int_t loLpt[2]={0,0}; Int_t loHpt[2]={0,0}; Int_t loApt[2]={0,0};
      decision->GetLutOutput(icirc, loLpt, loHpt, loApt);
      localtr[0] = icirc;
      localtr[1] = decision->GetStripX11(icirc);
      localtr[2] = decision->GetDev(icirc);
      localtr[3] = decision->GetStripY11(icirc);
      for (Int_t i=0; i<2; i++) {    // convert the Lut output in 1 digit 
	localtr[4] = localtr[4]+Int_t(loLpt[i]*pow(2,i));
	localtr[5] = localtr[5]+Int_t(loHpt[i]*pow(2,i));
	localtr[6] = localtr[6]+Int_t(loApt[i]*pow(2,i));
      }
      //      cout << loApt[0] << " , " << loApt[1] << " , " << localtr[6] << "\n";
      AddLocalTrigger(localtr);  // add a local trigger in the list
    }
  }
  delete decision;

  gAlice->TreeR()->Fill();
  ResetTrigger();
  char hname[30];
  sprintf(hname,"TreeR%d",nev);
  gAlice->TreeR()->Write(hname);
  gAlice->TreeR()->Reset();
  printf("\n End of trigger for event %d", nev);
}


//____________________________________________
void AliMUON::FindClusters(Int_t nev,Int_t lastEntry)
{
    TClonesArray *dig1, *dig2;
    Int_t ndig, k;
    dig1 = new TClonesArray("AliMUONDigit",1000);
    dig2 = new TClonesArray("AliMUONDigit",1000);
    AliMUONDigit *digit;
//
// Loop on chambers and on cathode planes
//
    
    for (Int_t ich=0;ich<10;ich++) {
	AliMUONChamber* iChamber=(AliMUONChamber*) (*fChambers)[ich];
	AliMUONClusterFinder* rec = iChamber->ReconstructionModel();    
	gAlice->ResetDigits();
	gAlice->TreeD()->GetEvent(lastEntry);
	TClonesArray *muonDigits  = this->DigitsAddress(ich);
	ndig=muonDigits->GetEntriesFast();
	printf("\n 1 Found %d digits in %p %d", ndig, muonDigits,ich);
	TClonesArray &lhits1 = *dig1;
	Int_t n=0;
	for (k=0; k<ndig; k++) {
	    digit=(AliMUONDigit*) muonDigits->UncheckedAt(k);
	    if (rec->TestTrack(digit->fTracks[0]))
		new(lhits1[n++]) AliMUONDigit(*digit);
	}
	gAlice->ResetDigits();
	gAlice->TreeD()->GetEvent(lastEntry+1);
	muonDigits  = this->DigitsAddress(ich);
	ndig=muonDigits->GetEntriesFast();
	printf("\n 2 Found %d digits in %p %d", ndig, muonDigits, ich);
	TClonesArray &lhits2 = *dig2;
	n=0;
	
	for (k=0; k<ndig; k++) {
	    digit= (AliMUONDigit*) muonDigits->UncheckedAt(k);
	    if (rec->TestTrack(digit->fTracks[0]))
	    new(lhits2[n++]) AliMUONDigit(*digit);
	}

	if (rec) {	  
	    rec->SetDigits(dig1, dig2);
	    rec->SetChamber(ich);
	    rec->FindRawClusters();
	}
	dig1->Delete();
	dig2->Delete();
    } // for ich
    gAlice->TreeR()->Fill();
    ResetRawClusters();
    char hname[30];
    sprintf(hname,"TreeR%d",nev);
    gAlice->TreeR()->Write(hname);
    gAlice->TreeR()->Reset();
    printf("\n End of cluster finding for event %d", nev);
    
    delete dig1;
    delete dig2;
    //gObjectTable->Print();
}
 

void AliMUON::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliMUON.
      AliMUONChamber       *iChamber;
      AliMUONTriggerCircuit *iTriggerCircuit;
      AliMUONSegmentation  *segmentation;
      AliMUONResponse      *response;
      TClonesArray         *digitsaddress;
      TClonesArray         *rawcladdress;
      
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliDetector::Streamer(R__b);
      R__b >> fNPadHits;
      R__b >> fPadHits; // diff
      R__b >> fNLocalTrigger;       
      R__b >> fLocalTrigger;       
      R__b >> fNGlobalTrigger;       
      R__b >> fGlobalTrigger;   
      R__b >> fDchambers;
      R__b >> fRawClusters;
      R__b.ReadArray(fNdch);
      R__b.ReadArray(fNrawch);
      R__b >> fAccCut;
      R__b >> fAccMin;
      R__b >> fAccMax; 
      R__b >> fChambers;
      R__b >> fTriggerCircuits;
      for (Int_t i =0; i<kNTriggerCircuit; i++) {
	iTriggerCircuit=(AliMUONTriggerCircuit*) (*fTriggerCircuits)[i];
	iTriggerCircuit->Streamer(R__b);
      }
// Stream chamber related information
      for (Int_t i =0; i<kNCH; i++) {
	  iChamber=(AliMUONChamber*) (*fChambers)[i];
	  iChamber->Streamer(R__b);
	  if (iChamber->Nsec()==1) {
	      segmentation=iChamber->SegmentationModel(1);
	      if (segmentation)
	      segmentation->Streamer(R__b);
	  } else {
	      segmentation=iChamber->SegmentationModel(1);
	      if (segmentation)
	      segmentation->Streamer(R__b);
	      if (segmentation)
	      segmentation=iChamber->SegmentationModel(2);
	      segmentation->Streamer(R__b);
	  }
          response=iChamber->ResponseModel();
	  if (response)
	  response->Streamer(R__b);	  
	  digitsaddress=(TClonesArray*) (*fDchambers)[i];
	  digitsaddress->Streamer(R__b);
	  if (i < kNTrackingCh) {
	      rawcladdress=(TClonesArray*) (*fRawClusters)[i];
	      rawcladdress->Streamer(R__b);
	  }
      }
      
   } else {
      R__b.WriteVersion(AliMUON::IsA());
      AliDetector::Streamer(R__b);
      R__b << fNPadHits;
      R__b << fPadHits; // diff
      R__b << fNLocalTrigger;       
      R__b << fLocalTrigger;       
      R__b << fNGlobalTrigger;       
      R__b << fGlobalTrigger; 
      R__b << fDchambers;
      R__b << fRawClusters;
      R__b.WriteArray(fNdch, kNCH);
      R__b.WriteArray(fNrawch, kNTrackingCh);

      R__b << fAccCut;
      R__b << fAccMin;
      R__b << fAccMax; 

      R__b << fChambers;
      R__b << fTriggerCircuits;
      for (Int_t i =0; i<kNTriggerCircuit; i++) {
	iTriggerCircuit=(AliMUONTriggerCircuit*) (*fTriggerCircuits)[i];
	iTriggerCircuit->Streamer(R__b);
      }
      for (Int_t i =0; i<kNCH; i++) {
	  iChamber=(AliMUONChamber*) (*fChambers)[i];
	  iChamber->Streamer(R__b);
	  if (iChamber->Nsec()==1) {
	      segmentation=iChamber->SegmentationModel(1);
	      if (segmentation)
	      segmentation->Streamer(R__b);
	  } else {
	      segmentation=iChamber->SegmentationModel(1);
	      if (segmentation)
	      segmentation->Streamer(R__b);
	      segmentation=iChamber->SegmentationModel(2);
	      if (segmentation)
	      segmentation->Streamer(R__b);
	  }
          response=iChamber->ResponseModel();
	  if (response)
	  response->Streamer(R__b);
	  digitsaddress=(TClonesArray*) (*fDchambers)[i];
	  digitsaddress->Streamer(R__b);
	  if (i < kNTrackingCh) {
	      rawcladdress=(TClonesArray*) (*fRawClusters)[i];
	      rawcladdress->Streamer(R__b);
	  }
      }
   }
}
AliMUONPadHit* AliMUON::FirstPad(AliMUONHit*  hit, TClonesArray *clusters) 
{
//
    // Initialise the pad iterator
    // Return the address of the first padhit for hit
    TClonesArray *theClusters = clusters;
    Int_t nclust = theClusters->GetEntriesFast();
    if (nclust && hit->fPHlast > 0) {
	AliMUON::fMaxIterPad=hit->fPHlast;
	AliMUON::fCurIterPad=hit->fPHfirst;
	return (AliMUONPadHit*) clusters->UncheckedAt(AliMUON::fCurIterPad-1);
    } else {
	return 0;
    }
}

AliMUONPadHit* AliMUON::NextPad(TClonesArray *clusters) 
{
    AliMUON::fCurIterPad++;
    if (AliMUON::fCurIterPad <= AliMUON::fMaxIterPad) {
	return (AliMUONPadHit*) clusters->UncheckedAt(AliMUON::fCurIterPad-1);
    } else {
	return 0;
    }
}


AliMUONRawCluster *AliMUON::RawCluster(Int_t ichamber, Int_t icathod, Int_t icluster)
{
    TClonesArray *muonRawCluster  = RawClustAddress(ichamber);
    ResetRawClusters();
    TTree *treeR = gAlice->TreeR();
    Int_t nent=(Int_t)treeR->GetEntries();
    treeR->GetEvent(nent-2+icathod-1);
    //treeR->GetEvent(icathod);
    //Int_t nrawcl = (Int_t)muonRawCluster->GetEntriesFast();

    AliMUONRawCluster * mRaw = (AliMUONRawCluster*)muonRawCluster->UncheckedAt(icluster);
    //printf("RawCluster _ nent nrawcl icluster mRaw %d %d %d%p\n",nent,nrawcl,icluster,mRaw);
    
    return  mRaw;
}

AliMUON& AliMUON::operator = (const AliMUON& rhs)
{
// copy operator
// dummy version
    return *this;
}



















