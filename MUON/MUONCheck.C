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

//
// Macro for checking aliroot output and associated files contents
// Gines Martinez, Subatech June 2003
//
#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include "TBranch.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TParticle.h"
#include "TTree.h"

// STEER includes
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliStack.h"

// MUON includes
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONHit.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTrack.h"
#endif


void MUONkine(char * filename="galice.root",Int_t event2Check=0)
{
  // Stack of particle for each event
  AliStack* stack;
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }

  RunLoader->LoadKinematics("READ");
  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();

  for(ievent=0; ievent<nevents; ievent++) {  // Event loop
    if (event2Check!=0) ievent=event2Check;
    Int_t iparticle, nparticles;
    // Getting event ievent
    RunLoader->GetEvent(ievent); 
    stack = RunLoader->Stack();
    nparticles = (Int_t) stack->GetNtrack();
    printf(">>> Event %d, Number of particles is %d \n",ievent, nparticles);
    for(iparticle=0; iparticle<nparticles; iparticle++) {
      stack->Particle(iparticle)->Print("");  
    }
    if (event2Check!=0) ievent=nevents;
  }
  RunLoader->UnloadKinematics();
}


void MUONhits(char * filename="galice.root", Int_t event2Check=0)
{
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  // Loading MUON subsystem
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadHits("READ");  // Loading Tree of hits for MUON
  AliMUONData muondata(MUONLoader,"MUON","MUON");  // Creating MUON data container
  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();

  for(ievent=0; ievent<nevents; ievent++) {  // Event loop
    if (event2Check!=0) ievent=event2Check;
    printf(">>> Event %d \n",ievent);
    // Getting event ievent
    RunLoader->GetEvent(ievent); 
    muondata.SetTreeAddress("H");
    Int_t itrack, ntracks;
    ntracks = (Int_t) muondata.GetNtracks();
    for (itrack=0; itrack<ntracks; itrack++) { // Track loop
      printf(">>> Track %d \n",itrack);

      //Getting List of Hits of Track itrack
      muondata.GetTrack(itrack); 

      Int_t ihit, nhits;
      nhits = (Int_t) muondata.Hits()->GetEntriesFast();
      printf(">>> Number of hits  %d \n",nhits);
      AliMUONHit* mHit;
      for(ihit=0; ihit<nhits; ihit++) {
	mHit = static_cast<AliMUONHit*>(muondata.Hits()->At(ihit));
  	Int_t Nch      = mHit->Chamber();  // chamber number
	Int_t detele   = mHit-> DetElemId(); // Detection element if defined
	Int_t hittrack = mHit->Track();
	Float_t x      = mHit->X();
  	Float_t y      = mHit->Y();
  	Float_t z      = mHit->Z();
  	Float_t elos   = mHit->Eloss();
  	Float_t theta  = mHit->Theta();
  	Float_t phi    = mHit->Phi();
  	Float_t momentum = mHit->Momentum();
  	printf(">>> Hit %2d Chamber %2d DetEle %4d Track %4d x %6.3f y %6.3f z %7.3f elos %g  momentum %5.3f\n",
	       ihit, Nch, detele, hittrack,x,y,z,elos,momentum);
      }
      muondata.ResetHits();
    } // end track loop
    if (event2Check!=0) ievent=nevents;
  }  // end event loop
  MUONLoader->UnloadHits();
}


void MUONdigits(char * filename="galice.root", Int_t event2Check=0)
{
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  // Loading MUON subsystem
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadDigits("READ");
  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");

  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();
  AliMUONDigit * mDigit;
  
  for(ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    printf(">>> Event %d \n",ievent);
    RunLoader->GetEvent(ievent);
  
    // Addressing
    Int_t ichamber, nchambers;
    nchambers = AliMUONConstants::NCh(); ;
    muondata.SetTreeAddress("D"); 
    //    char branchname[30];    
 
    Int_t icathode, ncathodes;
    ncathodes=2;
    //Loop on cathodes 
    for(icathode=0; icathode<ncathodes; icathode++) {
      printf(">>> Cathode %d\n",icathode);
      muondata.GetCathode(icathode);
      // Loop on chambers
      for( ichamber=0; ichamber<nchambers; ichamber++) {
	printf(">>> Chamber %d\n",ichamber+1);
	
	Int_t idigit, ndigits;
	ndigits = (Int_t) muondata.Digits(ichamber)->GetEntriesFast();
	
	for(idigit=0; idigit<ndigits; idigit++) {
	  mDigit = static_cast<AliMUONDigit*>(muondata.Digits(ichamber)->At(idigit));
	  Int_t PadX   = mDigit->PadX();     // Pad X number
	  Int_t PadY   = mDigit->PadY();     // Pad Y number
	  Int_t Signal = mDigit->Signal();   // Physics Signal
	  Int_t Physics= mDigit->Physics();  // Physics contribution to signal
	  Int_t Hit    = mDigit->Hit();      // iHit
	  Int_t Cathode= mDigit->Cathode();  // Cathode
	  Int_t Track0 = mDigit->Track(0);
	  Int_t Track1 = mDigit->Track(1); 
	  Int_t Track2 = mDigit->Track(2);
	  Int_t TCharges0 = mDigit->TrackCharge(0);  //charge per track making this digit (up to 10)
	  Int_t TCharges1 = mDigit->TrackCharge(1);
	  Int_t TCharges2 = mDigit->TrackCharge(2);
	  Int_t idDE = mDigit->DetElemId();
		  
	  printf(">>>IdDE %d Digit %4d cathode %1d hit %4d PadX %3d PadY %3d Signal %4d Physics %4d Track0 %4d TrackCharge0 %4d Track1 %4d TrackCharge1 %4d Track2 %4d TrackCharge2 %4d \n",
		 idDE, idigit, Cathode,Hit, PadX, PadY, Signal, Physics, Track0, 
		 TCharges0, Track1, TCharges1, Track2, TCharges2);
	} // end digit loop
      } // end chamber loop
      muondata.ResetDigits();
    } // end cathode loop
    if (event2Check!=0) ievent=nevents;
  }  // end event loop
  MUONLoader->UnloadDigits();
}

void MUONrecpoints(char * filename="galice.root", Int_t event2Check=0) {

  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  // Getting MUONloader
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadRecPoints("READ");
  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");

  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();
  AliMUONRawCluster * mRecPoint = 0;
  
  for(ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    printf(">>> Event %d \n",ievent);
    RunLoader->GetEvent(ievent);
    // Addressing
    Int_t ichamber, nchambers;
    nchambers = AliMUONConstants::NTrackingCh();
    muondata.SetTreeAddress("RC"); 
    char branchname[30];    
    muondata.GetRawClusters();
    // Loop on chambers
    for( ichamber=0; ichamber<nchambers; ichamber++) {
      printf(">>> Chamber %d\n",ichamber);
      sprintf(branchname,"MUONRawClusters%d",ichamber+1);
      //printf(">>>  branchname %s\n",branchname);
      Int_t irecpoint, nrecpoints;
      nrecpoints = (Int_t) muondata.RawClusters(ichamber)->GetEntriesFast();
      printf("number of recpoints = %6d \n",nrecpoints);
      for(irecpoint=0; irecpoint<nrecpoints; irecpoint++) {
	mRecPoint = static_cast<AliMUONRawCluster*>(muondata.RawClusters(ichamber)->At(irecpoint));
//     Int_t       fTracks[3];        //labels of overlapped tracks
//     Int_t       fQ[2]  ;           // Q of cluster (in ADC counts)     
//     Float_t     fX[2]  ;           // X of cluster
//     Float_t     fY[2]  ;           // Y of cluster
//     Float_t     fZ[2]  ;           // Z of cluster
//     Int_t       fPeakSignal[2];    // Peak signal 
//     Int_t       fIndexMap[50][2];  // indeces of digits
//     Int_t       fOffsetMap[50][2]; // Emmanuel special
//     Float_t     fContMap[50][2];   // Contribution from digit
//     Int_t       fPhysicsMap[50];   // Distinguish signal and background contr.
//     Int_t       fMultiplicity[2];  // Cluster multiplicity
//     Int_t       fNcluster[2];      // Number of clusters
//     Int_t       fClusterType;      // Cluster type
//     Float_t     fChi2[2];          // Chi**2 of fit
//     Int_t       fGhost;            // 0 if not a ghost or ghost problem solved
//                                    // >0 if ghost problem remains because
//                                    // 1 both (true and ghost) satify 
//                                    //   charge chi2 compatibility
//                                    // 2 none give satisfactory chi2

	Int_t Track0 = mRecPoint->GetTrack(0);
	Int_t Track1 = mRecPoint->GetTrack(1); 
	Int_t Track2 = mRecPoint->GetTrack(2);
	Int_t Q0 = mRecPoint->GetCharge(0);
	Int_t Q1 = mRecPoint->GetCharge(1);
	Float_t x0 = mRecPoint->GetX(0);
	Float_t x1 = mRecPoint->GetX(1);
	Float_t y0 = mRecPoint->GetY(0);
	Float_t y1 = mRecPoint->GetY(1);
	Float_t z0 = mRecPoint->GetZ(0);
	Float_t z1 = mRecPoint->GetZ(1);
	Float_t chi2_0 =  mRecPoint->GetChi2(0);
	Float_t chi2_1 =  mRecPoint->GetChi2(1);
	//	Int_t de = mRecPoint->GetDetElementID();
	Int_t de = mRecPoint->DetElemId();
	//	printf(">>> RecPoint %4d x %6.3f %6.3f y %6.3f %6.3f z %6.3f %6.3f Q0 %4d  Q1 %4d Hit %4d Track1 %4d Track2 %4d Chi2 %6.3f %6.3f \n",
	//irecpoint, x0, x1, y0, y1, z0, z1, Q0, Q1, Track0, Track1, Track2, chi2_0, chi2_1);
	cout << mRecPoint->DetElemId() << endl;
	printf(">>> RecPoint %4d x %6.3f y %6.3f z %6.3f DetElem %2d \n",irecpoint,x0,y0,z0,de);
      } // end recpoint loop
    } // end chamber loop
    muondata.ResetRawClusters();
    if (event2Check!=0) ievent=nevents;
  }  // end event loop
  MUONLoader->UnloadRecPoints();
}

void MUONTestTrigger (char * filename="galice.root", Int_t event2Check=0){
// reads and dumps trigger objects from MUON.RecPoints.root
  TClonesArray * globalTrigger;
  TClonesArray * localTrigger;
  
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadDigits("READ");
  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");
  
  
  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();
  
  AliMUONGlobalTrigger *gloTrg;
  AliMUONLocalTrigger *locTrg;
  
  for (ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    RunLoader->GetEvent(ievent);
    
    muondata.SetTreeAddress("GLT"); 
    muondata.GetTriggerD();
    
    globalTrigger = muondata.GlobalTrigger();
    localTrigger = muondata.LocalTrigger();
    
    Int_t nglobals = (Int_t) globalTrigger->GetEntriesFast(); // should be 1
    Int_t nlocals  = (Int_t) localTrigger->GetEntriesFast(); // up to 234
    printf("###################################################\n");
    cout << " event " << ievent 
	 << " nglobals nlocals: " << nglobals << " " << nlocals << "\n"; 
    
    for (Int_t iglobal=0; iglobal<nglobals; iglobal++) { // Global Trigger
      gloTrg = static_cast<AliMUONGlobalTrigger*>(globalTrigger->At(iglobal));
      
      printf("===================================================\n");
      printf(" Global Trigger output       Low pt  High pt   All\n");
      printf(" number of Single Plus      :\t");
      printf("%i\t%i\t%i\t",gloTrg->SinglePlusLpt(),
	     gloTrg->SinglePlusHpt(),gloTrg->SinglePlusApt());
      printf("\n");
      printf(" number of Single Minus     :\t");
      printf("%i\t%i\t%i\t",gloTrg->SingleMinusLpt(),
	     gloTrg->SingleMinusHpt(),gloTrg->SingleMinusApt());
      printf("\n");
      printf(" number of Single Undefined :\t"); 
      printf("%i\t%i\t%i\t",gloTrg->SingleUndefLpt(),
	     gloTrg->SingleUndefHpt(),gloTrg->SingleUndefApt());
      printf("\n");
      printf(" number of UnlikeSign pair  :\t"); 
      printf("%i\t%i\t%i\t",gloTrg->PairUnlikeLpt(),
	     gloTrg->PairUnlikeHpt(),gloTrg->PairUnlikeApt());
      printf("\n");
      printf(" number of LikeSign pair    :\t");  
      printf("%i\t%i\t%i\t",gloTrg->PairLikeLpt(),
	     gloTrg->PairLikeHpt(),gloTrg->PairLikeApt());
      printf("\n");
      printf("===================================================\n");
      
    } // end of loop on Global Trigger
    
    for (Int_t ilocal=0; ilocal<nlocals; ilocal++) { // Local Trigger
      cout << " >>> Output for Local Trigger " << ilocal << "\n";
      
      locTrg = static_cast<AliMUONLocalTrigger*>(localTrigger->At(ilocal));
      
      cout << "Circuit StripX Dev StripY: " 
	   << locTrg->LoCircuit() << " " 
	   << locTrg->LoStripX() << " " 
	   << locTrg->LoDev() << " " 
	   << locTrg->LoStripY() 
	   << "\n";
      cout << "Lpt Hpt Apt: "     
	   << locTrg->LoLpt() << " "   
	   << locTrg->LoHpt() << " "  
	   << locTrg->LoApt() << "\n";
      
    } // end of loop on Local Trigger
    muondata.ResetTrigger();
    if (event2Check!=0) ievent=nevents;
  } // end loop on event  
  MUONLoader->UnloadRecPoints();
}



void MUONRecTracks (char * filename="galice.root", Int_t event2Check=0 ){
// reads and dumps trigger objects from MUON.RecPoints.root
  TClonesArray * RecTracks;
  
  // Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }
  
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadTracks("READ");
  // Creating MUON data container
  AliMUONData muondata(MUONLoader,"MUON","MUON");
  
    Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();
  
  //  AliMUONTrack * rectrack;
  
  for (ievent=0; ievent<nevents; ievent++) {
    if (event2Check!=0) ievent=event2Check;
    RunLoader->GetEvent(ievent);
    
    muondata.SetTreeAddress("RT");
    muondata.GetRecTracks();
    RecTracks = muondata.RecTracks();
    
    
    Int_t nrectracks = (Int_t) RecTracks->GetEntriesFast(); //

    printf(">>> Event %d Number of Recconstructed tracks %d \n",ievent, nrectracks);
   
    muondata.ResetRecTracks();
    if (event2Check!=0) ievent=nevents;
  } // end loop on event  
  MUONLoader->UnloadTracks();
}










