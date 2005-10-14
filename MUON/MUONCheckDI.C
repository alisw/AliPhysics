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
// with using AliMUONDataInterface
// By Bruce Becker, DAPNIA/SPhN/CEA Saclay
// According to MUONCheck.C
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
#include "AliMUONDataInterface.h"
#endif


void MUONkine(char * filename="galice.root")
{
  AliMUONDataInterface amdi;
  amdi.SetFile(filename);
  Int_t nevents = amdi.NumberOfEvents();
  
  for(Int_t ievent=0; ievent<nevents; ievent++) 
    {  // Event loop
      Int_t iparticle, nparticles;
      // Getting event ievent
      amdi.GetEvent(ievent); 
      nparticles = amdi.NumberOfParticles();
      printf(">>> Event %d, Number of particles is %d \n",ievent, nparticles);
      for(iparticle=0; iparticle<nparticles; iparticle++) 
	{
	amdi.Particle(iparticle)->Print("");
	}
    }
}


void MUONhits(char * filename="galice.root")
{
  
  AliMUONDataInterface amdi;
  amdi.SetFile(filename);
  Int_t nevents = amdi.NumberOfEvents();
  
  for(Int_t ievent=0; ievent<nevents; ievent++) 
    {  // Event loop
      printf(">>> Event %d \n",ievent);
      // Getting event ievent
      amdi.GetEvent(ievent); 
      Int_t ntracks = amdi.NumberOfTracks();
      for (Int_t itrack=0; itrack<ntracks; itrack++) 
	{ // Track loop
	  printf(">>> Track %d \n",itrack);
	  //Getting List of Hits of Track itrack
	  // amdi.GetTrack(itrack); 
	  
	  Int_t ihit, nhits;
	  nhits = amdi.NumberOfHits(itrack);
	  printf(">>> Number of hits  %d \n",nhits);
	  AliMUONHit* mHit;
	  for(ihit=0; ihit<nhits; ihit++) 
	    {
	      mHit = amdi.Hit(itrack,ihit);
	      Int_t Nch      = mHit->Chamber();  // chamber number
	      Int_t detele   = mHit-> DetElemId(); // Detection element if defined
	      Int_t hittrack = mHit->Track();
	      Float_t x      = mHit->X();
	      Float_t y      = mHit->Y();
	      Float_t z      = mHit->Z();
	      Float_t elos   = mHit->Eloss();
	      // Float_t theta  = mHit->Theta();
	      // Float_t phi    = mHit->Phi();
	      Float_t momentum = mHit->Momentum();
	      printf(">>> Hit %2d Chamber %2d DetEle %4d Track %4d x %6.3f y %6.3f z %7.3f elos %g  momentum %5.3f\n",
		     ihit, Nch, detele, hittrack,x,y,z,elos,momentum);
	    }
	} // end track loop
    }  // end event loop
}


void MUONdigits(char * filename="galice.root")
{
  
  AliMUONDataInterface amdi;
  amdi.SetFile(filename);
  Int_t ievent, nevents;
  nevents = amdi.NumberOfEvents();
  AliMUONDigit * mDigit;
  
  for(ievent=0; ievent<nevents; ievent++) 
    {
      printf(">>> Event %d \n",ievent);
      amdi.GetEvent(ievent);
      
      // Addressing
      Int_t ichamber, nchambers;
      nchambers = AliMUONConstants::NCh(); ;
      
      //      Int_t icathode, ncathodes;
      //      ncathodes=2;
      for( ichamber=0; ichamber<nchambers; ichamber++) 
	{
	  printf(">>> Chamber %d\n",ichamber+1);
	  
	  Int_t idigit, ndigits;
	  ndigits = amdi.NumberOfDigits(ichamber,0); // second parameter here is cathode...
	  
	  for(idigit=0; idigit<ndigits; idigit++) 
	    {
	      mDigit = amdi.Digit(ichamber,0,idigit);
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
	      //	  printf(">>> Cathode %d\n",Cathode);
	      
	      printf(">>>IdDE %d Digit %4d cathode %1d hit %4d PadX %3d PadY %3d Signal %4d Physics %4d Track0 %4d TrackCharge0 %4d Track1 %4d TrackCharge1 %4d Track2 %4d TrackCharge2 %4d \n",
		     idDE, idigit, Cathode,Hit, PadX, PadY, Signal, Physics, Track0, 
		     TCharges0, Track1, TCharges1, Track2, TCharges2);
	    } // end digit loop
	} // end chamber loop
    }  // end event loop
}

void MUONrecpoints(char * filename="galice.root") 
{
  AliMUONDataInterface amdi;
  amdi.SetFile(filename);
  
  Int_t ievent, nevents;
  nevents = amdi.NumberOfEvents();
  AliMUONRawCluster * mRecPoint = 0;
  
  for(ievent=0; ievent<nevents; ievent++) 
    {
      printf(">>> Event %d \n",ievent);
      amdi.GetEvent(ievent);
      Int_t ichamber, nchambers;
      nchambers = AliMUONConstants::NTrackingCh();
      // Loop on chambers
      for( ichamber=0; ichamber<nchambers; ichamber++) 
	{
	  printf(">>> Chamber %d\n",ichamber);
	  Int_t irecpoint, nrecpoints;
	  nrecpoints = amdi.NumberOfRawClusters(ichamber);
	  printf("number of recpoints = %6d \n",nrecpoints);
	  for(irecpoint=0; irecpoint<nrecpoints; irecpoint++) 
	    {
	      mRecPoint = amdi.RawCluster(ichamber,irecpoint);
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
	      Int_t de = mRecPoint->GetDetElemId();
	      printf(">>> RecPoint %4d  DetElem %2d x %6.3f %6.3f y %6.3f %6.3f z %6.3f %6.3f Q0 %4d  Q1 %4d Hit %4d Track1 %4d Track2 %4d Chi2 %6.3f %6.3f \n",
		     irecpoint,de,x0,x1,y0,y1,z0,z1,Q0,Q1,Track0, Track1, Track2, chi2_0, chi2_1);
	    } // end recpoint loop
	} // end chamber loop
    }  // end event loop
}

void MUONTestTrigger (char * filename="galice.root")
{
  // reads and dumps trigger objects from MUON.RecPoints.root
  
  AliMUONDataInterface amdi;
  amdi.SetFile(filename);
  
  Int_t ievent, nevents;
  nevents = amdi.NumberOfEvents();
  
  AliMUONGlobalTrigger *gloTrg;
  AliMUONLocalTrigger *locTrg;
  
  for (ievent=0; ievent<nevents; ievent++) 
    {
      amdi.GetEvent(ievent);
      Int_t nglobals = amdi.NumberOfGlobalTriggers(); // should be 1
      Int_t nlocals  = amdi.NumberOfLocalTriggers(); // up to 234
      printf("###################################################\n");
      cout << " event " << ievent 
	   << " nglobals nlocals: " << nglobals << " " << nlocals << "\n"; 
      
      for (Int_t iglobal=0; iglobal<nglobals; iglobal++) 
	{ // Global Trigger
	  gloTrg = amdi.GlobalTrigger(iglobal);
	  
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
      
      for (Int_t ilocal=0; ilocal<nlocals; ilocal++) 
	{ // Local Trigger
	  cout << " >>> Output for Local Trigger " << ilocal << "\n";
	  locTrg = amdi.LocalTrigger(ilocal);
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
    } // end loop on event  
}

void MUONRecTracks (char * filename="galice.root")
{
  
  // reads and dumps trigger objects from MUON.RecPoints.root
  AliMUONDataInterface amdi;
  amdi.SetFile(filename);
  AliMUONTrack* rectrack;
  AliMUONTrackParam* trackparam;
  Double_t x,y,z;
  Int_t ievent, nevents;
  nevents = amdi.NumberOfEvents();
  
  
  for (ievent=0; ievent<nevents; ievent++) 
    {
      amdi.GetEvent(ievent);
      Int_t nrectracks = amdi.NumberOfRecTracks(); 
      printf(">>> Event %d Number of Recconstructed tracks %d \n",ievent, nrectracks);
      // loop over rec tracks and print vertex parameters
      for(Int_t rectracki=0; rectracki < nrectracks;rectracki++)
	{
	  rectrack  = amdi.RecTrack(rectracki);
	  trackparam = rectrack->GetTrackParamAtVertex();
	  x = trackparam->GetNonBendingCoor();
	  y = trackparam->GetBendingCoor();
	  z = trackparam->GetZ();
	  printf("Track Vertex : (x,y,z) = (%f,%f,%f \n",x,y,z);
	} // end of loop over tracks
    } // end loop on event  
}
