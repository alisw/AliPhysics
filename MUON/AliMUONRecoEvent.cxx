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

//Authors: Mihaela Gheata, Andrei Gheata 09/10/00
////////////////////////////////////////////////////////////////////
//                                                                //
// AliMUONRecoEvent, AliMUONRecoTrack (and AliMUONRecoDisplay)    //
//                                                                //
// Theses classes are used to store and retrieve                  //
// MUON reconstructed events.                                     //
// The corresponding tree is constructed and filled               //
// during the FillEvent() method of AliMUONEventReconstructor,    //
// when all reconstruction information is available.              //
//                                                                //
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
//                                                                //
// AliMUONRecoEvent                                               //
//                                                                //
// This class handles an array of reconstructed tracks.           //
// It provides :                                                  //
//	- filling the tracks array according to the information   //
//        stored in AliMUONEventReconstructor class ;             //
//	- printing event and track informations : event numer,    //
//	  number of tracks, hits positions, reconstr. mometum.    //
//                                                                //
////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <AliRun.h>
#include <TClonesArray.h>
#include <TClass.h>

#include <TFile.h>
#include <TMatrixD.h>
#include <TParticle.h>

#include "AliMUONRecoEvent.h"
#include "AliMUONRecoTrack.h"
#include "AliMUONEventReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackK.h"
#include "AliMUONTrackParam.h"
#include "AliMUONHitForRec.h"
#include "AliMUONTrackHit.h"
#include "AliHeader.h"

ClassImp(AliMUONRecoEvent)

//-------------------------------------------------------------------
AliMUONRecoEvent::AliMUONRecoEvent(Int_t eventNo) 
{
// Reconstructed event constructor
   fTracks 	= new TClonesArray("AliMUONRecoTrack",200);
   fNevr 	= eventNo;
   fNtracks = 0;
}

//-------------------------------------------------------------------
AliMUONRecoEvent::~AliMUONRecoEvent() 
{
// Destructor of AliMUONRecoEvent
   fTracks->Delete();
   delete fTracks;
   fTracks = 0;
}

//-------------------------------------------------------------------
AliMUONRecoTrack* AliMUONRecoEvent::AddEmptyTrack()
{
// Add a empty AliMUONRecoTrackObject to the track list
   TClonesArray &dumptracks = *fTracks;
   return (new(dumptracks[fNtracks++])AliMUONRecoTrack(kTRUE));
}

//-------------------------------------------------------------------
void AliMUONRecoEvent::Clear(Option_t * /*option*/)
{
// Clears all track pointers from the list
//   fTracks->Clear(option);
   fTracks->Delete();
   fNtracks=0;
}

//-------------------------------------------------------------------
void AliMUONRecoEvent::EventInfo()
{
// Prints reconstructed event information
   cout << "*********************Reco Dumper**************************" << endl;
   cout << "Event number : " << fNevr << endl;
   cout << "   Number of tracks : " << fNtracks << endl;
   AliMUONRecoTrack *currentTrack =0;
   Int_t trackIndex = 0;
   for(trackIndex=0; trackIndex<fNtracks; trackIndex++) {
      currentTrack = (AliMUONRecoTrack*)fTracks->UncheckedAt(trackIndex);
      cout << "Track : " << trackIndex << endl;
      cout << "   Sign : " << currentTrack->GetSign() << endl;
      cout << "   Vertex position    : " << currentTrack->GetVertexPos() << endl;
      Double_t momreco[3];
      for (Int_t axis=0; axis<3; axis++) {
         momreco[axis] = currentTrack->GetMomReconstr(axis);
      }
      cout << "   Reconstructed mom. : " << "Px=" << momreco[0] << "Py=" << momreco[1] << "Pz=" << momreco[2] << endl;
      cout << "   Chi squared        : " << currentTrack->GetChi2r() << endl;
      cout << "   Hits positions     : " << endl;
      Double_t xhit, yhit, zhit;
      for (Int_t chamber=0; chamber<10; chamber++) {
         xhit = currentTrack->GetPosX(chamber);
         yhit = currentTrack->GetPosY(chamber);
         zhit = currentTrack->GetPosZ(chamber);
//         cout <<"      chamber" << chamber << " X=" << xhit << " Y=" << yhit << " Z=" << zhit << endl;
      }	 
   }
   cout << "**********************************************************" << endl;
}

//-------------------------------------------------------------------
Bool_t AliMUONRecoEvent::MakeDumpTracks(Int_t muons, TClonesArray *tracksPtr, 
  AliMUONEventReconstructor *EventReco)
{
// This method takes the pointer of the list of reconstructed tracks from
// AliMUONEventReconstructor and fill the reconstructed AliMUONRecoEvent
// fields.

	cout << "Enter MakeDumpTracks..." << endl;
   Int_t nTracks = tracksPtr->GetEntriesFast();
   cout << "nTracks = "<< nTracks << endl;
   if (nTracks == 0) {
      cout << "AliMUONRecoEvent::MakeDumpTracks: Number of tracks is zero !" << endl;
      //AZ return kFALSE;
   }
   cout << tracksPtr << endl;
   if (!tracksPtr) {
      cout << "AliMUONRecoEvent::MakeDumpTracks() : You should call SetRecoTracksPtr() first..." << endl;
      return kFALSE;
   }
	// Get event number
   Int_t noEvent = gAlice->GetHeader()->GetEvent();
   cout << "noEvent = "<< nTracks << endl;
   tracksPtr->Compress();  // simple loop
   AliMUONRecoTrack *currentTrack;
   Int_t trackIndex, nTrackHits = 0;
   Double_t z, pYZ, bendingSlope, nonBendingSlope;
   Double_t pX, pY, pZ;			// reconstructed momentum components
   Int_t isign, flag=0;       	// charge sign, flag of reconstructed track
   Double_t alpha, beta;
   TObjArray *hitsOnTrack = 0;
   AliMUONTrackHit *trackHit = 0;
   AliMUONTrack *track = 0;
   AliMUONTrackK *trackK = 0;
   TMatrixD *trackParamK; //AZnon
   AliMUONTrackParam *trackParam = 0;
   // Fill event number and number of tracks
   fNevr = noEvent;
   fMuons = muons; //AZ - number of muons within acceptance
   // Loop over reconstructed tracks
   for (trackIndex=0; trackIndex<nTracks; trackIndex++) {
      cout << " trackIndex = " << trackIndex << endl;
      currentTrack = AddEmptyTrack();
      cout << " currentTrack = " << currentTrack << endl;

      if (EventReco->GetTrackMethod() == 2) { // Kalman

        trackK = (AliMUONTrackK*) ((*tracksPtr)[trackIndex]);
	nTrackHits = trackK->GetNTrackHits();
	trackParamK = trackK->GetTrackParameters();
	isign = Int_t(TMath::Sign(1., (*trackParamK)(4,0)));
	z = trackK->GetZ();
	alpha = (*trackParamK)(2,0);
	beta = (*trackParamK)(3,0);
	pYZ = TMath::Cos(beta)/TMath::Abs((*trackParamK)(4,0));
	pZ = pYZ/TMath::Cos(alpha);
	pX = TMath::Sin(beta)/TMath::Abs((*trackParamK)(4,0));
	pY = pYZ*TMath::Sin(alpha);

	currentTrack->SetVertexPos(z);
	currentTrack->SetMomReconstr(pX,pY,pZ);
	currentTrack->SetSign(isign);
	currentTrack->SetChi2r(trackK->GetTrackQuality());

	// Check hits on the track
	hitsOnTrack = trackK->GetHitOnTrack();
	Float_t signal = 0;
	Float_t tht = 0;
	for (int ihit = 0; ihit < nTrackHits; ihit++) {
	  signal += ((AliMUONHitForRec*)((*hitsOnTrack)[ihit]))->GetGeantSignal();
	  tht += TMath::Min (1,((AliMUONHitForRec*)((*hitsOnTrack)[ihit]))->GetTHTrack());
	}
	signal /= nTrackHits;
	tht /= nTrackHits;
	flag = 0;
	if (TMath::Nint(signal) > 0) { // signal muon
	  for (int ihit = 0; ihit < nTrackHits ; ihit++) {
	    if (((AliMUONHitForRec*)((*hitsOnTrack)[ihit]))->GetTHTrack() != TMath::Nint(tht)) flag++;
	  }
	} else flag = -9; // background track
	//cout << TMath::Nint(signal) << " " << TMath::Nint(tht) << " " << recTrackNt->fFlag << endl;
      	currentTrack->SetFlag(flag);
      } else { // default tracking

        track = (AliMUONTrack*) ((*tracksPtr)[trackIndex]);
	nTrackHits = track->GetNTrackHits();
	// track parameters at Vertex
	trackParam = track->GetTrackParamAtVertex();
	bendingSlope = trackParam->GetBendingSlope();
	nonBendingSlope = trackParam->GetNonBendingSlope();

	z = trackParam->GetZ();
	pYZ = 1/TMath::Abs(trackParam->GetInverseBendingMomentum());
	pZ = pYZ/TMath::Sqrt(1+bendingSlope*bendingSlope);
	pX = pZ * nonBendingSlope;
	pY = pZ * bendingSlope;
	
	if (trackParam->GetInverseBendingMomentum()<0) isign=-1; else isign=1;
	currentTrack->SetVertexPos(z);
	currentTrack->SetMomReconstr(pX,pY,pZ);
	currentTrack->SetSign(isign);
	//         currentTrack->SetChi2r(trackParam->GetChi2());
	currentTrack->SetChi2r(0);

	// Check hits on the track
	hitsOnTrack = track->GetTrackHitsPtr();
	Float_t signal = 0;
	Float_t tht = 0;
	AliMUONHitForRec *hitForRec = 0;
	for (int ihit = 0; ihit < nTrackHits; ihit++) {
	  hitForRec = ((AliMUONTrackHit*)(*hitsOnTrack)[ihit])->GetHitForRecPtr();
	  signal += hitForRec->GetGeantSignal();
	  tht += TMath::Min (1,hitForRec->GetTHTrack());
	}
	signal /= nTrackHits;
	tht /= nTrackHits;
	flag = 0;
	if (TMath::Nint(signal) > 0) { // signal muon
	  for (int ihit = 0; ihit < nTrackHits ; ihit++) {
	    hitForRec = ((AliMUONTrackHit*)(*hitsOnTrack)[ihit])->GetHitForRecPtr();
	    if (hitForRec->GetTHTrack() != TMath::Nint(tht)) flag++;
	  }
	} else flag = -9; // background track
	//cout << TMath::Nint(signal) << " " << TMath::Nint(tht) << " " << recTrackNt->fFlag << endl;
      	currentTrack->SetFlag(flag);
      }
     
      Double_t xhit,yhit,zhit;
      // Loop over track hits
      for (Int_t trackHitIndex = 0; trackHitIndex < nTrackHits; trackHitIndex++) {
	if (EventReco->GetTrackMethod() == 2) { // Kalman
	  xhit = ((AliMUONHitForRec*)((*hitsOnTrack)[trackHitIndex]))->GetNonBendingCoor();
	  yhit = ((AliMUONHitForRec*)((*hitsOnTrack)[trackHitIndex]))->GetBendingCoor();
	  zhit = ((AliMUONHitForRec*)((*hitsOnTrack)[trackHitIndex]))->GetZ();    
	} else {
	  trackHit = (AliMUONTrackHit*) (*(track->GetTrackHitsPtr()))[trackHitIndex];
	  xhit = trackHit->GetHitForRecPtr()->GetNonBendingCoor();
	  yhit = trackHit->GetHitForRecPtr()->GetBendingCoor();
	  zhit = trackHit->GetHitForRecPtr()->GetZ();
	}
	if (trackHitIndex >= 0 && trackHitIndex < 10) {
	  currentTrack->SetHitPosition(trackHitIndex,xhit,yhit,zhit);
	} else { cout << "track " << trackIndex << " hit out of range" << endl;} 
      }
   
   }
   cout << "Leave MakeDumpTracks..." << endl;
   return kTRUE;
}

//-------------------------------------------------------------------
void AliMUONRecoEvent::Streamer(TBuffer &R__b)
{
// Streams an object of class AliMUONRecoEvent
   if (R__b.IsReading()) {
      fTracks->Clear();
      AliMUONRecoEvent::Class()->ReadBuffer(R__b, this);
   } else {
      cout << "...writing event to file...\n";
      AliMUONRecoEvent::Class()->WriteBuffer(R__b, this);
   }
}
