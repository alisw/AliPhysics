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
Revision 1.3  2000/12/21 17:51:54  morsch
RN3 violations corrected

Revision 1.2  2000/11/23 10:09:38  gosset
Bug correction in AliMUONRecoDisplay.
Copyright, $Log$
Copyright, Revision 1.3  2000/12/21 17:51:54  morsch
Copyright, RN3 violations corrected
Copyright,, $Id$, comments at the right place for automatic documentation,
in AliMUONRecoEvent and AliMUONRecoDisplay

*/

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
////////////////////////////////////////////////////////////////////
//                                                                //
// AliMUONRecoTrack                                               //
//                                                                //
// This class represents a reconstructed muon track               //
// in the tree of reconstructed events.                           //
//                                                                //
////////////////////////////////////////////////////////////////////

#include <iostream.h>
#include <AliRun.h>
#include <TClonesArray.h>
#include <TClass.h>
#include "AliMUONRecoEvent.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONHitForRec.h"
#include "AliMUONTrackHit.h"

ClassImp(AliMUONRecoTrack)
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
void AliMUONRecoEvent::Clear(Option_t *option)
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
Bool_t AliMUONRecoEvent::MakeDumpTracks(TClonesArray *tracksPtr)
{
// This method takes the pointer of the list of reconstructed tracks from
// AliMUONEventReconstructor and fill the reconstructed AliMUONRecoEvent
// fields.
	cout << "Enter MakeDumpTracks..." << endl;
   Int_t nTracks = tracksPtr->GetEntriesFast();
   if (nTracks == 0) {
      cout << "AliMUONRecoEvent::MakeDumpTracks: Number of tracks is zero !" << endl;
      return kFALSE;
   }
   cout << tracksPtr << endl;
   if (!tracksPtr) {
      cout << "AliMUONRecoEvent::MakeDumpTracks() : You should call SetRecoTracksPtr() first..." << endl;
      return kFALSE;
   }
	// Get event number
   Int_t noEvent = gAlice->GetHeader()->GetEvent();
   tracksPtr->Compress();  // simple loop
   AliMUONRecoTrack *currentTrack;
   Int_t trackIndex = 0, nTrackHits = 0;
   Double_t z,bendingSlope, nonBendingSlope, pYZ;
   Double_t pX, pY, pZ;			// reconstructed momentum components
   Int_t isign;					// charge sign
   AliMUONTrack *track = 0;
   AliMUONTrackParam *trackParam = 0;
   // Fill event number and number of tracks
   fNevr = noEvent;
   // Loop over reconstructed tracks
   for (trackIndex=0; trackIndex<nTracks; trackIndex++) {
      currentTrack = AddEmptyTrack();
      track = (AliMUONTrack*) ((*tracksPtr)[trackIndex]);
      nTrackHits = track->GetNTrackHits();
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
//	   currentTrack->SetChi2r(trackParam->GetChi2());
      currentTrack->SetChi2r(0);
      AliMUONTrackHit *trackHit;
      Double_t xhit,yhit,zhit;
	  // Loop over track hits
      for (Int_t trackHitIndex = 0; trackHitIndex < nTrackHits; trackHitIndex++) {
         trackHit = (AliMUONTrackHit*) (*(track->GetTrackHitsPtr()))[trackHitIndex];
         xhit = trackHit->GetHitForRecPtr()->GetNonBendingCoor();
         yhit = trackHit->GetHitForRecPtr()->GetBendingCoor();
         zhit = trackHit->GetHitForRecPtr()->GetZ();
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

//-------------------------------------------------------------------
AliMUONRecoTrack::AliMUONRecoTrack(Bool_t active)
{
//Constructor of AliMUONRecoTrack
   fSign  = 0;
   fZvr   = 0.0;
   fChi2r = 0.0;
   if (active) {
   	for (Int_t axis=0; axis<3; axis++) {
      	fPr[axis] = 0.0;
   	}
   	for (Int_t chamber=0; chamber<10; chamber++) {
      	fPosX[chamber] = 0.0;
      	fPosY[chamber] = 0.0;
      	fPosZ[chamber] = 0.0;
   	}
   }
}

//-------------------------------------------------------------------
const Double_t AliMUONRecoTrack::Phi()
{
// Return trach phi angle
	return TMath::ATan2(fPr[2], fPr[1]);
}

//-------------------------------------------------------------------
const Double_t AliMUONRecoTrack::Theta()
{
// Return trach theta angle
   return TMath::ACos(fPr[3] / P());
}

//-------------------------------------------------------------------
void AliMUONRecoTrack::SetMomReconstr(Double_t px, Double_t py, Double_t pz)
{
// Set the track reconstructed momentum 
   fPr[0] = px;
   fPr[1] = py;
   fPr[2] = pz;            
} 
   
//-------------------------------------------------------------------		
void AliMUONRecoTrack::SetHitPosition(Int_t chamber, Double_t x, Double_t y, Double_t z)
{
// Set hit coordinates in chamber[0..9]
   fPosX[chamber] = x;
   fPosY[chamber] = y;
   fPosZ[chamber] = z;
}
//-------------------------------------------------------------------		
void AliMUONRecoTrack::TrackInfo()
{
// Prints momentum info for this track
   cout << "Px=" << GetMomReconstr(0) << " Py=" << GetMomReconstr(1) <<
          " Pz=" << GetMomReconstr(2) << " P=" << P() << endl;
}
