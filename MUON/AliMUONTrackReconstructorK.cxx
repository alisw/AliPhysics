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

////////////////////////////////////
//
// MUON track reconstructor using the kalman method
//
// This class contains as data:
// * the parameters for the track reconstruction
//
// It contains as methods, among others:
// * MakeTracks to build the tracks
//
////////////////////////////////////

#include <stdlib.h>
#include <Riostream.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TMatrixD.h>

#include "AliMUONVTrackReconstructor.h"
#include "AliMUONTrackReconstructorK.h"
#include "AliMUONData.h"
#include "AliMUONConstants.h"
#include "AliMUONHitForRec.h"
#include "AliMUONRawCluster.h"
#include "AliMUONSegment.h"
#include "AliMUONTrackK.h" 
#include "AliLog.h"

ClassImp(AliMUONTrackReconstructorK) // Class implementation in ROOT context

//__________________________________________________________________________
AliMUONTrackReconstructorK::AliMUONTrackReconstructorK(AliMUONData* data, const Option_t* TrackMethod)
  : AliMUONVTrackReconstructor(data),
    fTrackMethod(2), //tracking method (2-Kalman 3-Combination-Kalman/Clustering)
    fMuons(0)
{
  /// Constructor for class AliMUONTrackReconstructorK

  if (strstr(TrackMethod,"Kalman")) {
    cout << " *** Tracking with the Kalman filter *** " << endl;
    fTrackMethod = 2;
  } else if (strstr(TrackMethod,"Combi")) {
    cout << " *** Combined cluster / track finder ***" << endl;
    fTrackMethod = 3;
  } else AliFatal(Form("Tracking method %s not available",TrackMethod));
  
  // Memory allocation for the TClonesArray of reconstructed tracks
  fRecTracksPtr = new TClonesArray("AliMUONTrackK", 10);
}

  //__________________________________________________________________________
AliMUONTrackReconstructorK::~AliMUONTrackReconstructorK(void)
{
  /// Destructor for class AliMUONTrackReconstructorK
  delete fRecTracksPtr;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::AddHitsForRecFromRawClusters()
{
  /// To add to the list of hits for reconstruction all the raw clusters
  TTree *treeR;
  AliMUONHitForRec *hitForRec;
  AliMUONRawCluster *clus;
  Int_t iclus, nclus, nTRentries;
  TClonesArray *rawclusters;
  AliDebug(1,"Enter AddHitsForRecFromRawClusters");

  treeR = fMUONData->TreeR();
  if (!treeR) {
    AliError("TreeR must be loaded");
    exit(0);
  }
  
  if (fTrackMethod != 3) { //AZ
    nTRentries = Int_t(treeR->GetEntries());
    
    if (!(fMUONData->IsRawClusterBranchesInTree())) {
      AliError(Form("RawCluster information is not avalaible, nTRentries = %d not equal to 1",nTRentries));
      exit(0);
    }

    fMUONData->SetTreeAddress("RC"); //AZ
    fMUONData->GetRawClusters(); // only one entry  
  }

  // Loop over tracking chambers
  for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
    // number of HitsForRec to 0 for the chamber
    fNHitsForRecPerChamber[ch] = 0;
    // index of first HitForRec for the chamber
    if (ch == 0) fIndexOfFirstHitForRecPerChamber[ch] = 0;
    else fIndexOfFirstHitForRecPerChamber[ch] = fNHitsForRec;
    rawclusters =fMUONData->RawClusters(ch);
    nclus = (Int_t) (rawclusters->GetEntries());
    // Loop over (cathode correlated) raw clusters
    for (iclus = 0; iclus < nclus; iclus++) {
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(iclus);
      // new AliMUONHitForRec from raw cluster
      // and increment number of AliMUONHitForRec's (total and in chamber)
      hitForRec = new ((*fHitsForRecPtr)[fNHitsForRec]) AliMUONHitForRec(clus);
      fNHitsForRec++;
      (fNHitsForRecPerChamber[ch])++;
      // more information into HitForRec
      hitForRec->SetBendingReso2(clus->GetErrY() * clus->GetErrY());
      hitForRec->SetNonBendingReso2(clus->GetErrX() * clus->GetErrX());
      //  original raw cluster
      hitForRec->SetChamberNumber(ch);
      hitForRec->SetHitNumber(iclus);
      // Z coordinate of the raw cluster (cm)
      hitForRec->SetZ(clus->GetZ(0));
      
      StdoutToAliDebug(3,
                       cout << "Chamber " << ch <<
                       " raw cluster  " << iclus << " : " << endl;
                       clus->Print("full");
                       cout << "AliMUONHitForRec number (1...): " << fNHitsForRec << endl;
                       hitForRec->Print("full");
                       );
    } // end of cluster loop
  } // end of chamber loop
  SortHitsForRecWithIncreasingChamber(); 
  
  AliDebug(1,"End of AddHitsForRecFromRawClusters");
    if (AliLog::GetGlobalDebugLevel() > 0) {
      AliDebug(1, Form("NHitsForRec: %d",fNHitsForRec));
      for (Int_t ch = 0; ch < AliMUONConstants::NTrackingCh(); ch++) {
	AliDebug(1, Form("Chamber(0...): %d",ch));
	AliDebug(1, Form("NHitsForRec: %d", fNHitsForRecPerChamber[ch]));
	AliDebug(1, Form("Index(first HitForRec): %d", fIndexOfFirstHitForRecPerChamber[ch]));
	for (Int_t hit = fIndexOfFirstHitForRecPerChamber[ch];
	     hit < fIndexOfFirstHitForRecPerChamber[ch] + fNHitsForRecPerChamber[ch];
	     hit++) {
	  AliDebug(1, Form("HitForRec index(0...): %d",hit));
	  ((*fHitsForRecPtr)[hit])->Dump();
      }
    }
  }
  
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::MakeSegments(void)
{
  /// To make the list of segments in all stations,
  /// from the list of hits to be reconstructed
  AliDebug(1,"Enter MakeSegments");
  //AZ ResetSegments();
  // Loop over stations
  for (Int_t st = 3; st < AliMUONConstants::NTrackingCh()/2; st++) MakeSegmentsPerStation(st); 
  
  StdoutToAliDebug(3,
    cout << "end of MakeSegments" << endl;
    for (Int_t st = 0; st < AliMUONConstants::NTrackingCh()/2; st++) 
    {
      cout << "station " << st
	    << "  has " << fNSegments[st] << " segments:"
	    << endl;
      for (Int_t seg = 0; seg < fNSegments[st]; seg++) 
      {
	      ((*fSegmentsPtr[st])[seg])->Print();
      }
    }
                   );
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::MakeTracks(void)
{
  /// To make the tracks from the list of segments and points in all stations
  AliDebug(1,"Enter MakeTracks");
  // The order may be important for the following Reset's
  //AZ ResetTracks();
  MakeTrackCandidates();
  if (fRecTracksPtr->GetEntriesFast() == 0) return;
  // Follow tracks in stations(1..) 3, 2 and 1
  FollowTracks();
  // Remove double tracks
  RemoveDoubleTracks();
  // Propagate tracks to the vertex thru absorber
  GoToVertex();
  // Fill AliMUONTrack data members
  FillMUONTrack();
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::MakeTrackCandidates(void)
{
  /// To make initial tracks for Kalman filter from the list of segments
  Int_t istat, iseg;
  AliMUONSegment *segment;
  AliMUONTrackK *trackK;

  AliDebug(1,"Enter MakeTrackCandidatesK");

  AliMUONTrackK a(this, fHitsForRecPtr);
  // Loop over stations(1...) 5 and 4
  for (istat=4; istat>=3; istat--) {
    // Loop over segments in the station
    for (iseg=0; iseg<fNSegments[istat]; iseg++) {
      // Transform segments to tracks and evaluate covariance matrix
      segment = (AliMUONSegment*) ((*fSegmentsPtr[istat])[iseg]);
      trackK = new ((*fRecTracksPtr)[fNRecTracks++]) AliMUONTrackK(segment);
    } // for (iseg=0;...)
  } // for (istat=4;...)
  return;
}

//__________________________________________________________________________
void AliMUONTrackReconstructorK::FollowTracks(void)
{
  /// Follow tracks using Kalman filter
  Bool_t ok;
  Int_t icand, ichamBeg = 0, ichamEnd, chamBits;
  Double_t zDipole1, zDipole2;
  AliMUONTrackK *trackK;
  AliMUONHitForRec *hit;
  AliMUONRawCluster *clus;
  TClonesArray *rawclusters;
  clus = 0; rawclusters = 0;

  zDipole1 = fSimpleBPosition + fSimpleBLength/2;
  zDipole2 = zDipole1 - fSimpleBLength;

  // Print hits
  trackK = (AliMUONTrackK*) ((*fRecTracksPtr)[0]);

  if (trackK->DebugLevel() > 0) {
    for (Int_t i1=0; i1<fNHitsForRec; i1++) {
      hit = (AliMUONHitForRec*) ((*fHitsForRecPtr)[i1]);
      printf(" Hit # %d %10.4f %10.4f %10.4f",
             hit->GetChamberNumber(), hit->GetBendingCoor(),
             hit->GetNonBendingCoor(), hit->GetZ());
 
      // from raw clusters
      rawclusters = fMUONData->RawClusters(hit->GetChamberNumber());
      clus = (AliMUONRawCluster*) rawclusters->UncheckedAt(hit->
							   GetHitNumber());
      printf(" %d", clus->GetTrack(1));
      if (clus->GetTrack(2) != -1) printf(" %d \n", clus->GetTrack(2));
      else printf("\n");
     
    }
  } // if (trackK->DebugLevel() > 0)

  icand = -1;
  Int_t nSeeds;
  nSeeds = fNRecTracks; // starting number of seeds
  // Loop over track candidates
  while (icand < fNRecTracks-1) {
    icand ++;
    if (trackK->DebugLevel()>0) cout << " *** Kalman track candidate No. " << icand << endl;
    trackK = (AliMUONTrackK*) ((*fRecTracksPtr)[icand]);
    if (trackK->GetRecover() < 0) continue; // failed track

    // Discard candidate which will produce the double track
    /*
    if (icand > 0) {
      ok = CheckCandidate(icand,nSeeds);
      if (!ok) {
        trackK->SetRecover(-1); // mark candidate to be removed
        continue;
      }
    }
    */

    ok = kTRUE;
    if (trackK->GetRecover() == 0) 
      hit = (AliMUONHitForRec*) trackK->GetTrackHits()->Last(); // last hit
    else 
      hit = trackK->GetHitLastOk(); // hit where track stopped

    if (hit) ichamBeg = hit->GetChamberNumber();
    ichamEnd = 0;
    // Check propagation direction
    if (!hit) { ichamBeg = ichamEnd; AliFatal(" ??? "); }
    else if (trackK->GetTrackDir() < 0) {
      ichamEnd = 9; // forward propagation
      ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kFALSE,zDipole1,zDipole2);
      if (ok) {
        ichamBeg = ichamEnd;
        ichamEnd = 6; // backward propagation
	// Change weight matrix and zero fChi2 for backpropagation
        trackK->StartBack();
	trackK->SetTrackDir(1);
        ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kTRUE,zDipole1,zDipole2);
        ichamBeg = ichamEnd;
        ichamEnd = 0;
      }
    } else {
      if (trackK->GetBPFlag()) {
	// backpropagation
        ichamEnd = 6; // backward propagation
	// Change weight matrix and zero fChi2 for backpropagation
        trackK->StartBack();
        ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kTRUE,zDipole1,zDipole2);
        ichamBeg = ichamEnd;
        ichamEnd = 0;
      }
    }

    if (ok) {
      trackK->SetTrackDir(1);
      trackK->SetBPFlag(kFALSE);
      ok = trackK->KalmanFilter(ichamBeg,ichamEnd,kFALSE,zDipole1,zDipole2);
    }
    if (!ok) { trackK->SetRecover(-1); continue; } // mark candidate to be removed

    // Apply smoother
    if (trackK->GetRecover() >= 0) {
      ok = trackK->Smooth();
      if (!ok) trackK->SetRecover(-1); // mark candidate to be removed
    }

    // Majority 3 of 4 in first 2 stations
    if (!ok) continue;
    chamBits = 0;
    Double_t chi2max = 0;
    for (Int_t i=0; i<trackK->GetNTrackHits(); i++) {
      hit = (AliMUONHitForRec*) (*trackK->GetTrackHits())[i];
      chamBits |= BIT(hit->GetChamberNumber());
      if (trackK->GetChi2PerPoint(i) > chi2max) chi2max = trackK->GetChi2PerPoint(i);
    }
    if (!((chamBits&3)==3 || (chamBits>>2&3)==3) && chi2max > 25) {
      //trackK->Recover();
      trackK->SetRecover(-1); //mark candidate to be removed
      continue;
    }
    if (ok) trackK->SetTrackQuality(0); // compute "track quality"
  } // while

  for (Int_t i=0; i<fNRecTracks; i++) {
    trackK = (AliMUONTrackK*) ((*fRecTracksPtr)[i]);
    if (trackK->GetRecover() < 0) fRecTracksPtr->RemoveAt(i);
  }

  // Compress TClonesArray
  fRecTracksPtr->Compress();
  fNRecTracks = fRecTracksPtr->GetEntriesFast();
  return;
}

//__________________________________________________________________________
Bool_t AliMUONTrackReconstructorK::CheckCandidate(Int_t icand, Int_t nSeeds) const
{
  /// Discards track candidate if it will produce the double track (having
  /// the same seed segment hits as hits of a good track found before)
  AliMUONTrackK *track1, *track2;
  AliMUONHitForRec *hit1, *hit2, *hit;

  track1 = (AliMUONTrackK*) ((*fRecTracksPtr)[icand]);
  hit1 = (AliMUONHitForRec*) (*track1->GetTrackHits())[0]; // 1'st hit
  hit2 = (AliMUONHitForRec*) (*track1->GetTrackHits())[1]; // 2'nd hit

  for (Int_t i=0; i<icand; i++) {
    track2 = (AliMUONTrackK*) ((*fRecTracksPtr)[i]);
    //if (track2->GetRecover() < 0) continue;
    if (track2->GetRecover() < 0 && icand >= nSeeds) continue;

    if (track1->GetStartSegment() == track2->GetStartSegment()) {
      return kFALSE;
    } else {
      Int_t nSame = 0;
      for (Int_t j=0; j<track2->GetNTrackHits(); j++) {
        hit = (AliMUONHitForRec*) (*track2->GetTrackHits())[j];
        if (hit == hit1 || hit == hit2) {
          nSame++;
          if (nSame == 2) return kFALSE;
        }
      } // for (Int_t j=0;
    }
  } // for (Int_t i=0;
  return kTRUE;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::RemoveDoubleTracks(void)
{
  /// Removes double tracks (sharing more than half of their hits).
  /// Keeps the track with higher quality
  AliMUONTrackK *track1, *track2, *trackToKill;

  // Sort tracks according to their quality
  fRecTracksPtr->Sort();

  // Loop over first track of the pair
  track1 = (AliMUONTrackK*) fRecTracksPtr->First();
  Int_t debug = track1->DebugLevel();
  while (track1) {
    // Loop over second track of the pair
    track2 = (AliMUONTrackK*) fRecTracksPtr->After(track1);
    while (track2) {
      // Check whether or not to keep track2
      if (!track2->KeepTrack(track1)) {
        if (debug >= 0) cout << " Killed track: " << 1/(*track2->GetTrackParameters())(4,0) <<
	  " " << track2->GetTrackQuality() << endl;
        trackToKill = track2;
        track2 = (AliMUONTrackK*) fRecTracksPtr->After(track2);
        trackToKill->Kill();
        fRecTracksPtr->Compress();
      } else track2 = (AliMUONTrackK*) fRecTracksPtr->After(track2);
    } // track2
    track1 = (AliMUONTrackK*) fRecTracksPtr->After(track1);
  } // track1

  fNRecTracks = fRecTracksPtr->GetEntriesFast();
  if (debug >= 0) cout << " Number of Kalman tracks: " << fNRecTracks << endl;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::GoToVertex(void)
{
  /// Propagates track to the vertex thru absorber
  /// (using Branson correction for now)
  Double_t zVertex;
  zVertex = 0;
  for (Int_t i=0; i<fNRecTracks; i++) {
    //((AliMUONTrackK*)(*fRecTracksPtr)[i])->Branson();
    ((AliMUONTrackK*)(*fRecTracksPtr)[i])->SetTrackQuality(1); // compute Chi2
    //((AliMUONTrackK*)(*fRecTracksPtr)[i])->GoToZ(zVertex); // w/out absorber
    ((AliMUONTrackK*)(*fRecTracksPtr)[i])->GoToVertex(1); // with absorber
  }
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::FillMUONTrack()
{
  // Set track parameters at hits for Kalman track. Fill fTrackParamAtHit of AliMUONTrack's
  AliMUONTrackK *track;
  track = (AliMUONTrackK*) fRecTracksPtr->First();
  while (track) {
    track->FillMUONTrack();
    track = (AliMUONTrackK*) fRecTracksPtr->After(track);
  } 
  return;
}

  //__________________________________________________________________________
void AliMUONTrackReconstructorK::EventDump(void)
{
  /// Dump reconstructed event (track parameters at vertex and at first hit),
  /// and the particle parameters
  AliDebug(1,"****** enter EventDump ******");
  AliDebug(1, Form("Number of Reconstructed tracks : %d", fNRecTracks)); 
  
  AliWarning("AliMUONTrackReconstructorK::EventDump: not implemented"); 
  
  // informations about generated particles NO !!!!!!!!
  
//    for (Int_t iPart = 0; iPart < np; iPart++) {
//      p = gAlice->Particle(iPart);
//      printf(" particle %d: type= %d px= %f py= %f pz= %f pdg= %d\n",
//  	   iPart, p->GetPdgCode(), p->Px(), p->Py(), p->Pz(), p->GetPdgCode());    
//    }
  return;
}


