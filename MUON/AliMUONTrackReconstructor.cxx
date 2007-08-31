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

//-----------------------------------------------------------------------------
/// \class AliMUONTrackReconstructor
/// MUON track reconstructor using the original method
///
/// This class contains as data:
/// - the parameters for the track reconstruction
///
/// It contains as methods, among others:
/// - MakeTracks to build the tracks
//-----------------------------------------------------------------------------

#include "AliMUONTrackReconstructor.h"

#include "AliMUONConstants.h"
#include "AliMUONHitForRec.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliLog.h"

#include <TMinuit.h>
#include <Riostream.h>
#include <TMath.h>
#include <TMatrixD.h>

// Functions to be minimized with Minuit
void TrackChi2(Int_t &nParam, Double_t *gradient, Double_t &chi2, Double_t *param, Int_t flag);

/// \cond CLASSIMP
ClassImp(AliMUONTrackReconstructor) // Class implementation in ROOT context
/// \endcond

//************* Parameters for reconstruction
const Double_t AliMUONTrackReconstructor::fgkBendingVertexDispersion = 10.;
const Double_t AliMUONTrackReconstructor::fgkNonBendingVertexDispersion = 10.;


  //__________________________________________________________________________
AliMUONTrackReconstructor::AliMUONTrackReconstructor()
  : AliMUONVTrackReconstructor()
{
  /// Constructor for class AliMUONTrackReconstructor
  AliInfo("*** Original tracking ***");
}

  //__________________________________________________________________________
AliMUONTrackReconstructor::~AliMUONTrackReconstructor()
{
/// Destructor
} 

  //__________________________________________________________________________
void AliMUONTrackReconstructor::MakeTrackCandidates()
{
  /// To make track candidates (assuming linear propagation if the flag fgkMakeTrackCandidatesFast is set to kTRUE):
  /// Start with segments station(1..) 4 or 5 then follow track in station 5 or 4.
  /// Good candidates are made of at least three hitForRec's.
  /// Keep only best candidates or all of them according to the flag fgkTrackAllTracks.
  TClonesArray *segments;
  AliMUONTrack *track;
  Int_t iCandidate = 0;
  Bool_t hitFound;

  AliDebug(1,"Enter MakeTrackCandidates");

  // Loop over stations(1..) 5 and 4 and make track candidates
  for (Int_t istat=4; istat>=3; istat--) {
    
    // Make segments in the station
    segments = MakeSegmentsInStation(istat);
    
    // Loop over segments
    for (Int_t iseg=0; iseg<segments->GetEntriesFast(); iseg++) 
    {
      AliDebug(1,Form("Making primary candidate(1..) %d",++iCandidate));
      
      // Transform segments to tracks and put them at the end of fRecTracksPtr
      track = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack((AliMUONObjectPair*)((*segments)[iseg]));
      fNRecTracks++;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2))
      {
        cout<<endl<<"Track parameter covariances at first hit with multiple Coulomb scattering effects:"<<endl;
        ((AliMUONTrackParam*) track->GetTrackParamAtHit()->First())->GetCovariances().Print();
      }
      
      // Look for compatible hitForRec(s) in the other station
      if (fgkMakeTrackCandidatesFast) hitFound = FollowLinearTrackInStation(*track,7-istat);
      else hitFound = FollowTrackInStation(*track,7-istat);
      
      // Remove track if no hit found
      if (!hitFound) {
        fRecTracksPtr->Remove(track);
  	fNRecTracks--;
      }
      
    }
    
    // delete the array of segments
    delete segments;
  }
  
  fRecTracksPtr->Compress(); // this is essential before checking tracks
  
  // Keep all different tracks or only the best ones as required
  if (fgkTrackAllTracks) RemoveIdenticalTracks();
  else RemoveDoubleTracks();
  
  AliDebug(1,Form("Number of good candidates = %d",fNRecTracks));
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::FollowTracks()
{
  /// Follow tracks in stations(1..) 3, 2 and 1
  AliDebug(1,"Enter FollowTracks");
  
  AliMUONTrack *track, *nextTrack;
  Int_t currentNRecTracks;
  Bool_t hitFound;
  
  for (Int_t station = 2; station >= 0; station--) {
    
    // Save the actual number of reconstructed track in case of
    // tracks are added or suppressed during the tracking procedure
    // !! Do not compress fRecTracksPtr until the end of the loop over tracks !!
    currentNRecTracks = fNRecTracks;
    
    for (Int_t iRecTrack = 0; iRecTrack <currentNRecTracks; iRecTrack++) {
      AliDebug(1,Form("FollowTracks: track candidate(1..) %d", iRecTrack+1));
      
      track = (AliMUONTrack*) fRecTracksPtr->UncheckedAt(iRecTrack);
      
      // Fit the track:
      // Do not take into account the multiple scattering to speed up the fit
      // Calculate the track parameter covariance matrix
      // If "station" is station(1..) 3 then use the vertex to better constrain the fit
      if (station==2) {
        SetVertexForFit(*track);
        track->SetFitWithVertex(kTRUE);
      } else track->SetFitWithVertex(kFALSE);
      Fit(*track, kFALSE, kTRUE);
      
      // Remove the track if the normalized chi2 is too high
      if (track->GetNormalizedChi2() > fgkSigmaToCutForTracking * fgkSigmaToCutForTracking) {
  	fRecTracksPtr->Remove(track);
  	fNRecTracks--;
  	continue;
      }
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
        cout<<endl<<"Track parameter covariances at first hit with multiple Coulomb scattering effects:"<<endl;
        ((AliMUONTrackParam*) track->GetTrackParamAtHit()->First())->GetCovariances().Print();
      }
      
      // Look for compatible hitForRec in station(0..) "station"
      hitFound = FollowTrackInStation(*track,station);
      
      // Try to recover track if required
      if (!hitFound && fgkRecoverTracks) hitFound = RecoverTrack(*track,station);
      
      // remove track if no hit found
      if (!hitFound) {
	fRecTracksPtr->Remove(track);
  	fNRecTracks--;
      }
      
    }
    
    // Compress fRecTracksPtr for the next step
    fRecTracksPtr->Compress();
    
    // Keep only the best tracks if required
    if (!fgkTrackAllTracks) RemoveDoubleTracks();
    
  }
  
  // Last fit of track candidates with all station
  // Take into account the multiple scattering and remove bad tracks
  Int_t trackIndex = -1;
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    
    trackIndex++;
    nextTrack = (AliMUONTrack*) fRecTracksPtr->After(track); // prepare next track
    
    track->SetFitWithVertex(kFALSE); // just to be sure
    Fit(*track, kTRUE, kTRUE);
    
    // Printout for debuging
    if (AliLog::GetGlobalDebugLevel() >= 3) {
      cout << "FollowTracks: track candidate(0..) " << trackIndex << " after final fit" << endl;
      track->RecursiveDump();
    } 
    
    // Remove the track if the normalized chi2 is too high
    if (track->GetNormalizedChi2() > fgkSigmaToCutForTracking * fgkSigmaToCutForTracking) {
      fRecTracksPtr->Remove(track);
      fNRecTracks--;
    }
    
    track = nextTrack;
    
  }
  
  fRecTracksPtr->Compress();
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::FollowTrackInStation(AliMUONTrack &trackCandidate, Int_t nextStation)
{
  /// Follow trackCandidate in station(0..) nextStation and search for compatible HitForRec(s)
  /// Keep all possibilities or only the best one(s) according to the flag fgkTrackAllTracks:
  /// kTRUE:  duplicate "trackCandidate" if there are several possibilities and add the new tracks at the end of
  ///         fRecTracksPtr to avoid conficts with other track candidates at this current stage of the tracking procedure.
  ///         Remove the obsolete "trackCandidate" at the end.
  /// kFALSE: add only the best hit(s) to the "trackCandidate". Try to add a couple of hits in priority.
  AliDebug(1,Form("Enter FollowTrackInStation(1..) %d", nextStation+1));
  
  // Order the chamber according to the propagation direction (tracking starts with chamber 2):
  // - nextStation == station(1...) 5 => forward propagation
  // - nextStation < station(1...) 5 => backward propagation
  Int_t ch1, ch2;
  if (nextStation==4) {
    ch1 = 2*nextStation+1;
    ch2 = 2*nextStation;
  } else {
    ch1 = 2*nextStation;
    ch2 = 2*nextStation+1;
  }
  
  Double_t chi2WithOneHitForRec = 1.e10;
  Double_t chi2WithTwoHitForRec = 1.e10;
  Double_t maxChi2WithOneHitForRec = 2. * fgkSigmaToCutForTracking * fgkSigmaToCutForTracking; // 2 because 2 quantities in chi2
  Double_t maxChi2WithTwoHitForRec = 4. * fgkSigmaToCutForTracking * fgkSigmaToCutForTracking; // 4 because 4 quantities in chi2
  Double_t bestChi2WithOneHitForRec = maxChi2WithOneHitForRec;
  Double_t bestChi2WithTwoHitForRec = maxChi2WithTwoHitForRec;
  Bool_t foundOneHit = kFALSE;
  Bool_t foundTwoHits = kFALSE;
  AliMUONTrack *newTrack = 0x0;
  AliMUONHitForRec *hitForRecCh1, *hitForRecCh2;
  AliMUONTrackParam extrapTrackParam;
  AliMUONTrackParam extrapTrackParamAtHit1;
  AliMUONTrackParam extrapTrackParamAtHit2;
  AliMUONTrackParam bestTrackParamAtHit1;
  AliMUONTrackParam bestTrackParamAtHit2;
  Bool_t *hitForRecCh1Used = new Bool_t[fNHitsForRecPerChamber[ch1]];
  for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) hitForRecCh1Used[hit1] = kFALSE;
  
  // Get track parameters
  AliMUONTrackParam extrapTrackParamAtCh(*(AliMUONTrackParam*)trackCandidate.GetTrackParamAtHit()->First());
  
  // Add MCS effect
  AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(),1.);
  
  // Add MCS in the missing chamber if any (only 1 chamber can be missing according to tracking criteria)
  if (ch1 < ch2 && extrapTrackParamAtCh.GetHitForRecPtr()->GetChamberNumber() > ch2 + 1) {
    // extrapolation to the missing chamber
    AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(ch2 + 1));
    // add MCS effect
    AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(),1.);
  }
  
  //Extrapolate trackCandidate to chamber "ch2"
  AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(ch2));
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    cout<<endl<<"Track parameter covariances at first hit extrapolated to z = "<<AliMUONConstants::DefaultChamberZ(ch2)<<":"<<endl;
    extrapTrackParamAtCh.GetCovariances().Print();
  }
  
  // Printout for debuging
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
    cout << "FollowTrackInStation: look for hits in chamber(1..): " << ch2+1 << endl;
  }
  
  // look for candidates in chamber 2 
  for (Int_t hit2 = 0; hit2 < fNHitsForRecPerChamber[ch2]; hit2++) {
    
    hitForRecCh2 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch2]+hit2);
    
    // try to add the current hit fast
    if (!TryOneHitForRecFast(extrapTrackParamAtCh, hitForRecCh2)) continue;
    
    // try to add the current hit accuratly
    chi2WithOneHitForRec = TryOneHitForRec(extrapTrackParamAtCh, hitForRecCh2, extrapTrackParamAtHit2);
    
    // if good chi2 then try to attach a hitForRec in the other chamber too
    if (chi2WithOneHitForRec < maxChi2WithOneHitForRec) {
      Bool_t foundSecondHit = kFALSE;
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: found one hit in chamber(1..): " << ch2+1
	     << " (Chi2 = " << chi2WithOneHitForRec << ")" << endl;
        cout << "                      look for second hits in chamber(1..): " << ch1+1 << " ..." << endl;
      }
      
      // add MCS effect for next step
      AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtHit2,AliMUONConstants::ChamberThicknessInX0(),1.);
      
      // copy new track parameters for next step
      extrapTrackParam = extrapTrackParamAtHit2;
      
      //Extrapolate track parameters to chamber "ch1"
      AliMUONTrackExtrap::ExtrapToZ(&extrapTrackParam, AliMUONConstants::DefaultChamberZ(ch1));
      
      for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) {
        
	hitForRecCh1 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch1]+hit1);
	
    	// try to add the current hit fast
    	if (!TryOneHitForRecFast(extrapTrackParam, hitForRecCh1)) continue;
    	
    	// try to add the current hit accuratly
	chi2WithTwoHitForRec = TryTwoHitForRec(extrapTrackParamAtHit2, hitForRecCh1, extrapTrackParamAtHit1);
        
	// if good chi2 then create a new track by adding the 2 hitForRec to the "trackCandidate"
	if (chi2WithTwoHitForRec < maxChi2WithTwoHitForRec) {
	  foundSecondHit = kTRUE;
          foundTwoHits = kTRUE;
          
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowTrackInStation: found second hit in chamber(1..): " << ch1+1
	  	 << " (Global Chi2 = " << chi2WithTwoHitForRec << ")" << endl;
	  }
	  
	  if (fgkTrackAllTracks) {
	    // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
            newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	    UpdateTrack(*newTrack,extrapTrackParamAtHit1,extrapTrackParamAtHit2);
	    fNRecTracks++;
	    
	    // Tag hitForRecCh1 as used
	    hitForRecCh1Used[hit1] = kTRUE;
	    
	    // Printout for debuging
	    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	      cout << "FollowTrackInStation: added two hits in station(1..): " << nextStation+1 << endl;
	      if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	    }
	    
          } else if (chi2WithTwoHitForRec < bestChi2WithTwoHitForRec) {
	    // keep track of the best couple of hits
	    bestChi2WithTwoHitForRec = chi2WithTwoHitForRec;
	    bestTrackParamAtHit1 = extrapTrackParamAtHit1;
	    bestTrackParamAtHit2 = extrapTrackParamAtHit2;
          }
	  
	}
	
      }
      
      // if no hitForRecCh1 found then consider to add hitForRecCh2 only
      if (!foundSecondHit) {
        foundOneHit = kTRUE;
        
	if (fgkTrackAllTracks) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
          newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  UpdateTrack(*newTrack,extrapTrackParamAtHit2);
	  fNRecTracks++;
	  
	  // Printout for debuging
	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
	    cout << "FollowTrackInStation: added one hit in chamber(1..): " << ch2+1 << endl;
	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
	  }
	  
	} else if (!foundTwoHits && chi2WithOneHitForRec < bestChi2WithOneHitForRec) {
	  // keep track of the best single hitForRec except if a couple of hits has already been found
	  bestChi2WithOneHitForRec = chi2WithOneHitForRec;
	  bestTrackParamAtHit1 = extrapTrackParamAtHit2;
        }
	
      }
      
    }
    
  }
  
  // look for candidates in chamber 1 not already attached to a track
  // if we want to keep all possible tracks or if no good couple of hitForRec has been found
  if (fgkTrackAllTracks || !foundTwoHits) {
    
    // Printout for debuging
    if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
      cout << "FollowTrackInStation: look for single hits in chamber(1..): " << ch1+1 << endl;
    }
    
    // add MCS effect for next step
    AliMUONTrackExtrap::AddMCSEffect(&extrapTrackParamAtCh,AliMUONConstants::ChamberThicknessInX0(),1.);
    
    //Extrapolate trackCandidate to chamber "ch1"
    AliMUONTrackExtrap::ExtrapToZCov(&extrapTrackParamAtCh, AliMUONConstants::DefaultChamberZ(ch1));
    
    for (Int_t hit1 = 0; hit1 < fNHitsForRecPerChamber[ch1]; hit1++) {
      
      hitForRecCh1 = (AliMUONHitForRec*) fHitsForRecPtr->UncheckedAt(fIndexOfFirstHitForRecPerChamber[ch1]+hit1);
      
      if (hitForRecCh1Used[hit1]) continue; // Skip hitForRec already used
      
      // try to add the current hit fast
      if (!TryOneHitForRecFast(extrapTrackParamAtCh, hitForRecCh1)) continue;
      
      // try to add the current hit accuratly
      chi2WithOneHitForRec = TryOneHitForRec(extrapTrackParamAtCh, hitForRecCh1, extrapTrackParamAtHit1);
    
      // if good chi2 then consider to add hitForRecCh1
      // We do not try to attach a hitForRec in the other chamber too since it has already been done above
      if (chi2WithOneHitForRec < maxChi2WithOneHitForRec) {
        foundOneHit = kTRUE;
	  
	// Printout for debuging
  	if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	  cout << "FollowTrackInStation: found one hit in chamber(1..): " << ch1+1
  	       << " (Chi2 = " << chi2WithOneHitForRec << ")" << endl;
  	}
	
  	if (fgkTrackAllTracks) {
	  // copy trackCandidate into a new track put at the end of fRecTracksPtr and add the new hitForRec's
  	  newTrack = new ((*fRecTracksPtr)[fRecTracksPtr->GetLast()+1]) AliMUONTrack(trackCandidate);
	  UpdateTrack(*newTrack,extrapTrackParamAtHit1);
	  fNRecTracks++;
  	  
	  // Printout for debuging
  	  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
  	    cout << "FollowTrackInStation: added one hit in chamber(1..): " << ch1+1 << endl;
  	    if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
  	  }
	  
  	} else if (chi2WithOneHitForRec < bestChi2WithOneHitForRec) {
	  // keep track of the best single hitForRec except if a couple of hits has already been found
  	  bestChi2WithOneHitForRec = chi2WithOneHitForRec;
	  bestTrackParamAtHit1 = extrapTrackParamAtHit1;
  	}
	
      }
      
    }
    
  }
  
  // fill out the best track if required else clean up the fRecTracksPtr array
  if (!fgkTrackAllTracks) {
    if (foundTwoHits) {
      UpdateTrack(trackCandidate,bestTrackParamAtHit1,bestTrackParamAtHit2);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: added the two best hits in station(1..): " << nextStation+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else if (foundOneHit) {
      UpdateTrack(trackCandidate,bestTrackParamAtHit1);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "FollowTrackInStation: added the best hit in chamber(1..): " << bestTrackParamAtHit1.GetHitForRecPtr()->GetChamberNumber()+1 << endl;
        if (AliLog::GetGlobalDebugLevel() >= 3) newTrack->RecursiveDump();
      }
      
    } else {
      delete [] hitForRecCh1Used;
      return kFALSE;
    }
    
  } else if (foundOneHit || foundTwoHits) {
    
    // remove obsolete track
    fRecTracksPtr->Remove(&trackCandidate);
    fNRecTracks--;
    
  } else {
    delete [] hitForRecCh1Used;
    return kFALSE;
  }
  
  delete [] hitForRecCh1Used;
  return kTRUE;
  
}

  //__________________________________________________________________________
Double_t AliMUONTrackReconstructor::TryTwoHitForRec(const AliMUONTrackParam &trackParamAtHit1, AliMUONHitForRec* hitForRec2,
						    AliMUONTrackParam &trackParamAtHit2)
{
/// Test the compatibility between the track and the 2 hitForRec together (using trackParam's covariance matrix):
/// return the corresponding Chi2 accounting for covariances between the 2 hitForRec
/// return trackParamAtHit1 & 2
  
  // extrapolate track parameters at the z position of the second hit
  trackParamAtHit2.SetParameters(trackParamAtHit1.GetParameters());
  trackParamAtHit2.SetZ(trackParamAtHit1.GetZ());
  AliMUONTrackExtrap::ExtrapToZ(&trackParamAtHit2, hitForRec2->GetZ());
  
  // set pointer to hit2 into trackParamAtHit2
  trackParamAtHit2.SetHitForRecPtr(hitForRec2);
  
  // Set differences between track and the 2 hitForRec in the bending and non bending directions
  AliMUONHitForRec* hitForRec1 = trackParamAtHit1.GetHitForRecPtr();
  TMatrixD dPos(4,1);
  dPos(0,0) = hitForRec1->GetNonBendingCoor() - trackParamAtHit1.GetNonBendingCoor();
  dPos(1,0) = hitForRec1->GetBendingCoor() - trackParamAtHit1.GetBendingCoor();
  dPos(2,0) = hitForRec2->GetNonBendingCoor() - trackParamAtHit2.GetNonBendingCoor();
  dPos(3,0) = hitForRec2->GetBendingCoor() - trackParamAtHit2.GetBendingCoor();
  
  // Calculate the error matrix from the track parameter covariances at first hitForRec
  TMatrixD error(4,4);
  error.Zero();
  if (trackParamAtHit1.CovariancesExist()) {
    // Save track parameters at first hitForRec
    AliMUONTrackParam trackParamAtHit1Save(trackParamAtHit1);
    TMatrixD paramAtHit1Save(trackParamAtHit1Save.GetParameters());
    Double_t z1 = trackParamAtHit1Save.GetZ();
    
    // Save track coordinates at second hitForRec
    Double_t nonBendingCoor2	     = trackParamAtHit2.GetNonBendingCoor();
    Double_t bendingCoor2  	     = trackParamAtHit2.GetBendingCoor();
    
    // add MCS effect at first hitForRec
    AliMUONTrackExtrap::AddMCSEffect(&trackParamAtHit1Save,AliMUONConstants::ChamberThicknessInX0(),1.);
    
    // Get the pointer to the parameter covariance matrix at first hitForRec
    const TMatrixD& kParamCov = trackParamAtHit1Save.GetCovariances();
    
    // Calculate the jacobian related to the transformation between track parameters
    // at first hitForRec and track coordinates at the 2 hitForRec z-position
    TMatrixD jacob(4,5);
    jacob.Zero();
    // first derivative at the first hitForRec:
    jacob(0,0) = 1.; // dx1/dx
    jacob(1,2) = 1.; // dy1/dy
    // first derivative at the second hitForRec:
    TMatrixD dParam(5,1);
    for (Int_t i=0; i<5; i++) {
      // Skip jacobian calculation for parameters with no associated error
      if (kParamCov(i,i) == 0.) continue;
      // Small variation of parameter i only
      for (Int_t j=0; j<5; j++) {
        if (j==i) {
          dParam(j,0) = TMath::Sqrt(kParamCov(i,i));
	  if (j == 4) dParam(j,0) *= TMath::Sign(1.,-paramAtHit1Save(4,0)); // variation always in the same direction
        } else dParam(j,0) = 0.;
      }
      
      // Set new track parameters at first hitForRec
      trackParamAtHit1Save.SetParameters(paramAtHit1Save);
      trackParamAtHit1Save.AddParameters(dParam);
      trackParamAtHit1Save.SetZ(z1);
      
      // Extrapolate new track parameters to the z position of the second hitForRec
      AliMUONTrackExtrap::ExtrapToZ(&trackParamAtHit1Save,hitForRec2->GetZ());
      
      // Calculate the jacobian
      jacob(2,i) = (trackParamAtHit1Save.GetNonBendingCoor()  - nonBendingCoor2) / dParam(i,0); // dx2/dParami
      jacob(3,i) = (trackParamAtHit1Save.GetBendingCoor()     - bendingCoor2   ) / dParam(i,0); // dy2/dParami
    }
    
    // Calculate the error matrix
    TMatrixD tmp(jacob,TMatrixD::kMult,kParamCov);
    error = TMatrixD(tmp,TMatrixD::kMultTranspose,jacob);
  }
  
  // Add hitForRec resolution to the error matrix
  error(0,0) += hitForRec1->GetNonBendingReso2();
  error(1,1) += hitForRec1->GetBendingReso2();
  error(2,2) += hitForRec2->GetNonBendingReso2();
  error(3,3) += hitForRec2->GetBendingReso2();
  
  // invert the error matrix for Chi2 calculation
  if (error.Determinant() != 0) {
    error.Invert();
  } else {
    AliWarning(" Determinant error=0");
    return 1.e10;
  }
  
  // Compute the Chi2 value
  TMatrixD tmp2(dPos,TMatrixD::kTransposeMult,error);
  TMatrixD result(tmp2,TMatrixD::kMult,dPos);
  
  return result(0,0);
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtHit)
{
  /// Add 1 hit to the track candidate
  /// Update chi2 of the track 
  
  // Compute local chi2
  AliMUONHitForRec* hit = trackParamAtHit.GetHitForRecPtr();
  Double_t deltaX = trackParamAtHit.GetNonBendingCoor() - hit->GetNonBendingCoor();
  Double_t deltaY = trackParamAtHit.GetBendingCoor() - hit->GetBendingCoor();
  Double_t localChi2 = deltaX*deltaX / hit->GetNonBendingReso2() +
  		       deltaY*deltaY / hit->GetBendingReso2();
  
  // Flag hit as being not removable
  trackParamAtHit.SetRemovable(kFALSE);
  trackParamAtHit.SetLocalChi2(0.); // --> Local chi2 not used
  
  // Update the chi2 of the new track
  track.SetFitFMin(track.GetFitFMin() + localChi2);
  
  // Update TrackParamAtHit
  track.AddTrackParamAtHit(&trackParamAtHit,trackParamAtHit.GetHitForRecPtr());
  track.GetTrackParamAtHit()->Sort();
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::UpdateTrack(AliMUONTrack &track, AliMUONTrackParam &trackParamAtHit1, AliMUONTrackParam &trackParamAtHit2)
{
  /// Add 2 hits to the track candidate
  /// Update track and local chi2
  
  // Update local chi2 at first hit
  AliMUONHitForRec* hit1 = trackParamAtHit1.GetHitForRecPtr();
  Double_t deltaX = trackParamAtHit1.GetNonBendingCoor() - hit1->GetNonBendingCoor();
  Double_t deltaY = trackParamAtHit1.GetBendingCoor() - hit1->GetBendingCoor();
  Double_t localChi2AtHit1 = deltaX*deltaX / hit1->GetNonBendingReso2() +
  			     deltaY*deltaY / hit1->GetBendingReso2();
  trackParamAtHit1.SetLocalChi2(localChi2AtHit1);
  
  // Flag first hit as being removable
  trackParamAtHit1.SetRemovable(kTRUE);
  
  // Update local chi2 at second hit
  AliMUONHitForRec* hit2 = trackParamAtHit2.GetHitForRecPtr();
  deltaX = trackParamAtHit2.GetNonBendingCoor() - hit2->GetNonBendingCoor();
  deltaY = trackParamAtHit2.GetBendingCoor() - hit2->GetBendingCoor();
  Double_t localChi2AtHit2 = deltaX*deltaX / hit2->GetNonBendingReso2() +
  			     deltaY*deltaY / hit2->GetBendingReso2();
  trackParamAtHit2.SetLocalChi2(localChi2AtHit2);
  
  // Flag first hit as being removable
  trackParamAtHit2.SetRemovable(kTRUE);
  
  // Update the chi2 of the new track
  track.SetFitFMin(track.GetFitFMin() + localChi2AtHit1 + localChi2AtHit2);
  
  // Update TrackParamAtHit
  track.AddTrackParamAtHit(&trackParamAtHit1,trackParamAtHit1.GetHitForRecPtr());
  track.AddTrackParamAtHit(&trackParamAtHit2,trackParamAtHit2.GetHitForRecPtr());
  track.GetTrackParamAtHit()->Sort();
  
}

  //__________________________________________________________________________
Bool_t AliMUONTrackReconstructor::RecoverTrack(AliMUONTrack &trackCandidate, Int_t nextStation)
{
  /// Try to recover the track candidate in the next station
  /// by removing the worst of the two hits attached in the current station
  /// Return kTRUE if recovering succeeds
  AliDebug(1,"Enter RecoverTrack");
  
  // Do not try to recover track until we have attached hit(s) on station(1..) 3
  if (nextStation > 1) return kFALSE;
  
  Int_t worstHitNumber = -1;
  Double_t localChi2, worstLocalChi2 = 0.;
  
  // Look for the hit to remove
  for (Int_t hitNumber = 0; hitNumber < 2; hitNumber++) {
    AliMUONTrackParam *trackParamAtHit = (AliMUONTrackParam*)trackCandidate.GetTrackParamAtHit()->UncheckedAt(hitNumber);
    
    // check if current hit is removable
    if (!trackParamAtHit->IsRemovable()) return kFALSE;
    
    // Pick up hit with the worst chi2
    localChi2 = trackParamAtHit->GetLocalChi2();
    if (localChi2 > worstLocalChi2) {
      worstLocalChi2 = localChi2;
      worstHitNumber = hitNumber;
    }
  }
  
  // Reset best hit as being NOT removable
  ((AliMUONTrackParam*)trackCandidate.GetTrackParamAtHit()->UncheckedAt((worstHitNumber+1)%2))->SetRemovable(kFALSE);
  
  // Remove the worst hit
  trackCandidate.RemoveTrackParamAtHit((AliMUONTrackParam*)trackCandidate.GetTrackParamAtHit()->UncheckedAt(worstHitNumber));
  
  // Re-fit the track:
  // Do not take into account the multiple scattering to speed up the fit
  // Calculate the track parameter covariance matrix
  trackCandidate.SetFitWithVertex(kFALSE); // To be sure
  Fit(trackCandidate, kFALSE, kTRUE);
  
  // Look for new hit(s) in next station
  return FollowTrackInStation(trackCandidate,nextStation);
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::SetVertexForFit(AliMUONTrack &trackCandidate)
{
  /// Add the vertex as a measured hit to constrain the fit of the "trackCandidate"
  /// Compute the vertex resolution from natural vertex dispersion and
  /// multiple scattering effets according to trackCandidate path in absorber
  /// It is necessary to account for multiple scattering effects here instead of during the fit of
  /// the "trackCandidate" to do not influence the result by changing track resolution at vertex
  AliDebug(1,"Enter SetVertexForFit");
  
  Double_t nonBendingReso2 = fgkNonBendingVertexDispersion * fgkNonBendingVertexDispersion;
  Double_t bendingReso2 = fgkBendingVertexDispersion * fgkBendingVertexDispersion;
  // add multiple scattering effets
  AliMUONTrackParam paramAtVertex(*((AliMUONTrackParam*)(trackCandidate.GetTrackParamAtHit()->First())));
  paramAtVertex.DeleteCovariances(); // to be sure to account only for multiple scattering
  AliMUONTrackExtrap::ExtrapToVertexUncorrected(&paramAtVertex,0.);
  const TMatrixD& kParamCov = paramAtVertex.GetCovariances();
  nonBendingReso2 += kParamCov(0,0);
  bendingReso2 += kParamCov(2,2);
  // Set the vertex
  AliMUONHitForRec vertex; // Coordinates set to (0.,0.,0.) by default
  vertex.SetNonBendingReso2(nonBendingReso2);
  vertex.SetBendingReso2(bendingReso2);
  trackCandidate.SetVertex(&vertex);
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::Fit(AliMUONTrack &track, Bool_t includeMCS, Bool_t calcCov)
{
  /// Fit the track "track" w/wo multiple Coulomb scattering according to "includeMCS".
  
  Double_t benC, errorParam, invBenP, nonBenC, x, y;
  AliMUONTrackParam *trackParam;
  Double_t arg[1], fedm, errdef, fitFMin;
  Int_t npari, nparx;
  Int_t status, covStatus;
  
  // Instantiate gMinuit if not already done
  if (!gMinuit) gMinuit = new TMinuit(6);
  // Clear MINUIT parameters
  gMinuit->mncler();
  // Give the fitted track to MINUIT
  gMinuit->SetObjectFit(&track);
  if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 2) || (AliLog::GetGlobalDebugLevel() >= 2)) {
    // Define print level
    arg[0] = 1;
    gMinuit->mnexcm("SET PRI", arg, 1, status);
    // Print covariance matrix
    gMinuit->mnexcm("SHO COV", arg, 0, status);
  } else {
    arg[0] = -1;
    gMinuit->mnexcm("SET PRI", arg, 1, status);
  }
  // No warnings
  gMinuit->mnexcm("SET NOW", arg, 0, status);
  // Define strategy
  //arg[0] = 2;
  //gMinuit->mnexcm("SET STR", arg, 1, status);
  
  // set flag w/wo multiple scattering according to "includeMCS"
  track.SetFitWithMCS(includeMCS);
  if (includeMCS) {
    // compute hit weights only once
    if (!track.ComputeHitWeights()) {
      AliWarning("cannot take into account the multiple scattering effects");
      track.SetFitWithMCS(kFALSE);
    }
  }
  
  // Set fitting function
  gMinuit->SetFCN(TrackChi2);
  
  // Set fitted parameters (!! The order is very important for the covariance matrix !!)
  trackParam = (AliMUONTrackParam*) (track.GetTrackParamAtHit()->First());
  // could be tried with no limits for the search (min=max=0) ????
  // mandatory limits in non Bending to avoid NaN values of parameters
  gMinuit->mnparm(0, "X", trackParam->GetNonBendingCoor(), 0.03, -500.0, 500.0, status);
  gMinuit->mnparm(1, "NonBenS", trackParam->GetNonBendingSlope(), 0.001, -0.5, 0.5, status);
  // mandatory limits in Bending to avoid NaN values of parameters
  gMinuit->mnparm(2, "Y", trackParam->GetBendingCoor(), 0.10, -500.0, 500.0, status);
  gMinuit->mnparm(3, "BenS", trackParam->GetBendingSlope(), 0.001, -0.5, 0.5, status);
  gMinuit->mnparm(4, "InvBenP", trackParam->GetInverseBendingMomentum(), 0.003, -0.4, 0.4, status);
  
  // minimization
  gMinuit->mnexcm("MIGRAD", arg, 0, status);
  
  // Calculate the covariance matrix more accurately if required
  if (calcCov) gMinuit->mnexcm("HESSE", arg, 0, status);
   
  // get results into "invBenP", "benC", "nonBenC" ("x", "y")
  gMinuit->GetParameter(0, x, errorParam);
  trackParam->SetNonBendingCoor(x);
  gMinuit->GetParameter(1, nonBenC, errorParam);
  trackParam->SetNonBendingSlope(nonBenC);
  gMinuit->GetParameter(2, y, errorParam);
  trackParam->SetBendingCoor(y);
  gMinuit->GetParameter(3, benC, errorParam);
  trackParam->SetBendingSlope(benC);
  gMinuit->GetParameter(4, invBenP, errorParam);
  trackParam->SetInverseBendingMomentum(invBenP);
  
  // global result of the fit
  gMinuit->mnstat(fitFMin, fedm, errdef, npari, nparx, covStatus);
  track.SetFitFMin(fitFMin);
  
  // Get the covariance matrix if required
  if (calcCov) {
    // Covariance matrix according to HESSE status
    // If problem then keep only the diagonal terms (variances)
    Double_t matrix[5][5];
    gMinuit->mnemat(&matrix[0][0],5);
    if (covStatus == 3) trackParam->SetCovariances(matrix);
    else trackParam->SetVariances(matrix);
  } else trackParam->DeleteCovariances();
  
}

  //__________________________________________________________________________
void TrackChi2(Int_t & /*nParam*/, Double_t * /*gradient*/, Double_t &chi2, Double_t *param, Int_t /*flag*/)
{
  /// Return the "Chi2" to be minimized with Minuit for track fitting.
  /// Assumes that the track hits are sorted according to increasing Z.
  /// Track parameters at each TrackHit are updated accordingly.
  /// Vertex is used according to the flag "trackBeingFitted->GetFitWithVertex()".
  /// Multiple Coulomb scattering is taken into account according to the flag "trackBeingFitted->GetFitWithMCS()".
  
  AliMUONTrack *trackBeingFitted = (AliMUONTrack*) gMinuit->GetObjectFit();
  AliMUONTrackParam* trackParamAtHit = (AliMUONTrackParam*) trackBeingFitted->GetTrackParamAtHit()->First();
  Double_t dX, dY;
  chi2 = 0.; // initialize chi2
  
  // update track parameters
  trackParamAtHit->SetNonBendingCoor(param[0]);
  trackParamAtHit->SetNonBendingSlope(param[1]);
  trackParamAtHit->SetBendingCoor(param[2]);
  trackParamAtHit->SetBendingSlope(param[3]);
  trackParamAtHit->SetInverseBendingMomentum(param[4]);
  trackBeingFitted->UpdateTrackParamAtHit();
  
  // Take the vertex into account in the fit if required
  if (trackBeingFitted->GetFitWithVertex()) {
    AliMUONTrackParam paramAtVertex(*trackParamAtHit);
    AliMUONTrackExtrap::ExtrapToZ(&paramAtVertex, 0.);
    AliMUONHitForRec *vertex = trackBeingFitted->GetVertex();
    if (!vertex) {
      cout<<"Error in TrackChi2: Want to use the vertex in tracking but it has not been created!!"<<endl;
      exit(-1);
    }
    dX = vertex->GetNonBendingCoor() - paramAtVertex.GetNonBendingCoor();
    dY = vertex->GetBendingCoor() - paramAtVertex.GetBendingCoor();
    chi2 += dX * dX / vertex->GetNonBendingReso2() + dY * dY / vertex->GetBendingReso2();
  }
  
  // compute chi2 w/wo multiple scattering
  if (trackBeingFitted->GetFitWithMCS()) chi2 += trackBeingFitted->ComputeGlobalChi2(kTRUE);
  else chi2 += trackBeingFitted->ComputeGlobalChi2(kFALSE);
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::ImproveTracks()
{
  /// Improve tracks by removing clusters with local chi2 highter than the defined cut
  /// Recompute track parameters and covariances at the remaining clusters
  AliDebug(1,"Enter ImproveTracks");
  
  Double_t localChi2, worstLocalChi2;
  Int_t worstChamber;
  AliMUONTrackParam *trackParamAtHit, *worstTrackParamAtHit;
  
  // Remove double track to improve only "good" tracks
  RemoveDoubleTracks();
  
  AliMUONTrack *track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    
    while (!track->IsImproved()) {
      
      // Update track parameters and covariances
      track->UpdateCovTrackParamAtHit();
      
      // Compute local chi2 of each hits
      track->ComputeLocalChi2(kTRUE);
      
      // Look for the hit to remove
      worstTrackParamAtHit = 0;
      worstLocalChi2 = 0.;
      trackParamAtHit = (AliMUONTrackParam*) track->GetTrackParamAtHit()->First();
      while (trackParamAtHit) {
        
        // Pick up hit with the worst chi2
        localChi2 = trackParamAtHit->GetLocalChi2();
        if (localChi2 > worstLocalChi2) {
          worstLocalChi2 = localChi2;
          worstTrackParamAtHit = trackParamAtHit;
        }
        
      trackParamAtHit = (AliMUONTrackParam*) track->GetTrackParamAtHit()->After(trackParamAtHit);
      }
      
      // Check if bad hit found
      if (!worstTrackParamAtHit) {
        track->SetImproved(kTRUE);
        break;
      }
      
      // check whether the worst hit is removable or not
      if (!worstTrackParamAtHit->IsRemovable()) {
        track->SetImproved(kTRUE);
        break;
      }
      
      // Check whether the worst chi2 is under requirement or not
      if (worstLocalChi2 < 2. * fgkSigmaToCutForImprovement * fgkSigmaToCutForImprovement) { // 2 because 2 quantities in chi2
        track->SetImproved(kTRUE);
        break;
      }
      
      // Reset the second hit in the same station as the bad one as being NOT removable
      worstChamber = worstTrackParamAtHit->GetHitForRecPtr()->GetChamberNumber();
      if (worstChamber%2 == 0) ((AliMUONTrackParam*)track->GetTrackParamAtHit()->After(worstTrackParamAtHit))->SetRemovable(kFALSE);
      else ((AliMUONTrackParam*)track->GetTrackParamAtHit()->Before(worstTrackParamAtHit))->SetRemovable(kFALSE);
      
      // Remove the worst hit
      track->RemoveTrackParamAtHit(worstTrackParamAtHit);
      
      // Re-fit the track:
      // Take into account the multiple scattering
      // Calculate the track parameter covariance matrix
      track->SetFitWithVertex(kFALSE); // To be sure
      Fit(*track, kTRUE, kTRUE);
      
      // Printout for debuging
      if ((AliLog::GetDebugLevel("MUON","AliMUONTrackReconstructor") >= 1) || (AliLog::GetGlobalDebugLevel() >= 1)) {
        cout << "ImproveTracks: track " << fRecTracksPtr->IndexOf(track)+1 << " improved " << endl;
      }
      
    }
    
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
  }
  
}

  //__________________________________________________________________________
void AliMUONTrackReconstructor::Finalize()
{
  /// Fill AliMUONTrack's fHitForRecAtHit array
  /// Recompute track parameters and covariances at each attached cluster from those at the first one
  AliMUONTrack *track;
  AliMUONTrackParam *trackParamAtHit;
  
  track = (AliMUONTrack*) fRecTracksPtr->First();
  while (track) {
    
    // update track parameters if not already done
    if (!track->IsImproved()) track->UpdateCovTrackParamAtHit();
    
    trackParamAtHit = (AliMUONTrackParam*) (track->GetTrackParamAtHit()->First());
    while (trackParamAtHit) {
      
      // update array of track hit
      track->AddHitForRecAtHit(trackParamAtHit->GetHitForRecPtr());
      
      trackParamAtHit = (AliMUONTrackParam*) (track->GetTrackParamAtHit()->After(trackParamAtHit));
    }
    
    track = (AliMUONTrack*) fRecTracksPtr->After(track);
    
  }
    
}


