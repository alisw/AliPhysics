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

// ROOT includes
#include <TArrayI.h>
#include <TClonesArray.h>

// STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliCDBManager.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONRecoParam.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONVCluster.h"

#include "AliAnalysisTaskESDMCLabelAddition.h"

ClassImp(AliAnalysisTaskESDMCLabelAddition)

//----------------------------------------------------------------------
AliAnalysisTaskESDMCLabelAddition::AliAnalysisTaskESDMCLabelAddition():
AliAnalysisTaskSE(),
fDefaultStorage(""),
fRequestedStationMask(0),
fRequest2ChInSameSt45(kFALSE),
fExternalTrkSigmaCut(-1.),
fSigmaCut(-1.),
fExternalTrgSigmaCut(-1.),
fSigmaCutTrig(-1.)
{
  /// Default constructor
}


//----------------------------------------------------------------------
AliAnalysisTaskESDMCLabelAddition::AliAnalysisTaskESDMCLabelAddition(const char* name):
AliAnalysisTaskSE(name),
fDefaultStorage("raw://"),
fRequestedStationMask(0),
fRequest2ChInSameSt45(kFALSE),
fExternalTrkSigmaCut(-1.),
fSigmaCut(-1.),
fExternalTrgSigmaCut(-1.),
fSigmaCutTrig(-1.)
{
  /// Constructor
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::UserCreateOutputObjects()
{
  /// Create output objects
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::NotifyRun()
{
  /// Load OCDB inputs
  
  // load OCDB objects only once
  if (fSigmaCut > 0) return;
  
  // set OCDB location
  AliCDBManager* cdbm = AliCDBManager::Instance();
  if (cdbm->IsDefaultStorageSet()) printf("MCLabelAddition: CDB default storage already set!\n");
  else cdbm->SetDefaultStorage(fDefaultStorage.Data());
  if (cdbm->GetRun() > -1) printf("MCLabelAddition: run number already set!\n");
  else cdbm->SetRun(fCurrentRunNumber);
  
  // load recoParam
  const AliMUONRecoParam* recoParam = (AliMUONESDInterface::GetTracker())
  ? AliMUONESDInterface::GetTracker()->GetRecoParam()
  : AliMUONCDB::LoadRecoParam();
  
  if (!recoParam) {
    fRequestedStationMask = 0;
    fRequest2ChInSameSt45 = kFALSE;
    fSigmaCut = -1.;
    fSigmaCutTrig = -1.;
    return;
  }
  
  // compute the mask of requested stations from recoParam
  fRequestedStationMask = 0;
  for (Int_t i = 0; i < 5; i++) if (recoParam->RequestStation(i)) fRequestedStationMask |= ( 1 << i );
  
  // get from recoParam whether a track need 2 chambers hit in the same station (4 or 5) or not to be reconstructible
  fRequest2ChInSameSt45 = !recoParam->MakeMoreTrackCandidates();
  
  // get sigma cut to associate clusters with TrackRefs from recoParam if not already set manually
  if (fExternalTrkSigmaCut > 0) fSigmaCut = fExternalTrkSigmaCut;
  else if (recoParam->ImproveTracks()) fSigmaCut = recoParam->GetSigmaCutForImprovement();
  else fSigmaCut = recoParam->GetSigmaCutForTracking();
  
  // get sigma cut to associate trigger to triggerable track from recoParam if not already set manually
  if (fExternalTrgSigmaCut > 0) fSigmaCutTrig = fExternalTrgSigmaCut;
  else fSigmaCutTrig = recoParam->GetSigmaCutForTrigger();
  
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::UserExec(Option_t */*option*/)
{
  /// Execute analysis for current event				
  
  AliDebug(1, Form("MCLabel Addition: Analysing event # %5d\n",(Int_t) Entry())); 
  
  // make sure necessary information from OCDB have been loaded
  if (fSigmaCut < 0) return;
  
  /// Load ESD event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) {
    AliError("Cannot get input event");
    return;
  }      
  
  // Load MC event 
  AliMCEventHandler* mcH = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if ( ! mcH ) {
    AliError ("MCH event handler not found. Nothing done!");
    return;
  }
  
  
  // Get reference tracks
  AliMUONRecoCheck rc(esd,mcH);
  AliMUONVTrackStore* trackRefStore = rc.TrackRefs(-1);
  AliMUONVTriggerTrackStore* triggerTrackRefStore = rc.TriggerableTracks(-1);
  
  // Loop over reconstructed tracks
  AliESDMuonTrack *esdTrack = 0x0;
  Int_t nMuTracks = esd->GetNumberOfMuonTracks();
  for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack) {
    
    esdTrack = esd->GetMuonTrack(nMuTrack);
    
    // tracker tracks
    if (esdTrack->ContainTrackerData()) {
      
      // convert ESD track to MUON track (without recomputing track parameters at each clusters)
      AliMUONTrack muonTrack;
      AliMUONESDInterface::ESDToMUON(*esdTrack, muonTrack, kFALSE);
      
      // try to match, by position, the reconstructed track with a simulated one
      Int_t nMatchClustersByPosition = 0;
      AliMUONTrack* matchedTrackRefByPosition = rc.FindCompatibleTrack(muonTrack, *trackRefStore, nMatchClustersByPosition, kFALSE, fSigmaCut);
      Bool_t isMatchedYetByPosition = kFALSE;
      Bool_t isRecoDecayByPosition = kFALSE;
      Int_t decayLabelByPosition = -1, lastChDecayByPosition = 0;
      if (!matchedTrackRefByPosition || !matchedTrackRefByPosition->IsValid(fRequestedStationMask, fRequest2ChInSameSt45)) {
	decayLabelByPosition = IsDecayByPosition(muonTrack, *trackRefStore, isRecoDecayByPosition, lastChDecayByPosition);
	if (decayLabelByPosition >= 0) matchedTrackRefByPosition = 0x0;
	else if (matchedTrackRefByPosition) isMatchedYetByPosition = kTRUE;
      }
      Bool_t isFakeByPosition = (!matchedTrackRefByPosition && decayLabelByPosition < 0);
      
      // try to match, by using MC labels, the reconstructed track with a simulated one
      Int_t nMatchClustersByLabel = 0;
      AliMUONTrack* matchedTrackRefByLabel = rc.FindCompatibleTrack(muonTrack, *trackRefStore, nMatchClustersByLabel, kTRUE, fSigmaCut);
      Bool_t isMatchedYetByLabel = kFALSE;
      Bool_t isRecoDecayByLabel = kFALSE;
      Int_t decayLabelByLabel = -1, lastChDecayByLabel = 0;
      if (!matchedTrackRefByLabel || !matchedTrackRefByLabel->IsValid(fRequestedStationMask, fRequest2ChInSameSt45)) {
	decayLabelByLabel = IsDecayByLabel(muonTrack, isRecoDecayByLabel, lastChDecayByLabel);
	if (decayLabelByLabel >= 0) matchedTrackRefByLabel = 0x0;
	else if (matchedTrackRefByLabel) isMatchedYetByLabel = kTRUE;
      }
      Bool_t isFakeByLabel = (!matchedTrackRefByLabel && decayLabelByLabel < 0);
      
      // choose the best, or the only available, matched track
      AliMUONTrack* matchedTrackRef = 0x0;
      Bool_t isMatchedYet = kFALSE, isRecoDecay = kFALSE;
      Int_t decayLabel = -1;
      if (matchedTrackRefByPosition && matchedTrackRefByLabel && ((!isMatchedYetByPosition && !isMatchedYetByLabel) ||
								  (isMatchedYetByPosition && isMatchedYetByLabel))) {
	
	Int_t nMatchClusters = TMath::Max(nMatchClustersByPosition, nMatchClustersByLabel);
	matchedTrackRef = (nMatchClusters == nMatchClustersByPosition) ? matchedTrackRefByPosition : matchedTrackRefByLabel;
	isMatchedYet = isMatchedYetByPosition;
	
      } else if (matchedTrackRefByPosition && (!isMatchedYetByPosition || isFakeByLabel)) {
	
	matchedTrackRef = matchedTrackRefByPosition;
	isMatchedYet = isMatchedYetByPosition;
	
      } else if (matchedTrackRefByLabel && (!isMatchedYetByLabel || isFakeByPosition)) {
	
	matchedTrackRef = matchedTrackRefByLabel;
	isMatchedYet = isMatchedYetByLabel;
	
	// choose the best, or the only available, decay chain
      } else if (decayLabelByPosition >= 0 && decayLabelByLabel >= 0 && ((isRecoDecayByPosition && isRecoDecayByLabel) ||
									 (!isRecoDecayByPosition && !isRecoDecayByLabel))) {
	
	decayLabel = (lastChDecayByLabel > lastChDecayByPosition) ? decayLabelByLabel : decayLabelByPosition;
	isRecoDecay = isRecoDecayByPosition;
	
      } else if (decayLabelByPosition >= 0 && (isRecoDecayByPosition || decayLabelByLabel < 0)) {
	
	decayLabel = decayLabelByPosition;
	isRecoDecay = isRecoDecayByPosition;
	
      } else if (decayLabelByLabel >= 0) {
	
	decayLabel = decayLabelByLabel;
	isRecoDecay = isRecoDecayByLabel;
	
      }
      
      // set the MC label, the decay flag (bit 22) and the not-reconstructible flag (bit 23)
      if (matchedTrackRef) {
	
	esdTrack->SetLabel(matchedTrackRef->GetUniqueID());
	esdTrack->SetBit(BIT(22), kFALSE);
	esdTrack->SetBit(BIT(23), isMatchedYet);
	
      } else if (decayLabel >= 0) {
	
	esdTrack->SetLabel(decayLabel);
	esdTrack->SetBit(BIT(22), kTRUE);
	esdTrack->SetBit(BIT(23), !isRecoDecay);
	
      } else {
	
	esdTrack->SetLabel(-1);
	esdTrack->SetBit(BIT(22), kFALSE);
	esdTrack->SetBit(BIT(23), kFALSE);
	
      }
      
    } else { // ghosts
      
      // Convert ESD track to trigger track
      AliMUONLocalTrigger locTrg;
      AliMUONESDInterface::ESDToMUON(*esdTrack, locTrg);
      AliMUONTriggerTrack trigTrack;
      rc.TriggerToTrack(locTrg, trigTrack);
      
      // try to match the reconstructed track with a simulated one
      AliMUONTriggerTrack* matchedTrigTrackRef = rc.FindCompatibleTrack(trigTrack, *triggerTrackRefStore, fSigmaCutTrig);
      
      // set the MC label
      if (matchedTrigTrackRef) esdTrack->SetLabel(matchedTrigTrackRef->GetUniqueID());
      else esdTrack->SetLabel(-1);
      
    }
    
  }
  
}


//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::Terminate(Option_t */*option*/)
{
  /// Terminate analysis
  AliDebug(2, "Terminate()");
}

//----------------------------------------------------------------------
Int_t AliAnalysisTaskESDMCLabelAddition::IsDecay(Int_t nClusters, Int_t *chId, Int_t *labels,
						 Bool_t &isReconstructible, Int_t &lastCh) const
{
  /// Check whether this combination of clusters correspond to a decaying particle or not:
  /// More than 50% of clusters, including 1 before and 1 after the dipole, must be connected.
  /// - Return the MC label of the most downstream decay product or -1 if not a decay.
  /// - "isReconstructible" tells if the combination of matched clusters fulfil the reconstruction criteria.
  /// - As soon as we realized the decay chain cannot be tagged as reconstructible, we reject any chain ending
  ///   on a chamber equal to or upstream "lastCh" (used to select the best chain in case of multiple choices).
  /// - "lastCh" is reset the most downstream chamber of the found decay chain if any.
  
  Int_t halfCluster = nClusters/2;
  
  // loop over last clusters (if nClusters left < halfCluster the conditions cannot be fulfilled)
  Int_t firstLabel = -1, decayLabel = -1;
  isReconstructible = kFALSE;
  for (Int_t iCluster1 = nClusters-1; iCluster1 >= halfCluster; iCluster1--) {
    
    // if the last cluster is not on station 4 or 5 the conditions cannot be fulfilled
    if (chId[iCluster1] < 6) break;
    
    // skip clusters with no label or same label as at the begining of the previous step (already tested)
    if (labels[iCluster1] < 0 || labels[iCluster1] == firstLabel) continue;
    
    // is there any chance the hypothetical decay chain can be tagged reconstructible?
    Int_t stationId = chId[iCluster1]/2;
    Int_t stationMask = 1 << stationId;
    Int_t requestedStations = fRequestedStationMask >> stationId;
    Bool_t isValid = ((1 & requestedStations) == requestedStations);
    
    // if not: check whether we can find a better chain than already found
    if (!isValid && chId[iCluster1] <= lastCh) break;
    
    // count the number of fired chambers on stations 4 & 5
    Int_t nChHitInSt45[2] = {0, 0};
    nChHitInSt45[stationId-3] = 1;
    Int_t currentCh = chId[iCluster1];
    
    // get the ancestors
    TArrayI chainLabels(100);
    Int_t nParticles = 0;
    Int_t currentLabel = labels[iCluster1];
    do {
      chainLabels[nParticles++] = currentLabel;
      if (nParticles >= chainLabels.GetSize()) chainLabels.Set(2*chainLabels.GetSize());
      AliMCParticle* currentParticle = static_cast<AliMCParticle*>(fMCEvent->GetTrack(currentLabel));
      currentLabel = (currentParticle) ? currentParticle->GetMother() : -1;
    } while (currentLabel >= 0);
    
    // Loop over prior clusters
    firstLabel = labels[iCluster1];
    Int_t nCompatibleLabel = 1;
    Int_t currentParticle = 0;
    for (Int_t iCluster2 = iCluster1-1; iCluster2 >= 0; iCluster2--) {
      
      // if the number of clusters left is not enough the conditions cannot be fulfilled
      if (iCluster2 < halfCluster-nCompatibleLabel) break;
      
      if (labels[iCluster2] < 0) continue;
      
      // check if the cluster belong to the same particle or one of its ancestors
      Bool_t matchFound = kFALSE;
      for (Int_t iParticle = currentParticle; iParticle < nParticles; iParticle++) {
	if (labels[iCluster2] == chainLabels[iParticle]) {
	  currentParticle = iParticle;
	  matchFound = kTRUE;
	  break;
	}
      }
      if (matchFound) nCompatibleLabel++;
      else continue;
      
      // add this station to the mask
      stationId = chId[iCluster2]/2;
      stationMask |= 1 << stationId;
      
      // count the number of fired chamber on stations 4 & 5
      if (stationId > 2 && chId[iCluster2] < currentCh) {
	nChHitInSt45[stationId-3]++;
	currentCh = chId[iCluster2];
      }
      
      // check if we matched enough clusters to tag the track as a decay
      if (nCompatibleLabel <= halfCluster || chId[iCluster2] > 3 || chainLabels[currentParticle] == firstLabel) continue;
      
      // check if this chain is better than already found
      if (chId[iCluster1] > lastCh) {
	decayLabel = firstLabel;
	lastCh = chId[iCluster1];
      }
      
      // is there enough matched clusters on station 4 & 5 to make the track reconstructible?
      Bool_t isEnoughClOnSt45 = fRequest2ChInSameSt45 ? (nChHitInSt45[0] == 2 || nChHitInSt45[1] == 2)
      : (nChHitInSt45[0]+nChHitInSt45[1] >= 2);
      
      // is there any chance the current decay chain can still be tagged reconstructible?
      requestedStations = fRequestedStationMask >> stationId;
      isValid = (((stationMask >> stationId) & requestedStations) == requestedStations &&
		 (chId[iCluster2] > 5 || isEnoughClOnSt45));
      
      // if not then we cannot do better with this trial
      if (!isValid) break;
      
      // take in priority the decay chain that can be tagged reconstructible
      if (((stationMask & fRequestedStationMask) == fRequestedStationMask) && isEnoughClOnSt45) {
	lastCh = chId[iCluster1];
	isReconstructible = kTRUE;
	return firstLabel;
      }
      
    }
    
  }
  
  return decayLabel;
}

//----------------------------------------------------------------------
void AliAnalysisTaskESDMCLabelAddition::AddCompatibleClusters(const AliMUONTrack &track, const AliMUONTrack &trackRef,
							      TArrayI *labels, Int_t *nLabels) const
{
  /// Try to match clusters between track and trackRef and add the corresponding MC labels to the arrays
  
  Double_t chi2Max = 2. * fSigmaCut * fSigmaCut; // 2 because 2 quantities in chi2
  
  // Loop over clusters of first track
  Int_t nCl1 = track.GetNClusters();
  for(Int_t iCl1 = 0; iCl1 < nCl1; iCl1++) {
    AliMUONVCluster *cluster1 = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCl1))->GetClusterPtr();
    
    // Loop over clusters of second track
    Int_t nCl2 = trackRef.GetNClusters();
    for(Int_t iCl2 = 0; iCl2 < nCl2; iCl2++) {
      AliMUONVCluster *cluster2 = static_cast<AliMUONTrackParam*>(trackRef.GetTrackParamAtCluster()->UncheckedAt(iCl2))->GetClusterPtr();
      
      // check DE Id
      if (cluster1->GetDetElemId() != cluster2->GetDetElemId()) continue;
      
      // check local chi2
      Double_t dX = cluster1->GetX() - cluster2->GetX();
      Double_t dY = cluster1->GetY() - cluster2->GetY();
      Double_t chi2 = dX * dX / (cluster1->GetErrX2() + cluster2->GetErrX2()) + dY * dY / (cluster1->GetErrY2() + cluster2->GetErrY2());
      if (chi2 > chi2Max) continue;
      
      // expand array if needed
      if (nLabels[iCl1] >= labels[iCl1].GetSize()) labels[iCl1].Set(2*labels[iCl1].GetSize());
      
      // save label
      labels[iCl1][nLabels[iCl1]] = static_cast<Int_t>(trackRef.GetUniqueID());
      nLabels[iCl1]++;
      break;
      
    }
    
  }
  
}

//----------------------------------------------------------------------
Int_t AliAnalysisTaskESDMCLabelAddition::IsDecayByLabel(const AliMUONTrack &track, Bool_t &isReconstructible,
							Int_t &lastCh) const
{
  /// Check whether this track correspond to a decaying particle by using cluster MC labels.
  /// "lastCh" contains the chamber Id of the most downstream chamber hit by the decay chain
  
  Int_t nClusters = track.GetNClusters();
  Int_t *chId = new Int_t[nClusters];
  Int_t *labels = new Int_t[nClusters];
  
  // copy labels and chamber Ids
  for (Int_t iCluster = 0; iCluster < nClusters; iCluster++) {
    AliMUONVCluster* cluster = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCluster))->GetClusterPtr();
    chId[iCluster] = cluster->GetChamberId();
    labels[iCluster] = cluster->GetMCLabel();
  }
  
  // look for decay
  lastCh = 0;
  Int_t decayLabel = IsDecay(nClusters, chId, labels, isReconstructible, lastCh);
  
  delete[] chId;
  delete[] labels;
  
  return decayLabel;
}

//----------------------------------------------------------------------
Int_t AliAnalysisTaskESDMCLabelAddition::IsDecayByPosition(const AliMUONTrack &track, const AliMUONVTrackStore &trackRefStore,
							   Bool_t &isReconstructible, Int_t &lastCh) const
{
  /// Check whether this track correspond to a decaying particle by comparing clusters position
  /// All possible combinations of compatible clusters from every trackRefs are considered
  
  Int_t nClusters = track.GetNClusters();
  Int_t *chId = new Int_t[nClusters];
  Int_t *nLabels = new Int_t[nClusters];
  TArrayI *labels = new TArrayI[nClusters];
  
  // copy chamber Ids
  for (Int_t iCluster = 0; iCluster < nClusters; iCluster++) {
    AliMUONVCluster* cluster = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCluster))->GetClusterPtr();
    chId[iCluster] = cluster->GetChamberId();
    nLabels[iCluster] = 0;
    labels[iCluster].Set(100);
  }
  
  // loop over trackRef store and add label of compatible clusters
  TIter next1(trackRefStore.CreateIterator());
  AliMUONTrack* trackRef;
  while ( ( trackRef = static_cast<AliMUONTrack*>(next1()) ) )
    AddCompatibleClusters(track, *trackRef, labels, nLabels);
  
  // complete the arrays of labels with "-1" if no label was found for a given cluster
  for (Int_t iCluster = 0; iCluster < nClusters; iCluster++) {
    if (nLabels[iCluster] == 0) {
      labels[iCluster][0] = -1;
      nLabels[iCluster]++;
    }
  }
  
  // loop over all possible combinations
  Int_t *iLabel = new Int_t[nClusters];
  memset(iLabel,0,nClusters*sizeof(Int_t));
  iLabel[nClusters-1] = -1;
  Int_t *currentLabels = new Int_t[nClusters];
  Int_t decayLabel = -1;
  lastCh = 0;
  isReconstructible = kFALSE;
  while (kTRUE) {
    
    // go to the next combination
    Int_t iCl = nClusters-1;
    while (++iLabel[iCl] >= nLabels[iCl] && iCl > 0) iLabel[iCl--] = 0;
    if (iLabel[iCl] >= nLabels[iCl]) break; // no more combination
    
    // copy labels
    for (Int_t iCluster = 0; iCluster < nClusters; iCluster++)
      currentLabels[iCluster] = labels[iCluster][iLabel[iCluster]];
    
    // look for decay
    Int_t currentDecayLabel = IsDecay(nClusters, chId, currentLabels, isReconstructible, lastCh);
    if (currentDecayLabel >= 0) {
      decayLabel = currentDecayLabel;
      if (isReconstructible) break;
    }
    
  }
  
  delete[] chId;
  delete[] nLabels;
  delete[] labels;
  delete[] iLabel;
  delete[] currentLabels;
  
  return decayLabel;  
}

