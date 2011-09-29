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
/// \class AliMUONTrackHitPattern
///
/// This class propagates tracks to trigger chambers 
/// searching for matching trigger tracks and fired strips.
///
/// To each track, a hit pattern for trigger chambers is set.
/// 
/// The main method is:
/// * ExecuteValidation
///
///  \author Diego Stocco
//-----------------------------------------------------------------------------


#include "AliMUONTrackHitPattern.h"

#include "AliMUONConstants.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONTriggerUtilities.h"

#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"
#include "AliMpArea.h"
#include "AliMpConstants.h"

#include "AliLog.h"
#include "AliESDMuonTrack.h"
#include "AliCodeTimer.h"

#include <Riostream.h>
#include <TArrayI.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TArrayF.h>
#include <TRandom.h>

/// \cond CLASSIMP
ClassImp(AliMUONTrackHitPattern) // Class implementation in ROOT context
/// \endcond


//______________________________________________________________________________
AliMUONTrackHitPattern::AliMUONTrackHitPattern(const AliMUONRecoParam* recoParam,
                                               const AliMUONGeometryTransformer& transformer,
                                               const AliMUONVDigitStore& digitStore,
                                               const AliMUONTriggerUtilities* triggerUtilities)
: TObject(),
fkRecoParam(recoParam),
fkTransformer(transformer),
fkDigitStore(digitStore),
fkTriggerUtilities(triggerUtilities),
fkMaxDistance(99999.)
{
  /// Default constructor
  AliMUONTrackExtrap::SetField();
}


//______________________________________________________________________________
AliMUONTrackHitPattern::~AliMUONTrackHitPattern(void)
{
  /// Destructor
}


//______________________________________________________________________________
void AliMUONTrackHitPattern::ExecuteValidation(const AliMUONVTrackStore& trackStore,
					       const AliMUONVTriggerTrackStore& triggerTrackStore,
					       const AliMUONVTriggerStore& triggerStore) const
{
  //
  /// Main method:
  /// Loops on reco tracks, extrapolates them to trigger chambers
  /// and searches for matching trigger tracks and digits
  //

  // Get the hit pattern for all trigger tracks
  AliMUONTriggerTrack* triggerTrack;
  TIter itTriggerTrack(triggerTrackStore.CreateIterator());
  while ( ( triggerTrack = static_cast<AliMUONTriggerTrack*>(itTriggerTrack() ) ) ){
    UShort_t pattern = GetHitPattern(triggerTrack);
    triggerTrack->SetHitsPatternInTrigCh(pattern);
    AliDebug(1, Form("Hit pattern: hits 0x%x  slat %2i  board %3i  effFlag %i",
		     pattern & 0xFF, AliESDMuonTrack::GetSlatOrInfo(pattern),
		     triggerTrack->GetLoTrgNum(), AliESDMuonTrack::GetEffFlag(pattern)));
  }

  // Match tracker tracks with trigger tracks.
  TIter itTrack(trackStore.CreateIterator());
  AliMUONTrack* track;

  const Int_t kFirstTrigCh = AliMUONConstants::NTrackingCh();

  while ( ( track = static_cast<AliMUONTrack*>(itTrack()) ) )
  { 
    AliMUONTrackParam trackParam(*((AliMUONTrackParam*) (track->GetTrackParamAtCluster()->Last())));

    ApplyMCSCorrections(trackParam);
    AliMUONTrackExtrap::ExtrapToZCov(&trackParam, AliMUONConstants::DefaultChamberZ(kFirstTrigCh)); // extrap to 1st trigger chamber

    AliMUONTriggerTrack *matchedTriggerTrack = MatchTriggerTrack(track, trackParam, triggerTrackStore, triggerStore);

    // Copy trigger tracks hit pattern if there is matching,
    // otherwise calculate the hit pattern directly from tracker track:
    // the obtained pattern is good for check, but not good for efficiency determination.
    UShort_t pattern = matchedTriggerTrack ?
      matchedTriggerTrack->GetHitsPatternInTrigCh() : 
      GetHitPattern(&trackParam);

    track->SetHitsPatternInTrigCh(pattern);
  }
}


//______________________________________________________________________________
AliMUONTriggerTrack * 
AliMUONTrackHitPattern::MatchTriggerTrack(AliMUONTrack* track,
					  AliMUONTrackParam& trackParam,
					  const AliMUONVTriggerTrackStore& triggerTrackStore,
					  const AliMUONVTriggerStore& triggerStore) const
{
  //
  /// Match track with trigger track
  //

  Int_t matchTrigger = 0;
  Int_t loTrgNum(-1);
  TMatrixD paramDiff(3,1);
  Double_t chi2;
  Double_t chi2MatchTrigger = 0., minChi2MatchTrigger = 999.;
  Int_t doubleMatch = -1; // Check if track matches 2 trigger tracks
  Double_t doubleChi2 = -1.;
  AliMUONTriggerTrack* doubleTriggerTrack = 0x0;
  AliMUONTriggerTrack* matchedTriggerTrack = 0x0;
    

  // Covariance matrix 3x3 (X,Y,slopeY) for trigger tracks
  TMatrixD trackCov(3,3);

  AliMUONTriggerTrack *triggerTrack;
  TIter itTriggerTrack(triggerTrackStore.CreateIterator());
  while ( ( triggerTrack = static_cast<AliMUONTriggerTrack*>(itTriggerTrack() ) ) )
  {
    AliMUONTrackExtrap::LinearExtrapToZCov(&trackParam, triggerTrack->GetZ11());
    const TMatrixD& kParamCov = trackParam.GetCovariances();
    
    Double_t xTrack = trackParam.GetNonBendingCoor();
    Double_t yTrack = trackParam.GetBendingCoor();
    Double_t ySlopeTrack = trackParam.GetBendingSlope();

    paramDiff(0,0) = triggerTrack->GetX11() - xTrack;
    paramDiff(1,0) = triggerTrack->GetY11() - yTrack;
    paramDiff(2,0) = triggerTrack->GetSlopeY() - ySlopeTrack;

    // Covariance matrix 3x3 (X,Y,slopeY) for tracker tracks
    trackCov.Zero();
    trackCov(0,0) = kParamCov(0,0);
    trackCov(1,1) = kParamCov(2,2);
    trackCov(2,2) = kParamCov(3,3);
    trackCov(1,2) = kParamCov(2,3);
    trackCov(2,1) = kParamCov(3,2);

    // Covariance matrix 3x3 (X,Y,slopeY) for trigger tracks
    TMatrixD trigCov(triggerTrack->GetCovariances());

    TMatrixD sumCov(trackCov,TMatrixD::kPlus,trigCov);
    if (sumCov.Determinant() != 0) {
      sumCov.Invert();
      
      TMatrixD tmp(sumCov,TMatrixD::kMult,paramDiff);
      TMatrixD chi2M(paramDiff,TMatrixD::kTransposeMult,tmp);
      chi2 = chi2M(0,0);      
    } else {
      AliWarning(" Determinant = 0");
      Double_t sigma2 = 0.;
      chi2 = 0.;
      for (Int_t iVar = 0; iVar < 3; iVar++) {
	sigma2 = trackCov(iVar,iVar) + trigCov(iVar,iVar);
	chi2 += paramDiff(iVar,0) * paramDiff(iVar,0) / sigma2;
      }
    }

    chi2 /= 3.; // Normalized Chi2: 3 degrees of freedom (X,Y,slopeY)
    if (chi2 < GetRecoParam()->GetMaxNormChi2MatchTrigger()) 
    {
      Bool_t isDoubleTrack = (TMath::Abs(chi2 - minChi2MatchTrigger)<1.);
      if (chi2 < minChi2MatchTrigger && chi2 < GetRecoParam()->GetMaxNormChi2MatchTrigger()) 
      {
	if(isDoubleTrack)
	{
	  doubleMatch = loTrgNum;
	  doubleChi2 = chi2MatchTrigger;
	  doubleTriggerTrack = matchedTriggerTrack;
	}
	minChi2MatchTrigger = chi2;
	chi2MatchTrigger = chi2;
	loTrgNum = triggerTrack->GetLoTrgNum();
	matchedTriggerTrack = triggerTrack;
	AliMUONLocalTrigger* locTrg = triggerStore.FindLocal(loTrgNum);
	matchTrigger = 1;
	if(locTrg->LoLpt()>0) matchTrigger = 2;
	if(locTrg->LoHpt()>0) matchTrigger = 3;
      }
      else if(isDoubleTrack) 
      {
	doubleMatch = triggerTrack->GetLoTrgNum();
	doubleChi2 = chi2;
      }
    }
  }
  if(doubleMatch>=0)
  { // If two trigger tracks match, select the one passing more trigger cuts
    AliDebug(1, Form("Two candidates found: %i and %i",loTrgNum,doubleMatch));
    AliMUONLocalTrigger* locTrg1 = triggerStore.FindLocal(doubleMatch);
    if((locTrg1->LoLpt()>0 && matchTrigger<2) || (locTrg1->LoHpt() && matchTrigger<3))
    {
      if(locTrg1->LoHpt()>0) matchTrigger=3;
      else matchTrigger = 2;
      loTrgNum = doubleMatch;
      chi2MatchTrigger = doubleChi2;
      matchedTriggerTrack = doubleTriggerTrack;
    }
  }
    
  track->SetMatchTrigger(matchTrigger);
  track->SetChi2MatchTrigger(chi2MatchTrigger);

  AliMUONLocalTrigger* locTrg = static_cast<AliMUONLocalTrigger*>(triggerStore.FindLocal(loTrgNum));

  if (locTrg)
  {    
    track->SetLocalTrigger(locTrg->LoCircuit(),
			   locTrg->LoStripX(),
			   locTrg->LoStripY(),
			   locTrg->GetDeviation(),
			   locTrg->LoLpt(),
			   locTrg->LoHpt(),
			   locTrg->GetTriggerWithoutChamber());
  }

  return matchedTriggerTrack;
}


//______________________________________________________________________________
UShort_t AliMUONTrackHitPattern::GetHitPattern(const AliMUONTriggerTrack* matchedTriggerTrack) const
{
  //
  /// Get hit pattern on trigger chambers for the current trigger track
  //
  UShort_t pattern = 0;
  PerformTrigTrackMatch(pattern, matchedTriggerTrack);
  return pattern;
}


//______________________________________________________________________________
UShort_t AliMUONTrackHitPattern::GetHitPattern(AliMUONTrackParam* trackParam) const
{
  //
  /// Get hit pattern on trigger chambers for the current tracker track
  //
  UShort_t pattern = 0;
  Bool_t isMatch[2];
  const Int_t kNTrackingCh = AliMUONConstants::NTrackingCh();

  for(Int_t ch=0; ch<4; ++ch)
  {
    Int_t iChamber = kNTrackingCh+ch;
    AliMUONTrackExtrap::ExtrapToZCov(trackParam, AliMUONConstants::DefaultChamberZ(iChamber));
    FindPadMatchingTrack(*trackParam, isMatch, iChamber);
    for(Int_t cath=0; cath<2; ++cath)
    {
      if(isMatch[cath]) AliESDMuonTrack::SetFiredChamber(pattern, cath, ch);
    }
  }

  // pattern obtained by propagation of tracker track
  // when it does not match the trigger.
  AliESDMuonTrack::AddEffInfo(pattern, AliESDMuonTrack::kTrackerTrackPattern);

  return pattern;
}

//______________________________________________________________________________
void 
AliMUONTrackHitPattern::ApplyMCSCorrections(AliMUONTrackParam& trackParam) const
{
  //
  /// Returns uncertainties on extrapolated position.
  /// Takes into account Branson plane corrections in the iron wall.
  //

  const Float_t kZFilterOut = AliMUONConstants::MuonFilterZEnd();
  const Float_t kFilterThickness = kZFilterOut-AliMUONConstants::MuonFilterZBeg(); // cm

  AliMUONTrackExtrap::ExtrapToZCov(&trackParam, kZFilterOut); // Extrap to muon filter end
  AliMUONTrackExtrap::AddMCSEffect(&trackParam, kFilterThickness, AliMUONConstants::MuonFilterX0()); // Add MCS effects
  return;
}


//______________________________________________________________________________
void 
AliMUONTrackHitPattern::FindPadMatchingTrack(const AliMUONTrackParam& trackParam,
                                             Bool_t isMatch[2], Int_t iChamber) const
{
    //
    /// Given the tracker track position, searches for matching digits.
    //

    Float_t minMatchDist[2];

    for(Int_t cath=0; cath<2; ++cath)
    {
      isMatch[cath]=kFALSE;
      minMatchDist[cath]=fkMaxDistance/10.;
    }

    TIter next(fkDigitStore.CreateTriggerIterator());
    AliMUONVDigit* mDigit;

    while ( ( mDigit = static_cast<AliMUONVDigit*>(next()) ) )
    {
      Int_t currDetElemId = mDigit->DetElemId();
      Int_t currCh = AliMpDEManager::GetChamberId(currDetElemId);
      if(currCh!=iChamber) continue;
      Int_t cathode = mDigit->Cathode();
      Int_t ix = mDigit->PadX();
      Int_t iy = mDigit->PadY();
      Float_t xpad, ypad, zpad;
      const AliMpVSegmentation* seg = AliMpSegmentation::Instance()
        ->GetMpSegmentation(currDetElemId,AliMp::GetCathodType(cathode));
      
      AliMpPad pad = seg->PadByIndices(ix,iy,kTRUE);
      Float_t xlocal1 = pad.GetPositionX();
      Float_t ylocal1 = pad.GetPositionY();
      Float_t dpx = pad.GetDimensionX();
      Float_t dpy = pad.GetDimensionY();

      fkTransformer.Local2Global(currDetElemId, xlocal1, ylocal1, 0, xpad, ypad, zpad);
      Float_t matchDist = MinDistanceFromPad(xpad, ypad, zpad, dpx, dpy, trackParam);
      if(matchDist>minMatchDist[cathode])continue;
      isMatch[cathode] = kTRUE;
      if(isMatch[0] && isMatch[1]) break;
      minMatchDist[cathode] = matchDist;
    }
}


//______________________________________________________________________________
Float_t 
AliMUONTrackHitPattern::MinDistanceFromPad(Float_t xPad, Float_t yPad, Float_t zPad,
                                           Float_t dpx, Float_t dpy, 
                                           const AliMUONTrackParam& trackParam) const
{
    //
    /// Decides if the digit belongs to the tracker track.
    //

    AliMUONTrackParam trackParamAtPadZ(trackParam);
    AliMUONTrackExtrap::LinearExtrapToZCov(&trackParamAtPadZ, zPad);

    Float_t xTrackAtPad = trackParamAtPadZ.GetNonBendingCoor();
    Float_t yTrackAtPad = trackParamAtPadZ.GetBendingCoor();

    const Float_t kNSigma = GetRecoParam()->GetSigmaCutForTrigger();

    const TMatrixD& kCovParam = trackParamAtPadZ.GetCovariances();
    
    Float_t sigmaX = TMath::Sqrt(kCovParam(0,0));
    Float_t sigmaY = TMath::Sqrt(kCovParam(2,2));

    Float_t maxDistX = kNSigma * sigmaX; // in cm
    Float_t maxDistY = kNSigma * sigmaY; // in cm

    Float_t deltaX = TMath::Abs(xPad-xTrackAtPad)-dpx;
    Float_t deltaY = TMath::Abs(yPad-yTrackAtPad)-dpy;

    Float_t matchDist = fkMaxDistance;
    if(deltaX<=maxDistX && deltaY<=maxDistY) matchDist = TMath::Max(deltaX, deltaY);

    return matchDist;
}


//_____________________________________________________________________________
Bool_t AliMUONTrackHitPattern::FindPadMatchingTrig(const TVector3& vec11, const TVector3& vec21,
                                                  Int_t matchedDetElemId[2], TObjArray& pads) const
{
  //
  /// Check slat and board number of digit matching trigger track
  //
  
  enum {kBending, kNonBending};
  
  Float_t minMatchDist[2];
  Int_t padsInCheckArea[2];
  
  Int_t inputDetElemId = matchedDetElemId[0];
  
  for(Int_t cath=0; cath<2; cath++){
    minMatchDist[cath] = fkMaxDistance/10.;
    padsInCheckArea[cath] = 0;
  }
  
  Int_t iChamber = AliMpDEManager::GetChamberId(inputDetElemId);  
  Int_t iSlat = inputDetElemId%100;
  
  TIter next(fkDigitStore.CreateTriggerIterator());
  AliMUONVDigit* mDigit;
  TVector3 localExtrap;
  
  Int_t previousDetElemId = -1;
  
  while ( ( mDigit = static_cast<AliMUONVDigit*>(next()) ) )
  {
    Int_t currDetElemId = mDigit->DetElemId();
    Int_t currCh = AliMpDEManager::GetChamberId(currDetElemId);
    if ( currCh != iChamber ) continue;
    Int_t currSlat = currDetElemId%100;
    Int_t slatDiff = TMath::Abs(currSlat-iSlat);
    if ( slatDiff>1 && slatDiff<17 ) continue; // Check neighbour slats
    Int_t cathode = mDigit->Cathode();
    Int_t ix = mDigit->PadX();
    Int_t iy = mDigit->PadY();
    const AliMpVSegmentation* seg = AliMpSegmentation::Instance()
    ->GetMpSegmentation(currDetElemId,AliMp::GetCathodType(cathode));
    
    AliMpPad pad = seg->PadByIndices(ix,iy,kTRUE);
    
    if ( currDetElemId != previousDetElemId ) {
      PosInDetElemIdLocal(localExtrap, vec11, vec21, currDetElemId);
      previousDetElemId = currDetElemId;
    }
    
    AliDebug(2, Form("\nDetElemId = %i  Cathode = %i  Pad = (%i,%i) = (%.2f,%.2f)  Dim = (%.2f,%.2f)  Track = (%.2f,%.2f)\n",
                     currDetElemId,cathode,ix,iy,pad.GetPositionX(),pad.GetPositionY(),pad.GetDimensionX(),pad.GetDimensionY(),localExtrap.X(),localExtrap.Y()));
    Float_t matchDist = PadMatchTrack(pad, localExtrap);
    if ( matchDist < fkMaxDistance/2. ) padsInCheckArea[cathode]++;
    if ( matchDist > minMatchDist[cathode] ) continue;
    if ( pads.At(cathode) ) delete pads.RemoveAt(cathode);
    pads.AddAt((AliMpPad*)pad.Clone(),cathode);
    minMatchDist[cathode] = matchDist;
    matchedDetElemId[cathode] = currDetElemId; // Set the input detection element id to the matched one
  } // loop on digits
  
  // If track matches many pads, it is not good for effciency determination.
  // However we still want to calculate the hit pattern.
  for ( Int_t cath=0; cath<2; cath++ ){
    if ( padsInCheckArea[cath] > 2 ) {
      AliDebug(1, Form("padsInCheckArea[%i] = %i\n",cath,padsInCheckArea[cath]));
      return kFALSE;
    }
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
Float_t AliMUONTrackHitPattern::PadMatchTrack(const AliMpPad& pad, const TVector3& trackPosAtPad) const
{
  //
  /// Decides if the digit belongs to the trigger track.
  //
  
  Float_t xPad = pad.GetPositionX();
  Float_t yPad = pad.GetPositionY();
  Float_t dpx = pad.GetDimensionX();
  Float_t dpy = pad.GetDimensionY();  
  

  Float_t maxDist = GetRecoParam()->GetStripCutForTrigger() * 2. * TMath::Min(dpx,dpy); // cm
  if ( maxDist<2. ) maxDist = 2.;
  Float_t maxDistCheckArea = GetRecoParam()->GetMaxStripAreaForTrigger() * 2. *  TMath::Min(dpx,dpy); // cm

  Float_t matchDist = fkMaxDistance;

  Float_t deltaX = TMath::Abs(xPad-trackPosAtPad.X())-dpx;
  Float_t deltaY = TMath::Abs(yPad-trackPosAtPad.Y())-dpy;
  Float_t maxDistX = maxDist;
  Float_t maxDistY = maxDist;

  if(deltaX<=maxDistX && deltaY<=maxDistY) matchDist = TMath::Max(deltaX, deltaY);
  else if(deltaX<=maxDistCheckArea && deltaY<=maxDistCheckArea) matchDist = fkMaxDistance/5.;
  return matchDist;
}


//_____________________________________________________________________________
Int_t AliMUONTrackHitPattern::DetElemIdFromPos(Float_t x, Float_t y, 
                                               Int_t chamber, Int_t foundDetElemId[2]) const
{
  //
  /// Given the (x,y) position in the chamber,
  /// it returns the corresponding slats
  /// Caveat: at the border the result is not univoque
  //
  
  Int_t nFound = 0;
  
  foundDetElemId[0] = foundDetElemId[1] = 0;
  AliMpArea pointArea(x, y, 2.*AliMpConstants::LengthTolerance(), 2.*AliMpConstants::LengthTolerance());
  AliMpDEIterator it;
  for ( it.First(chamber-1); ! it.IsDone(); it.Next() ){
    Int_t detElemId = it.CurrentDEId();
    AliMpArea* deArea = fkTransformer.GetDEArea(detElemId);
    
    if ( deArea->Contains(pointArea) ) {
      foundDetElemId[nFound++] = detElemId;
      if ( nFound == 2 ) break;
    }
    else if ( nFound > 0 ) break;
  } // loop on detElemId
  
  return nFound;
}



//_____________________________________________________________________________
Bool_t AliMUONTrackHitPattern::PadsFromPos(const TVector3& vec11, const TVector3& vec21,
                                           Int_t detElemId, TObjArray& pads) const
{
  //
  /// Given the (x,y) position in the chamber,
  /// it returns the corresponding local board
  //
  
  pads.Delete();
  
  TVector3 localCoor;
  PosInDetElemIdLocal(localCoor, vec11, vec21, detElemId);
  for ( Int_t icath=0; icath<2; icath++ ) {
    const AliMpVSegmentation* seg = 
    AliMpSegmentation::Instance()
    ->GetMpSegmentation(detElemId,AliMp::GetCathodType(icath));
    AliMpPad pad = seg->PadByPosition(localCoor.X(),localCoor.Y(),kFALSE);
    if ( pad.IsValid() ) {
      pads.AddAt(pad.Clone(), icath);
    }
  }
  
  return pads.GetEntries();
}


//_____________________________________________________________________________
Bool_t AliMUONTrackHitPattern::PosInDetElemIdLocal(TVector3& localCoor, const TVector3& globalPoint1,
                                                   const TVector3& globalPoint2, Int_t detElemId) const
{
  /// Given two points belonging to a line (global coordinates)
  /// it returns the intersection point with the detElemId (local coordinates)
  
  Double_t xloc, yloc, zloc;
  fkTransformer.Global2Local(detElemId, globalPoint1.X(), globalPoint1.Y(), globalPoint1.Z(), xloc, yloc, zloc);
  TVector3 localPoint1(xloc, yloc, zloc);
  fkTransformer.Global2Local(detElemId, globalPoint2.X(), globalPoint2.Y(), globalPoint2.Z(), xloc, yloc, zloc);
  TVector3 localPoint2(xloc, yloc, zloc);
  localCoor = localPoint1 - ( localPoint1.Z() / ( localPoint2.Z() - localPoint1.Z() ) ) * ( localPoint2 - localPoint1 );
  
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliMUONTrackHitPattern::IsCloseToAccEdge(TObjArray& pads, Int_t detElemId, Float_t coor[2]) const
{
  AliMpArea* deArea = fkTransformer.GetDEArea(detElemId);
  Float_t dpx = ((AliMpPad*)pads.At(1))->GetDimensionX();
  Float_t dpy = ((AliMpPad*)pads.At(0))->GetDimensionY();
  Float_t resolution[2] = { 0.75*dpx, 0.75*dpy };
  Float_t sign[3] = {0., 1., -1.};
  for ( Int_t ineighx=0; ineighx<3; ineighx++ ) {
    for ( Int_t ineighy=0; ineighy<3; ineighy++ ) {
      AliMpArea pointArea(coor[0]+sign[ineighx]*resolution[0],
                          coor[1]+sign[ineighy]*resolution[1],
                          2.*AliMpConstants::LengthTolerance(),
                          2.*AliMpConstants::LengthTolerance());
      if ( ! deArea->Contains(pointArea) ) return kTRUE;
    } // loop on y neighbours
  } // loop on x neighbours
  return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliMUONTrackHitPattern::IsMasked(const AliMpPad& pad, Int_t detElemId, Int_t cathode, const TVector3& vec11, const TVector3& vec21) const
{
  //
  /// Check if pad or its neighbours are masked
  //
  
  Int_t nMasked = 0;
  
  if ( fkTriggerUtilities->IsMasked(pad, detElemId, cathode) ) ++nMasked;
  
  // Check closest neighbour
  
  TVector3 localCoor;
  PosInDetElemIdLocal(localCoor, vec11, vec21, detElemId);
  
  Float_t padPos[2] = { pad.GetPositionX(), pad.GetPositionY()};
  Float_t padDim[2] = { pad.GetDimensionX(), pad.GetDimensionY()};
  Float_t inpactPos[2] = { localCoor.X(), localCoor.Y()};
  Float_t sign[2] = {-1., 1.};
  Int_t icoor = 1-cathode;
  Int_t addSlatSign[3] = {0,1,-1};
  Int_t inputCh = AliMpDEManager::GetChamberId(detElemId)+1;
  Int_t inputSlat = detElemId%100;
  
  Float_t newPos[2];
  for ( Int_t ineigh=0; ineigh<2; ineigh++ ) {
    newPos[1-icoor] = inpactPos[1-icoor];
    newPos[icoor] = inpactPos[icoor] + 1.05 * sign[ineigh] * padDim[icoor];
    if ( TMath::Abs(newPos[icoor] - padPos[icoor]) < padDim[icoor] ) continue;
    const AliMpVSegmentation* seg = 
      AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));
    AliMpPad neighPad = seg->PadByPosition(newPos[0],newPos[1],kFALSE);
    if ( neighPad.IsValid() ) {
      if ( fkTriggerUtilities->IsMasked(neighPad, detElemId, cathode) ) ++nMasked;
    }
    else {
      TVector3 deltaVec(newPos[0]-inpactPos[0],newPos[1]-inpactPos[1],0.);
      TVector3 transVec11 = vec11+deltaVec;
      TVector3 transVec21 = vec21+deltaVec;
      TObjArray padsFromPos;
      padsFromPos.SetOwner();
      for ( Int_t iAddSlat=0; iAddSlat<2; iAddSlat++ ) {
        Int_t currSlat = (inputSlat + addSlatSign[iAddSlat])%18;
        Int_t currDetElemId = 100 * inputCh + currSlat;
        PadsFromPos(transVec11,transVec21,currDetElemId,padsFromPos);
        AliMpPad* currPad = (AliMpPad*)padsFromPos.UncheckedAt(cathode);
        Bool_t isMasked = ( currPad ) ? fkTriggerUtilities->IsMasked(*currPad, currDetElemId, cathode) : kFALSE;
        padsFromPos.Delete();
        if ( isMasked ) ++nMasked;
      }
    } // loop on slats
  } // loop on neigbours
  
  Double_t maskedProb = ((Double_t)nMasked)/3.;
  if ( gRandom->Rndm() < maskedProb ) return kTRUE;
  
  return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliMUONTrackHitPattern::PerformTrigTrackMatch(UShort_t &pattern,
						     const AliMUONTriggerTrack* matchedTrigTrack) const
{
  //
  /// It searches for matching digits around the trigger track.
  //
  
  AliCodeTimerAuto("",0);

  enum {kBending, kNonBending};

  TArrayF zMeanChamber(AliMUONConstants::NTriggerCh());
  zMeanChamber[0] = matchedTrigTrack->GetZ11();
  zMeanChamber[1] = matchedTrigTrack->GetZ11() + AliMUONConstants::DefaultChamberZ(11) - AliMUONConstants::DefaultChamberZ(10);
  zMeanChamber[2] = matchedTrigTrack->GetZ21();
  zMeanChamber[3] = matchedTrigTrack->GetZ21() + AliMUONConstants::DefaultChamberZ(13) - AliMUONConstants::DefaultChamberZ(12);

  TArrayI digitPerTrack(2);
  digitPerTrack.Reset();

  Float_t trackIntersectCh[2];

  Float_t slopeX = matchedTrigTrack->GetSlopeX();
  Float_t slopeY = matchedTrigTrack->GetSlopeY();
  
  Float_t z11 = matchedTrigTrack->GetZ11();
  Float_t x11 = slopeX * z11;
  Float_t y11 = matchedTrigTrack->GetY11();
  Float_t z21 =  matchedTrigTrack->GetZ21();
  Float_t x21 = slopeX * z21;
  Float_t y21 = y11 + slopeY * (z21-z11);
  TVector3 vec11(x11, y11, z11), vec21(x21, y21, z21);
  
  
  Int_t firstSlat = -1, firstBoard = -1;
  AliESDMuonTrack::EAliTriggerChPatternFlag goodForEff = AliESDMuonTrack::kBoardEff;
  TObjArray matchedPads(2), padsFromPos(2), validPads(2);
  matchedPads.SetOwner(); padsFromPos.SetOwner();
  Int_t matchedDetElemId[2], detElemIdFromTrack[2];
  
  for(Int_t ich=0; ich<AliMUONConstants::NTriggerCh(); ich++) { // chamber loop
    
    // searching track intersection with chambers (first approximation)
    Float_t deltaZ = zMeanChamber[ich] - zMeanChamber[0];
    trackIntersectCh[0] = zMeanChamber[ich] * slopeX;
    trackIntersectCh[1] = y11 + deltaZ * slopeY;
    Int_t nFound = DetElemIdFromPos(trackIntersectCh[0], trackIntersectCh[1], 11+ich, detElemIdFromTrack);
    if ( nFound == 0 ) {
      // track is rejected since the extrapolated track
      // does not match a slat (border effects)
      AliESDMuonTrack::AddEffInfo(pattern, AliESDMuonTrack::kTrackOutsideGeometry);
      goodForEff = AliESDMuonTrack::kNoEff;
      AliDebug(1, "Warning: trigger track outside trigger chamber\n");
      continue;
    }
    
    matchedDetElemId[0] = matchedDetElemId[1] = detElemIdFromTrack[0];
    
    if ( ! FindPadMatchingTrig(vec11, vec21, matchedDetElemId, matchedPads) ) {
      // if ! FindPadMatchingTrig => too many digits matching pad =>
      //                          => Event not clear => Do not use for efficiency calculation
      AliESDMuonTrack::AddEffInfo(pattern, AliESDMuonTrack::kTrackMatchesManyPads);
      goodForEff = AliESDMuonTrack::kNoEff;
      AliDebug(1, Form("Warning: track = %p (%i) matches many pads. Rejected!\n",(void *)matchedTrigTrack, matchedDetElemId[0]));
    }
    
    Int_t nMatched = 0;

    Int_t mostProbDEmatched = detElemIdFromTrack[0];
    for ( Int_t icath=0; icath<2; icath++ ) {
      if ( matchedPads.UncheckedAt(icath) ) {
        nMatched++;
        // Fill pattern anyway
        AliESDMuonTrack::SetFiredChamber(pattern, icath, ich);
        digitPerTrack[icath]++;
        mostProbDEmatched = matchedDetElemId[icath];
      }
    }
    Int_t mostProbDEfromTrack = detElemIdFromTrack[0];
    for ( Int_t ifound=0; ifound<nFound; ifound++ ) {
      if ( detElemIdFromTrack[ifound] == mostProbDEmatched ) {
        mostProbDEfromTrack = mostProbDEmatched;
        break;
      }
    }
    
    if ( goodForEff == AliESDMuonTrack::kNoEff ) continue;

    if ( nMatched < 2 ) PadsFromPos(vec11, vec21, mostProbDEfromTrack, padsFromPos);

    for ( Int_t cath=0; cath<2; cath++ ) {
      if ( matchedPads.UncheckedAt(cath) ) validPads.AddAt(matchedPads.UncheckedAt(cath),cath);
      else if ( padsFromPos.UncheckedAt(cath) ) {
        AliMpPad* currPad = (AliMpPad*)padsFromPos.UncheckedAt(cath);
        validPads.AddAt(currPad,cath);
        if ( IsMasked(*currPad, mostProbDEfromTrack, cath, vec11, vec21) ) {
          // Check if strip was masked (if inefficient strip is found)
          AliESDMuonTrack::AddEffInfo(pattern,25,AliESDMuonTrack::kNoEff); // pad is masked
          AliDebug(1,Form("DetElemId %i  cath %i  strip %i is masked: effFlag 0", mostProbDEfromTrack, cath, currPad->GetLocalBoardId(0)));
          goodForEff = AliESDMuonTrack::kNoEff;
        }
      }
      else goodForEff = AliESDMuonTrack::kNoEff;
    } // loop on cathodes
        
    if ( goodForEff != AliESDMuonTrack::kNoEff ) {
      
      if ( nMatched == 0 && IsCloseToAccEdge(padsFromPos, mostProbDEfromTrack, trackIntersectCh) ) {
        // Non of cathodes is fired.
        // If we are close to the edge of the RPC
        // it could be a problem of acceptance 
        AliESDMuonTrack::AddEffInfo(pattern, AliESDMuonTrack::kTrackOutsideGeometry);
        goodForEff = AliESDMuonTrack::kNoEff;
        AliDebug(1, "Warning: trigger track at the edge of the chamber\n");
      }
      
      Int_t currSlat = mostProbDEmatched%100;
      if ( firstSlat < 0 ) firstSlat = currSlat;
      else if ( currSlat != firstSlat ) {
        goodForEff = AliESDMuonTrack::kChEff;
        firstSlat = AliESDMuonTrack::kCrossDifferentSlats;
      }
    
      if ( firstBoard < 0 ) firstBoard = ((AliMpPad*)validPads[kBending])->GetLocalBoardId(0);
        
      for ( Int_t cath=0; cath<2; cath++ ){      
      
        if ( goodForEff == AliESDMuonTrack::kBoardEff) {
          Bool_t atLeastOneLoc = kFALSE;
          AliMpPad* currPad = (AliMpPad*)validPads.UncheckedAt(cath);
          for ( Int_t iloc=0; iloc<currPad->GetNofLocations(); iloc++) {
            if ( currPad->GetLocalBoardId(iloc) == firstBoard ) {
              atLeastOneLoc = kTRUE;
              break;
            }
          } // loop on locations
          if ( ! atLeastOneLoc ) goodForEff = AliESDMuonTrack::kSlatEff;
        }
      } // loop on cathodes
    } // if track good for efficiency
    matchedPads.Delete();
    padsFromPos.Delete();
  } // end chamber loop
  
  if ( goodForEff == AliESDMuonTrack::kNoEff ) return kFALSE;

  for(Int_t cath=0; cath<2; cath++){
    if(digitPerTrack[cath]<3) {
      // track is rejected since the number of associated
      // digits found is less than 3.
      AliESDMuonTrack::AddEffInfo(pattern, AliESDMuonTrack::kTrackMatchesFewPads);
      goodForEff = AliESDMuonTrack::kNoEff;
      AliDebug(1, Form("Warning: found %i digits for trigger track cathode %i.\nRejecting event\n", digitPerTrack[cath],cath));
    }
  } // loop on cathodes 

  if ( goodForEff == AliESDMuonTrack::kNoEff ) return kFALSE;
  
  AliESDMuonTrack::AddEffInfo(pattern, firstSlat, goodForEff);
  return kTRUE;
}
