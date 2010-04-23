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
#include "AliMUONDigitMaker.h"
#include "AliMUONDigitStoreV1.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackParam.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"
#include "AliMUONReconstructor.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONVTriggerTrackStore.h"

#include "AliMpConstants.h"

#include "AliLog.h"
#include "AliTracker.h"
#include "AliESDMuonTrack.h"

#include <Riostream.h>
#include <TArrayS.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TSystem.h>
#include <TVector2.h>

#include <cassert>

/// \cond CLASSIMP
ClassImp(AliMUONTrackHitPattern) // Class implementation in ROOT context
/// \endcond


//______________________________________________________________________________
AliMUONTrackHitPattern::AliMUONTrackHitPattern(const AliMUONRecoParam* recoParam,
                                               const AliMUONGeometryTransformer& transformer,
                                               const AliMUONDigitMaker& digitMaker)
: TObject(),
fkRecoParam(recoParam),
fkTransformer(transformer),
fkDigitMaker(digitMaker),
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

//_____________________________________________________________________________
void AliMUONTrackHitPattern::CheckConstants() const
{
/// Check consistence of redefined constants 

  assert(fgkNcathodes == AliMpConstants::NofCathodes());    
  assert(fgkNchambers == AliMpConstants::NofTriggerChambers());    
  assert(fgkNplanes == AliMpConstants::NofTriggerChambers() * fgkNcathodes);    
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

  AliMUONDigitStoreV1 digitStore;
  fkDigitMaker.TriggerToDigitsStore(triggerStore,digitStore);

  // Get the hit pattern for all trigger tracks
  AliMUONTriggerTrack* triggerTrack;
  TIter itTriggerTrack(triggerTrackStore.CreateIterator());
  while ( ( triggerTrack = static_cast<AliMUONTriggerTrack*>(itTriggerTrack() ) ) ){
    UShort_t pattern = GetHitPattern(triggerTrack, digitStore);
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
      GetHitPattern(&trackParam, digitStore);

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
UShort_t AliMUONTrackHitPattern::GetHitPattern(AliMUONTriggerTrack* matchedTriggerTrack,
					       AliMUONVDigitStore& digitStore) const
{
  //
  /// Get hit pattern on trigger chambers for the current trigger track
  //
  UShort_t pattern = 0;
  PerformTrigTrackMatch(pattern, matchedTriggerTrack, digitStore);
  return pattern;
}


//______________________________________________________________________________
UShort_t AliMUONTrackHitPattern::GetHitPattern(AliMUONTrackParam* trackParam,
					       AliMUONVDigitStore& digitStore) const
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
    FindPadMatchingTrack(digitStore, *trackParam, isMatch, iChamber);
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
AliMUONTrackHitPattern::FindPadMatchingTrack(const AliMUONVDigitStore& digitStore,
                                             const AliMUONTrackParam& trackParam,
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

    TIter next(digitStore.CreateIterator());
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
Int_t AliMUONTrackHitPattern::FindPadMatchingTrig(const AliMUONVDigitStore& digitStore, Int_t &detElemId,
						  Float_t coor[2], Bool_t isMatch[2],
						  TArrayI nboard[2], TArrayF &zRealMatch, Float_t y11) const
{
    //
    /// Check slat and board number of digit matching trigger track
    //

    enum {kBending, kNonBending};

    Float_t minMatchDist[fgkNcathodes];
    Int_t padsInCheckArea[fgkNcathodes];

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	isMatch[cath] = kFALSE;
	minMatchDist[cath] = fkMaxDistance/10.;
	padsInCheckArea[cath] = 0;
    }
    Int_t iChamber = AliMpDEManager::GetChamberId(detElemId);
    Int_t ch = iChamber-10;
    Float_t oldDeltaZ = AliMUONConstants::DefaultChamberZ(iChamber) - AliMUONConstants::DefaultChamberZ(10);
    Float_t y = coor[1];
    Int_t iSlat = detElemId%100;
    Int_t trigDigitBendPlane = -1;
    Int_t foundDetElemId = detElemId;
    Float_t foundZmatch=999.;
    Float_t yCoorAtPadZ=999.;

    TIter next(digitStore.CreateIterator());
    AliMUONVDigit* mDigit;
    Int_t idigit=0;
    
    while ( ( mDigit = static_cast<AliMUONVDigit*>(next()) ) )
    {
	idigit++;
	Int_t currDetElemId = mDigit->DetElemId();
	Int_t currCh = AliMpDEManager::GetChamberId(currDetElemId);
	if(currCh!=iChamber) continue;
	Int_t currSlat = currDetElemId%100;
	Int_t slatDiff = TMath::Abs(currSlat-iSlat);
	if(slatDiff>1 && slatDiff<17) continue; // Check neighbour slats
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
	AliDebug(2, Form("\nDetElemId = %i  Cathode = %i  Pad = (%i,%i) = (%.2f,%.2f)  Dim = (%.2f,%.2f)  Track = (%.2f,%.2f)\n",currDetElemId,cathode,ix,iy,xpad,ypad,dpx,dpy,coor[0],coor[1]));
	// searching track intersection with chambers (second approximation)
	if(ch%2==1){
	    //if(iChamber%2==1){
	    Float_t deltaZ = zpad - zRealMatch[0];
	    y = (coor[1]-y11)*deltaZ/oldDeltaZ + y11;
	    if(TMath::Abs(y-coor[1])>0.1) AliDebug(3, Form("oldDeltaZ = %7.2f   newDeltaZ = %7.2f\toldY = %7.2f   new y = %7.2f\n",oldDeltaZ,deltaZ,coor[1],y));
	}
	Float_t matchDist = PadMatchTrack(xpad, ypad, dpx, dpy, coor[0], y);
	if(matchDist<fkMaxDistance/2.) padsInCheckArea[cathode]++;
	if(matchDist>minMatchDist[cathode])continue;
	isMatch[cathode] = kTRUE;
	minMatchDist[cathode] = matchDist;
	foundDetElemId = currDetElemId;
	foundZmatch=zpad;
	yCoorAtPadZ=y;
	if(cathode==kBending) trigDigitBendPlane = idigit;
	for (Int_t loc=0; loc<pad.GetNofLocations(); loc++){
	    nboard[cathode][loc] = pad.GetLocalBoardId(loc);
	}
	for(Int_t loc=pad.GetNofLocations(); loc<fgkNlocations; loc++){
	    nboard[cathode][loc]=-1;
	}
    }

    if(isMatch[kBending] || isMatch[kNonBending]){
	detElemId = foundDetElemId;
	zRealMatch[ch] = foundZmatch;
	coor[1] = yCoorAtPadZ;
    }

    // If track matches many pads, it is not good for effciency determination.
    // However we still want to calculate the hit pattern.
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
      if(padsInCheckArea[cath]>2) {
	AliDebug(1, Form("padsInCheckArea[%i] = %i\n",cath,padsInCheckArea[cath]));
	return -500;
      }
    }

    return trigDigitBendPlane;
}

//_____________________________________________________________________________
Float_t AliMUONTrackHitPattern::PadMatchTrack(Float_t xPad, Float_t yPad,
						Float_t dpx, Float_t dpy, 
						Float_t xTrackAtPad, Float_t yTrackAtPad) const
{
    //
    /// Decides if the digit belongs to the trigger track.
    //

  Float_t maxDist = GetRecoParam()->GetStripCutForTrigger() * 2. * TMath::Min(dpx,dpy); // cm
  if(maxDist<2.) maxDist = 2.;
  Float_t maxDistCheckArea = GetRecoParam()->GetMaxStripAreaForTrigger() * 2. *  TMath::Min(dpx,dpy); // cm

    Float_t matchDist = fkMaxDistance;

    Float_t deltaX = TMath::Abs(xPad-xTrackAtPad)-dpx;
    Float_t deltaY = TMath::Abs(yPad-yTrackAtPad)-dpy;
    Float_t maxDistX = maxDist;
    Float_t maxDistY = maxDist;

    if(deltaX<=maxDistX && deltaY<=maxDistY) matchDist = TMath::Max(deltaX, deltaY);
    else if(deltaX<=maxDistCheckArea && deltaY<=maxDistCheckArea) matchDist = fkMaxDistance/5.;
    return matchDist;
}


//_____________________________________________________________________________
Int_t AliMUONTrackHitPattern::DetElemIdFromPos(Float_t x, Float_t y, 
						 Int_t chamber, Int_t cathode) const
{
    //
    /// Given the (x,y) position in the chamber,
    /// it returns the corresponding slat
    //

    Int_t resultingDetElemId = -1;
    AliMpDEIterator it;
    Float_t minDist = 999.;
    for ( it.First(chamber-1); ! it.IsDone(); it.Next() ){
	Int_t detElemId = it.CurrentDEId();
	Int_t ich = detElemId/100-10;
	Float_t tolerance=0.2*((Float_t)ich);
	Float_t currDist=9999.;

	const AliMpVSegmentation* seg = 
	    AliMpSegmentation::Instance()
	    ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));
	if (!seg) continue;

	Float_t deltax = seg->GetDimensionX();
	Float_t deltay = seg->GetDimensionY();
	Float_t xlocal1 =  -deltax;
	Float_t ylocal1 =  -deltay;
	Float_t xlocal2 =  +deltax;
	Float_t ylocal2 =  +deltay;
	Float_t xg01, yg01, zg1, xg02, yg02, zg2;
	fkTransformer.Local2Global(detElemId, xlocal1, ylocal1, 0, xg01, yg01, zg1);
	fkTransformer.Local2Global(detElemId, xlocal2, ylocal2, 0, xg02, yg02, zg2);

	Float_t xg1 = xg01, xg2 = xg02, yg1 = yg01, yg2 = yg02;

	if(xg01>xg02){
	    xg1 = xg02;
	    xg2 = xg01;
	}
	if(yg01>yg02){
	    yg1 = yg02;
	    yg2 = yg01;
	}

	if(x>=xg1-tolerance && x<=xg2+tolerance && y>=yg1-tolerance && y<=yg2+tolerance){ // takes into account errors in extrapolation
	    if(y<yg1) currDist = yg1-y;
	    else if(y>yg2) currDist = y-yg2;
	    if(currDist<minDist) {
		resultingDetElemId = detElemId;
		minDist=currDist;
		continue;
	    }
	    resultingDetElemId = detElemId;
	    break;
	}
    } // loop on detElemId
    return resultingDetElemId;
}


//_____________________________________________________________________________
void AliMUONTrackHitPattern::LocalBoardFromPos(Float_t x, Float_t y,
						 Int_t detElemId, Int_t cathode,
						 Int_t localBoard[4]) const
{
    //
    /// Given the (x,y) position in the chamber,
    /// it returns the corresponding local board
    //

    for(Int_t loc=0; loc<fgkNlocations; loc++){
	localBoard[loc]=-1;
    }
    Float_t xg, yg, zg;
    fkTransformer.Local2Global(detElemId, 0., 0., 0., xg, yg, zg);
    Float_t z = zg - (y-yg)*TMath::Tan(TMath::DegToRad()*AliMUONConstants::St345Inclination());
    Float_t xl, yl, zl;
    fkTransformer.Global2Local(detElemId, x, y, z, xl, yl, zl);
    const AliMpVSegmentation* seg = 
	AliMpSegmentation::Instance()
          ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));
    if ( ! seg ) return;
    AliMpPad pad = seg->PadByPosition(xl,yl,kFALSE);
    if ( ! pad.IsValid() ){
      AliWarning(Form("Pad not found! DetElemId %i  position global (%f, %f, %f)  local (%f, %f, %f)\n", detElemId, x, y, z, xl, yl, zl));
      return;
    }

    AliDebug(2, Form("\nWould match: DetElemId = %i  Pad = (%i,%i) = (%.2f,%.2f,%.2f)   Board %3i\n",detElemId,pad.GetIx(),pad.GetIy(),xl,yl,zl,pad.GetLocalBoardId(0)));

    for (Int_t loc=0; loc<pad.GetNofLocations(); loc++){
      localBoard[loc] = pad.GetLocalBoardId(loc);
    }
}


//_____________________________________________________________________________
Bool_t AliMUONTrackHitPattern::PerformTrigTrackMatch(UShort_t &pattern,
						     const AliMUONTriggerTrack* matchedTrigTrack,
						     AliMUONVDigitStore& digitStore) const
{
  //
  /// It searches for matching digits around the trigger track.
  //

  enum {kBending, kNonBending};

  Int_t chOrder[fgkNchambers] = {0,2,1,3};

  TArrayF zRealMatch(fgkNchambers);

  Bool_t isMatch[fgkNcathodes];
  for(Int_t cath=0; cath<fgkNcathodes; cath++){
    isMatch[cath] = kFALSE;
  }

  TArrayF zMeanChamber(fgkNchambers);
  zMeanChamber[0] = matchedTrigTrack->GetZ11();
  zMeanChamber[1] = matchedTrigTrack->GetZ11() + AliMUONConstants::DefaultChamberZ(11) - AliMUONConstants::DefaultChamberZ(10);
  zMeanChamber[2] = matchedTrigTrack->GetZ21();
  zMeanChamber[3] = matchedTrigTrack->GetZ21() + AliMUONConstants::DefaultChamberZ(13) - AliMUONConstants::DefaultChamberZ(12);

  TArrayI digitPerTrack(fgkNcathodes);

  Float_t trackIntersectCh[fgkNchambers][fgkNcathodes];

  TArrayI triggeredDigits;
  triggeredDigits.Set(fgkNchambers);
  triggeredDigits.Reset(-1);

  TArrayI trigScheme[fgkNcathodes];
  TArrayI slatThatTriggered[fgkNcathodes];
  for(Int_t cath=0; cath<fgkNcathodes; cath++){
    trigScheme[cath].Set(fgkNchambers);
    slatThatTriggered[cath].Set(fgkNchambers);
  }

  Int_t boardThatTriggered[fgkNchambers][fgkNcathodes][fgkNlocations];
  TArrayI nboard[fgkNcathodes];
  for(Int_t cath=0; cath<fgkNcathodes; cath++){
    nboard[cath].Set(fgkNlocations);
  }
  Int_t ineffBoard[fgkNlocations];
  for(Int_t loc=0; loc<fgkNlocations; loc++){
    ineffBoard[loc] = -1;
  }

  digitPerTrack.Reset();
  for(Int_t ch=0; ch<fgkNchambers; ch++){
    zRealMatch[ch] = zMeanChamber[ch];
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
      for(Int_t loc=0; loc<fgkNlocations; loc++){
	boardThatTriggered[ch][cath][loc]=-1;
      }
    }
  }

  for(Int_t cath=0; cath<fgkNcathodes; cath++){
    slatThatTriggered[cath].Reset(-1);
    trigScheme[cath].Reset();
  }

  Bool_t isClearEvent = kTRUE;

  Float_t y11 = matchedTrigTrack->GetY11();// y position (info from bending plane)
  Float_t slopeX = matchedTrigTrack->GetSlopeX();
  Float_t slopeY = matchedTrigTrack->GetSlopeY();

  for(Int_t ch=0; ch<fgkNchambers; ch++) { // chamber loop
    Int_t currCh = chOrder[ch];
    AliDebug(3, Form("zMeanChamber[%i] = %.2f\tzRealMatch[0] = %.2f\n",currCh,zMeanChamber[currCh],zRealMatch[0]));

    // searching track intersection with chambers (first approximation)
    Float_t deltaZ = zMeanChamber[currCh] - zMeanChamber[0];
    trackIntersectCh[currCh][0] = zMeanChamber[currCh] * slopeX;
    trackIntersectCh[currCh][1] = y11 + deltaZ * slopeY;
    Int_t detElemIdFromTrack = DetElemIdFromPos(trackIntersectCh[currCh][0], trackIntersectCh[currCh][1], 11+currCh, 0);
    if(detElemIdFromTrack<0) {
      AliDebug(1, "Warning: trigger track outside trigger chamber\n");
      isClearEvent = kFALSE;

      // track is rejected since the extrapolated track
      // does not match a slat (border effects)
      AliESDMuonTrack::AddEffInfo(pattern, AliESDMuonTrack::kTrackOutsideGeometry);
      continue;
    }
		
    triggeredDigits[currCh] = FindPadMatchingTrig(digitStore, detElemIdFromTrack, trackIntersectCh[currCh], isMatch, nboard, zRealMatch, y11);

    // if FindPadMatchingTrig = -500 => too many digits matching pad =>
    //                               => Event not clear => Do not use for efficiency calculation
    if(triggeredDigits[currCh]<-100){
      isClearEvent = kFALSE;
      // track is rejected since it matches many pads
      AliESDMuonTrack::AddEffInfo(pattern, AliESDMuonTrack::kTrackMatchesManyPads);
      AliDebug(1, Form("Warning: track = %p (%i) matches many pads. Rejected!\n",(void *)matchedTrigTrack, detElemIdFromTrack));
    }

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
      if(!isMatch[cath]) continue;
      AliESDMuonTrack::SetFiredChamber(pattern, cath, currCh);
      digitPerTrack[cath]++;
      trigScheme[cath][currCh]++;
      slatThatTriggered[cath][currCh] = detElemIdFromTrack;
      for(Int_t loc=0; loc<fgkNlocations; loc++){
	boardThatTriggered[currCh][cath][loc] = nboard[cath][loc];
      }
    }
  } // end chamber loop

  for(Int_t cath=0; cath<fgkNcathodes; cath++){
    if(digitPerTrack[cath]<3) {
      isClearEvent = kFALSE;
      // track is rejected since the number of associated
      // digits found is less than 3.
      AliESDMuonTrack::AddEffInfo(pattern, AliESDMuonTrack::kTrackMatchesFewPads);
      AliDebug(1, Form("Warning: found %i digits for trigger track cathode %i.\nRejecting event\n", digitPerTrack[cath],cath));
    }
  }

  if(!isClearEvent) return kFALSE;

  AliESDMuonTrack::EAliTriggerChPatternFlag goodForEff = AliESDMuonTrack::kBoardEff;

  Int_t ineffSlat = -1;
  Int_t ineffDetElId = -1;
  Int_t firstSlat = slatThatTriggered[kBending][0]%100;
  if(firstSlat<0) firstSlat = slatThatTriggered[kBending][1]%100;
  Int_t firstBoard = boardThatTriggered[0][kBending][0];
  if(firstBoard<0) firstBoard = boardThatTriggered[1][kBending][0];
  for(Int_t ch=0; ch<fgkNchambers; ch++){
    Bool_t isCurrChIneff = kFALSE;
    Int_t currSlat = slatThatTriggered[kBending][ch]%100;
    if(currSlat<0){
      ineffDetElId = DetElemIdFromPos(trackIntersectCh[ch][0], trackIntersectCh[ch][1], 11+ch, kBending);
      currSlat = ineffDetElId%100;
      ineffSlat = currSlat;
      isCurrChIneff = kTRUE;
    }
    if(currSlat!=firstSlat) {
      AliESDMuonTrack::AddEffInfo(pattern,
				  AliESDMuonTrack::kCrossDifferentSlats,
				  AliESDMuonTrack::kChEff); // track crosses different slats
      return kTRUE;
    }
    Bool_t atLeastOneLoc=kFALSE;
    if(isCurrChIneff) LocalBoardFromPos(trackIntersectCh[ch][0], trackIntersectCh[ch][1], ineffDetElId, kBending, ineffBoard);
    for(Int_t loc=0; loc<fgkNlocations; loc++){
      Int_t currBoard = boardThatTriggered[ch][kBending][loc];
      if(isCurrChIneff) currBoard = ineffBoard[loc];
      if(currBoard==firstBoard){
	atLeastOneLoc=kTRUE;
	break;
      }
    }
    if(!atLeastOneLoc) goodForEff = AliESDMuonTrack::kSlatEff;
  } // end chamber loop
  
  AliESDMuonTrack::AddEffInfo(pattern, firstSlat, goodForEff);
  return kTRUE;
}
