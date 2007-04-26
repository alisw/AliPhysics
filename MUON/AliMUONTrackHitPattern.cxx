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

////////////////////////////////////
///
/// \class AliMUONTrackHitPattern
///
/// This class propagates tracks to trigger chambers 
/// searching for fired strips.
///
/// To each track, a hit pattern for trigger chambers is set.
/// The hit pattern is a UShort_t with 8 bits used:
///
///            1  1  0  1    1  1  0  1
///           |           |            |
///            ----------- ------------
/// chamber:  11 12 13 14 | 11 12 13 14
/// cathode:    bending   | non-bending
///
/// The main method is:
/// * GetHitPattern
///
///  \author Diego Stocco
///
////////////////////////////////////


#include "AliMUONTrackHitPattern.h"
#include "AliMUONData.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONDigit.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONLocalTriggerBoard.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONDigitMaker.h"

#include "AliMpPad.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"

#include "AliLog.h"
#include "AliTracker.h"
#include "AliMagF.h"

#include <Riostream.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TArrayS.h>

/// \cond CLASSIMP
ClassImp(AliMUONTrackHitPattern) // Class implementation in ROOT context
/// \endcond


//______________________________________________________________________________
AliMUONTrackHitPattern::AliMUONTrackHitPattern(AliMUONData *data)
    : TObject(),
      fMUONData(data),
      fTransformer(new AliMUONGeometryTransformer(kTRUE)),
      fCrateManager(new AliMUONTriggerCrateStore()),
      fDigitMaker(new AliMUONDigitMaker())
{
    /// Default constructor

    // Set magnetic field
    const AliMagF* kField = AliTracker::GetFieldMap();
    if (!kField) AliFatal("No field available");
    AliMUONTrackExtrap::SetField(kField);

    // Geometry transformer
    fTransformer->ReadGeometryData("volpath.dat", "geometry.root");

    // Crate manager to retrieve local boards
    fCrateManager->ReadFromFile();

    // set to digit maker
    fDigitMaker->SetCrateManager(fCrateManager);

    for(Int_t ch=0; ch<4; ch++){
	fTriggerDigitsList[ch].Clear();
    }
}


//______________________________________________________________________________
AliMUONTrackHitPattern::~AliMUONTrackHitPattern(void)
{
/// Destructor
    for(Int_t ch=0; ch<4; ch++){
	fTriggerDigitsList[ch].Delete();
    }
    delete fCrateManager;
}


//______________________________________________________________________________
void AliMUONTrackHitPattern::GetHitPattern(TClonesArray *recTracksPtr)
{
    //
    /// Main method:
    /// Loops on reco tracks, extrapolates them to trigger chambers
    /// and searches for matching digits
    //
    
    const Int_t mask[2][4]={{0x80, 0x40, 0x20, 0x10},
			    {0x08, 0x04, 0x02, 0x01}};
    Bool_t isMatch[2];
    UShort_t pattern=0;
    TriggerDigits();
    Int_t nRecTracks = (Int_t)recTracksPtr->GetEntriesFast();
    for(Int_t iTrack=0; iTrack<nRecTracks; iTrack++){
	pattern = 0;
	AliMUONTrack *muonTrack = (AliMUONTrack*) recTracksPtr->At(iTrack);
	AliMUONTrackParam *trackParam = (AliMUONTrackParam*) ((muonTrack->GetTrackParamAtHit())->Last());
	for(Int_t ch=0; ch<4; ch++){
	    AliMUONTrackExtrap::ExtrapToZCov(trackParam, AliMUONConstants::DefaultChamberZ(10+ch));
	    FindPadMatchingTrack(trackParam, isMatch, ch);
	    for(Int_t cath=0; cath<2; cath++){
		if(isMatch[cath]) pattern |= mask[cath][ch];
	    }
	}
	muonTrack->SetHitsPatternInTrigCh(pattern);
    }
    return;
}


//______________________________________________________________________________
void AliMUONTrackHitPattern::FindPadMatchingTrack(AliMUONTrackParam *trackParam,
						  Bool_t isMatch[2], Int_t iChamber)
{
    //
    /// Given track position, searches for matching digits.
    //

    Float_t minMatchDist[2];

    for(Int_t cath=0; cath<2; cath++){
	isMatch[cath]=kFALSE;
	minMatchDist[cath]=9999.;
    }

    Int_t ndigits = (Int_t)fTriggerDigitsList[iChamber].GetEntries();
    AliMUONDigit * mDigit = 0x0;
    for(Int_t idigit=0; idigit<ndigits; idigit++) { // digit loop
	mDigit = (AliMUONDigit*)fTriggerDigitsList[iChamber].At(idigit);
	Int_t currDetElemId = mDigit->DetElemId();

	Int_t cathode = mDigit->Cathode();
	Int_t ix = mDigit->PadX();
	Int_t iy = mDigit->PadY();
	Float_t xpad, ypad, zpad;
	const AliMpVSegmentation* seg = AliMpSegmentation::Instance()
	    ->GetMpSegmentation(currDetElemId,AliMp::GetCathodType(cathode));

	AliMpPad pad = seg->PadByIndices(AliMpIntPair(ix,iy),kTRUE);
	Float_t xlocal1 = pad.Position().X();
	Float_t ylocal1 = pad.Position().Y();
	Float_t dpx = pad.Dimensions().X();
	Float_t dpy = pad.Dimensions().Y();
	fTransformer->Local2Global(currDetElemId, xlocal1, ylocal1, 0, xpad, ypad, zpad);
	Float_t matchDist = MinDistanceFromPad(xpad, ypad, zpad, dpx, dpy, trackParam);
	if(matchDist>minMatchDist[cathode])continue;
	isMatch[cathode] = kTRUE;
	minMatchDist[cathode] = matchDist;
    }

    return;
}


//______________________________________________________________________________
Float_t AliMUONTrackHitPattern::MinDistanceFromPad(Float_t xPad, Float_t yPad, Float_t zPad,
						   Float_t dpx, Float_t dpy, AliMUONTrackParam *trackParam)
{
    //
    /// Decides if the digit belongs to the track.
    //
    Float_t xTrackAtPad = trackParam->GetNonBendingCoor();
    Float_t yTrackAtPad = trackParam->GetBendingCoor();

    Float_t sigmaX, sigmaY, sigmaMS;
    GetPosUncertainty(trackParam, zPad, sigmaX, sigmaY, sigmaMS);

    Float_t maxDistX = 3.*(sigmaX + sigmaMS); // in cm
    Float_t maxDistY = 3.*(sigmaY + sigmaMS); // in cm

    Float_t deltaX = TMath::Abs(xPad-xTrackAtPad)-dpx;
    Float_t deltaY = TMath::Abs(yPad-yTrackAtPad)-dpy;

    Float_t matchDist = 99999.;
    if(deltaX<=maxDistX && deltaY<=maxDistY) matchDist = TMath::Max(deltaX, deltaY);
    return matchDist;
}


//______________________________________________________________________________
void AliMUONTrackHitPattern::GetPosUncertainty(AliMUONTrackParam *trackParam, Float_t zChamber, 
					       Float_t &sigmaX, Float_t &sigmaY, Float_t &sigmaMS)
{
    //
    /// Returns uncertainties on extrapolated position.
    /// Takes into account Branson plane corrections in the iron wall.
    //

    const Float_t alpha = 0.1123; // GeV/c
    
    // Find a better way to get such parameters ???
    const Float_t kZFilterIn = 1471.; // From STRUCT/SHILConst2.h
    const Float_t kZFilterOut = kZFilterIn + 120.; // From STRUCT/SHILConst2.h
    
    const Float_t zBranson = - (kZFilterIn + (kZFilterOut - kZFilterIn)*2./3. ); // - sign because distance are positive
    Float_t zDistFromWall = TMath::Abs(zChamber - zBranson);
    Float_t zDistFromLastTrackCh = TMath::Abs(zChamber - AliMUONConstants::DefaultChamberZ(9));

    TMatrixD *covParam = trackParam->GetCovariances();
    
    sigmaX = (*covParam)(0,0);
    sigmaY = (*covParam)(2,2);

    // If covariance matrix is not extrapolated, use "reasonable" errors
    // (To be removed as soon as covariance matrix is correctly propagated).
    if (sigmaX==0.)sigmaX = 0.003 * zDistFromLastTrackCh;
    if (sigmaY==0.)sigmaY = 0.004 * zDistFromLastTrackCh;

    Float_t p = trackParam->P();
    Float_t thetaMS = alpha/p;
    sigmaMS = zDistFromWall * TMath::Tan(thetaMS);

    return;
}


//____________________________________________________________________
Bool_t AliMUONTrackHitPattern::TriggerDigits()
{
    //
    /// make (S)Digit for trigger
    //


    Int_t nBoard;

    TList digitList;

    digitList.Clear();

    AliMUONLocalTrigger *locTrg = 0x0;

    fMUONData->SetTreeAddress("RC,TC");
    fMUONData->GetTrigger();

    TClonesArray *localTrigger = fMUONData->LocalTrigger();
    Int_t nLocTrig = (Int_t) localTrigger->GetEntriesFast();

    for(Int_t iLoc=0; iLoc<nLocTrig; iLoc++) {
      locTrg = (AliMUONLocalTrigger*)localTrigger->UncheckedAt(iLoc);

      TArrayS xyPattern[2];
      xyPattern[0].Set(4);
      xyPattern[1].Set(4);

      xyPattern[0].AddAt(locTrg->GetX1Pattern(),0);
      xyPattern[0].AddAt(locTrg->GetX2Pattern(),1);
      xyPattern[0].AddAt(locTrg->GetX3Pattern(),2);
      xyPattern[0].AddAt(locTrg->GetX4Pattern(),3);

      xyPattern[1].AddAt(locTrg->GetY1Pattern(),0);
      xyPattern[1].AddAt(locTrg->GetY2Pattern(),1);
      xyPattern[1].AddAt(locTrg->GetY3Pattern(),2);
      xyPattern[1].AddAt(locTrg->GetY4Pattern(),3);

      for(Int_t cath=0; cath<2; cath++){
	  for(Int_t ch=0; ch<4; ch++){
	      if(xyPattern[cath][ch]==0) continue;
	  }
      }

      nBoard    = locTrg->LoCircuit();
      fDigitMaker->TriggerDigits(nBoard, xyPattern, digitList);


    } // loop on localTriggers

    for (Int_t iEntry = 0; iEntry < digitList.GetEntries(); ++iEntry) {
      AliMUONDigit* digit = (AliMUONDigit*)digitList.At(iEntry);
      Int_t detElemId = digit->DetElemId();
      Int_t iChamber  = detElemId/100 - 11; //FIXEME should be given by mapping
      fTriggerDigitsList[iChamber].Add(digit);

    }

    return kTRUE;
}
