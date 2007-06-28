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

/// \class AliMUONTriggerChamberEff
/// implementation of the trigger chamber efficiency determination from
/// data, and returns the
/// efficiencyCells.dat with the calculated efficiencies
///
/// \author Diego Stocco (Torino)

#include "AliMUONTriggerChamberEff.h"
#include "AliMUONVDigit.h"
#include "AliMUONConstants.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONDigitMaker.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGeometryTransformer.h"

#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"

#include "AliMUONDigitStoreV1.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVTriggerStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONVTrackStore.h"

#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpPad.h"
#include "AliMpDEIterator.h"
#include "AliMpPlaneType.h"
#include "AliMpDEManager.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>

#include <TSeqCollection.h>
#include <TTree.h>
#include <TROOT.h>


/// \cond CLASSIMP
ClassImp(AliMUONTriggerChamberEff)
/// \endcond


//_____________________________________________________________________________
AliMUONTriggerChamberEff::AliMUONTriggerChamberEff()
: TObject(),
  fTransformer(0x0),
  fDigitMaker(0x0),
  fReproduceTrigResponse(kFALSE),
  fPrintInfo(kFALSE),
  fWriteOnESD(kFALSE),
  fDebugLevel(0),
  fkMaxDistance(99999.)
{
/// Standard constructor
    ResetArrays();
}


//_____________________________________________________________________________
AliMUONTriggerChamberEff::AliMUONTriggerChamberEff(const AliMUONGeometryTransformer* transformer,
						   const AliMUONDigitMaker* digitMaker,
						   Bool_t writeOnESD)
: TObject(),
  fTransformer(transformer),
  fDigitMaker(digitMaker),
  fReproduceTrigResponse(kFALSE),
  fPrintInfo(kFALSE),
  fWriteOnESD(writeOnESD),
  fDebugLevel(0),
  fkMaxDistance(99999.)
{
/// Standard constructor
    ResetArrays();
}


//_____________________________________________________________________________
AliMUONTriggerChamberEff::~AliMUONTriggerChamberEff()
{
/// Destructor
    Bool_t writeOnESD=fWriteOnESD;
    fWriteOnESD=kFALSE;
    if(writeOnESD) SaveInESDFile();
}


//_____________________________________________________________________________
AliMUONTriggerChamberEff::AliMUONTriggerChamberEff(const AliMUONTriggerChamberEff& other)
    :TObject(other),
     fTransformer(0x0),
     fDigitMaker(0x0),
     fReproduceTrigResponse(other.fReproduceTrigResponse),
     fPrintInfo(other.fPrintInfo),
     fWriteOnESD(other.fWriteOnESD),
     fDebugLevel(other.fDebugLevel),
     fkMaxDistance(other.fkMaxDistance)
{
/// copy constructor

    for(Int_t ch=0; ch<fgkNchambers; ch++){
	for(Int_t cath=0; cath<fgkNcathodes; cath++){
	    fTrigger34[ch][cath] = other.fTrigger34[ch][cath];
	    fTrigger44[cath] = other.fTrigger44[cath];
	    for(Int_t slat=0; slat<fgkNslats; slat++){
		fInefficientSlat[ch][cath][slat] = other.fInefficientSlat[ch][cath][slat];
		fHitPerSlat[ch][cath][slat] = other.fHitPerSlat[ch][cath][slat];
	    }
	    for(Int_t board=0; board<fgkNboards; board++){
		fInefficientBoard[ch][cath][board] = other.fInefficientBoard[ch][cath][board];
		fHitPerBoard[ch][cath][board] = other.fHitPerBoard[ch][cath][board];
	    }
	}
    }
}


//_____________________________________________________________________________
AliMUONTriggerChamberEff& AliMUONTriggerChamberEff::operator=(const AliMUONTriggerChamberEff& other)
{
    /// Asignment operator
    // check assignement to self
    if (this == &other)
	return *this;

    // base class assignement
    TObject::operator=(other);

    fTransformer = 0x0;
    fDigitMaker = 0x0;
    fReproduceTrigResponse = other.fReproduceTrigResponse;
    fPrintInfo = other.fPrintInfo;
    fWriteOnESD = other.fWriteOnESD;
    fDebugLevel = other.fDebugLevel;
    //fkMaxDistance = other.fkMaxDistance;
    return *this;
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::ResetArrays()
{
    //
    /// Sets the data member counters to 0.
    //

    for(Int_t ch=0; ch<fgkNchambers; ch++){
	for(Int_t cath=0; cath<fgkNcathodes; cath++){
	    fTrigger34[ch][cath] = 0;
	    fTrigger44[cath] = 0;
	    for(Int_t slat=0; slat<fgkNslats; slat++){
		fInefficientSlat[ch][cath][slat] = 0;
		fHitPerSlat[ch][cath][slat] = 0;
	    }
	    for(Int_t board=0; board<fgkNboards; board++){
		fInefficientBoard[ch][cath][board] = 0;
		fHitPerBoard[ch][cath][board] = 0;
	    }
	}
    }
}


//______________________________________________________________________________
Bool_t 
AliMUONTriggerChamberEff::TriggerDigits(const AliMUONVTriggerStore& triggerStore,
					AliMUONVDigitStore& digitStore) const
{
    //
    /// make (S)Digit for trigger
    //
    digitStore.Clear();

    AliMUONLocalTrigger* locTrg;
    TIter next(triggerStore.CreateLocalIterator());
  
    while ( ( locTrg = static_cast<AliMUONLocalTrigger*>(next()) ) ) 
    {
	if (locTrg->IsNull()) continue;
   
	TArrayS xyPattern[2];
	locTrg->GetXPattern(xyPattern[0]);
	locTrg->GetYPattern(xyPattern[1]);
    
	Int_t nBoard = locTrg->LoCircuit();
	fDigitMaker->TriggerDigits(nBoard, xyPattern, digitStore);
    }
    return kTRUE;
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::InfoDigit(AliMUONVDigitStore& digitStore)
{
    //
    /// Prints information on digits (for debugging)
    //
    TIter next(digitStore.CreateIterator());
    AliMUONVDigit* mDigit=0x0;
    
    while ( ( mDigit = static_cast<AliMUONVDigit*>(next()) ) )
    {
	mDigit->Print();
    } // end digit loop
    printf("\n");
}


//_____________________________________________________________________________
Int_t AliMUONTriggerChamberEff::MatchingPad(AliMUONVDigitStore& digitStore, Int_t &detElemId,
					    Float_t coor[2], Bool_t isMatch[fgkNcathodes],
					    Int_t nboard[fgkNcathodes][4],
					    Float_t zRealMatch[fgkNchambers], Float_t y11)
{
    //
    /// Check slat and board number of digit matching track
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

	AliMpPad pad = seg->PadByIndices(AliMpIntPair(ix,iy),kTRUE);
	Float_t xlocal1 = pad.Position().X();
	Float_t ylocal1 = pad.Position().Y();
	Float_t dpx = pad.Dimensions().X();
	Float_t dpy = pad.Dimensions().Y();
	fTransformer->Local2Global(currDetElemId, xlocal1, ylocal1, 0, xpad, ypad, zpad);
	if(fDebugLevel>2)printf("DetElemId = %i\tCathode = %i\t(x,y) Pad = (%i,%i) = (%.2f,%.2f)\tDim = (%.2f,%.2f)\tTrack = (%.2f,%.2f)\n",currDetElemId,cathode,ix,iy,xpad,ypad,dpx,dpy,coor[0],coor[1]);
	// searching track intersection with chambers (second approximation)
	if(ch%2==1){
	    //if(iChamber%2==1){
	    Float_t deltaZ = zpad - zRealMatch[0];
	    y = (coor[1]-y11)*deltaZ/oldDeltaZ + y11;
	    if(fDebugLevel>=3 && TMath::Abs(y-coor[1])>0.1)printf("oldDeltaZ = %7.2f   newDeltaZ = %7.2f\toldY = %7.2f   new y = %7.2f\n",oldDeltaZ,deltaZ,coor[1],y);
	}
	Float_t matchDist = PadMatchTrack(xpad, ypad, dpx, dpy, coor[0], y, ch);
	if(matchDist<fkMaxDistance/2.) padsInCheckArea[cathode]++;
	if(matchDist>minMatchDist[cathode])continue;
	isMatch[cathode] = kTRUE;
	minMatchDist[cathode] = matchDist;
	foundDetElemId = currDetElemId;
	foundZmatch=zpad;
	yCoorAtPadZ=y;
	if(cathode==kBending) trigDigitBendPlane = idigit;
	for (Int_t loc=0; loc<pad.GetNofLocations(); loc++){
	    AliMpIntPair location = pad.GetLocation(loc);
	    nboard[cathode][loc] = location.GetFirst();
	}
	for(Int_t loc=pad.GetNofLocations(); loc<4; loc++){
	    nboard[cathode][loc]=-1;
	}
    }

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	if(padsInCheckArea[cath]>2) {
	    if(fDebugLevel>=1) printf("padsInCheckArea[%i] = %i\n",cath,padsInCheckArea[cath]);
	    return -500;
	}
    }

    if(isMatch[kBending] || isMatch[kNonBending]){
	detElemId = foundDetElemId;
	zRealMatch[ch] = foundZmatch;
	coor[1] = yCoorAtPadZ;
	if(fDebugLevel>2){
	    Int_t whichCathode=kBending;
	    if(!isMatch[kBending])whichCathode=kNonBending;
	}
    }
    return trigDigitBendPlane;
}

//_____________________________________________________________________________
Float_t AliMUONTriggerChamberEff::PadMatchTrack(Float_t xPad, Float_t yPad,
						Float_t dpx, Float_t dpy, 
						Float_t xTrackAtPad, Float_t yTrackAtPad,
						Int_t chamber)
{
    //
    /// Decides if the digit belongs to the trigger track.
    //

    Float_t maxDist = 2.;//3. // cm
    Float_t maxDistCheckArea = 6.; // cm

    Float_t matchDist = fkMaxDistance;

    Float_t deltaX = TMath::Abs(xPad-xTrackAtPad)-dpx;
    Float_t deltaY = TMath::Abs(yPad-yTrackAtPad)-dpy;
    Float_t maxDistX = maxDist;
    Float_t maxDistY = maxDist;
    
    if(fReproduceTrigResponse){
	maxDistX = dpx;
	maxDistY = dpy;
	deltaX = TMath::Abs(xPad-xTrackAtPad);
	deltaY = TMath::Abs(yPad-yTrackAtPad);
	if(dpx<dpy && chamber>=2) maxDistX = 3.*dpx;// Non-bending plane: check the +- 1 strip between stations
	if(dpy<dpx && chamber%2) maxDistY = 3.*dpy;// bending plane: check the +- 1 strip between planes in the same station
    }

    if(deltaX<=maxDistX && deltaY<=maxDistY) matchDist = TMath::Max(deltaX, deltaY);
    else if(deltaX<=maxDistCheckArea && deltaY<=maxDistCheckArea) matchDist = fkMaxDistance/5.;
    return matchDist;
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::CalculateEfficiency(Int_t trigger44, Int_t trigger34,
						   Float_t &efficiency, Float_t &error,
						   Bool_t failuresAsInput)
{
    //
    /// Returns the efficiency.
    //

    efficiency=-9.;
    error=0.;
    if(trigger34>0){
	efficiency=(Double_t)trigger44/((Double_t)trigger34);
	if(failuresAsInput)efficiency=1.-(Double_t)trigger44/((Double_t)trigger34);
    }
    Double_t q = TMath::Abs(1-efficiency);
    if(efficiency<0)error=0.0;
    else error = TMath::Sqrt(efficiency*q/((Double_t)trigger34));
}


//_____________________________________________________________________________
Int_t AliMUONTriggerChamberEff::DetElemIdFromPos(Float_t x, Float_t y, 
						 Int_t chamber, Int_t cathode)
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

	Float_t deltax = seg->Dimensions().X();
	Float_t deltay = seg->Dimensions().Y();
	Float_t xlocal1 =  -deltax;
	Float_t ylocal1 =  -deltay;
	Float_t xlocal2 =  +deltax;
	Float_t ylocal2 =  +deltay;
	Float_t xg01, yg01, zg1, xg02, yg02, zg2;
	fTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg01, yg01, zg1);
	fTransformer->Local2Global(detElemId, xlocal2, ylocal2, 0, xg02, yg02, zg2);

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
void AliMUONTriggerChamberEff::LocalBoardFromPos(Float_t x, Float_t y,
						 Int_t detElemId, Int_t cathode,
						 Int_t localBoard[4])
{
    //
    /// Given the (x,y) position in the chamber,
    /// it returns the corresponding local board
    //

    for(Int_t loc=0; loc<4; loc++){
	localBoard[loc]=-1;
    }
    Float_t xl, yl, zl;
    fTransformer->Global2Local(detElemId, x, y, 0, xl, yl, zl);
    TVector2 pos(xl,yl);
    const AliMpVSegmentation* seg = 
	AliMpSegmentation::Instance()
          ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));
    if (seg){
	AliMpPad pad = seg->PadByPosition(pos,kFALSE);
	for (Int_t loc=0; loc<pad.GetNofLocations(); loc++){
	    AliMpIntPair location = pad.GetLocation(loc);
	    localBoard[loc] = location.GetFirst();
	}
    }
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::EventChamberEff(const AliMUONVTriggerStore& triggerStore,
					       const AliMUONVTriggerTrackStore& trigTrackStore,
					       const AliMUONVTrackStore& trackStore)
{
    //
    /// Main method.
    /// It loops over the the trigger rec. tracks in the event.
    /// Then it search for matching digits around the track.
    /// Finally it calculates the efficiency for each trigger board.
    /// Files with calculated efficiency are placed in the user defined outputDir.
    //

    if(!fTransformer || ! fDigitMaker) {
	AliError(Form("AliMUONGeometryTransformer or AliMUONDigitMaker not properly initialized!!"));
	return;
    }

    enum {kBending, kNonBending};
    Float_t rad2deg = 180./TMath::Pi();

    Int_t chOrder[fgkNchambers] = {0,2,1,3};
    Float_t zRealMatch[fgkNchambers] = {0.0};
    Float_t correctFactor[fgkNcathodes] = {1.};

    Bool_t match[fgkNchambers][fgkNcathodes] = {{kFALSE}};
    Bool_t matchPad[fgkNcathodes]={kFALSE};


    Float_t zMeanChamber[fgkNchambers];
    for(Int_t ch=0; ch<fgkNchambers; ch++){
	zMeanChamber[ch] = AliMUONConstants::DefaultChamberZ(10+ch);
    }

    Int_t digitPerTrack[fgkNcathodes] = {0};

    Float_t trackIntersectCh[fgkNchambers][2]={{0.0}};

    Int_t triggeredDigits[2][fgkNchambers] = {{-1}};

    Int_t trigScheme[fgkNchambers][fgkNcathodes]={{0}};
    Int_t slatThatTriggered[fgkNchambers][fgkNcathodes]={{-1}};
    Int_t boardThatTriggered[fgkNchambers][fgkNcathodes][4]={{{-1}}};
    Int_t nboard[fgkNcathodes][4]={{-1}};
    Int_t ineffBoard[4]={-1};

    AliMUONDigitStoreV1 digitStore;   
    TriggerDigits(triggerStore,digitStore);

    for(Int_t ch=0; ch<fgkNchambers; ch++){
	for(Int_t itrack=0; itrack<2; itrack++){
	    triggeredDigits[itrack][ch]=-1;
	}
    }

    AliMUONTriggerTrack *recTrigTrack = 0x0;

    TIter next(trigTrackStore.CreateIterator());
    
    while ( ( recTrigTrack = static_cast<AliMUONTriggerTrack*>(next()) ) )
    {
	for(Int_t cath=0; cath<fgkNcathodes; cath++){
	    digitPerTrack[cath]=0;
	}
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		match[ch][cath]=kFALSE;
		slatThatTriggered[ch][cath]=-1;
		for(Int_t loc=0; loc<4; loc++){
		    boardThatTriggered[ch][cath][loc]=-1;
		}
	    }
	}

	Bool_t isClearEvent = kTRUE;
	Bool_t doubleCountTrack = kFALSE;

	if(!IsCleanTrack(recTrigTrack, trackStore)) {
	    if(fDebugLevel>=1) printf("\tTrack %p (%f, %f) don't match tracker track: rejected!\n",recTrigTrack,recTrigTrack->GetX11(),recTrigTrack->GetY11());
	    continue;
	}

	Float_t x11 = recTrigTrack->GetX11();// x position (info from non-bending plane)
	Float_t y11 = recTrigTrack->GetY11();// y position (info from bending plane)
	Float_t thetaX = recTrigTrack->GetThetax();
	Float_t thetaY = recTrigTrack->GetThetay();

	if(fDebugLevel>=3) printf("\tTrack = %p\npos from track: (x,y) = (%f, %f), (thetaX, thetaY) = (%f, %f)\n",recTrigTrack,x11,y11,thetaX*rad2deg,thetaY*rad2deg);

	for(Int_t ch=0; ch<fgkNchambers; ch++) {
	    zRealMatch[ch] = zMeanChamber[ch];
	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		trigScheme[ch][cath] = 0;
	    }
	}

	for(Int_t ch=0; ch<fgkNchambers; ch++) { // chamber loop
	    Int_t currCh = chOrder[ch];
	    if(fDebugLevel>=2)
		printf("zMeanChamber[%i] = %.2f\tzRealMatch[0] = %.2f\n",currCh,zMeanChamber[currCh],zRealMatch[0]);

	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		correctFactor[cath]=1.;
	    }
	    // calculate corrections to trigger track theta
	    if(ch>=1)correctFactor[kNonBending] = zMeanChamber[0]/zRealMatch[0];// corrects x position
	    if(ch>=2)correctFactor[kBending] = (zMeanChamber[2] - zMeanChamber[0]) / (zRealMatch[2] - zRealMatch[0]);// corrects y position

	    // searching track intersection with chambers (first approximation)
	    Float_t deltaZ = zMeanChamber[currCh] - zMeanChamber[0];
	    trackIntersectCh[currCh][0] = zMeanChamber[currCh] * TMath::Tan(thetaX) * correctFactor[kNonBending];// x position (info from non-bending plane) 
	    trackIntersectCh[currCh][1] = y11 + deltaZ * TMath::Tan(thetaY) * correctFactor[kBending];// y position (info from bending plane)
	    Int_t detElemIdFromTrack = DetElemIdFromPos(trackIntersectCh[currCh][0], trackIntersectCh[currCh][1], 11+currCh, 0);
	    if(detElemIdFromTrack<0) {
		if(fDebugLevel>1) printf("Warning: trigger track outside trigger chamber\n");
		continue;
	    }
		
	    triggeredDigits[1][currCh] = MatchingPad(digitStore, detElemIdFromTrack, trackIntersectCh[currCh], matchPad, nboard, zRealMatch, y11);

	    // if MatchingPad = -500 => too many digits matching pad =>
	    //                       => Event not clear => Reject track
	    if(triggeredDigits[1][currCh]<-100){
		isClearEvent = kFALSE;
		if(fDebugLevel>=1) printf("Warning: track = %p (%i) matches many pads. Rejected!\n",recTrigTrack, detElemIdFromTrack);
		break;
	    }

	    // deciding if digit matches track
	    Bool_t isDiffLocBoard = kFALSE;
	    if(fReproduceTrigResponse && ch>2){
		for(Int_t cath=0; cath<fgkNcathodes; cath++){
		    if(boardThatTriggered[currCh][cath][0]>=0){
			if(boardThatTriggered[currCh][cath][0]!=boardThatTriggered[currCh-1][cath][0]) isDiffLocBoard = kTRUE;
		    }
		}
	    }

	    if(isDiffLocBoard && fDebugLevel>=1)printf("\tDifferent local board\n");

	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		match[currCh][cath] = (matchPad[cath] && !isDiffLocBoard);
		if(!match[currCh][cath]) continue;
		digitPerTrack[cath]++;
		trigScheme[currCh][cath]++;
		slatThatTriggered[currCh][cath] = detElemIdFromTrack;
		for(Int_t loc=0; loc<4; loc++){
		    boardThatTriggered[currCh][cath][loc] = nboard[cath][loc];
		}
	    }
	} // end chamber loop

	for(Int_t cath=0; cath<fgkNcathodes; cath++){
	    if(digitPerTrack[cath]<3)isClearEvent = kFALSE;
	    if(fDebugLevel>=1 && !isClearEvent)printf("Warning: found %i digits for trigger track cathode %i.\nRejecting event\n", digitPerTrack[cath],cath);
	}

	if(!isClearEvent && !fReproduceTrigResponse) continue;

	Int_t commonDigits = 0;
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    if(triggeredDigits[1][ch]==triggeredDigits[0][ch]) commonDigits++; // Compare with previous track
	    triggeredDigits[0][ch] = triggeredDigits[1][ch]; // Store this track parameters for comparison with next one
	}
	if(commonDigits>=2){
	    doubleCountTrack=kTRUE;
	}

	if(!doubleCountTrack || fReproduceTrigResponse){
	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		Int_t is44 = 1;
		Bool_t goodForSlatEff = kTRUE;
		Bool_t goodForBoardEff = kTRUE;
		Int_t ineffSlat = -1;
		Int_t ineffDetElId = -1;
		Int_t firstSlat = slatThatTriggered[0][cath]%100;
		if(firstSlat<0) firstSlat=slatThatTriggered[1][cath]%100;
		Int_t firstBoard = boardThatTriggered[0][kBending][0];
		if(firstBoard<0) firstBoard=boardThatTriggered[1][kBending][0];
		for(Int_t ch=0; ch<fgkNchambers; ch++){
		    Bool_t isCurrChIneff = kFALSE;
		    is44 *= trigScheme[ch][cath];
		    Int_t currSlat = slatThatTriggered[ch][cath]%100;
		    if(currSlat<0){
			ineffDetElId = DetElemIdFromPos(trackIntersectCh[ch][0], trackIntersectCh[ch][1], 11+ch, cath);
			currSlat = ineffDetElId%100;
			ineffSlat = currSlat;
			isCurrChIneff = kTRUE;
		    }
		    if(currSlat!=firstSlat)goodForSlatEff=kFALSE;
		    Bool_t atLeastOneLoc=kFALSE;
		    if(isCurrChIneff) LocalBoardFromPos(trackIntersectCh[ch][0], trackIntersectCh[ch][1], ineffDetElId, cath, ineffBoard);
		    for(Int_t loc=0; loc<4; loc++){
			Int_t currBoard = boardThatTriggered[ch][cath][loc];
			if(isCurrChIneff) currBoard = ineffBoard[loc];
			if(currBoard==firstBoard){
			    atLeastOneLoc=kTRUE;
			    break;
			}
		    }
		    if(!atLeastOneLoc)goodForBoardEff=kFALSE;
		} // end chamber loop

		// Trigger 4/4
		if(is44==1){
		    fTrigger44[cath]++;
		    if(fDebugLevel>=1)printf("Trigger44[%i] = %i\n",cath,fTrigger44[cath]);
		    if(goodForSlatEff){
			for(Int_t ch=0; ch<fgkNchambers; ch++){
			    fHitPerSlat[ch][cath][firstSlat]++;
			    if(fDebugLevel>=1)printf("Slat that triggered = %i\n",slatThatTriggered[ch][cath]);
			    if(goodForBoardEff && firstBoard>0){
				fHitPerBoard[ch][cath][firstBoard-1]++;
				if(fDebugLevel>=1)printf("Board that triggered = %i\n",firstBoard);
			    }
			    else if(fDebugLevel>=1) printf("Track = %p: Particle crossed different boards: rejected!\n",recTrigTrack);
			}
		    }
		    else if(fDebugLevel>=1) printf("Track = %p: Particle crossed different slats: rejected!\n",recTrigTrack);
		    //cout<<"fTrigger44["<<cath<<"] = "<<fTrigger44[cath]<<"\tfHitPerSlat["<<0<<"]["<<cath<<"]["<<firstSlat<<"] = "<<fHitPerSlat[0][cath][firstSlat]<<"\tfHitPerBoard["<<0<<"]["<<cath<<"]["<<firstBoard-1<<"] = "<<fHitPerBoard[0][cath][firstBoard-1]<<endl; //REMEMBER TO CUT
		}

		// Trigger 3/4
		if(ineffDetElId>0){
		    Int_t ineffCh = ineffDetElId/100-11;
		    fTrigger34[ineffCh][cath]++;
		    if(fDebugLevel>=1) printf("Trigger34[%i][%i] = %i\n",ineffCh,cath,fTrigger34[ineffCh][cath]);
		    if(goodForSlatEff){
			if(fDebugLevel>=1) printf("Slat non efficient = %i\n",ineffDetElId);
			fInefficientSlat[ineffCh][cath][ineffSlat]++;

			if(goodForBoardEff && firstBoard>0){
			    if(fDebugLevel>=1) printf("Board non efficient = %i\n",firstBoard);
			    fInefficientBoard[ineffCh][cath][firstBoard-1]++;
			}
			else if(fDebugLevel>=1) printf("Track = %p: Particle crossed different boards: rejected!\n",recTrigTrack);
		    }
		    else if(fDebugLevel>=1) printf("Track = %p: Particle crossed different slats: rejected!\n",recTrigTrack);
		    //cout<<"fTrigger34["<<ineffCh<<"]["<<cath<<"] = "<<fTrigger34[ineffCh][cath]<<"\tfInefficientSlat["<<ineffCh<<"]["<<cath<<"]["<<ineffSlat<<"] = "<<fInefficientSlat[ineffCh][cath][ineffSlat]<<"\tfInefficientBoard["<<ineffCh<<"]["<<cath<<"]["<<firstBoard-1<<"] = "<<fInefficientBoard[ineffCh][cath][firstBoard-1]<<endl; //REMEMBER TO CUT
		}
	    } // end loop on cathodes
	}
	else if(doubleCountTrack){
	    if(fDebugLevel>=1)
		printf("\n\tTrack = %p: \nDouble Count Track: Track rejected!\n",recTrigTrack);
	}
    } // end trigger tracks loop

    if(fPrintInfo) InfoDigit(digitStore);
}

//_____________________________________________________________________________
void AliMUONTriggerChamberEff::WriteEfficiencyMap(const char* outputDir)
{
    //
    /// Writes information on calculated efficiency.
    /// It writes: triggerChamberEff.root file containing efficiency histograms.
    //

    char *cathCode[fgkNcathodes] = {"bendPlane", "nonBendPlane"};

    char outFileName[100];

    sprintf(outFileName, "%s/MUON.TriggerEfficiencyMap.root",outputDir);
    TFile *outputHistoFile = new TFile(outFileName,"RECREATE");
    TDirectory *dir = gDirectory;

    enum {kSlatIn11, kSlatIn12, kSlatIn13, kSlatIn14, kChamberEff};
    char *yAxisTitle = "trigger efficiency (a.u.)";
    char *xAxisTitle = "chamber";

    TH1F *histo[fgkNcathodes][fgkNchambers+1];
    TH1F *histoBoard[fgkNcathodes][fgkNchambers];

    // ADDED for check
    enum {allChEff, chNonEff, numOfHistoTypes};
    char *histoTypeName[numOfHistoTypes] = {"CountInCh", "NonCountInCh"};
    char *histoTypeTitle[numOfHistoTypes] = {"counted", "non counted"};
    TH1F *histoCheckSlat[fgkNcathodes][fgkNchambers][numOfHistoTypes];
    TH1F *histoCheckBoard[fgkNcathodes][fgkNchambers][numOfHistoTypes];
    // end ADDED for check

    char histoName[40];
    char histoTitle[90];

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	for(Int_t ch=0; ch<fgkNchambers+1; ch++){
	    if(ch==kChamberEff){
		sprintf(histoName, "%sChamberEff", cathCode[cath]);
		sprintf(histoTitle, "Chamber efficiency %s", cathCode[cath]);
		histo[cath][ch] = new TH1F(histoName, histoTitle, fgkNchambers, 11-0.5, 15-0.5);
		histo[cath][ch]->SetXTitle(xAxisTitle);
		histo[cath][ch]->SetYTitle(yAxisTitle);
		histo[cath][ch]->GetXaxis()->SetNdivisions(fgkNchambers);
	    }
	    else {
		sprintf(histoName, "%sSlatEffChamber%i", cathCode[cath], 11+ch);
		sprintf(histoTitle, "Chamber %i: slat efficiency %s", 11+ch, cathCode[cath]);
		histo[cath][ch] = new TH1F(histoName, histoTitle, fgkNslats, 0-0.5, fgkNslats-0.5);
		histo[cath][ch]->SetXTitle("slat");
		histo[cath][ch]->SetYTitle(yAxisTitle);
		histo[cath][ch]->GetXaxis()->SetNdivisions(fgkNslats);
		
		sprintf(histoName, "%sBoardEffChamber%i", cathCode[cath], 11+ch);
		sprintf(histoTitle, "Chamber %i: board efficiency %s", 11+ch, cathCode[cath]);
		histoBoard[cath][ch] = new TH1F(histoName, histoTitle, fgkNboards, 1-0.5, fgkNboards+1.-0.5);
		histoBoard[cath][ch]->SetXTitle("boards");
		histoBoard[cath][ch]->SetYTitle(yAxisTitle);
		histoBoard[cath][ch]->GetXaxis()->SetNdivisions(fgkNboards);

		// ADDED for check
		for(Int_t hType=0; hType<numOfHistoTypes; hType++){
		    sprintf(histoName, "%sSlat%s%i", cathCode[cath], histoTypeName[hType], 11+ch);
		    sprintf(histoTitle, "Chamber %i: slat %s %s", 11+ch, histoTypeTitle[hType], cathCode[cath]);
		    histoCheckSlat[cath][ch][hType] = new TH1F(histoName, histoTitle, fgkNslats, 0-0.5, fgkNslats-0.5);
		    histoCheckSlat[cath][ch][hType]->SetXTitle("slat");
		    histoCheckSlat[cath][ch][hType]->SetYTitle(yAxisTitle);
		    histoCheckSlat[cath][ch][hType]->GetXaxis()->SetNdivisions(fgkNslats);

		    sprintf(histoName, "%sBoard%s%i", cathCode[cath], histoTypeName[hType], 11+ch);
		    sprintf(histoTitle, "Chamber %i: board %s %s", 11+ch, histoTypeTitle[hType], cathCode[cath]);
		    histoCheckBoard[cath][ch][hType] = new TH1F(histoName, histoTitle, fgkNboards, 1-0.5, fgkNboards+1.-0.5);
		    histoCheckBoard[cath][ch][hType]->SetXTitle("boards");
		    histoCheckBoard[cath][ch][hType]->SetYTitle(yAxisTitle);
		    histoCheckBoard[cath][ch][hType]->GetXaxis()->SetNdivisions(fgkNboards);
		}
		// end ADDED for check
	    }
	}
    }

    Float_t efficiency, efficiencyError;
    Int_t bin;

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    for(Int_t slat=0; slat<fgkNslats; slat++){
		CalculateEfficiency(fHitPerSlat[ch][cath][slat], fHitPerSlat[ch][cath][slat]+fInefficientSlat[ch][cath][slat], efficiency, efficiencyError, kFALSE);
		bin = histo[cath][ch]->FindBin(slat);
		histo[cath][ch]->SetBinContent(bin, efficiency);
		histo[cath][ch]->SetBinError(bin, efficiencyError);

		// ADDED for check
		histoCheckSlat[cath][ch][allChEff]->SetBinContent(bin, fHitPerSlat[ch][cath][slat]);
		histoCheckSlat[cath][ch][chNonEff]->SetBinContent(bin, fInefficientSlat[ch][cath][slat]);
	    }
	    CalculateEfficiency(fTrigger44[cath], fTrigger34[ch][cath]+fTrigger44[cath], efficiency, efficiencyError, kFALSE);
	    bin = histo[cath][ch]->FindBin(11+ch);
	    histo[cath][kChamberEff]->SetBinContent(bin, efficiency);
	    histo[cath][kChamberEff]->SetBinError(bin, efficiencyError);

	    for(Int_t board=0; board<fgkNboards; board++){
		CalculateEfficiency(fHitPerBoard[ch][cath][board], fHitPerBoard[ch][cath][board]+fInefficientBoard[ch][cath][board], efficiency, efficiencyError, kFALSE);
		bin = histoBoard[cath][ch]->FindBin(board+1);
		histoBoard[cath][ch]->SetBinContent(bin, efficiency);
		histoBoard[cath][ch]->SetBinError(bin, efficiencyError);

		// ADDED for check
		histoCheckBoard[cath][ch][allChEff]->SetBinContent(bin, fHitPerBoard[ch][cath][board]);
		histoCheckBoard[cath][ch][chNonEff]->SetBinContent(bin, fInefficientBoard[ch][cath][board]);
	    }
	}
    }

// write all histos
    outputHistoFile->cd();
    dir->GetList()->Write();
    outputHistoFile->Close();
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::WriteEfficiencyMapTxt(const char* outputDir)
{
    //
    /// Writes the calculated efficiency in the text file efficiencyCells.dat
    ///
    /// The file can be further put in $ALICE_ROOT/MUON/data
    /// and used to run simulations with measured trigger chamber efficiencies.
    //

    Int_t effOutWidth=4;

    Float_t efficiency, efficiencyError;

    Int_t aCapo[] = {16, 38, 60, 76, 92, 108, 117, 133, 155, 177, 193, 209, 225, 234};

    char filename[70];
    sprintf(filename, "%s/efficiencyCells.dat", outputDir);
    ofstream outFile(filename);
    outFile << "localBoards" << endl;
    for(Int_t ch=0; ch<fgkNchambers; ch++){
	//Print information
	outFile << "\n\ndetElemId:\t" << 11+ch;
	outFile << "00";
	for(Int_t cath=0; cath<fgkNcathodes; cath++){
	    outFile << "\n cathode:\t" << cath << endl;
	    Int_t currLine=0;
	    for(Int_t board=0; board<fgkNboards; board++){

		if(board==aCapo[currLine]){
		    outFile << endl;
		    currLine++;
		}
		CalculateEfficiency(fHitPerBoard[ch][cath][board], fHitPerBoard[ch][cath][board]+fInefficientBoard[ch][cath][board], efficiency, efficiencyError, kFALSE);
		outFile << " " << setw(effOutWidth) << efficiency;
	    }// loop on boards
	    outFile << endl;
	}// loop on cathodes
    }// loop on chambers    
}


//_____________________________________________________________________________
Bool_t AliMUONTriggerChamberEff::IsCleanTrack(AliMUONTriggerTrack *triggerTrack,
					      const AliMUONVTrackStore& trackStore)
{
    //
    /// Try to match track from tracking system with trigger track
    //
    const Double_t kDistSigma[3]={1,1,0.02}; // sigma of distributions (trigger-track) X,Y,slopeY
    const Double_t kMaxChi2MatchTrigger = 16.0;
  
    AliMUONTrackParam trackParam; 

    Double_t distTriggerTrack[3];
    Double_t xTrack, yTrack, ySlopeTrack, chi2;
  
    AliMUONTrack* track;
    TIter next(trackStore.CreateIterator());
    
    while ( ( track = static_cast<AliMUONTrack*>(next()) ) )
    {
	trackParam = *((AliMUONTrackParam*) (track->GetTrackParamAtHit()->Last()));
	AliMUONTrackExtrap::ExtrapToZ(&trackParam, AliMUONConstants::DefaultChamberZ(10)); // extrap to 1st trigger chamber
    
	xTrack = trackParam.GetNonBendingCoor();
	yTrack = trackParam.GetBendingCoor();
	ySlopeTrack = trackParam.GetBendingSlope();
  
	distTriggerTrack[0] = (triggerTrack->GetX11()-xTrack)/kDistSigma[0];
	distTriggerTrack[1] = (triggerTrack->GetY11()-yTrack)/kDistSigma[1];
	distTriggerTrack[2] = (TMath::Tan(triggerTrack->GetThetay())-ySlopeTrack)/kDistSigma[2];
	chi2 = 0.;
	for (Int_t iVar = 0; iVar < 3; iVar++) chi2 += distTriggerTrack[iVar]*distTriggerTrack[iVar];
	chi2 /= 3.; // Normalized Chi2: 3 degrees of freedom (X,Y,slopeY)
	if (chi2 < kMaxChi2MatchTrigger) return kTRUE;
    }

    return kFALSE;
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::SaveInESDFile()
{
    //
    /// Store AliMUONTriggerChamberEff in esd file
    //
    TDirectory *dir = gDirectory;
    TFile *logFile = 0x0;
    TSeqCollection *list = gROOT->GetListOfFiles();
    Int_t n = list->GetEntries();
    for(Int_t i=0; i<n; i++) {
	logFile = (TFile*)list->At(i);
	if (strstr(logFile->GetName(), "AliESDs.root")) break;
    }
    if(logFile){
	TTree *esdTree = (TTree*)logFile->Get("esdTree");
	if(esdTree){
	    if(!esdTree->GetUserInfo()->FindObject("AliMUONTriggerChamberEff")){
		AliInfo(Form("Adding AliMUONTrigChamberEff in %s",logFile->GetName()));
		esdTree->GetUserInfo()->Add(this->Clone());
		esdTree->Write("",TObject::kOverwrite);
	    }
	}
    }
    dir->cd();
}
