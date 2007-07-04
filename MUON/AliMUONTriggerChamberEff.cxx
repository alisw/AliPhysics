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
#include "AliMpConstants.h"

#include "AliMpDDLStore.h"
#include "AliMpDDL.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"

#include "AliLog.h"

#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>

#include <TSeqCollection.h>
#include <TTree.h>
#include <TROOT.h>

#include <TStyle.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TSystem.h>
#include <TPaveLabel.h>

#include <cassert>

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
/// Default constructor

    CheckConstants();
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

    CheckConstants();
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
     fkMaxDistance(other.fkMaxDistance),
     fTrigger44(other.fTrigger44),
     fTrigger34(other.fTrigger34)
{
/// Copy constructor

    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
	fInefficientSlat[chCath] = other.fInefficientSlat[chCath];
	fHitPerSlat[chCath] = other.fHitPerSlat[chCath];
	fInefficientBoard[chCath] = other.fInefficientBoard[chCath];
	fHitPerBoard[chCath] = other.fHitPerBoard[chCath];
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

    fTrigger44 = other.fTrigger44;
    fTrigger34 = other.fTrigger34;

    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
	fInefficientSlat[chCath] = other.fInefficientSlat[chCath];
	fHitPerSlat[chCath] = other.fHitPerSlat[chCath];
	fInefficientBoard[chCath] = other.fInefficientBoard[chCath];
	fHitPerBoard[chCath] = other.fHitPerBoard[chCath];
    }
    return *this;
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::ResetArrays()
{
    //
    /// Sets the data member counters to 0.
    //

    fTrigger44.Set(fgkNcathodes);
    fTrigger34.Set(fgkNplanes);

    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
	fInefficientSlat[chCath].Set(fgkNslats);
	fHitPerSlat[chCath].Set(fgkNslats);
	fInefficientBoard[chCath].Set(AliMpConstants::NofLocalBoards());
	fHitPerBoard[chCath].Set(AliMpConstants::NofLocalBoards());
    }

    fTrigger44.Reset();
    fTrigger34.Reset();

    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
	fInefficientSlat[chCath].Reset();
	fHitPerSlat[chCath].Reset();
	fInefficientBoard[chCath].Reset();
	fHitPerBoard[chCath].Reset();
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
   
	TArrayS xyPattern[fgkNcathodes];
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
					    Float_t coor[2], Bool_t isMatch[2],
					    TArrayI nboard[2], TArrayF &zRealMatch, Float_t y11)
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
	for(Int_t loc=pad.GetNofLocations(); loc<fgkNlocations; loc++){
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

    for(Int_t loc=0; loc<fgkNlocations; loc++){
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

    TArrayF zRealMatch(fgkNchambers);
    TArrayF correctFactor(fgkNcathodes);

    Bool_t match[fgkNchambers][fgkNcathodes];
    Bool_t matchPad[fgkNcathodes];
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	matchPad[cath] = kFALSE;
    }

    TArrayF zMeanChamber(fgkNchambers);
    for(Int_t ch=0; ch<fgkNchambers; ch++){
	zMeanChamber[ch] = AliMUONConstants::DefaultChamberZ(10+ch);
    }

    TArrayI digitPerTrack(fgkNcathodes);

    Float_t trackIntersectCh[fgkNchambers][fgkNcathodes];

    TArrayI triggeredDigits[2];
    for(Int_t itrack=0; itrack<2; itrack++){
	triggeredDigits[itrack].Set(fgkNchambers);
	triggeredDigits[itrack].Reset(-1);
    }

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

    AliMUONDigitStoreV1 digitStore;   
    TriggerDigits(triggerStore,digitStore);

    AliMUONTriggerTrack *recTrigTrack = 0x0;

    TIter next(trigTrackStore.CreateIterator());
    
    while ( ( recTrigTrack = static_cast<AliMUONTriggerTrack*>(next()) ) )
    {
	if(!IsCleanTrack(recTrigTrack, trackStore)) {
	    if(fDebugLevel>=1) printf("\tTrack %p (%f, %f) don't match tracker track: rejected!\n",(void *)recTrigTrack,recTrigTrack->GetX11(),recTrigTrack->GetY11());
	    continue;
	}

	digitPerTrack.Reset();
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    zRealMatch[ch] = zMeanChamber[ch];
	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		match[ch][cath]=kFALSE;
		//slatThatTriggered[ch][cath]=-1;
		//trigScheme[ch][cath] = 0;
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
	Bool_t doubleCountTrack = kFALSE;

	Float_t x11 = recTrigTrack->GetX11();// x position (info from non-bending plane)
	Float_t y11 = recTrigTrack->GetY11();// y position (info from bending plane)
	Float_t thetaX = recTrigTrack->GetThetax();
	Float_t thetaY = recTrigTrack->GetThetay();

	if(fDebugLevel>=3) printf("\tTrack = %p\npos from track: (x,y) = (%f, %f), (thetaX, thetaY) = (%f, %f)\n",(void *)recTrigTrack,x11,y11,thetaX*rad2deg,thetaY*rad2deg);

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
		if(fDebugLevel>=1) printf("Warning: track = %p (%i) matches many pads. Rejected!\n",(void *)recTrigTrack, detElemIdFromTrack);
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
		trigScheme[cath][currCh]++;
		slatThatTriggered[cath][currCh] = detElemIdFromTrack;
		for(Int_t loc=0; loc<fgkNlocations; loc++){
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
		Int_t firstSlat = slatThatTriggered[cath][0]%100;
		if(firstSlat<0) firstSlat=slatThatTriggered[cath][1]%100;
		Int_t firstBoard = boardThatTriggered[0][kBending][0];
		if(firstBoard<0) firstBoard=boardThatTriggered[1][kBending][0];
		for(Int_t ch=0; ch<fgkNchambers; ch++){
		    Bool_t isCurrChIneff = kFALSE;
		    is44 *= trigScheme[cath][ch];
		    Int_t currSlat = slatThatTriggered[cath][ch]%100;
		    if(currSlat<0){
			ineffDetElId = DetElemIdFromPos(trackIntersectCh[ch][0], trackIntersectCh[ch][1], 11+ch, cath);
			currSlat = ineffDetElId%100;
			ineffSlat = currSlat;
			isCurrChIneff = kTRUE;
		    }
		    if(currSlat!=firstSlat)goodForSlatEff=kFALSE;
		    Bool_t atLeastOneLoc=kFALSE;
		    if(isCurrChIneff) LocalBoardFromPos(trackIntersectCh[ch][0], trackIntersectCh[ch][1], ineffDetElId, cath, ineffBoard);
		    for(Int_t loc=0; loc<fgkNlocations; loc++){
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
			    Int_t chCath = fgkNchambers*cath + ch;
			    fHitPerSlat[chCath][firstSlat]++;
			    if(fDebugLevel>=1)printf("Slat that triggered = %i\n",slatThatTriggered[cath][ch]);
			    if(goodForBoardEff && firstBoard>0){
				fHitPerBoard[chCath][firstBoard-1]++;
				if(fDebugLevel>=1)printf("Board that triggered = %i\n",firstBoard);
			    }
			    else if(fDebugLevel>=1) printf("Track = %p: Particle crossed different boards: rejected!\n",(void *)recTrigTrack);
			}
		    }
		    else if(fDebugLevel>=1) printf("Track = %p: Particle crossed different slats: rejected!\n",(void *)recTrigTrack);
		}

		// Trigger 3/4
		if(ineffDetElId>0){
		    Int_t ineffCh = ineffDetElId/100-11;
		    Int_t chCath = fgkNchambers*cath + ineffCh;
		    fTrigger34[chCath]++;
		    if(fDebugLevel>=1) printf("Trigger34[%i] = %i\n",chCath,fTrigger34[chCath]);
		    if(goodForSlatEff){
			if(fDebugLevel>=1) printf("Slat non efficient = %i\n",ineffDetElId);
			fInefficientSlat[chCath][ineffSlat]++;

			if(goodForBoardEff && firstBoard>0){
			    if(fDebugLevel>=1) printf("Board non efficient = %i\n",firstBoard);
			    fInefficientBoard[chCath][firstBoard-1]++;
			}
			else if(fDebugLevel>=1) printf("Track = %p: Particle crossed different boards: rejected!\n",(void *)recTrigTrack);
		    }
		    else if(fDebugLevel>=1) printf("Track = %p: Particle crossed different slats: rejected!\n",(void *)recTrigTrack);
		}
	    } // end loop on cathodes
	}
	else if(doubleCountTrack){
	    if(fDebugLevel>=1)
		printf("\n\tTrack = %p: \nDouble Count Track: Track rejected!\n",(void *)recTrigTrack);
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

    char *yAxisTitle = "trigger efficiency (a.u.)";
    char *xAxisTitle = "chamber";

    const Int_t kNumOfBoards = AliMpConstants::NofLocalBoards();

    TH1F *histoChamber[fgkNcathodes];
    TH1F *histoSlat[8];
    TH1F *histoBoard[8];

    // ADDED for check
    enum {kAllChEff, kChNonEff, kNumOfHistoTypes};
    char *histoTypeName[kNumOfHistoTypes] = {"CountInCh", "NonCountInCh"};
    char *histoTypeTitle[kNumOfHistoTypes] = {"counted", "non counted"};
    TH1F *histoCheckSlat[8][kNumOfHistoTypes];
    TH1F *histoCheckBoard[8][kNumOfHistoTypes];
    // end ADDED for check

    char histoName[40];
    char histoTitle[90];

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	sprintf(histoName, "%sChamberEff", cathCode[cath]);
	sprintf(histoTitle, "Chamber efficiency %s", cathCode[cath]);
	histoChamber[cath] = new TH1F(histoName, histoTitle, fgkNchambers, 11-0.5, 15-0.5);
	histoChamber[cath]->SetXTitle(xAxisTitle);
	histoChamber[cath]->SetYTitle(yAxisTitle);
	histoChamber[cath]->GetXaxis()->SetNdivisions(fgkNchambers);
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    Int_t chCath = fgkNchambers*cath + ch;
	    sprintf(histoName, "%sSlatEffChamber%i", cathCode[cath], 11+ch);
	    sprintf(histoTitle, "Chamber %i: slat efficiency %s", 11+ch, cathCode[cath]);
	    histoSlat[chCath] = new TH1F(histoName, histoTitle, fgkNslats, 0-0.5, fgkNslats-0.5);
	    histoSlat[chCath]->SetXTitle("slat");
	    histoSlat[chCath]->SetYTitle(yAxisTitle);
	    histoSlat[chCath]->GetXaxis()->SetNdivisions(fgkNslats);
		
	    sprintf(histoName, "%sBoardEffChamber%i", cathCode[cath], 11+ch);
	    sprintf(histoTitle, "Chamber %i: board efficiency %s", 11+ch, cathCode[cath]);
	    histoBoard[chCath] = new TH1F(histoName, histoTitle, kNumOfBoards, 1-0.5, kNumOfBoards+1.-0.5);
	    histoBoard[chCath]->SetXTitle("boards");
	    histoBoard[chCath]->SetYTitle(yAxisTitle);
	    histoBoard[chCath]->GetXaxis()->SetNdivisions(kNumOfBoards);

	    // ADDED for check
	    for(Int_t hType=0; hType<kNumOfHistoTypes; hType++){
		sprintf(histoName, "%sSlat%s%i", cathCode[cath], histoTypeName[hType], 11+ch);
		sprintf(histoTitle, "Chamber %i: slat %s %s", 11+ch, histoTypeTitle[hType], cathCode[cath]);
		histoCheckSlat[chCath][hType] = new TH1F(histoName, histoTitle, fgkNslats, 0-0.5, fgkNslats-0.5);
		histoCheckSlat[chCath][hType]->SetXTitle("slat");
		histoCheckSlat[chCath][hType]->SetYTitle(yAxisTitle);
		histoCheckSlat[chCath][hType]->GetXaxis()->SetNdivisions(fgkNslats);

		sprintf(histoName, "%sBoard%s%i", cathCode[cath], histoTypeName[hType], 11+ch);
		sprintf(histoTitle, "Chamber %i: board %s %s", 11+ch, histoTypeTitle[hType], cathCode[cath]);
		histoCheckBoard[chCath][hType] = new TH1F(histoName, histoTitle, kNumOfBoards, 1-0.5, kNumOfBoards+1.-0.5);
		histoCheckBoard[chCath][hType]->SetXTitle("boards");
		histoCheckBoard[chCath][hType]->SetYTitle(yAxisTitle);
		histoCheckBoard[chCath][hType]->GetXaxis()->SetNdivisions(kNumOfBoards);
	    }
	    // end ADDED for check
	}
    }

    Float_t efficiency, efficiencyError;
    Int_t bin;

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    Int_t chCath = fgkNchambers*cath + ch;
	    for(Int_t slat=0; slat<fgkNslats; slat++){
		CalculateEfficiency(fHitPerSlat[chCath][slat], fHitPerSlat[chCath][slat]+fInefficientSlat[chCath][slat], efficiency, efficiencyError, kFALSE);
		bin = histoSlat[chCath]->FindBin(slat);
		histoSlat[chCath]->SetBinContent(bin, efficiency);
		histoSlat[chCath]->SetBinError(bin, efficiencyError);

		// ADDED for check
		histoCheckSlat[chCath][kAllChEff]->SetBinContent(bin, fHitPerSlat[chCath][slat]);
		histoCheckSlat[chCath][kChNonEff]->SetBinContent(bin, fInefficientSlat[chCath][slat]);
	    }
	    CalculateEfficiency(fTrigger44[cath], fTrigger34[chCath]+fTrigger44[cath], efficiency, efficiencyError, kFALSE);
	    bin = histoChamber[cath]->FindBin(11+ch);
	    histoChamber[cath]->SetBinContent(bin, efficiency);
	    histoChamber[cath]->SetBinError(bin, efficiencyError);

	    for(Int_t board=0; board<kNumOfBoards; board++){
		CalculateEfficiency(fHitPerBoard[chCath][board], fHitPerBoard[chCath][board]+fInefficientBoard[chCath][board], efficiency, efficiencyError, kFALSE);
		bin = histoBoard[chCath]->FindBin(board+1);
		histoBoard[chCath]->SetBinContent(bin, efficiency);
		histoBoard[chCath]->SetBinError(bin, efficiencyError);

		// ADDED for check
		histoCheckBoard[chCath][kAllChEff]->SetBinContent(bin, fHitPerBoard[chCath][board]);
		histoCheckBoard[chCath][kChNonEff]->SetBinContent(bin, fInefficientBoard[chCath][board]);
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
	    Int_t chCath = fgkNchambers*cath + ch;
	    Int_t currLine=0;
	    for(Int_t board=0; board<AliMpConstants::NofLocalBoards(); board++){

		if(board==aCapo[currLine]){
		    outFile << endl;
		    currLine++;
		}
		CalculateEfficiency(fHitPerBoard[chCath][board], fHitPerBoard[chCath][board]+fInefficientBoard[chCath][board], efficiency, efficiencyError, kFALSE);
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


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::DisplayEfficiency(Bool_t perSlat)
{
    //
    /// Display calculated efficiency.
    //

    const Int_t kNumOfBoards = AliMpConstants::NofLocalBoards();

    Int_t side, col, line, nbx, slat;
    Float_t xCenter, yCenter, zCenter, xWidth, yWidth;
    Int_t x1, y1, x2, y2, board=0;
    char name[8], text[200];

    gStyle->SetPalette(1);

    TString boardName[234];

    // instanciate the elec. mapping class
    AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();

    // loop over the trigger DDL (Right: 20, Left: 21)
    for (Int_t iDDL = 20; iDDL <= 21; ++iDDL) {

      // get ddl object
      AliMpDDL* ddl = ddlStore->GetDDL(iDDL);

      Int_t nCrate = ddl->GetNofTriggerCrates();
    
      // loop over the number of crate in DDL
      for (Int_t index = 0; index < nCrate; ++index) {

	// get crate object
	AliMpTriggerCrate* crate = ddlStore->GetTriggerCrate(iDDL, index);

	Int_t nLocal = crate->GetNofLocalBoards();

	for (Int_t iLocal = 0; iLocal  < nLocal; ++iLocal) {

	  // get local board Id from crate object
	  board = crate->GetLocalBoardId(iLocal);
	  if(board>kNumOfBoards || board<=0) continue;

	  AliMpLocalBoard* localBoard = ddlStore->GetLocalBoard(board);
	  boardName[board-1] = localBoard->GetName();
	}
      }
    }

    char *cathCode[fgkNcathodes] = {"bendPlane", "nonBendPlane"};

    Float_t boardsX = 257.00;  // cm
    Float_t boardsY = 307.00;  // cm

    TH2F *histo[8];
    TPaveLabel *boardLabel[8][234];

    char histoName[40];
    char histoTitle[90];

    Float_t efficiency, efficiencyError;

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    Int_t chCath = fgkNchambers*cath + ch;
	    sprintf(histoName, "%sChamber%i", cathCode[cath], 11+ch);
	    sprintf(histoTitle, "Chamber %i: efficiency %s", 11+ch, cathCode[cath]);
	    histo[chCath] = new TH2F(histoName, histoTitle, (Int_t)boardsX, -boardsX, boardsX, (Int_t)boardsY, -boardsY, boardsY);
	    histo[chCath]->SetXTitle("X (cm)");
	    histo[chCath]->SetYTitle("Y (cm)");
	}
    }

    TString mapspath = gSystem->Getenv("ALICE_ROOT");
    mapspath.Append("/MUON/data");

    sprintf(text,"%s/guimapp11.txt",mapspath.Data());
    FILE *fmap = fopen(text,"r");

    for (Int_t ib = 0; ib < kNumOfBoards; ib++) {
	fscanf(fmap,"%d   %d   %d   %d   %f   %f   %f   %f   %f   %s   \n",&side,&col,&line,&nbx,&xCenter,&yCenter,&xWidth,&yWidth,&zCenter,&name[0]);

	slat = (line-5)%fgkNslats;
	for(Int_t iBoard=0; iBoard<kNumOfBoards; iBoard++){
	    if(!boardName[iBoard].Contains(name)) continue;
	    board = iBoard;
	    break;
	}

	for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
	    x1 = histo[chCath]->GetXaxis()->FindBin(xCenter-xWidth/2.)+1;
	    y1 = histo[chCath]->GetYaxis()->FindBin(yCenter-yWidth/2.)+1;
	    x2 = histo[chCath]->GetXaxis()->FindBin(xCenter+xWidth/2.)-1;
	    y2 = histo[chCath]->GetYaxis()->FindBin(yCenter+yWidth/2.)-1;
	    
	    boardLabel[chCath][board] = new TPaveLabel(xCenter-xWidth/2., yCenter-yWidth/2., xCenter+xWidth/2., yCenter+yWidth/2., Form("%3d",board+1));
	    boardLabel[chCath][board]->SetFillStyle(0);
	    boardLabel[chCath][board]->SetBorderSize(0);

	    if(!perSlat){
		CalculateEfficiency(fHitPerBoard[chCath][board], fHitPerBoard[chCath][board]+fInefficientBoard[chCath][board], efficiency, efficiencyError, kFALSE);
	    }
	    else{
		CalculateEfficiency(fHitPerSlat[chCath][slat], fHitPerSlat[chCath][slat]+fInefficientSlat[chCath][slat], efficiency, efficiencyError, kFALSE);
	    }
	    
	    for(Int_t binX=x1; binX<=x2; binX++){
		for(Int_t binY=y1; binY<=y2; binY++){
		    histo[chCath]->SetBinContent(binX, binY, efficiency);
		    histo[chCath]->SetBinError(binX, binY, efficiencyError);
		}
	    }
	}
    }

    TCanvas *can[8];
    for(Int_t chCath=0; chCath<fgkNplanes; chCath++){
	sprintf(histoName, "%sCan", histo[chCath]->GetName());
	sprintf(histoTitle, "%s", histo[chCath]->GetTitle());
	can[chCath] = new TCanvas(histoName, histoTitle, 100+10*chCath, 10*chCath, 700, 700);
	can[chCath]->SetRightMargin(0.14);
	can[chCath]->SetLeftMargin(0.12);
	histo[chCath]->GetZaxis()->SetRangeUser(0.,1.);
	histo[chCath]->GetYaxis()->SetTitleOffset(1.4);
	histo[chCath]->SetStats(kFALSE);
	histo[chCath]->Draw("COLZ");
	for (Int_t board = 0; board < kNumOfBoards; board++) {
	    boardLabel[chCath][board]->Draw("same");
	}
    }
}

//_____________________________________________________________________________
void AliMUONTriggerChamberEff::CheckConstants() const
{
/// Check consistence of redefined constants 

  assert(fgkNcathodes == AliMpConstants::NofCathodes());    
  assert(fgkNchambers == AliMpConstants::NofTriggerChambers());    
  assert(fgkNplanes == AliMpConstants::NofTriggerChambers() * fgkNcathodes);    
}
