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
/// Author: Diego Stocco (Torino)

#include "AliMUONTriggerChamberEff.h"
#include "AliMUONDigit.h"
#include "AliMUONConstants.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONSegmentation.h"
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONTriggerTrack.h"

#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpPad.h"
#include "AliMpDEIterator.h"
#include "AliMpPlaneType.h"

#include "AliRunLoader.h"
#include "AliRun.h"

#include <Riostream.h>
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>

/// \cond CLASSIMP
ClassImp(AliMUONTriggerChamberEff)
/// \endcond

//_____________________________________________________________________________
AliMUONTriggerChamberEff::AliMUONTriggerChamberEff(const char* galiceFile, 
						   Int_t firstEvent, Int_t lastEvent) 
: TObject(),
  fFirstEvent(firstEvent),
  fLastEvent(lastEvent),
  fFirstRun(-1),
  fLastRun(-1),
  fRunLoader(0x0),
  fData(0x0),
  fReproduceTrigResponse(kFALSE),
  fPrintInfo(kFALSE),
  fMUON(0x0),
  fDebugLevel(0),
  fGaliceDir(0x0)
{
/// Standard constructor
    SetGaliceFile(galiceFile);
    ResetArrays();
}

//_____________________________________________________________________________
AliMUONTriggerChamberEff::AliMUONTriggerChamberEff(Int_t firstRun, Int_t lastRun,
						   const char* galiceRunDir, 
						   Int_t firstEvent, Int_t lastEvent) 
: TObject(),
  fFirstEvent(firstEvent),
  fLastEvent(lastEvent),
  fFirstRun(firstRun),
  fLastRun(lastRun),
  fRunLoader(0x0),
  fData(0x0),
  fReproduceTrigResponse(kFALSE),
  fPrintInfo(kFALSE),
  fMUON(0x0),
  fDebugLevel(0),
  fGaliceDir(galiceRunDir)
{
    ResetArrays();
    delete gAlice;
}

//_____________________________________________________________________________
AliMUONTriggerChamberEff::~AliMUONTriggerChamberEff()
{
/// Destructor
    fRunLoader->UnloadAll();
    delete fRunLoader;
    delete fData;
}

//_____________________________________________________________________________
void AliMUONTriggerChamberEff::SetGaliceFile(const char *galiceFile)
{
    //
    /// Opens the galice.root and loads tracks and digits.
    //

    fRunLoader = AliRunLoader::Open(galiceFile,"MUONFolder","READ");
    if (!fRunLoader) 
    {
	AliError(Form("Error opening %s file \n",galiceFile));
    }  
    else
    {
	fRunLoader->LoadgAlice();
	gAlice = fRunLoader->GetAliRun();
	fMUON = (AliMUON*)gAlice->GetModule("MUON");

	if(fLastEvent<=0 || fLastEvent>fRunLoader->GetNumberOfEvents())fLastEvent = fRunLoader->GetNumberOfEvents()-1;
	if(fFirstEvent<0)fFirstEvent=0;


	AliLoader* loader = fRunLoader->GetLoader("MUONLoader");
	if ( loader )
	{
	    fData = new AliMUONData(loader,"MUON","MUON");
	    loader->LoadTracks("READ");
	    loader->LoadDigits("READ");
	}
	else
	{
	    AliError(Form("Could get MUONLoader"));
	}
    }
}

//_____________________________________________________________________________
void AliMUONTriggerChamberEff::CleanGalice()
{
    //
    /// Unload all loaded data
    //
    
    fRunLoader->UnloadAll();
    delete fRunLoader;
    fRunLoader = 0;
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
	}
    }
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::InfoDigit()
{
    //
    /// Prints information on digits (for debugging)
    //

    AliMUONDigit * mDigit=0x0;
    Int_t firstTrigCh = AliMUONConstants::NTrackingCh();
    // Addressing
    Int_t nchambers = AliMUONConstants::NCh();
    fData->SetTreeAddress("D,GLT");
    
    fData->GetDigits();
    // Loop on chambers
    for(Int_t ichamber=firstTrigCh; ichamber<nchambers; ichamber++) {
      TClonesArray* digits = fData->Digits(ichamber);
      digits->Sort();
      Int_t ndigits = (Int_t)digits->GetEntriesFast();
      for(Int_t idigit=0; idigit<ndigits; idigit++) {
	  mDigit = (AliMUONDigit*)digits->At(idigit);
	  mDigit->Print();
      } // end digit loop
    } // end chamber loop
    fData->ResetDigits();
}


//_____________________________________________________________________________
Bool_t AliMUONTriggerChamberEff::PadMatchTrack(Float_t xPad, Float_t yPad, Float_t dpx, Float_t dpy, 
					       Float_t xTrackAtPad, Float_t yTrackAtPad, Int_t chamber)
{
    //
    /// Decides if the digit belongs to the trigger track.
    //

    Float_t numOfHalfWidth = 5.;
    Bool_t match = kFALSE;
    Float_t maxDistX = dpx;
    if(fReproduceTrigResponse && chamber>=2) maxDistX = 3.*dpx;// Non-bending plane: check the +- 1 strip between stations
    if(!fReproduceTrigResponse)maxDistX = numOfHalfWidth*dpx;
    Float_t maxDistY = dpy;
    if(fReproduceTrigResponse && chamber%2) maxDistY = 3*dpy;// bending plane: check the +- 1 strip between planes in the same station
    if(!fReproduceTrigResponse) maxDistY = numOfHalfWidth*dpy;
    Float_t deltaX = TMath::Abs(xPad-xTrackAtPad);
    Float_t deltaY = TMath::Abs(yPad-yTrackAtPad);
    if(deltaX<=maxDistX && deltaY<=maxDistY)match = kTRUE;
    return match;
}


//_____________________________________________________________________________
Bool_t AliMUONTriggerChamberEff::IsDiffLocalBoard(Int_t currDetElemId, Int_t iy, Int_t detElemIdP1, Int_t iyDigitP1) const
{
    //
    /// Determins if the digits belong to the same local board.
    /// Used only if one wants to reproduce the trigger algorithm result.
    /// (fReproduceTrigResponse = kTRUE).
    //

    Bool_t isDiff = kTRUE;
    if(detElemIdP1<0 || iyDigitP1<0)return kFALSE;
    Int_t currSlat = currDetElemId%100;
    Int_t slatP1 = detElemIdP1%100;
    Int_t currLoc = iy/16;
    Int_t locP1 = iyDigitP1/16;
    if(currSlat==slatP1 && currLoc==locP1)isDiff = kFALSE;
    return isDiff;
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::CalculateEfficiency(Int_t trigger44, Int_t trigger34,
						   Float_t &efficiency, Float_t &error, Bool_t failuresAsInput)
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
Int_t AliMUONTriggerChamberEff::DetElemIdFromPos(Float_t x, Float_t y, Int_t chamber, Int_t cathode)
{
    //
    /// Given the (x,y) position in the chamber,
    /// it returns the corresponding slat
    //

    Int_t resultingDetElemId = -1;
    AliMpDEIterator it;
    const AliMUONGeometryTransformer *kGeomTransformer = fMUON->GetGeometryTransformer();
    AliMUONSegmentation *segmentation = fMUON->GetSegmentation();
    for ( it.First(chamber-1); ! it.IsDone(); it.Next() ){
	Int_t detElemId = it.CurrentDE();

	if (  segmentation->HasDE(detElemId) ){
	    const AliMpVSegmentation* seg = 
		AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,cathode);
	    if (seg){
		Float_t deltax = seg->Dimensions().X();
		Float_t deltay = seg->Dimensions().Y();
		Float_t xlocal1 =  -deltax;
		Float_t ylocal1 =  -deltay;
		Float_t xlocal2 =  +deltax;
		Float_t ylocal2 =  +deltay;
		Float_t xg01, yg01, zg1, xg02, yg02, zg2;
		kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg01, yg01, zg1);
		kGeomTransformer->Local2Global(detElemId, xlocal2, ylocal2, 0, xg02, yg02, zg2);

		Float_t xg1 = xg01, xg2 = xg02, yg1 = yg01, yg2 = yg02;

		if(xg01>xg02){
		    xg1 = xg02;
		    xg2 = xg01;
		}
		if(yg01>yg02){
		    yg1 = yg02;
		    yg2 = yg01;
		}

		if(x>=xg1 && x<=xg2 && y>=yg1 && y<=yg2){
		    resultingDetElemId = detElemId;
		    break;
		}
	    }
	}
    }
    return resultingDetElemId;
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::LocalBoardFromPos(Float_t x, Float_t y, Int_t detElemId, Int_t cathode, Int_t localBoard[4])
{
    //
    /// Given the (x,y) position in the chamber,
    /// it returns the corresponding local board
    //

    for(Int_t loc=0; loc<4; loc++){
	localBoard[loc]=-1;
    }
    const AliMUONGeometryTransformer *kGeomTransformer = fMUON->GetGeometryTransformer();
    Float_t xl, yl, zl;
    kGeomTransformer->Global2Local(detElemId, x, y, 0, xl, yl, zl);
    TVector2 pos(xl,yl);
    const AliMpVSegmentation* seg = 
	AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,cathode);
    if (seg){
	AliMpPad pad = seg->PadByPosition(pos,kFALSE);
	for (Int_t loc=0; loc<pad.GetNofLocations(); loc++){
	    AliMpIntPair location = pad.GetLocation(loc);
	    localBoard[loc] = location.GetFirst();
	}
    }
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::PrintTrigger(AliMUONGlobalTrigger *globalTrig)
{
    //
    /// Print trigger response.
    //

    printf("===================================================\n");
    printf(" Global Trigger output\t \tLow pt\tHigh pt\n");

    printf(" number of Single:\t \t"); 
    printf("%i\t",globalTrig->SingleLpt());
    printf("%i\t",globalTrig->SingleHpt());
    printf("\n");

    printf(" number of UnlikeSign pair:\t"); 
    printf("%i\t",globalTrig->PairUnlikeLpt());
    printf("%i\t",globalTrig->PairUnlikeHpt());
    printf("\n");

    printf(" number of LikeSign pair:\t");  
    printf("%i\t",globalTrig->PairLikeLpt());
    printf("%i\t",globalTrig->PairLikeHpt());
    printf("\n");
    printf("===================================================\n");
    printf("\n");
}

void AliMUONTriggerChamberEff::PerformTriggerChamberEff(const char* outputDir)
{
    //
    /// Main method.
    /// It loops over the the trigger rec. tracks in the event.
    /// Then it search for matching digits around the track.
    /// Finally it calculates the efficiency for each trigger board.
    /// Files with calculated efficiency are placed in the user defined outputDir.
    //

    enum {kBending, kNonBending};
    Int_t evtBeforePrint = 1000;
    Float_t rad2deg = 180./TMath::Pi();

    Int_t chOrder[] = {0,2,1,3};
    Int_t station[] = {0,0,1,1};
    Float_t zRealMatch[fgkNchambers] = {0.0};
    Float_t correctFactor[fgkNcathodes] = {1.};

    Bool_t match[fgkNchambers][fgkNcathodes] = {{kFALSE}};

    TClonesArray *recTrigTracksArray = 0x0;
    AliMUONTriggerTrack *recTrigTrack = 0x0;
    AliMUONDigit * mDigit = 0x0;

    Float_t zMeanChamber[fgkNchambers];
    for(Int_t ch=0; ch<fgkNchambers; ch++){
	zMeanChamber[ch] = AliMUONConstants::DefaultChamberZ(10+ch);
    }

    TClonesArray * globalTrigger = 0x0;
    AliMUONGlobalTrigger * gloTrg = 0x0; 

    Int_t partNumOfTrig[fgkNchambers][fgkNcathodes] = {{0}};
    Int_t totNumOfTrig[fgkNchambers][fgkNcathodes] = {{0}};
    Int_t atLeast1MuPerEv[fgkNchambers][fgkNcathodes] = {{0}};
    Int_t digitPerTrack[fgkNcathodes] = {0};

    Float_t trackIntersectCh[2][fgkNchambers]={{0.0}};

    Int_t slatInPlane1[2][fgkNcathodes];
    Int_t iyDigitInPlane1[2][fgkNcathodes];

    const Int_t kMaxNumOfTracks = 10;
    Int_t trigScheme[kMaxNumOfTracks][fgkNchambers][fgkNcathodes]={{{0}}};
    Int_t triggeredDigits[kMaxNumOfTracks][fgkNchambers][fgkNcathodes] = {{{-1}}};
    Int_t slatThatTriggered[kMaxNumOfTracks][fgkNchambers][fgkNcathodes]={{{-1}}};
    Int_t boardThatTriggered[kMaxNumOfTracks][fgkNchambers][fgkNcathodes][4]={{{{-1}}}};
    Int_t nboard[4]={-1};
    Int_t ineffBoard[4]={-1};

    const Int_t kMaxNumOfDigits = 20;
    Int_t detElOfDigitsInData[kMaxNumOfDigits][fgkNchambers][fgkNcathodes] = {{{-1}}};

    char filename[150];
    FileStat_t fs;

    if(fFirstRun<0)fFirstRun=fLastRun=-1;

    for(Int_t iRun = fFirstRun; iRun <= fLastRun; iRun++){// Loop over runs
    // open run loader and load gAlice
    if(fFirstRun>=0){
	cout<<"\n\nRun = "<<iRun<<endl;
	sprintf(filename, "%s/run%i/galice.root", fGaliceDir.Data(), iRun);
	if(gSystem->GetPathInfo(filename,fs)){
	    cout<<"Warning: "<<filename<<" not found. Skip to next one"<<endl;
	    continue;
	}
	cout<<"Opening file "<<filename<<endl;
	SetGaliceFile(filename);
    }

    for (Int_t ievent=fFirstEvent; ievent<=fLastEvent; ievent++) { // event loop
	Bool_t isClearEvent = kTRUE;

	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		partNumOfTrig[ch][cath]=0;
		match[ch][cath]=kFALSE;
		for(Int_t itrack=0; itrack<kMaxNumOfTracks; itrack++){
		    triggeredDigits[itrack][ch][cath]=-1;
		    slatThatTriggered[itrack][ch][cath]=-1;
		    for(Int_t loc=0; loc<4; loc++){
			boardThatTriggered[itrack][ch][cath][loc]=-1;
		    }
		}
		for(Int_t idig=0; idig<kMaxNumOfDigits; idig++){
		    detElOfDigitsInData[idig][ch][cath]=-1;
		}
	    }
	}

	fRunLoader->GetEvent(ievent);
	if (ievent%evtBeforePrint==0) printf("\t Event = %d\n",ievent);

	fData->SetTreeAddress("RL");
	fData->GetRecTriggerTracks();
	recTrigTracksArray = fData->RecTriggerTracks();
	Int_t nRecTrigTracks = (Int_t) recTrigTracksArray->GetEntriesFast();

	fData->SetTreeAddress("D,GLT");
	fData->GetDigits();

	const AliMUONGeometryTransformer* kGeomTransformer = fMUON->GetGeometryTransformer();

	for (Int_t iRecTrigTrack=0; iRecTrigTrack<nRecTrigTracks; iRecTrigTrack++) {
	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		digitPerTrack[cath]=0;
		for(Int_t sta=0; sta<2; sta++){
		    slatInPlane1[sta][cath] = -9999;
		    iyDigitInPlane1[sta][cath] = -9999;
		}
	    }

	    Bool_t doubleCountTrack = kFALSE;

	    // reading info from tracks
	    recTrigTrack = (AliMUONTriggerTrack *)recTrigTracksArray->At(iRecTrigTrack);
	    Float_t x11 = recTrigTrack->GetX11();// x position (info from non-bending plane)
	    Float_t y11 = recTrigTrack->GetY11();// y position (info from bending plane)
	    Float_t thetaX = recTrigTrack->GetThetax();
	    Float_t thetaY = recTrigTrack->GetThetay();

	    if(fDebugLevel>=3)printf("\tEvent = %i, Track = %i\npos from track: (x,y) = (%f, %f), (thetaX, thetaY) = (%f, %f)\n",ievent,iRecTrigTrack,x11,y11,thetaX*rad2deg,thetaY*rad2deg);

	    for(Int_t ch=0; ch<fgkNchambers; ch++) {
		zRealMatch[ch] = zMeanChamber[ch];
		for(Int_t cath=0; cath<fgkNcathodes; cath++){
		    trigScheme[iRecTrigTrack][ch][cath] = 0;
		}
	    }

	    for(Int_t ch=0; ch<fgkNchambers; ch++) { // chamber loop
		Int_t currCh = chOrder[ch];
		Int_t ichamber = 10+currCh;
		Int_t currStation = station[currCh];
		TClonesArray* digits = fData->Digits(ichamber);
		digits->Sort();
		Int_t ndigits = (Int_t)digits->GetEntriesFast();
		if(fDebugLevel>=2){
		    if(fDebugLevel<3)printf("\tEvent = %i, Track = %i\n", ievent, iRecTrigTrack);
		    printf("DigitNum: %i digits detected\n",ndigits);
		}

		for(Int_t cath=0; cath<fgkNcathodes; cath++){
		    correctFactor[cath]=1.;
		}
		// calculate corrections to trigger track theta
		if(ch>=1)correctFactor[kNonBending] = zMeanChamber[0]/zRealMatch[0];// corrects x position
		if(ch>=2)correctFactor[kBending] = (zMeanChamber[2] - zMeanChamber[0]) / (zRealMatch[2] - zRealMatch[0]);// corrects y position

		// searching track intersection with chambers (first approximation)
		Float_t deltaZ = zMeanChamber[currCh] - zMeanChamber[0];
		trackIntersectCh[0][currCh] = zMeanChamber[currCh] * TMath::Tan(thetaX) * correctFactor[kNonBending];// x position (info from non-bending plane) 
		trackIntersectCh[1][currCh] = y11 + deltaZ * TMath::Tan(thetaY) * correctFactor[kBending];// y position (info from bending plane)

		for(Int_t idigit=0; idigit<ndigits; idigit++) { // digit loop
		    mDigit = (AliMUONDigit*)digits->At(idigit);
		    for(Int_t loc=0; loc<4; loc++){
			nboard[loc]=-1;
		    }

		    // searching loaded digit global position and dimension
		    Int_t detElemId = mDigit->DetElemId();
		    Int_t cathode = mDigit->Cathode();
		    Int_t ix = mDigit->PadX();
		    Int_t iy = mDigit->PadY();
		    Float_t xpad, ypad, zpad;
		    if(detElOfDigitsInData[idigit][ch][cathode]==-1)detElOfDigitsInData[idigit][ch][cathode] = detElemId;

		    if(fDebugLevel>=2)printf("cathode = %i\n",cathode);
		    const AliMpVSegmentation* seg = 
			AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,cathode);

		    AliMpPad pad = seg->PadByIndices(AliMpIntPair(ix,iy),kTRUE);
		    for (Int_t loc=0; loc<pad.GetNofLocations(); loc++){
			AliMpIntPair location = pad.GetLocation(loc);
			nboard[loc] = location.GetFirst();
		    }

		    // get the pad position and dimensions
		    Float_t xlocal1 = pad.Position().X();
		    Float_t ylocal1 = pad.Position().Y();
		    Float_t dpx = pad.Dimensions().X();
		    Float_t dpy = pad.Dimensions().Y();

		    kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xpad, ypad, zpad);

		    if(fDebugLevel>=3)printf("ch = %i\t cath = %i\tpad = (%4.1f, %4.1f, %4.1f)\tsize = (%3.1f, %3.1f)\n",currCh,cathode,xpad,ypad,zpad,dpx,dpy);

		    // searching track intersection with chambers (second approximation)
		    if(ch>=2){
			deltaZ = zpad - zRealMatch[0];
			trackIntersectCh[1][currCh] = y11 + deltaZ * TMath::Tan(thetaY) * correctFactor[kBending];// y position (info from bending plane)
		    }

		    // deciding if digit matches track
		    Bool_t isDiffLocBoard = kFALSE;
		    if(fReproduceTrigResponse)isDiffLocBoard = IsDiffLocalBoard(detElemId, iy, slatInPlane1[currStation][cathode], iyDigitInPlane1[currStation][cathode]);
		    Bool_t matchPad = PadMatchTrack(xpad, ypad, dpx, dpy, trackIntersectCh[0][currCh], trackIntersectCh[1][currCh], currCh);

		    if(matchPad && ch<2){
			slatInPlane1[currStation][cathode] = detElemId;
			iyDigitInPlane1[currStation][cathode] = iy;
			if(fDebugLevel>=3)printf("slatInPlane1[%i][%i] = %i\tiyDigitInPlane1[%i][%i] = %i\n",currStation,cathode,slatInPlane1[currStation][cathode],currStation,cathode,iyDigitInPlane1[currStation][cathode]);
		    }

		    if(isDiffLocBoard && fDebugLevel>=1)printf("\tDifferent local board\n");

		    match[currCh][cathode] = (matchPad && !isDiffLocBoard);

		    if(match[currCh][cathode]){
			digitPerTrack[cathode]++;
			trigScheme[iRecTrigTrack][currCh][cathode]++;
			triggeredDigits[iRecTrigTrack][currCh][cathode] = idigit;
			slatThatTriggered[iRecTrigTrack][currCh][cathode] = detElemId;
			for(Int_t loc=0; loc<4; loc++){
			    boardThatTriggered[iRecTrigTrack][currCh][cathode][loc] = nboard[loc];
			}
			if(digitPerTrack[cathode]>4 && !fReproduceTrigResponse)isClearEvent = kFALSE;
		    }

		    // in case of match, store real z position of the chamber
		    if(matchPad)zRealMatch[currCh] = zpad;

		} // end digit loop
	    } // end chamber loop

	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		if(digitPerTrack[cath]<3 && !fReproduceTrigResponse)isClearEvent = kFALSE;
	    }

	    if(!isClearEvent && !fReproduceTrigResponse){ 
		fData->ResetDigits();
		continue;
	    }

	    Int_t commonDigits = 0;
	    Int_t doubleTrack = -1;
	    for(Int_t itrack=0; itrack<iRecTrigTrack; itrack++){
		for(Int_t ch=0; ch<fgkNchambers; ch++){
		    if(triggeredDigits[itrack][ch][kBending]==triggeredDigits[iRecTrigTrack][ch][kBending])commonDigits++;
		}
		if(commonDigits>=2){
		    doubleCountTrack=kTRUE;
		    doubleTrack = itrack;
		    break;
		}
	    }

	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		if(!doubleCountTrack || fReproduceTrigResponse){
		    Int_t is44 = 1;
		    Bool_t goodForSlatEff=kTRUE;
		    Bool_t goodForBoardEff=kTRUE;
		    Int_t firstSlat=slatThatTriggered[iRecTrigTrack][0][cath]%100;
		    if(firstSlat<0)firstSlat=slatThatTriggered[iRecTrigTrack][1][cath]%100;
		    Int_t firstBoard=boardThatTriggered[iRecTrigTrack][0][0][0];
		    if(firstBoard<0)firstBoard=boardThatTriggered[iRecTrigTrack][1][0][0];
		    for(Int_t ch=0; ch<fgkNchambers; ch++){
			is44 *= trigScheme[iRecTrigTrack][ch][cath];
			Int_t currSlat=slatThatTriggered[iRecTrigTrack][ch][cath]%100;
			if(currSlat<0)continue;
			if(currSlat!=firstSlat)goodForSlatEff=kFALSE;
			Bool_t atLeastOneLoc=kFALSE;
			for(Int_t loc=0; loc<4; loc++){
			    Int_t currBoard = boardThatTriggered[iRecTrigTrack][ch][cath][loc];
			    if(currBoard==firstBoard){
				atLeastOneLoc=kTRUE;
				break;
			    }
			}
			if(!atLeastOneLoc)goodForBoardEff=kFALSE;
		    }
		    if(fDebugLevel==1)printf("\tEvent = %i, Track = %i\n", ievent, iRecTrigTrack);
		    if(is44==1){
			fTrigger44[cath]++;
			if(fDebugLevel>=1)printf("Trigger44[%i] = %i\n",cath,fTrigger44[cath]);
			if(goodForSlatEff){
			    for(Int_t ch=0; ch<fgkNchambers; ch++){
				for(Int_t slat=0; slat<fgkNslats; slat++){
				    Int_t corrDetEl = (ch+11)*100 + slat;
				    if(corrDetEl==slatThatTriggered[iRecTrigTrack][ch][cath]){
					fHitPerSlat[ch][cath][slat]++;
					if(fDebugLevel>=1)printf("Slat that triggered = %i\n",corrDetEl);
					if(goodForBoardEff && firstBoard>0){
					    fHitPerBoard[ch][cath][firstBoard-1]++;
					    if(fDebugLevel>=1)printf("Board that triggered = %i\n",firstBoard);
					}
					else if(fDebugLevel>=1)printf("Event = %i, Track = %i: Particle crossed different boards: rejected!\n",ievent,iRecTrigTrack);
				    }
				}
			    }
			}
			else printf("Event = %i, Track = %i: Particle crossed different slats: rejected!\n",ievent,iRecTrigTrack);
		    }
		    if(digitPerTrack[cath]==3){
			for(Int_t ch=0; ch<fgkNchambers; ch++){
			    if(match[ch][cath])partNumOfTrig[ch][cath]++;
			    if(trigScheme[iRecTrigTrack][ch][cath]==0){
				fTrigger34[ch][cath]++;
				if(fDebugLevel>=1)printf("Trigger34[%i][%i] = %i\n",ch,cath,fTrigger34[ch][cath]);
				if(!goodForSlatEff){
				    printf("Event %i, Track = %i: Particle crossed different slats: rejected!\n",ievent,iRecTrigTrack);
				    continue;
				}
				Int_t ineffSlat = DetElemIdFromPos(trackIntersectCh[0][ch], trackIntersectCh[1][ch], 11+ch, cath);
				if(fDebugLevel>=1)printf("Slat non efficient = %i\n",ineffSlat);
				if(ineffSlat>0){
				    Int_t slatInCh = ineffSlat%100;
				    fInefficientSlat[ch][cath][slatInCh]++;
				    for(Int_t idig=0; idig<kMaxNumOfDigits; idig++){
					if(ineffSlat==detElOfDigitsInData[idig][ch][cath])cout<<"Warning: "<<ineffSlat<<" is not inefficient!!!"<<endl;
				    }
				    LocalBoardFromPos(trackIntersectCh[0][ch], trackIntersectCh[1][ch], ineffSlat, cath, ineffBoard);
				    Int_t boardNonEff=-1;
				    for(Int_t loc=0; loc<4; loc++){
					if(ineffBoard[loc]==firstBoard){
					    boardNonEff=ineffBoard[loc];
					    break;
					}
				    }
				    if(fDebugLevel>=1)printf("Board non efficient = %i\n",boardNonEff);
				    if(boardNonEff>0)fInefficientBoard[ch][cath][boardNonEff-1]++;
				    else if(fDebugLevel>=1){
					printf("Inefficient board should be %i.\tBoards found:\n", firstBoard);
					for(Int_t loc=0; loc<4; loc++){
					    printf("%i\t",ineffBoard[loc]);
					}
					printf("\n");
				    }
				    
				}
			    }
			}
		    }
		}
		else if(doubleCountTrack){
		    if(fDebugLevel<=1)printf("\n\tEvent = %i, Track = %i: ", ievent,iRecTrigTrack);
		    printf("Double Count Track: %i similar to %i. Track rejected!\n",iRecTrigTrack, doubleTrack);
		}
	    }
	}// end trigger tracks loop
	if(nRecTrigTracks<=0){ 
	    fData->ResetDigits();
	    continue;
	}
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    for(Int_t cath=0; cath<fgkNcathodes; cath++){
		totNumOfTrig[ch][cath] += partNumOfTrig[ch][cath];
		if(partNumOfTrig[ch][cath]>0)atLeast1MuPerEv[ch][cath]++;
	    }
	}

	if(fPrintInfo){
	    //Global trigger
	    globalTrigger = fData->GlobalTrigger();
	    Int_t nglobals = (Int_t) globalTrigger->GetEntriesFast(); // should be 1

	    for (Int_t iglobal=0; iglobal<nglobals; iglobal++) { // Global Trigger
		gloTrg = (AliMUONGlobalTrigger*)globalTrigger->At(iglobal);
	    }
	    PrintTrigger(gloTrg);
	    InfoDigit();
	    cout<<"\n"<<endl;
	}

	fData->ResetDigits();
    }// end event loop
    if(fFirstRun>=0)CleanGalice();
    } //end loop over run
    
    // Write output data
    WriteEfficiencyMap(outputDir);

    WriteOutput(outputDir, totNumOfTrig, atLeast1MuPerEv);
}

//_____________________________________________________________________________
void AliMUONTriggerChamberEff::WriteOutput(const char* outputDir, Int_t totNumOfTrig[4][2], Int_t atLeast1MuPerEv[4][2])
{
    //
    /// Writes information on calculated efficiency.
    /// It writes: triggerChamberEff.root file containing efficiency histograms.
    ///
    /// In addition a text file triggerChamberEff.out is created,
    /// with further informations on efficiencies.
    //

    char *cathodeName[fgkNcathodes]={"Bending plane", "Non-Bending plane"};
    char *cathCode[fgkNcathodes] = {"bendPlane", "nonBendPlane"};

    char outFileName[100];
    sprintf(outFileName, "%s/triggerChamberEff.out",outputDir);
    FILE *outfile = fopen(outFileName, "w");
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	fprintf(outfile,"%s:\n",cathodeName[cath]);
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    fprintf(outfile,"Total number of muon triggers chamber 1%i = %i\n",ch+1,totNumOfTrig[ch][cath]);
	}
	fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	fprintf(outfile,"%s:\n",cathodeName[cath]);
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    fprintf(outfile,"At least 1 muon triggered chamber 1%i = %i\n",ch+1, atLeast1MuPerEv[ch][cath]);
	}
	fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n\n");
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	fprintf(outfile,"%s:\n",cathodeName[cath]);
	fprintf(outfile,"Number of triggers where all chambers counted = %i\n",fTrigger44[cath]);
	fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	fprintf(outfile,"%s:\n",cathodeName[cath]);
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    fprintf(outfile,"Number of triggers where chamber 1%i did not count = %i\n",ch+1,fTrigger34[ch][cath]);
	}
	fprintf(outfile,"\n");
    }

    fprintf(outfile,"\n");
    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	fprintf(outfile,"%s:\n",cathodeName[cath]);
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    Int_t sumIneff = 0, sumHits = 0;
	    fprintf(outfile,"\n Chamber %1i\n", ch+1);
	    for(Int_t slat=0; slat<fgkNslats; slat++){
		fprintf(outfile,"Number of triggers where slat %2i - did not count = %5i - was hit (hit%sCh%iSlat%i) = %5i\n",slat,fInefficientSlat[ch][cath][slat],cathCode[cath],11+ch,slat,fHitPerSlat[ch][cath][slat]);
		sumIneff += fInefficientSlat[ch][cath][slat];
		sumHits += fHitPerSlat[ch][cath][slat];
	    }
	    fprintf(outfile,"Number of triggers where chamber %1i - did not count = %5i - was hit (hit%sCh%i) = %5i\n",ch+1,sumIneff,cathCode[cath],11+ch,sumHits);
	}
	fprintf(outfile,"\n");
    }
    fclose(outfile);

    sprintf(outFileName, "%s/triggerChamberEff.root",outputDir);
    TFile *outputHistoFile = new TFile(outFileName,"RECREATE");
    TDirectory *dir = gDirectory;

    enum {kSlatIn11, kSlatIn12, kSlatIn13, kSlatIn14, kChamberEff};
    char *yAxisTitle = "trigger efficiency (a.u.)";
    char *xAxisTitle = "chamber";

    TH1F *histo[fgkNcathodes][fgkNchambers+1];
    TH1F *histoBoard[fgkNcathodes][fgkNchambers];

    char histoName[30];
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
	    }
	}
    }

    Float_t efficiency, efficiencyError;

    for(Int_t cath=0; cath<fgkNcathodes; cath++){
	for(Int_t ch=0; ch<fgkNchambers; ch++){
	    for(Int_t slat=0; slat<fgkNslats; slat++){
		CalculateEfficiency(fHitPerSlat[ch][cath][slat], fHitPerSlat[ch][cath][slat]+fInefficientSlat[ch][cath][slat], efficiency, efficiencyError, kFALSE);
		histo[cath][ch]->SetBinContent(slat+1, efficiency);
		histo[cath][ch]->SetBinError(slat+1, efficiencyError);
	    }
	    CalculateEfficiency(fTrigger44[cath], fTrigger34[ch][cath]+fTrigger44[cath], efficiency, efficiencyError, kFALSE);
	    histo[cath][kChamberEff]->SetBinContent(ch+1, efficiency);
	    histo[cath][kChamberEff]->SetBinError(ch+1, efficiencyError);

	    for(Int_t board=0; board<fgkNboards; board++){
		CalculateEfficiency(fHitPerBoard[ch][cath][board], fHitPerBoard[ch][cath][board]+fInefficientBoard[ch][cath][board], efficiency, efficiencyError, kFALSE);
		histoBoard[cath][ch]->SetBinContent(board+1, efficiency);
		histoBoard[cath][ch]->SetBinError(board+1, efficiencyError);
	    }
	}
    }

// write all histos
    outputHistoFile->cd();
    dir->GetList()->Write();
    outputHistoFile->Close();
}


//_____________________________________________________________________________
void AliMUONTriggerChamberEff::WriteEfficiencyMap(const char* outputDir)
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

