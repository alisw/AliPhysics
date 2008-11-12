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

//Class to check the response of the detection elements of the  MUON tracking chambers 
//in function of the position in the detection element.
//Author:  Nicolas LE BRIS - SUBATECH Nantes

//PWG3/muon:
#include "AliAnalysisTaskMuonTrackingEff.h"
#include "AliCheckMuonDetEltResponse.h"

//include STEER:
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliTracker.h"
#include "AliESDMuonTrack.h"

//include MUON:
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONESDInterface.h"

//include MUON/mapping:
#include "mapping/AliMpDEManager.h"

//include ROOT:
#include <Riostream.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TClonesArray.h>
#include <TPaveLabel.h>
#include <TList.h>

/// \cond CLASSIMP
ClassImp(AliCheckMuonDetEltResponse)
/// \endcond

const Int_t AliCheckMuonDetEltResponse::fNbrOfChamber          = 10;
const Int_t AliCheckMuonDetEltResponse::fNbrOfStation          = 5;
const Int_t AliCheckMuonDetEltResponse::fNbrOfDetectionElt[10] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
const Int_t AliCheckMuonDetEltResponse::fFirstDetectionElt[10] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
const Int_t AliCheckMuonDetEltResponse::fOffset                = 100;
const Int_t AliCheckMuonDetEltResponse::fOverlapSize           = 15;
const Int_t AliCheckMuonDetEltResponse::fYSlatSize             = 40;

//_____________________________________________________________________________
AliCheckMuonDetEltResponse::AliCheckMuonDetEltResponse() 
: TObject(),
  fNCh(0),
  fNSt(0),
  fNDE(0),
  fTransformer(0x0),
  fESD(0x0),
  fTracksTotalNbr(0x0),
  fTrackParams(0x0),
  fTrackParam(0x0),
  fCluster(0x0),
  fDetEltTDHistList(0x0),
  fDetEltTTHistList(0x0)
{
/// Default constructor

    fNCh = AliCheckMuonDetEltResponse::fNbrOfChamber;
    fNSt = AliCheckMuonDetEltResponse::fNbrOfStation;
    fNDE = AliAnalysisTaskMuonTrackingEff::fTotNbrOfDetectionElt;

    for (Int_t iCluster = 0; iCluster<fNCh; ++iCluster)
    {
      fNbrClustersCh[iCluster] = 0;
    }  

    for (Int_t i=0; i<2; ++i)
    {
      fGetDetElt[i] = 0;
    }
    for (Int_t i=0; i<fNCh; ++i)
    {
      fTrackFilter[i] = 0;
    } 
}

//_____________________________________________________________________________
AliCheckMuonDetEltResponse::AliCheckMuonDetEltResponse(const AliCheckMuonDetEltResponse& src) 
: TObject(src),
  fNCh(0),
  fNSt(0),
  fNDE(0),
  fTransformer(0x0),
  fESD(0x0),
  fTracksTotalNbr(0x0),
  fTrackParams(0x0),
  fTrackParam(0x0),
  fCluster(0x0),
  fDetEltTDHistList(0x0),
  fDetEltTTHistList(0x0)
{
 src.Copy(*this);
}
//_____________________________________________________________________________
AliCheckMuonDetEltResponse& AliCheckMuonDetEltResponse::operator=(const AliCheckMuonDetEltResponse& src) 
{
  /// assignement operator
  if ( this != &src ) 
  {
    src.Copy(*this);
  }
  return *this;
}

//_____________________________________________________________________________
AliCheckMuonDetEltResponse::AliCheckMuonDetEltResponse(const AliMUONGeometryTransformer* transformer,
						       AliESDEvent* esd,
						       TClonesArray* detEltTDHistList,
						       TClonesArray* detEltTTHistList) 
: TObject(),
  fNCh(0),
  fNSt(0),
  fNDE(0),
  fTransformer(transformer),
  fESD(esd),
  fTracksTotalNbr(0),
  fTrackParams(0x0),
  fTrackParam(0),
  fCluster(0), 
  fDetEltTDHistList(detEltTDHistList),
  fDetEltTTHistList(detEltTTHistList)
{
/// Constructor

    fNCh = AliCheckMuonDetEltResponse::fNbrOfChamber;
    fNSt = AliCheckMuonDetEltResponse::fNbrOfStation;
    fNDE = AliAnalysisTaskMuonTrackingEff::fTotNbrOfDetectionElt;

    for (Int_t iCluster = 0; iCluster<fNCh; ++iCluster)
    {
      fNbrClustersCh[iCluster] = 0;
    }

    for (Int_t i = 0; i <2 ; ++i)
    {
      fGetDetElt[i] = 0;
    }
    for (Int_t i=0; i<fNCh; ++i)
    {
      fTrackFilter[i] = 0;
    } 
}


//_____________________________________________________________________________
AliCheckMuonDetEltResponse::~AliCheckMuonDetEltResponse()

{
/// Destructor
    delete fTrackParams;
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::CheckDetEltResponse()
{
//
//Cataloging positions (X,Y) of the clusters detected in the detection elements
//(fDetEltTDHistList), and positions of crossing points between all the
//tracks and the detection elements (fDetEltTTHistList).
//Efficiency = 100 * fDetEltTDHistList / fDetEltTTHistList.

//Loop on tracks
//--------------
    TrackLoop();
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::TrackLoop()
{
    AliESDMuonTrack* esdTrack;
    AliMUONTrack track;
    Int_t nTracks, iTrack;

    nTracks = (Int_t)fESD -> GetNumberOfMuonTracks();
    fTrackParams = new TClonesArray();
 ///Begininig of the loop:
    for (iTrack = 0; iTrack < nTracks; iTrack++)
    {
      esdTrack   = fESD -> GetMuonTrack(iTrack);
  
      if( esdTrack->ContainTrackerData() )
      {
	AliMUONESDInterface::ESDToMUON(*esdTrack, track);
        fTrackParams = track.GetTrackParamAtCluster();
	TrackParamLoop(); //!<Loop on trackParam.
      }
    }
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::TrackParamLoop()
{
    Int_t nTrackParams = (Int_t) fTrackParams->GetEntriesFast();  //!<Number of trackParams in the track.
    Int_t iTrackParam = 0;                                        //!<Number of the trackParam of the track.
    Int_t oldChamber = -1, newChamber = 0; //!<To check if there is 0, 1 or 2 (overlap cases) clusters in the same chamber for a track.                                      
    Int_t detElt;                          //!<Detection element Id.
  
    for (Int_t ch = 0; ch < fNCh; ++ch)
    {
      fTrackFilter[ch] = 0;
    }

    Double_t posXL, posYL, posZL;          //!<Local positions.
    Double_t posXG, posYG, posZG;          //!<Global. positions.
    Int_t chamberResponse [10] = {0};      //!<1 if the chamber has responded; 0 if not

    for (iTrackParam = 0; iTrackParam < nTrackParams; ++iTrackParam)
    { 
      fTrackParam = (AliMUONTrackParam*) fTrackParams->At(iTrackParam);
      fCluster    = (AliMUONVCluster*  ) fTrackParam ->GetClusterPtr();    
      fTrackFilter   [fCluster->GetChamberId()] = 1;
      chamberResponse[fCluster->GetChamberId()] = 1;
    }

    for (Int_t station = 0; station < fNSt-1; ++station)
    {
      Int_t filter;                                                  //<!
      Int_t ch1, ch2, ch3, ch4;                                      //<!
      ch1 = 2*station;                                               //<!
      ch2 = 2*station + 1;                                           //<!
      ch3 = 2*station + 2;                                           //<!
      ch4 = 2*station + 3;                                           //<!
                                                                     //<!For the efficiency calculation the tracks
      if (station < 3 )                                              //<!reconstructed must have responded to the
      {                                                              //<!criteria of the tracking. 
	filter            = fTrackFilter[ch1];                       //<!And that's why the tracks usable for the 
	fTrackFilter[ch1] = fTrackFilter[ch2];                       //<!intrinsic efficiency calculation are
	fTrackFilter[ch2] = filter;                                  //<!the tracks which have one or two clusters
      }                                                              //<!in each station. So the case where a track
                                                                     //<!hasn't a cluster in a station is not
      else                                                           //<!taking into account.
      {                                                              //<!This part solves the problem. See the ALICE 
	if (chamberResponse[ch3]*chamberResponse[ch4] != 0)          //<!note of Diego STOCCO on the trigger efficiency
	{                                                            //<!
	filter            = fTrackFilter[ch1];                       //<!
	fTrackFilter[ch1] = fTrackFilter[ch2];                       //<!
	fTrackFilter[ch2] = filter;                                  //<!
	}                                                            //<!
	else                                                         //<!
	{                                                            //<!
	fTrackFilter[ch1] = 0;                                       //<!
	fTrackFilter[ch2] = 0;                                       //<!
	}                                                            //<!
                                                                     //<!
	if (chamberResponse[ch1]*chamberResponse[ch2] != 0)          //<!
	{                                                            //<!
	filter            = fTrackFilter[ch3];                       //<!
	fTrackFilter[ch3] = fTrackFilter[ch4];                       //<!
	fTrackFilter[ch4] = filter;                                  //<!
	}                                                            //<!
	else                                                         //<!
	{                                                            //<!
	fTrackFilter[ch3] = 0;                                       //<!
	fTrackFilter[ch4] = 0;                                       //<!
	}                                                            //<!
      }                                                              //<!
    }                                                                //<!

    for (Int_t ch = 0; ch < fNCh; ++ch)
    {
      if (fTrackFilter[ch] == 1)
      {
	if ( chamberResponse[ch] != 0) ((TH2F*) fDetEltTDHistList->UncheckedAt(fNDE))->Fill(ch, 0);
	if (!fDetEltTTHistList) ((TH2F*) fDetEltTTHistList->UncheckedAt(fNDE))->Fill(ch, 0);
      }
    } 

 ///Begining of the loop:
    for (iTrackParam = 0; iTrackParam < nTrackParams; ++iTrackParam)
    {
      fTrackParam = (AliMUONTrackParam*) fTrackParams->At(iTrackParam);
      fCluster    = (AliMUONVCluster*  ) fTrackParam ->GetClusterPtr(); 
      
      newChamber  = fCluster->GetChamberId();
      detElt      = fCluster->GetDetElemId();

   ///Global and local positions calculation:
      posXG = fTrackParam->GetNonBendingCoor(); 
      posYG = fTrackParam->GetBendingCoor(); 
      posZG = fTrackParam->GetZ(); 

      fTransformer->Global2Local(detElt, posXG, posYG, posZG, posXL, posYL, posZL);  //!<Transfomation from global to local positions.
   
   ///Filling histograms of the cluster positions on the detection element of the TRACKS DETECTED (TD):
      FillDetEltTDHisto(newChamber, detElt, posXL, posYL);

   ///Filling histograms of the cluster positions on the detection element of ALL THE TRACKS (TT):
      FillDetEltTTHisto(newChamber, detElt, posXG, posYG, posZG, posXL, posYL, posZL);

      if (newChamber != oldChamber) 
      {
	if (newChamber > oldChamber + 1)                                 //!<Check if it doesn't miss a chamber.
	{
	  Int_t nbrOfMissCh = newChamber - (oldChamber+1);                //!<Number of missing chambers.
	  CalculMissClusterParam(fTrackParam, oldChamber+1, nbrOfMissCh); //!<Calculation of the parameters of the missing cluster(s).
	}
	if ( iTrackParam == nTrackParams - 1 && newChamber != fNCh-1)           //!<Check if the last chamber, chamber 9 (from 0 to 9) has responded.
	{
	  CalculMissClusterParam(fTrackParam, fNCh-1, 1);                      //!<Calculation of the parameters of the missing cluster(s) in the last chamber.
	}
      }
      oldChamber = newChamber; 
    }
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::FillDetEltTDHisto(Int_t chamber,
						   Int_t detElt,
						   Double_t posXL,
						   Double_t posYL)
{
    if(fTrackFilter[chamber]== 1)
    {
      Int_t iDet = 0; //!<Position of the detection element in the histograms' list.
      iDet = FromDetElt2iDet(chamber, detElt);
      ((TH2F*) fDetEltTDHistList->UncheckedAt(iDet))->Fill(posXL, posYL);
    }
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::FillDetEltTTHisto(Int_t chamber,
						   Int_t detElt,
						   Double_t posXG,
						   Double_t posYG,
						   Double_t posZG,
						   Double_t posXL,
						   Double_t posYL,
						   Double_t posZL)
{
    if(fTrackFilter[chamber] == 1)
    {
      Int_t iDet = 0; //!<Position of the detection element in the histograms' list.
      iDet = FromDetElt2iDet(chamber, detElt);
      if (!fDetEltTTHistList)  ((TH2F*) fDetEltTTHistList->UncheckedAt(iDet)) -> Fill(posXL, posYL);
      
      if(TMath::Abs(posYL) > fOverlapSize) //!<It is an overlap area. It can have two clusters for this track in this chamber.
      {
	GetDetEltFromPosition(chamber, posXG, posYG, posZG);
	if(fGetDetElt[1] != 0) //<!There is a second detection element for the same (X,Y) in this chamber (overlap).
	{
	  if(fGetDetElt[1] != detElt)
	  {
	    fTransformer->Global2Local(fGetDetElt[1], posXG, posYG, posZG, posXL, posYL, posZL);  //!<Transfomation from global to local positions.
	    iDet = FromDetElt2iDet(chamber, fGetDetElt[1]);
	    if (!fDetEltTTHistList) ((TH2F*) fDetEltTTHistList->UncheckedAt(iDet)) -> Fill(posXL, posYL);
	  }
	  else
	  {
	    fTransformer->Global2Local(fGetDetElt[0], posXG, posYG, posZG, posXL, posYL, posZL);  //!<Transfomation from global to local positions.
	    iDet = FromDetElt2iDet(chamber, fGetDetElt[0]);
	    if (!fDetEltTTHistList) ((TH2F*) fDetEltTTHistList->UncheckedAt(iDet)) -> Fill(posXL, posYL);
	  }
	}
      }
    }
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::CalculMissClusterParam(AliMUONTrackParam* extrapTrackParam,
							Int_t firstMissCh,
							Int_t nbrOfMissCh)
{
//Calculation of the cluster parameter which was not detect by
//the chamber.

    for (Int_t iCh = 0; iCh<nbrOfMissCh; ++iCh)
    {
      Double_t exTraXL, exTraYL, exTraZL;   //!<Extrapolated local positions.
      Double_t exTraXG, exTraYG, exTraZG;   //!<Extrapolated global positions.
      Int_t missChamber= firstMissCh + iCh; //!<The missing chamber.
      Int_t missDetElt = 0;

      exTraZG = AliMUONConstants::DefaultChamberZ(missChamber);
      AliMUONTrackExtrap::ExtrapToZ(extrapTrackParam, exTraZG);      //!<Extrapolation to the missing chamber.
      exTraXG = extrapTrackParam->GetNonBendingCoor();               //!<Global X position extrapolated to the missing chamber.
      exTraYG = extrapTrackParam->GetBendingCoor();                  //!<Global Y position extrapolated to the missing chamber.

      GetDetEltFromPosition(missChamber, exTraXG, exTraYG, exTraZG); //!<Gives the detection element (fGetDetElt) with the position.
      missDetElt = fGetDetElt[0];

      if( missDetElt != 0 ) //!<Check if the track has passed trough a detection area
      {
	fTransformer->Global2Local(missDetElt, exTraXG, exTraYG, exTraZG, exTraXL, exTraYL, exTraZL);  //!<Transfomation from global to local positions.
	FillDetEltTTHisto(missChamber, missDetElt, exTraXG, exTraYG, exTraZG, exTraXL, exTraYL, exTraZL);
      }
    }
}



//_____________________________________________________________________________
Int_t AliCheckMuonDetEltResponse::FromDetElt2iDet(Int_t chamber, 
						  Int_t detElt)
{
///
///Connexion between the detection element X and its position in the list of histograms iX.
///

    Int_t iDet = 0; //!<Position of the detection element (detElt) in the histograms' list.

    if (chamber<4)             iDet = detElt-fOffset*(chamber+1)+ 4* chamber      ; 
    if (chamber>3 && chamber<6) iDet = detElt-fOffset*(chamber+1)+18*(chamber-4)+16;
    if (chamber>5)             iDet = detElt-fOffset*(chamber+1)+26*(chamber-6)+52;

    return iDet;    
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::GetDetEltFromPosition(Int_t chamber,
						       Double_t posX,
						       Double_t posY,
						       Double_t posZ)
{
///
///Gives the detetection element fGetDetElt[0] corresponding to the position. In the case 
///of an overlap (two detection element with the same (X,Y)) fGetDetElt[1] is calculated.
///

    Int_t nbrOfDetElt =  AliMpDEManager::GetNofDEInChamber(chamber, kTRUE); //!<Number of detection elements in the chamber.
    Int_t detElt = 0, lastDetElt = 0;

    Double_t posXL, posYL, posZL; //!<Local positions.
    posXL = posYL = posZL = 1000.;
    fGetDetElt[0] = fGetDetElt[1] = 0;

    if(chamber>3)
    {
      Int_t shift  = 0;                                                    //!<|
      if(posX<0) shift =    nbrOfDetElt /4 + 1;                            //!<|Method to avoid the loop on all elements of
      if(posX>0) shift = (3*nbrOfDetElt)/4 + 1;                            //!<|detection, and just loop through fourth chamber.
                                                                           //!<|
      detElt = fOffset*(chamber+1) + shift;                                //!<|
      lastDetElt = fOffset*(chamber+1) +(shift+nbrOfDetElt/2)%nbrOfDetElt; //!<|
      
      Int_t nbr = 0;
      while(detElt != lastDetElt) //!<Loop on fourth chamber.
      {

	fTransformer->Global2Local(detElt, posX, posY, posZ, posXL, posYL, posZL);  //!<Transfomation from global to local positions.

	if(TMath::Abs(posYL)< fYSlatSize) //!<If |posYL|<20 => the cluster is in the detection element (-20<Ylocal<20 for each slat). 	
	{
	  fGetDetElt[nbr] = detElt;  
	  ++nbr;
	}	
	shift  = (shift + 1)%nbrOfDetElt;
	detElt = fOffset*(chamber+1) + shift;
      }
    }
    
    else //!<For the station 1 & 2 (4 detection elements in each chamber). 
    {
      if(posX>0 && posY>0) fGetDetElt[0] = fOffset*(chamber+1)    ;
      if(posX<0 && posY>0) fGetDetElt[0] = fOffset*(chamber+1) + 1;
      if(posX<0 && posY<0) fGetDetElt[0] = fOffset*(chamber+1) + 2;
      if(posX>0 && posY<0) fGetDetElt[0] = fOffset*(chamber+1) + 3; 
    }
}
