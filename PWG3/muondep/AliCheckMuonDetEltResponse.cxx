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
//Modified by Matthieu LENHARDT - SUBATECH Nantes

//PWG3/muondep:
#include "AliCheckMuonDetEltResponse.h"

//include STEER:
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"

//include MUON:
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMUONConstants.h"
#include "AliMUONESDInterface.h"
#include "AliMUONGeometryTransformer.h"

//include MUON/mapping:
#include "mapping/AliMpDEManager.h"
#include "mapping/AliMpSegmentation.h"
#include "mapping/AliMpSlatSegmentation.h"
#include "mapping/AliMpSectorSegmentation.h"
#include "mapping/AliMpPad.h"

//include ROOT:
#include <TH2F.h>
#include <TList.h>
#include <TClonesArray.h>

/// \cond CLASSIMP
ClassImp(AliCheckMuonDetEltResponse)
/// \endcond

const Int_t AliCheckMuonDetEltResponse::fgkNCh                   = AliMUONConstants::NTrackingCh();
const Int_t AliCheckMuonDetEltResponse::fgkNSt                   = AliMUONConstants::NTrackingSt();
const Int_t AliCheckMuonDetEltResponse::fgkNDE                   = 156;
const Int_t AliCheckMuonDetEltResponse::fgkNbrOfDetectionElt[10] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
const Int_t AliCheckMuonDetEltResponse::fgkFirstDetectionElt[10] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
const Int_t AliCheckMuonDetEltResponse::fgkOffset                = 100;

//_____________________________________________________________________________
AliCheckMuonDetEltResponse::AliCheckMuonDetEltResponse() 
: TObject(),
  fkTransformer(0x0),
  fESD(0x0),
  fTracksTotalNbr(0x0),
  fNbrUsableTracks(0),
  fTrackParams(0x0),
  fTrackParam(0x0),
  fCluster(0x0),
  fDetEltTDHistList(0),
  fDetEltTTHistList(0),
  fChamberTDHistList(0),
  fChamberTTHistList(0)
{
/// Default constructor

    fNbrUsableTracks = 0;

    for (Int_t iCluster = 0; iCluster<fgkNCh; ++iCluster)
      fNbrClustersCh[iCluster] = 0;

    for (Int_t i=0; i<fgkNCh; ++i)
      fTrackFilter[i] = 0;
}

//_____________________________________________________________________________
AliCheckMuonDetEltResponse::AliCheckMuonDetEltResponse(const AliCheckMuonDetEltResponse& src) 
: TObject(src),
  fkTransformer(0x0),
  fESD(0x0),
  fTracksTotalNbr(0x0),
  fNbrUsableTracks(0),
  fTrackParams(0x0),
  fTrackParam(0x0),
  fCluster(0x0),
  fDetEltTDHistList(0),
  fDetEltTTHistList(0),
  fChamberTDHistList(0),
  fChamberTTHistList(0)
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
						       TList* detEltTDHistList,
						       TList* detEltTTHistList,
						       TList* chamberTDHistList,
						       TList* chamberTTHistList) 
: TObject(),
  fkTransformer(transformer),
  fESD(esd),
  fTracksTotalNbr(0),
  fNbrUsableTracks(0),
  fTrackParams(0x0),
  fTrackParam(0),
  fCluster(0),
  fDetEltTDHistList(detEltTDHistList),
  fDetEltTTHistList(detEltTTHistList),
  fChamberTDHistList(chamberTDHistList),
  fChamberTTHistList(chamberTTHistList)
{
/// Constructor

    fNbrUsableTracks = 0;

    for (Int_t iCluster = 0; iCluster<fgkNCh; ++iCluster)
      fNbrClustersCh[iCluster] = 0;
    
    for (Int_t i=0; i<fgkNCh; ++i)
      fTrackFilter[i] = 0;
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
// Check if the track is kept
  AliESDMuonTrack* esdTrack;
  AliMUONTrack track;
  Int_t nTracks, iTrack;

  nTracks = (Int_t)fESD -> GetNumberOfMuonTracks();
  fTrackParams = new TClonesArray();
  ///Begininig of the loop:
  //if (fESD->IsTriggerClassFired("CINT1B-ABCE-NOPF-ALL"))
  {
    for (iTrack = 0; iTrack < nTracks; iTrack++)
      {
	esdTrack   = fESD -> GetMuonTrack(iTrack);
	
	if(esdTrack->ContainTrackerData() && esdTrack->GetMatchTrigger() > 0)
	  {
	    AliMUONESDInterface::ESDToMUON(*esdTrack, track);
	    fTrackParams = track.GetTrackParamAtCluster();
	    TrackParamLoop(); //!<Loop on trackParam.
	    fNbrUsableTracks += 1;
	  }
      }
  }
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::TrackParamLoop()
{
// Loop on all the track params and fill the histos
  Int_t nTrackParams = (Int_t) fTrackParams->GetEntriesFast();  //!<Number of trackParams in the track.
  Int_t iTrackParam = 0;                                        //!<Number of the trackParam of the track.
  Int_t oldChamber = -1, newChamber = 0; //!<To check if there is 0, 1 or 2 (overlap cases) clusters in the same chamber for a track.                                      
  Int_t detElt;                          //!<Detection element Id.
  
  for (Int_t ch = 0; ch < fgkNCh; ++ch)
    fTrackFilter[ch] = 0;

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

  for (Int_t station = 0; station < fgkNSt-1; ++station)
    {
      Int_t filter;                                                       //<!
      Int_t ch1, ch2, ch3, ch4;                                           //<!
      ch1 = 2*station;                                                    //<!
      ch2 = 2*station + 1;                                                //<!
      ch3 = 2*station + 2;                                                //<!
      ch4 = 2*station + 3;                                                //<!
                                                                          //<!For the efficiency calculation the tracks
      if (station < 3 )                                                   //<!reconstructed must have responded to the
	{                                                                 //<!criteria of the tracking. 
	  filter            = fTrackFilter[ch1];                          //<!And that's why the tracks usable for the 
	  fTrackFilter[ch1] = fTrackFilter[ch2];                          //<!intrinsic efficiency calculation are
	  fTrackFilter[ch2] = filter;                                     //<!the tracks which have one or two clusters
	}                                                                 //<!in each station. So the case where a track
                                                                          //<!hasn't a cluster in a station is not
      else                                                                //<!taking into account.
	{                                                                 //<!This part solves the problem. See the ALICE 
	  if (chamberResponse[ch3]*chamberResponse[ch4] != 0)             //<!note of Diego STOCCO on the trigger efficiency
	    {                                                             //<!
	      filter            = fTrackFilter[ch1];                      //<!
	      fTrackFilter[ch1] = fTrackFilter[ch2];                      //<!
	      fTrackFilter[ch2] = filter;                                 //<!
	    }                                                             //<!
	  else                                                            //<!
	    {                                                             //<!
	      fTrackFilter[ch1] = 0;                                      //<!
	      fTrackFilter[ch2] = 0;                                      //<!
	    }                                                             //<!

	  if (chamberResponse[ch1]*chamberResponse[ch2] != 0)
	    {
	      filter            = fTrackFilter[ch3];
	      fTrackFilter[ch3] = fTrackFilter[ch4];
	      fTrackFilter[ch4] = filter;
	    }
	  else
	    {
	      fTrackFilter[ch3] = 0;
	      fTrackFilter[ch4] = 0;
	    }
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
      
      fkTransformer->Global2Local(detElt, posXG, posYG, posZG, posXL, posYL, posZL);  //!<Transfomation from global to local positions.
      
      ///Filling histograms of the cluster positions on the detection element of the TRACKS DETECTED (TD):
      FillTDHistos(newChamber, detElt, posXL, posYL);
    
      ///Filling histograms of the cluster positions on the detection element of ALL THE TRACKS (TT):
      FillTTHistos(newChamber, detElt, posXL, posYL);

      if (newChamber != oldChamber) 
	{
	  if (newChamber > oldChamber + 1)                                 //!<Check if it doesn't miss a chamber.
	    {
	      Int_t nbrMissChamber = newChamber - (oldChamber + 1);
	      FindAndFillMissedDetElt(fTrackParam, oldChamber+1, nbrMissChamber); //!<Calculation of the parameters of the missing cluster(s).
	    }
	    
	  if ( iTrackParam == nTrackParams - 1 && newChamber != fgkNCh-1)           //!<Check if the last chamber, chamber 9 (from 0 to 9) has responded.
	    FindAndFillMissedDetElt(fTrackParam, fgkNCh-1, 1);                      //!<Calculation of the parameters of the missing cluster(s) in the last chamber.
	    
	}
      oldChamber = newChamber; 
    } 
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::FillTDHistos(Int_t chamber,
					      Int_t detElt,
					      Double_t posXL,
					      Double_t posYL)
{
// Fill the histo for tracks detected
  if(fTrackFilter[chamber]== 1)
    {
      Int_t iDet = 0; //!<Position of the detection element in the histograms' list.
      iDet = FromDetElt2iDet(chamber, detElt);
      ((TH2F*) fDetEltTDHistList->At(iDet))->Fill(posXL, posYL);

      Int_t detEltLocalId = 0;  //!<Id of the detection element in the station
      detEltLocalId =  FromDetElt2LocalId(chamber, detElt);
      ((TH1F*) fChamberTDHistList->At(chamber))->Fill(detEltLocalId);
      ((TH1F*) fChamberTDHistList->At(fgkNCh))->Fill(chamber + 1);
   }
}




//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::FillTTHistos(Int_t chamber,
					      Int_t detElt,
					      Double_t posXL,
					      Double_t posYL)
{
// Fill the histo for total number of tracks
  if(fTrackFilter[chamber] == 1)
    {
      Int_t iDet = 0; //!<Position of the detection element in the histograms' list.
      iDet = FromDetElt2iDet(chamber, detElt);
      ((TH2F*) fDetEltTTHistList->At(iDet)) -> Fill(posXL, posYL);
     
      Int_t detEltLocalId = 0;  //!<Id of the detection element in the station
      detEltLocalId =  FromDetElt2LocalId(chamber, detElt);
      ((TH1F*) fChamberTTHistList->At(chamber))->Fill(detEltLocalId);
      ((TH1F*) fChamberTTHistList->At(fgkNCh))->Fill(chamber + 1);
    }
}





//_____________________________________________________________________________
Int_t AliCheckMuonDetEltResponse::FromDetElt2iDet(Int_t chamber, 
						  Int_t detElt) const
{
///
///Connexion between the detection element X and its position in the list of histograms iX.
///

    Int_t iDet = 0; //!<Position of the detection element (detElt) in the histograms' list.

    if (chamber<4)             iDet = detElt-fgkOffset*(chamber+1)+ 4* chamber      ; 
    if (chamber>3 && chamber<6) iDet = detElt-fgkOffset*(chamber+1)+18*(chamber-4)+16;
    if (chamber>5)             iDet = detElt-fgkOffset*(chamber+1)+26*(chamber-6)+52;

    return iDet;    
}



//_____________________________________________________________________________
Int_t AliCheckMuonDetEltResponse::FromDetElt2LocalId(Int_t chamber, 
						     Int_t detElt) const
{
///
///Connexion between the detection element X and its number in the station.
///

    Int_t localId = 0; //!<Position of the detection element (detElt) in the histograms' list.
    localId = detElt - (chamber+1) * 100;

    return localId;    
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::FindAndFillMissedDetElt(AliMUONTrackParam* extrapTrackParam, 
							 Int_t firstMissCh,
							 Int_t nbrMissCh)
{
///
///Find which detection elements should have been hit but were missed, 
///and fill the TT histos appropriately
///
  for (Int_t iCh = 0; iCh < nbrMissCh; ++iCh)
    {
      Int_t chamber = firstMissCh + iCh;
      Int_t nbrOfDetElt =  AliMpDEManager::GetNofDEInChamber(chamber, kTRUE); //!<Number of detection elements in the chamber.
      
      Double_t pos1[6] = {0, 0, 0, 0, 0, 0};        //!<First point used to compute the extrapolated point (first 3 for global coordinates, last 3 for local).
      Double_t pos2[6] = {0, 0, 0, 0, 0, 0};        //!<Second point used to compute the extrapolated point (first 3 for global coordinates, last 3 for local).
      Double_t posMiss[2] = {0, 0};                 //!<(X, Y) local coordinates of the missing cluster.
            
      pos1[2] = AliMUONConstants::DefaultChamberZ(chamber);           //!<Z of point 1, defined by being the Z of the chamber in "perfect" position.
      AliMUONTrackExtrap::ExtrapToZ(extrapTrackParam, pos1[2]);
      pos1[0] = extrapTrackParam->GetNonBendingCoor();                //!<X of point 1, extrapolated by following the Track.
      pos1[1] = extrapTrackParam->GetBendingCoor();                   //!<Y of point 1, extrapolated by following the Track.
      
      pos2[2] = AliMUONConstants::DefaultChamberZ(chamber) + AliMUONConstants::DzCh();   //!<Z of point 2, defined by being the Z of the chamber in "perfect" position 
      AliMUONTrackExtrap::ExtrapToZ(extrapTrackParam, pos2[2]);                           //!< + plus a small shift (the distance between two stations in a same chamber).
      pos2[0] = extrapTrackParam->GetNonBendingCoor();                                   //!<X of point 2, extrapolated by following the Track.                        
      pos2[1] = extrapTrackParam->GetBendingCoor();                                      //!<Y of point 2, extrapolated by following the Track.                          
      
      
      
	for (Int_t iDE = 0; iDE < nbrOfDetElt; iDE++)                    //!<Loop on all the detection element of the chamber
	  {
	    Int_t deId = (chamber + 1)*fgkOffset + iDE;                        //!<detection element Id 
	    
	    fkTransformer->Global2Local(deId, pos1[0], pos1[1], pos1[2], pos1[3], pos1[4], pos1[5]);      //!<convesrion of point 1 and 2 in the local coordinates
	    fkTransformer->Global2Local(deId, pos2[0], pos2[1], pos2[2], pos2[3], pos2[4], pos2[5]);
	    
	    CoordinatesOfMissingCluster(pos1[3], pos1[4], pos1[5], pos2[3], pos2[4], pos2[5], posMiss[0], posMiss[1]);

	    Bool_t isMissed = kFALSE;
	    if (chamber < 4)
	      isMissed = CoordinatesInDetEltSt12(deId, posMiss[0], posMiss[1]);
	    else
	      isMissed = CoordinatesInDetEltSt345(deId, posMiss[0], posMiss[1]);

	    if (isMissed)
	      FillTTHistos(chamber, deId, posMiss[0], posMiss[1]);
	  }
    }
}



//_____________________________________________________________________________
void AliCheckMuonDetEltResponse::CoordinatesOfMissingCluster(Double_t x1, Double_t y1, Double_t z1,
							     Double_t x2, Double_t y2, Double_t z2,
							     Double_t& x, Double_t& y) const
{
  //
  //Compute the coordinates of the missing cluster.
  //There are defined by the intersection between the straigth line joining two extrapolated points (1 and 2) and the detection element plane.
  //In the local coordinates, this means Z=0 in the parametric equation of the line.
  //

  Double_t t = 0;
  t = - z1 / (z2 - z1);
  
  x = t * (x2 - x1) + x1;
  y = t * (y2 - y1) + y1;
}


//_____________________________________________________________________________
Bool_t AliCheckMuonDetEltResponse::CoordinatesInDetEltSt345(Int_t DeId, Double_t x, Double_t y)
{
  //
  //Return kTRUE if the coordinates are in the Detection Element, for station 3, 4 and 5.
  //This is done by checking if a pad correspond to the (x, y) position.
  //  

  AliMpPad pad1;
  AliMpPad pad2;

  AliMpSlatSegmentation *segm1 = new AliMpSlatSegmentation(AliMpSegmentation::Instance(kFALSE)->GetSlat(DeId, AliMp::kCath0));
  AliMpSlatSegmentation *segm2 = new AliMpSlatSegmentation(AliMpSegmentation::Instance(kFALSE)->GetSlat(DeId, AliMp::kCath1));
  pad1 = segm1->PadByPosition(x, y, kFALSE);
  pad2 = segm2->PadByPosition(x, y, kFALSE);
 
  if (pad1.IsValid() && pad2.IsValid())
    return kTRUE;
  else
    return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliCheckMuonDetEltResponse::CoordinatesInDetEltSt12(Int_t DeId, Double_t x, Double_t y)
{
  //Return kTRUE if the coordinates are in the Detection Element, for station 1 and 2.
  //This is done by checking if a pad correspond to the (x, y) position.
  
  AliMpPad pad1;
  AliMpPad pad2;

  AliMpSectorSegmentation *segm1 = new AliMpSectorSegmentation(AliMpSegmentation::Instance(kFALSE)->GetSector(DeId, AliMp::kCath0));
  AliMpSectorSegmentation *segm2 = new AliMpSectorSegmentation(AliMpSegmentation::Instance(kFALSE)->GetSector(DeId, AliMp::kCath1));
  pad1 = segm1->PadByPosition(x, y, kFALSE);
  pad2 = segm2->PadByPosition(x, y, kFALSE);
 
  if (pad1.IsValid() && pad2.IsValid())
    return kTRUE;
  else
    return kFALSE;
}
