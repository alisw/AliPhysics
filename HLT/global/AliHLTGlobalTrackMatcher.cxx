//$Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jacek Otwinowski (Jacek.Otwinowski@gsi.de)            *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTGlobalTrackMatcher.cxx
    @author Svein Lindal
    @date   
    @brief  The HLT class Matching Calorimeter clusters to TPC tracks
*/

//#include "AliTPCReconstructor.h"

#include "AliHLTGlobalTrackMatcher.h"
#include "AliExternalTrackParam.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "AliPHOSGeoUtils.h"
#include "AliESDEvent.h"

//struct AliHLTCaloClusterHeaderStruct;

#if __GNUC__>= 3
using namespace std;
#endif

ClassImp(AliHLTGlobalTrackMatcher)

AliHLTGlobalTrackMatcher::AliHLTGlobalTrackMatcher() :
  fClusterReader(NULL),
  fPHOSGeom(NULL), 
  fMaxZ(64.+10.),
  fMaxX(72. + 72.*TMath::Sin(20) + 72.*TMath::Sin(40) +10. ), 
  fMinX(-72.-10.),

  fDetRadius(-99999),
  fMatchDistanceSq(400),
  fBestMatchesArray(NULL),
  fTrackDistanceArray(NULL)
{
  fClusterReader= new AliHLTCaloClusterReader();
  fPHOSGeom = new AliPHOSGeoUtils("PHOS", "noCPV");
}

//_____________________________________________________________________________
AliHLTGlobalTrackMatcher::~AliHLTGlobalTrackMatcher()
{
  //Destructor
  if (fPHOSGeom)
    delete fPHOSGeom;
  fPHOSGeom = NULL;

  if (fBestMatchesArray)
    delete[] fBestMatchesArray;
  fBestMatchesArray = NULL;

  if (fTrackDistanceArray)
    delete[] fTrackDistanceArray;
  fTrackDistanceArray = NULL;

}

//_____________________________________________________________________________
Bool_t AliHLTGlobalTrackMatcher::Match(AliESDEvent* esdEvent, AliHLTCaloClusterHeaderStruct * clusterHeader)
{

  Double_t fNormVector[5][3] = {
    {0, 0, 0}, 
    {0, 0, 0}, 
    {0, 10, 0}, 
    {10*TMath::Cos(20), 10*TMath::Sin(20), 0}, 
    {10*TMath::Cos(40), 10*TMath::Sin(40), 0}, 
  };
  
  
  fClusterReader->SetMemory(clusterHeader);
  
  int nClusters = clusterHeader->fNClusters;
  
  if( !(nClusters> 0) )
    return false;

  fBestMatchesArray = new int[nClusters];
  fTrackDistanceArray = new Float_t[nClusters];
  
  //Position of track intersection with cylindrical plane defined by PHOS radius
  Double_t trackPos[3] = {0,0,0};
  
  Double_t bz = esdEvent->GetMagneticField();
  int nt = esdEvent->GetNumberOfTracks();
  for (int it = 0; it < nt; it++ ) {
    
    AliExternalTrackParam * track = static_cast<AliExternalTrackParam*> (esdEvent->GetTrack(it)) ;
 
    //See if track is even close to detector volume
    if (! (track->GetXYZAt(fDetRadius, bz, trackPos)) )
      continue;

    if (trackPos[3] > fMaxZ || trackPos[3] < -fMaxZ)
      continue;
  
    if (trackPos[0] > fMaxX || trackPos[0] < fMinX)
      continue;
 
    //Track is close to Detector volume
    //loop over clusters to find clusters that are fairly close to track
    
    //Set cluster reader to first cluster
    fClusterReader->SetMemory(clusterHeader);
    //Loop over clusters
    int clusterIndex = -1;
    AliHLTCaloClusterDataStruct* cluster;
    while( (cluster = fClusterReader->NextCluster()) ) {
            
      clusterIndex++;
      //Get approximate distance between cluster and track
      Double_t clusterPos[3] = {cluster->fGlobalPos[0], cluster->fGlobalPos[1], cluster->fGlobalPos[2]};
      Double_t distanceSq = 0;
      for(int i = 0; i < 3; i++) {
	Float_t distance = trackPos[i] - clusterPos[1];
	distanceSq += distance*distance;
      }
    
      if (distanceSq >( fMatchDistanceSq + 60) ) 
	continue;
      
      //We have a cluster in relatively close proximity to the track, do more accurate matching

      //Get the module where the cluster is located
      //Need CellId of one of the cells of the cluster;
      UShort_t cellId = -1;
      Double_t cellAmp = -1;
    
      //Get a (any) cell in the cluster
      fClusterReader->GetCell(cluster, cellId, cellAmp, (UInt_t) (0) );
    
      //Use cellId to find module nr;
      Int_t relNumbering[4];
      fPHOSGeom->AbsToRelNumbering(cellId, relNumbering);
    
      //Project the track to the plane defined by the module geometry and find intersection point
      track->Intersect(clusterPos, fNormVector[relNumbering[0]], bz);
    
      //Is the track closer to the cluster than previous matches?
      distanceSq = clusterPos[0]*clusterPos[0] + clusterPos[2]* clusterPos[2];
      if ( distanceSq  <  fTrackDistanceArray[clusterIndex] )  {
	fBestMatchesArray[clusterIndex] = it;
	fTrackDistanceArray[clusterIndex] = distanceSq;
      }
    }  
  }//done looping over clusters 
  
  //now set the values in the ESD, or where? Return them as what?
  

  //Reset Matching arrays
  if(fBestMatchesArray)
    delete[] fBestMatchesArray;
  fBestMatchesArray = NULL;

  if(fTrackDistanceArray)
    delete[] fTrackDistanceArray;
  fTrackDistanceArray = NULL;

}  
