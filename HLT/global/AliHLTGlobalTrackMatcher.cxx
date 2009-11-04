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
    @author Jacek Otwinowski
    @date   
    @brief  The HLT global merger base class
*/

//#include "AliTPCReconstructor.h"

#include "AliHLTGlobalTrackMatcher.h"
#include "AliExternalTrackParam.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliPHOSGeoUtils.h"

#if __GNUC__>= 3
using namespace std;
#endif

ClassImp(AliHLTGlobalTrackMatcher)

AliHLTGlobalTrackMatcher::AliHLTGlobalTrackMatcher() :
  fMaxZ(64.+10.),
  fMaxX(72. + 72.*TMath::Sin(20) + 72.*TMath::Sin(40) +10. ), 
  fMinX(-72.-10.),
  fPHOSGeom(NULL), 
  fDetRadius(99999999),
  fMaxSqDistance(99999),
  fMatchingDistanceSq(9999999)
{
  //Default constructor
  
  fPHOSGeom = new AliPHOSGeoUtils("PHOS", "noCPV");
  //  fMaxX= ();
  //fMinX=-72.-10.;
}

//_____________________________________________________________________________
AliHLTGlobalTrackMatcher::~AliHLTGlobalTrackMatcher()
{
  //Destructor
}

  



//_____________________________________________________________________________
Bool_t AliHLTGlobalTrackMatcher::Match(AliHLTComponentBlockData* pBlock)
{
 
  Float_t shortestDistance = fMatchingDistance;

  //Should be moved to constructor, fPHOSAngles should be a
  Float_t fPHOSAngles[5] = {40, 20 , 0, -20 , 40};
  Double_t fNormVector[1][3] = {{0, 0, 0}};
  
  //Helper variables, must be set dynamically;
  int nt = 0;
  AliExternalTrackParam * track = NULL;
  AliHLTCaloClusterDataStruct * cluster = NULL;
  Double_t bz = 0.0;

  Double_t trackPos[3] = {0,0,0};
  
  for (int it = 0; it < nt; it++ ) {
  
 
    //Get track
 
    //See it track is even close to detector volume
    if (! (track->GetXYZAt(fDetRadius, bz, trackPos)) )
      continue;

    if (trackPos[3] > fMaxZ || trackPos[3] < -fMaxZ)
      continue;
  
    if (trackPos[0] > fMaxX || trackPos[0] < fMinX)
      continue;
 
    //Track is close to Detector volume
    //loop over clusters and find clusters that are fairly near track
    //for(int ic = 0; ic < nClusters; ic++) {
    Double_t clusterPos[3] = {cluster->fGlobalPos[0], cluster->fGlobalPos[1], cluster->fGlobalPos[2]};
    Double_t distanceSq = 0;
    
    //Find approximate distance between track and cluster
    for(int i = 0; i < 3; i++) {
      Float_t distance = trackPos[i] - clusterPos[1];
      distanceSq += distance*distance;
    }
    
    
    if (distanceSq > fMaxSqDistance ) 
      continue;

    //We have a track in relatively close proximity to a cluster, do more accurate matching

    //Get the module where the cluster is located
    //Need CellId of one of the cells of the cluster;
    int cellId = 0;
    
    Int_t relNumbering[4];
    //BALLE check formart of relnumbering!!!!
    fPHOSGeom->AbsToRelNumbering(cellId, relNumbering);
    
    //Project the track to the plane defined by the module geometry.
    track->Intersect(clusterPos, fNormVector[(int) relNumbering[0]], bz);
    
    //Is the track closer to the cluster than previous matches.
    if ( (Float_t distanceSq = clusterPos[0]*clusterPos[0] + clusterPos[2]* clusterPos[2]) < fTrackDistance[clusterIndex]*fTrackDistance[clusterIndex] ) {
      
      fBestMatches[clusterIndex] = trackIndex;
      fTrackDistance[clusterIndex] = distanceSq;
      
    }
    
  }//done looping over clusters 
  
}

  



