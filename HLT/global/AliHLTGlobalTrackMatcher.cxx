//$Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Svein Lindal     (slindal@fys.uio.no)                 *
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

#include "AliHLTGlobalTrackMatcher.h"
#include "AliESDCaloCluster.h"

#if __GNUC__>= 3
using namespace std;
#endif

ClassImp(AliHLTGlobalTrackMatcher)


AliHLTGlobalTrackMatcher::AliHLTGlobalTrackMatcher() :
  fPhosMaxZ(0),
  fPhosMaxX(0),
  fEmcalMaxZ(0),
  fEmcalMaxX(0),
  fPhosRadius(460),
  fEmcalRadius(448),
  fMatchDistance(0),
  fPhosClustersArray(NULL), 
  fEmcalClustersArray(NULL)
{
  //Default constructor

  DoInit();
}

//_____________________________________________________________________________
AliHLTGlobalTrackMatcher::~AliHLTGlobalTrackMatcher()
{
  //Destructor

}

void AliHLTGlobalTrackMatcher::DoInit( ) {

  fMatchDistance = 40*40;

  fPhosMaxX = 355 + TMath::Sqrt(fMatchDistance) + 30;
  fPhosMaxZ = 64.+ TMath::Sqrt(fMatchDistance) + 30;

  fEmcalMaxZ = 350 + TMath::Sqrt(fMatchDistance) + 30;
  fEmcalMaxX = 3000;

  fPhosClustersArray = new TRefArray();
  fEmcalClustersArray = new TRefArray();

}

Bool_t AliHLTGlobalTrackMatcher::Match( AliESDEvent* event ){

  Int_t nTracks = event->GetNumberOfTracks();

  //Fill arrays with clusters
  Int_t nPhosClusters = event->GetPHOSClusters(fPhosClustersArray);
  Int_t nEmcalClusters = event->GetEMCALClusters(fEmcalClustersArray);

  if ( ( nTracks <= 0 ) || ( (nEmcalClusters <= 0) && (nPhosClusters <= 0)) ) {
    //printf(Form("No tracks or clusters in event"));
    return 0;
  }
  

  Float_t bestMatchPhos[nPhosClusters];   
  for(int ic = 0; ic < nPhosClusters; ic++) {
    bestMatchPhos[ic] = 999999;
  }

  Float_t bestMatchEmcal[nEmcalClusters];    
  for(int ic = 0; ic < nEmcalClusters; ic++) {
    bestMatchEmcal[ic] = 999999;
  }


  Double_t trackPosition[3];
  
  //Loop over tracks
  for (int it = 0; it < nTracks; it++ ) {

    AliESDtrack * track = event->GetTrack(it);

    Double_t bz = event->GetMagneticField();
    
    //Does track pass close to either calorimeter volume
    if ( IsTrackCloseToDetector(track, bz, trackPosition, fPhosMaxX, kFALSE, fPhosMaxZ ) ) {
      MatchTrackToClusters( track, fPhosClustersArray, nPhosClusters, bestMatchPhos, bz);
    
    } else if ( IsTrackCloseToDetector(track, bz, trackPosition, fEmcalMaxX, kTRUE, fEmcalMaxZ ) ) {
      MatchTrackToClusters( track, fEmcalClustersArray, nEmcalClusters, bestMatchEmcal, bz);
    } 
  } // track loop 


  return 0;
} 



void AliHLTGlobalTrackMatcher::MatchTrackToClusters( AliESDtrack * track, TRefArray * clustersArray, Int_t nClusters, Float_t * bestMatch, Double_t bz) {
    
  //loop over clusters to find matches with track
 
  Float_t clusterPosition[3];
  Double_t trackPosition[3];
 
  for(int ic = 0; ic < nClusters; ic++) {
       
    AliESDCaloCluster * cluster = static_cast<AliESDCaloCluster*>(clustersArray->At(ic));
      
    //Get Cluster Global Coordinates
    cluster->GetPosition(clusterPosition);

    //Get track postion at radius of cluster
    Double_t rCluster = TMath::Sqrt(clusterPosition[0]*clusterPosition[0] + clusterPosition[1]*clusterPosition[1]);      
    if (! (track->GetXYZAt(rCluster, bz, trackPosition)) ) {
      cout << "XXXXXXXXXXXXX  BALLE Track reached detector but not cluster, needs to be handled no? !!!!!!" << endl;
    }
      
    Float_t match = 0;
    Float_t dxyz = 0;
    for(int i = 0; i < 3; i++) {
      dxyz = trackPosition[i] - clusterPosition[i];
      match += dxyz*dxyz;
    }
      
      
    if( match > fMatchDistance  )  {     
 	   continue;
    }


    //Track is close to cluster, add to cluster's array of matching tracks
    TArrayI *matchedTracksArray = cluster->GetTracksMatched();
    matchedTracksArray->Set(matchedTracksArray->GetSize() + 1);
    if ( match  < bestMatch[ic] )  {
      bestMatch[ic] = match;
      matchedTracksArray->AddAt(matchedTracksArray->At(0), matchedTracksArray->GetSize() - 1);
      matchedTracksArray->AddAt(track->GetID(), 0);
    } else {
      matchedTracksArray->AddAt(track->GetID(), matchedTracksArray->GetSize() - 1);
    }
  } //cluster loop
}


Bool_t AliHLTGlobalTrackMatcher::IsTrackCloseToDetector(AliESDtrack * track, Double_t bz, Double_t * trackPosition, Double_t fMaxX, Bool_t ySign, Double_t fMaxZ) {

  //See header file for documentation
  
  
  //Get track instersection with cylinder defined by detector radius
  if (! (track->GetXYZAt(fPhosRadius, bz, trackPosition)) ) 
    return kFALSE;
  
 //Positive y for EMCAL, negative for PHOS
  if(ySign) {
    if (trackPosition[1] < 0 ) 
      return kFALSE;
  } else {
    if (trackPosition[1] > 0 ) 
      return kFALSE;
  }
  
  
  if ( (TMath::Abs(trackPosition[2]) > fMaxZ) ) 
    return kFALSE;
  
  if (TMath::Abs(trackPosition[0]) > fMaxX )
    return kFALSE;
  
  
  
  return kTRUE;  
}


// Bool_t AliHLTGlobalTrackMatcher::IsTrackCloseToEmcal(AliESDtrack * track, Double_t bz, Double_t * trackPosition) {
  
//   //Get track instersection with cylinder defined by detector radius
//   if (! (track->GetXYZAt(fEmcalRadius, bz, trackPosition)) ) 
//     return kFALSE;
    

//   //See header file for documentation



//   return kTRUE;

//   if (trackPosition[1] < 0 ) 
//     return kFALSE;
  

//   if ( TMath::Abs(trackPosition[2]) > fEmcalMaxZ ) 
//     return kFALSE;
  
//   if ( TMath::Abs(trackPosition[0]) > fEmcalMaxX )
//     return kFALSE;
  
  
//   //
//   return kTRUE;
// }

