//$Id$

#ifndef ALIHLTGLOBALTRACKMATCHER_H
#define ALIHLTGLOBALTRACKMATCHER_H


//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTGlobalTrackMatcher.h
    @author Svein Lindal (svein.lindal@fys.uio.no)
    @date   
    @brief  The HLT class matching TPC tracks to calorimeter clusters
*/



#include "AliHLTLogging.h"
#include "AliESDtrack.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TVector3.h"
#include "AliTrackerBase.h" // Marcel: for EMCal track-matching

class AliHLTGlobalTrackMatcher : public AliHLTLogging{

public:
  AliHLTGlobalTrackMatcher();

  /** destructor */
  virtual ~AliHLTGlobalTrackMatcher();

  //Main function, loops over tracks and calls appropriate functions to establish matches
  template <class T>
  Int_t Match( TObjArray * trackArray, vector<T*>  &phosClustersVector, vector<T*>  &emcalClustersVector,  Double_t bz ); 

private:
  
  void DoInit();

 
  //Loops over clusters and decides if track is a good match to any of these
  template <class T>
  Int_t MatchTrackToClusters( AliExternalTrackParam * track, vector<T*>  &clustersVector, Int_t nClusters, Float_t * bestMatch, Double_t bz); 

  template <class T>
  Int_t MatchTrackToEMCalClusters( AliExternalTrackParam * track, vector<T*>  &clustersVector, Int_t nClusters, Float_t * bestMatch, Double_t bz); // EMCal Track-Matching from recoUtils::ExtrapolateTracktoCluster

  //Add track Id to cluster's list of matching tracks
  Int_t AddTrackToCluster(Int_t tId, Int_t* clustersArray, Bool_t bestMatch, Int_t nMatches);
  Int_t AddTrackToCluster(Int_t tId, TArrayI* clustersArray, Bool_t bestMatch, Int_t nMatches);

  //Projects track to detector volume and decides if it's anywhere near calorimeter volume
  Bool_t IsTrackCloseToDetector(AliExternalTrackParam * track, Double_t bz, Double_t fMaxX, Bool_t ySign, Double_t fMaxZ, Double_t dRadius);

  // Geometrical cut off values used to decide whether track is anywhere near calorimeter volumes
  Float_t fPhosMaxZ;              // max Z track    (cm)
  Float_t fPhosMaxX;              // max X track    (cm)
  Float_t fEmcalMaxZ;             // max Z track    (cm)
  Float_t fEmcalMaxX;             // max X track    (cm)

  Float_t fMatchDistance;        // Square of maximum distance where track is considered a match to cluster (cm^2)
  Float_t fMatchDistanceEMCal; // Square of maximum distance where track is considered a match to cluster (EtaxPhi space) 

  const Double_t fPhosRadius;          // Radial position of PHOS 
  const Double_t fEmcalRadius;         // Radial position of EMCAL

  AliHLTGlobalTrackMatcher(const AliHLTGlobalTrackMatcher & );
  AliHLTGlobalTrackMatcher & operator = (const AliHLTGlobalTrackMatcher &);

  ClassDef(AliHLTGlobalTrackMatcher,1) 
};


template <class T>
Int_t AliHLTGlobalTrackMatcher::Match( TObjArray * trackArray, vector<T*>  &phosClustersVector, vector<T*> &emcalClustersVector,  Double_t bz ) {
  //See above for documentation

  Int_t nTracks = trackArray->GetEntriesFast();
  Int_t nPhosClusters = phosClustersVector.size();
  Int_t nEmcalClusters = emcalClustersVector.size(); 
  

  //HLTError("tracks phos emcal %d %d %d", nTracks, nPhosClusters, nEmcalClusters);

  //See if there are tracks and clusters to match
  if ( nTracks <= 0 ) {
    return 0;
  } else if ( (nEmcalClusters <= 0) && (nPhosClusters <= 0))  {
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

  //Loop over tracks
  for (int it = 0; it < nTracks; it++ ) {
    AliExternalTrackParam * track = static_cast<AliExternalTrackParam*>(trackArray->At(it));

    if ( IsTrackCloseToDetector(track, bz, fPhosMaxX, kFALSE, fPhosMaxZ, fPhosRadius ) ) {
      continue;//marcel test
      MatchTrackToClusters( track, phosClustersVector, nPhosClusters, bestMatchPhos, bz);

    } else if ( IsTrackCloseToDetector(track, bz, fEmcalMaxX, kTRUE, fEmcalMaxZ, fEmcalRadius ) ) {
//        MatchTrackToClusters( track, phosClustersVector, nPhosClusters, bestMatchPhos, bz);
       MatchTrackToEMCalClusters( track, emcalClustersVector, nEmcalClusters, bestMatchEmcal, bz);
    } 

  }   
    
  return 0;
} 

//MARCEL 
template <class T>
Int_t AliHLTGlobalTrackMatcher::MatchTrackToEMCalClusters( AliExternalTrackParam * track,  vector<T*> &clustersVector, Int_t nClusters, Float_t * bestMatch, Double_t  /*bz */) {
  
  //See header file for documentation
  Int_t iResult = 0;
 
  Float_t clusterPosition[3];
  //Double_t trackPosition[3];
  
  for(int ic = 0; ic < nClusters; ic++) {
    
     T * cluster = clustersVector.at(ic);

/* Comment:July 2011
   The lines below correspond to the method for Track-matching from RecoUtils:ExtrapolatetoCluster()
   In principle, the method ExtrapolateToCluster should be called directly from RecoUtils. The problems is that this method requires AliVCluster
   which would have to be created inside the GlobalMatcher, since the information that this class obtains from the Cluster comes from the
   data struct. In order to avoid the whole creation of a AliVCluster object, the code from RecoUtils
   was brought here in the same way it is written there.
*/ 
    if(cluster->E()<1.)continue;  
    if(track->Pt()<1.)continue;

    // This is a rough cut in acceptance, since it is pointless to try matching with tracks far from the EMCal
    if(track->Eta()<-.9)continue;
    if(track->Eta()>.9)continue;
    if(track->Phi()>4.)continue; 
    if(track->Phi()<.9)continue;
    
    cluster->GetPosition(clusterPosition);
    TVector3 VClsPos(clusterPosition[0], clusterPosition[1], clusterPosition[2]); //MARCEL  
    AliExternalTrackParam *trkParam = new AliExternalTrackParam(*track);//Retrieve the starting point every time before the extrapolation
    Double_t trkPos[3];
    TVector3 vec(clusterPosition[0],clusterPosition[1],clusterPosition[2]);
    Double_t alpha =  ((int)(vec.Phi()*TMath::RadToDeg()/20)+0.5)*20*TMath::DegToRad();
    vec.RotateZ(-alpha); //Rotate the cluster to the local extrapolation coordinate system
    trkParam->Rotate(alpha); //Rotate the track to the same local extrapolation system
    Float_t fMass=0.139;
    Float_t fStep=100.;
    if(!AliTrackerBase::PropagateTrackToBxByBz(trkParam, vec.X(), fMass, fStep,kFALSE, 0.8, -1)) return kFALSE; 

    trkParam->GetXYZ(trkPos); //Get the extrapolated global position
    TVector3 clsPosVec(clusterPosition[0],clusterPosition[1],clusterPosition[2]);
    TVector3 trkPosVec(trkPos[0],trkPos[1],trkPos[2]);
    Float_t clsPhi = (Float_t)clsPosVec.Phi();
    if(clsPhi<0) clsPhi+=2*TMath::Pi();
    Float_t trkPhi = (Float_t)trkPosVec.Phi();
    if(trkPhi<0) trkPhi+=2*TMath::Pi();
    Double_t tmpPhi = clsPhi-trkPhi;  // track cluster matching
    Double_t tmpEta = clsPosVec.Eta()-trkPosVec.Eta();  // track cluster matching
    Double_t tmpR=TMath::Sqrt(tmpEta*tmpEta + tmpPhi*tmpPhi);//MARCEL
    Double_t match = tmpR;

    if( match > fMatchDistanceEMCal )continue;
    
    if (match < bestMatch[ic]) {
      bestMatch[ic] = match;
      cluster->SetEmcCpvDistance(TMath::Sqrt(match));
      Double_t dx=tmpPhi;
      Double_t dz=tmpEta;  
      cluster->SetTrackDistance(dx,dz);
    }
    //Add track to cluster's array of matching tracks
    Int_t nTracksMatched = cluster->GetNTracksMatched();
    iResult = AddTrackToCluster(track->GetID(), cluster->GetTracksMatched(), match < bestMatch[ic], nTracksMatched);
  }
  
  return iResult;
}
//end of MARCEL TEST


template <class T>
Int_t AliHLTGlobalTrackMatcher::MatchTrackToClusters( AliExternalTrackParam * track, vector<T*>  &clustersVector, Int_t nClusters, Float_t * bestMatch, Double_t bz) {
  
  //See header file for documentation
  Int_t iResult = 0;
 
  Float_t clusterPosition[3];
  Double_t trackPosition[3];
  
  for(int ic = 0; ic < nClusters; ic++) {
    
    T * cluster = clustersVector.at(ic);
    
    //Get cluster global coordinates
    cluster->GetPosition(clusterPosition);

    Double_t rCluster = TMath::Sqrt(clusterPosition[0]*clusterPosition[0] + clusterPosition[1]*clusterPosition[1]);      

    //Rotate tracking system to the angle of the cluster
    TVector3 cVec(clusterPosition);
    if (! (track->Rotate(cVec.Phi())) ) {
      continue;
    }
   
    if(! (track->GetXYZAt(rCluster, bz, trackPosition)) ) {
      continue;
    }

    //Calculate track - cluster residuals
    Double_t match = 0;
    for(int i = 0; i < 3; i++) {
      Double_t dd = trackPosition[i] - clusterPosition[i];
      match += dd*dd;
    }
 

    if( match > fMatchDistance  )  {     
      continue;
    }

    if (match < bestMatch[ic]) {
      
      bestMatch[ic] = match;
      cluster->SetEmcCpvDistance(TMath::Sqrt(match));
      
      Double_t dx = trackPosition[0] - clusterPosition[0];
      Double_t dy = trackPosition[1] - clusterPosition[1];
      Double_t dz = trackPosition[2] - clusterPosition[2];
      cluster->SetTrackDistance( ((dx > 0) ? 1 : -1 )*TMath::Sqrt(dx*dx + dy*dy), dz);
    }
    
    //Add track to cluster's array of matching tracks
    Int_t nTracksMatched = cluster->GetNTracksMatched();
    iResult = AddTrackToCluster(track->GetID(), cluster->GetTracksMatched(), match < bestMatch[ic], nTracksMatched);
  }
  
  return iResult;
}

#endif
