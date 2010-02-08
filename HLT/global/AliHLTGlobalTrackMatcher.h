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
#include "TRefArray.h"
#include "AliESDEvent.h"

class AliHLTGlobalTrackMatcher {

public:
  AliHLTGlobalTrackMatcher();

  /** destructor */
  virtual ~AliHLTGlobalTrackMatcher();


  // Main function, loops over tracks, clusters
  // Finds best match for cluster
  Bool_t Match( AliESDEvent* event );


private:

  void DoInit();

  //Projects track to detector volume and decides if it passes in the vinity or not
  Bool_t IsTrackCloseToDetector(AliESDtrack * track, Double_t bz, Double_t * trackPosition, Double_t fMaxX, Bool_t ySign, Double_t fMaxZ);

  //Fills fClustersArray with refs to the clusters, returns number of clusters
  void MatchTrackToClusters( AliESDtrack * track, TRefArray * clustersArray, Int_t nClusters, Float_t * bestMatch, Double_t bz );
    
  // Geometrical paramaters of detector volume
  Float_t fPhosMaxZ;              // max Z track    (cm)
  Float_t fPhosMaxX;              // max X track    (cm)
  Float_t fEmcalMaxZ;
  Float_t fEmcalMaxX;

  Float_t fMatchDistance;

  const Double_t fPhosRadius;
  const Double_t fEmcalRadius;         // Radial position of detector volume
  //TRefArrays that will be filled with referenced to the PHOS and EMCAL Clusters
  TRefArray * fPhosClustersArray;
  TRefArray * fEmcalClustersArray;

  AliHLTGlobalTrackMatcher(const AliHLTGlobalTrackMatcher & );
  AliHLTGlobalTrackMatcher & operator = (const AliHLTGlobalTrackMatcher &);

  ClassDef(AliHLTGlobalTrackMatcher,1) //Merging base class
};

#endif
