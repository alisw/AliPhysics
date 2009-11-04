// $Id$
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


class AliPHOSGeoUtils;

class TClonesArray;
class TTreeStream;
class TTreeSRedirector;
class AliESDEvent;
class AliHLTCaloClusterReader;
struct AliHLTCaloClusterDataStruct;
struct AliHLTCaloClusterHeaderStruct;

#include "AliHLTLogging.h"
#include "AliESDtrack.h"

/** 
 * @class AliHLTGlobalTrackMatcher
 * Global track merger for the barrel section.
 *
 * @ingroup alihlt_global
 * @author Jacek.Otwinowski@gsi.de
 */
class AliHLTGlobalTrackMatcher : public AliHLTLogging {

public:
  AliHLTGlobalTrackMatcher();

  /** destructor */
  virtual ~AliHLTGlobalTrackMatcher();

  // set matching parameters
  void SetParameter(Double_t maxy=1., Double_t maxz=1., Double_t maxsnp=0.05, Double_t maxtgl=0.1, Double_t signed1Pt=0.001);

  // match tracks
  Bool_t Match(AliESDEvent *esdEvent, AliHLTCaloClusterHeaderStruct * clusterHeaderStruct);

private:
  
  //Helper class reading calocluster structs.
  AliHLTCaloClusterReader * fClusterReader;

  //PHOS Geometry
  AliPHOSGeoUtils* fPHOSGeom;


  // PHOS Geometry boundaries matching parameters
  const Double_t fMaxZ;    //! max Z track    (cm)
  const Double_t fMaxX;    //! max X track    (cm)
  const Double_t fMinX;    //  min X of track (cm)

  const Double_t fDetRadius;
  const Double_t fMatchDistanceSq;
  
  //Angle of PHOS Modules to Y 
  //Float_t fPHOSAngles[5];
  int fNModules;

  Int_t *fBestMatchesArray;
  Float_t *fTrackDistanceArray;

  AliHLTGlobalTrackMatcher(const AliHLTGlobalTrackMatcher & );
  AliHLTGlobalTrackMatcher & operator = (const AliHLTGlobalTrackMatcher &);

  Int_t *fBestMatchesArray;
  Float_t *fTrackDistanceArray;

  ClassDef(AliHLTGlobalTrackMatcher,1) //Merging base class
};

#endif
