// @(#) $Id$

#ifndef AliL3Models_H
#define AliL3Models_H

#include "AliL3RootTypes.h"

const Int_t MaxNClusters = 32;

struct AliL3ClusterModel {
  Byte_t fPresent;
  Float_t fDTime;
  Float_t fDPad;
  Float_t fDCharge;
  Float_t fDSigmaY;
  Float_t fDSigmaZ;
  UInt_t fNPads;
  Short_t fSlice;
#ifdef do_mc
  Int_t fTrackID[3];
#endif
};
typedef struct AliL3ClusterModel AliL3ClusterModel;

struct AliL3RemainingCluster {
  Float_t fPad;   
  Float_t fTime;  
  Float_t fSigmaY2;  
  Float_t fSigmaZ2;  
  UShort_t fCharge;
};
typedef struct AliL3RemainingCluster AliL3RemainingCluster;

struct AliL3RemainingRow {
  Byte_t fPadRow;       //1 byte
  UShort_t fNClusters;  //2 bytes
  AliL3RemainingCluster fClusters[0];
};
typedef struct AliL3RemainingRow AliL3RemainingRow;

struct AliL3TrackModel {//5 independent parameters is needed to encode the helix:
  Float_t fKappa; //Curvature
  Float_t fPhi;   //Azimuthal angle of DCAO (distance of closest approach to origo)
  Float_t fD;     //radius of DCA0
  Float_t fZ0;    //z-coordinate of DCA0
  Float_t fTgl;   //tan of dipangle
};
typedef struct AliL3TrackModel AliL3TrackModel;

#endif
