// @(#) $Id$

#ifndef AliHLTModels_H
#define AliHLTModels_H

#include "AliHLTRootTypes.h"

const Int_t MaxNClusters = 32;

struct AliHLTClusterModel {
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
typedef struct AliHLTClusterModel AliHLTClusterModel;
typedef AliHLTClusterModel AliL3ClusterModel;

struct AliHLTRemainingCluster {
  Float_t fPad;   
  Float_t fTime;  
  Float_t fSigmaY2;  
  Float_t fSigmaZ2;  
  UShort_t fCharge;
};
typedef struct AliHLTRemainingCluster AliHLTRemainingCluster;
typedef AliHLTRemainingCluster AliL3RemainingCluster;

struct AliHLTRemainingRow {
  Byte_t fPadRow;       //1 byte
  UShort_t fNClusters;  //2 bytes
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTRemainingCluster fClusters[1];
#else
  AliHLTRemainingCluster fClusters[0];
#endif
};
typedef struct AliHLTRemainingRow AliHLTRemainingRow;
typedef AliHLTRemainingRow AliL3RemainingRow;

struct AliHLTTrackModel {//5 independent parameters is needed to encode the helix:
  Float_t fKappa; //Curvature
  Float_t fPhi;   //Azimuthal angle of DCAO (distance of closest approach to origo)
  Float_t fD;     //radius of DCA0
  Float_t fZ0;    //z-coordinate of DCA0
  Float_t fTgl;   //tan of dipangle
};
typedef struct AliHLTTrackModel AliHLTTrackModel;
typedef AliHLTTrackModel AliL3TrackModel;

#endif
