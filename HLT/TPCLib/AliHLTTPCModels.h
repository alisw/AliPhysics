// @(#) $Id$
// Original: AliHLTModels.h,v 1.11 2004/05/17 16:37:19 hristov 

#ifndef AliHLTTPCModels_H
#define AliHLTTPCModels_H

#include "AliHLTTPCRootTypes.h"

const Int_t MaxNClusters = 32;

struct AliHLTTPCClusterModel {
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
typedef struct AliHLTTPCClusterModel AliHLTTPCClusterModel;

struct AliHLTTPCRemainingCluster {
  Float_t fPad;   
  Float_t fTime;  
  Float_t fSigmaY2;  
  Float_t fSigmaZ2;  
  UShort_t fCharge;
};
typedef struct AliHLTTPCRemainingCluster AliHLTTPCRemainingCluster;

struct AliHLTTPCRemainingRow {
  Byte_t fPadRow;       //1 byte
  UShort_t fNClusters;  //2 bytes
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTPCRemainingCluster fClusters[1];
#else
  AliHLTTPCRemainingCluster fClusters[0];
#endif
};
typedef struct AliHLTTPCRemainingRow AliHLTTPCRemainingRow;

struct AliHLTTPCTrackModel {//5 independent parameters is needed to encode the helix:
  Float_t fKappa; //Curvature
  Float_t fPhi;   //Azimuthal angle of DCAO (distance of closest approach to origo)
  Float_t fD;     //radius of DCA0
  Float_t fZ0;    //z-coordinate of DCA0
  Float_t fTgl;   //tan of dipangle
};
typedef struct AliHLTTPCTrackModel AliHLTTPCTrackModel;

#endif
