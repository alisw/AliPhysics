// @(#) $Id$

#ifndef AliHLTTPCModels_H
#define AliHLTTPCModels_H

#include "AliHLTTPCRootTypes.h"

const Int_t MaxNClusters = 32;

struct AliHLTTPCClusterModel {
  Byte_t fPresent;
  Float_t fDTime;
  Float_t fDPad;
  Float_t fDCharge;
  Float_t fDSigmaY2;
  Float_t fDSigmaZ2;
  UInt_t fNPads;
  Short_t fSlice;
#ifdef do_mc
  Int_t fTrackID[3];
#endif
};
typedef struct AliHLTTPCClusterModel AliHLTTPCClusterModel;

struct AliHLTTPCRemainingCluster {
  Float_t fY;
  Float_t fZ;
  UShort_t fCharge;
  Float_t fSigmaY2;
  Float_t fSigmaZ2;
};
typedef struct AliHLTTPCRemainingCluster AliHLTTPCRemainingCluster;

struct AliHLTTPCRemainingRow {
  Byte_t fPadRow;
  UShort_t fNClusters;
  AliHLTTPCRemainingCluster fClusters[0];
};
typedef struct AliHLTTPCRemainingRow AliHLTTPCRemainingRow;

struct AliHLTTPCTrackModel {
  Float_t fKappa;
  Float_t fFirstPointX;
  Float_t fFirstPointY;
  Float_t fFirstPointZ;
  Float_t fTgl;
  Float_t fPsi;
};
typedef struct AliHLTTPCTrackModel AliHLTTPCTrackModel;

#endif
