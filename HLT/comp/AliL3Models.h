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
  Float_t fDSigmaY2;
  Float_t fDSigmaZ2;
  UInt_t fNPads;
  Short_t fSlice;
#ifdef do_mc
  Int_t fTrackID[3];
#endif
};
typedef struct AliL3ClusterModel AliL3ClusterModel;

struct AliL3RemainingCluster {
  Float_t fY;
  Float_t fZ;
  UShort_t fCharge;
  Float_t fSigmaY2;
  Float_t fSigmaZ2;
};
typedef struct AliL3RemainingCluster AliL3RemainingCluster;

struct AliL3RemainingRow {
  Byte_t fPadRow;
  UShort_t fNClusters;
  AliL3RemainingCluster fClusters[0];
};
typedef struct AliL3RemainingRow AliL3RemainingRow;

struct AliL3TrackModel {
  Float_t fKappa;
  Float_t fFirstPointX;
  Float_t fFirstPointY;
  Float_t fFirstPointZ;
  Float_t fTgl;
  Float_t fPsi;
};
typedef struct AliL3TrackModel AliL3TrackModel;

#endif
