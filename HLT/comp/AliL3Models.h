#ifndef AliL3Models_H
#define AliL3Models_H

#include "AliL3RootTypes.h"

const Int_t MaxNClusters = 32;

struct AliL3ClusterModel {
  Bool_t fPresent;
  Float_t fDTime;
  Float_t fDPad;
  Float_t fDCharge;
  Float_t fDSigmaY2;
  Float_t fDSigmaZ2;
  UInt_t fNPads;
#ifdef do_mc
  Int_t fTrackID[3];
#endif
};
typedef struct AliL3ClusterModel AliL3ClusterModel;

struct AliL3TrackModel {
  Float_t fKappa;
  Float_t fFirstPointX;
  Float_t fFirstPointY;
  Float_t fFirstPointZ;
  Float_t fTgl;
  Float_t fPsi;
  Short_t fLength;
  Short_t fClusterCharge;
  Char_t fNClusters;
};
typedef struct AliL3TrackModel AliL3TrackModel;

#endif
