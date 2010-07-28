// @(#) $Id$
// Original: AliHLTModels.h,v 1.11 2004/05/17 16:37:19 hristov 

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#ifndef AliHLTTPCModels_H
#define AliHLTTPCModels_H

/// \author Anders Vestbo or Constantin Loizides
/// \date 7 Sep 2005
/// \brief Data structures used for TPC clusters.
///
/// Basic structures used by the AliHLTTPCModelTrack class.

#include "AliHLTTPCRootTypes.h"

const Int_t kMaxNClusters = 32;

struct AliHLTTPCClusterModel {
  Byte_t fPresent;  // Boolean flag, is the cluster present or not?
  Float_t fDTime;   // Digit time.
  Float_t fDPad;    // Digit pad.
  Float_t fDCharge; // Digit charge.
  Float_t fDSigmaY; // Sigma in Y.
  Float_t fDSigmaZ; // Sigma in X.
  UInt_t fNPads;    // Total number of pafs.
  Short_t fSlice;   // TPC slice.
#ifdef do_mc
  Int_t fTrackID[3];  // Track IDs from MC.
#endif
};
typedef struct AliHLTTPCClusterModel AliHLTTPCClusterModel;

struct AliHLTTPCRemainingCluster {
  Float_t fPad;     // Pad value.
  Float_t fTime;    // Time value.
  Float_t fSigmaY2; // Square of sigma in Y.
  Float_t fSigmaZ2; // Square of sigma in X.
  UShort_t fCharge; // Charge.
};
typedef struct AliHLTTPCRemainingCluster AliHLTTPCRemainingCluster;

struct AliHLTTPCRemainingRow {
  Byte_t fPadRow;       //1 byte
  UShort_t fNClusters;  //2 bytes
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTPCRemainingCluster fClusters[1];  // Array of clusters.
#else
  AliHLTTPCRemainingCluster fClusters[0];  // Array of clusters.
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
