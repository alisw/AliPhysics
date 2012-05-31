//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDTRACKLET_H
#define ALIHLTTRDTRACKLET_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliTRDseedV1.h"
#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTTRDCluster.h"

class AliHLTTRDTracklet
{
 public:
  AliHLTTRDTracklet();
  AliHLTTRDTracklet(const AliTRDseedV1* const inTracklet);
  
  void CopyDataMembers(const AliTRDseedV1* const inTracklet);
  void ExportTRDTracklet(AliTRDseedV1* const outTracklet) const;
  AliHLTUInt8_t *GetEndPointer() const // Returns pointer to the end of the tracklet
    { return ((AliHLTUInt8_t *)this + fSize); };
  AliHLTUInt32_t GetSize() const { return fSize; };
  void Print(Bool_t printClusters = kTRUE) const;
  static AliHLTUInt32_t SaveAt(AliHLTUInt8_t *const block, const AliTRDseedV1* const inTracklet);
  static AliHLTUInt32_t LoadFrom(AliTRDseedV1 *const outTracklet, const AliHLTUInt8_t *const block);
  
 private:
  AliHLTTRDTracklet(const AliHLTTRDTracklet&);
  AliHLTTRDTracklet& operator=(const AliHLTTRDTracklet&);
  void InitArrays();

  /* Defenitely need */
  UInt_t         fN;                     // number of clusters attached/used/shared
  Float_t        fdX;                    // length of time bin
  Float_t        fYref[2];               // Reference y
  Float_t        fZref[2];               // Reference z
  Float_t        fS2Y;                   // estimated resolution in the r-phi direction
  Float_t        fPt;                    // Momentum estimate for tracklet [GeV/c]
 
  /* Probably need */
  Float_t        fPad[3];                //  local pad definition : length/width/tilt 
  Float_t        fX0;                    //  X0 position
  Float_t        fYfit[2];               //  Y fit position +derivation
  Float_t        fZfit[2];               //  Z fit position
  Float_t        fC[2];                  //  Curvature Rieman[0] Vertex[1]
  Float_t        fChi2;                  //  Global chi2
  Float_t        fProb[AliPID::kSPECIES];// PID probabilities

  /* Not needed */
  // Float_t        fExB;                    // tg(a_L) @ tracklet location
  // Float_t        fVD;                     // drift velocity @ tracklet location
  // Float_t        fT0;                     // time 0 @ tracklet location
  // Float_t        fS2PRF;                  // sigma^2 PRF for xd->0 and phi=a_L 
  // Float_t        fDiffL;                  // longitudinal diffusion coefficient
  // Float_t        fDiffT;                  // transversal diffusion coefficient
  // Float_t        fX;                      // radial position of the tracklet
  // Float_t        fY;                      // r-phi position of the tracklet
  // Float_t        fZ;                      // z position of the tracklet
  // Float_t        fS2Z;                    // estimated resolution in the z direction 
  // Float_t        fdEdx[AliTRDseedV1::kNslices];         // dE/dx measurements for tracklet
  // Float_t        fRefCov[7];              // covariance matrix of the track in the yz plane + the rest of the diagonal elements
  // Float_t        fCov[3];                 // covariance matrix of the tracklet in the xy plane

  UChar_t        fPos[AliTRDseedV1::kNclusters]; // position of the cluster in the original array

  Short_t        fDet;                   // TRD detector
  UChar_t        fBits;                  // Bits of the tracklet
  AliHLTUInt8_t  fCount;                 // Number of clusters saved in the open array
  AliHLTUInt32_t fSize;                  // Size of the tracklet with clusters in the memory

#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTRDExtCluster fClusters[1];                         // Open array of clusters
#else
  AliHLTTRDExtCluster fClusters[0];                         // Open array of clusters
#endif

};

#endif
