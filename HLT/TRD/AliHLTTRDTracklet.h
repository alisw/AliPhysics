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
  // void ReadClustersFromMemory(void *input);
  
 private:
  AliHLTTRDTracklet(const AliHLTTRDTracklet&);
  AliHLTTRDTracklet& operator=(const AliHLTTRDTracklet&);
  void InitArrays();

  /* Defenitely need */
  UInt_t         fN;                     // number of clusters attached/used/shared
  Float_t        fdX;                    // length of time bin
  Float_t        fYref[2];               //  Reference y
  Float_t        fZref[2];               //  Reference z
  Float_t        fS2Y;                   //  "Robust" sigma in Y - line fit
  Float_t        fPt;                    //  Momentum estimate for tracklet [GeV/c]
 
  /* Probably need */
  Float_t        fPad[3];                //  local pad definition : length/width/tilt 
  Float_t        fX0;                    //  X0 position
  Float_t        fYfit[2];               //  Y fit position +derivation
  Float_t        fZfit[2];               //  Z fit position
  Float_t        fC;                     //  Curvature
  Float_t        fChi2;                  //  Global chi2
  Float_t        fProb[AliPID::kSPECIES];// PID probabilities
  Short_t        fDet;                   //  TRD detector

  /* Not needed */
  //Int_t          fLabels[3];            //  Labels
  //Float_t        fX[knTimebins];        //! X position
  //Float_t        fY[knTimebins];        //! Y position
  //Float_t        fZ[knTimebins];        //! Z position
  //Float_t        fYfitR[2];             //  Y fit position +derivation
  //Float_t        fZfitR[2];             //  Z fit position
  //Float_t        fMeanz;                //  Mean vaue of z
  //Float_t        fZProb;                //  Max probbable z
  //Int_t          fFreq;                 //  Frequency
  //Int_t          fNChange;              //  Change z counter
  //Float_t        fMPads;                //  Mean number of pads per cluster

  //Float_t        fCC;                   //  Curvature with constrain
  //Float_t        fChi2Z;                //  Global chi2

  AliHLTUInt16_t fCount;                  // Number of clusters saved in the open array
  AliHLTUInt32_t fSize;                   // Size of the tracklet with clusters in the memory

#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTRDCluster fClusters[1];                         // Open array of clusters and their index
#else
  AliHLTTRDCluster fClusters[0];                         // Open array of clusters and their index
#endif

};

#endif
