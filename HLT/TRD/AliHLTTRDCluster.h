//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDCLUSTER_H
#define ALIHLTTRDCLUSTER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

#include "AliTRDcluster.h"
#include "AliHLTDataTypes.h"

class AliHLTTRDCluster
{
 public:
  AliHLTTRDCluster();
  AliHLTTRDCluster(const AliTRDcluster* const inCluster);
  virtual ~AliHLTTRDCluster() {};
  void ExportTRDCluster(AliTRDcluster* const outCluster) const;
  void Print() const;
  
 private:
  // From AliCluster
  Float_t  fX;             // X of the cluster in the tracking c.s.
  Float_t  fY;             // Y of the cluster in the tracking c.s.
  Float_t  fZ;             // Z of the cluster in the tracking c.s.

  // From AliTRDcluster
  Short_t  fSignals[3];     // Signals in the cluster
  Short_t  fDetector;       // TRD detector number 
  Char_t   fLocalTimeBin;   // T0-calibrated time bin number
  UChar_t  fClusterMasking; // Bit field containing cluster status information;

  UChar_t  fPadCol;         // Central pad number in column direction 
  UChar_t  fPadRow;         // Central pad number in row direction 
  UChar_t  fPadTime;        // Uncalibrated time bin number 

  UChar_t  fBits;           // Bits of the cluster

  // UChar_t  fNPads;          //  Number of pads in cluster 
  // Float_t  fCenter;         //  Center of the cluster relative to the pad  

};

#endif
