//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDCLUSTER_H
#define ALIHLTTRDCLUSTER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTRDCluster.h
    @author Thedoor Rascanu
    @date   
    @brief  A datacontainer for clusters fitting component for the HLT. 
*/

#include "AliTRDcluster.h"
#include "AliHLTDataTypes.h"

class AliHLTTRDCluster
{
public:
  AliHLTTRDCluster();
  AliHLTTRDCluster(const AliTRDcluster* const inCluster);
  void ExportTRDCluster(AliTRDcluster* const outCluster) const;
  static AliHLTUInt32_t SaveAt(AliHLTUInt8_t *const block, const AliTRDcluster* const inClust);
  static AliHLTUInt32_t LoadFrom(AliTRDcluster *const outClust, const AliHLTUInt8_t *const block);
  
private:
  // From AliTRDcluster
  UInt_t   fSignals;        // Signals in the cluster

  UChar_t  fPadCol;         // Central pad number in column direction 
  UChar_t  fPadRow;         // Central pad number in row direction 
  UChar_t  fPadTime;        // Uncalibrated time bin number 
  UChar_t  fBits;           // Bits of the cluster
};

// disable warnings to avoid
// warning: base class ‘class ...’ has a non-virtual destructor
#if defined __GNUC__
#if __GNUC__ == 4 && __GNUC_MINOR__ > 3
#pragma GCC diagnostic ignored "-Weffc++"
#else
#pragma GCC system_header 
#endif
#elif defined __SUNPRO_CC
#pragma disable_warn
#elif defined _MSC_VER
#pragma warning(push, 1)
#endif

class AliHLTTRDExtCluster: public AliHLTTRDCluster
{
 public:
  AliHLTTRDExtCluster();
  AliHLTTRDExtCluster(const AliTRDcluster* const inCluster);
  void ExportTRDCluster(AliTRDcluster* const outCluster) const;
  void Print() const;

 private:
  // From AliCluster
  Float_t  fX;             // X of the cluster in the tracking c.s.
  Float_t  fY;             // Y of the cluster in the tracking c.s.
  Float_t  fZ;             // Z of the cluster in the tracking c.s.

  // UChar_t  fClusterMasking; // Bit field containing cluster status information
  // Char_t   fLocalTimeBin;   // T0-calibrated time bin number
  // UChar_t  fNPads;          //  Number of pads in cluster 
  // Float_t  fCenter;         //  Center of the cluster relative to the pad  

};

#if defined __GNUC__
#if __GNUC__ == 4 && __GNUC_MINOR__ > 3
#pragma GCC diagnostic warning "-Weffc++"
#endif
#elif defined __SUNPRO_CC
#pragma enable_warn
#elif defined _MSC_VER
#pragma warning(pop)
#endif

struct AliHLTTRDClustersArray {
#ifdef HAVE_NOT_ALITRD_CLUSTERIZER_r42837
  typedef AliHLTTRDExtCluster cluster_type;
#else
  typedef AliHLTTRDCluster cluster_type;
#endif
  AliHLTTRDClustersArray(Int_t det):fDetector(det),fCount(0){}
  Short_t  fDetector;
  UShort_t fCount;
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC) || defined(__clang__)
  cluster_type fCluster[1];
#else
  cluster_type fCluster[];
#endif
};

#endif
