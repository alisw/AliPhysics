//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCCLUSTERACCESSHLTOUT_H
#define ALIHLTTPCCLUSTERACCESSHLTOUT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCClusterAccessHLTOUT.h
/// @author Matthias Richter
/// @date   2011-06-06
/// @brief  Interface to HLT TPC clusters
///

#include "TObject.h"

class AliTPCClustersRow;
class AliHLTOUT;
class TClonesArray;

/**
 * @class AliHLTTPCClusterAccessHLTOUT
 * Generator for TPC cluster array from HLT TPC clusters in the HLTOUT
 * data stream.
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCClusterAccessHLTOUT : public TObject
{
 public:
  /** standard constructor */
  AliHLTTPCClusterAccessHLTOUT();
  /** destructor */
  ~AliHLTTPCClusterAccessHLTOUT();

  /// inherited from TObject: abstract command interface
  virtual void        Execute(const char *method,  const char *params, Int_t *error=0);

  /// inherited from TObject: return the cluster array if name id "clusterarray"
  virtual TObject    *FindObject(const char *name) const;

  /// inherited from TObject: cleanup
  virtual void        Clear(Option_t * option ="");

  /// inherited from TObject
  virtual void        Print(Option_t *option="") const;

  /// process the cluster data block of various formats from HLTOUT
  int ProcessClusters();

  /// process the cluster data block {CLUSTERS:TPC } from HLTOUT
  int ReadAliHLTTPCClusterData(AliHLTOUT* pHLTOUT, TClonesArray* pClusters) const;

 private:
  /// copy constructor prohibited
  AliHLTTPCClusterAccessHLTOUT(const AliHLTTPCClusterAccessHLTOUT&);
  /// assignment operator prohibited
  AliHLTTPCClusterAccessHLTOUT& operator=(const AliHLTTPCClusterAccessHLTOUT&);

  TClonesArray* fClusters;

  ClassDef(AliHLTTPCClusterAccessHLTOUT, 0)
};
#endif
