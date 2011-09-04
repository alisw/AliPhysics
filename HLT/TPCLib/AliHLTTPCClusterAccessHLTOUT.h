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
#include "AliHLTDataTypes.h"
#include "AliHLTTPCClusterMCData.h"
#include <map>

class AliTPCClustersRow;
class AliHLTOUT;
class TClonesArray;

typedef std::map<AliHLTUInt32_t, AliHLTTPCClusterMCLabel> AliHLTTPCClusterMCDataList;

/**
 * @class AliHLTTPCClusterAccessHLTOUT
 * Generator for TPC cluster array from HLT TPC clusters in the HLTOUT
 * data stream. It uses the TObject method interface. Combined with dynamic
 * loading, any cross dependency between the TPC code and HLT libraries
 * can be avoided.
 *
 * Creating the instance: 
 * The class is implemented in libAliHLTTPC.so and can be loaded
 * through the TClass interface (you might want to add
 * further error messages on the various error conditions).
 * <pre>
 *     TObject* pClusterAccess=NULL;
 *     TClass* pCl=NULL;
 *     ROOT::NewFunc_t pNewFunc=NULL;
 *     do {
 *       pCl=TClass::GetClass("AliHLTTPCClusterAccessHLTOUT");
 *     } while (!pCl && gSystem->Load("libAliHLTTPC.so")==0);
 *     if (pCl && (pNewFunc=pCl->GetNew())!=NULL) {
 *       void* p=(*pNewFunc)(NULL);
 *       if (p) {
 *         pClusterAccess=reinterpret_cast<TObject*>(p);
 *       }
 *     }
 * </pre>
 * 
 * Usage:
 * TObject::Execute can be used to execute commands. Command 'read'
 * will get hold on the HLTOUT data and read the clusters. The const char*
 * parameter 'param' is used to select the region.
 * - param=NULL: read all
 * - param="sector=sectorno"
 * 'sectorno' specifies sector number in the offline code, range 0 and 71,
 * enumerating first the 36 inner (partitions 0+1)  and then 36 outer sectors
 * (partitions 2-5).
 * 
 * Command 'verbosity=level' sets the verbositylevel which is default 0
 * (no info output).
 * <pre>
 *     pClusterAccess->Clear();
 *     pClusterAccess->Execute("read", param);
 *     TObject* pClusterAccess->FindObject("clusterarray");
 * </pre>
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
  int ProcessClusters(const char* params);

  /// process the cluster mc data block {CLMCINFO:TPC } from HLTOUT
  int ReadAliHLTTPCClusterMCData(AliHLTOUT* pHLTOUT, AliHLTTPCClusterMCDataList &tpcClusterLabels) const;

  /// process the cluster data block {CLUSTERS:TPC } from HLTOUT
  int ReadAliHLTTPCClusterData(AliHLTOUT* pHLTOUT, TClonesArray* pClusters, const AliHLTTPCClusterMCDataList *tpcClusterLabels=NULL) const;

  /// process the cluster data block {CLUSTRAW:TPC } from HLTOUT
  int ReadAliHLTTPCRawClusterData(AliHLTOUT* pHLTOUT, TClonesArray* pClusters, const AliHLTTPCClusterMCDataList *tpcClusterLabels);

  /// process the clusters of type {REMCLSCM:TPC } from HLTOUT
  int ReadRemainingClustersCompressed(AliHLTOUT* pHLTOUT, TClonesArray* pClusters, const AliHLTTPCClusterMCDataList *tpcClusterLabels);

  /// process clusters encoded by AliHLTDataDeflaterSimple
  int ReadAliHLTTPCRawClusterDataDeflateSimple(const AliHLTUInt8_t* pData, int dataSize,
					       int nofClusters, AliHLTUInt32_t specification,
					       TClonesArray* pClusters, const AliHLTTPCClusterMCDataList *tpcClusterLabels);

 private:
  /// copy constructor prohibited
  AliHLTTPCClusterAccessHLTOUT(const AliHLTTPCClusterAccessHLTOUT&);
  /// assignment operator prohibited
  AliHLTTPCClusterAccessHLTOUT& operator=(const AliHLTTPCClusterAccessHLTOUT&);

  int fVerbosity; //! verbosity level
  TClonesArray* fClusters; //! cluster array

  ClassDef(AliHLTTPCClusterAccessHLTOUT, 0)
};
#endif
