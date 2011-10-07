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
#include "AliTPCclusterMI.h"
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
 * - param="sector=sectorno"
 * 'sectorno' specifies sector number in the offline code, range 0 and 71,
 * enumerating first the 36 inner (partitions 0+1)  and then 36 outer sectors
 * (partitions 2-5).<br>
 * If the error pointer parameter is provided the result code is returned
 * - >=0 success, number of clusters
 * - -ENODATA  no data in HLTOUT
 * - -EINVAL   invalid parameter/argument
 * - -ENOMEM   memory allocation failed
 * - -EACCESS  no access to HLTOUT
 * - -NODEV    internal error, can not get AliHLTSystem
 * - -ENOBUFS  internal error, can not get cluster array
 * 
 * Command 'verbosity=level' sets the verbositylevel which is default 0
 * (no info output).
 *
 * <pre>
 *     pClusterAccess->Execute("read", param);
 *     TObject* pClusterAccess->FindObject("clusterarray");
 * </pre>
 *
 * After processing the loop of sectors, the instance should be cleaned.
 * <pre>
 *     pClusterAccess->Clear("event");
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

  /**
   * @class AliTPCclusterMIContainer
   * Cluster read interface for offline.
   * The class implements the interface to be used in the decoding
   * of compressed TPC data.
   */
  class AliTPCclusterMIContainer {
  public:
    AliTPCclusterMIContainer();
    virtual ~AliTPCclusterMIContainer();

    struct AliClusterIdBlock {
      AliClusterIdBlock() : fIds(NULL), fSize(0) {}
      AliHLTUInt32_t* fIds; //!
      AliHLTUInt32_t  fSize; //!
    };

    class iterator {
    public:
      iterator() : fClusterNo(-1), fData(NULL), fCluster(NULL), fClusterId(kAliHLTVoidDataSpec), fRowOffset(0) {}
      iterator(AliTPCclusterMIContainer* pData) : fClusterNo(-1), fData(pData), fCluster(NULL), fClusterId(fData?fData->GetClusterId(fClusterNo):kAliHLTVoidDataSpec), fRowOffset(0) {}
      iterator(const iterator& other) : fClusterNo(other.fClusterNo), fData(other.fData), fCluster(other.fCluster), fClusterId(other.fClusterId), fRowOffset(other.fRowOffset) {}
      iterator& operator=(const iterator& other) {
	fClusterNo=other.fClusterNo; fData=other.fData; fCluster=other.fCluster, fClusterId=other.fClusterId; fRowOffset=other.fRowOffset; return *this;
      }
      ~iterator() {}

      void SetPadRow(int row)          {if (fCluster) fCluster->SetRow(row-fRowOffset);}
      void SetPad(float pad) 	       {if (fCluster) fCluster->SetPad(pad);}
      void SetTime(float time) 	       {if (fCluster) fCluster->SetTimeBin(time);}
      void SetSigmaY2(float sigmaY2)   {if (fCluster) fCluster->SetSigmaY2(sigmaY2);}
      void SetSigmaZ2(float sigmaZ2)   {if (fCluster) fCluster->SetSigmaZ2(sigmaZ2);}
      void SetCharge(unsigned charge)  {if (fCluster) fCluster->SetQ(charge);}
      void SetQMax(unsigned qmax)      {if (fCluster) fCluster->SetMax(qmax);}

      // switch to next cluster
      iterator& Next(int slice, int partition);

    private:
      int fClusterNo; //! cluster no in the current block
      AliTPCclusterMIContainer* fData; //! pointer to actual data
      AliTPCclusterMI* fCluster; //! pointer to current cluster
      AliHLTUInt32_t fClusterId; //! id of the cluster, from optional cluster id blocks
      int fRowOffset;  //! row offset for current partition
    };

    /// iterator of remaining clusters block of specification
    iterator& BeginRemainingClusterBlock(int count, AliHLTUInt32_t specification);
    /// iterator of track model clusters
    iterator& BeginTrackModelClusterBlock(int count);

    /// add cluster mc data block
    int AddClusterMCData(const AliHLTComponentBlockData* pDesc);
    /// add cluster id block for remaining or track model clusters
    int AddClusterIds(const AliHLTComponentBlockData* pDesc);
    /// get the cluster id from the current cluster id block (optional)
    AliHLTUInt32_t GetClusterId(int clusterNo) const;

    /// internal cleanup
    virtual void  Clear(Option_t * option="");
    /// get the cluster array for a sector
    TObjArray* GetSectorArray(unsigned sector) const;
    /// print info
    virtual void Print(Option_t *option=NULL) const;

  protected:
    /// load next cluster from array of the sepcific sector
    AliTPCclusterMI* NextCluster(int slice, int partition);
    /// set MC data for the cluster
    int SetMC(AliTPCclusterMI* cluster, AliHLTUInt32_t clusterId);

  private:
    AliTPCclusterMIContainer(const AliTPCclusterMIContainer&);
    AliTPCclusterMIContainer& operator=(const AliTPCclusterMIContainer&);

    vector<TClonesArray*> fClusterArrays; //! cluster arrays per sector (offline notation 0-71)
    vector<AliClusterIdBlock> fRemainingClusterIds; //! clusters ids for remaining cluster ids
    AliClusterIdBlock fTrackModelClusterIds; //! cluster ids for track model clusters
    AliClusterIdBlock* fCurrentClusterIds; //! id block currently active in the iteration
    vector<const AliHLTTPCClusterMCData*> fClusterMCData; //! references to MC data blocks
    iterator fIterator; //!
  };

 private:
  /// copy constructor prohibited
  AliHLTTPCClusterAccessHLTOUT(const AliHLTTPCClusterAccessHLTOUT&);
  /// assignment operator prohibited
  AliHLTTPCClusterAccessHLTOUT& operator=(const AliHLTTPCClusterAccessHLTOUT&);

  enum EOptions {
    // skip the track clusters
    kSkipTrackClusters = BIT(15),
    // skip the partition (remaining) clusters
    kSkipPartitionClusters = BIT(16)
  };

  int fVerbosity; //! verbosity level
  AliTPCclusterMIContainer* fClusters; //! cluster container
  int fCurrentSector; //! current sector

  ClassDef(AliHLTTPCClusterAccessHLTOUT, 0)
};
#endif
