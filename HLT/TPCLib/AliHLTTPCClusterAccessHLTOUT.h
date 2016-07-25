//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCCLUSTERACCESSHLTOUT_H
#define ALIHLTTPCCLUSTERACCESSHLTOUT_H
//* This file is property of and copyright by the ALICE Project            * 
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
#include "AliHLTTPCRawCluster.h"
#include <vector>

class AliTPCParam;
class TClonesArray;
class AliHLTTPCDataCompressionDecoder;

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
 *     } while (!pCl && gSystem->Load("libAliHLTTPC")==0);
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
  virtual ~AliHLTTPCClusterAccessHLTOUT();

  /// inherited from TObject: abstract command interface
  virtual void        Execute(const char *method,  const char *params, Int_t *error=0);

  /// inherited from TObject: return the cluster array if name id "clusterarray"
  virtual TObject    *FindObject(const char *name) const;

  /// inherited from TObject: supports writing of data to AliTPCClustersRow
  virtual void Copy(TObject &object) const;

  /// inherited from TObject: cleanup
  virtual void        Clear(Option_t * option ="");

  /// inherited from TObject
  virtual void        Print(Option_t *option="") const;

  /// process the cluster data block of various formats from HLTOUT
  int ProcessClusters(const char* params);

  /// scan parameters
  /// known: sector=<n> row=<n>
  int ScanParameters(const char* params);

  /// helper struct to store cluster in a map together with MC info
  struct AliRawClusterEntry {
    AliRawClusterEntry() : fCluster(), fMC() {}
    AliRawClusterEntry(const AliRawClusterEntry& other) : fCluster(other.fCluster), fMC(other.fMC) {}
    AliRawClusterEntry& operator=(const AliRawClusterEntry& other) {
      if (&other==this) return *this;
      fCluster=other.fCluster; fMC=other.fMC;
      return *this;
    }
    AliHLTTPCRawCluster fCluster; //! cluster
    AliHLTTPCClusterMCLabel fMC; //! MC labels 
  };
  
  typedef vector<AliRawClusterEntry> AliRawClusterEntryVector;

  /**
   * @class AliRawClusterContainer
   * Cluster read interface for offline.
   * The class implements the interface to be used in the decoding
   * of compressed TPC data. The container handles 
   */
  class AliRawClusterContainer {
  public:
    AliRawClusterContainer();
    virtual ~AliRawClusterContainer();

    struct AliClusterIdBlock {
      AliClusterIdBlock() : fIds(NULL), fSize(0) {}
      AliHLTUInt32_t* fIds; //!
      AliHLTUInt32_t  fSize; //!
    };

    class iterator {
    public:
      iterator() : fClusterNo(-1), fData(NULL), fEntry(NULL), fRowOffset(0) {}
      iterator(AliRawClusterContainer* pData) : fClusterNo(-1), fData(pData), fEntry(NULL), fRowOffset(0) {}
      iterator(const iterator& other) : fClusterNo(other.fClusterNo), fData(other.fData), fEntry(other.fEntry), fRowOffset(other.fRowOffset) {}
      iterator& operator=(const iterator& other) {
	if (this==&other) return *this;
	fClusterNo=other.fClusterNo; fData=other.fData; fEntry=other.fEntry; fRowOffset=other.fRowOffset; return *this;
      }
      virtual ~iterator() {}

      void SetPadRow(int row)          {if (fEntry ) fEntry->fCluster.SetPadRow(row-fRowOffset);}
      void SetPad(float pad) 	       {if (fEntry ) fEntry->fCluster.SetPad(pad);}
      void SetTime(float time) 	       {if (fEntry ) fEntry->fCluster.SetTime(time);}
      void SetSigmaY2(float sigmaY2)   {if (fEntry ) fEntry->fCluster.SetSigmaPad2(sigmaY2);}
      void SetSigmaZ2(float sigmaZ2)   {if (fEntry ) fEntry->fCluster.SetSigmaTime2(sigmaZ2);}
      void SetCharge(unsigned charge)  {if (fEntry ) fEntry->fCluster.SetCharge(charge);}
      void SetQMax(unsigned qmax)      {if (fEntry ) fEntry->fCluster.SetQMax(qmax);}
      void SetFlags(unsigned short flags) {if (fEntry ) fEntry->fCluster.SetFlags(flags);}
      iterator& operator=(const AliHLTTPCRawCluster& rawcluster) {if (fEntry ) {
	  memcpy(&fEntry->fCluster, &rawcluster, sizeof(AliHLTTPCRawCluster));
	  // Note: offline code uses a different convention for row offset than the online code
	  // Online: first row of readout partition
	  // Offline: first row of readout chamber(inner: partition 0-1; outer: 2-5 
	  fEntry->fCluster.fPadRow-=fRowOffset;
	} return *this;}
      void SetMC(const AliHLTTPCClusterMCLabel* pMC) {
	if (fEntry && pMC ) fEntry->fMC=*pMC;
      }

      // switch to next cluster
      iterator& Next(int slice, int partition);

    private:

      static const Int_t fkRowOffsetOuterSector; //! transient

      int fClusterNo; //! cluster no in the current block
      AliRawClusterContainer* fData; //! pointer to actual data
      AliRawClusterEntry* fEntry; //! pointer to current cluster
      int fRowOffset;  //! row offset for current partition      
    };

    /// legacy, to be removed later
    iterator& BeginRemainingClusterBlock(int count, AliHLTUInt32_t specification) {
      return BeginPartitionClusterBlock(count, specification);
    }
    /// iterator of partition clusters block of specification
    iterator& BeginPartitionClusterBlock(int count, AliHLTUInt32_t specification);
    /// iterator of track model clusters
    iterator& BeginTrackModelClusterBlock(int count);

    /// check if the container is filled
    bool HaveData() const {return fHaveData;}
    /// mark the container valid
    void MarkValid() {fHaveData=true;}
    /// internal cleanup
    virtual void  Clear(Option_t * option="");
    /// get the cluster array for a sector
    TObjArray* GetSectorArray(unsigned sector, int propagateSplitClusterFlag = 0, int markEdgeClusters = 0) const;
    /// fill the cluster array for a sector and specific row if specified
    int FillSectorArray(TClonesArray* pSectorArray, unsigned sector, int row=-1, int propagateSplitClusterFlag = 0, int markEdgeClusters = 0) const;
    /// print info
    virtual void Print(Option_t *option=NULL) const;

  protected:
    /// load next cluster from array of the sepcific sector
    AliRawClusterEntry* NextCluster(int slice, int partition);

  private:
    /// copy constructor prohibited
    AliRawClusterContainer(const AliRawClusterContainer&);
    /// assignment operator prohibited
    AliRawClusterContainer& operator=(const AliRawClusterContainer&);

    vector<AliRawClusterEntryVector*> fClusterMaps; //! cluster pointer vectors per sector (offline notation 0-71)
    TClonesArray* fSectorArray; //! current sector array of clusters provided to caller
    iterator fIterator; //!
    bool fHaveData; //! indivate that the container is filled
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
  AliRawClusterContainer* fClusters; //! cluster container
  int fCurrentSector; //! current sector
  int fCurrentRow; //! current row
  int fPropagateSplitClusterFlag; //! propagate flags for split clusters
  int fMarkEdgeClusters; //! mark edge clusters during decoding
  AliHLTTPCDataCompressionDecoder* fpDecoder; //! decoder instance
  AliTPCParam* fTPCParam; //! pointer to TPC param

  ClassDef(AliHLTTPCClusterAccessHLTOUT, 0)
};
#endif
