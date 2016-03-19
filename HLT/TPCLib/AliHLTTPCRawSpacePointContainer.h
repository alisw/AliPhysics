//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCRAWSPACEPOINTCONTAINER_H
#define ALIHLTTPCRAWSPACEPOINTCONTAINER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCRawSpacePointContainer.h
/// @author Matthias Richter, Sergey Gorbunov
/// @date   2011-08-08
/// @brief  Helper class for handling of HLT TPC cluster data blocks from the
///         HW ClusterFinder
/// @note   Class is a duplicate of AliHLTTPCRawSpacePointContainer and should
///         be merged with it in a generic way

#include "AliHLTSpacePointContainer.h"
#include "AliHLTTPCRawCluster.h"
#include <map>
#include <vector>
using namespace std;

/**
 * @class AliHLTTPCRawSpacePointContainer
 * Handler class for HLT TPC hardware ClusterFinder space point data blocks.
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCRawSpacePointContainer : public AliHLTSpacePointContainer
{
 public:
  /// standard constructor
  AliHLTTPCRawSpacePointContainer(int mode=0);
  /// copy constructor
  AliHLTTPCRawSpacePointContainer(const AliHLTTPCRawSpacePointContainer& c);
  /// assignment operator
  AliHLTTPCRawSpacePointContainer& operator=(const AliHLTTPCRawSpacePointContainer& c);
  /// destructor
  ~AliHLTTPCRawSpacePointContainer();

  enum {
    kModeSingle = 0x1,
    kModeCreateMap = 0x2,
    kModeDifferentialPadTime = 0x4
  };

  virtual bool Check(AliHLTUInt32_t clusterID) const;
  virtual int GetClusterIDs(vector<AliHLTUInt32_t>& tgt) const;
  virtual const vector<AliHLTUInt32_t>* GetClusterIDs(AliHLTUInt32_t mask);
  virtual float GetX(AliHLTUInt32_t clusterID) const;
  virtual float GetXWidth(AliHLTUInt32_t clusterID) const;
  virtual float GetY(AliHLTUInt32_t clusterID) const;
  virtual float GetYWidth(AliHLTUInt32_t clusterID) const;
  virtual float GetZ(AliHLTUInt32_t clusterID) const;
  virtual float GetZWidth(AliHLTUInt32_t clusterID) const;
  virtual float GetCharge(AliHLTUInt32_t clusterID) const;
  virtual float GetQMax(AliHLTUInt32_t clusterID) const;
  virtual float GetMaxSignal(AliHLTUInt32_t clusterID) const {return GetQMax(clusterID);}
  virtual float GetPhi(AliHLTUInt32_t clusterID) const;

  /// add input block to the collection
  virtual int AddInputBlock(const AliHLTComponentBlockData* pDesc);

  virtual int PopulateAccessGrid(AliHLTSpacePointPropertyGrid* pGrid, AliHLTUInt32_t mask) const;
  int PopulateAccessGrid(AliHLTSpacePointPropertyGrid* pGrid, AliHLTTPCRawClusterData * pDecoder, int slice, int partition) const;
  virtual const AliHLTSpacePointPropertyGrid* GetSpacePointPropertyGrid(AliHLTUInt32_t mask) const;
  virtual int SetSpacePointPropertyGrid(AliHLTUInt32_t mask, AliHLTSpacePointPropertyGrid*);

  /// clear the object and reset pointer references
  virtual void Clear(Option_t * option ="");

  /// print information
  virtual void Print(ostream& out, Option_t *option="") const;

  /// create a collection of clusters for a space point mask
  virtual AliHLTSpacePointContainer* SelectByMask(AliHLTUInt32_t mask, bool bAlloc=false) const;

  /// create a collection of clusters for a specific track
  virtual AliHLTSpacePointContainer* SelectByTrack(int trackId, bool bAlloc=false) const;

  /// create a collection of clusters for a specific MC track
  virtual AliHLTSpacePointContainer* SelectByMC(int mcId, bool bAlloc=false) const;

  /// create a collection of all used clusters
  virtual AliHLTSpacePointContainer* UsedClusters(bool bAlloc=false) const;

  /// create a collection of all unused clusters
  virtual AliHLTSpacePointContainer* UnusedClusters(bool bAlloc=false) const;

  virtual int MarkUsed(const AliHLTUInt32_t* clusterIDs, int arraySize);
  virtual int SetTrackID(int trackID, const AliHLTUInt32_t* clusterIDs, int arraySize);
  virtual int GetTrackID(AliHLTUInt32_t clusterID) const;
  virtual int SetMCID(int clusterID, const AliHLTUInt32_t* clusterIDs, int arraySize);

  virtual int Write(AliHLTUInt8_t* outputPtr, AliHLTUInt32_t size,
		    vector<AliHLTComponentBlockData>& outputBlocks,
		    AliHLTDataDeflater* pDeflater,
		    const char* option="") const;
  virtual int Write(AliHLTUInt8_t* outputPtr, AliHLTUInt32_t size, AliHLTUInt32_t offset,
		    vector<AliHLTComponentBlockData>& outputBlocks,
		    AliHLTDataDeflater* pDeflater,
		    const char* option="") const;

  int WriteSorted(AliHLTUInt8_t* outputPtr, AliHLTUInt32_t size, AliHLTUInt32_t offset,
		  vector<AliHLTComponentBlockData>& outputBlocks,
		  AliHLTDataDeflater* pDeflater,
		  const char* option="") const;

  int WriteSorted(AliHLTUInt8_t* outputPtr, AliHLTUInt32_t size, AliHLTUInt32_t offset,
		  AliHLTTPCRawClusterData* pDecoder, AliHLTSpacePointPropertyGrid* pGrid,
		  AliHLTUInt32_t mask,
		  vector<AliHLTComponentBlockData>&  outputBlocks,
		  AliHLTDataDeflater* pDeflater,
		  const char* option) const;

  /// allocate index grid, one single point to define the dimensions
  static AliHLTSpacePointPropertyGrid* AllocateIndexGrid();

  class AliHLTTPCRawSpacePointProperties {
  public:
    AliHLTTPCRawSpacePointProperties();
    AliHLTTPCRawSpacePointProperties(const AliHLTTPCRawCluster* pCluster);
    AliHLTTPCRawSpacePointProperties(const AliHLTTPCRawSpacePointProperties& src);
    AliHLTTPCRawSpacePointProperties& operator=(const AliHLTTPCRawSpacePointProperties& src);

    ~AliHLTTPCRawSpacePointProperties() {}

    const AliHLTTPCRawCluster* GetCluster() const {return fpCluster;}
    bool IsUsed() const {return fUsed;}
    void MarkUsed(bool used=true) {fUsed=used;}
    int GetTrackId() const {return fTrackId;}
    void SetTrackId(int trackId) {fTrackId=trackId;}
    int GetMCId() const {return fMCId;}
    void SetMCId(int mcId) {fMCId=mcId;}

    void Print(ostream& out, Option_t *option="") const;

  private:
    const AliHLTTPCRawCluster* fpCluster; //! decoder for data block
    bool fUsed; //! transient
    int fTrackId; //! track id from reconstruction
    int fMCId; //! MC id
  };

  class AliHLTTPCRawSpacePointBlock {
  public:
    AliHLTTPCRawSpacePointBlock(AliHLTUInt32_t id=0, AliHLTTPCRawClusterData *pDecoder=NULL, AliHLTSpacePointPropertyGrid* pGrid=NULL)
      : fDecoder(pDecoder), fGrid(pGrid), fId(id) {}
    AliHLTTPCRawSpacePointBlock(const AliHLTTPCRawSpacePointBlock& s) 
      : fDecoder(s.fDecoder), fGrid(s.fGrid), fId(s.fId) {}
    AliHLTTPCRawSpacePointBlock& operator=(const AliHLTTPCRawSpacePointBlock& s) {
      if (this==&s) return *this;
      fDecoder=s.fDecoder; fGrid=s.fGrid; fId=s.fId; return *this;
    }
    ~AliHLTTPCRawSpacePointBlock() {}

    int GetNofSpacepoints() const {return fDecoder?fDecoder->fCount:0;}
    AliHLTUInt32_t GetId() const {return fId;}
    void SetId(AliHLTUInt32_t id) {fId=id;}
    AliHLTTPCRawClusterData* GetDecoder() const {return fDecoder;} 
    void SetDecoder(AliHLTTPCRawClusterData* pDecoder) {fDecoder=pDecoder;}
    AliHLTSpacePointPropertyGrid* GetGrid() const {return fGrid;}
    void SetGrid(AliHLTSpacePointPropertyGrid* pGrid) {fGrid=pGrid;}

  protected:
  private:
    AliHLTTPCRawClusterData* fDecoder; //!
    AliHLTSpacePointPropertyGrid* fGrid; //!
    AliHLTUInt32_t fId; //!
  };

 protected:

 private:
  /// map of clusters
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties> fClusters; //!

  /// map of cluster id collection for different masks
  std::map<AliHLTUInt32_t, vector<AliHLTUInt32_t>*> fSelections; //!

  /// array of decoders
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointBlock> fBlocks; //!

  /// the one instance for mode single (=1)
  AliHLTTPCRawSpacePointBlock fSingleBlock;

  /// mode
  int fMode; //!

  /// vector of cluster ids for writing
  vector<AliHLTUInt32_t>* fWrittenClusterIds; //!

  ClassDef(AliHLTTPCRawSpacePointContainer, 0)
};

ostream& operator<<(ostream &out, const AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties& p);

#endif
