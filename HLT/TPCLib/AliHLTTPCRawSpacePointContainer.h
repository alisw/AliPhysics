//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCRAWSPACEPOINTCONTAINER_H
#define ALIHLTTPCRAWSPACEPOINTCONTAINER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCRawSpacePointContainer.h
/// @author Matthias Richter
/// @date   2011-08-01
/// @brief  Helper class for handling of HLT TPC raw cluster data blocks
/// @note   Class is a duplicate of AliHLTTPCRawSpacePointContainer and should
///         be merged with it in a generic way

#include "AliHLTSpacePointContainer.h"
#include "AliHLTTPCRawCluster.h"
#include <map>
using namespace std;

/**
 * @class AliHLTTPCRawSpacePointContainer
 * Handler class for HLT TPCS raw space point data blocks.
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCRawSpacePointContainer : public AliHLTSpacePointContainer
{
 public:
  /// standard constructor
  AliHLTTPCRawSpacePointContainer();
  /// copy constructor
  AliHLTTPCRawSpacePointContainer(const AliHLTTPCRawSpacePointContainer& c);
  /// assignment operator
  AliHLTTPCRawSpacePointContainer& operator=(const AliHLTTPCRawSpacePointContainer& c);
  /// destructor
  ~AliHLTTPCRawSpacePointContainer();

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
  virtual float GetMaxSignal(AliHLTUInt32_t clusterID) const;
  virtual float GetPhi(AliHLTUInt32_t clusterID) const;

  /// add input block to the collection
  virtual int AddInputBlock(const AliHLTComponentBlockData* pDesc);

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

  class AliHLTTPCRawSpacePointProperties {
  public:
    AliHLTTPCRawSpacePointProperties();
    AliHLTTPCRawSpacePointProperties(const AliHLTTPCRawCluster* pCluster);
    AliHLTTPCRawSpacePointProperties(const AliHLTTPCRawSpacePointProperties& src);
    AliHLTTPCRawSpacePointProperties& operator=(const AliHLTTPCRawSpacePointProperties& src);

    ~AliHLTTPCRawSpacePointProperties() {}

    const AliHLTTPCRawCluster* Data() const {return fCluster;}
    bool IsUsed() const {return fUsed;}
    void MarkUsed(bool used=true) {fUsed=used;}
    int GetTrackId() const {return fTrackId;}
    void SetTrackId(int trackId) {fTrackId=trackId;}
    int GetMCId() const {return fMCId;}
    void SetMCId(int mcId) {fMCId=mcId;}

    void Print(ostream& out, Option_t *option="") const;

  private:
    const AliHLTTPCRawCluster* fCluster; //! transient
    bool fUsed; //! transient
    int fTrackId; //! track id from reconstruction
    int fMCId; //! MC id
  };

 protected:

 private:
  /// map of clusters
  std::map<AliHLTUInt32_t, AliHLTTPCRawSpacePointProperties> fClusters; //!

  /// map of cluster id collection for different masks
  std::map<AliHLTUInt32_t, vector<AliHLTUInt32_t>*> fSelections; //!

  ClassDef(AliHLTTPCRawSpacePointContainer, 0)
};

ostream& operator<<(ostream &out, const AliHLTTPCRawSpacePointContainer::AliHLTTPCRawSpacePointProperties& p);

#endif
