//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTPCSPACEPOINTCONTAINER_H
#define ALIHLTTPCSPACEPOINTCONTAINER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCSpacePointContainer.h
/// @author Matthias Richter
/// @date   2011-04-29
/// @brief  Helper class for handling of HLT TPC space point data blocks
///

#include "AliHLTSpacePointContainer.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCSpacePointData.h"
#include <map>
using namespace std;

/**
 * @class AliHLTTPCSpacePointContainer
 * Handler class for HLT TPCS space point data blocks.
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCSpacePointContainer : public AliHLTSpacePointContainer
{
 public:
  /// standard constructor
  AliHLTTPCSpacePointContainer();
  /// copy constructor
  AliHLTTPCSpacePointContainer(const AliHLTTPCSpacePointContainer& c);
  /// assignment operator
  AliHLTTPCSpacePointContainer& operator=(const AliHLTTPCSpacePointContainer& c);
  /// destructor
  ~AliHLTTPCSpacePointContainer();

  virtual int GetClusterIDs(vector<AliHLTUInt32_t>& tgt) const;
  virtual float GetX(AliHLTUInt32_t clusterID) const;
  virtual float GetXWidth(AliHLTUInt32_t clusterID) const;
  virtual float GetY(AliHLTUInt32_t clusterID) const;
  virtual float GetYWidth(AliHLTUInt32_t clusterID) const;
  virtual float GetZ(AliHLTUInt32_t clusterID) const;
  virtual float GetZWidth(AliHLTUInt32_t clusterID) const;
  virtual float GetCharge(AliHLTUInt32_t clusterID) const;
  virtual float GetPhi(AliHLTUInt32_t clusterID) const;

  /// add input block to the collection
  virtual int AddInputBlock(const AliHLTComponentBlockData* pDesc);

  /// clear the object and reset pointer references
  virtual void Clear(Option_t * option ="");

  /// print information
  virtual void Print(ostream& out, Option_t *option="") const;

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

  class AliHLTTPCSpacePointProperties {
  public:
    AliHLTTPCSpacePointProperties();
    AliHLTTPCSpacePointProperties(const AliHLTTPCSpacePointData* pCluster);
    AliHLTTPCSpacePointProperties(const AliHLTTPCSpacePointProperties& src);
    AliHLTTPCSpacePointProperties& operator=(const AliHLTTPCSpacePointProperties& src);

    ~AliHLTTPCSpacePointProperties() {}

    const AliHLTTPCSpacePointData* Data() const {return fCluster;}
    bool IsUsed() const {return fUsed;}
    void MarkUsed(bool used=true) {fUsed=used;}
    int GetTrackId() const {return fTrackId;}
    void SetTrackId(int trackId) {fTrackId=trackId;}
    int GetMCId() const {return fMCId;}
    void SetMCId(int mcId) {fMCId=mcId;}

    void Print(ostream& out, Option_t *option="") const;

  private:
    const AliHLTTPCSpacePointData* fCluster; //! transient
    bool fUsed; //! transient
    int fTrackId; //! track id from reconstruction
    int fMCId; //! MC id
  };

 protected:

 private:
  /// map of clusters
  std::map<AliHLTUInt32_t, AliHLTTPCSpacePointProperties> fClusters; //!

  ClassDef(AliHLTTPCSpacePointContainer, 0)
};

ostream& operator<<(ostream &out, const AliHLTTPCSpacePointContainer::AliHLTTPCSpacePointProperties& p);

#endif
