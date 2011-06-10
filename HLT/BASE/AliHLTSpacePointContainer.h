//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTSPACEPOINTCONTAINER_H
#define ALIHLTSPACEPOINTCONTAINER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTSpacePointContainer.h
/// @author Matthias Richter
/// @date   2011-04-29
/// @brief  Base helper class for handling of HLT space point data blocks
///

#include <vector>
#include <TObject.h>
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"
#include "AliHLTStdIncludes.h"

class TArrayC;
class TH1;

/**
 * @class AliHLTSpacePointContainer
 * Base class of helper classes for space point data blocks.
 * The class implements a couple of interface methods to be commonly used
 * for space point data blocks.
 *
 * @ingroup alihlt_base
 */
class AliHLTSpacePointContainer : public TObject, public AliHLTLogging
{
 public:
  /// standard constructor
  AliHLTSpacePointContainer();
  /// copy constructor
  AliHLTSpacePointContainer(const AliHLTSpacePointContainer&);
  /// assignment operator
  AliHLTSpacePointContainer& operator=(const AliHLTSpacePointContainer&);

  /// destructor
  ~AliHLTSpacePointContainer();

  /// add input block to the collection
  virtual int AddInputBlock(const AliHLTComponentBlockData* pDesc)=0;

  virtual int GetNumberOfSpacePoints() const;
  virtual bool Check(AliHLTUInt32_t clusterID) const;
  virtual int GetClusterIDs(vector<AliHLTUInt32_t>& tgt) const = 0;
  virtual float GetX(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetXWidth(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetY(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetYWidth(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetZ(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetZWidth(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetCharge(AliHLTUInt32_t clusterID) const = 0;
  virtual float GetPhi(AliHLTUInt32_t /*clusterID*/) const {return 0.0;}

  /// create a collection of clusters for a specific track
  virtual AliHLTSpacePointContainer* SelectByTrack(int trackId, bool bAlloc=false) const;

  /// create a collection of clusters for a specific MC track
  virtual AliHLTSpacePointContainer* SelectByMC(int mcId, bool bAlloc=false) const;

  /// create a collection of all used clusters
  virtual AliHLTSpacePointContainer* UsedClusters(bool bAlloc=false) const;

  /// create a collection of all unused clusters
  virtual AliHLTSpacePointContainer* UnusedClusters(bool bAlloc=false) const;

  int MarkUsed(AliHLTUInt32_t clusterID) {return MarkUsed(&clusterID, sizeof(clusterID));}
  virtual int MarkUsed(const AliHLTUInt32_t* clusterIDs, int arraySize);

  int SetTrackID(int trackID, AliHLTUInt32_t clusterID) {
    return SetTrackID(trackID, &clusterID, sizeof(clusterID));
  }
  virtual int SetTrackID(int trackID, const AliHLTUInt32_t* clusterIDs, int arraySize);

  virtual int GetTrackID(AliHLTUInt32_t /*clusterID*/) const {return -1;}

  int SetMCID(int mcID, AliHLTUInt32_t clusterID) {
    return SetMCID(mcID, &clusterID, sizeof(clusterID));
  }
  virtual int SetMCID(int clusterID, const AliHLTUInt32_t* clusterIDs, int arraySize);

  /// add input block from file to collection
  int AddInputBlock(const char* filename, AliHLTComponentDataType dt, unsigned specification);

  /// add input block from list of blank separated files to collection
  int AddInputBlocks(const char* filenames, AliHLTComponentDataType dt);

  /// inherited from TObject: clear the object and reset pointer references
  virtual void Clear(Option_t * /*option*/ ="");

  /// inherited from TObject
  virtual void Print(Option_t *option="") const;

  virtual void Print(ostream& out, Option_t *option="") const;

  void Draw(Option_t *option);

  TH1* DrawProjection(const char* plane) const {
    vector<AliHLTUInt32_t> selection; // empty list -> no selection
    return DrawProjection(plane, selection);
  }

  TH1* DrawProjection(const char* plane, AliHLTUInt32_t specification) const {
    vector<AliHLTUInt32_t> selection; selection.push_back(specification);
    return DrawProjection(plane, selection);
  }

  TH1* DrawProjection(const char* plane, const vector<AliHLTUInt32_t>& selection) const;

 protected:

 private:
  vector<TArrayC*> fBuffers; //! buffers of loaded files

  ClassDef(AliHLTSpacePointContainer, 0)
};

ostream& operator<<(ostream &out, const AliHLTSpacePointContainer& c);

#endif
