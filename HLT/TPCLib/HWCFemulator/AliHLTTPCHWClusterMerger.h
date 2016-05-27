//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCHWCLUSTERMERGER_H
#define ALIHLTTPCHWCLUSTERMERGER_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTPCHWClusterMerger.h
//  @author Matthias Richter, Sergey Gorbunov
//  @date   2011-11-25
//  @brief  Merger class for HLT TPC Hardware clusters
//          Handles merging of branch border clusters

#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCClusterMCData.h"
#include "AliHLTLogging.h"
#include <vector>
#include "TObject.h"

/**
 * @class AliHLTTPCHWClusterMerger
 *
 * @ingroup alihlt_base
 */
class AliHLTTPCHWClusterMerger : public AliHLTLogging
{
 public:
  /// standard constructor
  AliHLTTPCHWClusterMerger();
  /// destructor
  ~AliHLTTPCHWClusterMerger();

  void Init();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// check if a cluster is a candidate for merging
  template<typename T> 
  bool CheckCandidate(int slice, int partition, const T& c) {
    return CheckCandidate(slice, partition, c.GetPadRow(), c.GetPad(), c.GetTime(), c.GetSigmaPad2() );
  }

  /// check if a cluster is a candidate for merging
  bool CheckCandidate(int slice,
		      int partition,
		      int partitionrow, // local row in partition
		      float pad,
		      float time,
		      float sigmaPad2) ;

  /// cache cluster for later merging
  template<typename T> 
  int AddCandidate(int slice,
		   int partition,
		   AliHLTUInt32_t id,
		   const T& c) {
    return AddCandidate(slice,
			partition,
			c.GetPadRow(),
			c.GetPad(),
			c.GetTime(),
			c.GetSigmaY2(),
			c.GetSigmaZ2(),
			c.GetCharge(),
			c.GetQMax(),
			c.GetFlags(),
			id
			);
  }

  // cache cluster for later merging
  template<typename T> 
  int AddCandidate(int slice,
		   int partition,
		   AliHLTUInt32_t id,
		   const T& c,
		   const AliHLTTPCClusterMCLabel *mc) {
    return AddCandidate(slice,
			partition,
			c.GetPadRow(),
			c.GetPad(),
			c.GetTime(),
			c.GetSigmaPad2(),
			c.GetSigmaTime2(),
			c.GetCharge(),
			c.GetQMax(),
			c.GetFlags(),
			id,
			mc
			);
  }

  /// cache cluster for later merging
  int AddCandidate(int slice,
		   int partition,
		   short partitionrow, // local row in the partition
		   float pad,
		   float time,
		   float sigmaY2,
		   float sigmaZ2,
		   unsigned short charge,
		   unsigned short qmax,
		   unsigned short flags,
		   AliHLTUInt32_t id=~AliHLTUInt32_t(0),
		   const AliHLTTPCClusterMCLabel *mc=NULL
		   );

  /// merge clusters
  int Merge();

  /// cleanup
  void Clear();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// helper class to store relevant data for a cluster at border
  struct AliBorderRecord {
    AliHLTInt64_t fClusterRecordID;
    AliHLTUInt32_t fTimeBin;    
  };
 
  
  /// helper class to store relevant data for a cluster candidate
  class AliClusterRecord {
  public:
    AliClusterRecord()
      : fSlice(-1), fPartition(-1), fBorder(-1), fMergedFlag(-1), fId(~AliHLTUInt32_t(0)), fCluster(), fMC() {}
    AliClusterRecord(int slice, int partition, int border,int merged, AliHLTUInt32_t id, const AliHLTTPCRawCluster &cluster)
      : fSlice(slice), fPartition(partition), fBorder(border), fMergedFlag(merged), fId(id), fCluster(cluster), fMC() {}
    AliClusterRecord(int slice, int partition, int border,int merged, AliHLTUInt32_t id, const AliHLTTPCRawCluster &cluster, const AliHLTTPCClusterMCLabel &mc)
      : fSlice(slice), fPartition(partition), fBorder(border), fMergedFlag(merged), fId(id), fCluster(cluster), fMC(mc) {}

    AliClusterRecord(const AliClusterRecord& other)
      : fSlice(other.fSlice), fPartition(other.fPartition), fBorder(other.fBorder), fMergedFlag(other.fMergedFlag), fId(other.fId), fCluster(other.fCluster), fMC(other.fMC) {}
    AliClusterRecord& operator=(const AliClusterRecord& other) {
      if (this==&other) return *this;
      this->~AliClusterRecord();
      new (this) AliClusterRecord(other);
      return *this;
    }

    ~AliClusterRecord() {}    
    
    AliClusterRecord& operator=(const AliHLTTPCRawCluster& c) {
      fCluster=c;
      return *this;
    }

    int IsMergedTo() const { return fMergedFlag; }
    int GetSlice() const {return fSlice;}
    int GetBorder() const {return fBorder;}
    int GetPartition() const {return fPartition;}
    AliHLTUInt32_t GetId() const {return fId;}
    operator AliHLTTPCRawCluster() const {return fCluster;}
    const AliHLTTPCRawCluster& GetCluster() const {return fCluster;}
    const AliHLTTPCClusterMCLabel& GetMCLabel() const {return fMC;}
    void SetMergedTo( int ind ){ fMergedFlag = ind;}
    AliHLTTPCRawCluster &Cluster(){ return fCluster; }
    AliHLTTPCClusterMCLabel& MCLabel(){ return fMC; }
  private:
    int fSlice; //!
    int fPartition; //!
    int fBorder; //!
    int fMergedFlag; //!
    AliHLTUInt32_t fId; //!
    AliHLTTPCRawCluster fCluster; //!
    AliHLTTPCClusterMCLabel fMC; //!
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// iterator class to access merged and remaining clusters
  class iterator {
  public:
    iterator() : fArray(NULL), fIter() {}
    iterator(vector<AliClusterRecord>* pArray) : fArray(pArray), fIter() {if (fArray) fIter=fArray->begin();}
    iterator(const iterator& other) : fArray(other.fArray), fIter(other.fIter) {}
    iterator& operator=(const iterator& other) {
      if (this==&other) return *this;
      fArray=other.fArray; fIter=other.fIter; return *this;
    }
    ~iterator() {}

    AliClusterRecord& operator*() {return *fIter;}

    // prefix operators
    iterator& operator++() {
      if (!fArray || fIter==fArray->end()) return *this;
      while ((++fIter)!=fArray->end()) {
	if ( fIter->IsMergedTo()<0 ) break;	
      }	     
      return *this;
    }
    iterator& operator--() {
      if (!fArray) return *this;
      while (fIter!=fArray->begin()) {
	--fIter;
	if ( fIter->IsMergedTo()<0 ) break;	
      }
      return *this;
    }

    // postfix operators
    iterator operator++(int) {iterator i(*this); this->operator++(); return i;}
    iterator operator--(int) {iterator i(*this); this->operator--(); return i;}

    iterator& operator+=(int step) {
      if (!fArray) return *this;
      while ((fIter)!=fArray->end() && step-->0) {++fIter;}
      return *this;
    }

    bool operator==(const iterator& other) {
      return (other.fArray!=NULL && fArray!=NULL && other.fIter==fIter);
    }

    bool operator!=(const iterator& other) {
      return (other.fArray!=NULL && fArray!=NULL && other.fIter!=fIter);
    }

  protected:
  private:
    vector<AliClusterRecord>* fArray; //!
    vector<AliClusterRecord>::iterator fIter; //!
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// iterator function, start iteration
  iterator& begin() {
    fIter.~iterator();
    new (&fIter) iterator(&fClusters);
    fEnd=fIter;  fEnd+=fClusters.size();    
    // skip empty (merged) clusters
    while (fIter!=fEnd && ( (*fIter).IsMergedTo()>=0) ) {
      fIter++;
    }
    return fIter;
  }

  /// iterator function, end marker
  iterator& end() {
    return fEnd;
  }

  const vector<AliHLTTPCHWClusterMerger::AliClusterRecord> &GetRecords(){ return fClusters; }
  static int GetNSlices(){ return fkNSlices; }
  int GetNBorders() const { return fNBorders; }
  int GetBorderNClusters( int ib ) const { return fBorderNClusters[ib]; }
  int GetBorderFirstCluster( int ib ) const { return fBorderFirstCluster[ib]; }
  const AliHLTTPCHWClusterMerger::AliBorderRecord *GetBorderClusters() const { return fBorderClusters;}

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

 protected:

 private:
  /// copy constructor
  AliHLTTPCHWClusterMerger(const AliHLTTPCHWClusterMerger&);
  /// assignment operator
  AliHLTTPCHWClusterMerger& operator=(const AliHLTTPCHWClusterMerger&);

  int FillIndex();
  static bool CompareTime( const AliBorderRecord &b1, const AliBorderRecord &b2){
    return b1.fTimeBin > b2.fTimeBin;
  }
 
  static bool CompareMCWeights(const AliHLTTPCClusterMCWeight &a, const AliHLTTPCClusterMCWeight &b){
    return a.fWeight > b.fWeight;
  }
  static bool CompareMCLabels(const AliHLTTPCClusterMCWeight &a, const AliHLTTPCClusterMCWeight &b){
    return a.fMCID < b.fMCID;
  }

  AliHLTInt16_t *fMapping;//!
  int fNRows;//!
  int fNRowPads;//!
  int fNBorders;//!
  AliHLTFloat32_t *fBorders; //!
  int *fBorderNClusters; //!
  int *fBorderFirstCluster; //!
  AliBorderRecord *fBorderClusters;
  int fBorderNClustersTotal; //!

  vector<AliClusterRecord> fClusters; //! array of candidates
  vector<AliHLTUInt32_t> fRemovedClusterIds; //! array of removed clusters by id
  iterator fIter; //!
  iterator fEnd; //!
  static const int fkMergeWidth = 3;
  static const int fkNSlices = 36;
  static const int fkMergeTimeWindow = 3;
  ClassDef(AliHLTTPCHWClusterMerger, 0)
};

#endif //ALIHLTTPCHWCLUSTERMERGER_H
