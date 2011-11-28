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
#include "AliHLTIndexGrid.h"
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

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// check if a cluster is a candidate for merging
  template<typename T> 
  bool CheckCandidate(int slice, int partition, const T& c) const {
    return CheckCandidate(slice, partition, c.GetPadRow(), c.GetPad(), c.GetTime());
  }

  /// check if a cluster is a candidate for merging
  bool CheckCandidate(int slice,
		      int partition,
		      int partitionrow, // local row in partition
		      float pad,
		      float time) const;

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
			id
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
		   AliHLTUInt32_t id=~AliHLTUInt32_t(0)
		   );

  /// merge clusters
  int Merge();

  /// cleanup
  void Clear();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// helper class to store relevant data for a cluster candidate
  class AliClusterRecord {
  public:
    AliClusterRecord()
      : fSlice(-1), fPartition(-1), fId(~AliHLTUInt32_t(0)), fCluster() {}
    AliClusterRecord(int slice, int partition, AliHLTUInt32_t id, AliHLTTPCRawCluster cluster)
      : fSlice(slice), fPartition(partition), fId(id), fCluster(cluster) {}
    AliClusterRecord(const AliClusterRecord& other)
      : fSlice(other.fSlice), fPartition(other.fPartition), fId(other.fId), fCluster(other.fCluster) {}
    AliClusterRecord& operator=(const AliClusterRecord& other) {
      if (this==&other) return *this;
      this->~AliClusterRecord();
      new (this) AliClusterRecord(other);
      return *this;
    }

    ~AliClusterRecord() {}

    void Clear() {
      fSlice=-1; fPartition=-1; fId=~AliHLTUInt32_t(0);
      fCluster.~AliHLTTPCRawCluster();
      new (&fCluster) AliHLTTPCRawCluster;
    }

    AliClusterRecord& operator=(const AliHLTTPCRawCluster& c) {
      fCluster=c;
      return *this;
    }


    int GetSlice() const {return fSlice;}
    int GetPartition() const {return fPartition;}
    AliHLTUInt32_t GetId() const {return fId;}
    operator AliHLTTPCRawCluster() const {return fCluster;}
    const AliHLTTPCRawCluster& GetCluster() const {return fCluster;}

  private:
    int fSlice; //!
    int fPartition; //!
    AliHLTUInt32_t fId; //!
    AliHLTTPCRawCluster fCluster; //!
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

    AliClusterRecord operator*() {return *fIter;}

    // prefix operators
    iterator& operator++() {
      if (!fArray || fIter==fArray->end()) return *this;
      while ((++fIter)!=fArray->end()) {
	if (fIter->GetCluster().GetCharge()!=0 ||
	    fIter->GetCluster().GetQMax()!=0) {
	  break;
	} else {
	  continue;
	}
      }	     
      return *this;
    }
    iterator& operator--() {
      if (!fArray) return *this;
      while (fIter!=fArray->begin()) {
	--fIter;
	if (fIter->GetCluster().GetCharge()!=0 ||
	    fIter->GetCluster().GetQMax()!=0) {
	  break;
	}
      }
      return *this;
    }

    // postfix operators
    iterator operator++(int) {iterator i(*this); this->operator++(); return i;}
    iterator operator--(int) {iterator i(*this); this->operator--(); return i;}

    iterator& operator+=(int step) {
      if (!fArray) return *this;
      while ((++fIter)!=fArray->end() && step-->0) {}
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
    fEnd=fIter; fEnd+=fClusters.size();
    while (fIter!=fEnd &&
	   (*fIter).GetCluster().GetCharge()==0 &&
	   (*fIter).GetCluster().GetQMax()==0) {
      // skip empty (merged) clusters
      fIter++;
    }

    return fIter;
  }

  /// iterator function, end marker
  iterator& end() {
    return fEnd;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// index grid with coordinates slice, slicerow and time for position of cluster in array
  typedef AliHLTIndexGrid<int, AliHLTUInt32_t> AliSortedClusters;

 protected:

 private:
  /// copy constructor
  AliHLTTPCHWClusterMerger(const AliHLTTPCHWClusterMerger&);
  /// assignment operator
  AliHLTTPCHWClusterMerger& operator=(const AliHLTTPCHWClusterMerger&);

  int FillIndex();

  vector<AliClusterRecord> fClusters; //! array of candidates
  vector<AliHLTUInt32_t> fRemovedClusterIds; //! array of removed clusters by id
  AliSortedClusters fIndex; //! cluster index in slice, row and time
  iterator fIter; //!
  iterator fEnd; //!

  ClassDef(AliHLTTPCHWClusterMerger, 0)
};

#endif //ALIHLTTPCHWCLUSTERMERGER_H
