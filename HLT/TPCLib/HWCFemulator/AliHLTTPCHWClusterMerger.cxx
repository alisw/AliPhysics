// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sergey Gorbunov                                       *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

//  @file   AliHLTTPCHWClusterMerger.cxx
//  @author Matthias Richter, Sergey Gorbunov
//  @date   2011-11-25
//  @brief  Merger class for HLT TPC Hardware clusters
//          Handles merging of branch border clusters

#include "AliHLTTPCHWClusterMerger.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCSpacePointData.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCHWClusterMerger)

const int gkIndexGridTimeStep=10;
AliHLTTPCHWClusterMerger::AliHLTTPCHWClusterMerger()
  : AliHLTLogging()
  , fClusters()
  , fRemovedClusterIds()
  , fIndex(AliHLTTPCTransform::GetNSlice(), 1,
	   AliHLTTPCTransform::GetNRows(), 1,
	   AliHLTTPCTransform::GetNTimeBins(), gkIndexGridTimeStep)
  , fIter()
  , fEnd()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  // constructor
}

AliHLTTPCHWClusterMerger::~AliHLTTPCHWClusterMerger()
{
  // destructor
}

bool AliHLTTPCHWClusterMerger::CheckCandidate(int /*slice*/, int partition, int partitionrow, float pad, float /*time*/) const
{
  /// check cluster if it is a candidate for merging
  int slicerow=partitionrow+AliHLTTPCTransform::GetFirstRow(partition);

  // FIXME: implement the logic here
  if (TMath::Abs(pad-AliHLTTPCTransform::GetNPads(slicerow)/2)<2) {
    return true;
  }
  return false;
}

int AliHLTTPCHWClusterMerger::AddCandidate(int slice,
					   int partition,
					   short partitionrow,
					   float pad,
					   float time,
					   float sigmaY2,
					   float sigmaZ2,
					   unsigned short charge,
					   unsigned short qmax,
					   AliHLTUInt32_t id
					   )
{
  /// add a candidate for merging and register in the index grid
  int slicerow=partitionrow+AliHLTTPCTransform::GetFirstRow(partition);
  if (id!=~AliHLTUInt32_t(0)) {
    if (slice<0) {
      slice=AliHLTTPCSpacePointData::GetSlice(id);
    } else if ((unsigned)slice!=AliHLTTPCSpacePointData::GetSlice(id)) {
      HLTError("cluster id 0x%08x is not consistent with specified slice %d", id, slice);
    }
    if (partition<0) {
      partition=AliHLTTPCSpacePointData::GetPatch(id);
    } else if ((unsigned)partition!=AliHLTTPCSpacePointData::GetPatch(id)) {
      HLTError("cluster id 0x%08x is not consistent with specified partition %d", id, partition);
    }
  }
  fClusters.push_back(AliClusterRecord(slice, partition, id,
				       AliHLTTPCRawCluster(partitionrow, pad, time, sigmaY2, sigmaZ2, charge, qmax)));
  fIndex.CountSpacePoint(slice, slicerow, (int)time);

  return fClusters.size();
}

int AliHLTTPCHWClusterMerger::Merge()
{
  /// merge clusters
  /// first all candidates are filled into the index grid, looping over all clusters
  /// in the grid automatically results in time ordered clusters per row (ordering with
  /// respect to cells, within a cell two clusters can be reversed)
  int iResult=0;
  int count=0;
  if ((iResult=FillIndex())<0) {
    return iResult;
  }

  vector<int> cellClusters;
  for (AliSortedClusters::iterator& cli=fIndex.begin();
       cli!=fIndex.end(); cli++) {
    if (cellClusters.size()==0) {
      // no partner for merging, put this one on the stack and continue
      cellClusters.push_back(*cli);
      continue;
    }
    AliHLTTPCRawCluster c1=fClusters[*cli].GetCluster();
    int slicerow1=c1.GetPadRow()+AliHLTTPCTransform::GetFirstRow(fClusters[*cli].GetPartition());
    vector<int>::iterator partner=cellClusters.begin();
    bool bMerged=false;
    while (partner!=cellClusters.end()) {
      AliHLTTPCRawCluster c2=fClusters[*partner].GetCluster();
      int slicerow2=c2.GetPadRow()+AliHLTTPCTransform::GetFirstRow(fClusters[*partner].GetPartition());
      if ((fClusters[*partner].GetSlice()!=fClusters[*cli].GetSlice()) ||
	  (slicerow2!=slicerow1) ||
	  (c2.GetTime()+gkIndexGridTimeStep<c1.GetTime())) {
	// already out of range to be a partner for merging
	// remove and go to next
	partner=cellClusters.erase(partner);
	continue;
      }

      // check if two clusters need to be merged
      // FIXME: implement correct logic
      // - clusters on different sides of the branch border
      // - pad difference within limit
      // - time difference within limit (specify limit, studies needed)
      // - cluster shape using charge, qmax and sigma
      // - single pad cluster on either one or the other side or on both
      if (false && // the condition below was just for testing
	  TMath::Abs(c2.GetPad()-c1.GetPad())<2 &&
	  TMath::Abs(c2.GetTime()-c1.GetTime())<2) {
	// merge c1 and c2 into c1
        // FIXME: implement merging

	// store c1
	fClusters[*cli]=c1;

	// wipe c2
	fRemovedClusterIds.push_back(fClusters[*partner].GetId());
	HLTDebug("merging %d into %d", *partner, *cli);
	fClusters[*partner].Clear();
	cellClusters.erase(partner);
	bMerged=true;
	count++;
	break;
      }
      partner++;
    }
    if (!bMerged) {
      // did not find any partner, put on stack
      cellClusters.push_back(*cli);
    }
  }
  if (iResult<0) return iResult;
  return count;
}

int AliHLTTPCHWClusterMerger::FillIndex()
{
  /// loop over cached raw clusters and fill index
  for (AliHLTUInt32_t pos=0; pos<fClusters.size(); pos++) {
    int slicerow=fClusters[pos].GetCluster().GetPadRow()+AliHLTTPCTransform::GetFirstRow(fClusters[pos].GetPartition());
    fIndex.AddSpacePoint(pos, fClusters[pos].GetSlice(), slicerow, (int)fClusters[pos].GetCluster().GetTime());
  }
  return 0;
}

void AliHLTTPCHWClusterMerger::Clear()
{
  /// cleanup
  fClusters.clear();
  fRemovedClusterIds.clear();
  fIndex.Clear();
}
