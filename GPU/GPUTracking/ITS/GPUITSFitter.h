//**************************************************************************\
//* This file is property of and copyright by the ALICE Project            *\
//* ALICE Experiment at CERN, All rights reserved.                         *\
//*                                                                        *\
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *\
//*                  for The ALICE HLT Project.                            *\
//*                                                                        *\
//* Permission to use, copy, modify and distribute this software and its   *\
//* documentation strictly for non-commercial purposes is hereby granted   *\
//* without fee, provided that the above copyright notice appears in all   *\
//* copies and that both the copyright notice and this permission notice   *\
//* appear in the supporting documentation. The authors make no claims     *\
//* about the suitability of this software for any purpose. It is          *\
//* provided "as is" without express or implied warranty.                  *\
//**************************************************************************

/// \file GPUITSFitter.h
/// \author David Rohr, Maximiliano Puccio

#ifndef GPUITSFITTER_H
#define GPUITSFITTER_H

#include "GPUProcessor.h"
#include "GPUITSTrack.h"

namespace o2
{
namespace its
{
class Road;
struct TrackingFrameInfo;
struct Cluster;
class Cell;
} // namespace its
} // namespace o2

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUITSTrack;

class GPUITSFitter : public GPUProcessor
{
 public:
#ifndef GPUCA_GPUCODE
  void InitializeProcessor();
  void RegisterMemoryAllocation();
  void SetMaxData();

  void* SetPointersInput(void* mem);
  void* SetPointersTracks(void* mem);
  void* SetPointersMemory(void* mem);
#endif

  GPUd() o2::its::Road* roads()
  {
    return mRoads;
  }
  GPUd() void SetNumberOfRoads(int v) { mNumberOfRoads = v; }
  GPUd() int NumberOfRoads() { return mNumberOfRoads; }
  GPUd() GPUITSTrack* tracks()
  {
    return mTracks;
  }
  GPUd() GPUAtomic(unsigned int) & NumberOfTracks()
  {
    return mMemory->mNumberOfTracks;
  }
  GPUd() void SetNumberTF(int i, int v) { mNTF[i] = v; }
  GPUd() o2::its::TrackingFrameInfo** trackingFrame()
  {
    return mTF;
  }
  GPUd() const o2::its::Cluster** clusters()
  {
    return mClusterPtrs;
  }
  GPUd() const o2::its::Cell** cells()
  {
    return mCellPtrs;
  }

  void clearMemory();

  struct Memory {
    GPUAtomic(unsigned int) mNumberOfTracks = 0;
  };

 protected:
  int mNumberOfRoads = 0;
  int mNMaxTracks = 0;
  int mNTF[7] = {};
  Memory* mMemory = nullptr;
  o2::its::Road* mRoads = nullptr;
  o2::its::TrackingFrameInfo* mTF[7] = {};
  GPUITSTrack* mTracks = nullptr;

  const o2::its::Cluster* mClusterPtrs[7];
  const o2::its::Cell* mCellPtrs[5];

  short mMemoryResInput = -1;
  short mMemoryResTracks = -1;
  short mMemoryResMemory = -1;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
