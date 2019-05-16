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

/// \file GPUTPCSliceOutput.h
/// \author Sergey Gorbunov, Ivan Kisel, David Rohr

#ifndef GPUTPCSLICEOUTPUT_H
#define GPUTPCSLICEOUTPUT_H

#include "GPUTPCDef.h"
#ifndef GPUCA_GPUCODE
#include "GPUTPCSliceOutTrack.h" //Breaks OpenCL 1.2 since GPUTPCSliceOutput does not know the address space of the track
#endif

#if !defined(__OPENCL__)
#include <cstdlib>
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct GPUOutputControl;

/**
 * @class GPUTPCSliceOutput
 *
 * GPUTPCSliceOutput class is used to store the output of GPUTPCTracker{Component}
 * and transport the output to GPUTPCGBMerger{Component}
 *
 * The class contains all the necessary information about TPC tracks, reconstructed in one slice.
 * This includes the reconstructed track parameters and some compressed information
 * about the assigned clusters: clusterId, position and amplitude.
 *
 */
class GPUTPCSliceOutput
{
 public:
#if !defined(GPUCA_GPUCODE_DEVICE)
  GPUhd() unsigned int NTracks() const
  {
    return mNTracks;
  }
  GPUhd() unsigned int NLocalTracks() const { return mNLocalTracks; }
  GPUhd() unsigned int NTrackClusters() const { return mNTrackClusters; }
#ifndef GPUCA_GPUCODE
  GPUhd() const GPUTPCSliceOutTrack* GetFirstTrack() const
  {
    return mMemory;
  }
  GPUhd() GPUTPCSliceOutTrack* FirstTrack()
  {
    return mMemory;
  }
#endif
  GPUhd() size_t Size() const
  {
    return (mMemorySize);
  }

  static unsigned int EstimateSize(unsigned int nOfTracks, unsigned int nOfTrackClusters);
  static void Allocate(GPUTPCSliceOutput*& ptrOutput, int nTracks, int nTrackHits, GPUOutputControl* outputControl, void*& internalMemory);

  GPUhd() void SetNTracks(unsigned int v) { mNTracks = v; }
  GPUhd() void SetNLocalTracks(unsigned int v) { mNLocalTracks = v; }
  GPUhd() void SetNTrackClusters(unsigned int v) { mNTrackClusters = v; }

 private:
  GPUTPCSliceOutput() CON_DELETE;                                    // NOLINT: Must be private or ROOT tries to use them!
  ~GPUTPCSliceOutput() CON_DELETE;                                   // NOLINT
  GPUTPCSliceOutput(const GPUTPCSliceOutput&) CON_DELETE;            // NOLINT
  GPUTPCSliceOutput& operator=(const GPUTPCSliceOutput&) CON_DELETE; // NOLINT

  GPUh() void SetMemorySize(size_t val) { mMemorySize = val; }

  unsigned int mNTracks; // number of reconstructed tracks
  unsigned int mNLocalTracks;
  unsigned int mNTrackClusters; // total number of track clusters
  size_t mMemorySize;           // Amount of memory really used

// Must be last element of this class, user has to make sure to allocate anough memory consecutive to class memory!
// This way the whole Slice Output is one consecutive Memory Segment
#ifndef GPUCA_GPUCODE
#ifdef __OPENCL__
  GPUTPCSliceOutTrack mMemory[1]; // the memory where the pointers above point into
#else
  GPUTPCSliceOutTrack mMemory[0]; // the memory where the pointers above point into
#endif
#endif
#endif
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE
#endif
