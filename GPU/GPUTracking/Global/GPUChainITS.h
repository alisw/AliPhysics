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

/// \file GPUChainITS.h
/// \author David Rohr

#ifndef GPUCHAINITS_H
#define GPUCHAINITS_H

#include "GPUChain.h"
namespace o2
{
namespace its
{
class Cluster;
class Road;
class Cell;
class TrackingFrameInfo;
class TrackITSExt;
} // namespace its
} // namespace o2

namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUChainITS : public GPUChain
{
  friend class GPUReconstruction;

 public:
  ~GPUChainITS() override;
  void RegisterPermanentMemoryAndProcessors() override;
  void RegisterGPUProcessors() override;
  int Init() override;
  int PrepareEvent() override;
  int Finalize() override;
  int RunChain() override;
  void MemorySize(size_t& gpuMem, size_t& pageLockedHostMem) override;

  int PrepareAndRunITSTrackFit(std::vector<o2::its::Road>& roads, std::array<const o2::its::Cluster*, 7> clusters, std::array<const o2::its::Cell*, 5> cells, const std::array<std::vector<o2::its::TrackingFrameInfo>, 7>& tf, std::vector<o2::its::TrackITSExt>& tracks);
  int RunITSTrackFit(std::vector<o2::its::Road>& roads, std::array<const o2::its::Cluster*, 7> clusters, std::array<const o2::its::Cell*, 5> cells, const std::array<std::vector<o2::its::TrackingFrameInfo>, 7>& tf, std::vector<o2::its::TrackITSExt>& tracks);

  o2::its::TrackerTraits* GetITSTrackerTraits() { return mITSTrackerTraits.get(); }
  o2::its::VertexerTraits* GetITSVertexerTraits() { return mITSVertexerTraits.get(); }

 protected:
  GPUChainITS(GPUReconstruction* rec);
  std::unique_ptr<o2::its::TrackerTraits> mITSTrackerTraits;
  std::unique_ptr<o2::its::VertexerTraits> mITSVertexerTraits;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
