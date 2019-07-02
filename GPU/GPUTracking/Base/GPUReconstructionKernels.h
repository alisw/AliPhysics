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

/// \file GPUReconstructionKernels.h
/// \author David Rohr

#ifndef GPURECONSTRUCTIONKERNELS_H
#define GPURECONSTRUCTIONKERNELS_H

namespace GPUCA_NAMESPACE
{
namespace gpu
{
#ifdef GPUCA_GPURECONSTRUCTIONCPU_DECLONLY
template <>
class GPUReconstructionKernels<GPUReconstructionCPUBackend> : public GPUReconstructionCPUBackend
#define GPUCA_KRNL(...) ;
#define GPUCA_KRNL_CLASS GPUReconstructionCPUBackend
#else
template <class T>
class GPUReconstructionKernels : public T
#define GPUCA_EXPAND(...) __VA_ARGS__
#define GPUCA_KRNL(X) GPUCA_EXPAND X
#define GPUCA_KRNL_CLASS T
#endif
{
 public:
  using krnlRunRange = GPUReconstruction::krnlRunRange;
  using krnlExec = GPUReconstruction::krnlExec;
  using krnlEvent = GPUReconstruction::krnlEvent;
  template <class X, int Y = 0>
  using classArgument = GPUReconstruction::classArgument<X, Y>;

  virtual ~GPUReconstructionKernels() = default; // NOLINT: Do not declare override in template class! AMD hcc will not create the destructor otherwise.
  GPUReconstructionKernels(const GPUSettingsProcessing& cfg) : GPUCA_KRNL_CLASS(cfg) {}

 protected:
  virtual int runKernelImpl(classArgument<GPUTPCNeighboursFinder>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCNeighboursFinder>(x, y, z); }));                         // NOLINT
  virtual int runKernelImpl(classArgument<GPUTPCNeighboursCleaner>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCNeighboursCleaner>(x, y, z); }));                       // NOLINT
  virtual int runKernelImpl(classArgument<GPUTPCStartHitsFinder>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCStartHitsFinder>(x, y, z); }));                           // NOLINT
  virtual int runKernelImpl(classArgument<GPUTPCStartHitsSorter>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCStartHitsSorter>(x, y, z); }));                           // NOLINT
  virtual int runKernelImpl(classArgument<GPUTPCTrackletConstructor>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCTrackletConstructor>(x, y, z); }));                   // NOLINT
  virtual int runKernelImpl(classArgument<GPUTPCTrackletConstructor, 1>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCTrackletConstructor, 1>(x, y, z); }));             // NOLINT
  virtual int runKernelImpl(classArgument<GPUTPCTrackletSelector>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCTrackletSelector>(x, y, z); }));                         // NOLINT
  virtual int runKernelImpl(classArgument<GPUMemClean16>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z, void* ptr, unsigned long size) GPUCA_KRNL(({ return T::template runKernelBackend<GPUMemClean16>(x, y, z, ptr, size); })); // NOLINT
  virtual int runKernelImpl(classArgument<GPUTPCGMMergerTrackFit>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCGMMergerTrackFit>(x, y, z); }));                         // NOLINT
  virtual int runKernelImpl(classArgument<GPUTRDTrackerGPU>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTRDTrackerGPU>(x, y, z); }));                                     // NOLINT
  virtual int runKernelImpl(classArgument<GPUITSFitterKernel>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUITSFitterKernel>(x, y, z); }));                                 // NOLINT
  virtual int runKernelImpl(classArgument<GPUTPCConvertKernel>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCConvertKernel>(x, y, z); }));                               // NOLINT
  virtual int runKernelImpl(classArgument<GPUTPCCompressionKernels, 0>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCCompressionKernels, 0>(x, y, z); }));               // NOLINT
  virtual int runKernelImpl(classArgument<GPUTPCCompressionKernels, 1>, const krnlExec& x, const krnlRunRange& y, const krnlEvent& z) GPUCA_KRNL(({ return T::template runKernelBackend<GPUTPCCompressionKernels, 1>(x, y, z); }));               // NOLINT
};

#undef GPUCA_KRNL
#undef GPUCA_KRNL_CLASS
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
