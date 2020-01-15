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

/// \file GPUCommonAlgorithmCUDA.cuh
/// \author David Rohr

#ifndef GPUCOMMONALGORITHMCUDA_CUH
#define GPUCOMMONALGORITHMCUDA_CUH

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <thrust/device_ptr.h>
#pragma GCC diagnostic pop

#include <cuda.h>

namespace GPUCA_NAMESPACE
{
namespace gpu
{

template <class T>
GPUdi() void GPUCommonAlgorithm::sort(T* begin, T* end)
{
  thrust::device_ptr<T> thrustBegin(begin);
  thrust::device_ptr<T> thrustEnd(end);
  thrust::sort(thrust::seq, thrustBegin, thrustEnd);
}

template <class T, class S>
GPUdi() void GPUCommonAlgorithm::sort(T* begin, T* end, const S& comp)
{
  thrust::device_ptr<T> thrustBegin(begin);
  thrust::device_ptr<T> thrustEnd(end);
  thrust::sort(thrust::seq, thrustBegin, thrustEnd, comp);
}

template <class T>
GPUdi() void GPUCommonAlgorithm::sortInBlock(T* begin, T* end)
{
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
    thrust::device_ptr<T> thrustBegin(begin);
    thrust::device_ptr<T> thrustEnd(end);
    thrust::sort(thrust::cuda::par, thrustBegin, thrustEnd);
  }
}

template <class T, class S>
GPUdi() void GPUCommonAlgorithm::sortInBlock(T* begin, T* end, const S& comp)
{
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
    thrust::device_ptr<T> thrustBegin(begin);
    thrust::device_ptr<T> thrustEnd(end);
    thrust::sort(thrust::cuda::par, thrustBegin, thrustEnd, comp);
  }
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
