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

/// \file GPUGeneralKernels.cxx
/// \author David Rohr

#include "GPUGeneralKernels.h"
using namespace GPUCA_NAMESPACE::gpu;

template <>
GPUdii() void GPUMemClean16::Thread<0>(int nBlocks, int nThreads, int iBlock, int iThread, GPUsharedref() MEM_LOCAL(GPUSharedMemory) & GPUrestrict() smem, processorType& GPUrestrict() processors, GPUglobalref() void* ptr, unsigned long size)
{
  const unsigned long stride = get_global_size(0);
  int4 i0;
  i0.x = i0.y = i0.z = i0.w = 0;
  int4* ptra = (int4*)ptr;
  unsigned long len = (size + sizeof(int4) - 1) / sizeof(int4);
  for (unsigned long i = get_global_id(0); i < len; i += stride) {
    ptra[i] = i0;
  }
}
