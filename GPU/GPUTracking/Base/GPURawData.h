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

/// \file GPURawData.h
/// \author David Rohr

#ifndef O2_GPU_RAW_DATA_H
#define O2_GPU_RAW_DATA_H

// Raw data parser is not accessible from GPU, therefore we use this header to wrap direct access to the current RDH
// Since OpenCL currently doesn't support bit fields, we have to access the members directly

#include "GPUCommonDef.h"
#ifndef __OPENCL__
#include "Headers/RAWDataHeader.h"
#include "DetectorsRaw/RDHUtils.h"
#else
namespace o2
{
namespace header
{
struct RAWDataHeader {
  union {
    unsigned int words[8];
  };
};
} // namespace header
} // namespace o2
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
typedef o2::header::RAWDataHeader RAWDataHeaderGPU;

class GPURawDataUtils
{
 public:
  static GPUd() unsigned int getOrbit(const RAWDataHeaderGPU* rdh);
  static GPUd() unsigned int getBC(const RAWDataHeaderGPU* rdh);
  static GPUd() unsigned int getSize(const RAWDataHeaderGPU* rdh);
};

GPUdi() unsigned int GPURawDataUtils::getOrbit(const RAWDataHeaderGPU* rdh)
{
#ifndef __OPENCL__
  return o2::raw::RDHUtils::getHeartBeatOrbit(*rdh);
#else
  return (rdh->words[2] >> 32);           // TODO: Ad-hoc implementation for OpenCL, RDHV4, to be moved to RDHUtils
#endif
}

GPUdi() unsigned int GPURawDataUtils::getBC(const RAWDataHeaderGPU* rdh)
{
#ifndef __OPENCL__
  return o2::raw::RDHUtils::getHeartBeatBC(*rdh);
#else
  return (rdh->words[2] & 0xFFF);         // TODO: Ad-hoc implementation for OpenCL, RDHV4, to be moved to RDHUtils
#endif
}

GPUdi() unsigned int GPURawDataUtils::getSize(const RAWDataHeaderGPU* rdh)
{
#ifndef __OPENCL__
  return o2::raw::RDHUtils::getMemorySize(*rdh);
#else
  return ((rdh->words[1] >> 16) & 0xFFFF); // TODO: Ad-hoc implementation for OpenCL, RDHV4, to be moved to RDHUtils
#endif
}

} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
