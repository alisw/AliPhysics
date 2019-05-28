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

/// \file GPUConstantMem.h
/// \author David Rohr

#ifndef GPUCONSTANTMEM_H
#define GPUCONSTANTMEM_H

#include "GPUTPCTracker.h"
#include "GPUParam.h"

#if defined(GPUCA_NOCOMPAT_ALLCINT) && (!defined(GPUCA_GPULIBRARY) || !defined(GPUCA_ALIROOT_LIB))
#include "GPUTPCConvert.h"
#include "GPUTPCCompression.h"
#include "GPUTPCGMMerger.h"
#include "GPUITSFitter.h"
#include "GPUTRDTracker.h"
#else
namespace GPUCA_NAMESPACE
{
namespace gpu
{
class GPUTPCGMMerger
{
};
class GPUITSFitter
{
};
class GPUTRDTracker
{
  void SetMaxData() {}
};
class GPUTPCConvert
{
};
class GPUTPCCompression
{
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE
#endif

namespace GPUCA_NAMESPACE
{
namespace gpu
{
MEM_CLASS_PRE()
struct GPUConstantMem {
  MEM_CONSTANT(GPUParam)
  param;
  MEM_GLOBAL(GPUTPCTracker)
  tpcTrackers[GPUCA_NSLICES];
  GPUTPCConvert tpcConverter;
  GPUTPCCompression tpcCompressor;
  GPUTPCGMMerger tpcMerger;
  GPUTRDTracker trdTracker;
  GPUITSFitter itsFitter;
};
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
