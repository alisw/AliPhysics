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

/// \file GPUReconstructionIncludes.h
/// \author David Rohr

#ifndef GPURECONSTRUCTIONINCLUDES_H
#define GPURECONSTRUCTIONINCLUDES_H

// Disable assertions since they produce errors in GPU Code
#ifdef assert
#undef assert
#endif
#define assert(param)

#ifndef WIN32
#include <sys/syscall.h>
#include <semaphore.h>
#include <fcntl.h>
#include <sched.h>
#endif

#include "GPUDef.h"
#include "GPULogging.h"

#include <iostream>
#include <fstream>

#if defined(GPUCA_ALIROOT_LIB) && !defined(GPUCA_GPULIBRARY)
#include "AliHLTDefinitions.h"
#include "AliHLTSystem.h"
#endif

#include "GPUReconstructionIncludesITS.h"

#define RANDOM_ERROR
//#define RANDOM_ERROR || rand() % 500 == 1

#define GPUCA_GPUReconstructionUpdateDefailts()                                        \
  if (mProcessingSettings.trackletConstructorInPipeline < 0) {                         \
    mProcessingSettings.trackletConstructorInPipeline = GPUCA_CONSTRUCTOR_IN_PIPELINE; \
  }                                                                                    \
  if (mProcessingSettings.trackletSelectorInPipeline < 0) {                            \
    mProcessingSettings.trackletSelectorInPipeline = GPUCA_SELECTOR_IN_PIPELINE;       \
  }                                                                                    \
  if (mProcessingSettings.trackletSelectorSlices < 0) {                                \
    mProcessingSettings.trackletSelectorSlices = GPUCA_TRACKLET_SELECTOR_SLICE_COUNT;  \
  }                                                                                    \
  if (mProcessingSettings.alternateBorderSort < 0) {                                   \
    mProcessingSettings.alternateBorderSort = GPUCA_ALTERNATE_BORDER_SORT;             \
  }                                                                                    \
  if (mProcessingSettings.mergerSortTracks < 0) {                                      \
    mProcessingSettings.mergerSortTracks = GPUCA_SORT_BEFORE_FIT;                      \
  }                                                                                    \
  if (param().rec.loopInterpolationInExtraPass < 0) {                                  \
    param().rec.loopInterpolationInExtraPass = GPUCA_MERGER_SPLIT_LOOP_INTERPOLATION;  \
  }                                                                                    \
  if (mProcessingSettings.tpcCompressionGatherModeKernel < 0) {                        \
    mProcessingSettings.tpcCompressionGatherModeKernel = GPUCA_COMP_GATHER_KERNEL;     \
  }

#endif
